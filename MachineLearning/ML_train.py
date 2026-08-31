'''
MAE - mean absolute error
RMSE - root mean squared error
R2 - score of fit quality. Best is 1.0

The model learns:
[emitx_norm, beta_x0, alpha_x0, emit_y_norm, beta_y0, alpha_y0, K1L]
->
[sigx_OTR0X, sigx_OTR1X, sigx_OTR2X, sigx_OTR3X,
sigy_OTR0X, sigy_OTR1X, sigy_OTR2X, sigy_OTR3X]

'''

import sys
import torch
from pathlib import Path
import numpy as np
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupShuffleSplit
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

def screens_folder_name(screens):
    #screens = [str(screen).strip() for screen in (screens or []) if str(screen).strip()]
    #if not screens:
        return "all_screens"
    #return "_".join(screens)

def get_ml_model_file(machine_name, quad_name, screens):
    current_file = Path(__file__).resolve()
    project_root = current_file.parents[1]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    return project_root / "MachineLearning" / str(machine_name) / str(quad_name) / screens_folder_name(screens) / "EM_model.pt"

def get_ml_dataset_file(machine_name, quad_name, screens):
    current_file = Path(__file__).resolve()
    project_root = current_file.parents[1]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    return project_root / "MachineLearning" / str(machine_name) / str(quad_name) / screens_folder_name(screens) / "EM_dataset.npz"

try:
    torch.set_num_threads(1)
except Exception:
    pass

class NeuralNetwork(nn.Module):

    def __init__(self, n_inputs, n_outputs):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_inputs, 128),
            nn.LayerNorm(128),
            nn.SiLU(),
            nn.Dropout(0.05),
            nn.Linear(128, 128),
            nn.SiLU(),
            nn.Dropout(0.05),
            nn.Linear(128, 64),
            nn.SiLU(),
            #nn.Dropout(0.05),
            nn.Linear(64, n_outputs),
        )

    def forward(self, x):
        return self.net(x) # gets x through each layer of NN

class TrainModel:
    def __init__(self, screens, quad_name, machine_name, dataset_file=None, model_file=None, log_callback=None, progress_callback=None, stop_checker=None):
        self.screens = [str(screen) for screen in screens]
        self.quad_name = str(quad_name)
        self.machine_name = str(machine_name)
        if dataset_file is None:
            dataset_file = get_ml_dataset_file(self.machine_name, self.quad_name, self.screens)

        if model_file is None:
            model_file = get_ml_model_file(self.machine_name, self.quad_name, self.screens)

        self.dataset_file = Path(dataset_file)
        self.model_file = Path(model_file)
        self.random_seed = 2137
        self.test_size = 0.2 # 20% goes into test, 80% to training. with 3000 samples, 600 go into test to see if a NN can predict them
        self.batch_size = 256 # packets with 128 samples at the same time
        self.max_epochs = 2000 # maximum numbers of iterations though the dataset
        self.learning_rate = 1e-3 # size of a step during learning
        self.weight_decay = 1e-3 # small penalty for too big weights of a neural network
        self.patience = 30 # if through 120 epochs result doesn't get better, triggers early stopping
        self.log_callback = log_callback
        self.progress_callback = progress_callback
        self.stop_checker = stop_checker
        self.model = None
        self.x_scaler = None
        self.y_scaler = None
        self.param_names = []
        self.sigma_names = []
        self.target_observable = "sigma"
        self.metrics = {}
        self.sample_groups=None

    def log(self, message):
        if callable(self.log_callback):
            self.log_callback(str(message))
        else:
            print(message)

    def should_stop(self):
        return callable(self.stop_checker) and bool(self.stop_checker())

    def emit_progress(self, epoch, total_epochs):
        if callable(self.progress_callback):
            self.progress_callback(int(epoch), int(total_epochs))

    def set_random_seed(self):
        np.random.seed(self.random_seed)
        torch.manual_seed(self.random_seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(self.random_seed) # randomly asigns which sample is for test, which for training

    def get_device(self):
        return torch.device("cpu")

    @staticmethod
    def make_loader(X, Y, batch_size, shuffle):
        X_tensor = torch.tensor(X, dtype=torch.float32)
        Y_tensor = torch.tensor(Y, dtype=torch.float32)
        dataset = TensorDataset(X_tensor, Y_tensor) # makes pairs input -> output
        return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle) # shuffle shuffles samples

    def load_dataset(self):
        if not self.dataset_file.exists():
            raise FileNotFoundError(f"Dataset file not found: {self.dataset_file}")
        data = np.load(self.dataset_file, allow_pickle=True)
        X = np.asarray(data["X"], dtype=float)
        Y = np.asarray(data["Y"], dtype=float)

        if "param_names" in data:
            self.param_names = [str(v) for v in data["param_names"]]
        if "sigma_names" in data:
            self.sigma_names = [str(v) for v in data["sigma_names"]]
        self.target_observable = str(data["target_observable"].item()) if "target_observable" in data.files else "sigma2"
        if self.target_observable not in {"sigma", "sigma2"}:
            raise RuntimeError(f"Unsupported ML dataset target: {self.target_observable}")
        if "screens" in data:
            self.screens = [str(v) for v in data["screens"]]
        if "quad_name" in data:
            self.quad_name = str(data["quad_name"])

        finite = np.all(np.isfinite(X), axis=1) & np.all(np.isfinite(Y), axis=1) # deletes samples with nan, -inf, +inf
        n_k1l_per_twiss_set = 7
        groups = np.arange(X.shape[0],dtype = int) // max(1, n_k1l_per_twiss_set)
        X = X[finite]
        Y = Y[finite]
        self.sample_groups = groups[finite]

        if X.size ==0 or Y.size ==0:
            raise RuntimeError("Dataset contains no valid data samples.")
        return X, Y

    def train(self):
        self.set_random_seed()
        X, Y = self.load_dataset()
        if self.sample_groups is not None and len(np.unique(self.sample_groups)) > 1:
            splitter = GroupShuffleSplit(n_splits=1, test_size=self.test_size, random_state=self.random_seed)
            train_idx, test_idx = next(splitter.split(X, Y, groups=self.sample_groups))
            X_train_raw, X_test_raw = X[train_idx], X[test_idx]
            Y_train_raw, Y_test_raw = Y[train_idx], Y[test_idx]
            self.log(
                f"Grouped train/test split: train rows={len(train_idx)}, test rows={len(test_idx)}, "
                f"train Twiss groups={len(np.unique(self.sample_groups[train_idx]))}, "
                f"test Twiss groups={len(np.unique(self.sample_groups[test_idx]))}"
            )
        else:
            X_train_raw, X_test_raw, Y_train_raw, Y_test_raw = train_test_split(X, Y, test_size=self.test_size, random_state=self.random_seed, shuffle=True)

        self.x_scaler = StandardScaler() # scales data
        self.y_scaler = StandardScaler()
        X_train = self.x_scaler.fit_transform(X_train_raw)
        X_test = self.x_scaler.transform(X_test_raw)
        Y_train_log = np.log(np.maximum(Y_train_raw, 1e-30))
        Y_test_log = np.log(np.maximum(Y_test_raw, 1e-30))
        self.y_scaler = StandardScaler()
        Y_train = self.y_scaler.fit_transform(Y_train_log)
        Y_test = self.y_scaler.transform(Y_test_log)
        device = self.get_device()
        self.model = NeuralNetwork(X_train.shape[1], Y_train.shape[1]).float().to(device) # chooses gpu or cpu
                                                                                    # X_train.shape[1] = 7
                                                                                    # Y_train.shape[1] = 2 * n_screens
        optimizer = torch.optim.AdamW(self.model.parameters(), lr=self.learning_rate, weight_decay=self.weight_decay) # creates an algorithm, that fixes network weights, AdamW is a standard optimizer for
                                                                                                                    # neural networks, W stands for weight decay, that is optimized well
        loss_function = nn.MSELoss() # penalty for difference between prediction-true_sigma
        train_loader = self.make_loader(X_train, Y_train, batch_size=self.batch_size, shuffle=True)
        test_loader = self.make_loader(X_test, Y_test, batch_size=self.batch_size, shuffle=False)

        best_state = None
        best_test_loss = np.inf
        bad_epochs = 0

        for epoch in range(1, self.max_epochs+1):
            if self.should_stop():
                break
            self.model.train()
            train_losses = []

            for xb, yb in train_loader: # one data batch
                xb, yb = xb.to(device), yb.to(device)
                optimizer.zero_grad(set_to_none=True)
                prediction = self.model(xb) # predicts sigmas
                loss = loss_function(prediction, yb) # calculates difference
                loss.backward() # how to change weights to minimize error
                optimizer.step() # updates network's weights
                train_losses.append(float(loss.detach().cpu()))

            self.model.eval() # checks model on test data, without learning
            test_losses = []

            with torch.no_grad(): # model checks samples that were not seen during learning and we check how well it reconstructs Y
                for xb, yb in test_loader:
                    xb, yb = xb.to(device), yb.to(device)
                    prediction = self.model(xb)
                    loss = loss_function(prediction, yb)
                    test_losses.append(float(loss.detach().cpu()))

            test_loss = float(np.mean(test_losses)) if test_losses else np.inf
            if np.isfinite(test_loss) and test_loss < best_test_loss:
                best_test_loss = test_loss
                best_state = {k: v.detach().cpu().clone() for k,v in self.model.state_dict().items()}
                bad_epochs = 0
            else:
                bad_epochs += 1 # saves best model, if test loss got better, it saves current weights
            self.emit_progress(epoch, self.max_epochs)

            if epoch == 1 or epoch % 25 ==0: # prints result from 1st epoch and then every 25th, so doesn't create a spam
                train_loss = float(np.mean(train_losses)) if train_losses else np.nan
                self.log(f"Epoch {epoch}/{self.max_epochs}, train loss: {train_loss}, test loss: {test_loss}")
            if bad_epochs >= self.patience:
                self.log("Early stopping")
                break
        if best_state is not None:
            self.model.load_state_dict(best_state)

        Y_pred = self.predict_array(X_test_raw) # does prediction on unscaled test data
        r2_per_output = r2_score(Y_test_raw, Y_pred, multioutput="raw_values")
        self.metrics = {
            "mae": float(mean_absolute_error(Y_test_raw, Y_pred)),
            "rmse": float(np.sqrt(mean_squared_error(Y_test_raw, Y_pred))),
            "r2": float(r2_score(Y_test_raw, Y_pred, multioutput="variance_weighted")), # near 1 = perfect, it tells how well model is corresponding to data
            "r2_per_output": [float(value) for value in r2_per_output]
        }

        self.log(f"MAE: {self.metrics['mae']}, RMSE: {self.metrics['rmse']}, R2: {self.metrics['r2']}")
        if self.sigma_names and len(self.sigma_names) == len(r2_per_output):
            for name, value in zip(self.sigma_names, r2_per_output):
                self.log(f"R2[{name}] = {float(value)}")
        self.save_model()
        return self.metrics

    def predict_array(self, X): # for already trained model
        if self.model is None or self.x_scaler is None or self.y_scaler is None:
            raise RuntimeError("Model and scaler not initialized.")
        # model wants data in a format n_samples * n_params
        X = np.asarray(X, dtype=float)
        if X.ndim == 1:
            X = X.reshape(1, -1)

        model_inputs = int(self.model.net[0].in_features)
        model_outputs = int(self.model.net[-1].out_features)
        actual_inputs = int(X.shape[1])
        if actual_inputs != model_inputs:
            raise RuntimeError(
                f"ML input shape mismatch before PyTorch forward pass: X.shape={X.shape}, "
                f"model expects {model_inputs} inputs and returns {model_outputs} outputs. "
                f"Model file={self.model_file}."
            )
        if not np.all(np.isfinite(X)):
            bad_rows = np.where(~np.all(np.isfinite(X), axis=1))[0]
            raise RuntimeError(
                f"ML input contains NaN or Inf before prediction. Bad rows={bad_rows[:10].tolist()}, "
                f"X.shape={X.shape}, Model file={self.model_file}."
            )

        device = self.get_device()
        self.model.float().to(device)
        self.model.eval() # prediction mode, not training
        X_scaled = self.x_scaler.transform(X) # scales the same way, as during learning
        if not np.all(np.isfinite(X_scaled)):
            bad_rows = np.where(~np.all(np.isfinite(X_scaled), axis=1))[0]
            raise RuntimeError(
                f"Scaled ML input contains NaN or Inf. Bad rows={bad_rows[:10].tolist()}, "
                f"X_scaled.shape={X_scaled.shape}, Model file={self.model_file}."
            )

        try:
            with torch.no_grad(): # don't calculate gradients (what direction and how much should we change weights, so the error is smaller), because it's not learning
                X_tensor = torch.tensor(X_scaled, dtype=torch.float32, device=device)
                Y_scaled = self.model(X_tensor).detach().cpu().numpy()
        except Exception as e:
            raise RuntimeError(
                f"PyTorch forward pass failed inside ML surrogate. "
                f"X.shape={X.shape}, X_scaled.shape={X_scaled.shape}, X_tensor.shape={tuple(X_tensor.shape)}, "
                f"model_inputs={model_inputs}, model_outputs={model_outputs}, device={device}, "
                f"Model file={self.model_file}. Original error: {type(e).__name__}: {e}"
            ) from e

        if not np.all(np.isfinite(Y_scaled)):
            raise RuntimeError(
                f"ML model returned NaN or Inf. Y_scaled.shape={Y_scaled.shape}, Model file={self.model_file}."
            )
        Y_log = self.y_scaler.inverse_transform(Y_scaled) # unscales
        Y = np.exp(Y_log)

        return Y

    def save_model(self):
        if self.model is None or self.x_scaler is None or self.y_scaler is None:
            raise RuntimeError("Nothing to save. Train model first.")
        self.model_file.parent.mkdir(parents=True, exist_ok=True)
        torch.save({
            "state_dict": self.model.state_dict(), # trained network weights
            "n_inputs": int(self.model.net[0].in_features),
            "n_outputs": int(self.model.net[-1].out_features),
            "x_mean": self.x_scaler.mean_,
            "x_scale": self.x_scaler.scale_,
            "y_mean": self.y_scaler.mean_,
            "y_scale": self.y_scaler.scale_,
            "param_names": self.param_names,
            "sigma_names": self.sigma_names,
            "target_observable": self.target_observable,
            "screens": self.screens,
            "quad_name": self.quad_name,
            "machine_name": self.machine_name,
            "metrics": self.metrics,
        }, self.model_file)
        self.log(f"Saved ML model: {self.model_file}")

    def load_model(self):
        if not self.model_file.exists():
            raise FileNotFoundError(f"ML model file not found: {self.model_file}")
        checkpoint = torch.load(self.model_file, map_location="cpu", weights_only = False)
        self.model = NeuralNetwork(int(checkpoint["n_inputs"]), int(checkpoint["n_outputs"])).float() # takes the same neural network architecture
        self.model.load_state_dict(checkpoint["state_dict"]) # loads nn weights
        self.model.float()
        self.model.eval() # prediction mode
        self.x_scaler = StandardScaler()
        self.x_scaler.mean_ = np.asarray(checkpoint["x_mean"], dtype=float)
        self.x_scaler.scale_ = np.asarray(checkpoint["x_scale"], dtype=float)
        self.x_scaler.var_ = self.x_scaler.scale_ ** 2
        self.x_scaler.n_features_in_ = self.x_scaler.mean_.size

        self.y_scaler = StandardScaler()
        self.y_scaler.mean_ = np.asarray(checkpoint["y_mean"], dtype=float)
        self.y_scaler.scale_ = np.asarray(checkpoint["y_scale"], dtype=float)
        self.y_scaler.var_ = self.y_scaler.scale_ ** 2
        self.y_scaler.n_features_in_ = self.y_scaler.mean_.size

        self.param_names = [str(v) for v in checkpoint.get("param_names", [])]
        self.sigma_names = [str(v) for v in checkpoint.get("sigma_names", [])]
        self.target_observable = str(checkpoint.get("target_observable", "sigma2"))
        if self.target_observable not in {"sigma", "sigma2"}:
            raise RuntimeError(f"Unsupported ML model target: {self.target_observable}")
        self.screens = [str(v) for v in checkpoint.get("screens", self.screens)]
        self.quad_name = str(checkpoint.get("quad_name", self.quad_name))
        self.machine_name = str(checkpoint.get("machine_name", self.machine_name))
        self.metrics = dict(checkpoint.get("metrics", {}))
        return self # allows to do trainer = TrainModel(...).load_model(), it allows chaining

class MLInterface:
    def __init__(self, interface, quad_name, screens, machine_name):
        self.interface = interface
        self.quad_name = str(quad_name)
        self.screens = [str(s) for s in screens]
        self.machine_name = str(machine_name)
        self.model_file = get_ml_model_file(self.machine_name, self.quad_name, self.screens)
        self.trainer = TrainModel(screens=self.screens, quad_name=self.quad_name, machine_name=self.machine_name, model_file=self.model_file)
        self.trainer.load_model()
        if self.trainer.target_observable == "sigma2":
            print(f"Legacy ML model predicts sigma^2; retrain {self.model_file} to use direct sigma predictions.")

        self._screen_indices = [self.trainer.screens.index(screen) for screen in self.screens]

    def __getattr__(self, name):
        return getattr(self.interface, name) # if there is a function that MLInterface doesn't have (almost all of them), it gets them form the
                                            # interface, but has predict_emittance_scan_response, so uses that

    def predict_array(self, X):
        X = np.asarray(X, dtype=float)
        model_inputs = int(self.trainer.model.net[0].in_features)
        actual_inputs = int(X.size if X.ndim == 1 else X.shape[1])
        if actual_inputs != model_inputs:
            raise RuntimeError(
                f"MLInterface input shape mismatch: got {actual_inputs} input columns, "
                f"but model expects {model_inputs}. X.shape={X.shape}, Model file={self.model_file}."
            )
        return self.trainer.predict_array(X)

    def predict_emittance_scan_response(self, quad_name, screens, K1L_values, emit_x, emit_y, beta_x0, beta_y0, alpha_x0, alpha_y0, reference_screen = None, stop_checker = None):

        if callable(stop_checker) and stop_checker():
            raise RuntimeError("__OPTIMIZATION_STOP__")

        screens = [str(screen) for screen in screens]
        if str(quad_name) != self.quad_name:
            raise RuntimeError(f"ML model was loaded for a different quadrupole than optimizer requests.")
        if screens != self.screens:
            raise RuntimeError(f"ML model was loaded for different screens than optimizer requests.")

        X = []

        for K1L in K1L_values:
            X.append([emit_x, beta_x0, alpha_x0, emit_y, beta_y0, alpha_y0, K1L])
            # scan_points * 7

        X = np.asarray(X, dtype=float)
        Y = self.predict_array(X) # pytorch model
        n_trained_screens = len(self.trainer.screens)

        if Y.shape[1] != 2*n_trained_screens:
            raise RuntimeError(f"ML model output has wrong size.")

        # Y holds predictions for every trained screen; pick out only the requested ones, in order.
        prediction_sigx = Y[:, :n_trained_screens][:, self._screen_indices]
        prediction_sigy = Y[:, n_trained_screens:][:, self._screen_indices]

        if self.trainer.target_observable == "sigma2":
            prediction_sigx = np.sqrt(np.maximum(prediction_sigx, 0.0))
            prediction_sigy = np.sqrt(np.maximum(prediction_sigy, 0.0))
        else:
            prediction_sigx = np.maximum(prediction_sigx, 0.0)
            prediction_sigy = np.maximum(prediction_sigy, 0.0)

        return prediction_sigx, prediction_sigy
