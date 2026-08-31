# Installation
Requires: Python >= 3.11 and Poetry >= 2.0

### Install Poetry (if not yet installed)
```bash
curl -sSL https://install.python-poetry.org | python3 -
```
### Install project dependencies
```bash
poetry install --no-root
````
### Run

```bash
make BBA
make CRM
make SysID
```

### Alternatively:
```bash
poetry run python BBA_GUI.py
poetry run python ComputeResponseMatrix_GUI.py
poetry run python SysID_GUI.py
```

### While running Emittance_Measurement_GUI.py for the first time on a system:
#### Linux:
```bash
cd "Machine specifics, user implementations/ATF2/libmagnet"
make clean
make
gfortran -shared -fPIC -o libmagnet.so *.o
```
#### MacOS:
```bash
cd "Machine specifics, user implementations/ATF2/libmagnet"
make clean
make
gfortran -dynamiclib -o libmagnet.dylib *.o
```

# FlightSimulator's `Interface` description

This class defines how the user interacts with the real or simulated machine when reading beam instrumentation and steering devices and for writing corrector set-points in a control system for the specific accelerator facility. It wraps common operations on BPMs, ICTs, and corrector magnets.

Flight Simulator uses a common State / Interface architecture to make the same applications work across different accelerators and backends.

* `Interface` is the machine-dependent layer, it knows how to communicate with one specific machine or one specific simulation backend.
* `State` is a passive snapshot of the machine at one moment in time.
It stores device names, readings, settings, and other beamline-related information in a common structure.
---

## Constructor

### `__init__(nsamples: int = 1)`

Initializes the interface and builds internal lists of devices (sequence, BPMs, correctors) and EPICS-like indexes.

**Parameters**

* `nsamples` — number of consecutive samples to take when reading BPMs (with a 1 s pause between samples).

---

## Identification

### `get_name() -> str`

Returns a short identifier of this interface: e.g., `"ATF2_Ext"`.

---

## Lattice & device inventory

### `get_sequence() -> list[str]`

Returns the full **beamline-ordered** list of device names (BPMs and correctors).

### `get_hcorrectors_names() -> list[str]`

Returns the subset of correctors that are **horizontal** (prefixes need to be specified, e.g `ZH`).

### `get_vcorrectors_names() -> list[str]`

Returns the subset of correctors that are **vertical** (prefixes need to be specified, e.g `ZV`).

### `get_elements_indices(names: list[str] | str) -> list[int]`

Returns the **sequence indices** of the specified element name(s). Useful to map devices to their lattice order.

---

## Beam instrumentation reads

### `get_bpms() -> dict`

Reads BPM data from, e.g. `LINAC:monitors` for `nsamples` shots, one per second.
**Returns** a dictionary:

* `"names"` — list of BPM names
* `"x"` — `(nsamples × NBPM)` array of **x** positions in **mm**
* `"y"` — `(nsamples × NBPM)` array of **y** positions in **mm**
* `"tmit"` — `(nsamples × NBPM)` array of intensity/charge (TMIT), **zeroed where BPM status ≠ 1**

> Note: Internally selects only those rows corresponding to `self.bpms` using `self.bpm_indexes`.

### `get_icts() -> dict`

Reads charge from all configured ICT PVs.
**Returns** a dictionary:

* `"names"` — list of ICT PV names
* `"charge"` — 1-D NumPy array of the corresponding readings

### `get_screens() -> dict`

Reads screen measurements for all devices in `self.screens`.
**Returns** a dictionary:

* `"names"` — list of screen names
* `"hpixel"` — array of horizontal pixel sizes
* `"vpixel"` — array of vertical pixel sizes
* `"x"` — mean position of particles in x plane
* `"y"` — mean position of particles in y plane
* `"sigx"` — RMS beam size in x plane 
* `"sigy"` — RMS beam size in y plane 
* `"sum"` — number of particles on the screen; corresponds to intensity
* `"hedges"` — bin edges in x (nx + 1)
* `"vedges"` — bin edges in y (ny + 1)
* `"images"` — image[ i , j ] = nparticles in bin i on x axis and nparticles in bin j on y axis
* `"S"` — position of screen along the lattice
---

## Magnet reads & writes

### `get_correctors() -> dict`

Reads corrector strengths for all devices in `self.corrs`.
**Returns** a dictionary:

* `"names"` — list of corrector names
* `"bdes"` — desired/setpoint currents from `:<name>:currentWrite`
* `"bact"` — readback/actual currents from `:<name>:currentRead`

### `set_correcors(names: list[str] | str, corr_vals: list[float] | np.ndarray | float) -> None`

**Sets** the desired current for one or more correctors. If a single name/value is provided, it’s coerced to a list. Lengths must match. Writes to `:<name>:currentWrite` and sleeps 1 s.

### `vary_correctors(names: list[str] | str, corr_vals: list[float] | np.ndarray | float) -> None`

**Increments** the desired current for one or more correctors by the given value(s). Reads current value from `:<name>:currentWrite`, adds the increment, writes back, then sleeps 1 s.

### `get_quadrupoles() -> dict`

Reads quadrupole strengths for all devices in `self.quadrupoles`.
**Returns** a dictionary:

* `"names"` — list of quadrupole names
* `"bdes"` — desired/setpoint currents from `:<name>:currentWrite`
* `"bact"` — readback/actual currents from `:<name>:currentRead`

### `set_quadrupoles(names: list[str] | str, corr_vals: list[float] | np.ndarray | float) -> None`

**Sets** the desired current for one or more quadrupoles. If a single name/value is provided, it’s coerced to a list. Lengths must match. Writes to `:<name>:currentWrite` and sleeps 1 s.

### `vary_quadrupoles(names: list[str] | str, corr_vals: list[float] | np.ndarray | float) -> None`

**Increments** the desired strength for one or more quadrupoles by the given value(s). Reads current value from `:<name>:currentWrite`, adds the increment, writes back, then sleeps 1 s.

---

## DFS and WFS manipulations
For simulation-based interfaces, these may trigger retracking.

For online interfaces, they may correspond to control-system operations, like changing the klystron
phase or laser angle.

### `change_energy() -> None`


### `reset_energy() -> None`


### `change_intensity() -> None`


### `reset_intensity() -> None`


---

## State getter and other helpers

### `get_target_dispersion() -> target_disp_x: list[float], target_disp_y: list[float]`

Returns values of target dispersion in two planes.

### `get_state() -> State`

Builds and returns a State object containing the current machine snapshot.

### `restore_correctors_state() -> None`

Restores corrector settings from a previously saved State.

### `restore_quadrupoles_state() -> None`

Restores quadrupole settings from a previously saved State.

---

# How to add a new accelerator to the system

Flight Simulator uses a common State / Interface architecture and an interface_setup.py registry.
To add a new accelerator, follow the steps below.

### Create a new `Interface` class

Add a new Python file in Interfaces/`Accelerator name`/, for example:
```
Interfaces/AcceleratorName/Interface_AcceleratorName_DR.py
```

This class should inherit from AbstractMachineInterface and implement all needed methods, required by Flight Simulator tools.
Example:
```
class Interface_AcceleratorName_DR(AbstractMachineInterface)
```
### Implement the required methods:
Interface has to provide the machine state in the common format used by State.
At minimum, implement the methods needed in Flight Simulator tools:

At minimum, implement the methods needed by the applications:
*	get_name()
*	get_sequence()
*	get_correctors()
*	get_bpms()
*	get_icts()
*	get_quadrupoles()
*	get_screens()
*	get_hcorrectors_names()
*	get_vcorrectors_names()
*	set_correctors()
*	vary_correctors() if used by correction tools
*	set_quadrupoles() if quadrupole scans are used
*	change_energy() / reset_energy() if used
*	change_intensity() / reset_intensity() if used
*	get_state()
*	restore_correctors_state()
*	restore_quadrupoles_state()

Depending on the machine, some of these methods may be simple wrappers, but the external API should remain consistent.

The returned data must match the dictionaries expected by State, e.g.:
*	`correctors`: names, bdes, bact
*   `bpms`: names, x, y, tmit
*   `screens`: names, sigx, sigy, etc.
*   `quadrupoles`: names, bdes, bact

### Make sure State can be built:

State is the common snapshot container used by the GUIs and algorithms.
Your interface should gather the machine data and return a State object populated with:
*	correctors
*	BPMs
*	ICTs
*	screens
*	quadrupoles
*	sequence
*	horizontal / vertical corrector names

This is what allows the same applications to work across different accelerators.

### Add the new accelerator to interface_setup.py:
Register the new interface in interface_setup.py, located in `Interfaces` folder.

Example:
```
INTERFACE_SETUP = {
    "NewAccelerator": [
        {
            "display_name": "New Accelerator",
            "module": "Interfaces.NewAcc.InterfaceNewAcc",
            "class_name": "InterfaceNewAcc",
            "settings": {"nsamples": 3},
            "actions": [],
            "units": {
                "corrector_strength": "T*mm",
                "bpm_position": "mm",
                "sysid_corrector_kick": 0.01
            }
        }
    ]
}
```
In interface_setup.py, define the machine-specific data needed by the applications:

`corrector units`: e.g. T * mm

`BPM units`: e.g. mm

`default SysID kick`: e.g. 0.01

Remember, that `module` needs to be a Python-like path to the new Interface, `class_name` strictly a class specified in the interface and 
`display_name`, a desired name shown in the window used for initial selection. In `settings`, a user can provide `nsamples`
for the real machine interface, and additionally, `jitter` and `bpm_resolution` in the simulated one. 

If the machine requires special initialization steps, these can be added in the actions list in interface_setup.py and
will be executed automatically after the interface is created. For example, a simulated interface can have misaligned 
components:
```
"actions": ["align_everything","misalign_quadrupoles","misalign_bpms"],
```

---
## Notes & best practices

* **Sampling**: Increase `nsamples` in the constructor if you want automatic multi-shot acquisition with per-shot 1 s spacing.
* **Units**: BPM positions are returned in **mm**; intensities are passed through the BPM status gate.
* **Safety**: `set_correctors` and `vary_correctors` talk directly to setpoint PVs. Use with care in live machines and follow your facility’s MPS/permit procedures.
* **Extensibility**: The placeholder energy/intensity methods can be filled with site-specific PVs or finite-state-machine hooks if needed.
---


