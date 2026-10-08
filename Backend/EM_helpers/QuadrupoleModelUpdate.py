import os
from datetime import datetime, timezone
import numpy as np

def quadrupole_tracking_model(interface):
    model = getattr(interface, "tracking_interface", None)
    if model is None and getattr(interface, "is_simulation", False):
        model = interface
    if model is None or not getattr(model, "is_simulation", False):
        raise ValueError("A simulation interface is required to load model quadrupoles.")
    if model.get_name() != "CLEAR_RFT":
        raise ValueError("Loading quadrupole current snapshots is currently supported for CLEAR RF-Track.")
    return model


def validate_quadrupole_status(payload, model):
    unit = np.asarray(payload.get("value_unit", ""))
    if unit.size != 1 or str(unit.item()) != "A":
        raise ValueError("The quadrupole snapshot must declare value_unit='A'.")
    names_array = np.asarray(payload.get("names", []))
    if names_array.ndim != 1 or not names_array.size:
        raise ValueError("The snapshot must contain a nonempty one-dimensional names array.")
    names = [str(name) for name in names_array.tolist()]
    if len(set(names)) != len(names):
        raise ValueError("The quadrupole snapshot contains duplicate names.")
    unknown = sorted(set(names) - set(model.quadrupoles))
    if unknown:
        raise ValueError(f"Quadrupoles not present in the RF-Track model: {', '.join(unknown)}")
    readback_fields = [field for field in ("iact", "bact") if field in payload]
    if not readback_fields:
        raise ValueError("The snapshot needs measured currents in iact or bact.")
    for field in readback_fields:
        values = np.asarray(payload[field], dtype=float)
        if values.ndim != 1 or values.size != len(names) or not np.all(np.isfinite(values)):
            raise ValueError(f"{field} must contain one finite current per quadrupole.")
    field = readback_fields[0]
    return names, np.asarray(payload[field], dtype=float).copy(), field


def quadrupole_status_metadata(payload, source_path):
    captured_at = np.asarray(payload.get("captured_at_utc", ""))
    captured_at = str(captured_at.item()) if captured_at.size == 1 else ""
    momentum = np.asarray(payload.get("reference_momentum_MeV_c", []), dtype=float)
    momentum = float(momentum.item()) if momentum.size == 1 and np.all(np.isfinite(momentum)) else None
    return {
        "source_path": os.path.abspath(os.fspath(source_path)),
        "captured_at_utc": captured_at or None,
        "legacy_without_timestamp": not bool(captured_at),
        "reference_momentum_MeV_c": momentum,
        "value_unit": str(np.asarray(payload.get("value_unit", "")).item()),
        "quadrupole_count": int(np.asarray(payload.get("names", [])).size),
    }


def load_model_quadrupoles(interface, source_path):
    model = quadrupole_tracking_model(interface)
    with np.load(source_path, allow_pickle=False) as saved:
        payload = {key: saved[key].copy() for key in saved.files}
    names, currents, field = validate_quadrupole_status(payload, model)
    metadata = quadrupole_status_metadata(payload, source_path)
    metadata.update({
        "applied_at_utc": datetime.now(timezone.utc).isoformat(),
        "readback_field": field,
        "model_name": model.get_name(),
        "model_reference_momentum_MeV_c": float(model.Pref),
    })
    original_lattice = model.lattice
    model.lattice = original_lattice.clone()
    try:
        model.set_quadrupoles(names, currents, track=False)
    except Exception:
        model.lattice = original_lattice
        raise
    return payload, metadata


def save_quadrupole_status(payload, target_path):
    target_path = os.path.abspath(os.fspath(target_path))
    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    if os.path.exists(target_path):
        try:
            with np.load(target_path, allow_pickle=False) as saved:
                identical = set(saved.files) == set(payload)
                if identical:
                    for key, value in payload.items():
                        current = saved[key]
                        expected = np.asarray(value)
                        numeric = np.issubdtype(current.dtype, np.number) and np.issubdtype(expected.dtype, np.number)
                        if not np.array_equal(current, expected, equal_nan=numeric):
                            identical = False
                            break
                if identical:
                    return target_path
        except Exception:
            pass
        stem, extension = os.path.splitext(target_path)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        target_path = f"{stem}_{stamp}{extension}"
        index = 1
        while os.path.exists(target_path):
            target_path = f"{stem}_{stamp}_{index}{extension}"
            index += 1
    np.savez(target_path, **payload)
    return target_path
