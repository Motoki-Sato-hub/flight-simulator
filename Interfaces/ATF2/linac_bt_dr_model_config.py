"""Shared, explicitly offline configuration for Linac--BT--DR studies."""

# This voltage brings the current SAD-derived Linac+BT reference exit to the
# deterministic RF-Track synchronous momentum at the physical DR injection
# boundary KII.1.  It is a model-tuning value, not an ATF operating setting:
# RF phase, BT flight time and injection-kicker calibration remain external
# inputs before any real-machine interpretation.
MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV = 75.85950789
