# ATF Linac → BT → DR model-capture pipeline

This offline pipeline tracks a six-dimensional bunch from the SAD Linac
entrance marker `IPP1L`, through the accelerating Linac and BT, then applies
an explicit IPZT-to-`KII.1` handoff and tracks the DR turn by turn.  It never
connects to ATF controls.

```
IPP1L  →  SAD-derived accelerating Linac  →  BT / IPZT
       →  calibrated-or-design handoff  →  DR KII.1  →  N turns
```

`N`-turn survival is called **model capture** or **ring survival**.  It is the
end-to-end quantity that can be optimised now.  It is not physical final ATF
transmission until a surveyed handoff, pulsed-injection settings and a complete
aperture/loss map are supplied.  The current periodic DR model has no
extraction line, so it does not predict extracted-beam transmission.

The Linac+BT input preserves the SAD daihon element inventory, geometry and
design reference orbit, and re-scales magnets for the configured local energy.
Its standard RF-Track TFS representation does **not** yet reproduce the
historical BT SAD MARK Twiss at every downstream marker (notably SAD `COORD`
and BEND fringe maps).  Consequently, use the pipeline for staged model-capture
studies, not as a standalone precision BT beam-size/acceptance prediction.

The SAD `RING0` sequence contains two `KII` injection kickers and two `KIX`
extraction kickers.  Their daihon entries define only zero-angle geometry, not
an operating pulse.  The model now retains them as zero-strength RF-Track
correctors.  A study may give a *calibrated physical* first-turn kick in rad:

```python
machine = ATF2LinacBTDRRFTrack(
    dr_first_turn_kicks_rad={"KII": (theta_x_rad, theta_y_rad)}
)
```

The periodic lattice begins immediately upstream of `KII.1`, so the first
tracked passage applies the injection pulse before any downstream aperture.
`"KII"` applies the same pair to `KII.1` and `KII.2`; an instance name such as
`"KII.1"` selects one magnet.  The kick is cleared immediately after turn one,
including when tracking fails.  This supplies the time-gating mechanism; it
does not infer waveform, delay, amplitude calibration, or the IPZT-to-`IIN`
transport from the daihon.

## Reproducible software gate

From the repository root:

```bash
MPLCONFIGDIR=/tmp/mpl-rftrack PYTHONPATH=. \
  /home/motokisato/rftrack-env/bin/python \
  Interfaces/ATF2/validate_linac_bt_dr_end_to_end.py \
  --long-turns 100 --long-sample-every 20
```

This checks a deterministic finite bunch through Linac→BT→DR for 10 turns and
a reference particle for the selected long horizon.  It is an integration
gate, not a measured-transmission prediction.

To verify the direct measured/reconstructed bunch interface independently:

```bash
MPLCONFIGDIR=/tmp/mpl-rftrack PYTHONPATH=. \
  /home/motokisato/rftrack-env/bin/python \
  Interfaces/ATF2/validate_linac_bt_dr_direct_6d_input.py
```

## Inputs that replace study defaults

`benchmark_linac_bt_dr_capture_proxy.py` accepts a direct six-dimensional
input bunch, a fitted handoff and an explicit aperture screen:

```bash
MPLCONFIGDIR=/tmp/mpl-rftrack PYTHONPATH=. \
  /home/motokisato/rftrack-env/bin/python \
  Interfaces/ATF2/benchmark_linac_bt_dr_capture_proxy.py \
  --entrance-bunch-json /path/to/ipp1l_bunch.json \
  --handoff-json /path/to/ipzt_to_kii1_handoff.json \
  --dr-first-turn-kicks-json /path/to/kii_pulse_rad.json \
  --short-turns 10 --long-turns 1000 \
  --long-turn-history-sample-every 20
```

The direct bunch JSON must use exactly this coordinate order and boundary:

```json
{
  "coordinate_order": "x_mm,xp_mrad,y_mm,yp_mrad,t_mm_c,p_mev_c",
  "coordinates_mm_mrad_mm_c_mev_c": [[0, 0, 0, 0, 0, 80]],
  "charge_e": 1000000000.0,
  "provenance": "diagnostic and reconstruction description",
  "location": "IPP1L"
}
```

The handoff JSON is the `TransverseHandoff.as_dict()` payload produced by
`linac_bt_dr_handoff_fit.py`, or the nested `handoff` object in that tool's
result.  It maps `[x_mm, xp_mrad, y_mm, yp_mrad]` at the BT endpoint into the
DR `KII.1` upstream coordinates.  The JSON must contain
`"target_location": "KII.1"`; an old or differently referenced handoff is
rejected rather than silently used.  Supplying a bunch alone does not
calibrate this map.

The optional kicker JSON has physical kick radians, for example
`{"KII": [0.001, 0.0]}`. Omit it unless its amplitude and timing have been
calibrated independently; this model applies it to the first DR turn only.

An aperture table, when needed by a custom caller, has this form:

```json
{
  "linac_bt": {"IPZT": [10.0, 10.0, "rectangular"]},
  "dr": {"KIX.1": [5.0, 5.0, "circular"]}
}
```

Only entries explicitly listed are active.  Do not infer missing chambers
from magnet names or design apertures.

## Fast screen versus multi-turn evaluation

Use the same candidate bunch and machine settings in three stages:

1. `dr_turns=0`: rank many candidates by Linac+BT survival, energy error,
   DR-start orbit and dispersion-subtracted projected mismatch.
2. `dr_turns=5` or `10`: reject prompt dynamic-capture losses.
3. Run the selected candidates for the declared capture horizon with
   `record_turn_history=True`.  The reported `dr_survival` at the requested
   final turn is the model-capture objective.

The benchmark prints separate construction, endpoint and per-turn timings.
On the current offline host, a four-particle KII.1-start study took 0.078 s
for the endpoint and 0.25 s per DR turn (30 turns: 7.62 s including the
endpoint).  Treat those as capacity-planning figures, not portable
performance guarantees.

The endpoint values are deliberately not a replacement for the last stage:
the handoff sweep (`sweep_atf2_dr_injection_capture_proxy.py`) demonstrates
that a beam can traverse Linac+BT while being lost on a later DR turn.
`--turn-history-sample-every 1` records the exact first loss turn; a larger
value bounds report size but only brackets it.

The required capture horizon is an experiment-level choice.  A few turns test
prompt injection/capture; it cannot establish radiation damping or stored-beam
lifetime.  The DR equilibrium studies estimate a roughly 70,000-turn 1/e
damping time, so a physical storage/extraction study needs a separately
declared horizon and its corresponding RF, radiation, aperture and extraction
models.

## Corrector optimisation

`correct_linac_bt_dr_injection_orbit.py` is the local response-matrix
baseline.  `benchmark_linac_bt_dr_response_matrix_vs_rl.py` compares it with
a deliberately small model-free CEM baseline under equal synthetic errors and
actuator bounds.  Both optimise the cheap endpoint observation and validate
the selected setting with DR turns.  A method should not be called better than
the response matrix unless it also improves held-out multi-turn model capture
under identical actuator limits and, later, the calibrated physical model.

## Finite-turn dynamic acceptance and transmission

`scan_atf_dr_dynamic_aperture.py` scans signed local transverse amplitudes,
independent x/y betatron phases, and relative momentum offsets at `KII.1`.
Each probe starts from the same model state and is tracked independently.
It records the first rejection turn, all sampled survival islands, and a
bracket for the first radial loss. If all sampled radii survive, the reported
boundary is **beyond the scan**, not equal to its largest sampled radius.
The local mm amplitude parameter is not an invariant ring acceptance.

```bash
MPLCONFIGDIR=/tmp/mpl-rftrack PYTHONPATH=. \
  /home/motokisato/rftrack-env/bin/python \
  -m Interfaces.ATF2.scan_atf_dr_dynamic_aperture \
  --turns 100 --radii-mm 0,0.25,0.5,1,2,3,6 \
  --rays-deg 0,90,180,270 --phases-deg 0,90 \
  --delta-p=0,-0.005,0.005 --pipeline-particles 16 \
  --output-dir /tmp/atf-dr-acceptance
```

Outputs are `dr_acceptance.json` and `dr_acceptance.png`. The map's fraction
of surviving scan points is **not bunch transmission**: this grid is not a
sample of the incoming beam distribution. With `--pipeline-particles N`, a
separate synthetic equal-weight Linac bunch is transported through BT and the
existing handoff; its particles are directly tracked in the DR, without
interpolating or clipping against the coarse map. `--entrance-bunch-json`
can replace that synthetic input. The output separates conditional DR
survival from charge transmission through Linac+BT and their product.
The default assumes zero injection pulses. A periodic DA map does not replace
a calibrated first-turn injection trajectory/kicker calculation.

By default, RF and deterministic radiation are enabled: label the result
**damped RF capture / N-turn acceptance**, not conservative Hamiltonian DA.
`--rf-mode disabled` gives an RF/radiation-off transverse comparison, not
longitudinal capture. A conservative RF-on, radiation-off DA calculation is
not supplied by this command. Quantum excitation and scattering losses are
not included. Numerical guards of 100 mm and 100 mrad, evaluated only at turn
boundaries, are reported separately from RF-Track loss and are not physical
walls. Reaching the guard is not a measurement of chamber aperture.

Repeat with `--aperture-screen historical-kix` to compare the existing partial
5-mm KIX physical screen with the no-wall model. Its circular shape is an
assumption; this does not locate all masks or supply a complete chamber model.
The historical theoretical DA of +/-6 mm in [SLAC-R-559, printed p. 92,
April 2000](https://www.slac.stanford.edu/pubs/slacreports/reports08/slac-r-559.pdf)
is a literature comparison only: it is never imposed as a tracking wall or
assumed to be a limit at the current `KII.1` observation point.

The input's mean arrival time is aligned once to the model synchronous phase
by the existing handoff. `--time-offset-mm-c` then adds an explicit common
timing offset to both the scan and evaluated bunch, preserving bunch length
and time--momentum correlation. No per-particle or per-turn rephasing occurs.

Bookkeeping regression tests:

```bash
PYTHONPATH=. /home/motokisato/rftrack-env/bin/python \
  -m Interfaces.ATF2.validate_atf_dr_dynamic_aperture
```
