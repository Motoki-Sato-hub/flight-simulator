# ATF DR RF-Track lattice

`ATF_DR_RFTrack_lattice.py` builds the RF-Track damping-ring lattice from the
checked-in `ATF_DR_RFTrack_lattice.json` data.

The reference SAD daihon is
**`operation/daihon/atfdr-design-20111111b.sad`** (relative to the `sad`
directory).  The JSON stores its SHA-256, element definitions, and the order in
`LINE RING0`.  It was made without executing SAD.  To regenerate it after an
intentional daihon update:

```bash
python generate_atf_dr_rftrack_lattice.py \
  /path/to/sad/operation/daihon/atfdr-design-20111111b.sad
```

The current RF-Track model is for transverse correction studies.  The SAD
zero-length RF cavity is represented by a zero-length drift, and synchrotron
radiation is not enabled.  This keeps particle momentum constant.  Do not use
this version for longitudinal capture, radiation damping, or equilibrium
emittance studies.

The generated ring has a circumference of 138.55953941176458 m, 36 bends,
100 quadrupole instances, 98 monitor instances, 50 horizontal steerers and 51
vertical steerers.  `InterfaceATF2_DR_RFTrack` uses this builder instead of
inserting artificial steerers into every third quadrupole.

`ATF_DR_twiss_file.tws` remains the legacy MAD-X target-optics file and was
not made from this exact daihon revision.  The checked-in
`ATF_DR_20111111b_RFTrack_twiss.tfs` is instead generated from this RF-Track
lattice, linearized about its nominal periodic closed orbit.  Its filename
identifies the reference daihon revision.  It contains one row for every
RF-Track lattice element (including BPMs, magnets, drifts and the added thin
skew multipoles), with beta, alpha, closed orbit and dispersion.  It is the
nominal optics reference used by the correction-validation notebook.  The
table also contains the magnetic/import fields needed by `RF_Track.Lattice`,
so the notebook reads this TFS directly as its lattice instead of rebuilding
the lattice in memory.  Regenerate it with:

```bash
PYTHONPATH=. /path/to/rftrack-env/bin/python \
  Interfaces/ATF2/DR_ATF2/ATF_DR_RFTrack_twiss.py
```

A final SAD-to-RF-Track optics equivalence claim still requires reference
Twiss/COD/dispersion output from the exact SAD daihon.

## Periodic orbit correction

`ATF_DR_RFTrack_correction.py` implements the model-side equivalent of the SAD
workflow in `operation/lib/cod.n`:

1. solve the periodic fixed point `one_turn(z) = z`;
2. read the COD at the RF-Track BPMs;
3. perturb the real ZH or ZV steerers and build the periodic response matrix;
4. calculate an SVD least-squares correction;
5. optionally use SAD-like forward steerer selection with `max_correctors`;
6. keep calculation and application as separate operations.

`InterfaceATF2_DR_RFTrack` exposes this through
`get_closed_orbit`, `get_model_dispersion`, `get_nominal_model_dispersion`,
`compute_periodic_orbit_response`, `suggest_periodic_orbit_correction`, and
`apply_periodic_orbit_correction`.  For example:

```python
response = interface.compute_periodic_orbit_response("x")
suggestion = interface.suggest_periodic_orbit_correction(
    response,
    max_correctors=20,
    rcond=1e-4,
    gain=0.8,
)
print(suggestion.rms_before, suggestion.rms_predicted)
interface.apply_periodic_orbit_correction(suggestion)
```

The SAD `EtaCorrect` combination of dispersion and orbit constraints is
available through `compute_periodic_dispersion_response`,
`suggest_periodic_dispersion_correction`, and
`apply_periodic_dispersion_correction`.  The nominal target is captured from
the unperturbed lattice when the interface starts.  The dispersion and COD
weights are independent; their defaults match the 0.05 dispersion factor used
in `operation/lib/coddispersion.n`.

The interface response is in mm per interface corrector unit, consistent with
`get_correctors` and `set_correctors`.  The lower-level correction class can
instead use raw RF-Track corrector-strength units.

## Skew-coupling correction

The SAD procedure in `operation/lib/skewcor.n` uses two horizontal-steerer
changes and measures the induced vertical COD at the BPMs.  It then adjusts
the skew-quadrupole components at `SD1R.*` and `SF1R.*`.  The RF-Track builder
now adds an initially-zero thin skew K1L multipole directly after each of the
34 SD1R and 34 SF1R sextupoles; the normal sextupole and all daihon lengths
are unchanged.  The actuator names are, for example, `SF1R.1$SKEW` and
`SD1R.1$SKEW`.

`compute_periodic_coupling_response`,
`suggest_periodic_coupling_correction`, and
`apply_periodic_coupling_correction` provide the same response/SVD/apply
separation as the COD corrections.  The default SVD cutoff is the SAD value
of 0.4.  SAD's SD/SF current limits require calibration values from operation
data; because those values are not available here, use `max_abs_delta` to
provide a deliberate normalized K1L limit.

```python
response = interface.compute_periodic_coupling_response(
    probe_corrector_names=("ZH1R", "ZH2R"),
)
suggestion = interface.suggest_periodic_coupling_correction(
    response,
    max_skew_correctors=20,
    max_abs_delta=1e-3,
)
interface.apply_periodic_coupling_correction(suggestion)
```

`get_skew_correctors` and `set_skew_correctors` expose those normalized
integrated K1L values for simulation-error setup and inspection.

The DR DFS energy step is 0.5%.  This is close to the sub-percent momentum
change used by the SAD frequency-shift dispersion procedure and stays inside
the momentum range validated for this transverse-only model.
