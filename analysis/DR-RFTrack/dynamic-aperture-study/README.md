# Offline ATF DR acceptance study — 2026-10-07

RF-Track 2.6.3; current SAD-derived 2011 DR model, entry KII.1,
deterministic RF/radiation, zero pulsed kickers. These are finite-turn damped
capture results, not measured ATF dynamic aperture or final transmission.
The literature +/-6 mm value is recorded for context but never used as a cut.

## 100-turn comparison

Directories `no-physical-screen` and `historical-kix-screen` contain JSON
and PNG outputs. Each scan used radii 0, 0.25, 0.5, 1, 2, 3, 6 mm,
signed horizontal/vertical rays, all four combinations of x/y betatron
phases 0 and 90 degrees, and nominal momentum. The amplitude parameters
are local Courant--Snyder amplitudes, not simply instantaneous x and y.
There are 112 samples, including duplicate zero-amplitude cases; the grid
is not an incoming beam distribution and does not sample diagonal rays.

- Horizontal sampled first-loss brackets: 1 mm survives, 2 mm fails.
- Vertical sampled first-loss brackets: 3 mm survives, 6 mm fails.
- Each screen condition: 80 surviving samples, 30 RF-Track losses,
  and 2 numerical-guard rejections. These are not bunch-transmission counts.
- The separately tracked synthetic 16-particle Linac bunch had Linac+BT
  charge transmission 1.0 and conditional 100-turn DR survival 1.0.
- The two screen conditions produced identical sampled survival decisions.
  This does not imply physical apertures are irrelevant or fully represented.
- Runtime was about 84–85 seconds per scan while both jobs overlapped;
  this is not a controlled single-process performance benchmark.

The historical KIX screen uses 5-mm half aperture and an assumed circular
cross-section. It does not locate arc, wiggler or south-straight masks.
Numerical guards (100 mm / 100 mrad at the turn boundary) are not chamber walls.

## Momentum slices: 10 turns

`momentum-slices` uses radii 0, 1, 3 mm, positive horizontal/vertical rays,
four betatron-phase combinations and dp/p = 0, -0.005, +0.005. Surviving
sample counts were respectively 20/24, 16/24 and 20/24. The directly supplied
four-particle synthetic entrance fixture also survived the full pipeline.
Neither sample-count ratio is a measured transmission or an integrated
momentum acceptance. The horizon differs from the 100-turn comparison.

## Reproduction

From the flight-simulator repository root, with RF-Track available:

```bash
MPLCONFIGDIR=/tmp/mpl-rftrack PYTHONPATH=. \
  /home/motokisato/rftrack-env/bin/python \
  -m Interfaces.ATF2.scan_atf_dr_dynamic_aperture \
  --turns 100 --radii-mm 0,0.25,0.5,1,2,3,6 \
  --rays-deg 0,90,180,270 --phases-deg 0,90 \
  --pipeline-particles 16 --output-dir /tmp/atf-dr-acceptance
```

Repeat in another directory with `--aperture-screen historical-kix`.
For the momentum test use `--turns 10 --radii-mm 0,1,3 --rays-deg 0,90
--delta-p=0,-0.005,0.005` and replace `--pipeline-particles 16` with
`--entrance-bunch-json Interfaces/ATF2/linac_bt_dr_entrance_bunch_synthetic_example.json`.

Further validation needs denser diagonal/phase/momentum sampling and longer
horizons, a calibrated injection map/pulse, complete physical apertures,
field-error/nonlinearity checks and a realistic entrance distribution.
