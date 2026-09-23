from __future__ import annotations
import sys, os , argparse, json, math
from pathlib import Path
from pprint import pprint
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))
from Backend.State import State
filename = Path("/local/home/clearop/CERN-Flight_Simulator-Data/CLEAR_20260902_184342_Orbit/machine_status.pkl")

state = State(filename=filename)
quadrupoles = state.get_quadrupoles()
pprint(quadrupoles)


# p_files = sorted(folder.glob("DATA_*_p0001.pkl"))
#
# for p_file in p_files:
#     corrector_name = (p_file.name.removeprefix("DATA_").removesuffix("_p0001.pkl"))
#     print(f"\n{'=' * 60}\nCorrector: {corrector_name}")
#     for pm in ("p", "m"):
#         filename = folder / f"DATA_{corrector_name}_{pm}0001.pkl"
#         state = State(filename=str(filename))
#         corrector = state.get_correctors([corrector_name])
#
#         print(f"{pm}:")
#         print("  file:", filename.name)
#         print("  bdes:", corrector["bdes"][0])
#         print("  bact:", corrector["bact"][0])
# print("===============================================")

# p_files = sorted(folder.glob("DATA_*_p0001.pkl"))
# for p_file in p_files:
#     corrector_name = (p_file.name.removeprefix("DATA_").removesuffix("_p0001.pkl"))
#     print(f"\n{'=' * 60}\nCorrector: {corrector_name}")
#     for pm in ("p", "m"):
#         filename = folder / f"DATA_{corrector_name}_{pm}0001.pkl"
#         state = State(filename=str(filename))
#         corrector = state.get_correctors([corrector_name])
#
#         print(f"{pm}:")
#         print("  file:", filename.name)
#         print("  bdes:", corrector["bdes"][0])
#         print("  bact:", corrector["bact"][0])
# print("===============================================")
#
# import argparse
# import json
# import math
# import sys
# from collections import Counter
# from dataclasses import dataclass
# from pathlib import Path
#
# # Allow this helper to run directly from its subdirectory.
# REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
# if str(REPOSITORY_ROOT) not in sys.path:
#     sys.path.insert(0, str(REPOSITORY_ROOT))
#
# from Backend.State import State
#
#
# DEFAULT_ROOT = Path("/local/home/clearop/CERN-Flight_Simulator-Data")
#
#
# @dataclass
# class DirectoryPhaseSummary:
#     """RF-phase information read from the saved states in one SysID directory."""
#
#     directory: Path
#     data_files: int
#     phase_counts: Counter
#     unreadable_files: list[tuple[Path, str]]
#
#     @property
#     def phases(self) -> tuple[float, ...]:
#         return tuple(sorted(self.phase_counts))
#
#
# def _read_phase(state_file: Path) -> tuple[float | None, object | None]:
#     """Return the saved CLEAR RF phase and state timestamp, if present."""
#
#     state = State(filename=str(state_file))
#     phase = state.get_beam_settings().get("energy", {}).get("mks11_phase")
#     try:
#         phase = float(phase)
#     except (TypeError, ValueError):
#         return None, state.get_timestamp()
#
#     return (phase if math.isfinite(phase) else None), state.get_timestamp()
#
#
# def summarise_data_directory(directory: Path, show_files: bool = False) -> DirectoryPhaseSummary:
#     """Read all saved SysID states from *directory* and report their RF phases."""
#
#     files = sorted(directory.glob("DATA_*.pkl"))
#     phases: Counter = Counter()
#     unreadable: list[tuple[Path, str]] = []
#
#     print(f"\n{directory}")
#     if not files:
#         print("  No DATA_*.pkl files found.")
#         return DirectoryPhaseSummary(directory, 0, phases, unreadable)
#
#     for state_file in files:
#         try:
#             phase, timestamp = _read_phase(state_file)
#         except Exception as error:  # A damaged historical file should not stop the full audit.
#             unreadable.append((state_file, str(error)))
#             if show_files:
#                 print(f"  {state_file.name}: unreadable ({error})")
#             continue
#
#         if phase is None:
#             if show_files:
#                 print(f"  {state_file.name}: phase unavailable; timestamp={timestamp}")
#             continue
#
#         phases[phase] += 1
#         if show_files:
#             print(f"  {state_file.name}: phase={phase:g} degrees; timestamp={timestamp}")
#
#     phase_text = ", ".join(
#         f"{phase:g} degrees ({count} files)" for phase, count in sorted(phases.items())
#     )
#     print(f"  Phase: {phase_text or 'unavailable'}")
#     if len(phases) > 1:
#         print("  WARNING: this directory contains more than one saved RF phase.")
#     if unreadable:
#         print(f"  WARNING: could not read {len(unreadable)} file(s).")
#
#     return DirectoryPhaseSummary(directory, len(files), phases, unreadable)
#
#
# def _directories_with_sysid_data(root: Path) -> list[Path]:
#     return sorted({path.parent for path in root.rglob("DATA_*.pkl")})
#
#
# def _resolve_data_directory(path: Path, scan_root: Path) -> Path | None:
#     """Resolve the path saved in BBA settings, including archived AllCLEAR data."""
#
#     saved_path = path.expanduser().resolve()
#     candidates = [saved_path, scan_root / saved_path.name]
#     for candidate in candidates:
#         if candidate.is_dir():
#             return candidate.resolve()
#     return None
#
#
# def _summary_for_path(
#     path: Path,
#     scan_root: Path,
#     cache: dict[Path, DirectoryPhaseSummary],
# ) -> DirectoryPhaseSummary | None:
#     directory = _resolve_data_directory(path, scan_root)
#     if directory is None:
#         return None
#     if directory not in cache:
#         cache[directory] = summarise_data_directory(directory)
#     return cache[directory]
#
#
# def _format_phases(summary: DirectoryPhaseSummary | None) -> str:
#     if summary is None:
#         return "directory not found"
#     if not summary.phases:
#         return "phase unavailable"
#     return ", ".join(f"{phase:g} degrees" for phase in summary.phases)
#
#
# def report_bba_sessions(root: Path, cache: dict[Path, DirectoryPhaseSummary]) -> None:
#     """Report the energy states of the Orbit/DFS data chosen in BBA sessions."""
#
#     settings_files = sorted(root.rglob("correction_settings.json"))
#     if not settings_files:
#         print("\nNo BBA correction_settings.json files found.")
#         return
#
#     print("\nBBA Orbit/DFS phase check")
#     print("=" * 80)
#     for settings_file in settings_files:
#         try:
#             settings = json.loads(settings_file.read_text())
#             data_dirs = settings["data_dirs"]
#         except (OSError, ValueError, KeyError) as error:
#             print(f"{settings_file.parent.name}: cannot read settings ({error})")
#             continue
#
#         orbit_dir = data_dirs.get("traj")
#         dfs_dir = data_dirs.get("dfs")
#         if not orbit_dir or not dfs_dir:
#             continue
#
#         orbit = _summary_for_path(Path(orbit_dir), root, cache)
#         dfs = _summary_for_path(Path(dfs_dir), root, cache)
#         orbit_phases = set(orbit.phases) if orbit else set()
#         dfs_phases = set(dfs.phases) if dfs else set()
#         dfs_weight = float(settings.get("w2", 0.0))
#
#         if dfs_weight <= 0:
#             status = "DFS disabled for this BBA session"
#         elif len(orbit_phases) == len(dfs_phases) == 1:
#             status = (
#                 "OK: different recorded RF phases"
#                 if orbit_phases != dfs_phases
#                 else "INVALID FOR DFS: identical recorded RF phase"
#             )
#         else:
#             status = "CHECK MANUALLY: missing or mixed phase data"
#
#         print(f"\n{settings_file.parent.name}")
#         print(f"  DFS weight w2: {dfs_weight:g}")
#         print(f"  Orbit: {_format_phases(orbit)} ({orbit_dir})")
#         print(f"  DFS:   {_format_phases(dfs)} ({dfs_dir})")
#         print(f"  {status}")
#
#
# def main() -> None:
#     parser = argparse.ArgumentParser(
#         description="Read RF phase metadata from saved SysID pickle files."
#     )
#     parser.add_argument(
#         "root",
#         type=Path,
#         nargs="?",
#         default=DEFAULT_ROOT,
#         help=f"Directory to scan recursively (default: {DEFAULT_ROOT})",
#     )
#     parser.add_argument(
#         "--show-files",
#         action="store_true",
#         help="Print phase and timestamp for every DATA_*.pkl file, not only directory summaries.",
#     )
#     parser.add_argument(
#         "--bba-sessions",
#         action="store_true",
#         help="Also compare the Orbit and DFS directories referenced by saved BBA sessions.",
#     )
#     args = parser.parse_args()
#
#     root = args.root.expanduser().resolve()
#     if not root.is_dir():
#         parser.error(f"Directory does not exist: {root}")
#
#     directories = _directories_with_sysid_data(root)
#     if not directories:
#         print(f"No DATA_*.pkl files found below {root}")
#         return
#
#     print(f"Scanning {len(directories)} SysID data directory/directories below {root}")
#     cache: dict[Path, DirectoryPhaseSummary] = {}
#     for directory in directories:
#         cache[directory.resolve()] = summarise_data_directory(directory, args.show_files)
#
#     if args.bba_sessions:
#         report_bba_sessions(root, cache)
#
#
# if __name__ == "__main__":
#     main()
