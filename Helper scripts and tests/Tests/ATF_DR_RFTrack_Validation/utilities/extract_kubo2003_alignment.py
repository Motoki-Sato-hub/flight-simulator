"""Extract Kubo (2003) Fig. 1/2 main-magnet alignment points from its PDF.

The source article stores the two plots as vector circles.  Reading those
circles avoids an imprecise raster digitisation.  The points are associated in
longitudinal order with the 204 quadrupoles, sextupoles, and bends of a chosen
compatible SAD RING0 export.  This association is documented in the JSON and
must not be described as an original 2003 SAD alignment database.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import zlib

import matplotlib.pyplot as plt
import numpy as np

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import get_magnet_centres


PAPER_CITATION = "K. Kubo, Phys. Rev. ST Accel. Beams 6, 092801 (2003)"


def _circle_centres(pdf_path: Path) -> list[tuple[float, float]]:
    """Read filled four-Bézier circles from the PDF content stream of page 3."""
    payload = pdf_path.read_bytes()
    page_object = re.search(rb"(?m)^3 0 obj\s*(.*?)\s*endobj", payload, re.S)
    if page_object is None:
        raise ValueError("Expected Kubo paper page-3 PDF object was not found")
    stream = re.search(
        rb"stream\r?\n(.*?)\r?\nendstream", page_object.group(1), re.S
    )
    if stream is None:
        raise ValueError("Expected compressed page-3 content stream was not found")
    content = zlib.decompress(stream.group(1)).decode("latin1")
    circles = re.compile(
        r"(?m)^([\d.]+) ([\d.]+) m\n(?:[^\n]* c\n){4}f\n"
    )
    result = []
    for match in circles.finditer(content):
        lines = match.group(0).splitlines()
        right_x, y = map(float, lines[0].split()[:2])
        left_x = float(lines[2].split()[4])
        result.append(((left_x + right_x) / 2.0, y))
    return result


def _figure_points(circles, *, x_range, y_range):
    selected = [
        point
        for point in circles
        if x_range[0] <= point[0] <= x_range[1] and y_range[0] <= point[1] <= y_range[1]
    ]
    if len(selected) != 204:
        raise ValueError(f"Expected 204 Fig. points, found {len(selected)}")
    return sorted(selected)


def extract(pdf_path: Path, lattice_data_path: Path) -> dict:
    circles = _circle_centres(pdf_path)
    # PDF plot coordinates, calibrated from labelled axes.  PDF y rises upward.
    horizontal = _figure_points(circles, x_range=(90.0, 290.0), y_range=(115.0, 230.0))
    vertical = _figure_points(circles, x_range=(350.0, 550.0), y_range=(625.0, 745.0))
    x_s0, x_s_per_pdf = 94.080, 20.0 / (121.613 - 94.080)
    y_s0, y_s_per_pdf = 357.078, 20.0 / (384.611 - 357.078)
    x_zero, x_mm_per_pdf = 164.437, 0.1 / (179.182 - 164.437)
    y_zero, y_mm_per_pdf = 686.067, 0.1 / (703.347 - 686.067)
    horizontal_data = [
        ((x - x_s0) * x_s_per_pdf, (y - x_zero) * x_mm_per_pdf)
        for x, y in horizontal
    ]
    vertical_data = [
        ((x - y_s0) * y_s_per_pdf, (y - y_zero) * y_mm_per_pdf)
        for x, y in vertical
    ]
    magnets = get_magnet_centres(lattice_data_path)
    if len(magnets) != 204:
        raise ValueError(f"Expected 204 main magnets in SAD export, found {len(magnets)}")

    records = []
    for magnet, (sx, dx_mm), (sy, dy_mm) in zip(magnets, horizontal_data, vertical_data):
        records.append(
            {
                **magnet,
                "figure_s_x_m": sx,
                "figure_s_y_m": sy,
                "s_x_residual_m": float(magnet["s_m"] - sx),
                "s_y_residual_m": float(magnet["s_m"] - sy),
                "dx_mm": dx_mm,
                "dy_mm": dy_mm,
            }
        )
    return {
        "metadata": {
            "citation": PAPER_CITATION,
            "source": "Kubo 2003 Fig. 1 (horizontal) and Fig. 2 (vertical)",
            "pdf_sha256": hashlib.sha256(pdf_path.read_bytes()).hexdigest(),
            "extraction": "vector-PDF filled-circle centres; labelled-axis linear calibration",
            "association": "longitudinal-order association to compatible SAD RING0 main magnets",
            "lattice_data": str(lattice_data_path),
            "warning": "Digitised publication data, not an original 2003 SAD alignment file.",
        },
        "magnets": records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--lattice-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()
    result = extract(args.pdf, args.lattice_data)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    residuals = [abs(row["s_x_residual_m"]) for row in result["magnets"]]
    if args.figure:
        rows = result["magnets"]
        s = np.asarray([row["s_m"] for row in rows])
        figure, axes = plt.subplots(2, 1, figsize=(8.0, 5.2), sharex=True, constrained_layout=True)
        axes[0].plot(s, [row["dx_mm"] for row in rows], ".", label="dx")
        axes[0].plot(s, [row["dy_mm"] for row in rows], ".", label="dy")
        axes[0].set(ylabel="published offset [mm]", title="Kubo (2003) Fig. 1/2 vector digitisation")
        axes[0].legend(ncol=2)
        axes[1].plot(s, [row["s_x_residual_m"] for row in rows], ".", color="tab:purple")
        axes[1].axhline(0.0, color="0.5", linewidth=0.8)
        axes[1].set(xlabel="associated lattice s [m]", ylabel="s association residual [m]")
        for axis in axes:
            axis.grid(alpha=0.25)
        args.figure.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(args.figure, dpi=180)
    print(f"Wrote {args.output} ({len(result['magnets'])} magnets)")
    print(f"Maximum longitudinal association residual: {max(residuals):.4f} m")


if __name__ == "__main__":
    main()
