#!/usr/bin/env python3
"""3D surface plot of F(theta) = -H(theta) over the full ROI [0,400]x[0,400].

Reads the full-ROI F(theta) .npy file produced by the Fortran driver
and renders a 3D surface plot without markers.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    output_dir = Path(__file__).resolve().parent.parent / "outputs"
    npy_path = output_dir / "task_3b1_F_full.npy"
    meta_path = output_dir / "task_3b1_full_meta.dat"
    plot_path = output_dir / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    if not npy_path.exists():
        print(f"Error: {npy_path} not found. Run the Fortran driver first.")
        sys.exit(1)

    F_values = np.load(npy_path)  # shape: (nx, ny)
    nx, ny = F_values.shape

    if meta_path.exists():
        meta = np.loadtxt(meta_path)
        x_min, x_max, y_min, y_max = meta[0], meta[1], meta[2], meta[3]
    else:
        x_min, y_min = 0.0, 0.0
        x_max, y_max = 400.0, 400.0

    x_coords = np.linspace(x_min, x_max, nx)
    y_coords = np.linspace(y_min, y_max, ny)
    X, Y = np.meshgrid(x_coords, y_coords)
    Z = F_values.T  # transpose to (ny, nx) for meshgrid

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot_surface(
        X, Y, Z,
        cmap="viridis",
        edgecolor="none",
        alpha=0.9,
        antialiased=True,
    )

    ax.set_xlabel("x [m]", fontsize=11)
    ax.set_ylabel("y [m]", fontsize=11)
    ax.set_zlabel(r"$F(\theta) = -H(\theta)$", fontsize=11)
    ax.set_title(
        "Task 3-B.1: F(theta) over Full ROI [0, 400] x [0, 400]",
        fontsize=13,
    )

    out_file = plot_path / "task_3b1_F_surface_3d.pdf"
    fig.savefig(out_file, dpi=150, bbox_inches="tight")
    print(f"Saved {out_file}")

    png_file = plot_path / "task_3b1_F_surface_3d.png"
    fig.savefig(png_file, dpi=150, bbox_inches="tight")
    print(f"Saved {png_file}")

    plt.show()


if __name__ == "__main__":
    main()
