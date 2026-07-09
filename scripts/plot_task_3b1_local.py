#!/usr/bin/env python3
"""2D heatmap (top view) of F(theta) over a local region around the true emitter.

Reads the local-patch F(theta) .npy file and position data produced by
the Fortran driver.  Marks:
  - True emitter position      : red filled circle
  - Rough estimate (on-grid)   : black filled circle
  - Refined estimate (off-grid): blue hollow triangle

The visual effect shows the true position (red dot) sitting inside
the refined estimate's blue triangle.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    output_dir = Path(__file__).resolve().parent.parent / "outputs"
    npy_path = output_dir / "task_3b1_F_local.npy"
    meta_path = output_dir / "task_3b1_local_meta.dat"
    pos_path = output_dir / "task_3b1_positions.dat"
    plot_path = output_dir / "plots"
    plot_path.mkdir(parents=True, exist_ok=True)

    if not npy_path.exists():
        print(f"Error: {npy_path} not found. Run the Fortran driver first.")
        sys.exit(1)

    F_values = np.load(npy_path)  # shape: (nx, ny)
    nx, ny = F_values.shape

    # Load local-patch bounds
    if meta_path.exists():
        meta = np.loadtxt(meta_path)
        x_min, x_max, y_min, y_max = meta[0], meta[1], meta[2], meta[3]
    else:
        x_min, x_max = 161.0, 221.0
        y_min, y_max = 161.0, 221.0

    # Load positions: [true_x, true_y, rough_x, rough_y, refined_x, refined_y]
    if pos_path.exists():
        pos = np.loadtxt(pos_path)
        true_x, true_y = pos[0], pos[1]
        rough_x, rough_y = pos[2], pos[3]
        refined_x, refined_y = pos[4], pos[5]
    else:
        true_x = true_y = 191.0
        rough_x, rough_y = 180.0, 200.0
        refined_x, refined_y = 191.0, 191.0

    # Coordinate arrays (F_values stored x-fast, y-slow; transpose for imshow)
    x_coords = np.linspace(x_min, x_max, nx)
    y_coords = np.linspace(y_min, y_max, ny)
    Z = F_values.T  # (ny, nx)

    # --- 2D heatmap ---
    fig, ax = plt.subplots(figsize=(9, 8))

    extent = [x_min, x_max, y_min, y_max]
    im = ax.imshow(
        Z,
        origin="lower",
        extent=extent,
        cmap="viridis",
        aspect="equal",
        interpolation="bilinear",
    )

    # Markers
    ax.plot(
        true_x, true_y, "o",
        color="red", markersize=10, markeredgewidth=1.5,
        markeredgecolor="darkred",
        label="True emitter",
    )
    ax.plot(
        rough_x, rough_y, "o",
        color="black", markersize=10,
        label="Rough estimate (on-grid)",
    )
    ax.plot(
        refined_x, refined_y, "^",
        color="blue", markersize=14, markeredgewidth=2.0,
        markerfacecolor="none",
        label="Refined estimate (off-grid)",
    )

    ax.set_xlabel("x [m]", fontsize=12)
    ax.set_ylabel("y [m]", fontsize=12)
    ax.set_title(
        "Task 3-B.1: F(theta) — Local Region with Estimates",
        fontsize=13,
    )
    ax.legend(fontsize=10, loc="lower right")

    cbar = fig.colorbar(im, ax=ax, shrink=0.85, aspect=15)
    cbar.set_label(r"$F(\theta) = -H(\theta)$", fontsize=11)

    # Save
    out_pdf = plot_path / "task_3b1_F_local_2d.pdf"
    fig.savefig(out_pdf, dpi=150, bbox_inches="tight")
    print(f"Saved {out_pdf}")

    out_png = plot_path / "task_3b1_F_local_2d.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"Saved {out_png}")

    plt.show()


if __name__ == "__main__":
    main()
