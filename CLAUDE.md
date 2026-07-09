# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

AREM (Accurate Radio Environment Map) — a Fortran implementation of a sparse-recovery/dictionary-learning algorithm for emitter localization and power estimation from RSS measurements. Implements Algorithm 1 from `paper-notes.tex` and its 10 numerical tasks (Sections 3-B, 3-C, 4-A).

Built with **fortran-lang/fpm**. Dependencies: `blas` (mapped to system OpenBLAS), `openmp`. The `stdlib` dependency was originally declared but the line is currently commented in `fpm.toml`; the project currently uses `stdlib_linalg`, `stdlib_stats_distribution_uniform`, `stdlib_stats_distribution_normal` modules — if these fail to resolve, uncomment `stdlib` in `fpm.toml`.

## Build & Run Commands

```bash
fpm build                          # Build all modules and drivers
fpm run                            # Run the main smoke-test executable
fpm run task_3b1_single_emitter    # Run a specific app driver
fpm test                           # Run tests
fpm build --compiler gfortran      # Specify compiler
```

All app drivers are auto-discovered from `app/*.f90` (`auto-executables = true`). Output goes to `outputs/` for post-processing by Python scripts.

## Architecture

**Primary paradigm: Array + Modular** (no OOP — no component has 3+ interchangeable variants). Plain derived types for data containers; `elemental pure` functions for path-loss math; structured subroutines for algorithmic pipelines.

### Module Dependency Hierarchy (bottom-up)

```
mod_precision (wp = real64)
  ├── mod_constants   — PI, DEFAULT_D0, DEFAULT_ETA, DEFAULT_GRID_DELTA, DEFAULT_ROI_SIZE
  ├── mod_types       — Plain derived types: grid_spec, emitter_set, sensor_set, arem_config, arem_result, rss_map
  ├── mod_io          — write_array_2d, write_vector, write_sweep_csv, ensure_output_dir → outputs/
  ├── mod_pathloss    — elemental: path_loss(), path_loss_deriv_x/y(); build_path_loss_vector(), build_path_loss_dictionary()
  ├── mod_preprocess  — compute_preprocessing(A, T, Aprime): thin-SVD whitening T = diag(1/S)·Uᵀ, A' = T·A
  ├── mod_gradient    — compute_H(), compute_gradient() via stdlib_linalg::lstsq and stdlib_linalg::svd
  ├── mod_arem        — arem_solve(): outer loop (rough est → support update) + inner alternating loop (LS omega → prune → GD position → residual)
  ├── mod_simulation  — create_grid(), generate_sensors/emitters(), generate_measurements(), compute_rss_map(), add_noise_by_snr()
  └── mod_metrics     — position_error(), power_error(), success_probability(), awerr()
```

### Preprocessing: SVD-based whitening, not QR

The `compute_preprocessing` in `mod_preprocess` uses thin SVD (`stdlib_linalg::svd`): `T = diag(1/s_i)·U^T` (r×M) and `A' = V^T` (r×N) where r = min(M,N). This produces decorrelated dictionary columns. Unlike the QR-based approach in the paper, SVD-whitening is dimensionally correct when M < N (fewer sensors than grid points). `T` must be passed to `arem_solve()` so the residual can be preprocessed: `gamma' = T @ gamma`.

### AREM Algorithm Flow

`arem_solve(z, A_grid, T, Aprime, sensors, grid, config, result)`:
1. **Outer loop** — adds one emitter per iteration:
   - Preprocess residual: `gamma' = T @ gamma`
   - Rough estimation: `idx = maxloc(abs(matmul(transpose(Aprime), gamma')))`
   - Update support set with new grid point
2. **Inner alternating loop** — refines positions/powers:
   - `omega = lstsq(A_Tk, z)` (LS power update)
   - Prune by threshold and ROI bounds
   - Gradient descent on `H(θ) = ||z - A_Tk·A_Tk†·z||²` via `mod_gradient`
   - Residual update and convergence check

### Key stdlib Usage Notes

- **`lstsq(A, b)` requires `A` to be `intent(inout)`** — the LAPACK backend may overwrite it. Always pass an allocated copy: `allocate(A_work, source=A); x = lstsq(A_work, b)`. Do not pass `intent(in)` arrays directly.
- **`rvs_normal(loc, 0.0_wp, n)` with zero std-dev produces undefined results** — guard with `if (sigma > 0)` and use constant `0.0_wp` or `1.0_wp` otherwise.
- **`random_seed` from stdlib conflicts with the Fortran intrinsic** — use the intrinsic directly via a helper (`set_random_seed` in `mod_simulation`) that queries the seed array size.

## Coding Conventions

From `coding_standards.md` (enforced for all new code):

- **Every procedure and derived type must live inside a module.** No bare procedures.
- **`implicit none` + `private`** at module level; `public :: ` lists all exported symbols.
- **`real(wp)` for all floats; `_wp` suffix on all literals.** `wp = real64` from `mod_precision`.
- **Every dummy argument has explicit `intent(in|out|inout)`** and its own `!>` doc-comment line above.
- **Lowercase** for all Fortran keywords, identifiers, and intrinsics.
- **`write(*, fmt) ...`** never `print *`; use explicit format strings.
- **`!>` for procedure/argument doc-comments; `!` for in-body comments.** Formal English throughout.
- **No OpenMP in initial code** — optimization is a later phase. The `openmp` dependency in `fpm.toml` is declared but should not be used yet.
- **No magic numbers** — every physical/numerical constant is named with units.

## Data Flow

```
Fortran app/ drivers → call arem_solve() + mod_metrics → write to outputs/*.dat, *.csv
                                                              ↓
                                              Python scripts/ read, plot → outputs/plots/*.pdf
```

## Key Types (mod_types)

| Type | Purpose |
|---|---|
| `grid_spec` | ROI bounds, Δ, N, nx, ny, allocatable `points(N,2)` |
| `arem_config` | th_error, alpha (GD step), max_outer/inner_iter, prune_threshold, use_double_grid, eta, d0 |
| `arem_result` | allocatable positions(:,:), powers(:), n_found, n_outer_iter, converged |
| `rss_map` | allocatable values(:) [dBm], weights(:) [sum=1, proportional to linear RSS] |

## App Drivers

Each is a standalone `program` in `app/`. Currently implemented: `task_3b1_single_emitter`. Remaining tasks (3-B.2, 3-C, 4-A.1–4-A.7) follow the same pattern: configure parameters → create grid via `mod_simulation::create_grid()` → build dictionary → `compute_preprocessing` → Monte Carlo loop calling `arem_solve()` → compute metrics → write output.
