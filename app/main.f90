!> Compilation smoke test — verifies all core modules compile.
!> Individual numerical tasks are run via the dedicated drivers in app/.
program main
    use mod_precision, only: wp
    use mod_constants, only: PI, DEFAULT_D0, DEFAULT_ETA
    use mod_types, only: point_2d, grid_spec, arem_config
    use mod_pathloss, only: path_loss
    use mod_preprocess, only: compute_pinv_via_svd
    use mod_gradient, only: compute_H
    use mod_arem, only: arem_solve
    use mod_simulation, only: generate_sensors, compute_snr
    use mod_metrics, only: position_error, awerr
    use mod_io, only: ensure_output_dir
    implicit none

    call ensure_output_dir()

    write (*, "(A)") "AREM (Accurate Radio Environment Map) project"
    write (*, "(A, F8.4)") "  PI          = ", PI
    write (*, "(A, F8.4)") "  DEFAULT_D0  = ", DEFAULT_D0
    write (*, "(A, F8.4)") "  DEFAULT_ETA = ", DEFAULT_ETA
    write (*, "(A)") "All core modules compiled successfully."

end program main
