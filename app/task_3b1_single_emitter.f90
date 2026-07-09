!> Task 3-B.1: Single-emitter convergence verification.
!>
!> Places one emitter at [191, 191] (off-grid for Delta = 20 m),
!> runs AREM, verifies that F(theta) = -H(theta) attains its maximum
!> at the true emitter location, and outputs F(theta) over a dense
!> 2D patch for 3D visualization.
program task_3b1_single_emitter
    use mod_precision, only: wp
    use mod_constants, only: PI
    use mod_types, only: grid_spec, emitter_set, sensor_set, arem_config, arem_result
    use mod_pathloss, only: build_path_loss_dictionary, path_loss
    use mod_preprocess, only: compute_preprocessing
    use mod_arem, only: arem_solve
    use mod_simulation, only: create_grid, generate_sensors, &
        generate_measurements
    use mod_gradient, only: compute_H
    use mod_io, only: write_array_2d, write_vector, ensure_output_dir
    implicit none

    ! --- Parameters ---
    real(wp), parameter :: ROI_SIZE = 400.0_wp
    real(wp), parameter :: DELTA   = 20.0_wp
    integer,  parameter :: M       = 200
    real(wp), parameter :: ETA     = 2.0_wp
    real(wp), parameter :: D0      = 5.0_wp
    real(wp), parameter :: EMITTER_X = 191.0_wp
    real(wp), parameter :: EMITTER_Y = 191.0_wp
    real(wp), parameter :: POWER_DBM = 0.0_wp

    ! Local variables
    type(grid_spec)    :: grid
    type(sensor_set)   :: sensors
    type(emitter_set)  :: emitters
    type(arem_config)  :: config
    type(arem_result)  :: result
    real(wp), allocatable :: z(:), A_grid(:, :), T(:, :), Aprime(:, :)
    real(wp), allocatable :: F_values(:, :)
    integer :: i, j, idx, nx_patch, ny_patch
    real(wp) :: x, y, F_val
    type(grid_spec) :: patch_grid

    ! --- Initialization ---
    call ensure_output_dir()
    write (*, "(A)") "=== Task 3-B.1: Single-Emitter Convergence ==="

    ! Create grid
    call create_grid(0.0_wp, ROI_SIZE, 0.0_wp, ROI_SIZE, DELTA, grid)
    write (*, "(A, I0, A, I0, A, I0)") "Grid: ", grid%nx, " x ", grid%ny, &
        " = ", grid%N, " points"

    ! Generate sensors (deterministic seed)
    call generate_sensors(0.0_wp, ROI_SIZE, 0.0_wp, ROI_SIZE, M, 42, sensors)
    write (*, "(A, I0)") "Sensors: M = ", sensors%M

    ! Place single emitter
    emitters%K = 1
    allocate (emitters%positions(1, 2))
    allocate (emitters%powers(1))
    emitters%positions(1, 1) = EMITTER_X
    emitters%positions(1, 2) = EMITTER_Y
    emitters%powers(1) = POWER_DBM
    write (*, "(A, F7.1, A, F7.1)") "True emitter at: ", EMITTER_X, ", ", EMITTER_Y

    ! Generate noiseless measurements
    allocate (z(M))
    call generate_measurements(sensors, emitters, D0, ETA, &
                               0.0_wp, 0.0_wp, 123, z)
    write (*, "(A, ES12.5)") "||z||_2 = ", norm2(z)

    ! Build full dictionary and preprocessing
    A_grid = build_path_loss_dictionary(sensors%positions, grid, D0, ETA)
    write (*, "(A, I0, A, I0)") "Dictionary A: ", size(A_grid, 1), &
        " x ", size(A_grid, 2)

    call compute_preprocessing(A_grid, T, Aprime)
    write (*, "(A)") "Preprocessing T and A'' computed."

    ! --- Run AREM ---
    config%eta = ETA
    config%d0 = D0
    config%th_error = 1e-10_wp
    config%alpha = 0.5_wp
    config%max_outer_iter = 3
    config%max_inner_iter = 200
    config%prune_threshold = 1e-6_wp
    config%use_double_grid = .false.

    call arem_solve(z, A_grid, T, Aprime, sensors%positions, grid, config, result)

    write (*, "(A, I0)") "Emitters found: ", result%n_found
    write (*, "(A, L1)") "Converged: ", result%converged

    if (result%n_found > 0) then
        write (*, "(A, F10.2, A, F10.2)") "Estimated position: ", &
            result%positions(1, 1), ", ", result%positions(1, 2)
        write (*, "(A, F10.6)") "Estimated power: ", result%powers(1)
        write (*, "(A, F10.4)") "Position error [m]: ", &
            norm2(result%positions(1, :) - [EMITTER_X, EMITTER_Y])
    end if

    ! --- Compute F(theta) = -H(theta) over a dense patch around the emitter ---
    ! Patch: 40x40 m around emitter, 1 m resolution
    nx_patch = 41; ny_patch = 41
    allocate (F_values(nx_patch, ny_patch))

    write (*, "(A)") "Computing F(theta) surface..."
    do j = 1, ny_patch
        y = EMITTER_Y - 20.0_wp + (j - 1) * 1.0_wp
        do i = 1, nx_patch
            x = EMITTER_X - 20.0_wp + (i - 1) * 1.0_wp

            ! Build A_Tk for this single position
            block
                real(wp), allocatable :: A_Tk(:, :), support_pos(:, :)
                allocate (support_pos(1, 2))
                support_pos(1, 1) = x
                support_pos(1, 2) = y
                A_Tk = build_A_Tk_patch(sensors%positions, support_pos, D0, ETA)
                F_val = -compute_H(z, A_Tk)
            end block

            F_values(i, j) = F_val
        end do
    end do

    ! Write F(theta) surface data
    call write_array_2d("task_3b1_F_surface.dat", F_values)
    write (*, "(A)") "F(theta) surface written to outputs/task_3b1_F_surface.dat"

    ! Write patch metadata
    block
        real(wp) :: meta(4)
        meta = [EMITTER_X - 20.0_wp, EMITTER_X + 20.0_wp, &
                EMITTER_Y - 20.0_wp, EMITTER_Y + 20.0_wp]
        call write_vector("task_3b1_patch_bounds.dat", meta)
    end block

    write (*, "(A)") "=== Task 3-B.1 complete ==="

contains

    !> Build a support dictionary for a given set of positions.
    !> (Local copy to avoid circular module dependencies in the driver.)
    function build_A_Tk_patch(sens, supp, d0_loc, eta_loc) result(A_Tk)
        real(wp), intent(in) :: sens(:, :)
        real(wp), intent(in) :: supp(:, :)
        real(wp), intent(in) :: d0_loc, eta_loc
        real(wp), allocatable :: A_Tk(:, :)

        integer  :: ii, mm, M_sens, n_supp
        real(wp) :: dx, dy, dist

        M_sens = size(sens, 1)
        n_supp = size(supp, 1)
        allocate (A_Tk(M_sens, n_supp))

        do ii = 1, n_supp
            do mm = 1, M_sens
                dx = sens(mm, 1) - supp(ii, 1)
                dy = sens(mm, 2) - supp(ii, 2)
                dist = sqrt(dx**2 + dy**2)
                A_Tk(mm, ii) = path_loss(dist, d0_loc, eta_loc)
            end do
        end do
    end function build_A_Tk_patch

end program task_3b1_single_emitter
