!> Task 3-B.1: Single-emitter convergence verification.
!>
!> Places one emitter at [191, 191] (off-grid for Delta = 20 m),
!> runs AREM, verifies that F(theta) = -H(theta) attains its maximum
!> at the true emitter location, and outputs F(theta) surfaces for
!> 3D and 2D visualization.
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
    use mod_io, only: write_vector, ensure_output_dir
    use stdlib_io_npy, only: save_npy
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
    real(wp), allocatable :: F_full(:, :), F_local(:, :)
    real(wp), allocatable :: gamma_prime(:), corr(:)
    real(wp), allocatable :: A_Tk_patch(:, :), support_pos(:, :)
    integer  :: i, j, idx_max, nx_full, ny_full, nx_local, ny_local
    real(wp) :: x, y, F_val, rough_x, rough_y, refined_x, refined_y
    real(wp), parameter :: FULL_DELTA  = 10.0_wp  ! 10 m spacing for full ROI
    real(wp), parameter :: LOCAL_DELTA = 1.0_wp   ! 1 m spacing for local patch
    real(wp), parameter :: LOCAL_HALF  = 30.0_wp  ! half-width of local patch

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

    ! --- Rough estimate (Eq. 15) ---
    allocate (gamma_prime(size(T, 1)))
    allocate (corr(grid%N))
    gamma_prime = matmul(T, z)
    corr = abs(matmul(transpose(Aprime), gamma_prime))
    idx_max = maxloc(corr, dim=1)
    rough_x = grid%points(idx_max, 1)
    rough_y = grid%points(idx_max, 2)
    write (*, "(A, F7.1, A, F7.1)") "Rough estimate: ", rough_x, ", ", rough_y

    ! --- Run AREM ---
    config%eta = ETA
    config%d0 = D0
    config%th_error = 1e-10_wp
    config%alpha = 1.0_wp
    config%max_outer_iter = 3
    config%max_inner_iter = 200
    config%prune_threshold = 1e-6_wp
    config%use_double_grid = .false.

    call arem_solve(z, A_grid, T, Aprime, sensors%positions, grid, config, result)

    write (*, "(A, I0)") "Emitters found: ", result%n_found
    write (*, "(A, L1)") "Converged: ", result%converged

    if (result%n_found > 0) then
        refined_x = result%positions(1, 1)
        refined_y = result%positions(1, 2)
        write (*, "(A, F10.2, A, F10.2)") "Refined position:  ", refined_x, ", ", refined_y
        write (*, "(A, F10.6)") "Estimated power: ", result%powers(1)
        write (*, "(A, F10.4)") "Position error [m]: ", &
            norm2(result%positions(1, :) - [EMITTER_X, EMITTER_Y])
    else
        refined_x = rough_x
        refined_y = rough_y
    end if

    ! --- Save position data for plotting ---
    block
        real(wp) :: pos_data(6)
        pos_data = [EMITTER_X, EMITTER_Y, rough_x, rough_y, refined_x, refined_y]
        call write_vector("task_3b1_positions.dat", pos_data)
    end block

    allocate (support_pos(1, 2))

    ! --- Compute F(theta) over the FULL ROI at 10 m spacing ---
    nx_full = nint(ROI_SIZE / FULL_DELTA) + 1  ! 41
    ny_full = nint(ROI_SIZE / FULL_DELTA) + 1
    allocate (F_full(nx_full, ny_full))
    allocate (A_Tk_patch(M, 1))

    write (*, "(A, I0, A, I0, A)") "Computing F(theta) over full ROI (", &
        nx_full, " x ", ny_full, ")..."
    do j = 1, ny_full
        y = (j - 1) * FULL_DELTA
        do i = 1, nx_full
            x = (i - 1) * FULL_DELTA
            support_pos(1, 1) = x
            support_pos(1, 2) = y
            A_Tk_patch = build_A_Tk_scan(sensors%positions, support_pos, D0, ETA)
            F_full(i, j) = -compute_H(z, A_Tk_patch)
        end do
    end do

    call save_npy("outputs/task_3b1_F_full.npy", F_full)
    write (*, "(A)") "Full-ROI F(theta) saved to outputs/task_3b1_F_full.npy"

    ! Save full-ROI metadata: x_min, x_max, y_min, y_max
    block
        real(wp) :: full_meta(4)
        full_meta = [0.0_wp, ROI_SIZE, 0.0_wp, ROI_SIZE]
        call write_vector("task_3b1_full_meta.dat", full_meta)
    end block

    ! --- Compute F(theta) over LOCAL patch at 1 m spacing ---
    nx_local = nint(2.0_wp * LOCAL_HALF / LOCAL_DELTA) + 1  ! 61
    ny_local = nint(2.0_wp * LOCAL_HALF / LOCAL_DELTA) + 1
    allocate (F_local(nx_local, ny_local))

    write (*, "(A, I0, A, I0, A)") "Computing F(theta) over local patch (", &
        nx_local, " x ", ny_local, ")..."
    do j = 1, ny_local
        y = EMITTER_Y - LOCAL_HALF + (j - 1) * LOCAL_DELTA
        do i = 1, nx_local
            x = EMITTER_X - LOCAL_HALF + (i - 1) * LOCAL_DELTA
            support_pos(1, 1) = x
            support_pos(1, 2) = y
            A_Tk_patch = build_A_Tk_scan(sensors%positions, support_pos, D0, ETA)
            F_local(i, j) = -compute_H(z, A_Tk_patch)
        end do
    end do

    call save_npy("outputs/task_3b1_F_local.npy", F_local)
    write (*, "(A)") "Local-patch F(theta) saved to outputs/task_3b1_F_local.npy"

    ! Save local-patch metadata
    block
        real(wp) :: local_meta(4)
        local_meta = [EMITTER_X - LOCAL_HALF, EMITTER_X + LOCAL_HALF, &
                      EMITTER_Y - LOCAL_HALF, EMITTER_Y + LOCAL_HALF]
        call write_vector("task_3b1_local_meta.dat", local_meta)
    end block

    write (*, "(A)") "=== Task 3-B.1 complete ==="

contains

    !> Build a single-column dictionary for F(theta) scanning.
    function build_A_Tk_scan(sens, supp, d0_loc, eta_loc) result(A_Tk)
        real(wp), intent(in) :: sens(:, :)
        real(wp), intent(in) :: supp(:, :)
        real(wp), intent(in) :: d0_loc, eta_loc
        real(wp), allocatable :: A_Tk(:, :)

        integer  :: mm, M_sens
        real(wp) :: dx, dy, dist

        M_sens = size(sens, 1)
        allocate (A_Tk(M_sens, 1))

        do mm = 1, M_sens
            dx = sens(mm, 1) - supp(1, 1)
            dy = sens(mm, 2) - supp(1, 2)
            dist = sqrt(dx**2 + dy**2)
            A_Tk(mm, 1) = path_loss(dist, d0_loc, eta_loc)
        end do
    end function build_A_Tk_scan

end program task_3b1_single_emitter
