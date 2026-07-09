!> Simulation framework for the AREM numerical tasks.
!>
!> Provides procedures to generate random sensor/emitter placements,
!> synthetic RSS measurements (with shadowing and noise), SNR control,
!> and RSS map computation.
!>
!> Uses `stdlib` random number generators for reproducibility.
module mod_simulation
    use mod_precision, only: wp
    use mod_types, only: grid_spec, emitter_set, sensor_set, arem_result, rss_map
    use mod_pathloss, only: path_loss
    use stdlib_stats_distribution_uniform, only: rvs_uniform
    use stdlib_stats_distribution_normal, only: rvs_normal
    implicit none
    private

    public :: create_grid
    public :: generate_sensors
    public :: generate_emitters
    public :: generate_measurements
    public :: compute_rss_map
    public :: compute_rss_map_from_arem
    public :: compute_snr
    public :: add_noise_by_snr

contains

    !> Initialize a regular grid over a rectangular ROI.
    subroutine create_grid(x_min, x_max, y_min, y_max, delta, grid)
        !> ROI x bounds [m]
        real(wp), intent(in) :: x_min, x_max
        !> ROI y bounds [m]
        real(wp), intent(in) :: y_min, y_max
        !> Grid spacing [m]
        real(wp), intent(in) :: delta
        !> Output grid specification
        type(grid_spec), intent(out) :: grid

        integer :: i, j, idx

        grid%x_min = x_min
        grid%x_max = x_max
        grid%y_min = y_min
        grid%y_max = y_max
        grid%delta = delta

        ! Number of grid points in each direction
        grid%nx = nint((x_max - x_min) / delta) + 1
        grid%ny = nint((y_max - y_min) / delta) + 1
        grid%N = grid%nx * grid%ny

        allocate (grid%points(grid%N, 2))

        idx = 0
        do j = 0, grid%ny - 1
            do i = 0, grid%nx - 1
                idx = idx + 1
                grid%points(idx, 1) = x_min + i * delta
                grid%points(idx, 2) = y_min + j * delta
            end do
        end do
    end subroutine create_grid

    !> Generate M sensor positions uniformly distributed in the ROI.
    subroutine generate_sensors(x_min, x_max, y_min, y_max, M, seed, sensors)
        !> ROI x lower bound [m]
        real(wp), intent(in) :: x_min, x_max
        !> ROI y bounds [m]
        real(wp), intent(in) :: y_min, y_max
        !> Number of sensors
        integer,  intent(in) :: M
        !> Random seed for reproducibility
        integer,  intent(in) :: seed
        !> Output sensor set
        type(sensor_set), intent(out) :: sensors

        call set_random_seed(seed)

        sensors%M = M
        allocate (sensors%positions(M, 2))
        sensors%positions(:, 1) = rvs_uniform(x_min, x_max - x_min, M)
        sensors%positions(:, 2) = rvs_uniform(y_min, y_max - y_min, M)
    end subroutine generate_sensors

    !> Generate K emitter positions at 1m resolution with uniform random powers.
    !>
    !> Positions are rounded to integer coordinates (1 m resolution),
    !> matching the numerical task specification in Section 3-B.2.
    subroutine generate_emitters(x_min, x_max, y_min, y_max, K, seed, &
                                  power_min, power_max, emitters)
        !> ROI x bounds [m]
        real(wp), intent(in) :: x_min, x_max
        !> ROI y bounds [m]
        real(wp), intent(in) :: y_min, y_max
        !> Number of emitters
        integer,  intent(in) :: K
        !> Random seed
        integer,  intent(in) :: seed
        !> Minimum transmit power [dBm]
        real(wp), intent(in) :: power_min
        !> Maximum transmit power [dBm]
        real(wp), intent(in) :: power_max
        !> Output emitter set
        type(emitter_set), intent(out) :: emitters

        integer :: i

        call set_random_seed(seed)

        emitters%K = K
        allocate (emitters%positions(K, 2))
        allocate (emitters%powers(K))

        do i = 1, K
            ! Generate random position and round to 1 m grid
            emitters%positions(i, 1) = nint(rvs_uniform(x_min, x_max - x_min))
            emitters%positions(i, 2) = nint(rvs_uniform(y_min, y_max - y_min))
            ! Generate random power in [power_min, power_max] dBm
            emitters%powers(i) = rvs_uniform(power_min, power_max - power_min)
        end do
    end subroutine generate_emitters

    !> Generate synthetic RSS measurements (Eq. 1).
    !>
    !> Power values are in dBm; converted to linear for the measurement model:
    !> p_lin = 10^(p_dBm / 10).
    subroutine generate_measurements(sensors, emitters, d0, eta, &
                                      sigma_noise, sigma_shadowing, seed, z)
        !> Sensor positions and count
        type(sensor_set), intent(in) :: sensors
        !> Emitter positions and powers (in dBm)
        type(emitter_set), intent(in) :: emitters
        !> Reference distance [m]
        real(wp), intent(in) :: d0
        !> Path-loss exponent
        real(wp), intent(in) :: eta
        !> Standard deviation of measurement noise [linear]
        real(wp), intent(in) :: sigma_noise
        !> Standard deviation of log-normal shadowing [dB]
        real(wp), intent(in) :: sigma_shadowing
        !> Random seed
        integer,  intent(in) :: seed
        !> Output measurements z [M]
        real(wp), intent(out) :: z(:)

        integer  :: m, k, M_sensors, K_emitters
        real(wp) :: dx, dy, dist, p_linear, signal_sum
        real(wp), allocatable :: shadowing(:), noise(:)

        M_sensors = sensors%M
        K_emitters = emitters%K

        allocate (shadowing(M_sensors), noise(M_sensors))
        call set_random_seed(seed)

        ! Log-normal shadowing: exp(normal(0, sigma_shadowing))
        ! Guard against zero std dev (degenerate normal)
        if (sigma_shadowing > 0.0_wp) then
            shadowing = exp(rvs_normal(0.0_wp, sigma_shadowing, M_sensors))
        else
            shadowing = 1.0_wp
        end if

        ! Gaussian measurement noise
        if (sigma_noise > 0.0_wp) then
            noise = rvs_normal(0.0_wp, sigma_noise, M_sensors)
        else
            noise = 0.0_wp
        end if

        do m = 1, M_sensors
            signal_sum = 0.0_wp
            do k = 1, K_emitters
                dx = sensors%positions(m, 1) - emitters%positions(k, 1)
                dy = sensors%positions(m, 2) - emitters%positions(k, 2)
                dist = sqrt(dx**2 + dy**2)

                ! Convert dBm to linear power
                p_linear = 10.0_wp**(emitters%powers(k) / 10.0_wp)

                signal_sum = signal_sum + p_linear * path_loss(dist, d0, eta)
            end do
            z(m) = shadowing(m) * signal_sum + noise(m)
        end do
    end subroutine generate_measurements

    !> Compute the ground-truth RSS map (received power at each grid point in dBm).
    subroutine compute_rss_map(grid, emitters, d0, eta, P_map)
        !> Grid specification
        type(grid_spec), intent(in) :: grid
        !> True emitter parameters
        type(emitter_set), intent(in) :: emitters
        !> Reference distance [m]
        real(wp), intent(in) :: d0
        !> Path-loss exponent
        real(wp), intent(in) :: eta
        !> Output RSS map
        type(rss_map), intent(out) :: P_map

        integer  :: i, k, N_grid, K_emitters
        real(wp) :: dx, dy, dist, p_linear, signal_sum

        N_grid = grid%N
        K_emitters = emitters%K

        allocate (P_map%values(N_grid))
        allocate (P_map%weights(N_grid))

        do i = 1, N_grid
            signal_sum = 0.0_wp
            do k = 1, K_emitters
                dx = grid%points(i, 1) - emitters%positions(k, 1)
                dy = grid%points(i, 2) - emitters%positions(k, 2)
                dist = sqrt(dx**2 + dy**2)

                p_linear = 10.0_wp**(emitters%powers(k) / 10.0_wp)
                signal_sum = signal_sum + p_linear * path_loss(dist, d0, eta)
            end do

            if (signal_sum > 0.0_wp) then
                P_map%values(i) = 10.0_wp * log10(signal_sum)
            else
                P_map%values(i) = -200.0_wp  ! effectively -inf dBm
            end if
        end do

        ! Compute AWErr weights: lambda_i proportional to linear power
        call compute_awerr_weights(P_map)
    end subroutine compute_rss_map

    !> Compute RSS map from AREM estimation results.
    subroutine compute_rss_map_from_arem(grid, result, d0, eta, P_est)
        !> Grid specification
        type(grid_spec), intent(in) :: grid
        !> AREM estimation results
        type(arem_result), intent(in) :: result
        !> Reference distance [m]
        real(wp), intent(in) :: d0
        !> Path-loss exponent
        real(wp), intent(in) :: eta
        !> Output estimated RSS map
        type(rss_map), intent(out) :: P_est

        integer  :: i, k, N_grid, K_est
        real(wp) :: dx, dy, dist, p_linear, signal_sum

        N_grid = grid%N
        K_est = result%n_found

        allocate (P_est%values(N_grid))
        allocate (P_est%weights(N_grid))

        if (K_est == 0) then
            P_est%values = -200.0_wp
            P_est%weights = 1.0_wp / N_grid
            return
        end if

        do i = 1, N_grid
            signal_sum = 0.0_wp
            do k = 1, K_est
                dx = grid%points(i, 1) - result%positions(k, 1)
                dy = grid%points(i, 2) - result%positions(k, 2)
                dist = sqrt(dx**2 + dy**2)

                p_linear = 10.0_wp**(result%powers(k) / 10.0_wp)
                signal_sum = signal_sum + p_linear * path_loss(dist, d0, eta)
            end do

            if (signal_sum > 0.0_wp) then
                P_est%values(i) = 10.0_wp * log10(signal_sum)
            else
                P_est%values(i) = -200.0_wp
            end if
        end do
    end subroutine compute_rss_map_from_arem

    !> Compute SNR in dB from signal and noise powers.
    pure real(wp) function compute_snr(signal_power, noise_power) result(snr_dB)
        !> Signal power (linear)
        real(wp), intent(in) :: signal_power
        !> Noise power (linear)
        real(wp), intent(in) :: noise_power

        snr_dB = 10.0_wp * log10(signal_power / noise_power)
    end function compute_snr

    !> Add noise to achieve a target SNR.
    !>
    !> Computes signal power from z_true, then scales Gaussian noise
    !> to achieve target_snr in dB.
    subroutine add_noise_by_snr(z_true, target_snr_dB, seed, z_noisy)
        !> Clean measurement vector [M]
        real(wp), intent(in) :: z_true(:)
        !> Target SNR [dB]
        real(wp), intent(in) :: target_snr_dB
        !> Random seed
        integer,  intent(in) :: seed
        !> Noisy measurement vector [M]
        real(wp), intent(out) :: z_noisy(:)

        integer  :: M_sensors
        real(wp) :: signal_power, noise_power, noise_std
        real(wp), allocatable :: noise(:)

        M_sensors = size(z_true)

        call set_random_seed(seed)

        ! Signal power = mean of squared measurements
        signal_power = sum(z_true**2) / M_sensors

        ! Required noise power for target SNR:
        ! SNR_dB = 10*log10(P_signal / P_noise) => P_noise = P_signal / 10^(SNR/10)
        noise_power = signal_power / 10.0_wp**(target_snr_dB / 10.0_wp)
        noise_std = sqrt(noise_power)

        noise = rvs_normal(0.0_wp, noise_std, M_sensors)
        z_noisy = z_true + noise
    end subroutine add_noise_by_snr

    !> Compute AWErr weights: lambda_i proportional to linear RSS power,
    !> normalized to sum to 1.
    subroutine compute_awerr_weights(P_map)
        !> RSS map (weights field filled in-place)
        type(rss_map), intent(inout) :: P_map

        integer  :: i, N_grid
        real(wp) :: total_linear
        real(wp), allocatable :: p_linear(:)

        N_grid = size(P_map%values)
        allocate (p_linear(N_grid))

        do i = 1, N_grid
            if (P_map%values(i) > -199.0_wp) then
                p_linear(i) = 10.0_wp**(P_map%values(i) / 10.0_wp)
            else
                p_linear(i) = 0.0_wp
            end if
        end do

        total_linear = sum(p_linear)
        if (total_linear > 0.0_wp) then
            P_map%weights = p_linear / total_linear
        else
            P_map%weights = 1.0_wp / N_grid
        end if
    end subroutine compute_awerr_weights

    !> Set the Fortran intrinsic random seed from a scalar integer,
    !> ensuring reproducibility across runs.
    subroutine set_random_seed(seed)
        !> Scalar seed value
        integer, intent(in) :: seed

        integer              :: n, i
        integer, allocatable :: seed_arr(:)

        call random_seed(size=n)
        allocate (seed_arr(n))
        ! Fill the seed array with a deterministic sequence from the scalar
        do i = 1, n
            seed_arr(i) = seed + i - 1
        end do
        call random_seed(put=seed_arr)
    end subroutine set_random_seed

end module mod_simulation
