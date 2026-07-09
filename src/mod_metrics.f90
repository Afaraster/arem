!> Performance metrics for the AREM numerical tasks.
!>
!> Implements position error (Eq. 23), power error (Eq. 24),
!> success probability (Eq. 25), and AWErr (Eq. 26).
module mod_metrics
    use mod_precision, only: wp
    implicit none
    private

    public :: position_error
    public :: power_error
    public :: success_probability
    public :: awerr

contains

    !> Compute average position estimation error (Eq. 23).
    !>
    !> Uses greedy nearest-neighbor matching to pair estimated emitters
    !> with true emitters. Unmatched emitters contribute the maximum
    !> possible error (ROI diagonal).
    function position_error(est_pos, true_pos) result(err)
        !> Estimated emitter positions [K_est, 2]
        real(wp), intent(in) :: est_pos(:, :)
        !> True emitter positions [K_true, 2]
        real(wp), intent(in) :: true_pos(:, :)
        !> Average position error [m]
        real(wp) :: err

        integer  :: K_est, K_true, k, j, best_j, n_matched
        real(wp) :: best_dist, dist, total_err
        logical, allocatable :: matched(:)

        K_est = size(est_pos, 1)
        K_true = size(true_pos, 1)

        if (K_est == 0 .or. K_true == 0) then
            err = huge(1.0_wp)
            return
        end if

        allocate (matched(K_est))
        matched = .false.
        total_err = 0.0_wp
        n_matched = 0

        do k = 1, K_true
            best_dist = huge(1.0_wp)
            best_j = -1
            do j = 1, K_est
                if (matched(j)) cycle
                dist = norm2(est_pos(j, :) - true_pos(k, :))
                if (dist < best_dist) then
                    best_dist = dist
                    best_j = j
                end if
            end do

            if (best_j > 0) then
                matched(best_j) = .true.
                total_err = total_err + best_dist
                n_matched = n_matched + 1
            end if
        end do

        if (n_matched > 0) then
            err = total_err / n_matched
        else
            err = huge(1.0_wp)
        end if
    end function position_error

    !> Compute average power estimation error (Eq. 24).
    !>
    !> Uses the same greedy matching as position_error to pair emitters.
    function power_error(est_pos, est_pwr, true_pos, true_pwr) result(err)
        !> Estimated emitter positions [K_est, 2]
        real(wp), intent(in) :: est_pos(:, :)
        !> Estimated emitter powers [K_est] in dBm
        real(wp), intent(in) :: est_pwr(:)
        !> True emitter positions [K_true, 2]
        real(wp), intent(in) :: true_pos(:, :)
        !> True emitter powers [K_true] in dBm
        real(wp), intent(in) :: true_pwr(:)
        !> Average power error (relative squared)
        real(wp) :: err

        integer  :: K_est, K_true, k, j, best_j, n_matched
        real(wp) :: best_dist, dist, total_err
        logical, allocatable :: matched(:)

        K_est = size(est_pos, 1)
        K_true = size(true_pos, 1)

        if (K_est == 0 .or. K_true == 0) then
            err = huge(1.0_wp)
            return
        end if

        allocate (matched(K_est))
        matched = .false.
        total_err = 0.0_wp
        n_matched = 0

        do k = 1, K_true
            best_dist = huge(1.0_wp)
            best_j = -1
            do j = 1, K_est
                if (matched(j)) cycle
                dist = norm2(est_pos(j, :) - true_pos(k, :))
                if (dist < best_dist) then
                    best_dist = dist
                    best_j = j
                end if
            end do

            if (best_j > 0) then
                matched(best_j) = .true.
                ! Relative squared error: ((p_est - p_true) / p_true)^2
                if (abs(true_pwr(k)) > 1e-10_wp) then
                    total_err = total_err &
                        + ((est_pwr(best_j) - true_pwr(k)) / true_pwr(k))**2
                end if
                n_matched = n_matched + 1
            end if
        end do

        if (n_matched > 0) then
            err = total_err / n_matched
        else
            err = huge(1.0_wp)
        end if
    end function power_error

    !> Compute probability of correctly counting the number of emitters (Eq. 25).
    pure real(wp) function success_probability(n_correct, n_total) result(Pr)
        !> Number of trials where emitter count was correct
        integer, intent(in) :: n_correct
        !> Total number of trials
        integer, intent(in) :: n_total

        if (n_total > 0) then
            Pr = real(n_correct, wp) / real(n_total, wp)
        else
            Pr = 0.0_wp
        end if
    end function success_probability

    !> Compute the Average Weighted Error (AWErr) for REM reconstruction
    !> (Eq. 26).
    !>
    !> $$\mathrm{AWErr} = \sum_{i=1}^N \lambda_i |P_i - \dot{P}_i|$$
    !>
    !> where weights sum to 1 and are proportional to the true RSS power.
    function awerr(P_true, P_est, weights) result(err)
        !> True RSS values [N] in dBm
        real(wp), intent(in) :: P_true(:)
        !> Estimated RSS values [N] in dBm
        real(wp), intent(in) :: P_est(:)
        !> AWErr weights [N] (sum to 1)
        real(wp), intent(in) :: weights(:)
        !> Average weighted error
        real(wp) :: err

        err = sum(weights * abs(P_true - P_est))
    end function awerr

end module mod_metrics
