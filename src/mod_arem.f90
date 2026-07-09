!> AREM (Accurate Radio Environment Map) core algorithm.
!>
!> Implements Algorithm 1 from the paper: an alternating-minimization method
!> for solving the non-convex bilinear inverse problem of off-grid emitter
!> localization and power estimation from RSS measurements.
!>
!> The algorithm alternates between:
!> 1. Rough estimation via decorrelated dictionary correlation
!> 2. Least-squares power (omega) update
!> 3. Pruning to maintain sparsity
!> 4. Gradient-descent dictionary (position) refinement
!> 5. Residual update
module mod_arem
    use mod_precision, only: wp
    use mod_types, only: arem_config, arem_result, grid_spec
    use mod_pathloss, only: build_path_loss_vector, path_loss
    use mod_preprocess, only: compute_preprocessing
    use mod_gradient, only: compute_gradient
    use stdlib_linalg, only: lstsq, eye
    implicit none
    private

    public :: arem_solve

contains

    !> Run the AREM algorithm to estimate emitter positions and powers.
    !>
    !> @note `A_grid` is the original dictionary. `Aprime` is the preprocessed
    !> dictionary used for rough estimation. Both are computed once upfront
    !> (in the calling driver) and passed in. When `config%use_double_grid`
    !> is true, `Aprime` should be built from a finer grid.
    subroutine arem_solve(z, A_grid, T, Aprime, sensors, grid, config, result)
        !> RSS measurement vector $\mathbf{z}$ [M]
        real(wp), intent(in) :: z(:)
        !> Original path-loss dictionary $\mathbf{A}(\boldsymbol{\theta})$ [M, N]
        real(wp), intent(in) :: A_grid(:, :)
        !> Preprocessing matrix $\mathbf{T}$ [r, M] (from compute_preprocessing)
        real(wp), intent(in) :: T(:, :)
        !> Preprocessed dictionary $\mathbf{A}'$ [r, N]
        real(wp), intent(in) :: Aprime(:, :)
        !> Sensor positions [M, 2]
        real(wp), intent(in) :: sensors(:, :)
        !> Grid specification for the original grid
        type(grid_spec), intent(in) :: grid
        !> Algorithm configuration
        type(arem_config), intent(in) :: config
        !> Estimation results (output)
        type(arem_result), intent(out) :: result

        integer  :: M, N, n_support, k_outer, k_inner
        integer  :: i, idx_max
        real(wp) :: eps_prev, eps_curr
        real(wp) :: alpha

        ! Working arrays
        real(wp), allocatable :: gamma(:), gamma_prime(:)
        real(wp), allocatable :: support_pos(:, :)
        real(wp), allocatable :: omega(:)
        real(wp), allocatable :: A_Tk(:, :)
        real(wp), allocatable :: grad(:, :)
        real(wp), allocatable :: corr(:)
        integer,  allocatable :: support_grid_idx(:)
        logical,  allocatable :: keep_mask(:)

        M = size(z)
        N = size(A_grid, 2)

        ! --- Initialization ---
        gamma = z
        eps_prev = norm2(gamma)
        k_outer = 0
        allocate (support_pos(0, 2))
        allocate (omega(0))
        allocate (support_grid_idx(0))
        alpha = config%alpha

        ! --- Outer loop: add one emitter per iteration ---
        do while (k_outer < config%max_outer_iter)

            ! Data preprocessing: gamma' = T @ gamma
            gamma_prime = matmul(T, gamma)

            ! Rough estimation: find grid point with max correlation
            call rough_estimation(Aprime, gamma_prime, idx_max)

            ! Check if this grid point is already in the support set
            if (allocated(support_grid_idx)) then
                if (any(support_grid_idx == idx_max)) then
                    exit  ! No new emitter found
                end if
            end if

            ! Update support set
            n_support = size(support_pos, 1) + 1
            call expand_support(support_pos, support_grid_idx, &
                                grid%points(idx_max, 1), grid%points(idx_max, 2), idx_max)

            ! Build A_Tk for current support
            A_Tk = build_A_Tk(sensors, support_pos, config%d0, config%eta)

            ! --- Inner alternating update loop ---
            k_inner = 0
            do while (k_inner < config%max_inner_iter)
                k_inner = k_inner + 1

                ! Step 1: Update omega (least squares)
                omega = update_omega(A_Tk, z)

                ! Step 2: Pruning
                call prune_support(support_pos, support_grid_idx, omega, &
                                   A_Tk, grid, config%prune_threshold)
                n_support = size(support_pos, 1)

                if (n_support == 0) exit

                ! Rebuild A_Tk after pruning (dimensions may have changed)
                A_Tk = build_A_Tk(sensors, support_pos, config%d0, config%eta)

                ! Step 3: Gradient descent update of emitter positions
                allocate (grad(n_support, 2))
                call compute_gradient(z, A_Tk, sensors, support_pos, &
                                      config%eta, config%d0, grad)

                ! Gradient step: theta ← theta - alpha * grad_H^T
                ! Note: grad has shape [n_support, 2], where grad(i,1)=∂H/∂x,
                ! grad(i,2)=∂H/∂y. We step in the negative gradient direction.
                do i = 1, n_support
                    support_pos(i, 1) = support_pos(i, 1) &
                        - alpha * grad(i, 1)
                    support_pos(i, 2) = support_pos(i, 2) &
                        - alpha * grad(i, 2)
                end do
                deallocate (grad)

                ! Clamp positions to ROI bounds
                call clamp_to_roi(support_pos, grid)

                ! Rebuild A_Tk with updated positions
                A_Tk = build_A_Tk(sensors, support_pos, config%d0, config%eta)

                ! Step 4: Update residual and check convergence
                omega = update_omega(A_Tk, z)
                gamma = z - matmul(A_Tk, omega)
                eps_curr = norm2(gamma)

                if (abs(eps_curr - eps_prev) < config%th_error) exit
                eps_prev = eps_curr
            end do

            k_outer = k_outer + 1

            ! Check outer termination: residual small enough
            if (eps_curr < config%th_error .or. n_support == 0) exit
        end do

        ! --- Pack results ---
        result%n_found = size(support_pos, 1)
        result%n_outer_iter = k_outer
        result%converged = (k_outer < config%max_outer_iter)

        if (result%n_found > 0) then
            allocate (result%positions(result%n_found, 2))
            allocate (result%powers(result%n_found))
            result%positions = support_pos
            result%powers = omega
        end if
    end subroutine arem_solve

    !> Rough estimation: find grid index with maximum absolute correlation
    !> between the preprocessed dictionary columns and the preprocessed residual.
    !>
    !> Both `Aprime` and `gamma_prime` have been transformed by T:
    !>   Aprime = T @ A,  gamma_prime = T @ gamma
    subroutine rough_estimation(Aprime, gamma_prime, idx_max)
        !> Preprocessed dictionary $\mathbf{A}' = \mathbf{T}\mathbf{A}$ [r, N]
        real(wp), intent(in) :: Aprime(:, :)
        !> Preprocessed residual $\boldsymbol{\gamma}' = \mathbf{T}\boldsymbol{\gamma}$ [r]
        real(wp), intent(in) :: gamma_prime(:)
        !> Index of grid point with maximum correlation
        integer, intent(out) :: idx_max

        real(wp), allocatable :: corr(:)

        ! corr = A'^T @ gamma'  [N] (all inner products at once)
        corr = matmul(transpose(Aprime), gamma_prime)
        idx_max = maxloc(abs(corr), dim=1)
    end subroutine rough_estimation

    !> Update omega via least squares (Eq. 16):
    !> $$\boldsymbol{\omega}_k = (\mathbf{A}_{\mathcal{T}_k}^T
    !>   \mathbf{A}_{\mathcal{T}_k})^{-1} \mathbf{A}_{\mathcal{T}_k}^T
    !>   \mathbf{z}$$
    function update_omega(A_Tk, z) result(omega)
        !> Support-restricted dictionary [M, |T_k|]
        real(wp), intent(in) :: A_Tk(:, :)
        !> Measurement vector [M]
        real(wp), intent(in) :: z(:)
        !> Estimated powers $\boldsymbol{\omega}$ [|T_k|]
        real(wp), allocatable :: omega(:)

        real(wp), allocatable :: A_work(:, :)

        ! Pass a writable copy: lstsq requires intent(inout) for A
        allocate (A_work, source=A_Tk)
        omega = lstsq(A_work, z)
    end function update_omega

    !> Prune support set: remove emitters with small power estimates
    !> or positions outside the ROI (Eq. 17).
    subroutine prune_support(support_pos, support_idx, omega, A_Tk, grid, threshold)
        !> Support positions [n, 2] (inout — may shrink)
        real(wp), allocatable, intent(inout) :: support_pos(:, :)
        !> Grid indices [n] (inout)
        integer,  allocatable, intent(inout) :: support_idx(:)
        !> Power estimates [n] (inout)
        real(wp), allocatable, intent(inout) :: omega(:)
        !> Support dictionary [M, n] (inout — may shrink)
        real(wp), allocatable, intent(inout) :: A_Tk(:, :)
        !> Grid specification (for ROI bounds)
        type(grid_spec), intent(in) :: grid
        !> Pruning threshold for |omega|
        real(wp), intent(in) :: threshold

        integer  :: i, n, n_keep
        logical, allocatable :: keep(:)

        n = size(omega)
        allocate (keep(n))
        keep = .true.

        do i = 1, n
            ! Remove if power too small
            if (abs(omega(i)) < threshold) keep(i) = .false.

            ! Remove if outside ROI
            if (support_pos(i, 1) < grid%x_min .or. &
                support_pos(i, 1) > grid%x_max .or. &
                support_pos(i, 2) < grid%y_min .or. &
                support_pos(i, 2) > grid%y_max) keep(i) = .false.
        end do

        n_keep = count(keep)
        if (n_keep == n) return  ! Nothing to prune

        if (n_keep == 0) then
            deallocate (support_pos, support_idx, omega, A_Tk)
            allocate (support_pos(0, 2))
            allocate (support_idx(0))
            allocate (omega(0))
            allocate (A_Tk(0, 0))
            return
        end if

        ! Compact arrays (keep only non-pruned entries)
        call compact_support(support_pos, support_idx, omega, A_Tk, keep, n_keep)
    end subroutine prune_support

    !> Add a new emitter position to the support set.
    subroutine expand_support(support_pos, support_idx, x, y, idx)
        real(wp), allocatable, intent(inout) :: support_pos(:, :)
        integer,  allocatable, intent(inout) :: support_idx(:)
        real(wp), intent(in) :: x, y
        integer,  intent(in) :: idx

        real(wp), allocatable :: tmp_pos(:, :)
        integer,  allocatable :: tmp_idx(:)
        integer :: n_old

        if (.not. allocated(support_pos)) then
            allocate (support_pos(1, 2))
            allocate (support_idx(1))
            support_pos(1, 1) = x
            support_pos(1, 2) = y
            support_idx(1) = idx
            return
        end if

        n_old = size(support_pos, 1)
        allocate (tmp_pos(n_old, 2))
        allocate (tmp_idx(n_old))
        tmp_pos = support_pos
        tmp_idx = support_idx

        deallocate (support_pos, support_idx)
        allocate (support_pos(n_old + 1, 2))
        allocate (support_idx(n_old + 1))
        support_pos(1:n_old, :) = tmp_pos
        support_pos(n_old + 1, 1) = x
        support_pos(n_old + 1, 2) = y
        support_idx(1:n_old) = tmp_idx
        support_idx(n_old + 1) = idx
    end subroutine expand_support

    !> Compact support arrays by keeping only entries where `keep(i)` is true.
    subroutine compact_support(support_pos, support_idx, omega, A_Tk, keep, n_keep)
        real(wp), allocatable, intent(inout) :: support_pos(:, :)
        integer,  allocatable, intent(inout) :: support_idx(:)
        real(wp), allocatable, intent(inout) :: omega(:)
        real(wp), allocatable, intent(inout) :: A_Tk(:, :)
        logical,  intent(in) :: keep(:)
        integer,  intent(in) :: n_keep

        integer :: i, j, n_old, M_sensors

        n_old = size(support_pos, 1)
        M_sensors = size(A_Tk, 1)

        ! Compact positions, indices, omega
        call compact_1d_pos(support_pos, keep, n_keep, 2)
        call compact_1d_int(support_idx, keep, n_keep)
        call compact_1d_real(omega, keep, n_keep)

        ! Compact A_Tk columns
        block
            real(wp), allocatable :: tmp(:, :)
            integer :: col

            allocate (tmp(M_sensors, n_keep))
            col = 0
            do i = 1, n_old
                if (keep(i)) then
                    col = col + 1
                    tmp(:, col) = A_Tk(:, i)
                end if
            end do
            deallocate (A_Tk)
            allocate (A_Tk(M_sensors, n_keep))
            A_Tk = tmp
        end block
    end subroutine compact_support

    !> Helper: compact a 2D real array along the first dimension using a mask.
    subroutine compact_1d_pos(arr, keep, n_keep, ncol)
        real(wp), allocatable, intent(inout) :: arr(:, :)
        logical,  intent(in) :: keep(:)
        integer,  intent(in) :: n_keep, ncol

        real(wp), allocatable :: tmp(:, :)
        integer :: i, j, row

        allocate (tmp(n_keep, ncol))
        row = 0
        do i = 1, size(arr, 1)
            if (keep(i)) then
                row = row + 1
                tmp(row, :) = arr(i, :)
            end if
        end do
        deallocate (arr)
        allocate (arr(n_keep, ncol))
        arr = tmp
    end subroutine compact_1d_pos

    subroutine compact_1d_int(arr, keep, n_keep)
        integer, allocatable, intent(inout) :: arr(:)
        logical, intent(in) :: keep(:)
        integer, intent(in) :: n_keep

        integer, allocatable :: tmp(:)
        integer :: i, j

        allocate (tmp(n_keep))
        j = 0
        do i = 1, size(arr)
            if (keep(i)) then
                j = j + 1
                tmp(j) = arr(i)
            end if
        end do
        deallocate (arr)
        allocate (arr(n_keep))
        arr = tmp
    end subroutine compact_1d_int

    subroutine compact_1d_real(arr, keep, n_keep)
        real(wp), allocatable, intent(inout) :: arr(:)
        logical, intent(in) :: keep(:)
        integer, intent(in) :: n_keep

        real(wp), allocatable :: tmp(:)
        integer :: i, j

        allocate (tmp(n_keep))
        j = 0
        do i = 1, size(arr)
            if (keep(i)) then
                j = j + 1
                tmp(j) = arr(i)
            end if
        end do
        deallocate (arr)
        allocate (arr(n_keep))
        arr = tmp
    end subroutine compact_1d_real

    !> Build the support-restricted dictionary $\mathbf{A}_{\mathcal{T}_k}$
    !> from current support positions.
    function build_A_Tk(sensors, support_pos, d0, eta) result(A_Tk)
        !> Sensor positions [M, 2]
        real(wp), intent(in) :: sensors(:, :)
        !> Support emitter positions [|T_k|, 2]
        real(wp), intent(in) :: support_pos(:, :)
        !> Reference distance $d_0$
        real(wp), intent(in) :: d0
        !> Path-loss exponent $\eta$
        real(wp), intent(in) :: eta
        !> Support dictionary [M, |T_k|]
        real(wp), allocatable :: A_Tk(:, :)

        integer :: i, m, M_sensors, n_support
        real(wp) :: dx, dy, dist

        M_sensors = size(sensors, 1)
        n_support = size(support_pos, 1)

        if (n_support == 0) then
            allocate (A_Tk(0, 0))
            return
        end if

        allocate (A_Tk(M_sensors, n_support))

        do i = 1, n_support
            do m = 1, M_sensors
                dx = sensors(m, 1) - support_pos(i, 1)
                dy = sensors(m, 2) - support_pos(i, 2)
                dist = sqrt(dx**2 + dy**2)
                A_Tk(m, i) = path_loss(dist, d0, eta)
            end do
        end do
    end function build_A_Tk

    !> Clamp emitter positions to stay within ROI bounds.
    subroutine clamp_to_roi(positions, grid)
        !> Emitter positions [n, 2] (modified in-place)
        real(wp), intent(inout) :: positions(:, :)
        !> Grid specification (provides ROI bounds)
        type(grid_spec), intent(in) :: grid

        integer :: i

        do i = 1, size(positions, 1)
            positions(i, 1) = max(grid%x_min, min(grid%x_max, positions(i, 1)))
            positions(i, 2) = max(grid%y_min, min(grid%y_max, positions(i, 2)))
        end do
    end subroutine clamp_to_roi

end module mod_arem
