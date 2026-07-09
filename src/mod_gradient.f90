!> Gradient computation for the AREM dictionary refinement step.
!>
!> Implements the objective $H(\boldsymbol{\theta})$ (Eq. 18) and its gradient
!> $\nabla H(\boldsymbol{\theta})$ (Eqs. 19–21) used in the alternating update
!> for continuous-space emitter localization.
module mod_gradient
    use mod_precision, only: wp
    use mod_pathloss, only: path_loss_deriv_x, path_loss_deriv_y
    use mod_preprocess, only: compute_pinv_via_svd
    use stdlib_linalg, only: lstsq, svd, solve_lstsq
    implicit none
    private

    public :: compute_H
    public :: compute_gradient
    public :: compute_A_pinv_z

contains

    !> Compute the objective function $H(\boldsymbol{\theta})$ (Eq. 18):
    !> $$H = \mathbf{z}^T\mathbf{z} -
    !>     \mathbf{z}^T \mathbf{A}_{\mathcal{T}_k}
    !>     \mathbf{A}_{\mathcal{T}_k}^{\dagger} \mathbf{z}$$
    function compute_H(z, A_Tk) result(H)
        !> Measurement vector $\mathbf{z}$ [M]
        real(wp), intent(in) :: z(:)
        !> Support-restricted dictionary $\mathbf{A}_{\mathcal{T}_k}$ [M, |T_k|]
        real(wp), intent(in) :: A_Tk(:, :)
        !> Objective value
        real(wp) :: H

        real(wp), allocatable :: w(:)
        real(wp) :: ztz

        real(wp), allocatable :: A_work(:, :)

        ztz = dot_product(z, z)
        ! Pass a writable copy: lstsq requires intent(inout) for A
        allocate (A_work, source=A_Tk)
        w = lstsq(A_work, z)
        H = ztz - dot_product(z, matmul(A_Tk, w))
    end function compute_H

    !> Compute $\mathbf{A}_{\mathcal{T}_k}^{\dagger} \mathbf{z}$ via least squares.
    function compute_A_pinv_z(A_Tk, z) result(w)
        !> Support-restricted dictionary [M, |T_k|]
        real(wp), intent(in) :: A_Tk(:, :)
        !> Measurement vector [M]
        real(wp), intent(in) :: z(:)
        !> $\mathbf{w} = \mathbf{A}_{\mathcal{T}_k}^{\dagger} \mathbf{z}$
        real(wp), allocatable :: w(:)

        real(wp), allocatable :: A_work(:, :)

        ! Pass a writable copy: lstsq requires intent(inout) for A
        allocate (A_work, source=A_Tk)
        w = lstsq(A_work, z)
    end function compute_A_pinv_z

    !> Compute the gradient $\nabla H(\boldsymbol{\theta})$ w.r.t. all emitter
    !> positions in the support set (Eqs. 19–21).
    !>
    !> For each $\boldsymbol{\theta}_i \in \mathcal{T}_k$, computes
    !> $\partial H / \partial x_i^g$ and $\partial H / \partial y_i^g$ via:
    !> $$\mathbf{R}_i = [\mathbf{I} - \mathbf{A}\mathbf{A}^{\dagger}]
    !>   \frac{\partial \mathbf{A}}{\partial x_i^g} \mathbf{A}^{\dagger}$$
    !> $$\frac{\partial H}{\partial x_i^g} =
    !>   \mathbf{z}^T (\mathbf{R}_i + \mathbf{R}_i^T) \mathbf{z}$$
    subroutine compute_gradient(z, A_Tk, sensors, support, eta, d0, grad)
        !> Measurement vector $\mathbf{z}$ [M]
        real(wp), intent(in) :: z(:)
        !> Support-restricted dictionary $\mathbf{A}_{\mathcal{T}_k}$ [M, |T_k|]
        real(wp), intent(in) :: A_Tk(:, :)
        !> Sensor positions [M, 2]
        real(wp), intent(in) :: sensors(:, :)
        !> Support set emitter positions [|T_k|, 2]
        real(wp), intent(in) :: support(:, :)
        !> Path-loss exponent $\eta$
        real(wp), intent(in) :: eta
        !> Antenna far-field reference distance $d_0$ [m]
        real(wp), intent(in) :: d0
        !> Gradient $\nabla H$ [|T_k|, 2]: column 1 = $\partial H/\partial x$,
        !> column 2 = $\partial H/\partial y$
        real(wp), intent(out) :: grad(:, :)

        integer  :: i, m, M_sensors, n_support
        real(wp), allocatable :: I_minus_AApinv(:, :)
        real(wp), allocatable :: A_pinv(:, :)
        real(wp), allocatable :: dA_dxi(:, :), Ri(:, :)
        real(wp), allocatable :: temp(:, :)
        real(wp) :: tmp

        M_sensors = size(sensors, 1)
        n_support = size(support, 1)

        ! A_pinv = A_Tk^† [|T_k|, M]
        A_pinv = compute_pinv_full(A_Tk)

        ! I - A * A^†  [M, M]
        I_minus_AApinv = compute_projection_complement(A_Tk, A_pinv)

        allocate (dA_dxi(M_sensors, n_support))
        allocate (Ri(M_sensors, M_sensors))
        allocate (temp(M_sensors, M_sensors))

        do i = 1, n_support
            ! --- ∂H/∂x_i ---
            dA_dxi = 0.0_wp
            do m = 1, M_sensors
                dA_dxi(m, i) = path_loss_deriv_x( &
                                                  sensors(m, 1), sensors(m, 2), &
                                                  support(i, 1), support(i, 2), eta, d0)
            end do

            ! R_i = [I - A A^†] * (∂A/∂x_i) * A^†
            temp = matmul(dA_dxi, A_pinv)
            Ri = matmul(I_minus_AApinv, temp)

            ! ∂H/∂x_i = z^T (R_i + R_i^T) z
            tmp = dot_product(z, matmul(Ri + transpose(Ri), z))
            grad(i, 1) = tmp

            ! --- ∂H/∂y_i ---
            dA_dxi = 0.0_wp
            do m = 1, M_sensors
                dA_dxi(m, i) = path_loss_deriv_y( &
                                                  sensors(m, 1), sensors(m, 2), &
                                                  support(i, 1), support(i, 2), eta, d0)
            end do

            temp = matmul(dA_dxi, A_pinv)
            Ri = matmul(I_minus_AApinv, temp)

            tmp = dot_product(z, matmul(Ri + transpose(Ri), z))
            grad(i, 2) = tmp
        end do
    end subroutine compute_gradient

    !> Compute the full pseudoinverse $\mathbf{A}^{\dagger}$ [N, M] via thin SVD
    !> and thresholded reciprocal singular values.
    function compute_pinv_full(A) result(A_pinv)
        !> Input matrix [M, N]
        real(wp), intent(in) :: A(:, :)
        !> Pseudoinverse [N, M]
        real(wp), allocatable :: A_pinv(:, :)

        integer  :: M, N, r
        real(wp), allocatable :: S(:), U(:, :), Vt(:, :), A_copy(:, :)

        M = size(A, 1)
        N = size(A, 2)
        r = min(M, N)

        allocate (A_copy(M, N))
        A_copy = A
        allocate (S(r), U(M, r), Vt(r, N))
        call svd(A_copy, S, U=U, Vt=Vt, full_matrices=.false.)

        A_pinv = compute_pinv_via_svd(S, U, Vt)
    end function compute_pinv_full

    !> Compute the projection complement $\mathbf{I} - \mathbf{A}\mathbf{A}^{\dagger}$.
    function compute_projection_complement(A, A_pinv) result(Pc)
        !> Input matrix [M, N]
        real(wp), intent(in) :: A(:, :)
        !> Pseudoinverse [N, M]
        real(wp), intent(in) :: A_pinv(:, :)
        !> $\mathbf{I} - \mathbf{A}\mathbf{A}^{\dagger}$ [M, M]
        real(wp), allocatable :: Pc(:, :)

        integer :: M, i

        M = size(A, 1)
        allocate (Pc(M, M))

        Pc = -matmul(A, A_pinv)
        do i = 1, M
            Pc(i, i) = Pc(i, i) + 1.0_wp
        end do
    end function compute_projection_complement

end module mod_gradient
