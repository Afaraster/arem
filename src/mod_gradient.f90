!> Gradient computation for the AREM dictionary refinement step.
!>
!> Computes the gradient of the objective
!> $H(\boldsymbol{\theta}) = \|\mathbf{z} - \mathbf{A}_{\mathcal{T}_k}
!> \boldsymbol{\omega}_k\|_2^2$
!> with respect to emitter positions using the Envelope Theorem:
!> $$\nabla_{\boldsymbol{\theta}_i} H =
!>   -2 \omega_i \mathbf{J}_{\mathbf{a}_i}^T \boldsymbol{\gamma}$$
!>
!> where $\boldsymbol{\gamma} = \mathbf{z} - \mathbf{A}_{\mathcal{T}_k}
!> \boldsymbol{\omega}_k$ is the residual and $\mathbf{J}_{\mathbf{a}_i}$
!> is the Jacobian of the path-loss vector with respect to $(x_i, y_i)$.
!>
!> This formulation avoids the numerically unstable pseudo-inverse
!> manipulations in the paper's Eqs. (19-21).
module mod_gradient
    use mod_precision, only: wp
    use mod_pathloss, only: path_loss_deriv_x, path_loss_deriv_y
    use stdlib_linalg, only: lstsq
    implicit none
    private

    public :: compute_H
    public :: compute_gradient
    public :: compute_A_pinv_z

contains

    !> Compute the objective function $H(\boldsymbol{\theta})$ (Eq. 18):
    !> $$H = \|\mathbf{z} - \mathbf{A}_{\mathcal{T}_k}
    !>   \mathbf{A}_{\mathcal{T}_k}^{\dagger} \mathbf{z}\|_2^2$$
    function compute_H(z, A_Tk) result(H)
        !> Measurement vector $\mathbf{z}$ [M]
        real(wp), intent(in) :: z(:)
        !> Support-restricted dictionary $\mathbf{A}_{\mathcal{T}_k}$ [M, |T_k|]
        real(wp), intent(in) :: A_Tk(:, :)
        !> Objective value
        real(wp) :: H

        real(wp), allocatable :: w(:), A_work(:, :)
        real(wp) :: ztz

        ztz = dot_product(z, z)
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

        allocate (A_work, source=A_Tk)
        w = lstsq(A_work, z)
    end function compute_A_pinv_z

    !> Compute the gradient of $H$ w.r.t. each emitter position via the
    !> Envelope Theorem.
    !>
    !> For emitter $i$:
    !> $$\nabla_{\boldsymbol{\theta}_i} H =
    !>   -2 \omega_i \mathbf{J}_{\mathbf{a}_i}^T \boldsymbol{\gamma}$$
    !>
    !> where $\boldsymbol{\gamma} = \mathbf{z} - \mathbf{A} \boldsymbol{\omega}$
    !> and $\mathbf{J}_{\mathbf{a}_i} = [\partial\mathbf{a}/\partial x_i,\;
    !> \partial\mathbf{a}/\partial y_i] \in \mathbb{R}^{M \times 2}$.
    subroutine compute_gradient(z, A_Tk, omega, sensors, support, eta, d0, grad)
        !> Measurement vector $\mathbf{z}$ [M]
        real(wp), intent(in) :: z(:)
        !> Support-restricted dictionary $\mathbf{A}_{\mathcal{T}_k}$ [M, |T_k|]
        real(wp), intent(in) :: A_Tk(:, :)
        !> Current power estimates $\boldsymbol{\omega}$ [|T_k|]
        real(wp), intent(in) :: omega(:)
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
        real(wp), allocatable :: gamma(:)
        real(wp) :: df_dx, df_dy
        real(wp) :: dot_x, dot_y

        M_sensors = size(sensors, 1)
        n_support = size(support, 1)

        ! Residual: gamma = z - A_Tk * omega
        gamma = z - matmul(A_Tk, omega)

        do i = 1, n_support
            ! Compute dot products: J_{a_i}^T * gamma
            dot_x = 0.0_wp
            dot_y = 0.0_wp
            do m = 1, M_sensors
                df_dx = path_loss_deriv_x( &
                    sensors(m, 1), sensors(m, 2), &
                    support(i, 1), support(i, 2), eta, d0)
                df_dy = path_loss_deriv_y( &
                    sensors(m, 1), sensors(m, 2), &
                    support(i, 1), support(i, 2), eta, d0)
                dot_x = dot_x + df_dx * gamma(m)
                dot_y = dot_y + df_dy * gamma(m)
            end do

            ! Envelope Theorem: grad_i = -2 * omega_i * J_i^T * gamma
            grad(i, 1) = -2.0_wp * omega(i) * dot_x
            grad(i, 2) = -2.0_wp * omega(i) * dot_y
        end do
    end subroutine compute_gradient

end module mod_gradient
