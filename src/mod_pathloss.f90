!> Path-loss model and dictionary construction for the AREM algorithm.
!>
!> Implements the piecewise path-loss function (Eq. 4), its partial derivatives
!> (Eq. 21), and utilities to build path-loss vectors / dictionaries.
!>
!> Core functions are `elemental pure` for array broadcasting and compiler
!> optimization.
module mod_pathloss
    use mod_precision, only: wp
    use mod_types, only: grid_spec
    implicit none
    private

    public :: path_loss
    public :: path_loss_deriv_x, path_loss_deriv_y
    public :: build_path_loss_vector
    public :: build_path_loss_dictionary

contains

    !> Piecewise path-loss function (Eq. 4):
    !> $$f(d) = 1 \text{ for } d \leq d_0, \quad f(d) = (d_0/d)^\eta \text{ for } d > d_0$$
    elemental pure real(wp) function path_loss(d, d0, eta) result(f)
        !> Euclidean distance between sensor and emitter [m]
        real(wp), intent(in) :: d
        !> Antenna far-field reference distance [m]
        real(wp), intent(in) :: d0
        !> Path-loss exponent
        real(wp), intent(in) :: eta

        if (d <= d0) then
            f = 1.0_wp
        else
            f = (d0 / d)**eta
        end if
    end function path_loss

    !> Partial derivative of the path-loss function w.r.t. the emitter
    !> x-coordinate (Eq. 21):
    !> $$\frac{\partial f}{\partial x_i^g} = \eta d_0^\eta \frac{x_m^s - x_i^g}{\|\mathbf{s}_m - \boldsymbol{\theta}_i\|_2^{\eta+2}}$$
    elemental pure real(wp) function path_loss_deriv_x(sx, sy, tx, ty, eta, d0) result(df)
        !> Sensor x-coordinate $x_m^s$ [m]
        real(wp), intent(in) :: sx
        !> Sensor y-coordinate $y_m^s$ [m]
        real(wp), intent(in) :: sy
        !> Emitter/grid x-coordinate $x_i^g$ [m]
        real(wp), intent(in) :: tx
        !> Emitter/grid y-coordinate $y_i^g$ [m]
        real(wp), intent(in) :: ty
        !> Path-loss exponent $\eta$
        real(wp), intent(in) :: eta
        !> Antenna far-field reference distance $d_0$ [m]
        real(wp), intent(in) :: d0

        real(wp) :: dist, dx

        dx = sx - tx
        dist = sqrt(dx**2 + (sy - ty)**2)

        if (dist <= d0) then
            df = 0.0_wp
        else
            df = eta * d0**eta * dx / dist**(eta + 2.0_wp)
        end if
    end function path_loss_deriv_x

    !> Partial derivative of the path-loss function w.r.t. the emitter
    !> y-coordinate: analogous to `path_loss_deriv_x` with $(y_m^s - y_i^g)$.
    elemental pure real(wp) function path_loss_deriv_y(sx, sy, tx, ty, eta, d0) result(df)
        !> Sensor x-coordinate $x_m^s$ [m]
        real(wp), intent(in) :: sx
        !> Sensor y-coordinate $y_m^s$ [m]
        real(wp), intent(in) :: sy
        !> Emitter/grid x-coordinate $x_i^g$ [m]
        real(wp), intent(in) :: tx
        !> Emitter/grid y-coordinate $y_i^g$ [m]
        real(wp), intent(in) :: ty
        !> Path-loss exponent $\eta$
        real(wp), intent(in) :: eta
        !> Antenna far-field reference distance $d_0$ [m]
        real(wp), intent(in) :: d0

        real(wp) :: dist, dy

        dy = sy - ty
        dist = sqrt((sx - tx)**2 + dy**2)

        if (dist <= d0) then
            df = 0.0_wp
        else
            df = eta * d0**eta * dy / dist**(eta + 2.0_wp)
        end if
    end function path_loss_deriv_y

    !> Build the path-loss vector $\mathbf{a}(\mathbf{u}_k) \in \mathbb{R}^M$
    !> for a single emitter position (Eq. 5).
    function build_path_loss_vector(sensors, theta, d0, eta) result(a)
        !> Sensor positions [M, 2]
        real(wp), intent(in) :: sensors(:,:)
        !> Emitter/grid coordinates $[x, y]$
        real(wp), intent(in) :: theta(2)
        !> Antenna far-field reference distance [m]
        real(wp), intent(in) :: d0
        !> Path-loss exponent
        real(wp), intent(in) :: eta
        !> Path-loss vector [M]
        real(wp), allocatable :: a(:)

        integer  :: m, M_sensors
        real(wp) :: dx, dy, dist

        M_sensors = size(sensors, 1)
        allocate (a(M_sensors))

        do m = 1, M_sensors
            dx = sensors(m, 1) - theta(1)
            dy = sensors(m, 2) - theta(2)
            dist = sqrt(dx**2 + dy**2)
            a(m) = path_loss(dist, d0, eta)
        end do
    end function build_path_loss_vector

    !> Build the full path-loss dictionary $\mathbf{A}(\boldsymbol{\theta}) \in
    !> \mathbb{R}^{M \times N}$ for all grid points (Eq. 6).
    function build_path_loss_dictionary(sensors, grid, d0, eta) result(A)
        !> Sensor positions [M, 2]
        real(wp), intent(in) :: sensors(:,:)
        !> Grid specification (must have `points` allocated)
        type(grid_spec), intent(in) :: grid
        !> Antenna far-field reference distance [m]
        real(wp), intent(in) :: d0
        !> Path-loss exponent
        real(wp), intent(in) :: eta
        !> Path-loss dictionary [M, N]
        real(wp), allocatable :: A(:, :)

        integer :: i, m, M_sensors, N_grid
        real(wp) :: dx, dy, dist

        M_sensors = size(sensors, 1)
        N_grid = grid%N
        allocate (A(M_sensors, N_grid))

        do i = 1, N_grid
            do m = 1, M_sensors
                dx = sensors(m, 1) - grid%points(i, 1)
                dy = sensors(m, 2) - grid%points(i, 2)
                dist = sqrt(dx**2 + dy**2)
                A(m, i) = path_loss(dist, d0, eta)
            end do
        end do
    end function build_path_loss_dictionary

end module mod_pathloss
