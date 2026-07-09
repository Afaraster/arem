!> Physical and numerical constants used throughout the AREM project.
!> All values are in SI units unless otherwise noted.
module mod_constants
    use mod_precision, only: wp
    implicit none
    private

    public :: PI
    public :: DEFAULT_D0, DEFAULT_ETA
    public :: DEFAULT_GRID_DELTA, DEFAULT_ROI_SIZE

    !> $\pi$
    real(wp), parameter :: PI = acos(-1.0_wp)

    !> Default antenna far-field reference distance [m]
    real(wp), parameter :: DEFAULT_D0 = 5.0_wp

    !> Default path-loss exponent
    real(wp), parameter :: DEFAULT_ETA = 2.0_wp

    !> Default grid spacing [m]
    real(wp), parameter :: DEFAULT_GRID_DELTA = 20.0_wp

    !> Default ROI side length [m] (square region)
    real(wp), parameter :: DEFAULT_ROI_SIZE = 400.0_wp

end module mod_constants
