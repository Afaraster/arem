!> Plain derived types (Level 0 — passive data containers) for the AREM project.
!> No type-bound procedures; all behavior lives in separate modules per the
!> Array + Modular paradigm.
module mod_types
    use mod_precision, only: wp
    implicit none
    private

    public :: point_2d
    public :: grid_spec
    public :: emitter_set
    public :: sensor_set
    public :: arem_config
    public :: arem_result
    public :: rss_map

    !> A 2D point with x and y coordinates [m].
    type :: point_2d
        real(wp) :: x = 0.0_wp
        real(wp) :: y = 0.0_wp
    end type point_2d

    !> Regular grid specification over a rectangular ROI.
    type :: grid_spec
        !> ROI lower-left x [m]
        real(wp) :: x_min = 0.0_wp
        !> ROI upper-right x [m]
        real(wp) :: x_max = 400.0_wp
        !> ROI lower-left y [m]
        real(wp) :: y_min = 0.0_wp
        !> ROI upper-right y [m]
        real(wp) :: y_max = 400.0_wp
        !> Grid spacing [m]
        real(wp) :: delta = 20.0_wp
        !> Total number of grid points
        integer  :: N = 0
        !> Number of grid points in x-direction
        integer  :: nx = 0
        !> Number of grid points in y-direction
        integer  :: ny = 0
        !> Grid point coordinates [N, 2]: column 1 = x, column 2 = y
        real(wp), allocatable :: points(:,:)
    end type grid_spec

    !> Set of emitter (transmitter) parameters.
    type :: emitter_set
        !> Number of emitters
        integer  :: K = 0
        !> Emitter positions [K, 2]: column 1 = x, column 2 = y
        real(wp), allocatable :: positions(:,:)
        !> Emitter transmit powers [K] in dBm
        real(wp), allocatable :: powers(:)
    end type emitter_set

    !> Set of sensor (receiver) positions.
    type :: sensor_set
        !> Number of sensors
        integer  :: M = 0
        !> Sensor positions [M, 2]: column 1 = x, column 2 = y
        real(wp), allocatable :: positions(:,:)
    end type sensor_set

    !> Configuration parameters for the AREM algorithm.
    type :: arem_config
        !> Convergence threshold for residual change in inner alternating loop
        real(wp) :: th_error = 1e-6_wp
        !> Gradient descent step size $\alpha$
        real(wp) :: alpha = 1.0_wp
        !> Maximum number of outer iterations (new emitter additions)
        integer  :: max_outer_iter = 20
        !> Maximum number of inner iterations (alternating updates per outer iter)
        integer  :: max_inner_iter = 100
        !> Pruning threshold: remove $\omega_i$ with $|\omega_i| <$ prune_threshold
        real(wp) :: prune_threshold = 1e-6_wp
        !> Whether to double grid density in rough estimation step
        logical  :: use_double_grid = .false.
        !> Path-loss exponent $\eta$
        real(wp) :: eta = 2.0_wp
        !> Antenna far-field reference distance $d_0$ [m]
        real(wp) :: d0 = 5.0_wp
    end type arem_config

    !> Results produced by one run of the AREM algorithm.
    type :: arem_result
        !> Estimated emitter positions [n_found, 2]
        real(wp), allocatable :: positions(:,:)
        !> Estimated emitter powers [n_found]
        real(wp), allocatable :: powers(:)
        !> Number of emitters found
        integer  :: n_found = 0
        !> Number of outer iterations executed
        integer  :: n_outer_iter = 0
        !> Whether the algorithm converged
        logical  :: converged = .false.
    end type arem_result

    !> RSS (Radio Signal Strength) map over the grid.
    type :: rss_map
        !> True or estimated RSS values at each grid point [N] in dBm
        real(wp), allocatable :: values(:)
        !> AWErr weights $\lambda_i$ [N], normalized to sum to 1
        real(wp), allocatable :: weights(:)
    end type rss_map

end module mod_types
