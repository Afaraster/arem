module mod_precision
    use iso_fortran_env, only: real64
    implicit none
    private

    public :: wp

    integer, parameter :: wp = real64

end module mod_precision
