!> File I/O utilities for the AREM project.
!>
!> Writes numerical results to `outputs/` in whitespace-delimited or CSV format
!> for post-processing by Python scripts.
module mod_io
    use mod_precision, only: wp
    implicit none
    private

    public :: write_array_2d, write_vector
    public :: write_sweep_csv
    public :: ensure_output_dir

    !> Output directory path (relative to project root)
    character(len=*), parameter, public :: OUTPUT_DIR = "outputs"

contains

    !> Write a 2D array as whitespace-delimited text.
    subroutine write_array_2d(filename, data)
        !> File name (appended to output directory)
        character(len=*), intent(in) :: filename
        !> Data array [rows, cols]
        real(wp), intent(in) :: data(:, :)

        integer :: i, unit

        open (newunit=unit, file=OUTPUT_DIR//"/"//filename, &
              status="replace", action="write")
        do i = 1, size(data, 1)
            write (unit, "(*(ES15.7,1X))") data(i, :)
        end do
        close (unit)
    end subroutine write_array_2d

    !> Write a 1D array as a single-column text file.
    subroutine write_vector(filename, data)
        !> File name (appended to output directory)
        character(len=*), intent(in) :: filename
        !> Data array [n]
        real(wp), intent(in) :: data(:)

        integer :: i, unit

        open (newunit=unit, file=OUTPUT_DIR//"/"//filename, &
              status="replace", action="write")
        do i = 1, size(data)
            write (unit, "(ES15.7)") data(i)
        end do
        close (unit)
    end subroutine write_vector

    !> Write sweep results as a CSV file with a header row.
    !> Supports multiple y-columns (e.g., three SNR curves) via an optional
    !> labels array.
    subroutine write_sweep_csv(filename, header, x_values, y_values, labels)
        !> File name (appended to output directory)
        character(len=*), intent(in) :: filename
        !> Header string (e.g., "M,AWErr")
        character(len=*), intent(in) :: header
        !> x-axis sweep values [n_points]
        real(wp), intent(in) :: x_values(:)
        !> y-axis values [n_points, n_curves]
        real(wp), intent(in) :: y_values(:, :)
        !> Optional labels for each y-column
        character(len=*), intent(in), optional :: labels(:)

        integer            :: i, j, unit, n_points, n_curves
        character(len=256) :: fmt_string

        n_points = size(x_values)
        n_curves = size(y_values, 2)

        open (newunit=unit, file=OUTPUT_DIR//"/"//filename, &
              status="replace", action="write")

        ! Write header
        write (unit, "(A)") trim(header)

        ! Write data rows: x, y1, y2, ..., yn
        do i = 1, n_points
            write (fmt_string, "(""(*(ES15.7, '',''))"")")
            write (unit, "(*(ES15.7,','))") x_values(i), y_values(i, :)
        end do

        close (unit)
    end subroutine write_sweep_csv

    !> Ensure the output directory exists; create it if necessary.
    subroutine ensure_output_dir()
        logical :: exists

        inquire (file=OUTPUT_DIR, exist=exists)
        if (.not. exists) then
            call execute_command_line("mkdir -p "//OUTPUT_DIR)
        end if
    end subroutine ensure_output_dir

end module mod_io
