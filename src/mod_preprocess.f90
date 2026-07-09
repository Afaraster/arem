!> Decorrelation preprocessing for the AREM rough estimation step.
!>
!> Implements Eqs. (13-14): computes the preprocessing matrix T and the
!> transformed dictionary A' = T @ A to reduce column correlations.
!>
!> Uses thin SVD: A = U * diag(S) * V^T, then:
!>   T = diag(1/S) * U^T  [r x M]  (pseudoinverse of U*Sigma, i.e. whitening)
!>   A' = T * A = V^T      [r x N]  (perfectly decorrelated columns)
module mod_preprocess
    use mod_precision, only: wp
    use stdlib_linalg, only: svd
    implicit none
    private

    public :: compute_preprocessing
    public :: compute_pinv_via_svd

contains

    !> Compute the preprocessing matrix T and transformed dictionary A'.
    !>
    !> Uses thin SVD to build T = pinv(U*Sigma) = diag(1/s_i) * U^T.
    !> This yields A' = V^T, whose columns are orthogonal — the ideal
    !> decorrelation shown in Figure 2(b) of the paper.
    subroutine compute_preprocessing(A, T, Aprime)
        !> Original path-loss dictionary [M, N]
        real(wp), intent(in)  :: A(:, :)
        !> Preprocessing matrix [r, M] where r = min(M, N)
        real(wp), allocatable, intent(out) :: T(:, :)
        !> Transformed dictionary [r, N]
        real(wp), allocatable, intent(out) :: Aprime(:, :)

        integer  :: i, M, N, r
        real(wp) :: smax, tol
        real(wp), allocatable :: S(:), U(:, :), Vt(:, :), A_copy(:, :)

        M = size(A, 1)
        N = size(A, 2)
        r = min(M, N)

        ! Make a copy for SVD (stdlib may overwrite)
        allocate (A_copy(M, N))
        A_copy = A

        ! Thin SVD: A = U * diag(S) * V^T
        allocate (S(r), U(M, r), Vt(r, N))
        call svd(A_copy, S, U=U, Vt=Vt, full_matrices=.false.)

        ! Build T = diag(1/s_i) * U^T  [r x M]
        ! Only include singular values above tolerance
        smax = maxval(S)
        tol = max(M, N) * smax * epsilon(1.0_wp)

        allocate (T(r, M))
        T = 0.0_wp
        do i = 1, r
            if (S(i) > tol) then
                T(i, :) = U(:, i) / S(i)
            end if
        end do

        ! A' = T * A = V^T (up to tolerance)  [r x N]
        Aprime = matmul(T, A)
    end subroutine compute_preprocessing

    !> Build the Moore-Penrose pseudoinverse from thin SVD components.
    !>
    !> Given A = U * diag(S) * V^T, returns
    !> A† = V * diag(1/s_i) * U^T where s_i below tol are treated as zero.
    function compute_pinv_via_svd(S, U, Vt, tol) result(A_pinv)
        !> Singular values s_i [r]
        real(wp), intent(in) :: S(:)
        !> Left singular vectors [M, r] (thin)
        real(wp), intent(in) :: U(:, :)
        !> Right singular vectors transposed [r, N] (thin)
        real(wp), intent(in) :: Vt(:, :)
        !> Tolerance for singular value cutoff
        real(wp), intent(in), optional :: tol
        !> Pseudoinverse A† [N, M]
        real(wp), allocatable :: A_pinv(:, :)

        integer  :: i, r, M, N
        real(wp) :: tol_actual, smax
        real(wp), allocatable :: SinvUt(:, :)

        r = size(S)
        M = size(U, 1)
        N = size(Vt, 2)

        if (present(tol)) then
            tol_actual = tol
        else
            smax = maxval(S)
            tol_actual = max(M, N) * smax * epsilon(1.0_wp)
        end if

        ! Compute diag(1/S) * U^T  [r, M]
        allocate (SinvUt(r, M))
        SinvUt = 0.0_wp
        do i = 1, r
            if (S(i) > tol_actual) then
                SinvUt(i, :) = U(:, i) / S(i)
            end if
        end do

        ! A_pinv = V * SinvUt = Vt^T * SinvUt  [N, r] * [r, M] = [N, M]
        A_pinv = matmul(transpose(Vt), SinvUt)
    end function compute_pinv_via_svd

end module mod_preprocess
