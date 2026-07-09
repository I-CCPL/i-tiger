MODULE lin_mat
  USE kinds, ONLY: DP
  IMPLICIT NONE
CONTAINS
  FUNCTION inv3x3(mat) RESULT(mat_inv)
    REAL(DP), INTENT(IN) :: mat(3, 3)
    INTEGER::idx, jdx, kdx, ldx, mdx, ndx
    REAL(DP) :: mat_inv(3, 3), det

    det = 0
    DO idx = 1, 3
      jdx = MOD(idx, 3) + 1
      kdx = MOD(idx + 1, 3) + 1
      det = det + mat(1, idx)*(mat(2, jdx)*mat(3, kdx) - mat(2, kdx)*mat(3, jdx))
    END DO

    IF (ABS(det) <= 1.0D-14) CALL errore(1, 'inv3x3', 'Singular 3x3 matrix.')

    ! Inverse = transpose(cofactor)/determinant
    DO idx = 1, 3
      kdx = MOD(idx, 3) + 1
      ldx = MOD(idx + 1, 3) + 1
      DO jdx = 1, 3
        mdx = MOD(jdx, 3) + 1
        ndx = MOD(jdx + 1, 3) + 1
        mat_inv(jdx, idx) = (mat(kdx, mdx)*mat(ldx, ndx) - mat(kdx, ndx)*mat(ldx, mdx))/det
      END DO
    END DO
  END FUNCTION inv3x3
  !
  SUBROUTINE mat_mul(A, Ac, B, Bc, C)
    !< C = op(A) * op(B)
    !< op(X): N = X, T = X^T, C = X^dagger
    !< (ndim, ndim) Square matrices
    USE kinds, ONLY: DP
    USE constants, ONLY: cmplx_0, cmplx_1
    COMPLEX(DP), INTENT(IN)  :: A(:, :)
    !< : left matrix
    CHARACTER(len=1), INTENT(IN) :: Ac
    !< : [left matrix]
    !< N : normal / T : transpose / C : conjugate transpose
    COMPLEX(DP), INTENT(IN)  :: B(:, :)
    !< : right matrix
    CHARACTER(len=1), INTENT(IN) :: Bc
    !< : [right matrix]
    !< N : normal / T : transpose / C : conjugate transpose
    COMPLEX(DP), INTENT(OUT) :: C(:, :)
    !< output matrix
    INTEGER::m, n, k

    m = SIZE(C, 1)
    n = SIZE(C, 2)
    IF (Ac /= 'Z') THEN
      k = SIZE(A, 1)
    ELSE
      k = SIZE(A, 2)
    END IF

    CALL ZGEMM(Ac, Bc, m, n, k, cmplx_1, A, SIZE(A, 1), B, SIZE(B, 1), cmplx_0, C, m)
  END SUBROUTINE mat_mul
END MODULE lin_mat
