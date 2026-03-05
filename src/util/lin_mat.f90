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
  SUBROUTINE mat_mul(ndim, C, A, Ac, B, Bc)
    !< C_mn = \sum_k (A^Ac)_mk * (B^Bc)_kn \
    !< (ndim, ndim) Square matrices
    USE kinds, ONLY: DP
    INTEGER, INTENT(IN) :: ndim

    COMPLEX(DP), INTENT(IN)  :: A(ndim, ndim)
    !< : left matrix
    CHARACTER(len=1), INTENT(IN) :: Ac
    !< : [left matrix]
    !< N : normal / T : transpose / C : complex conjugate
    COMPLEX(DP), INTENT(IN)  :: B(ndim, ndim)
    !< : right matrix
    CHARACTER(len=1), INTENT(IN) :: Bc
    !< : [right matrix]
    !< N : normal / T : transpose / C : complex conjugate
    COMPLEX(DP), INTENT(OUT) :: C(ndim, ndim)
    !< output matrix
    COMPLEX(DP), PARAMETER::ALPHA = (1.0_DP, 0.0_DP)
    COMPLEX(DP), PARAMETER::BETA = (0.0_DP, 0.0_DP)

    CALL ZGEMM(Ac, Bc, ndim, ndim, ndim, ALPHA, A, ndim, B, ndim, BETA, C, ndim)
  END SUBROUTINE mat_mul
END MODULE lin_mat
