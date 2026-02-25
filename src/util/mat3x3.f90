MODULE mat3x3_util
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  PUBLIC::inv3x3
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
END MODULE mat3x3_util
