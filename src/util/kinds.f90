MODULE kinds
  USE ISO_FORTRAN_ENV, ONLY: INT8
  IMPLICIT NONE
  INTEGER, PARAMETER::DP = SELECTED_REAL_KIND(15, 307)
CONTAINS
  SUBROUTINE print_kinds_info()
    USE io_param, ONLY: stdout
    WRITE (stdout, '(A)') 'Data Kinds Information:'

    WRITE (stdout, 3001) 'REAL: ', 'DP'
    WRITE (stdout, 3002) 'Kind value:', KIND(0.0_DP)
    WRITE (stdout, 3003) 'Precision:', PRECISION(0.0_DP)
    WRITE (stdout, 3003) 'negligible relative to 1:', EPSILON(0.0_DP)
    WRITE (stdout, 3003) 'Smallest positive:', TINY(0.0_DP)
    WRITE (stdout, 3003) 'Largest number:', HUGE(0.0_DP)

    WRITE (stdout, 3001) 'INTEGER: ', '(default)'
    WRITE (stdout, 3002) 'Kind value:', KIND(0)
    WRITE (stdout, 3002) 'Bit size:', BIT_SIZE(0)
    WRITE (stdout, 3002) 'Largest number:', HUGE(0)

    WRITE (stdout, 3001) 'LOGICAL: ', '(default)'
    WRITE (stdout, 3002) 'Kind value:', KIND(.TRUE.)

    WRITE (stdout, 3001) 'CHARACTER: ', '(default)'
    WRITE (stdout, 3002) 'Kind value:', KIND('C')

3001 FORMAT('= ', A, T30, A)
3002 FORMAT(4X, A, T30, I0)
3003 FORMAT(4X, A, T30, 1PE12.3E3)

  END SUBROUTINE print_kinds_info
  !
  FUNCTION eq_real(A, B, tol) RESULT(retval)
    REAL(DP), INTENT(IN) :: A, B
    REAL(DP), INTENT(IN), OPTIONAL:: tol
    LOGICAL :: retval
    REAL(DP) :: tolerance
    IF (PRESENT(tol)) THEN
      tolerance = tol
    ELSE
      tolerance = 1.0D-12
    END IF
    retval = ABS(A - B) < tolerance
  END FUNCTION eq_real
  !
  FUNCTION eq_vec_real(A, B, tol) RESULT(retval)
    REAL(DP), INTENT(IN) :: A(3), B(3)
    REAL(DP), INTENT(IN), OPTIONAL:: tol
    LOGICAL :: retval
    REAL(DP) :: tolerance
    INTEGER :: i
    IF (PRESENT(tol)) THEN
      tolerance = tol
    ELSE
      tolerance = 1.0D-12
    END IF
    retval = ALL(ABS(A - B) < tolerance)
  END FUNCTION eq_vec_real
END MODULE kinds
