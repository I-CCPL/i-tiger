SUBROUTINE check_Hermiticity(nkpt, ldX, Ham, tol)
  USE kinds, ONLY: DP
  USE io_global, ONLY: stdout
  USE system, ONLY: Nw
  IMPLICIT NONE
  INTEGER, INTENT(IN)::nkpt
  INTEGER, INTENT(IN)::ldX
  COMPLEX(DP), INTENT(IN)::Ham(ldX, Nw, Nw, nkpt)
  REAL(DP), INTENT(IN)::tol
  REAL(DP)::herm_abs_max, h_abs_max, herm_rel
  INTEGER::ikpt, iw, jw
  COMPLEX(DP), ALLOCATABLE::antiherm(:)
  REAL(DP), PARAMETER::eps = 1.0D-14
  CHARACTER(LEN=256)::msg
  !
  herm_abs_max = 0.0_DP
  h_abs_max = 0.0_DP
  !
  ALLOCATE (antiherm(ldX))
  DO ikpt = 1, nkpt
    DO jw = 1, Nw
      DO iw = jw, Nw
        antiherm(:) = Ham(:, iw, jw, ikpt) - CONJG(Ham(:, jw, iw, ikpt))
        herm_abs_max = MAX(herm_abs_max, MAXVAL(ABS(antiherm)))
        h_abs_max = MAX(h_abs_max, MAXVAL(ABS(Ham(:, iw, jw, ikpt))), MAXVAL(ABS(Ham(:, jw, iw, ikpt))))
      END DO
    END DO
  END DO
  DEALLOCATE (antiherm)
  !
  herm_rel = herm_abs_max/MAX(h_abs_max, eps)

  WRITE (*, '(2X, A, 1X, ES12.4E3)') '- hermiticity |H-H^+|_max:', herm_abs_max
  WRITE (*, '(2X, A, 1X, ES12.4E3)') '- max element magnitude  :', h_abs_max
  WRITE (*, '(2X, A, 1X, ES12.4E3)') '- hermiticity relative   :', herm_rel
  IF (herm_rel > tol) THEN
    WRITE (msg, '(A,1X,ES12.4E3)') '- hermiticity check failed. relative=', herm_rel
    CALL errore(1, 'check_Hermiticity', TRIM(msg))
  END IF
END SUBROUTINE check_Hermiticity
