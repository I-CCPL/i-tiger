SUBROUTINE check_Hermiticity(nkpt, Ham, tol)
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  INTEGER, INTENT(IN)::nkpt
  COMPLEX(DP), INTENT(IN)::Ham(Nw, Nw, nkpt)
  REAL(DP), INTENT(IN)::tol
  REAL(DP)::herm_abs_max, h_abs_max, herm_rel
  INTEGER::ikpt
  COMPLEX(DP), ALLOCATABLE::antiherm(:, :)
  REAL(DP), PARAMETER::eps = 1.0D-14
  CHARACTER(LEN=256)::msg
  !
  herm_abs_max = 0.0_DP
  h_abs_max = 0.0_DP
  !
  ALLOCATE (antiherm(Nw, Nw))
  DO ikpt = 1, nkpt
    antiherm = Ham(:, :, ikpt) - TRANSPOSE(CONJG(Ham(:, :, ikpt)))
    herm_abs_max = MAX(herm_abs_max, MAXVAL(ABS(antiherm)))
    h_abs_max = MAX(h_abs_max, MAXVAL(ABS(Ham(:, :, ikpt))))
  END DO
  DEALLOCATE (antiherm)
  !
  herm_rel = herm_abs_max/MAX(h_abs_max, eps)

  WRITE (stdout, '(2X, A, 1X, ES12.4E3)') '- hermiticity |H-H^+|_max:', herm_abs_max
  WRITE (stdout, '(2X, A, 1X, ES12.4E3)') '- max element magnitude  :', h_abs_max
  WRITE (stdout, '(2X, A, 1X, ES12.4E3)') '- hermiticity relative   :', herm_rel
  IF (herm_rel > tol) THEN
    WRITE (msg, '(A,1X,ES12.4E3)') '- hermiticity check failed. relative=', herm_rel
    CALL errore(1, 'check_Hermiticity', TRIM(msg))
  END IF
END SUBROUTINE check_Hermiticity
