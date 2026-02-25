SUBROUTINE check_Hermiticity(Nw, nkpt, Ham, herm_abs_max, h_abs_max, herm_rel)
  USE kinds, ONLY: DP
  INTEGER, INTENT(IN)::Nw, nkpt
  COMPLEX(DP), INTENT(IN)::Ham(Nw, Nw, nkpt)
  REAL(DP), INTENT(OUT)::herm_abs_max, h_abs_max, herm_rel
  INTEGER::ikpt
  COMPLEX(DP), ALLOCATABLE::antiherm(:, :)
  REAL(DP), PARAMETER::eps = 1.0D-14
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
END SUBROUTINE check_Hermiticity
