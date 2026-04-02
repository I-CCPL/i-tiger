SUBROUTINE check_nan(msg, val, idx)
  USE kinds, ONLY: DP
  USE, INTRINSIC::IEEE_ARITHMETIC, ONLY: ieee_is_nan
  CHARACTER(len=*), INTENT(IN)::msg
  COMPLEX(DP), INTENT(IN)::val
  INTEGER, INTENT(IN)::idx(:)
  IF (ieee_is_nan(DBLE(val)) .OR. ieee_is_nan(AIMAG(val))) THEN
    WRITE (*, *) 'NaN ', msg, ieee_is_nan(DBLE(val)), ieee_is_nan(AIMAG(val)), idx(:)
    STOP
  END IF
END SUBROUTINE check_nan
