  SUBROUTINE errore(code, routine, msg)
    USE mp_global, ONLY: mp_abort
    USE io_global, ONLY: stdout
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: code
    CHARACTER(len=*), INTENT(IN)::routine
    CHARACTER(len=*), INTENT(IN)::msg
    IF (code == 0) RETURN
    CALL write_bold_line()
    WRITE (stdout, '(2X, "Error in routine: ", A, " (", I0, ")")') TRIM(routine), code
    WRITE (stdout, '(4X, a)') TRIM(msg)
    CALL write_bold_line()
    CALL mp_abort(code)
  END SUBROUTINE errore
