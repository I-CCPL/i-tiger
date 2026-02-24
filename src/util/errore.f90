  SUBROUTINE errore(code, routine, msg)
    USE mp_global, ONLY: mp_abort
    USE io_global, ONLY: stdout
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: code
    CHARACTER(len=*), INTENT(IN)::routine
    CHARACTER(len=*), INTENT(IN), OPTIONAL::msg
    IF (code == 0) RETURN
    WRITE (stdout, '(2X, "Error in routine: ", A, " (", I0, ")")') TRIM(routine), code
    IF (PRESENT(msg)) THEN
      WRITE (stdout, '(4X, a)') TRIM(msg)
    END IF
    CALL mp_abort(code)
  END SUBROUTINE errore
  !
