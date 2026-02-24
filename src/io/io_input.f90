MODULE io_input
  IMPLICIT NONE
  PRIVATE
  PUBLIC::read_input
  !
CONTAINS
  SUBROUTINE read_input()
    USE io_w90, ONLY: read_w90
    !
    CALL read_control()
    CALL read_w90()
  END SUBROUTINE read_input
  !
  SUBROUTINE read_control()
    USE io_global, ONLY: ionode, stdin, stdout, prefix
    NAMELIST /control/ prefix
    !
    WRITE (stdout, '(2X, A)') 'Reading &CONTROL Namelist...'
    IF (ionode) THEN
      READ (stdin, nml=control)
    END IF
  END SUBROUTINE read_control
END MODULE io_input
