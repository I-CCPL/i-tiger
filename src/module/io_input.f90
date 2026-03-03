MODULE io_input
  IMPLICIT NONE
  PRIVATE
  PUBLIC::read_input
  !
CONTAINS
  SUBROUTINE read_input()
    USE io_global, ONLY: stdout, prefix, write_sep_line
    USE wannier90, ONLY: w90data, read_w90
    !
    CALL read_control()
    w90data%prefix = TRIM(prefix)
    CALL read_w90()
    WRITE (stdout, '(2X, A)') 'Input reading completed.'
    CALL write_sep_line()
  END SUBROUTINE read_input
  !
  SUBROUTINE read_control()
    USE io_global, ONLY: ionode, stdin, stdout, &
                         prefix, debug
    NAMELIST /control/ prefix, debug
    !
    WRITE (stdout, '(2X, A)') 'Reading &CONTROL Namelist...'
    IF (ionode) THEN
      READ (stdin, nml=control)
      IF (debug) THEN
        WRITE (stdout, '(2X, A)') 'Debug mode is ON.'
        CALL read_debug()
      END IF
    END IF
  CONTAINS
    SUBROUTINE read_debug()
      USE wannier90, ONLY: chk_w90, Hq_band
      NAMELIST /debug/ chk_w90, Hq_band
      WRITE (stdout, '(2X, A)') 'Reading &DEBUG Namelist...'
      READ (stdin, nml=debug)
    END SUBROUTINE read_debug
  END SUBROUTINE read_control
END MODULE io_input
