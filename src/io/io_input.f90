MODULE io_input
  IMPLICIT NONE
  PRIVATE
  PUBLIC::read_input
  !
CONTAINS
  SUBROUTINE read_input()
    USE io_w90, ONLY: read_w90
    USE w90_chk_dump, ONLY: chk_dump, write_w90_chk_dump
    !
    CALL read_control()
    CALL read_w90()
    IF (chk_dump) CALL write_w90_chk_dump()
  END SUBROUTINE read_input
  !
  SUBROUTINE read_control()
    USE io_global, ONLY: ionode, stdin, stdout, prefix
    USE w90_chk_dump, ONLY: chk_dump, chk_dump_file
    NAMELIST /control/ prefix, chk_dump, chk_dump_file
    !
    WRITE (stdout, '(2X, A)') 'Reading &CONTROL Namelist...'
    IF (ionode) THEN
      READ (stdin, nml=control)
    END IF
  END SUBROUTINE read_control
END MODULE io_input
