MODULE io_global
  USE io_param, ONLY: stdin, stdout, ionode, debug
  IMPLICIT NONE
  CHARACTER(LEN=256)::prefix
CONTAINS
  FUNCTION get_free_unit() RESULT(io_unit)
    INTEGER :: io_unit, idx
    INTEGER, PARAMETER::min_unit = 10, max_unit = 99
    LOGICAL::is_free
    DO idx = min_unit, max_unit
      INQUIRE (unit=idx, opened=is_free)
      IF (.NOT. is_free) THEN
        io_unit = idx
        RETURN
      END IF
    END DO
    CALL errore(1, 'get_free_unit', ' No free unit available.')
  END FUNCTION get_free_unit
  !
  SUBROUTINE check_file(fname)
    USE mp_base, ONLY: mp_bcast
    CHARACTER(LEN=*), INTENT(IN)::fname
    LOGICAL::lexists
    !
    INQUIRE (file=TRIM(fname), exist=lexists)
    CALL mp_bcast(lexists)
    IF (.NOT. lexists) THEN
      CALL errore(1, 'check_file', ' File not found: '//TRIM(fname))
    END IF
  END SUBROUTINE check_file
  !
END MODULE io_global
!
SUBROUTINE write_sep_line()
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  WRITE (stdout, '(A)') REPEAT('-', 50)
END SUBROUTINE write_sep_line
SUBROUTINE write_bold_line()
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  WRITE (stdout, '(A)') REPEAT('=', 50)
END SUBROUTINE write_bold_line
