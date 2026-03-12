SUBROUTINE itg_warnings()
  USE io_global, ONLY: stdout, write_bold_line, write_sep_line
  IMPLICIT NONE
  CALL write_bold_line()
  WRITE (stdout, '(2X, A)') 'WARNING: This is a development version. Use with caution.'
  WRITE (stdout, '(2X, A)') '- Symmetric q-points are not implemented yet.'
  CALL write_bold_line()
END SUBROUTINE itg_warnings
