MODULE dump_vec_io
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  PUBLIC::dump_r, dump_c, dump_i, dump_l
CONTAINS
  SUBROUTINE dump_r(fname, vec)
    USE io_global, ONLY: get_free_unit
    CHARACTER(LEN=*), INTENT(IN) :: fname
    REAL(DP), INTENT(IN) :: vec(:)
    INTEGER :: iu, ios, i
    iu = get_free_unit()
    OPEN (unit=iu, file=TRIM(fname), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'dump_r', 'Failed to open '//TRIM(fname))
    WRITE (iu, '(I0)') SIZE(vec)
    DO i = 1, SIZE(vec)
      WRITE (iu, '(ES24.16E3)') vec(i)
    END DO
    CLOSE (iu)
  END SUBROUTINE dump_r

  SUBROUTINE dump_c(fname, vec)
    USE io_global, ONLY: get_free_unit
    CHARACTER(LEN=*), INTENT(IN) :: fname
    COMPLEX(DP), INTENT(IN) :: vec(:)
    INTEGER :: iu, ios, i
    iu = get_free_unit()
    OPEN (unit=iu, file=TRIM(fname), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'dump_c', 'Failed to open '//TRIM(fname))
    WRITE (iu, '(I0)') SIZE(vec)
    DO i = 1, SIZE(vec)
      WRITE (iu, '(2(1X,ES24.16E3))') REAL(vec(i), DP), AIMAG(vec(i))
    END DO
    CLOSE (iu)
  END SUBROUTINE dump_c

  SUBROUTINE dump_i(fname, vec)
    USE io_global, ONLY: get_free_unit
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: vec(:)
    INTEGER :: iu, ios, i
    iu = get_free_unit()
    OPEN (unit=iu, file=TRIM(fname), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'dump_i', 'Failed to open '//TRIM(fname))
    WRITE (iu, '(I0)') SIZE(vec)
    DO i = 1, SIZE(vec)
      WRITE (iu, '(I0)') vec(i)
    END DO
    CLOSE (iu)
  END SUBROUTINE dump_i

  SUBROUTINE dump_l(fname, vec)
    USE io_global, ONLY: get_free_unit
    CHARACTER(LEN=*), INTENT(IN) :: fname
    LOGICAL, INTENT(IN) :: vec(:)
    INTEGER :: iu, ios, i, iv
    iu = get_free_unit()
    OPEN (unit=iu, file=TRIM(fname), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'dump_l', 'Failed to open '//TRIM(fname))
    WRITE (iu, '(I0)') SIZE(vec)
    DO i = 1, SIZE(vec)
      IF (vec(i)) THEN
        iv = 1
      ELSE
        iv = 0
      END IF
      WRITE (iu, '(I0)') iv
    END DO
    CLOSE (iu)
  END SUBROUTINE dump_l
END MODULE dump_vec_io
