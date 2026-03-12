MODULE debug_data
  IMPLICIT NONE
CONTAINS
  SUBROUTINE write_diag_matrix(name, mat, ldim, ldX)
    USE kinds, ONLY: DP
    USE system, ONLY: Nw
    CHARACTER(LEN=*), INTENT(IN)::name
    REAL(DP), INTENT(IN)::mat(..)
    INTEGER, INTENT(IN)::ldim, ldX
    INTEGER::ldY, io_unit
    !
    OPEN (newunit=io_unit, file=TRIM(name), action='write')
    WRITE (io_unit, '(A)') '# idim, iw, idX, mat'
    ldY = SIZE(mat)/Nw/ldim/ldX
    CALL write_4D_matrix(io_unit, ldim, ldX, ldY, mat)
  END SUBROUTINE write_diag_matrix

  SUBROUTINE write_matrix(name, mat, ldim, ldX)
    USE kinds, ONLY: DP
    USE system, ONLY: Nw
    CHARACTER(LEN=*), INTENT(IN)::name
    INTEGER, INTENT(IN)::ldim, ldX
    COMPLEX(DP), INTENT(IN)::mat(..)
    INTEGER::ldY, io_unit
    !
    OPEN (newunit=io_unit, file=TRIM(name), action='write')
    WRITE (io_unit, '(A)') '# idim, iw, jw, idX, mat'
    ldY = SIZE(mat)/Nw/Nw/ldim/ldX
    CALL write_5D_matrix(io_unit, ldim, ldX, ldY, mat)
  END SUBROUTINE write_matrix
END MODULE debug_data

SUBROUTINE write_4D_matrix(io_unit, ldim, ldX, ldY, mat)
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  IMPLICIT NONE
  INTEGER, INTENT(IN)::io_unit
  INTEGER, INTENT(IN)::ldim, ldX, ldY
  REAL(DP), INTENT(IN)::mat(ldY, ldX, Nw, ldim)
  INTEGER::idim, iw, idX
  DO idim = 1, ldim
    DO iw = 1, Nw
      DO idX = 1, ldX
        WRITE (io_unit, '(3I4, *(SP, 1X, ES11.4))') idim, iw, idX, mat(:, idX, iw, idim)
      END DO
    END DO
  END DO
END SUBROUTINE write_4D_matrix

SUBROUTINE write_5D_matrix(io_unit, ldim, ldX, ldY, mat)
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  INTEGER, INTENT(IN)::io_unit
  INTEGER, INTENT(IN)::ldim, ldX, ldY
  COMPLEX(DP), INTENT(IN)::mat(ldY, ldX, Nw, Nw, ldim)
  INTEGER::idim, iw, jw, idX
  DO idim = 1, ldim
    DO iw = 1, Nw
      DO jw = 1, Nw
        DO idX = 1, ldX
          WRITE (io_unit, '(4I4, *(SP, 1X, 2ES11.4,"j"))') idim, iw, jw, idX, mat(:, idX, iw, jw, idim)
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE write_5D_matrix
