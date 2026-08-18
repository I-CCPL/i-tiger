MODULE degen_mod
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  LOGICAL, ALLOCATABLE, PRIVATE::l_degen(:, :, :)
CONTAINS
  SUBROUTINE degen_init()
    ALLOCATE (l_degen(Nw, Nw, t_kpt%nkpt))
  END SUBROUTINE degen_init
!
  SUBROUTINE degen_compute(iks, dE_thr, eigval)
    INTEGER, INTENT(IN)::iks
    REAL(DP), INTENT(IN)::dE_thr
    REAL(DP), INTENT(IN)::eigval(Nw)
    INTEGER::n, m
    l_degen(:, :, iks) = .FALSE.
    DO n = 1, Nw - 1
      DO m = n + 1, Nw
        IF (ABS(eigval(n) - eigval(m)) < dE_thr) THEN
          l_degen(n, m, iks) = .TRUE.
          l_degen(m, n, iks) = .TRUE.
        END IF
      END DO
    END DO
  END SUBROUTINE degen_compute
  !
  SUBROUTINE degen_write()
    USE io_global, ONLY: ionode, get_free_unit
    INTEGER::iks, n, m
    INTEGER::io_unit
    LOGICAL, ALLOCATABLE::l_degen_tot(:, :, :)
    IF (ionode) THEN
      ALLOCATE (l_degen_tot(Nw, Nw, t_kpt%nktot))
    ELSE
      ALLOCATE (l_degen_tot(0, 0, 0))
    END IF
    CALL t_kpt%gather_l(Nw*Nw, l_degen, l_degen_tot)
    DEALLOCATE (l_degen)

    IF (ionode) THEN
      io_unit = get_free_unit()
      OPEN (unit=io_unit, file='itg.degen.dat', &
            status='replace', action='write')

      WRITE (io_unit, '(A)') '# n, m, iks'
      DO n = 1, Nw - 1
        DO m = n + 1, Nw
          DO iks = 1, t_kpt%nktot
            IF (l_degen_tot(n, m, iks)) THEN
              WRITE (io_unit, '(I0, 1X, I0, 1X, I0)') n, m, iks
            END IF
          END DO
        END DO
      END DO

      CLOSE (io_unit)
    END IF
    DEALLOCATE (l_degen_tot)
  END SUBROUTINE degen_write
END MODULE degen_mod
