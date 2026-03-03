MODULE lin_eig_H
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  INTERFACE eig_H
    MODULE PROCEDURE eig_cH
  END INTERFACE eig_H
  PUBLIC :: eig_H, write_band
CONTAINS
  SUBROUTINE write_band(Hk, nkpt, eigval, eigvec)
    USE io_global, ONLY: get_free_unit
    USE system, ONLY: Nw
    COMPLEX(DP), INTENT(IN)::Hk(Nw, Nw, nkpt)
    INTEGER, INTENT(IN)::nkpt
    REAL(DP), INTENT(OUT)::eigval(Nw, nkpt)
    COMPLEX(DP), INTENT(OUT)::eigvec(Nw, Nw, nkpt)
    INTEGER::ikpt, iw, ibnd, io_unit
    !
    DO ikpt = 1, nkpt
      CALL eig_H(Nw, Hk(:, :, ikpt), &
                 eigval(:, ikpt), eigvec(:, :, ikpt))
    END DO

    io_unit = get_free_unit()
    OPEN (unit=io_unit, file='itg.eigval')

    WRITE (io_unit, '("#", A)') 'ibnd, ikpt, eigval'
    DO iw = 1, Nw
      DO ikpt = 1, nkpt
        WRITE (io_unit, '(I6, I6, ES13.4E3)') iw, ikpt, eigval(iw, ikpt)
      END DO
    END DO
    CLOSE (io_unit)
  END SUBROUTINE write_band
  !
  SUBROUTINE eig_cH(ld_cH, cH, eigval, eigvec)
    USE mp_base, ONLY: mp_bcast
    !< Calculate eigenvalues and eigenvectors of Complex H in Wannier gauge at each k.
    INTEGER, INTENT(IN) :: ld_cH
    COMPLEX(DP), INTENT(IN) :: cH(ld_cH, ld_cH)
    REAL(DP), INTENT(OUT) :: eigval(ld_cH)
    COMPLEX(DP), INTENT(OUT) :: eigvec(ld_cH, ld_cH)
    !
    INTEGER::lwork, nb, info
    REAL(DP), ALLOCATABLE :: rwork(:)
    COMPLEX(DP), ALLOCATABLE :: work(:)
    INTEGER, EXTERNAL :: ILAENV
    !
    ! optimal blocksize workspace
    nb = ILAENV(1, 'ZHETRD', 'U', ld_cH, -1, -1, -1)
    IF (nb < 1 .OR. nb >= ld_cH) THEN
      lwork = 2*ld_cH
    ELSE
      lwork = (nb + 1)*ld_cH
    END IF
    ALLOCATE (work(lwork))
    ALLOCATE (rwork(MAX(1, 3*ld_cH - 2)))
    !
    eigvec = cH
    CALL zheev('V', 'U', ld_cH, eigvec, ld_cH, eigval, work, lwork, rwork, info)
    !
    CALL mp_bcast(eigval)
    CALL mp_bcast(eigvec)
    !
    DEALLOCATE (work)
    DEALLOCATE (rwork)
  END SUBROUTINE eig_cH
END MODULE lin_eig_H
