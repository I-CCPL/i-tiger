MODULE lin_eig_H
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  INTERFACE eig_H
    MODULE PROCEDURE eig_cHk, eig_cH
  END INTERFACE eig_H
  PUBLIC :: eig_H, write_band
CONTAINS
  SUBROUTINE write_band(Hk, eigval, eigvec)
    USE io_global, ONLY: stdout, get_free_unit
    USE system, ONLY: Nw
    USE kpoints, ONLY: t_kpt, t_iks
    COMPLEX(DP), INTENT(IN)::Hk(Nw, Nw, t_kpt%nkpt)
    REAL(DP), INTENT(OUT)::eigval(Nw, t_kpt%nkpt)
    COMPLEX(DP), INTENT(OUT)::eigvec(Nw, Nw, t_kpt%nkpt)
    INTEGER::iw, ibnd, io_unit
    REAL(DP)::k_pos(t_kpt%nkpt)
    !
    k_pos(1) = 0.0_DP
    DO t_iks = 2, t_kpt%nkpt
      k_pos(t_iks) = k_pos(t_iks - 1) &
                     + SQRT(SUM((t_kpt%k_cart(:, t_iks) - t_kpt%k_cart(:, t_iks - 1))**2))
    END DO

    CALL eig_H(Nw, t_kpt%nkpt, Hk, eigval, eigvec)

    io_unit = get_free_unit()
    OPEN (unit=io_unit, file='itg.eigval')

    WRITE (stdout, '(2X, A)') '- Writing band structure data to "itg.eigval"...'
    WRITE (io_unit, '("#", A)') 'k_pos, eigval'
    DO iw = 1, Nw
      DO t_iks = 1, t_kpt%nkpt
        WRITE (io_unit, '(2F10.4)') k_pos(t_iks), eigval(iw, t_iks)
      END DO
      WRITE (io_unit, *) ! blank line
    END DO
    CLOSE (io_unit)
  END SUBROUTINE write_band
  !
  SUBROUTINE eig_cHk(ld_cH, nkpt, cHk, eigval, eigvec)
    !< Calculate eigenvalues and eigenvectors of Complex H in Wannier gauge at each k.
    USE io_global, ONLY: stdout
    INTEGER, INTENT(IN) :: ld_cH, nkpt
    COMPLEX(DP), INTENT(IN) :: cHk(ld_cH, ld_cH, nkpt)
    REAL(DP), INTENT(OUT) :: eigval(ld_cH, nkpt)
    COMPLEX(DP), INTENT(OUT) :: eigvec(ld_cH, ld_cH, nkpt)
    INTEGER::ikpt
    !
    WRITE (stdout, '(2X, A)') '- Diagonalizing matrix at each k...'
    DO ikpt = 1, nkpt
      CALL eig_cH(ld_cH, cHk(:, :, ikpt), eigval(:, ikpt), eigvec(:, :, ikpt))
    END DO
  END SUBROUTINE eig_cHk
  !
  SUBROUTINE eig_cH(ld_cH, cH, eigval, eigvec)
    !< Calculate eigenvalues and eigenvectors of Complex H in Wannier gauge at given k.
    USE mp_base, ONLY: mp_bcast
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
