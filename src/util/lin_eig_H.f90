MODULE lin_eig_H
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eig_H, write_band
CONTAINS
  SUBROUTINE write_band(Hk, eigval)
    USE io_global, ONLY: stdout, get_free_unit
    USE system, ONLY: Nw
    USE kpoints, ONLY: t_kpt, t_iks
    COMPLEX(DP), INTENT(IN)::Hk(Nw, Nw, t_kpt%nkpt)
    REAL(DP), INTENT(IN)::eigval(Nw, t_kpt%nkpt)
    INTEGER::iw, ibnd, io_unit
    REAL(DP)::k_pos(t_kpt%nkpt)
    !
    k_pos(1) = 0.0_DP
    DO t_iks = 2, t_kpt%nkpt
      k_pos(t_iks) = k_pos(t_iks - 1) &
                     + SQRT(SUM((t_kpt%k_cart(:, t_iks) - t_kpt%k_cart(:, t_iks - 1))**2))
    END DO

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
  SUBROUTINE eig_H(ld_cH, nkpt, cHk, eigval, eigvec)
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
      CALL eig_zheevd(ld_cH, cHk(:, :, ikpt), eigval(:, ikpt), eigvec(:, :, ikpt))
    END DO
  END SUBROUTINE eig_H
  !
  SUBROUTINE eig_zheev(ld_cH, cH, eigval, eigvec)
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
    nb = ILAENV(1, 'ZHETRD', 'L', ld_cH, -1, -1, -1)
    IF (nb < 1 .OR. nb >= ld_cH) THEN
      lwork = 2*ld_cH
    ELSE
      lwork = (nb + 1)*ld_cH
    END IF
    ALLOCATE (work(lwork))
    ALLOCATE (rwork(MAX(1, 3*ld_cH - 2)))
    !
    eigvec = cH
    CALL zheev('V', 'L', ld_cH, eigvec, ld_cH, eigval, work, lwork, rwork, info)
    IF (info /= 0) THEN
      CALL errore(info, 'eig_H', 'zheevd failed')
    END IF
    !
    DEALLOCATE (work)
    DEALLOCATE (rwork)
  END SUBROUTINE eig_zheev
  SUBROUTINE eig_zheevd(ld_cH, cH, eigval, eigvec)
    !< Calculate eigenvalues and eigenvectors of Complex H in Wannier gauge at given k.
    USE mp_base, ONLY: mp_bcast
    INTEGER, INTENT(IN) :: ld_cH
    COMPLEX(DP), INTENT(IN) :: cH(ld_cH, ld_cH)
    REAL(DP), INTENT(OUT) :: eigval(ld_cH)
    COMPLEX(DP), INTENT(OUT) :: eigvec(ld_cH, ld_cH)
    !
    INTEGER::lwork, lrwork, liwork, info
    INTEGER, ALLOCATABLE :: iwork(:)
    REAL(DP), ALLOCATABLE :: rwork(:)
    COMPLEX(DP), ALLOCATABLE :: work(:)
    COMPLEX(DP) :: work_q(1)
    REAL(DP) :: rwork_q(1)
    INTEGER :: iwork_q(1)
    !
    ! Workspace query for ZHEEVD (divide-and-conquer Hermitian eigensolver)
    eigvec = cH
    CALL zheevd('V', 'L', ld_cH, eigvec, ld_cH, eigval, work_q, -1, rwork_q, -1, iwork_q, -1, info)
    IF (info /= 0) THEN
      CALL errore(info, 'eig_H', 'zheevd query failed')
    END IF
    lwork = MAX(1, INT(REAL(work_q(1), DP)))
    lrwork = MAX(1, INT(rwork_q(1)))
    liwork = MAX(1, iwork_q(1))

    ALLOCATE (work(lwork))
    ALLOCATE (rwork(lrwork))
    ALLOCATE (iwork(liwork))
    !
    ! Eigenvalue decomposition
    eigvec = cH
    CALL zheevd('V', 'L', ld_cH, eigvec, ld_cH, eigval, work, lwork, rwork, lrwork, iwork, liwork, info)
    IF (info /= 0) THEN
      CALL errore(info, 'eig_H', 'zheevd failed')
    END IF
    !
    DEALLOCATE (work)
    DEALLOCATE (rwork)
    DEALLOCATE (iwork)
  END SUBROUTINE eig_zheevd
END MODULE lin_eig_H
