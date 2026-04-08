MODULE lin_eig_H
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eig_H
CONTAINS
  !
  SUBROUTINE eig_H(ld_cH, cHk, eigval, eigvec)
    !< Calculate eigenvalues and eigenvectors of Complex H in Wannier gauge at each k.
    USE io_global, ONLY: stdout
    INTEGER, INTENT(IN) :: ld_cH
    COMPLEX(DP), INTENT(IN) :: cHk(ld_cH, ld_cH)
    REAL(DP), INTENT(OUT) :: eigval(ld_cH)
    COMPLEX(DP), INTENT(OUT) :: eigvec(ld_cH, ld_cH)
    !
    CALL eig_zheevd(ld_cH, cHk(:, :), eigval(:), eigvec(:, :))
    ! CALL eig_zheev(ld_cH, cHk(:, :), eigval(:), eigvec(:, :))
    ! CALL eig_zhpevx(ld_cH, cHk(:, :), eigval(:), eigvec(:, :))
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
    CALL errore(info, 'eig_H', 'zheev failed')
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
    CALL errore(info, 'eig_H', 'zheevd failed')
    !
    DEALLOCATE (work)
    DEALLOCATE (rwork)
    DEALLOCATE (iwork)
  END SUBROUTINE eig_zheevd

  SUBROUTINE eig_zhpevx(ld_cH, mat, eig, rot)
    USE kinds, ONLY: DP
    USE constants, ONLY: zero

    INTEGER, INTENT(in) :: ld_cH
    COMPLEX(kind=dp), INTENT(in) :: mat(ld_cH, ld_cH)
    REAL(kind=dp), INTENT(out) :: eig(ld_cH)
    COMPLEX(kind=dp), INTENT(out) :: rot(ld_cH, ld_cH)

    COMPLEX(kind=dp) :: mat_pack((ld_cH*(ld_cH + 1))/2), cwork(2*ld_cH)
    REAL(kind=dp) :: rwork(7*ld_cH)
    INTEGER :: i, j, info, nfound, iwork(5*ld_cH), ifail(ld_cH)
    CHARACTER(len=120) :: errormsg

    DO j = 1, ld_cH
      DO i = 1, j
        mat_pack(i + ((j - 1)*j)/2) = mat(i, j)
      END DO
    END DO
    rot = zero; eig = 0.0_DP; cwork = zero; rwork = 0.0_DP; iwork = 0
    CALL ZHPEVX('V', 'A', 'U', ld_cH, mat_pack, 0.0_DP, 0.0_DP, 0, 0, -1.0_DP, &
                nfound, eig(1), rot, ld_cH, cwork, rwork, iwork, ifail, info)
    CALL errore(info, 'eig_H', 'zhpevx failed')
  END SUBROUTINE eig_zhpevx
END MODULE lin_eig_H
