MODULE lin_eig_H
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  INTERFACE eig_H
    MODULE PROCEDURE eig_cH
  END INTERFACE eig_H
  PUBLIC :: eig_H
CONTAINS
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
