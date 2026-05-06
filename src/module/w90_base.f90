SUBMODULE(wannier90) w90_base
  IMPLICIT NONE
CONTAINS
  MODULE SUBROUTINE read_w90_files(self)
    USE io_input, ONLY: lreq_mmn
    CLASS(w90data_type), INTENT(INOUT) :: self
    TYPE(chk_dum_type) :: chk_dum
    !
    CALL self%clear()
    CALL self%read_chk(chk_dum)
    CALL self%read_eig()

    IF (lreq_mmn) THEN
      CALL self%read_mmn()
    END IF

    IF (ionode .AND. chk_w90) THEN
      CALL write_chk_dump(self, chk_dum)
    END IF

    CALL chk_dum%clear()
  END SUBROUTINE read_w90_files

  MODULE SUBROUTINE clear_w90_data(self)
    CLASS(w90data_type), INTENT(INOUT) :: self
    self%nbnd = 0
    self%kpts%nkpt = 0
    self%nnb = 0

    IF (ALLOCATED(self%eigval)) DEALLOCATE (self%eigval)
    IF (ALLOCATED(self%kpts%k_cart)) DEALLOCATE (self%kpts%k_cart)
    IF (ALLOCATED(self%kpts%k_red)) DEALLOCATE (self%kpts%k_red)
    IF (ALLOCATED(self%v_matrix)) DEALLOCATE (self%v_matrix)
    IF (ALLOCATED(self%win_min)) DEALLOCATE (self%win_min)
    IF (ALLOCATED(self%ndimwin)) DEALLOCATE (self%ndimwin)
    IF (ALLOCATED(self%wannier_center_cart)) DEALLOCATE (self%wannier_center_cart)
    IF (ALLOCATED(self%wannier_spread)) DEALLOCATE (self%wannier_spread)
    IF (ALLOCATED(self%eigvec)) DEALLOCATE (self%eigvec)
    IF (ALLOCATED(self%Hq)) DEALLOCATE (self%Hq)
  END SUBROUTINE clear_w90_data

  MODULE FUNCTION wannier_gauge_diag(mat_H, v1, v2) RESULT(retval)
    USE constants, ONLY: cmplx_0
    REAL(DP), INTENT(IN) :: mat_H(:)
    COMPLEX(DP), INTENT(IN) :: v1(:), v2(:)
    COMPLEX(DP):: retval
    INTEGER::ibnd, nbnd
    nbnd = SIZE(v1)
    IF (SIZE(mat_H) /= nbnd .OR. SIZE(v2) /= nbnd) THEN
      CALL errore(1, 'wannier_gauge_diag', 'incompatible matrix and vector sizes')
    END IF

    retval = cmplx_0
    DO ibnd = 1, nbnd
      retval = retval + CONJG(v1(ibnd))*mat_H(ibnd)*v2(ibnd)
    END DO
  END FUNCTION wannier_gauge_diag

  MODULE FUNCTION wannier_gauge(mat_H, v1, v2) RESULT(retval)
    USE constants, ONLY: cmplx_0
    COMPLEX(DP), INTENT(IN) :: mat_H(:, :)
    COMPLEX(DP), INTENT(IN) :: v1(:), v2(:)
    COMPLEX(DP) :: retval
    INTEGER :: ibnd, jbnd, nbnd1, nbnd2
    nbnd1 = SIZE(v1)
    nbnd2 = SIZE(v2)
    IF (SIZE(mat_H, 1) /= nbnd1 .OR. SIZE(mat_H, 2) /= nbnd2) THEN
      CALL errore(1, 'wannier_gauge', 'incompatible matrix and vector sizes')
    END IF

    retval = cmplx_0
    DO ibnd = 1, nbnd1
      DO jbnd = 1, nbnd2
        retval = retval + CONJG(v1(ibnd))*mat_H(ibnd, jbnd)*v2(jbnd)
      END DO
    END DO
  END FUNCTION wannier_gauge

END SUBMODULE w90_base
