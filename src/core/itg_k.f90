MODULE itg_k
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  USE itg_R, ONLY: R_vec
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE::A_k_W(:, :, :), A_k_H(:, :, :)
  !< Berry connection (3, Nw, Nw)
  COMPLEX(DP), ALLOCATABLE::dH_k_W(:, :, :), dH_k_H(:, :, :)
  !< Derivative of Hamiltonian (3, Nw, Nw)
  COMPLEX(DP), ALLOCATABLE::v_k(:, :, :)
  !< Velocity matrix (3, Nw, Nw)
  REAL(DP), ALLOCATABLE::L_k(:, :, :)
  !< OAM matrix (3, Nw, nkpt)
CONTAINS
  SUBROUTINE make_k_data()
    USE itg_R, ONLY: H_R, A_R, dH_R
    USE fft_base, ONLY: fft_R2k
    USE io_input, ONLY: lOAM
    USE lin_eig_H, ONLY: eig_H
    USE kpoints, ONLY: t_iks
    INTEGER::iw, jw

    ! Eigenvalues and eigenvectors
    CALL fft_R2k(R_vec, H_R, t_kpt%H_k(:, :))
    DO iw = 1, Nw
      t_kpt%H_k(iw, iw) = REAL(t_kpt%H_k(iw, iw), DP)
    END DO
    CALL eig_H(Nw, t_kpt%H_k, t_kpt%eigval(:, t_iks), t_kpt%eigvec(:, :))

    IF (lOAM) THEN
      ! Derivative of Hamiltonian
      CALL fft_R2k(R_vec, dH_R, dH_k_W)
      CALL t_kpt%rotate(dH_k_W, dH_k_H)
      ! Berry connection
      CALL fft_R2k(R_vec, A_R, A_k_W)
      CALL t_kpt%rotate(A_k_W, A_k_H)
      CALL velocity(R_vec, A_k_H, dH_k_H, v_k)
      CALL OAM_mod_diag(t_kpt%eigval(:, t_iks), v_k, L_k(:, :, t_iks))
    END IF
  END SUBROUTINE make_k_data
  !
  SUBROUTINE write_k_data()
    USE io_global, ONLY: stdout, ionode
    USE io_input, ONLY: lBand, lOAM
    USE io_output, ONLY: io_output_init, write_band, write_OAM
    REAL(DP), ALLOCATABLE::eigval(:, :)
    REAL(DP), ALLOCATABLE::L_k_tot(:, :, :)
    !
    WRITE (stdout, '(2X, A)') 'Write k data...'
    CALL io_output_init()
    !
    IF (lBand) THEN
      IF (ionode) THEN
        ALLOCATE (eigval(Nw, t_kpt%nktot))
      ELSE
        ALLOCATE (eigval(0, 0))
      END IF
      CALL t_kpt%gather(Nw, t_kpt%eigval, eigval)
      CALL write_band('itg.band.dat', eigval)
      DEALLOCATE (eigval)
    END IF

    IF (lOAM) THEN
      IF (ionode) THEN
        ALLOCATE (L_k_tot(3, Nw, t_kpt%nktot))
      ELSE
        ALLOCATE (L_k_tot(0, 0, 0))
      END IF
      CALL t_kpt%gather(3*Nw, L_k, L_k_tot)
      CALL write_OAM('itg.OAM.dat', L_k_tot)
      DEALLOCATE (L_k_tot)
    END IF

    CALL write_sep_line()
  END SUBROUTINE write_k_data
  !
  SUBROUTINE allocate_k_data()
    USE system, ONLY: Nw
    CALL t_kpt%divide_k()
    ALLOCATE (t_kpt%H_k(Nw, Nw))
    ALLOCATE (t_kpt%eigval(Nw, t_kpt%nkpt))
    ALLOCATE (t_kpt%eigvec(Nw, Nw))
    ALLOCATE (A_k_W(3, Nw, Nw))
    ALLOCATE (A_k_H(3, Nw, Nw))
    ALLOCATE (dH_k_W(3, Nw, Nw))
    ALLOCATE (dH_k_H(3, Nw, Nw))
    ALLOCATE (v_k(3, Nw, Nw))
    ALLOCATE (L_k(3, Nw, t_kpt%nkpt))
    L_k = 0.0_DP
    t_kpt%eigval = 0.0_DP
  END SUBROUTINE allocate_k_data
  !
  SUBROUTINE clear_k_data()
    IF (ALLOCATED(A_k_W)) DEALLOCATE (A_k_W)
    IF (ALLOCATED(A_k_H)) DEALLOCATE (A_k_H)
    IF (ALLOCATED(dH_k_W)) DEALLOCATE (dH_k_W)
    IF (ALLOCATED(dH_k_H)) DEALLOCATE (dH_k_H)
    IF (ALLOCATED(v_k)) DEALLOCATE (v_k)
    IF (ALLOCATED(L_k)) DEALLOCATE (L_k)
  END SUBROUTINE clear_k_data
END MODULE itg_k
