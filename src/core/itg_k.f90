MODULE itg_k
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  USE itg_R, ONLY: R_vec
  USE kpoints, ONLY: t_kpt, kpoint_type
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE::A_k_W(:, :, :), A_k_H(:, :, :)
  !< Berry connection (3, Nw, Nw)
  COMPLEX(DP), ALLOCATABLE::dH_k_W(:, :, :), dH_k_H(:, :, :)
  !< Derivative of Hamiltonian (3, Nw, Nw)
  COMPLEX(DP), ALLOCATABLE::v_k(:, :, :)
  !< Velocity matrix (3, Nw, Nw)
  REAL(DP), ALLOCATABLE::L_k(:, :, :)
  !< OAM matrix (3, Nw, nktot)
CONTAINS
  SUBROUTINE make_k_data()
    USE itg_R, ONLY: H_R, A_R, dH_R
    USE fft_base, ONLY: fft_R2k
    USE lin_eig_H, ONLY: eig_H
    USE debug_data, ONLY: write_matrix, write_diag_matrix
    USE kpoints, ONLY: t_iks
    INTEGER::iw, jw

    ! TODO: k parallelization
    ALLOCATE (t_kpt%H_k(Nw, Nw))
    ALLOCATE (t_kpt%eigval(Nw, t_kpt%nkpt))
    ALLOCATE (t_kpt%eigvec(Nw, Nw))
    ALLOCATE (A_k_W(3, Nw, Nw))
    ALLOCATE (A_k_H(3, Nw, Nw))
    ALLOCATE (dH_k_W(3, Nw, Nw))
    ALLOCATE (dH_k_H(3, Nw, Nw))
    ALLOCATE (v_k(3, Nw, Nw))
    ALLOCATE (L_k(3, Nw, t_kpt%nkpt))
    DO t_iks = 1, t_kpt%nkpt
      ! Eigenvalues and eigenvectors
      CALL fft_R2k(R_vec, H_R, t_kpt%H_k(:, :))
      DO iw = 1, Nw
        t_kpt%H_k(iw, iw) = REAL(t_kpt%H_k(iw, iw), DP)
      END DO
      CALL eig_H(Nw, t_kpt%H_k, t_kpt%eigval(:, t_iks), t_kpt%eigvec(:, :))

      ! Derivative of Hamiltonian
      CALL fft_R2k(R_vec, dH_R, dH_k_W)
      CALL t_kpt%rotate(t_iks, dH_k_W, dH_k_H)
      ! Berry connection
      CALL fft_R2k(R_vec, A_R, A_k_W)
      CALL t_kpt%rotate(t_iks, A_k_W, A_k_H)
      CALL velocity(R_vec, A_k_H, dH_k_H, v_k)
      CALL OAM_mod_diag(t_kpt%eigval(:, t_iks), v_k, L_k(:, :, t_iks))
    END DO
  END SUBROUTINE make_k_data
  !
  SUBROUTINE write_k_data()
    USE io_global, ONLY: write_sep_line, stdout
    USE io_input, ONLY: lBand, lOAM
    USE io_output, ONLY: io_output_init, write_band, write_OAM
    CALL write_sep_line()
    WRITE (stdout, '(2X, A)') 'Write k data...'
    CALL io_output_init()
    IF (lBand) THEN
      CALL write_band('itg.band.dat', t_kpt%eigval)
    END IF
    IF (lOAM) THEN
      CALL write_OAM('itg.OAM.dat', L_k)
    END IF
  END SUBROUTINE write_k_data
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
