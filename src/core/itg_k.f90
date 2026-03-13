MODULE itg_k
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  USE itg_R, ONLY: R_vec
  USE kpoints, ONLY: t_kpt, kpoint_type
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE::H_k_H(:, :, :)
  COMPLEX(DP), ALLOCATABLE::A_k_W(:, :, :, :), A_k_H(:, :, :, :)
  COMPLEX(DP), ALLOCATABLE::v_k(:, :, :, :)
  REAL(DP), ALLOCATABLE::L_k(:, :, :)
CONTAINS
  SUBROUTINE make_k_data()
    USE itg_R, ONLY: H_R, A_R
    USE fft_base, ONLY: fft_R2k
    USE lin_eig_H, ONLY: eig_H
    USE debug_data, ONLY: write_matrix, write_diag_matrix
    USE kpoints, ONLY: t_iks
    INTEGER::iw, jw

    ! TODO: k parallelization
    ALLOCATE (t_kpt%H_k(Nw, Nw))
    ALLOCATE (t_kpt%eigval(Nw, t_kpt%nkpt))
    ALLOCATE (t_kpt%eigvec(Nw, Nw, t_kpt%nkpt))
    DO t_iks = 1, t_kpt%nkpt
      ! Eigenvalues and eigenvectors
      CALL fft_R2k(R_vec, H_R, t_kpt%H_k(:, :))
      DO iw = 1, Nw
        t_kpt%H_k(iw, iw) = REAL(t_kpt%H_k(iw, iw), DP)
      END DO
      CALL eig_H(Nw, t_kpt%H_k, t_kpt%eigval(:, t_iks), t_kpt%eigvec(:, :, t_iks))

    END DO

    ALLOCATE (A_k_W(3, Nw, Nw, t_kpt%nkpt))
    ALLOCATE (A_k_H(3, Nw, Nw, t_kpt%nkpt))
    DO t_iks = 1, t_kpt%nkpt
      CALL fft_R2k(R_vec, A_R, A_k_W(:, :, :, t_iks))
    END DO
    ! CALL write_matrix('A_k_W.itg', A_k_W, t_kpt%nkpt, 1)

    CALL t_kpt%rotate(A_k_W, A_k_H)
    ! CALL write_matrix('A_k_H.itg', A_k_H, t_kpt%nkpt, 1)

    ALLOCATE (v_k(3, Nw, Nw, t_kpt%nkpt))
    CALL velocity(R_vec, H_R, A_k_H, v_k)

    ALLOCATE (L_k(3, Nw, t_kpt%nkpt))
    CALL OAM_mod_diag(t_kpt%eigval, v_k, L_k)
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
    IF (ALLOCATED(H_k_H)) DEALLOCATE (H_k_H)
    IF (ALLOCATED(A_k_W)) DEALLOCATE (A_k_W)
    IF (ALLOCATED(A_k_H)) DEALLOCATE (A_k_H)
    IF (ALLOCATED(v_k)) DEALLOCATE (v_k)
    IF (ALLOCATED(L_k)) DEALLOCATE (L_k)
  END SUBROUTINE clear_k_data
END MODULE itg_k
