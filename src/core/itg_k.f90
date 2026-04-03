MODULE itg_k
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  USE itg_R, ONLY: R_vec
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  !... X_bar = U^+ X U
  !... X_k_H = X_bar only for Gauge-covariant X
  COMPLEX(DP), ALLOCATABLE::A_k_W(:, :, :), A_bar(:, :, :)
  !< Berry connection (3, Nw, Nw)
  COMPLEX(DP), ALLOCATABLE::dA_k_W(:, :, :, :), dA_bar(:, :, :, :)
  !< Derivative of Berry connection (3, 3, Nw, Nw)
  COMPLEX(DP), ALLOCATABLE::dH_k_W(:, :, :), dH_bar(:, :, :)
  !< Derivative of Hamiltonian (3, Nw, Nw)
  COMPLEX(DP), ALLOCATABLE::d2H_k_W(:, :, :, :), d2H_bar(:, :, :, :)
  !< Second derivative of Hamiltonian (3, 3, Nw, Nw)
  COMPLEX(DP), ALLOCATABLE::v_k_H(:, :, :)
  !< Velocity matrix (3, Nw, Nw)
  REAL(DP), ALLOCATABLE::L_k(:, :, :)
  !< OAM matrix (3, Nw, nkpt)
  REAL(DP), ALLOCATABLE::O_k(:, :)
  !< Berry curvature matrix (3, Nw)
  REAL(DP), ALLOCATABLE::berry(:, :)
  REAL(DP), ALLOCATABLE::berry_k(:, :, :)
  !< Berry curvature (3, nkpt)
CONTAINS
  SUBROUTINE make_k_data()
    USE itg_R, ONLY: H_R, A_R, dH_R, dA_R, d2H_R
    USE fft_base, ONLY: fft_R2k
    USE io_input, ONLY: lOAM, lBerry, lNLO
    USE lin_eig_H, ONLY: eig_H
    USE wannier90, ONLY: lreq_mmn
    USE kpoints, ONLY: t_iks, t_kpt
    USE NLO, ONLY: NLO_main
    INTEGER::iw
    CALL start_clock('make_k_data')

    ! Eigenvalues and eigenvectors
    CALL fft_R2k(R_vec, H_R, t_kpt%H_k(:, :))
    DO iw = 1, Nw
      t_kpt%H_k(iw, iw) = REAL(t_kpt%H_k(iw, iw), DP)
    END DO
    CALL eig_H(Nw, t_kpt%H_k, t_kpt%eigval(:, t_iks), t_kpt%eigvec(:, :))

    IF (lreq_mmn) THEN
      ! Derivative of Hamiltonian
      CALL fft_R2k(R_vec, dH_R, dH_k_W)
      CALL t_kpt%rotate(dH_k_W, dH_bar)
      ! Berry connection
      CALL fft_R2k(R_vec, A_R, A_k_W)
      CALL t_kpt%rotate(A_k_W, A_bar)
      CALL velocity(A_bar, dH_bar, v_k_H)
    END IF

    IF (lOAM) THEN
      CALL OAM_mod_diag(t_kpt%eigval(:, t_iks), v_k_H, L_k(:, :, t_iks))
    END IF
    IF (lBerry) THEN
      CALL Berry_mod(t_kpt%eigval(:, t_iks), v_k_H, O_k(:, :))
      berry_k(:, :, t_iks) = O_k(:, :)
      CALL Berry_sum(t_kpt%eigval(:, t_iks), O_k(:, :), berry(:, t_iks))
    END IF
    IF (lNLO) THEN
      CALL fft_R2k(R_vec, dA_R, dA_k_W)
      CALL t_kpt%rotate(dA_k_W, dA_bar)

      CALL fft_R2k(R_vec, d2H_R, d2H_k_W)
      CALL t_kpt%rotate(d2H_k_W, d2H_bar)

      CALL NLO_main(t_kpt, dH_bar, d2H_bar, A_bar, dA_bar, v_k_H)
    END IF
    CALL stop_clock('make_k_data')
  END SUBROUTINE make_k_data
  !
  SUBROUTINE write_k_data()
    USE io_global, ONLY: stdout, ionode
    USE io_input, ONLY: lBand, lOAM, lBerry, lNLO
    USE io_output, ONLY: io_output_init, write_band, write_OAM, &
                         write_Berry, write_Berry_k
    USE NLO, ONLY: NLO_write
    REAL(DP), ALLOCATABLE::eigval(:, :)
    REAL(DP), ALLOCATABLE::L_k_tot(:, :, :)
    REAL(DP), ALLOCATABLE::berry_tot(:, :)
    REAL(DP), ALLOCATABLE::berry_k_tot(:, :, :)
    !
    CALL start_clock('write_k_data')
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

    IF (lBerry) THEN
      IF (ionode) THEN
        ALLOCATE (berry_tot(3, t_kpt%nktot))
        ALLOCATE (berry_k_tot(3, Nw, t_kpt%nktot))
      ELSE
        ALLOCATE (berry_tot(0, 0))
        ALLOCATE (berry_k_tot(0, 0, 0))
      END IF
      CALL t_kpt%gather(3, berry, berry_tot)
      CALL t_kpt%gather(3*Nw, berry_k, berry_k_tot)
      CALL write_Berry('itg.Berry.dat', berry_tot)
      CALL write_Berry_k('itg.Berry_k.dat', berry_k_tot)
      DEALLOCATE (berry_tot)
      DEALLOCATE (berry_k_tot)
    END IF

    IF (lNLO) THEN
      CALL NLO_write(t_kpt)
    END IF

    CALL write_sep_line()
    CALL stop_clock('write_k_data')
  END SUBROUTINE write_k_data
  !
  SUBROUTINE allocate_k_data()
    USE system, ONLY: Nw
    USE io_input, ONLY: lOAM, lBerry, lNLO
    USE wannier90, ONLY: lreq_mmn
    USE NLO, ONLY: NLO_init
    INTEGER::i
    CALL t_kpt%divide_k()
    ALLOCATE (t_kpt%H_k(Nw, Nw))
    ALLOCATE (t_kpt%eigval(Nw, t_kpt%nkpt))
    ALLOCATE (t_kpt%eigvec(Nw, Nw))

    IF (lreq_mmn) THEN
      ALLOCATE (A_k_W(3, Nw, Nw))
      ALLOCATE (A_bar(3, Nw, Nw))
      ALLOCATE (dH_k_W(3, Nw, Nw))
      ALLOCATE (dH_bar(3, Nw, Nw))
      ALLOCATE (v_k_H(3, Nw, Nw))
    END IF
    IF (lNLO) THEN
      ALLOCATE (dA_k_W(3, 3, Nw, Nw))
      ALLOCATE (dA_bar(3, 3, Nw, Nw))
      ALLOCATE (d2H_k_W(3, 3, Nw, Nw))
      ALLOCATE (d2H_bar(3, 3, Nw, Nw))
      CALL NLO_init(t_kpt)
    END IF

    IF (lOAM) ALLOCATE (L_k(3, Nw, t_kpt%nkpt))
    IF (lBerry) ALLOCATE (O_k(3, Nw))
    IF (lBerry) ALLOCATE (berry(3, t_kpt%nkpt))
    IF (lBerry) ALLOCATE (berry_k(3, Nw, t_kpt%nkpt))
    IF (ALLOCATED(L_k)) L_k = 0.0_DP
    t_kpt%eigval = 0.0_DP
  END SUBROUTINE allocate_k_data
  !
  SUBROUTINE clear_k_data()
    USE NLO, ONLY: NLO_clear
    IF (ALLOCATED(A_k_W)) DEALLOCATE (A_k_W)
    IF (ALLOCATED(A_bar)) DEALLOCATE (A_bar)
    IF (ALLOCATED(dH_k_W)) DEALLOCATE (dH_k_W)
    IF (ALLOCATED(dH_bar)) DEALLOCATE (dH_bar)
    IF (ALLOCATED(v_k_H)) DEALLOCATE (v_k_H)
    IF (ALLOCATED(L_k)) DEALLOCATE (L_k)
    IF (ALLOCATED(O_k)) DEALLOCATE (O_k)
    IF (ALLOCATED(berry)) DEALLOCATE (berry)
    IF (ALLOCATED(berry_k)) DEALLOCATE (berry_k)
    CALL NLO_clear()
  END SUBROUTINE clear_k_data
END MODULE itg_k
