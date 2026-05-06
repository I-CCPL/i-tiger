MODULE itg_f
  !> Final properties to be calculated
  USE kinds, ONLY: DP
  USE f_params, ONLY: lBand, lOAM, lBerry, lBCD, lNLO
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE::v_k_H(:, :, :)
  !< Velocity matrix (Nw, Nw, 3)
  REAL(DP), ALLOCATABLE::L_k(:, :, :)
  !< OAM matrix (Nw, 3, nktot)

  REAL(DP), ALLOCATABLE::berry(:, :)
  !< Berry curvature (Nw, 3)
  REAL(DP), ALLOCATABLE::berry_k(:, :, :)
  !< Berry curvature (Nw, 3, nktot)

  REAL(DP), ALLOCATABLE::BCD_sea(:, :, :)
  !< Berry curvature dipole from Fermi sea (3, 3, NLO_nE)
  REAL(DP), ALLOCATABLE::BCD_surf(:, :, :)
  !< Berry curvature dipole from Fermi surface (3, 3, NLO_nE)
CONTAINS
  SUBROUTINE set_f_flag()
    USE itg_k, ONLY: k_data
    IF (lBerry) THEN
      k_data%bO_bar = .TRUE.
      k_data%bdH_bar = .TRUE.
      k_data%bA_bar = .TRUE.
    END IF

    IF (lBCD) THEN
      k_data%bO_bar = .TRUE.
      k_data%bA_bar = .TRUE.
      k_data%bdH_bar = .TRUE.
    END IF

    IF (lNLO) THEN
      k_data%bdH_bar = .TRUE.
      k_data%bd2H_bar = .TRUE.
      k_data%bA_bar = .TRUE.
      k_data%bdA_bar = .TRUE.
    END IF
  END SUBROUTINE set_f_flag
  !
  SUBROUTINE allocate_f()
    USE itg_R, ONLY: R_data
    USE kpoints, ONLY: t_kpt
    USE f_params, ONLY: Ef_nE
    USE system, ONLY: Nw
    USE NLO, ONLY: NLO_init
    !
    IF (R_data%bA_R) ALLOCATE (v_k_H(Nw, Nw, 3))
    IF (lOAM) ALLOCATE (L_k(Nw, 3, t_kpt%nkpt))
    IF (lBerry) THEN
      ALLOCATE (berry(3, t_kpt%nkpt))
      ALLOCATE (berry_k(Nw, 3, t_kpt%nkpt))
    END IF
    IF (lBCD) ALLOCATE (BCD_sea(3, 3, Ef_nE))
    IF (lNLO) CALL NLO_init(t_kpt)
  END SUBROUTINE allocate_f
  !
  SUBROUTINE clear_f()
    USE NLO, ONLY: NLO_clear
    IF (ALLOCATED(v_k_H)) DEALLOCATE (v_k_H)
    IF (ALLOCATED(L_k)) DEALLOCATE (L_k)
    IF (ALLOCATED(berry)) DEALLOCATE (berry)
    IF (ALLOCATED(berry_k)) DEALLOCATE (berry_k)
    IF (ALLOCATED(BCD_sea)) DEALLOCATE (BCD_sea)
    IF (ALLOCATED(BCD_surf)) DEALLOCATE (BCD_surf)
    CALL NLO_clear()
  END SUBROUTINE clear_f
  !
  SUBROUTINE make_f()
    USE itg_R, ONLY: R_data
    USE f_params, ONLY: Ef_nE
    USE itg_k, ONLY: k_data
    USE kpoints, ONLY: t_iks, t_kpt
    USE NLO, ONLY: NLO_main
    CALL start_clock('make_f')

    IF (R_data%bA_R) THEN
      ! CALL velocity(k_data%mA_bar, k_data%mdH_bar, v_k_H)
    END IF

    IF (lOAM) THEN
      CALL OAM_mod_diag(t_kpt%eigval(:, t_iks), v_k_H, L_k(:, :, t_iks))
    END IF

    IF (lBerry) THEN
      ! CALL Berry_mod(t_kpt%eigval(:, t_iks), v_k_H, O_k(:, :))
      ! berry_k(:, :, t_iks) = O_k(:, :)
      ! CALL Berry_sum(t_kpt%eigval(:, t_iks), O_k(:, :), berry(:, t_iks))
      CALL Berry_proj(t_kpt%eigval(:, t_iks), k_data%mO_bar, k_data%mA_bar, k_data%mdH_bar, berry(:, t_iks))
    END IF

    IF (lBCD) THEN
      !... Fermi surface

      IF (.NOT. ALLOCATED(BCD_surf)) THEN
        ALLOCATE (BCD_surf(3, 3, Ef_nE))
      END IF
      CALL compute_BCD_surf(k_data%mO_bar, k_data%mA_bar, k_data%mdH_bar, BCD_surf)

      !... Fermi sea
      ! CALL compute_BCD_sea(t_kpt%eigval(:, t_iks), k_data%mdH_bar, k_data%md2H_bar, k_data%mA_bar, k_data%mdA_bar, k_data%mdO_bar, BCD_sea)
    END IF

    IF (lNLO) THEN
      CALL NLO_main(t_kpt, k_data%mdH_bar, k_data%md2H_bar, k_data%mA_bar, k_data%mdA_bar, v_k_H)
    END IF

    CALL stop_clock('make_f')
  END SUBROUTINE make_f

  SUBROUTINE write_f()
    USE io_global, ONLY: stdout, ionode
    USE mp_base, ONLY: mp_sum
    USE f_params, ONLY: Ef_nE
    USE io_output, ONLY: io_output_init, write_band, write_OAM, &
                         write_Berry, write_Berry_k, write_BCD
    USE NLO, ONLY: NLO_write
    USE system, ONLY: Nw
    USE kpoints, ONLY: t_kpt
    REAL(DP), ALLOCATABLE::eigval(:, :)
    REAL(DP), ALLOCATABLE::L_k_tot(:, :, :)
    REAL(DP), ALLOCATABLE::berry_tot(:, :)
    REAL(DP), ALLOCATABLE::berry_k_tot(:, :, :)
    !
    CALL start_clock('write_f')
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
        ALLOCATE (L_k_tot(Nw, 3, t_kpt%nktot))
      ELSE
        ALLOCATE (L_k_tot(0, 0, 0))
      END IF
      CALL t_kpt%gather(Nw*3, L_k, L_k_tot)
      CALL write_OAM('itg.OAM.dat', L_k_tot)
      DEALLOCATE (L_k_tot)
    END IF

    IF (lBerry) THEN
      IF (ionode) THEN
        ALLOCATE (berry_tot(3, t_kpt%nktot))
        ALLOCATE (berry_k_tot(Nw, 3, t_kpt%nktot))
      ELSE
        ALLOCATE (berry_tot(0, 0))
        ALLOCATE (berry_k_tot(0, 0, 0))
      END IF
      CALL t_kpt%gather(3, berry, berry_tot)
      CALL t_kpt%gather(Nw*3, berry_k, berry_k_tot)
      CALL write_Berry('itg.Berry.dat', berry_tot)
      CALL write_Berry_k('itg.Berry_k.dat', berry_k_tot)
      DEALLOCATE (berry_tot)
      DEALLOCATE (berry_k_tot)
    END IF

    IF (lBCD) THEN
      ! CALL mp_sum(BCD_sea)
      ! CALL write_BCD('itg.BCD.sea.dat', BCD_sea)
      CALL mp_sum(BCD_surf)
      CALL write_BCD('itg.BCD.surf.dat', BCD_surf)
    END IF

    IF (lNLO) THEN
      CALL NLO_write(t_kpt)
    END IF

    CALL write_sep_line()
    CALL stop_clock('write_f')
  END SUBROUTINE write_f
END MODULE
