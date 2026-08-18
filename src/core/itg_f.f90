MODULE itg_f
  !> Final properties to be calculated
  USE kinds, ONLY: DP
  USE f_params
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE::v_k_H(:, :, :)
  !< Velocity matrix (Nw, Nw, 3)
  REAL(DP), ALLOCATABLE::L_k(:, :, :)
  !< OAM matrix (Nw, 3, nktot)

  REAL(DP), ALLOCATABLE::berry_k_p(:, :, :)
  !< Berry curvature (Nw, 3, nktot)
  REAL(DP), ALLOCATABLE::berry_sum_p(:, :)
  !< Berry curvature (Nw, 3)
  REAL(DP), ALLOCATABLE::berry_k_g(:, :, :)
  !< Berry curvature (Nw, 3, nktot)
  REAL(DP), ALLOCATABLE::berry_sum_g(:, :)
  !< Berry curvature (Nw, 3)

  REAL(DP), ALLOCATABLE::BCD_sea(:, :, :)
  !< Berry curvature dipole from Fermi sea (3, 3, NLO_nE)
  REAL(DP), ALLOCATABLE::BCD_surf(:, :, :)
  !< Berry curvature dipole from Fermi surface (3, 3, NLO_nE)
CONTAINS
  SUBROUTINE set_f_flag()
    USE itg_k, ONLY: k_data

    IF (lOAM) lOAM_g = .TRUE.
    IF (lBerry) lBerry_p = .TRUE.
    IF (lBCD) lBCD_p = .TRUE.
    IF (lNLO) lNLO_g = .TRUE.

    IF (lBerry_p) THEN
      k_data%bO_bar = .TRUE.
      k_data%bdH_bar = .TRUE.
      k_data%bA_bar = .TRUE.
      k_data%bD_bar = .TRUE.
    END IF

    IF (lBerry_g) THEN
      k_data%bdH_bar = .TRUE.
      k_data%bA_bar = .TRUE.
    END IF

    IF (lBCD_p) THEN
      k_data%bO_bar = .TRUE.
      k_data%bA_bar = .TRUE.
      ! k_data%bdH_bar = .TRUE.
      k_data%bD_bar = .TRUE.
    END IF

    IF (lNLO_g .OR. lshift_E_g .OR. ldielec_E_g .OR. &
        lshift_k_g .OR. ldielec_k_g .OR. lshift_vec_g) THEN
      k_data%bdH_bar = .TRUE.
      k_data%bd2H_bar = .TRUE.
      k_data%bA_bar = .TRUE.
      k_data%bdA_bar = .TRUE.
      k_data%bD_bar = .TRUE.
      IF (lshift_k_g .OR. ldielec_k_g .OR. lshift_vec_g) THEN
        k_data%binit_eigval = .TRUE.
      END IF
    END IF
  END SUBROUTINE set_f_flag
  !
  SUBROUTINE allocate_f()
    USE itg_R, ONLY: R_data
    USE kpoints, ONLY: t_kpt
    USE f_params, ONLY: Ef_nE
    USE system, ONLY: Nw
    USE NLO_g, ONLY: NLO_g_init
    !
    IF (R_data%bA_R) ALLOCATE (v_k_H(Nw, Nw, 3))
    IF (lOAM_g) ALLOCATE (L_k(Nw, 3, t_kpt%nkpt))
    IF (lBerry_p) THEN
      ALLOCATE (berry_sum_p(3, t_kpt%nkpt))
      ALLOCATE (berry_k_p(Nw, 3, t_kpt%nkpt))
    END IF
    IF (lBerry_g) THEN
      ALLOCATE (berry_sum_g(3, t_kpt%nkpt))
      ALLOCATE (berry_k_g(Nw, 3, t_kpt%nkpt))
    END IF
    IF (lBCD_p) ALLOCATE (BCD_sea(3, 3, Ef_nE))
    IF (lNLO_g .OR. lshift_E_g .OR. ldielec_E_g .OR. &
        lshift_k_g .OR. ldielec_k_g .OR. lshift_vec_g) CALL NLO_g_init(t_kpt)
  END SUBROUTINE allocate_f
  !
  SUBROUTINE clear_f()
    USE kpoints, ONLY: t_kpt
    USE NLO_g, ONLY: NLO_g_clear
    IF (ALLOCATED(t_kpt%eigval)) DEALLOCATE (t_kpt%eigval)
    IF (ALLOCATED(v_k_H)) DEALLOCATE (v_k_H)
    IF (ALLOCATED(L_k)) DEALLOCATE (L_k)
    IF (ALLOCATED(berry_sum_p)) DEALLOCATE (berry_sum_p)
    IF (ALLOCATED(berry_k_p)) DEALLOCATE (berry_k_p)
    IF (ALLOCATED(berry_sum_g)) DEALLOCATE (berry_sum_g)
    IF (ALLOCATED(berry_k_g)) DEALLOCATE (berry_k_g)
    IF (ALLOCATED(BCD_sea)) DEALLOCATE (BCD_sea)
    IF (ALLOCATED(BCD_surf)) DEALLOCATE (BCD_surf)
    CALL NLO_g_clear()
  END SUBROUTINE clear_f
  !
  SUBROUTINE make_f()
    USE itg_R, ONLY: R_data
    USE f_params, ONLY: Ef_nE
    USE itg_k, ONLY: k_data
    USE kpoints, ONLY: t_iks, t_kpt
    USE NLO_g, ONLY: NLO_g_main
    USE system, ONLY: Nw
    INTEGER::n
    CALL start_clock('make_f')

    IF (k_data%bA_bar .AND. k_data%bdH_bar) THEN
      CALL compute_v_k_H(k_data%mA_bar, k_data%mdH_bar, v_k_H)
    END IF

    IF (lOAM_g) THEN
      CALL OAM_g_mod_diag(t_kpt%eigval(:, t_iks), v_k_H, L_k(:, :, t_iks))
    END IF

    IF (lBerry_p) THEN
      DO n = 1, Nw
        CALL get_berry_nk_p(n, t_kpt%eigval(:, t_iks), k_data%mA_bar, k_data%mdH_bar, &
                            berry_k_p(n, :, t_iks))
      END DO
      CALL get_berry_sum_p(t_kpt%eigval(:, t_iks), berry_k_p(:, :, t_iks), berry_sum_p(:, t_iks))
      ! CALL get_berry_k_p(t_kpt%eigval(:, t_iks), k_data%mO_bar, k_data%mA_bar, k_data%mdH_bar, berry_sum_p(:, t_iks))
    END IF

    IF (lBerry_g) THEN
      CALL get_berry_nk_g(t_kpt%eigval(:, t_iks), v_k_H, berry_k_g(:, :, t_iks))
      CALL get_berry_sum_g(t_kpt%eigval(:, t_iks), berry_k_g(:, :, t_iks), berry_sum_g(:, t_iks))
    END IF

    IF (lBCD_p) THEN
      !... Fermi surface

      IF (.NOT. ALLOCATED(BCD_surf)) THEN
        ALLOCATE (BCD_surf(3, 3, Ef_nE))
        BCD_surf = 0.0_DP
      END IF
      CALL compute_BCD_surf_p(k_data%mO_bar, k_data%mA_bar, k_data%mdH_bar, BCD_surf)

      !... Fermi sea
    END IF

    IF (lNLO_g .OR. lshift_E_g .OR. ldielec_E_g .OR. &
        lshift_k_g .OR. ldielec_k_g .OR. lshift_vec_g) THEN
      CALL NLO_g_main(t_kpt, k_data%mdH_bar, k_data%md2H_bar, &
                      k_data%mA_bar, k_data%mdA_bar, v_k_H, k_data%mD_bar)
    END IF

    CALL stop_clock('make_f')
  END SUBROUTINE make_f

  SUBROUTINE write_f()
    USE io_global, ONLY: stdout, ionode
    USE mp_base, ONLY: mp_sum
    USE f_params, ONLY: Ef_nE
    USE io_output, ONLY: writing_info
    USE io_output, ONLY: io_output_init, write_band, write_OAM, &
                         write_Berry, write_Berry_k, write_BCD
    USE NLO_g, ONLY: NLO_g_write
    USE degen_mod, ONLY: degen_write
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
    CALL writing_info('degen. states', 'itg.degen.dat')
    CALL degen_write()
    !
    IF (lBand) THEN
      IF (ionode) THEN
        ALLOCATE (eigval(Nw, t_kpt%nktot))
      ELSE
        ALLOCATE (eigval(0, 0))
      END IF
      CALL t_kpt%gather_r(Nw, t_kpt%eigval, eigval)
      CALL write_band('itg.band.dat', eigval)
      DEALLOCATE (eigval)
    END IF

    IF (lOAM_g) THEN
      IF (ionode) THEN
        ALLOCATE (L_k_tot(Nw, 3, t_kpt%nktot))
      ELSE
        ALLOCATE (L_k_tot(0, 0, 0))
      END IF
      CALL t_kpt%gather_r(Nw*3, L_k, L_k_tot)
      CALL write_OAM('itg.OAM.dat', L_k_tot)
      DEALLOCATE (L_k_tot)
    END IF

    IF (lBerry_p) THEN
      IF (ionode) THEN
        ALLOCATE (berry_tot(3, t_kpt%nktot))
        ALLOCATE (berry_k_tot(Nw, 3, t_kpt%nktot))
      ELSE
        ALLOCATE (berry_tot(0, 0))
        ALLOCATE (berry_k_tot(0, 0, 0))
      END IF
      CALL t_kpt%gather_r(3, berry_sum_p, berry_tot)
      CALL t_kpt%gather_r(Nw*3, berry_k_p, berry_k_tot)
      CALL write_Berry('itg.Berry.p.dat', berry_tot)
      CALL write_Berry_k('itg.Berry_k.p.dat', berry_k_tot)
      DEALLOCATE (berry_tot)
      DEALLOCATE (berry_k_tot)
    END IF

    IF (lBerry_g) THEN
      IF (ionode) THEN
        ALLOCATE (berry_tot(3, t_kpt%nktot))
        ALLOCATE (berry_k_tot(Nw, 3, t_kpt%nktot))
      ELSE
        ALLOCATE (berry_tot(0, 0))
        ALLOCATE (berry_k_tot(0, 0, 0))
      END IF
      CALL t_kpt%gather_r(3, berry_sum_g, berry_tot)
      CALL t_kpt%gather_r(Nw*3, berry_k_g, berry_k_tot)
      CALL write_Berry('itg.Berry.g.dat', berry_tot)
      CALL write_Berry_k('itg.Berry_k.g.dat', berry_k_tot)
      DEALLOCATE (berry_tot)
      DEALLOCATE (berry_k_tot)
    END IF

    IF (lBCD_p) THEN
      ! CALL mp_sum(BCD_sea)
      ! CALL write_BCD('itg.BCD.sea.p.dat', BCD_sea)
      CALL mp_sum(BCD_surf)
      CALL write_BCD('itg.BCD.surf.p.dat', BCD_surf)
    END IF

    IF (lNLO_g .OR. lshift_E_g .OR. ldielec_E_g .OR. &
        lshift_k_g .OR. ldielec_k_g .OR. lshift_vec_g) THEN
      CALL NLO_g_write(t_kpt)
    END IF

    CALL write_sep_line()
    CALL stop_clock('write_f')
  END SUBROUTINE write_f
END MODULE
