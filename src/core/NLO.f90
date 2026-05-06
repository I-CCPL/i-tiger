MODULE NLO
  ! Nonlinear optics
  ! Ref. PRB 61, 5337 (2000)
  USE kinds, ONLY: DP
  USE io_input, ONLY: NLO_nE
  IMPLICIT NONE
  REAL(DP), ALLOCATABLE::NLO_hw(:)
  !> dielectric function (epsilon_r)
  COMPLEX(DP), ALLOCATABLE::epsilon_w(:, :)
  !> joint density of states
  REAL(DP), ALLOCATABLE::JDOS_w(:)
  !> shift current (sigma)
  REAL(DP), ALLOCATABLE::shift_w(:, :, :)
  !> injection current (eta)
  REAL(DP), ALLOCATABLE::injection_w(:, :, :)
  !... factors
  COMPLEX(DP)::fac_dielec
  REAL(DP)::fac_JDOS
  COMPLEX(DP)::fac_shift
  REAL(DP)::fac_injection

  INTEGER, PARAMETER :: bc2b(6) = (/1, 1, 2, 2, 3, 3/)
  INTEGER, PARAMETER :: bc2c(6) = (/1, 2, 2, 3, 3, 1/)

  INTEGER::is_mn, is_nm, ie_mn, ie_nm
CONTAINS
  SUBROUTINE NLO_init(t_kpt)
    USE constants, ONLY: cmplx_0, pi, cmplx_i, hbar_eVfs, &
                         e_chg_au, e_chg_si, FS2SEC, epsilon_0
    USE io_input, ONLY: NLO_Emin, NLO_dE
    USE system, ONLY: V_cell_3D
    USE kpoints, ONLY: kpoint_type
    TYPE(kpoint_type), INTENT(IN) :: t_kpt
    INTEGER::i

    ALLOCATE (NLO_hw(NLO_nE))
    ALLOCATE (epsilon_w(6, NLO_nE))
    ALLOCATE (JDOS_w(NLO_nE))
    ALLOCATE (shift_w(3, 6, NLO_nE))
    ALLOCATE (injection_w(3, 3, NLO_nE))

    IF (NLO_nE == 1) THEN
      NLO_hw(1) = NLO_Emin
    ELSE
      DO i = 1, NLO_nE
        NLO_hw(i) = NLO_Emin + REAL(i - 1, DP)*NLO_dE
      END DO
    END IF
    epsilon_w = cmplx_0
    JDOS_w = 0.0_DP
    shift_w = 0.0_DP
    injection_w = 0.0_DP

    !... dielectric function
    !> [e/V * 1/fs/Ang^3 * Ang*V/e] units
    !> kernel is [fs*Ang^2] units, so overall [1] units
    fac_dielec = cmplx_i*pi*e_chg_au**2/(hbar_eVfs*V_cell_3D) &
                 /epsilon_0*t_kpt%wk

    !... JDOS
    !> [1/eV * 1/fs] units
    !> kernel is [fs] units, so overall [1/eV] units
    fac_JDOS = 1/hbar_eVfs*t_kpt%wk

    !... shift current
    !> [e/fs * 1/V^2 * 1/fs/Ang^3] units
    !> kernel is [fs*Ang^3] units, so overall [e/fs * 1/V^2]
    fac_shift = -cmplx_i*pi*(e_chg_au**3)/(4.0_DP*(hbar_eVfs**2)*V_cell_3D) &
                *t_kpt%wk
    !> [e/fs] to [microA] units
    fac_shift = fac_shift &
                *1.0D6*e_chg_si/FS2SEC

    !... injection current
    !> [e/fs * 1/V^2 * 1/fs/Ang^3]
    !> kernel is [fs*Ang^2] units, so overall [1/fs * 1/V^2]
    fac_injection = pi*(e_chg_au**3)/((hbar_eVfs**2)*V_cell_3D)/2 &
                    *t_kpt%wk
    !> [e/fs] to [microA] units
    fac_injection = fac_injection &
                    *1.0D6*e_chg_si/FS2SEC

    !... convert delta_E to delta_w
    fac_dielec = fac_dielec*hbar_eVfs
    fac_JDOS = fac_JDOS*hbar_eVfs
    fac_shift = fac_shift*hbar_eVfs
    fac_injection = fac_injection*hbar_eVfs
  END SUBROUTINE NLO_init
  SUBROUTINE NLO_clear()
    IF (ALLOCATED(NLO_hw)) DEALLOCATE (NLO_hw)
    IF (ALLOCATED(epsilon_w)) DEALLOCATE (epsilon_w)
    IF (ALLOCATED(JDOS_w)) DEALLOCATE (JDOS_w)
    IF (ALLOCATED(shift_w)) DEALLOCATE (shift_w)
    IF (ALLOCATED(injection_w)) DEALLOCATE (injection_w)

  END SUBROUTINE NLO_clear
  SUBROUTINE NLO_write(t_kpt)
    USE mp_base, ONLY: mp_sum
    USE io_global, ONLY: ionode, get_free_unit
    USE io_output, ONLY: writing_info
    USE kpoints, ONLY: kpoint_type
    TYPE(kpoint_type), INTENT(IN) :: t_kpt
    INTEGER :: io_unit, iom, ia, ibc
    CHARACTER(LEN=256) :: fname_a
    CHARACTER(LEN=1), PARAMETER :: a_lab(3) = (/'x', 'y', 'z'/)
    CALL mp_sum(epsilon_w)
    CALL mp_sum(JDOS_w)
    CALL mp_sum(shift_w)
    CALL mp_sum(injection_w)

    IF (.NOT. ionode) RETURN
    io_unit = get_free_unit()

    !... Dielectric constant
    ! OPEN (unit=io_unit, file='itg.epsilon_r.dat')
    ! CALL writing_info('dielectric function', 'itg.epsilon_r.dat')
    ! WRITE (io_unit, 0947) 'dielectric function units: [1]'
    ! WRITE (io_unit, 0948) "xx", "xy", "yy", "yz", "zz", "zx"
    ! DO iom = 1, NLO_nE
    !   WRITE (io_unit, 0949) NLO_hw(iom), DBLE(epsilon_w(:, iom))
    ! END DO
    ! CLOSE (io_unit)

    OPEN (unit=io_unit, file='itg.epsilon_i.dat')
    CALL writing_info('dielectric function', 'itg.epsilon_i.dat')
    WRITE (io_unit, 0947) 'dielectric function units: [1]'
    WRITE (io_unit, 0948) "xx", "xy", "yy", "yz", "zz", "zx"
    DO iom = 1, NLO_nE
      WRITE (io_unit, 0949) NLO_hw(iom), AIMAG(epsilon_w(:, iom))
    END DO
    CLOSE (io_unit)

    !... Joint density of states
    OPEN (unit=io_unit, file='itg.JDOS.dat')
    CALL writing_info('joint density of states', 'itg.JDOS.dat')
    WRITE (io_unit, 0947) 'joint density of states units: [states/eV]'
    WRITE (io_unit, 0947) 'hw (eV), JDOS'
    DO iom = 1, NLO_nE
      WRITE (io_unit, 0949) NLO_hw(iom), JDOS_w(iom)
    END DO
    CLOSE (io_unit)

    !... Shift current
    DO ia = 1, 3
      fname_a = 'itg.shift_'//a_lab(ia)//'.dat'
      OPEN (unit=io_unit, file=fname_a)
      CALL writing_info('shift current', fname_a)
      WRITE (io_unit, 0947) 'shift current units: [microA/V^2]'
      WRITE (io_unit, 0948) "xx", "xy", "yy", "yz", "zz", "zx"

      DO iom = 1, NLO_nE
        WRITE (io_unit, 0949) NLO_hw(iom), &
          (shift_w(ia, ibc, iom), ibc=1, 6)
      END DO
      CLOSE (io_unit)
    END DO

    !... Injection current
    DO ia = 1, 3
      fname_a = 'itg.injection_'//a_lab(ia)//'.dat'
      OPEN (unit=io_unit, file=fname_a)
      CALL writing_info('injection current', fname_a)
      WRITE (io_unit, 0947) 'injection current units: [microA/V^2]'
      WRITE (io_unit, 0948) 'xy', 'yz', 'zx'

      DO iom = 1, NLO_nE
        WRITE (io_unit, 0949) NLO_hw(iom), &
          (injection_w(ia, ibc, iom), ibc=1, 3)
      END DO
      CLOSE (io_unit)
    END DO
0947 FORMAT("# ", A)
0948 FORMAT("# hw (eV)", 6(",", A16))
0949 FORMAT(F13.6, 6(1X, ES16.8E3))
  END SUBROUTINE NLO_write
  !
  SUBROUTINE NLO_main(t_kpt, dH_bar, d2H_bar, A_bar, dA_bar, v_k_H)
    USE constants, ONLY: hbar_eVfs, cmplx_0, cmplx_i
    USE io_input, ONLY: dE_thr, dE_eta, &
                        NLO_Emin, NLO_Emax, NLO_dE, NLO_nE, NLO_eta, NLO_w_thr
    USE system, ONLY: Nw
    USE delta_func, ONLY: dE_inv, w1gauss
    USE kpoints, ONLY: kpoint_type, t_iks
    TYPE(kpoint_type), INTENT(IN) :: t_kpt
    COMPLEX(DP), INTENT(IN) :: dH_bar(Nw, Nw, 3)
    COMPLEX(DP), INTENT(IN) :: d2H_bar(Nw, Nw, 3, 3)
    COMPLEX(DP), INTENT(IN) :: A_bar(Nw, Nw, 3)
    COMPLEX(DP), INTENT(IN) :: dA_bar(Nw, Nw, 3, 3)
    COMPLEX(DP), INTENT(IN) :: v_k_H(Nw, Nw, 3)
    INTEGER::n, m, p, a, b, iom
    REAL(DP)::inv_hbar, eig_n, eig_m, dE_nm, w_inv(Nw, Nw), occ(Nw), fmn, E_inv(Nw, Nw)
    COMPLEX(DP)::v_bar(Nw, Nw, 3), del_H_nm(3), del_bar_mn(3), dv_bar
    COMPLEX(DP)::psum, dr_mn, da_mn
    COMPLEX(DP)::gen_r(Nw, Nw, 3), gen_dr_mn(3, 3)
    REAL(DP), ALLOCATABLE::delta_E(:, :, :)
    CALL start_clock('NLO_main')
    inv_hbar = 1.0_DP/hbar_eVfs
    ALLOCATE (delta_E(NLO_nE, Nw, Nw))
    DO n = 1, Nw
      eig_n = t_kpt%eigval(n, t_iks)
      occ(n) = occ_T0(eig_n)
      DO m = 1, Nw
        eig_m = t_kpt%eigval(m, t_iks)
        dE_nm = eig_n - eig_m

        is_nm = MAX(INT((dE_nm - NLO_w_thr*NLO_eta - NLO_Emin)/NLO_dE + 1), 1)
        ie_nm = MIN(INT((dE_nm + NLO_w_thr*NLO_eta - NLO_Emin)/NLO_dE + 1), NLO_nE)
        delta_E(is_nm:ie_nm, n, m) = w1gauss(ie_nm - is_nm, (dE_nm - NLO_hw(is_nm:ie_nm)), NLO_eta, 0)
        IF (m == n) THEN
          w_inv(n, m) = cmplx_0
          E_inv(n, m) = cmplx_0
        ELSE
          w_inv(n, m) = dE_inv(dE_nm, dE_eta)*hbar_eVfs
          E_inv(n, m) = 1.0_DP/dE_nm*hbar_eVfs
        END IF

        v_bar(n, m, :) = dH_bar(n, m, :)*inv_hbar
        ! PRB 97, 245143 (2018) Eq. (22)
        gen_r(n, m, :) = -cmplx_i*v_bar(n, m, :)*E_inv(n, m) + A_bar(n, m, :)
      END DO
    END DO

    DO m = 1, Nw
      eig_m = t_kpt%eigval(m, t_iks)
      DO n = 1, Nw
        ! Reduce fmn = 0 case.
        IF (m == n) CYCLE
        fmn = occ(m) - occ(n)
        IF (ABS(fmn) <= 1.0D-14) CYCLE
        eig_n = t_kpt%eigval(n, t_iks)
        dE_nm = eig_n - eig_m
        ! IF (ABS(dE_nm) <= dE_thr) CYCLE

        del_bar_mn = (dH_bar(m, m, :) - dH_bar(n, n, :))*inv_hbar
        del_H_nm = v_k_H(n, n, :) - v_k_H(m, m, :)
        DO a = 1, 3
          DO b = 1, 3
            psum = cmplx_0
            dv_bar = d2H_bar(m, n, b, a)*inv_hbar
            DO p = 1, Nw
              IF (p == m .OR. p == n) CYCLE
              psum = psum &
                     + (v_bar(m, p, a)*v_bar(p, n, b))*w_inv(p, n) &
                     - (v_bar(m, p, b)*v_bar(p, n, a))*w_inv(m, p)
            END DO
            dr_mn = cmplx_i*E_inv(m, n) &
                    *( &
                    (v_bar(m, n, a)*del_bar_mn(b) &
                     + v_bar(m, n, b)*del_bar_mn(a))*E_inv(m, n) &
                    - dv_bar &
                    + psum &
                    )

            psum = cmplx_0
            DO p = 1, Nw
              IF (p == m .OR. p == n) CYCLE
              psum = psum &
                     + (v_bar(m, p, b)*A_bar(p, n, a))*w_inv(m, p) &
                     - (A_bar(m, p, a)*v_bar(p, n, b))*w_inv(p, n)
            END DO
            da_mn = dA_bar(m, n, a, b) &
                    - (A_bar(m, m, a) - A_bar(n, n, a))*v_bar(m, n, b)*E_inv(m, n) &
                    + psum
            ! PRB 97, 245143 (2018) Eq. (36)
            gen_dr_mn(b, a) = dr_mn + da_mn &
                              - (A_bar(m, m, b) - A_bar(n, n, b)) &
                              *(v_bar(m, n, a)*E_inv(m, n) + cmplx_i*A_bar(m, n, a))
          END DO
        END DO
        is_mn = MAX(INT((-dE_nm - NLO_w_thr*NLO_eta - NLO_Emin)/NLO_dE + 1), 1)
        ie_mn = MIN(INT((-dE_nm + NLO_w_thr*NLO_eta - NLO_Emin)/NLO_dE + 1), NLO_nE)
        is_nm = MAX(INT((dE_nm - NLO_w_thr*NLO_eta - NLO_Emin)/NLO_dE + 1), 1)
        ie_nm = MIN(INT((dE_nm + NLO_w_thr*NLO_eta - NLO_Emin)/NLO_dE + 1), NLO_nE)

        CALL dielectric(epsilon_w, delta_E(:, n, m), fmn, gen_r(n, m, :), gen_r(m, n, :))
        CALL Joint_DOS(JDOS_w, delta_E(:, n, m), fmn)
        CALL shift_current(shift_w, delta_E(:, n, m), delta_E(:, m, n), fmn, gen_r(n, m, :), gen_dr_mn)
        CALL injection_current(injection_w, delta_E(:, n, m), fmn, &
                               del_H_nm, gen_r(n, m, :), gen_r(m, n, :))
      END DO
    END DO
    CALL stop_clock('NLO_main')
  END SUBROUTINE NLO_main
  !
  SUBROUTINE dielectric(epsilon_w, delta_Enm, fmn, r_nm, r_mn)
    COMPLEX(DP), INTENT(INOUT) :: epsilon_w(6, NLO_nE)
    REAL(DP), INTENT(IN) :: delta_Enm(:), fmn
    COMPLEX(DP), INTENT(IN) :: r_nm(3), r_mn(3)
    INTEGER :: b, c, bc, iom
    COMPLEX(DP)::pref, kernel_mn(6)
    pref = fmn*fac_dielec
    DO bc = 1, 6
      b = bc2b(bc)
      c = bc2c(bc)
      kernel_mn(bc) = pref*r_mn(b)*r_nm(c)
    END DO

    DO iom = is_nm, ie_nm
      epsilon_w(:, iom) = epsilon_w(:, iom) &
                          + delta_Enm(iom)*kernel_mn
    END DO
  END SUBROUTINE dielectric
  !
  SUBROUTINE Joint_DOS(JDOS_w, delta_Enm, fmn)
    REAL(DP), INTENT(INOUT) :: JDOS_w(NLO_nE)
    REAL(DP), INTENT(IN) :: delta_Enm(:), fmn
    INTEGER :: iom
    REAL(DP) :: kernel_mn

    kernel_mn = fmn*fac_JDOS
    DO iom = is_nm, ie_nm
      JDOS_w(iom) = JDOS_w(iom) + delta_Enm(iom)*kernel_mn
    END DO
  END SUBROUTINE Joint_DOS
  !
  SUBROUTINE shift_current(shift_w, delta_Enm, delta_Emn, fmn, r_nm, dr_mn)
    REAL(DP), INTENT(INOUT) :: shift_w(3, 6, NLO_nE)
    REAL(DP), INTENT(IN) :: delta_Enm(:), delta_Emn(:), fmn
    COMPLEX(DP), INTENT(IN) :: r_nm(3), dr_mn(3, 3)
    INTEGER :: a, b, c, bc, iom
    COMPLEX(DP)::pref
    REAL(DP) :: kernel_mn(3, 6)
    !
    pref = fmn*fac_shift
    DO a = 1, 3
      DO bc = 1, 6
        b = bc2b(bc)
        c = bc2c(bc)
        kernel_mn(a, bc) = DBLE(pref*(r_nm(b)*dr_mn(a, c) + r_nm(c)*dr_mn(a, b)))
      END DO
    END DO

    DO iom = is_nm, ie_nm
      shift_w(:, :, iom) = shift_w(:, :, iom) &
                           + delta_Enm(iom)*kernel_mn(:, :)
    END DO
    DO iom = is_mn, ie_mn
      shift_w(:, :, iom) = shift_w(:, :, iom) &
                           + delta_Emn(iom)*kernel_mn(:, :)
    END DO
  END SUBROUTINE shift_current
  !
  SUBROUTINE injection_current(injection_w, delta_Enm, fmn, del_H_nm, r_nm, r_mn)
    REAL(DP), INTENT(INOUT) :: injection_w(3, 3, NLO_nE)
    REAL(DP), INTENT(IN) :: delta_Enm(:), fmn
    COMPLEX(DP), INTENT(IN) :: del_H_nm(3), r_nm(3), r_mn(3)
    INTEGER :: a, b, c, bc, iom
    COMPLEX(DP)::pref
    REAL(DP) :: kernel_mn(3, 3)
    !
    pref = fmn*fac_injection
    DO a = 1, 3
      DO bc = 1, 3
        b = bc2b(bc*2)
        c = bc2c(bc*2)
        kernel_mn(a, bc) = DBLE(pref*del_H_nm(a) &
                                *(r_nm(b)*r_mn(c) - r_nm(c)*r_mn(b)))
      END DO
    END DO

    DO iom = is_nm, ie_nm
      injection_w(:, :, iom) = injection_w(:, :, iom) &
                               + delta_Enm(iom)*kernel_mn(:, :)
    END DO
  END SUBROUTINE injection_current
  !
  PURE REAL(DP) FUNCTION occ_T0(en)
    USE io_input, ONLY: E_fermi
    REAL(DP), INTENT(IN) :: en
    IF (en <= E_fermi) THEN
      occ_T0 = 1.0_DP
    ELSE
      occ_T0 = 0.0_DP
    END IF
  END FUNCTION occ_T0
END MODULE NLO
