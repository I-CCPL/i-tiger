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
CONTAINS
  SUBROUTINE NLO_init(t_kpt)
    USE constants, ONLY: zero, pi, zi, hbar_eVfs, &
                         e_chg_au, e_chg_si, FS2SEC, epsilon_0
    USE io_input, ONLY: NLO_Emin, NLO_dE
    USE system, ONLY: V_cell_3D
    USE kpoints, ONLY: kpoint_type
    TYPE(kpoint_type), INTENT(IN) :: t_kpt
    INTEGER::i

    ALLOCATE (NLO_hw(NLO_nE))
    ALLOCATE (epsilon_w(6, NLO_nE))
    ALLOCATE (JDOS_w(NLO_nE))
    ALLOCATE (shift_w(6, 3, NLO_nE))
    ALLOCATE (injection_w(6, 3, NLO_nE))

    IF (NLO_nE == 1) THEN
      NLO_hw(1) = NLO_Emin
    ELSE
      DO i = 1, NLO_nE
        NLO_hw(i) = NLO_Emin + REAL(i - 1, DP)*NLO_dE
      END DO
    END IF
    epsilon_w = zero
    JDOS_w = 0.0_DP
    shift_w = 0.0_DP
    injection_w = 0.0_DP

    !... dielectric function
    !> [e/V * 1/fs/Ang^3 * Ang*V/e] units
    !> kernel is [fs*Ang^2] units, so overall [1] units
    fac_dielec = zi*pi*e_chg_au**2/(hbar_eVfs*V_cell_3D) &
                 /epsilon_0*t_kpt%wk

    !... JDOS
    !> [1/eV * 1/fs] units
    !> kernel is [fs] units, so overall [1/eV] units
    fac_JDOS = 1/hbar_eVfs*t_kpt%wk

    !... shift current
    !> [e/fs * 1/V^2 * 1/fs/Ang^3] units
    !> kernel is [fs*Ang^3] units, so overall [e/fs * 1/V^2]
    fac_shift = -zi*pi*(e_chg_au**3)/(4.0_DP*(hbar_eVfs**2)*V_cell_3D) &
                *t_kpt%wk
    !> [e/fs] to [microA] units
    fac_shift = fac_shift &
                *1.0D6*e_chg_si/FS2SEC

    !... injection current
    !> [e/fs * 1/V^2 * 1/fs/Ang^3]
    !> kernel is [fs*Ang^2] units, so overall [1/fs * 1/V^2]
    fac_injection = pi*(e_chg_au**3)/((hbar_eVfs**2)*(pi**3)*V_cell_3D) &
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
    INTEGER :: io_unit, iom, ia
    CHARACTER(LEN=256) :: fname_a
    CHARACTER(LEN=1), PARAMETER :: a_lab(3) = (/'x', 'y', 'z'/)
    CALL mp_sum(epsilon_w)
    CALL mp_sum(JDOS_w)
    CALL mp_sum(shift_w)
    CALL mp_sum(injection_w)

    IF (.NOT. ionode) RETURN
    io_unit = get_free_unit()

    !... Dielectric constant
    OPEN (unit=io_unit, file='itg.epsilon_r.dat')
    CALL writing_info('dielectric function', 'itg.epsilon_r.dat')
    WRITE (io_unit, 0947) 'dielectric function units: [1]'
    WRITE (io_unit, 0948)
    DO iom = 1, NLO_nE
      WRITE (io_unit, 0949) NLO_hw(iom), DBLE(epsilon_w(:, iom))
    END DO
    CLOSE (io_unit)

    OPEN (unit=io_unit, file='itg.epsilon_i.dat')
    CALL writing_info('dielectric function', 'itg.epsilon_i.dat')
    WRITE (io_unit, 0947) 'dielectric function units: [1]'
    WRITE (io_unit, 0948)
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
      WRITE (io_unit, 0948)

      DO iom = 1, NLO_nE
        WRITE (io_unit, 0949) NLO_hw(iom), shift_w(:, ia, iom)
      END DO
      CLOSE (io_unit)
    END DO

    !... Injection current
    DO ia = 1, 3
      fname_a = 'itg.injection_'//a_lab(ia)//'.dat'
      OPEN (unit=io_unit, file=fname_a)
      CALL writing_info('injection current', fname_a)
      WRITE (io_unit, 0947) 'injection current units: [microA/V^2]'
      WRITE (io_unit, 0948)

      DO iom = 1, NLO_nE
        WRITE (io_unit, 0949) NLO_hw(iom), injection_w(:, ia, iom)
      END DO
      CLOSE (io_unit)
    END DO
0947 FORMAT("# ", A)
0948 FORMAT("# hw (eV), xx, xy, yy, yz, zz, zx")
0949 FORMAT(F13.6, 6(1X, ES16.8E3))
  END SUBROUTINE NLO_write
  !
  SUBROUTINE NLO_main(t_kpt, dH_bar, d2H_bar, A_bar, dA_bar, v_k_H)
    USE constants, ONLY: hbar_eVfs, zero, zi
    USE io_input, ONLY: dE_thr, dE_eta
    USE system, ONLY: Nw
    USE delta_func, ONLY: dE_inv
    USE kpoints, ONLY: kpoint_type, t_iks
    TYPE(kpoint_type), INTENT(IN) :: t_kpt
    COMPLEX(DP), INTENT(IN) :: dH_bar(3, Nw, Nw)
    COMPLEX(DP), INTENT(IN) :: d2H_bar(3, 3, Nw, Nw)
    COMPLEX(DP), INTENT(IN) :: A_bar(3, Nw, Nw)
    COMPLEX(DP), INTENT(IN) :: dA_bar(3, 3, Nw, Nw)
    COMPLEX(DP), INTENT(IN) :: v_k_H(3, Nw, Nw)
    INTEGER::n, m, p, a, b
    REAL(DP)::inv_hbar, eig_n, eig_m, dE_nm, w_inv(Nw, Nw), occ(Nw), fmn
    COMPLEX(DP)::v_bar(3, Nw, Nw), del_H_nm(3), del_bar(3), dv_bar
    COMPLEX(DP)::psum, dr_mn, da_mn
    COMPLEX(DP)::gen_r(3, Nw, Nw), gen_dr_mn(3, 3)
    CALL start_clock('NLO_main')

    inv_hbar = 1.0_DP/hbar_eVfs
    DO n = 1, Nw
      eig_n = t_kpt%eigval(n, t_iks)
      occ(n) = occ_T0(eig_n)
      DO m = 1, Nw
        eig_m = t_kpt%eigval(m, t_iks)
        dE_nm = eig_n - eig_m
        w_inv(n, m) = dE_inv(dE_nm, dE_eta)*hbar_eVfs

        v_bar(:, n, m) = dH_bar(:, n, m)*inv_hbar
        ! PRB 97, 245143 (2018) Eq. (22)
        gen_r(:, n, m) = -zi*v_bar(:, n, m)*w_inv(n, m) + A_bar(:, n, m)
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

        del_bar = v_bar(:, m, m) - v_bar(:, n, n)
        del_H_nm = v_k_H(:, n, n) - v_k_H(:, m, m)
        DO a = 1, 3
          DO b = 1, 3
            psum = zero
            dv_bar = d2H_bar(a, b, m, n)*inv_hbar
            DO p = 1, Nw
              IF (p == m .OR. p == n) CYCLE
              psum = psum &
                     + (v_bar(a, m, p)*v_bar(b, p, n))*w_inv(p, n) &
                     - (v_bar(b, m, p)*v_bar(a, p, n))*w_inv(m, p)
            END DO
            dr_mn = zi*w_inv(m, n) &
                    *( &
                    (v_bar(a, m, n)*del_bar(b) &
                     + v_bar(b, m, n)*del_bar(a))*w_inv(m, n) &
                    - dv_bar &
                    + psum &
                    )

            psum = zero
            DO p = 1, Nw
              psum = psum &
                     + (v_bar(b, m, p)*A_bar(a, p, n))*w_inv(m, p) &
                     - (A_bar(a, m, p)*v_bar(b, p, n))*w_inv(p, n)
            END DO
            da_mn = dA_bar(b, a, m, n) &
                    + psum
            ! PRB 97, 245143 (2018) Eq. (36)
            gen_dr_mn(a, b) = dr_mn + da_mn &
                              - (A_bar(b, m, m) - A_bar(b, n, n)) &
                              *(v_bar(a, m, n)*w_inv(m, n) + zi*A_bar(a, m, n))
          END DO
        END DO

        CALL dielectric(NLO_hw, epsilon_w, dE_nm, fmn, gen_r(:, n, m), gen_r(:, m, n))
        CALL Joint_DOS(NLO_hw, JDOS_w, dE_nm, fmn)
        CALL shift_current(NLO_hw, shift_w, dE_nm, fmn, gen_r(:, n, m), gen_dr_mn)
        CALL injection_current(NLO_hw, injection_w, dE_nm, fmn, &
                               del_H_nm, gen_r(:, n, m), gen_r(:, m, n))
      END DO
    END DO
    CALL stop_clock('NLO_main')
  END SUBROUTINE NLO_main
  !
  SUBROUTINE dielectric(hw, epsilon_w, dE_nm, fmn, r_nm, r_mn)
    USE io_input, ONLY: dE_eta
    USE delta_func, ONLY: delta_gaussian
    REAL(DP), INTENT(IN) :: hw(:)
    COMPLEX(DP), INTENT(OUT) :: epsilon_w(6, NLO_nE)
    REAL(DP), INTENT(IN) :: dE_nm, fmn
    COMPLEX(DP), INTENT(IN) :: r_nm(3), r_mn(3)
    INTEGER :: b, c, bc, iom
    REAL(DP) :: delta_E
    COMPLEX(DP)::kernel_mn(6)
    DO bc = 1, 6
      b = bc2b(bc)
      c = bc2c(bc)
      kernel_mn(bc) = fmn*fac_dielec &
                      *r_nm(b)*r_mn(c)
    END DO

    DO iom = 1, NLO_nE
      delta_E = delta_gaussian(dE_nm - hw(iom), dE_eta)
      epsilon_w(:, iom) = epsilon_w(:, iom) &
                          + delta_E*kernel_mn
    END DO
  END SUBROUTINE dielectric
  !
  SUBROUTINE Joint_DOS(hw, JDOS_w, dE_nm, fmn)
    USE io_input, ONLY: dE_eta
    USE delta_func, ONLY: delta_gaussian
    REAL(DP), INTENT(IN) :: hw(:)
    REAL(DP), INTENT(OUT) :: JDOS_w(NLO_nE)
    REAL(DP), INTENT(IN) :: dE_nm, fmn
    INTEGER :: iom
    REAL(DP) :: delta_E
    REAL(DP) :: kernel_mn

    kernel_mn = fmn*fac_JDOS
    DO iom = 1, NLO_nE
      delta_E = delta_gaussian(dE_nm - hw(iom), dE_eta)
      JDOS_w(iom) = JDOS_w(iom) + delta_E*kernel_mn
    END DO
  END SUBROUTINE Joint_DOS
  !
  SUBROUTINE shift_current(hw, shift_w, dE_nm, fmn, r_nm, dr_mn)
    USE io_input, ONLY: NLO_eta
    USE delta_func, ONLY: delta_gaussian
    REAL(DP), INTENT(IN) :: hw(:)
    REAL(DP), INTENT(INOUT) :: shift_w(6, 3, NLO_nE)
    REAL(DP), INTENT(IN) :: dE_nm, fmn
    COMPLEX(DP), INTENT(IN) :: r_nm(3), dr_mn(3, 3)
    INTEGER :: a, b, c, bc, iom
    REAL(DP) :: delta_E
    REAL(DP) :: kernel_mn(6, 3)
    !
    DO a = 1, 3
      DO bc = 1, 6
        b = bc2b(bc)
        c = bc2c(bc)
        kernel_mn(bc, a) = DBLE(fmn*fac_shift &
                                *(r_nm(b)*dr_mn(c, a) + r_nm(c)*dr_mn(b, a)))
      END DO
    END DO

    DO iom = 1, NLO_nE
      delta_E = delta_gaussian(dE_nm - hw(iom), NLO_eta) &
                + delta_gaussian(-dE_nm - hw(iom), NLO_eta)
      shift_w(:, :, iom) = shift_w(:, :, iom) + delta_E*kernel_mn(:, :)
    END DO
  END SUBROUTINE shift_current
  !
  SUBROUTINE injection_current(hw, injection_w, dE_nm, fmn, del_H_nm, r_nm, r_mn)
    USE io_input, ONLY: NLO_eta
    USE delta_func, ONLY: delta_gaussian
    REAL(DP), INTENT(IN) :: hw(:)
    REAL(DP), INTENT(INOUT) :: injection_w(6, 3, NLO_nE)
    REAL(DP), INTENT(IN) :: dE_nm, fmn
    COMPLEX(DP), INTENT(IN) :: del_H_nm(3), r_nm(3), r_mn(3)
    INTEGER :: a, b, c, bc, iom
    REAL(DP) :: delta_E
    REAL(DP) :: kernel_mn(6, 3)
    !
    DO a = 1, 3
      DO bc = 1, 6
        b = bc2b(bc)
        c = bc2c(bc)
        kernel_mn(bc, a) = DBLE(fmn*fac_injection &
                                *del_H_nm(a)*r_nm(b)*r_mn(c))
      END DO
    END DO
    injection_w(:, :, iom) = injection_w(:, :, iom) + delta_E*kernel_mn(:, :)
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
