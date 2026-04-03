MODULE NLO
  ! Nonlinear optics
  USE kinds, ONLY: DP
  IMPLICIT NONE
  REAL(DP), ALLOCATABLE::shift_w(:, :, :)
  REAL(DP), ALLOCATABLE::shift_hw(:)
  !... factors
  COMPLEX(DP)::fac_dielec
  COMPLEX(DP)::fac_shift

CONTAINS
  SUBROUTINE NLO_init()
    USE constants, ONLY: pi, zi, hbar_evfs, e_chg_au, e_chg_si, FS2SEC
    USE io_input, ONLY: shift_nw, shift_wmin, shift_dw
    USE system, ONLY: V_cell_3D
    INTEGER::i
    REAL(DP)::conv_fact

    ALLOCATE (shift_w(3, 6, shift_nw))
    ALLOCATE (shift_hw(shift_nw))
    IF (shift_nw == 1) THEN
      shift_hw(1) = shift_wmin
    ELSE
      DO i = 1, shift_nw
        shift_hw(i) = shift_wmin + REAL(i - 1, DP)*shift_dw
      END DO
    END IF
    shift_w = 0.0_DP

    !... shift current
    !> [e/fs * 1/V^2 * 1/fs/Ang^3] units
    !> kernel is [fs*Ang^3] units, so overall [e/fs * 1/V^2]
    fac_shift = -zi*pi*(e_chg_au**3)/(4.0_DP*(hbar_evfs**2)*V_cell_3D)
    !> [e/fs] to [muA] units
    conv_fact = 1.0D6*e_chg_si/FS2SEC
    fac_shift = fac_shift*conv_fact

    !... convert delta_E to delta_w
    fac_shift = fac_shift*hbar_evfs
  END SUBROUTINE NLO_init
  SUBROUTINE NLO_clear()
    IF (ALLOCATED(shift_w)) DEALLOCATE (shift_w)
    IF (ALLOCATED(shift_hw)) DEALLOCATE (shift_hw)
  END SUBROUTINE NLO_clear
  SUBROUTINE NLO_write(t_kpt)
    USE mp_base, ONLY: mp_sum
    USE kpoints, ONLY: kpoint_type
    USE io_output, ONLY: write_shift
    TYPE(kpoint_type), INTENT(IN) :: t_kpt
    CALL mp_sum(shift_w)
    shift_w = shift_w*t_kpt%wk
    CALL write_shift('itg.shift', shift_hw, shift_w)
  END SUBROUTINE NLO_write
  !
  SUBROUTINE NLO_main(t_kpt, dH_bar, d2H_bar, A_bar, dA_bar)
    USE constants, ONLY: hbar_evfs, zero, zi
    USE io_input, ONLY: dE_thr, dE_eta
    USE system, ONLY: Nw
    USE delta_func, ONLY: dE_inv
    USE kpoints, ONLY: kpoint_type, t_iks
    TYPE(kpoint_type), INTENT(IN) :: t_kpt
    COMPLEX(DP), INTENT(IN) :: dH_bar(3, Nw, Nw)
    COMPLEX(DP), INTENT(IN) :: d2H_bar(3, 3, Nw, Nw)
    COMPLEX(DP), INTENT(IN) :: A_bar(3, Nw, Nw)
    COMPLEX(DP), INTENT(IN) :: dA_bar(3, 3, Nw, Nw)
    INTEGER::n, m, p, a, b
    REAL(DP)::dE, w_inv(Nw, Nw), occ(Nw), fmn
    COMPLEX(DP)::v_bar(3, Nw, Nw)
    COMPLEX(DP)::del_bar(3), dv_bar
    COMPLEX(DP)::psum, dr_mn, da_mn
    COMPLEX(DP)::gen_r_nm(3), gen_dr_mn(3, 3)
    CALL start_clock('NLO_main')

    DO n = 1, Nw
      occ(n) = occ_T0(t_kpt%eigval(n, t_iks))
      DO m = 1, Nw
        dE = t_kpt%eigval(n, t_iks) - t_kpt%eigval(m, t_iks)
        w_inv(n, m) = dE_inv(dE, dE_eta)*hbar_evfs
        v_bar(:, n, m) = dH_bar(:, n, m)/hbar_evfs
      END DO
    END DO

    DO m = 1, Nw
      DO n = 1, Nw
        IF (m == n) CYCLE
        fmn = occ(m) - occ(n)
        IF (ABS(fmn) <= 1.0D-14) CYCLE
        dE = t_kpt%eigval(m, t_iks) - t_kpt%eigval(n, t_iks)
        IF (ABS(dE) <= dE_thr) CYCLE

        del_bar(:) = v_bar(:, m, m) - v_bar(:, n, n)
        DO a = 1, 3
          ! PRB 97, 245143 (2018) Eq. (22)
          gen_r_nm(a) = v_bar(a, n, m)*w_inv(n, m)/zi + A_bar(a, n, m)
          DO b = 1, 3
            psum = zero
            dv_bar = d2H_bar(a, b, m, n)/hbar_evfs
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

        CALL shift_current(shift_hw, shift_w, dE, fmn, gen_r_nm, gen_dr_mn)
      END DO
    END DO
    CALL stop_clock('NLO_main')
  END SUBROUTINE NLO_main
  !
  SUBROUTINE shift_current(hw, sigma_w, dE, fmn, r_nm, dr_mn)
    USE io_input, ONLY: dE_thr, shift_eta
    USE delta_func, ONLY: delta_gaussian
    REAL(DP), INTENT(IN) :: hw(:)
    REAL(DP), INTENT(INOUT) :: sigma_w(3, 6, SIZE(hw))
    REAL(DP), INTENT(IN) :: dE, fmn
    COMPLEX(DP), INTENT(IN) :: r_nm(3), dr_mn(3, 3)
    INTEGER :: a, b, c, bc, iom
    INTEGER, PARAMETER :: bc2b(6) = (/1, 1, 2, 2, 3, 3/)
    INTEGER, PARAMETER :: bc2c(6) = (/1, 2, 2, 3, 3, 1/)
    REAL(DP) :: delta_E
    COMPLEX(DP) :: kernel_mn(3, 6)
    !
    DO a = 1, 3
      DO bc = 1, 6
        b = bc2b(bc)
        c = bc2c(bc)
        kernel_mn(a, bc) = (r_nm(b)*dr_mn(c, a) + r_nm(c)*dr_mn(b, a))
      END DO
    END DO

    DO iom = 1, SIZE(hw)
      delta_E = delta_gaussian(dE - hw(iom), shift_eta) &
                + delta_gaussian(-dE - hw(iom), shift_eta)
      sigma_w(:, :, iom) = sigma_w(:, :, iom) + fmn*delta_E*DBLE(fac_shift*kernel_mn(:, :))
    END DO
  END SUBROUTINE shift_current
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
