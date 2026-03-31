MODULE NLO
  ! Nonlinear optics
  USE kinds, ONLY: DP
  IMPLICIT NONE
CONTAINS
  SUBROUTINE shift_k_H(eigval, A_k_H, v_k_H, hw, sigma_w)
    USE constants, ONLY: pi, zi, hbar_evfs, e_chg_au, e_chg_si, ev2j, FS2SEC, ry2ev
    USE io_input, ONLY: dE_thr, shift_eta
    USE system, ONLY: Nw, V_cell_3D
    USE kpoints, ONLY: t_kpt
    REAL(DP), INTENT(IN) :: eigval(Nw)
    COMPLEX(DP), INTENT(IN) :: A_k_H(3, Nw, Nw)
    COMPLEX(DP), INTENT(IN) :: v_k_H(3, Nw, Nw)
    REAL(DP), INTENT(IN) :: hw(:)
    REAL(DP), INTENT(INOUT) :: sigma_w(3, 6, SIZE(hw))
    INTEGER :: n, m, iom
    REAL(DP) :: dE, fmn, delta_w, conv_fact
    COMPLEX(DP) :: pref_raw, factor
    COMPLEX(DP) :: kernel_mn(3, 6)
    !
    CALL start_clock('shift_k_H')
    !> (e/fs * 1/V^2/fs/Ang^3) units
    pref_raw = -zi*pi*(e_chg_au**3)/(2.0_DP*(hbar_evfs**2)*V_cell_3D)
    !> (e/fs to muA) units
    conv_fact = 1.0D6*e_chg_si/FS2SEC
    factor = pref_raw*conv_fact
    !
    DO m = 1, Nw
      DO n = 1, Nw
        ! Avoid minus photon resonance
        dE = (eigval(m) - eigval(n))
        IF (dE <= dE_thr) CYCLE
        fmn = occ_T0(eigval(m)) - occ_T0(eigval(n))
        IF (ABS(fmn) <= 1.0D-14) CYCLE
        CALL shift_kernel_mn(m, n, dE, A_k_H, v_k_H, kernel_mn)
        DO iom = 1, SIZE(hw)
          !> (fs) units
          delta_w = delta_gaussian(dE - hw(iom), shift_eta)*hbar_evfs
          sigma_w(:, :, iom) = sigma_w(:, :, iom) &
                               + DBLE(factor*kernel_mn(:, :))*fmn*delta_w
        END DO
      END DO
    END DO
    CALL stop_clock('shift_k_H')
  END SUBROUTINE shift_k_H
  !
  SUBROUTINE shift_kernel_mn(m, n, dE_mn, A_k_H, v_k_H, kernel_mn)
    USE constants, ONLY: hbar_evfs, zero
    USE io_input, ONLY: dE_thr
    USE system, ONLY: Nw
    INTEGER, INTENT(IN) :: m, n
    REAL(DP), INTENT(IN) :: dE_mn
    COMPLEX(DP), INTENT(IN) :: A_k_H(3, Nw, Nw)
    !< (Ang) units
    COMPLEX(DP), INTENT(IN) :: v_k_H(3, Nw, Nw)
    !< (Ang/fs) units
    COMPLEX(DP), INTENT(OUT) :: kernel_mn(3, 6)
    !< (Ang^3) units
    REAL(DP) :: pref_as
    COMPLEX(DP) :: rb, rc, psum, d_as(3, 3)
    !< (Ang^2) units
    INTEGER :: ia, ib, ic, p, ibc
    INTEGER, PARAMETER :: ibc2ib(6) = (/1, 1, 2, 2, 3, 3/)
    INTEGER, PARAMETER :: ibc2ic(6) = (/1, 2, 2, 3, 3, 1/)
    !
    kernel_mn = zero
    ! Avoid small energy denominator
    IF (ABS(dE_mn) <= dE_thr) RETURN
    !
    ! d_mn;b^a = pref_as [ r^a_mn d^b_mn + r^b_mn d^a_mn + psum ]
    ! d^a_mn = v^a_mm - v^a_nn
    pref_as = -hbar_evfs/dE_mn
    DO ib = 1, 3
      DO ia = 1, 3
        ! psum = sum_p [ v^a_mp r^b_pn - r^b_mp v^a_pn ]
        psum = zero
        DO p = 1, Nw
          IF (p == m .OR. p == n) CYCLE
          psum = psum &
                 + v_k_H(ia, m, p)*A_k_H(ib, p, n) &
                 - A_k_H(ib, m, p)*v_k_H(ia, p, n)
        END DO
        d_as(ib, ia) = pref_as &
                       *( &
                       A_k_H(ia, m, n)*(v_k_H(ib, m, m) - v_k_H(ib, n, n)) &
                       + A_k_H(ib, m, n)*(v_k_H(ia, m, m) - v_k_H(ia, n, n)) &
                       + psum &
                       )
      END DO
    END DO
    !
    ! k^abc_mn=(r^b_nm)*(d_mn;a^c) + (r^c_nm)*(d_mn;a^b)
    DO ibc = 1, 6
      ib = ibc2ib(ibc)
      ic = ibc2ic(ibc)
      rb = A_k_H(ib, n, m)
      rc = A_k_H(ic, n, m)
      DO ia = 1, 3
        kernel_mn(ia, ibc) = rb*d_as(ia, ic) + rc*d_as(ia, ib)
      END DO
    END DO
  END SUBROUTINE shift_kernel_mn
  !
  PURE REAL(DP) FUNCTION delta_lorentz(x, eta)
    USE constants, ONLY: pi
    REAL(DP), INTENT(IN) :: x, eta
    delta_lorentz = eta/(pi*(x*x + eta*eta))
  END FUNCTION delta_lorentz
  !
  PURE REAL(DP) FUNCTION delta_gaussian(x, eta)
    USE constants, ONLY: tpi
    REAL(DP), INTENT(IN) :: x, eta
    delta_gaussian = EXP(-(x/eta)**2)/(SQRT(tpi)*eta)
  END FUNCTION delta_gaussian
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
