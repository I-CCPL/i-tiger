MODULE NLO
  ! Nonlinear optics
  USE kinds, ONLY: DP
  IMPLICIT NONE
CONTAINS
  SUBROUTINE shift_current(eigval, A_bar, v_bar, hw, sigma_w)
    USE constants, ONLY: pi, zi, hbar_evfs, e_chg, ev2j, SEC2FS
    USE io_input, ONLY: dE_thr, shift_eta
    USE system, ONLY: Nw, V_cell_nD
    REAL(DP), INTENT(IN) :: eigval(Nw)
    COMPLEX(DP), INTENT(IN) :: A_bar(3, Nw, Nw)
    COMPLEX(DP), INTENT(IN) :: v_bar(3, Nw, Nw)
    REAL(DP), INTENT(IN) :: hw(:)
    REAL(DP), INTENT(INOUT) :: sigma_w(3, 6, SIZE(hw))
    INTEGER :: n, m, iom
    REAL(DP) :: dE, fmn, delta_w, conv_fact
    COMPLEX(DP) :: pref_raw, factor
    COMPLEX(DP) :: kernel_mn(3, 6)
    !
    pref_raw = -zi*pi/(2.0_DP*hbar_evfs*V_cell_nD)
    conv_fact = 1.0D6*e_chg**3*SEC2FS/(ev2j**2)
    factor = pref_raw*conv_fact
    !
    DO m = 1, Nw
      DO n = 1, Nw
        ! Avoid minus photon resonance
        dE = eigval(m) - eigval(n)
        IF (dE <= dE_thr) CYCLE
        fmn = occ_T0(eigval(m)) - occ_T0(eigval(n))
        IF (ABS(fmn) <= 1.0D-14) CYCLE
        CALL shift_kernel_mn(m, n, dE, A_bar, v_bar, kernel_mn)
        DO iom = 1, SIZE(hw)
          delta_w = delta_lorentz(dE - hw(iom), shift_eta)
          sigma_w(:, :, iom) = sigma_w(:, :, iom) &
                               + DBLE(factor*kernel_mn(:, :))*fmn*delta_w
        END DO
      END DO
    END DO
  END SUBROUTINE shift_current
  !
  SUBROUTINE shift_kernel_mn(m, n, dE_nm, A_bar, v_bar, kernel_mn)
    USE constants, ONLY: hbar_evfs, zero
    USE io_input, ONLY: dE_thr
    USE system, ONLY: Nw
    INTEGER, INTENT(IN) :: m, n
    REAL(DP), INTENT(IN) :: dE_nm
    COMPLEX(DP), INTENT(IN) :: A_bar(3, Nw, Nw)
    COMPLEX(DP), INTENT(IN) :: v_bar(3, Nw, Nw)
    COMPLEX(DP), INTENT(OUT) :: kernel_mn(3, 6)
    REAL(DP) :: pref_as
    COMPLEX(DP) :: rab, rac, psum, d_as(3, 3)
    INTEGER :: ia, ib, ic, p, ibc
    INTEGER, PARAMETER :: ibc2ib(6) = (/1, 1, 2, 2, 3, 3/)
    INTEGER, PARAMETER :: ibc2ic(6) = (/1, 2, 2, 3, 3, 1/)
    !
    kernel_mn = zero
    ! Avoid small energy denominator
    IF (ABS(dE_nm) <= dE_thr) RETURN
    !
    ! r^a_nm;b
    pref_as = -hbar_evfs/dE_nm
    DO ib = 1, 3
      DO ia = 1, 3
        ! psum = sum_p [ v^a_mp r^b_pn - r^b_mp v^a_pn ]
        psum = zero
        DO p = 1, Nw
          IF (p == m .OR. p == n) CYCLE
          psum = psum &
                 + v_bar(ia, m, p)*A_bar(ib, p, n) &
                 - A_bar(ib, m, p)*v_bar(ia, p, n)
        END DO
        d_as(ia, ib) = pref_as &
                       *( &
                       A_bar(ia, m, n)*(v_bar(ib, m, m) - v_bar(ib, n, n)) &
                       + A_bar(ib, m, n)*(v_bar(ia, m, m) - v_bar(ia, n, n)) &
                       + psum &
                       )
      END DO
    END DO
    !
    DO ibc = 1, 6
      ib = ibc2ib(ibc)
      ic = ibc2ic(ibc)
      rab = A_bar(ib, n, m)
      rac = A_bar(ic, n, m)
      DO ia = 1, 3
        kernel_mn(ia, ibc) = rab*d_as(ia, ic) + rac*d_as(ia, ib)
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
  PURE REAL(DP) FUNCTION occ_T0(en)
    USE io_input, ONLY: Ef
    REAL(DP), INTENT(IN) :: en
    IF (en <= Ef) THEN
      occ_T0 = 1.0_DP
    ELSE
      occ_T0 = 0.0_DP
    END IF
  END FUNCTION occ_T0
END MODULE NLO
