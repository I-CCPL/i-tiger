MODULE NLO
  ! Nonlinear optics
  USE kinds, ONLY: DP
  IMPLICIT NONE
CONTAINS
  SUBROUTINE shift_current(hw, sigma_w, dE, fmn, kernel_mn)
    USE constants, ONLY: pi, zi, hbar_evfs, e_chg_au, e_chg_si, FS2SEC
    USE io_input, ONLY: dE_thr, shift_eta
    USE system, ONLY: V_cell_3D
    USE delta_func, ONLY: delta_gaussian
    REAL(DP), INTENT(IN) :: hw(:)
    REAL(DP), INTENT(INOUT) :: sigma_w(3, 6, SIZE(hw))
    REAL(DP), INTENT(IN) :: dE, fmn
    COMPLEX(DP), INTENT(IN) :: kernel_mn(3, 6)
    INTEGER :: iom
    REAL(DP) :: delta_w, conv_fact
    REAL(DP) :: kernel_re(3, 6)
    COMPLEX(DP), SAVE :: factor = (0.0_DP, 0.0_DP)
    LOGICAL, SAVE :: binit = .FALSE.
    !
    IF (ABS(dE) <= dE_thr) RETURN
    IF (ABS(fmn) <= 1.0D-14) RETURN
    !
    IF (.NOT. binit) THEN
      !> (e/fs * 1/V^2/fs/Ang^3) units
      factor = -zi*pi*(e_chg_au**3)/(4.0_DP*(hbar_evfs**2)*V_cell_3D)
      !> (e/fs to muA) units
      conv_fact = 1.0D6*e_chg_si/FS2SEC
      factor = factor*conv_fact
      binit = .TRUE.
    END IF
    !
    kernel_re(:, :) = DBLE(factor*kernel_mn(:, :))*fmn
    DO iom = 1, SIZE(hw)
      !> (fs) units
      delta_w = (delta_gaussian(dE - hw(iom), shift_eta) &
                 + delta_gaussian(-dE - hw(iom), shift_eta))*hbar_evfs
      sigma_w(:, :, iom) = sigma_w(:, :, iom) + kernel_re(:, :)*delta_w
    END DO
  END SUBROUTINE shift_current
  !
  SUBROUTINE dr_gen(m, n, dE_mn, dr_mn, v_k_X, A_k_X, dA_k_X)
    USE constants, ONLY: hbar_evfs, zero, zi
    USE io_input, ONLY: dE_thr
    USE system, ONLY: Nw
    INTEGER, INTENT(IN) :: m, n
    REAL(DP), INTENT(IN) :: dE_mn
    COMPLEX(DP), INTENT(OUT) :: dr_mn(3, 3)
    !< (Ang^2) units
    COMPLEX(DP), INTENT(IN) :: v_k_X(3, Nw, Nw)
    !< (Ang/fs) units
    COMPLEX(DP), INTENT(IN) :: A_k_X(3, Nw, Nw)
    !< (Ang) units
    COMPLEX(DP), INTENT(IN), OPTIONAL :: dA_k_X(3, 3)
    !< (Ang^2) units
    REAL(DP) :: pref_as
    COMPLEX(DP) :: psum, dA_X(3, 3)
    INTEGER :: ia, ib, p
    !
    ! Avoid small energy denominator
    IF (ABS(dE_mn) <= dE_thr) THEN
      dr_mn = zero
      RETURN
    END IF
    !
    IF (PRESENT(dA_k_X)) THEN
      dA_X = dA_k_X(:, :)
    ELSE
      dA_X = zero
    END IF
    !
    ! d_mn;b^a = pref_as [ r^a_mn d^b_mn + r^b_mn d^a_mn + dAX^ab + psum ]
    ! d^a_mn = v^a_mm - v^a_nn
    pref_as = -hbar_evfs/dE_mn
    DO ib = 1, 3
      DO ia = 1, 3
        ! psum = sum_p [ v^a_mp r^b_pn - r^b_mp v^a_pn ]
        psum = zero
        DO p = 1, Nw
          IF (p == m .OR. p == n) CYCLE
          psum = psum &
                 + v_k_X(ia, m, p)*A_k_X(ib, p, n) &
                 - A_k_X(ib, m, p)*v_k_X(ia, p, n)
        END DO
        dr_mn(ib, ia) = pref_as &
                        *( &
                        A_k_X(ia, m, n)*(v_k_X(ib, m, m) - v_k_X(ib, n, n)) &
                        + A_k_X(ib, m, n)*(v_k_X(ia, m, m) - v_k_X(ia, n, n)) &
                        + zi*dA_X(ia, ib)/hbar_evfs &
                        + psum &
                        )
      END DO
    END DO
  END SUBROUTINE dr_gen
  !
  SUBROUTINE shift_kernel_mn(kernel_mn, A_nm, dr_mn)
    USE constants, ONLY: zero
    COMPLEX(DP), INTENT(OUT) :: kernel_mn(3, 6)
    !< (Ang^3) units
    COMPLEX(DP), INTENT(IN) :: A_nm(3)
    !< (Ang) units
    COMPLEX(DP) :: dr_mn(3, 3)
    !< (Ang^2) units
    INTEGER, PARAMETER :: ibc2ib(6) = (/1, 1, 2, 2, 3, 3/)
    INTEGER, PARAMETER :: ibc2ic(6) = (/1, 2, 2, 3, 3, 1/)
    INTEGER :: ia, ib, ic, ibc
    COMPLEX(DP) :: rb, rc
    !
    kernel_mn = zero
    ! k^abc_mn=(r^b_nm)*(d_mn;a^c) + (r^c_nm)*(d_mn;a^b)
    DO ibc = 1, 6
      ib = ibc2ib(ibc)
      ic = ibc2ic(ibc)
      rb = A_nm(ib)
      rc = A_nm(ic)
      DO ia = 1, 3
        kernel_mn(ia, ibc) = rb*dr_mn(ia, ic) + rc*dr_mn(ia, ib)
      END DO
    END DO
  END SUBROUTINE shift_kernel_mn
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
