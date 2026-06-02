SUBROUTINE get_berry_nk_g(eigval, v_k, O_k)
  USE kinds, ONLY: DP
  USE constants, ONLY: hbar_eVfs
  USE system, ONLY: Nw
  USE f_params, ONLY: dE_thr
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  COMPLEX(DP), INTENT(IN)::v_k(Nw, Nw, 3)
  REAL(DP), INTENT(OUT)::O_k(Nw, 3)
  INTEGER::iw1, iw2, ipol, jpol, kpol
  REAL(DP)::denom, factor
  ! v_k [Ang/fs], denom [eV] => O_k [Ang^2]
  factor = -2.0_DP*(hbar_eVfs**2)
  DO iw1 = 1, Nw
    O_k(iw1, :) = 0.0_DP
    DO iw2 = 1, Nw
      denom = eigval(iw1) - eigval(iw2)
      IF (ABS(denom) <= dE_thr) CYCLE
      DO kpol = 1, 3
        ipol = MOD(kpol, 3) + 1
        jpol = MOD(kpol + 1, 3) + 1
        O_k(iw1, kpol) = O_k(iw1, kpol) + factor*AIMAG(v_k(iw1, iw2, ipol)*v_k(iw2, iw1, jpol)/(denom**2))
      END DO
    END DO
  END DO
END SUBROUTINE get_berry_nk_g

SUBROUTINE get_berry_sum_g(eigval, O_k, berry)
  USE kinds, ONLY: DP
  USE f_params, ONLY: E_fermi
  USE system, ONLY: Nw
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  REAL(DP), INTENT(IN)::O_k(Nw, 3)
  REAL(DP), INTENT(OUT)::berry(3)
  INTEGER::iw, ipol
  berry = 0.0_DP
  DO iw = 1, Nw
    IF (eigval(iw) > E_fermi) CYCLE
    berry(:) = berry(:) + O_k(iw, :)
  END DO
END SUBROUTINE get_berry_sum_g

SUBROUTINE get_berry_k_g(eigval, v_k, berry)
  USE kinds, ONLY: DP
  USE constants, ONLY: hbar_eVfs
  USE system, ONLY: Nw
  USE f_params, ONLY: E_fermi, dE_thr
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  COMPLEX(DP), INTENT(IN)::v_k(Nw, Nw, 3)
  REAL(DP), INTENT(OUT)::berry(3)
  INTEGER::iw1, iw2, ipol, jpol, kpol
  REAL(DP)::denom, factor
  ! v_k [Ang/fs], denom [eV] => O_k [Ang^2]
  factor = -2.0_DP*(hbar_eVfs**2)
  berry = 0.0_DP
  DO iw1 = 1, Nw
    IF (eigval(iw1) > E_fermi) CYCLE
    DO iw2 = 1, Nw
      denom = eigval(iw1) - eigval(iw2)
      IF (ABS(denom) <= dE_thr) CYCLE
      DO kpol = 1, 3
        ipol = MOD(kpol, 3) + 1
        jpol = MOD(kpol + 1, 3) + 1
        berry(kpol) = berry(kpol) + factor*AIMAG(v_k(iw1, iw2, ipol)*v_k(iw2, iw1, jpol)/(denom**2))
      END DO
    END DO
  END DO
END SUBROUTINE get_berry_k_g
