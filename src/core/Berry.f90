SUBROUTINE Berry_mod(eigval, v_k, O_k)
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
  factor = -2.0_DP*(hbar_eVfs**2)
  DO iw1 = 1, Nw
    O_k(iw1, :) = 0.0_DP
    DO iw2 = 1, Nw
      denom = eigval(iw1) - eigval(iw2)
      IF (ABS(denom) <= dE_thr) CYCLE
      DO kpol = 1, 3
        ipol = MOD(kpol, 3) + 1
        jpol = MOD(kpol + 1, 3) + 1
        O_k(iw1, kpol) = O_k(iw1, kpol) + factor*DIMAG(v_k(iw1, iw2, ipol)*v_k(iw2, iw1, jpol)/(denom**2))
      END DO
    END DO
  END DO
END SUBROUTINE Berry_mod

SUBROUTINE Berry_sum(eigval, O_k, berry)
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
END SUBROUTINE Berry_sum

SUBROUTINE Berry_proj(eigval, O_bar, A_bar, dH_bar, berry)
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_0, cmplx_i
  USE f_params, ONLY: dE_thr, E_fermi
  USE system, ONLY: Nw
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  COMPLEX(DP), INTENT(IN)::O_bar(Nw, Nw, 3)
  COMPLEX(DP), INTENT(IN)::A_bar(Nw, Nw, 3)
  COMPLEX(DP), INTENT(IN)::dH_bar(Nw, Nw, 3)
  REAL(DP), INTENT(OUT)::berry(3)
  COMPLEX(DP):: D_val(Nw, Nw, 3), sum_val(3)
  INTEGER::m, n, a, b, c
  DO m = 1, Nw
    DO n = 1, Nw
      IF (ABS(eigval(m) - eigval(n)) <= dE_thr) THEN
        D_val(m, n, :) = 0.0_DP
      ELSE
        D_val(m, n, :) = dH_bar(m, n, :)/(eigval(n) - eigval(m))
      END IF
    END DO
  END DO

  berry = 0.0_DP
  DO n = 1, Nw
    IF (eigval(n) > E_fermi) CYCLE
    sum_val = cmplx_0
    DO m = 1, Nw
      IF (eigval(m) < E_fermi) CYCLE
      DO c = 1, 3
        a = MOD(c, 3) + 1
        b = MOD(a, 3) + 1
        sum_val(c) = sum_val(c) &
                     + D_val(n, m, a)*A_bar(m, n, b) &
                     - D_val(n, m, b)*A_bar(m, n, a) &
                     + cmplx_i*D_val(n, m, a)*D_val(m, n, b)
      END DO
    END DO
    berry = berry + REAL(O_bar(n, n, :) - 2*sum_val, DP)
  END DO
END SUBROUTINE Berry_proj
