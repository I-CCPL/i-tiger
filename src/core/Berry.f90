SUBROUTINE Berry_mod(eigval, v_k, O_k)
  USE constants, ONLY: DP, hbar_eVfs
  USE system, ONLY: Nw
  USE io_input, ONLY: dE_thr
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  COMPLEX(DP), INTENT(IN)::v_k(3, Nw, Nw)
  REAL(DP), INTENT(OUT)::O_k(3, Nw)
  INTEGER::iw1, iw2, ipol, jpol, kpol
  REAL(DP)::denom, factor
  factor = -2.0_DP*(hbar_eVfs**2)
  DO iw1 = 1, Nw
    O_k(:, iw1) = 0.0_DP
    DO iw2 = 1, Nw
      denom = eigval(iw1) - eigval(iw2)
      IF (ABS(denom) <= dE_thr) CYCLE
      DO kpol = 1, 3
        ipol = MOD(kpol, 3) + 1
        jpol = MOD(kpol + 1, 3) + 1
        O_k(kpol, iw1) = O_k(kpol, iw1) + factor*DIMAG(v_k(ipol, iw1, iw2)*v_k(jpol, iw2, iw1)/(denom**2))
      END DO
    END DO
  END DO
END SUBROUTINE Berry_mod

SUBROUTINE Berry_sum(eigval, O_k, berry)
  USE constants, ONLY: DP
  USE io_input, ONLY: Ef
  USE system, ONLY: Nw
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  REAL(DP), INTENT(IN)::O_k(3, Nw)
  REAL(DP), INTENT(OUT)::berry(3)
  INTEGER::iw, ipol
  berry = 0.0_DP
  DO iw = 1, Nw
    IF (eigval(iw) > Ef) CYCLE
    DO ipol = 1, 3
      berry(ipol) = berry(ipol) + O_k(ipol, iw)
    END DO
  END DO
END SUBROUTINE Berry_sum
