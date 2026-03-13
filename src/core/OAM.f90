SUBROUTINE OAM_mod(eigval, v_k, L_k)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, zi
  USE io_input, ONLY: OAM_thr
  USE system, ONLY: Nw
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  COMPLEX(DP), INTENT(IN)::v_k(3, Nw, Nw)
  COMPLEX(DP), INTENT(OUT)::L_k(3, Nw, Nw)
  REAL(DP)::dE_mk, dE_nk
  INTEGER::iw, jw, kw, ipol, jpol, kpol
  !
  DO jw = 1, Nw
    L_k(:, :, jw) = zero
    DO kw = 1, Nw
      dE_mk = eigval(jw) - eigval(kw)
      IF (ABS(dE_mk) < OAM_thr) CYCLE
      DO iw = 1, Nw
        dE_nk = eigval(iw) - eigval(kw)
        IF (ABS(dE_nk) < OAM_thr) CYCLE
        DO kpol = 1, 3
          ipol = MOD(kpol, 3) + 1
          jpol = MOD(kpol + 1, 3) + 1
          L_k(kpol, iw, jw) = L_k(kpol, iw, jw) &
                              - zi*(v_k(ipol, iw, kw)*v_k(jpol, kw, jw) - &
                                    v_k(jpol, iw, kw)*v_k(ipol, kw, jw)) &
                              *(1/dE_nk + 1/dE_mk)/2
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE OAM_mod

SUBROUTINE OAM_mod_diag(eigval, v_k, L_k)
  USE kinds, ONLY: DP
  USE io_input, ONLY: OAM_thr
  USE system, ONLY: Nw
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  COMPLEX(DP), INTENT(IN)::v_k(3, Nw, Nw)
  REAL(DP), INTENT(OUT)::L_k(3, Nw)
  INTEGER::iw, jw, ipol, jpol, kpol
  REAL(DP)::dE_mk
  !
  DO iw = 1, Nw
    L_k(:, iw) = 0.0_DP
    DO jw = 1, Nw
      dE_mk = eigval(iw) - eigval(jw)
      IF (ABS(dE_mk) < OAM_thr) CYCLE
      DO kpol = 1, 3
        ipol = MOD(kpol, 3) + 1
        jpol = MOD(kpol + 1, 3) + 1
        L_k(kpol, iw) = L_k(kpol, iw) &
                        + AIMAG(v_k(ipol, iw, jw)*v_k(jpol, jw, iw) - &
                                v_k(jpol, iw, jw)*v_k(ipol, jw, iw)) &
                        /dE_mk
      END DO
    END DO
  END DO
END SUBROUTINE OAM_mod_diag
