SUBROUTINE OAM_mod(eigval, v_k, L_k)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, zi
  USE io_input, ONLY: OAM_thr
  USE system, ONLY: Nw
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw, t_kpt%nkpt)
  COMPLEX(DP), INTENT(IN)::v_k(3, Nw, Nw, t_kpt%nkpt)
  COMPLEX(DP), INTENT(OUT)::L_k(3, Nw, Nw, t_kpt%nkpt)
  REAL(DP)::dE_mk, dE_nk
  INTEGER::ikpt, iw, jw, kw, ipol, jpol, kpol
  !
  DO ikpt = 1, t_kpt%nkpt
    L_k(:, :, :, ikpt) = zero
    DO jw = 1, Nw
      DO kw = 1, Nw
        dE_mk = eigval(jw, ikpt) - eigval(kw, ikpt)
        IF (ABS(dE_mk) < OAM_thr) CYCLE
        DO iw = 1, Nw
          dE_nk = eigval(iw, ikpt) - eigval(kw, ikpt)
          IF (ABS(dE_nk) < OAM_thr) CYCLE
          DO kpol = 1, 3
            ipol = MOD(kpol, 3) + 1
            jpol = MOD(kpol + 1, 3) + 1
            L_k(kpol, iw, jw, ikpt) = L_k(kpol, iw, jw, ikpt) &
                                      - zi*(v_k(ipol, iw, kw, ikpt)*v_k(jpol, kw, jw, ikpt) - &
                                            v_k(jpol, iw, kw, ikpt)*v_k(ipol, kw, jw, ikpt)) &
                                      *(1/dE_nk + 1/dE_mk)/2
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE OAM_mod

SUBROUTINE OAM_mod_diag(eigval, v_k, L_k)
  USE kinds, ONLY: DP
  USE io_input, ONLY: OAM_thr
  USE system, ONLY: Nw
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw, t_kpt%nkpt)
  COMPLEX(DP), INTENT(IN)::v_k(3, Nw, Nw, t_kpt%nkpt)
  REAL(DP), INTENT(OUT)::L_k(3, Nw, t_kpt%nkpt)
  INTEGER::ikpt, iw, jw, ipol, jpol, kpol
  REAL(DP)::dE_mk
  !
  DO ikpt = 1, t_kpt%nkpt
    DO iw = 1, Nw
      L_k(:, iw, ikpt) = 0.0_DP
      DO jw = 1, Nw
        dE_mk = eigval(iw, ikpt) - eigval(jw, ikpt)
        IF (ABS(dE_mk) < OAM_thr) CYCLE
        DO kpol = 1, 3
          ipol = MOD(kpol, 3) + 1
          jpol = MOD(kpol + 1, 3) + 1
          L_k(kpol, iw, ikpt) = L_k(kpol, iw, ikpt) &
                                + AIMAG(v_k(ipol, iw, jw, ikpt)*v_k(jpol, jw, iw, ikpt) - &
                                        v_k(jpol, iw, jw, ikpt)*v_k(ipol, jw, iw, ikpt)) &
                                /dE_mk
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE OAM_mod_diag
