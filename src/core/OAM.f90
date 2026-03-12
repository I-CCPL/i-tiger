SUBROUTINE OAM_mod(eigval, v_k, L_k)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, zi
  USE io_input, ONLY: OAM_thr
  USE system, ONLY: Nw
  USE kpoints, ONLY: t_kpt, t_iks
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw, t_kpt%nkpt)
  COMPLEX(DP), INTENT(IN)::v_k(3, Nw, Nw, t_kpt%nkpt)
  COMPLEX(DP), INTENT(OUT)::L_k(3, Nw, Nw, t_kpt%nkpt)
  REAL(DP)::dE_mk, dE_nk
  INTEGER::iw, jw, kw, ipol, jpol, kpol
  !
  DO t_iks = 1, t_kpt%nkpt
    L_k(:, :, :, t_iks) = zero
    DO jw = 1, Nw
      DO kw = 1, Nw
        dE_mk = eigval(jw, t_iks) - eigval(kw, t_iks)
        IF (ABS(dE_mk) < OAM_thr) CYCLE
        DO iw = 1, Nw
          dE_nk = eigval(iw, t_iks) - eigval(kw, t_iks)
          IF (ABS(dE_nk) < OAM_thr) CYCLE
          DO kpol = 1, 3
            ipol = MOD(kpol, 3) + 1
            jpol = MOD(kpol + 1, 3) + 1
            L_k(kpol, iw, jw, t_iks) = L_k(kpol, iw, jw, t_iks) &
                                       - zi*(v_k(ipol, iw, kw, t_iks)*v_k(jpol, kw, jw, t_iks) - &
                                             v_k(jpol, iw, kw, t_iks)*v_k(ipol, kw, jw, t_iks)) &
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
  USE kpoints, ONLY: t_kpt, t_iks
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw, t_kpt%nkpt)
  COMPLEX(DP), INTENT(IN)::v_k(3, Nw, Nw, t_kpt%nkpt)
  REAL(DP), INTENT(OUT)::L_k(3, Nw, t_kpt%nkpt)
  INTEGER::iw, jw, ipol, jpol, kpol
  REAL(DP)::dE_mk
  !
  DO t_iks = 1, t_kpt%nkpt
    DO iw = 1, Nw
      L_k(:, iw, t_iks) = 0.0_DP
      DO jw = 1, Nw
        dE_mk = eigval(iw, t_iks) - eigval(jw, t_iks)
        IF (ABS(dE_mk) < OAM_thr) CYCLE
        DO kpol = 1, 3
          ipol = MOD(kpol, 3) + 1
          jpol = MOD(kpol + 1, 3) + 1
          L_k(kpol, iw, t_iks) = L_k(kpol, iw, t_iks) &
                                 + AIMAG(v_k(ipol, iw, jw, t_iks)*v_k(jpol, jw, iw, t_iks) - &
                                         v_k(jpol, iw, jw, t_iks)*v_k(ipol, jw, iw, t_iks)) &
                                 /dE_mk
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE OAM_mod_diag
