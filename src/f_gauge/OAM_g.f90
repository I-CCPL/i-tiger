SUBROUTINE OAM_g_mod(eigval, v_k, L_k)
  !< OAM [hbar]
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_0, cmplx_i, m_e
  USE f_params, ONLY: dE_thr
  USE system, ONLY: Nw
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  COMPLEX(DP), INTENT(IN)::v_k(Nw, Nw, 3)
  COMPLEX(DP), INTENT(OUT)::L_k(Nw, Nw, 3)
  REAL(DP)::dE_mk, dE_nk, factor
  INTEGER::iw, jw, kw, ipol, jpol, kpol
  !
  factor = -cmplx_i*m_e/2
  DO jw = 1, Nw
    L_k(:, jw, :) = cmplx_0
    DO kw = 1, Nw
      dE_mk = eigval(jw) - eigval(kw)
      IF (ABS(dE_mk) < dE_thr) CYCLE
      DO iw = 1, Nw
        dE_nk = eigval(iw) - eigval(kw)
        IF (ABS(dE_nk) < dE_thr) CYCLE
        DO kpol = 1, 3
          ipol = MOD(kpol, 3) + 1
          jpol = MOD(kpol + 1, 3) + 1
          L_k(iw, jw, kpol) = L_k(iw, jw, kpol) &
                              + factor*(v_k(iw, kw, ipol)*v_k(kw, jw, jpol) - &
                                        v_k(iw, kw, jpol)*v_k(kw, jw, ipol)) &
                              *(1/dE_nk + 1/dE_mk)
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE OAM_g_mod

SUBROUTINE OAM_g_mod_diag(eigval, v_k, L_k)
  !< OAM [hbar] for diagonal elements
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_0, cmplx_i, m_e
  USE f_params, ONLY: dE_thr
  USE system, ONLY: Nw
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  COMPLEX(DP), INTENT(IN)::v_k(Nw, Nw, 3)
  REAL(DP), INTENT(OUT)::L_k(Nw, 3)
  INTEGER::iw, kw, ipol, jpol, kpol
  REAL(DP)::factor, dE_mk
  !
  factor = m_e
  DO iw = 1, Nw
    L_k(iw, :) = 0.0_DP
    DO kw = 1, Nw
      dE_mk = eigval(iw) - eigval(kw)
      IF (ABS(dE_mk) < dE_thr) CYCLE
      DO kpol = 1, 3
        ipol = MOD(kpol, 3) + 1
        jpol = MOD(kpol + 1, 3) + 1
        L_k(iw, kpol) = L_k(iw, kpol) &
                        + factor*AIMAG(v_k(iw, kw, ipol)*v_k(kw, iw, jpol) - &
                                       v_k(iw, kw, jpol)*v_k(kw, iw, ipol)) &
                        /dE_mk
      END DO
    END DO
  END DO
END SUBROUTINE OAM_g_mod_diag
