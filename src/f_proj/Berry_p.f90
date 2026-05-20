SUBROUTINE get_berry_p_nk(l, eigval, O_bar, A_bar, dH_bar, berry)
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_0, cmplx_1, cmplx_i
  USE f_params, ONLY: dE_thr, E_fermi
  USE system, ONLY: Nw
  USE kpoints, ONLY: t_kpt
  USE itg_k, ONLY: k_data
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l
  REAL(DP), INTENT(IN) :: eigval(Nw)
  COMPLEX(DP), INTENT(IN) :: O_bar(Nw, Nw, 3)
  COMPLEX(DP), INTENT(IN) :: A_bar(Nw, Nw, 3)
  COMPLEX(DP), INTENT(IN) :: dH_bar(Nw, Nw, 3)
  INTEGER::n, m, p, ief, a, b, c
  REAL(DP), INTENT(OUT) :: berry(3)
  REAL(DP)::occ(Nw)
  COMPLEX(DP)::f_list(Nw, Nw), g_list(Nw, Nw)
  COMPLEX(DP)::JJm(Nw, Nw, 3), JJp(Nw, Nw, 3)
  REAL(DP)::J0, J1, J1_1, J1_2, J2, J2_1
  !
  occ = 0.0_DP
  occ(l) = 1.0_DP
  CALL compute_JJ_list(occ, JJm, JJp)
  CALL compute_occ_mat(occ, f_list, g_list)

  DO c = 1, 3
    a = MOD(c, 3) + 1
    b = MOD(a, 3) + 1
    J0 = 0.0_DP
    J1 = 0.0_DP
    J2 = 0.0_DP
    DO n = 1, Nw
      DO m = 1, Nw
        J0 = J0 + DBLE(f_list(m, n)*k_data%mO_k_W(n, m, c))
        J1 = J1 - 2.0_DP*(AIMAG(k_data%mA_k_W(n, m, a)*JJp(m, n, b)) &
                          + AIMAG(k_data%mA_k_W(m, n, b)*JJm(n, m, a)))
        J2 = J2 - 2.0_DP*AIMAG(JJm(m, n, a)*JJp(n, m, b))
      END DO
    END DO
    berry(c) = J0 + J1 + J2
  END DO
END SUBROUTINE get_berry_p_nk

SUBROUTINE get_berry_p_k(eigval, O_bar, A_bar, dH_bar, berry)
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
END SUBROUTINE get_berry_p_k
