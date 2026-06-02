SUBROUTINE compute_D_bar(dH_bar, eigval, D_bar)
  USE kinds, ONLY: DP
  USE f_params, ONLY: dE_thr
  USE system, ONLY: Nw
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN)::dH_bar(Nw, Nw, 3)
  REAL(DP), INTENT(IN)::eigval(Nw)
  COMPLEX(DP), INTENT(OUT)::D_bar(Nw, Nw, 3)
  INTEGER::m, n
  ! dH_bar [eV*Ang], eigval [eV] => D_bar [Ang]
  DO m = 1, Nw
    DO n = 1, Nw
      IF (ABS(eigval(m) - eigval(n)) <= dE_thr) THEN
        D_bar(m, n, :) = 0.0_DP
      ELSE
        D_bar(m, n, :) = dH_bar(m, n, :)/(eigval(n) - eigval(m))
      END IF
    END DO
  END DO
END SUBROUTINE compute_D_bar

SUBROUTINE compute_v_k_H(A_bar, dH_bar, v_k_H)
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_i, hbar_eVfs
  USE system, ONLY: Nw
  USE kpoints, ONLY: t_kpt, t_iks
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN)::A_bar(Nw, Nw, 3)
  COMPLEX(DP), INTENT(IN)::dH_bar(Nw, Nw, 3)
  COMPLEX(DP), INTENT(OUT)::v_k_H(Nw, Nw, 3)
  INTEGER::iw, jw
  ! dH_bar [eV*Ang], A_bar [Ang], eigval [eV] => v_k_H [Ang/fs]
  DO jw = 1, Nw
    DO iw = 1, Nw
      v_k_H(iw, jw, :) = (dH_bar(iw, jw, :) + cmplx_i*A_bar(iw, jw, :) &
                          *(t_kpt%eigval(iw, t_iks) - t_kpt%eigval(jw, t_iks))) &
                         /hbar_eVfs
    END DO
  END DO
END SUBROUTINE compute_v_k_H

SUBROUTINE compute_occ_mat(occ, f_list, g_list)
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_0, cmplx_1
  USE system, ONLY: Nw
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::occ(Nw)
  COMPLEX(DP), INTENT(OUT)::f_list(Nw, Nw), g_list(Nw, Nw)
  INTEGER::m, n, p
  COMPLEX(DP)::tmp(Nw, Nw)
  !> TODO: Test performance
  f_list = cmplx_0
  DO n = 1, Nw
    DO m = 1, Nw
      DO p = 1, Nw
        f_list(n, m) = f_list(n, m) &
                       + t_kpt%eigvec(n, p)*occ(p)*CONJG(t_kpt%eigvec(m, p))
      END DO
      g_list(n, m) = -f_list(n, m)
      IF (m == n) g_list(n, n) = g_list(n, n) + cmplx_1
    END DO
  END DO

  ! ! occ * U
  ! DO p = 1, Nw
  !   tmp(:, p) = t_kpt%eigvec(:, p)*occ(p)
  ! END DO

  ! ! U^dag * occ * U
  ! CALL ZGEMM('N', 'C', Nw, Nw, Nw, cmplx_1, &
  !            tmp, Nw, t_kpt%eigvec, Nw, cmplx_0, f_list, Nw)

  ! ! g = 1-f
  ! g_list = -f_list
  ! DO p = 1, Nw
  !   g_list(p, p) = g_list(p, p) + cmplx_1
  ! END DO
END SUBROUTINE compute_occ_mat

SUBROUTINE compute_JJ_list(occ, JJm, JJp)
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_0, cmplx_i
  USE f_params, ONLY: E_fermi
  USE system, ONLY: Nw
  USE itg_k, ONLY: k_data, rotate_H2W
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::occ(Nw)
  COMPLEX(DP), INTENT(OUT)::JJm(Nw, Nw, 3), JJp(Nw, Nw, 3)
  INTEGER::m, n, a

  DO m = 1, Nw
    DO n = 1, Nw
      IF (occ(m) < 0.5_DP .AND. occ(n) > 0.5_DP) THEN
        JJm(n, m, :) = cmplx_i*k_data%mD_bar(n, m, :)
        JJp(m, n, :) = cmplx_i*k_data%mD_bar(m, n, :)
      ELSE
        JJm(n, m, :) = cmplx_0
        JJp(m, n, :) = cmplx_0
      END IF
    END DO
  END DO
  DO a = 1, 3
    CALL rotate_H2W(JJm(:, :, a))
    CALL rotate_H2W(JJp(:, :, a))
  END DO
END SUBROUTINE compute_JJ_list
