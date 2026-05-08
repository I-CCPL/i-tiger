SUBROUTINE velocity(A_bar, dH_bar, v_k_H)
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_i, hbar_eVfs
  USE system, ONLY: Nw
  USE R_vector, ONLY: R_vec_type
  USE kpoints, ONLY: t_kpt, t_iks
  USE fft_base, ONLY: fft_R2k
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN)::A_bar(Nw, Nw, 3)
  COMPLEX(DP), INTENT(IN)::dH_bar(Nw, Nw, 3)
  COMPLEX(DP), INTENT(OUT)::v_k_H(Nw, Nw, 3)
  INTEGER::iw, jw
  !
  DO jw = 1, Nw
    DO iw = 1, Nw
      v_k_H(iw, jw, :) = (dH_bar(iw, jw, :) + cmplx_i*A_bar(iw, jw, :) &
                          *(t_kpt%eigval(iw, t_iks) - t_kpt%eigval(jw, t_iks))) &
                         /hbar_eVfs
    END DO
  END DO
END SUBROUTINE velocity

SUBROUTINE vel_to_berry(eigval, v_k_H, A_k_H)
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_0, cmplx_i, hbar_eVfs
  USE f_params, ONLY: dE_thr
  USE system, ONLY: Nw
  IMPLICIT NONE
  REAL(DP), INTENT(IN)::eigval(Nw)
  COMPLEX(DP), INTENT(IN)::v_k_H(Nw, Nw, 3)
  COMPLEX(DP), INTENT(OUT)::A_k_H(Nw, Nw, 3)
  INTEGER::iw, jw
  REAL(DP)::dE
  COMPLEX(DP)::factor
  !
  factor = hbar_eVfs/cmplx_i
  DO jw = 1, Nw
    DO iw = 1, Nw
      IF (iw == jw) THEN
        A_k_H(iw, jw, :) = cmplx_0
        CYCLE
      END IF
      !
      dE = eigval(iw) - eigval(jw)
      IF (ABS(dE) <= dE_thr) THEN
        A_k_H(iw, jw, :) = cmplx_0
        CYCLE
      END IF
      !
      A_k_H(iw, jw, :) = factor*v_k_H(iw, jw, :)/dE
    END DO
  END DO
END SUBROUTINE vel_to_berry
