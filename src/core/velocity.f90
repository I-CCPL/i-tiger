SUBROUTINE velocity(R_vec, A_k_H, dH_k_H, v_k)
  USE kinds, ONLY: DP
  USE constants, ONLY: zi, hbar_evfs
  USE system, ONLY: Nw
  USE R_vector, ONLY: R_vec_type
  USE kpoints, ONLY: t_kpt, t_iks
  USE fft_base, ONLY: fft_R2k
  IMPLICIT NONE
  TYPE(R_vec_type), INTENT(INOUT)::R_vec
  COMPLEX(DP), INTENT(IN)::A_k_H(3, Nw, Nw)
  COMPLEX(DP), INTENT(IN)::dH_k_H(3, Nw, Nw)
  COMPLEX(DP), INTENT(OUT)::v_k(3, Nw, Nw)
  INTEGER::iw, jw
  !
  DO jw = 1, Nw
    DO iw = 1, Nw
      v_k(:, iw, jw) = (dH_k_H(:, iw, jw) + zi*A_k_H(:, iw, jw) &
                        *(t_kpt%eigval(iw, t_iks) - t_kpt%eigval(jw, t_iks))) &
                       /hbar_evfs
    END DO
  END DO
END SUBROUTINE velocity
