SUBROUTINE velocity(R_vec, H_R, A_k_H, v_k)
  USE kinds, ONLY: DP
  USE constants, ONLY: zi, zero
  USE system, ONLY: Nw
  USE debug_data, ONLY: write_matrix
  USE R_vector, ONLY: R_vec_type
  USE kpoints, ONLY: t_kpt, kpoint_type
  USE fft_base, ONLY: fft_R2k
  USE der_base, ONLY: der_R
  IMPLICIT NONE
  TYPE(R_vec_type), INTENT(INOUT)::R_vec
  COMPLEX(DP), INTENT(IN)::H_R(Nw, Nw, R_vec%nrpt)
  COMPLEX(DP), INTENT(IN)::A_k_H(3, Nw, Nw, t_kpt%nkpt)
  COMPLEX(DP), INTENT(OUT)::v_k(3, Nw, Nw, t_kpt%nkpt)
  COMPLEX(DP)::dH_R(3, Nw, Nw, R_vec%nrpt)
  COMPLEX(DP)::dH_k_W(3, Nw, Nw, t_kpt%nkpt)
  COMPLEX(DP)::dH_k_H(3, Nw, Nw, t_kpt%nkpt)
  INTEGER::ib, irpt, ikpt, iw, jw, ipol
  !
  CALL der_R(R_vec, H_R, dH_R)
  ! CALL write_matrix('dH_R.itg', dH_R, R_vec%nrpt, 1)
  CALL fft_R2k(R_vec, dH_R, dH_k_W)
  ! CALL write_matrix('dH_k_W.itg', dH_k_W, t_kpt%nkpt, 1)
  CALL t_kpt%rotate(dH_k_W, dH_k_H)
  ! CALL write_matrix('dH_k_H.itg', dH_k_H, t_kpt%nkpt, 1)
  DO ikpt = 1, t_kpt%nkpt
    DO iw = 1, Nw
      DO jw = 1, Nw
        v_k(:, iw, jw, ikpt) = dH_k_H(:, iw, jw, ikpt) + zi*A_k_H(:, iw, jw, ikpt) &
                               *(t_kpt%eigval(iw, ikpt) - t_kpt%eigval(jw, ikpt))
      END DO
    END DO
  END DO
  CALL check_Hermiticity(t_kpt%nkpt, 3, v_k, 1.0E-10_DP)
  ! CALL write_matrix('v_k.itg', v_k, t_kpt%nkpt, 1)
END SUBROUTINE velocity
