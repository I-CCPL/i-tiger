PROGRAM main
  USE kinds, ONLY: DP, eq_vec_real
  USE constants, ONLY: zero
  USE env, ONLY: env_start, env_end
  USE io_input, ONLY: read_input
  USE lin_eig_H, ONLY: write_band, eig_H
  USE fft_base, ONLY: fft_q2R, fft_R2k
  USE system, ONLY: Nw
  USE debug_data, ONLY: write_matrix, write_diag_matrix
  USE kpoints, ONLY: t_kpt, kpoint_type
  USE wannier90, ONLY: w90data
  USE R_vector, ONLY: R_vec_type
  IMPLICIT NONE
  TYPE(R_vec_type)::R_vec
  COMPLEX(DP), ALLOCATABLE::H_R(:, :, :), H_k_H(:, :, :)
  COMPLEX(DP), ALLOCATABLE::A_R(:, :, :, :), A_R_b(:, :, :, :)
  COMPLEX(DP), ALLOCATABLE::A_k_W(:, :, :, :), A_k_H(:, :, :, :)
  COMPLEX(DP), ALLOCATABLE::v_k(:, :, :, :)
  REAL(DP), ALLOCATABLE::L_k(:, :, :)
  INTEGER::inb, ikpt, irpt, jrpt, iw, jw
  !
  CALL env_start(__DATE__, __TIME__)
  CALL itg_warnings()
  CALL read_input()

  CALL w90data%build_Hq()
  ! CALL write_matrix('H_q.itg', w90data%Hq, w90data%kpts%nkpt, 1)
  CALL w90data%build_bvec()
  CALL w90data%build_Aq()
  !
  CALL R_vec%build_R(w90data)
  ALLOCATE (H_R(Nw, Nw, R_vec%nRpt))
  CALL fft_q2R(w90data, R_vec, w90data%Hq, H_R)
  ! CALL write_matrix('H_R.itg', H_R, R_vec%nRpt, 1)

  ALLOCATE (A_R(3, Nw, Nw, R_vec%nRpt))
  ALLOCATE (A_R_b(3, Nw, Nw, R_vec%nRpt))
  A_R = zero
  DO inb = 1, w90data%nnb
    CALL fft_q2R(w90data, R_vec, w90data%Aq(:, :, :, :, inb), A_R_b)
    ikpt = 0
    DO irpt = 1, R_vec%nRpt
      DO jrpt = 1, R_vec%nRpt
        IF (.NOT. eq_vec_real(R_Vec%R_red(:, irpt), &
                              -R_vec%R_red(:, jrpt), 1.0D-12)) THEN
          CYCLE
        END IF
        ikpt = ikpt + 1
        DO jw = 1, Nw
          DO iw = 1, Nw
            A_R(:, iw, jw, irpt) = A_R(:, iw, jw, irpt) + (A_R_b(:, iw, jw, irpt) + CONJG(A_R_b(:, jw, iw, jrpt)))/2
          END DO
        END DO
      END DO
    END DO
  END DO
  ! CALL write_matrix('A_R.itg', A_R, R_vec%nRpt, 1)

  ! TODO: k parallelization
  ALLOCATE (t_kpt%H_k(Nw, Nw, t_kpt%nkpt))
  CALL fft_R2k(R_vec, H_R, t_kpt%H_k)

  DO ikpt = 1, t_kpt%nkpt
    DO iw = 1, Nw
      t_kpt%H_k(iw, iw, ikpt) = REAL(t_kpt%H_k(iw, iw, ikpt), DP)
    END DO
  END DO
  ! CALL write_matrix('H_k_W.itg', t_kpt%H_k, t_kpt%nkpt, 1)

  ALLOCATE (t_kpt%eigval(Nw, t_kpt%nkpt))
  ALLOCATE (t_kpt%eigvec(Nw, Nw, t_kpt%nkpt))
  CALL eig_H(Nw, t_kpt%nkpt, t_kpt%H_k(:, :, :), t_kpt%eigval, t_kpt%eigvec)

  ! CALL write_matrix('eigvec.itg', t_kpt%eigvec, t_kpt%nkpt, 1)

  ALLOCATE (H_k_H(Nw, Nw, t_kpt%nkpt))
  CALL t_kpt%rotate(t_kpt%H_k, H_k_H)
  ! CALL write_matrix('H_k_H.itg', H_k_H, t_kpt%nkpt, 1)

  ALLOCATE (A_k_W(3, Nw, Nw, t_kpt%nkpt))
  ALLOCATE (A_k_H(3, Nw, Nw, t_kpt%nkpt))
  CALL fft_R2k(R_vec, A_R, A_k_W)
  ! CALL write_matrix('A_k_W.itg', A_k_W, t_kpt%nkpt, 1)

  CALL t_kpt%rotate(A_k_W, A_k_H)
  ! CALL write_matrix('A_k_H.itg', A_k_H, t_kpt%nkpt, 1)
  !

  ALLOCATE (v_k(3, Nw, Nw, t_kpt%nkpt))
  CALL velocity(R_vec, H_R, A_k_H, v_k)

  ALLOCATE (L_k(3, Nw, t_kpt%nkpt))
  CALL OAM_mod_diag(t_kpt%eigval, v_k, L_k)
  CALL write_diag_matrix('OAM_diag.itg', L_k, t_kpt%nkpt, 1)
  !
  ! TODO: write_band from input
  CALL write_band(t_kpt%H_k, t_kpt%eigval)
  !
  CALL w90data%clear()
  CALL R_vec%clear()
  CALL env_end()
END PROGRAM main
