PROGRAM main
  USE kinds, ONLY: DP
  USE env, ONLY: env_start, env_end
  USE io_input, ONLY: read_input
  USE lin_eig_H, ONLY: write_band
  USE system, ONLY: Nw
  USE kpoints, ONLY: kpoint_type
  USE wannier90, ONLY: w90data
  USE R_vector, ONLY: R_vec_type
  IMPLICIT NONE
  TYPE(R_vec_type)::R_vec
  COMPLEX(DP), ALLOCATABLE::H_R(:, :, :), H_k(:, :, :)
  REAL(DP), ALLOCATABLE::eigval(:, :)
  COMPLEX(DP), ALLOCATABLE::eigvec(:, :, :)
  TYPE(kpoint_type)::k_list
  !
  CALL env_start()
  CALL read_input()
  CALL w90data%build_Hq()
  !
  CALL R_vec%build_R(w90data)
  CALL R_vec%fft_q2R(w90data, w90data%Hq, H_R)
  !
  ! TODO: k-point from input
  CALL k_list%build_path(3, &
                         (/(/0.0_DP, 0.0_DP, -0.5_DP/), &
                           (/0.0_DP, 0.0_DP, 0.0_DP/), &
                           (/0.0_DP, 0.0_DP, 0.5_DP/)/), &
                         (/50, 50, 50/))
  ALLOCATE (H_k(Nw, Nw, k_list%nkpt))
  CALL R_vec%fft_R2k(k_list, H_R, H_k)
  !
  ! TODO: write_band from input
  ALLOCATE (eigval(Nw, k_list%nkpt))
  ALLOCATE (eigvec(Nw, Nw, k_list%nkpt))
  CALL write_band(H_k, k_list%nkpt, eigval, eigvec)
  !
  CALL w90data%clear()
  CALL R_vec%clear()
  CALL env_end()
END PROGRAM main
