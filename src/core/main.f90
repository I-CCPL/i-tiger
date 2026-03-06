PROGRAM main
  USE kinds, ONLY: DP
  USE env, ONLY: env_start, env_end
  USE io_input, ONLY: read_input
  USE lin_eig_H, ONLY: write_band
  USE system, ONLY: Nw
  USE kpoints, ONLY: t_kpt
  USE wannier90, ONLY: w90data
  USE R_vector, ONLY: R_vec_type
  IMPLICIT NONE
  TYPE(R_vec_type)::R_vec
  COMPLEX(DP), ALLOCATABLE::H_R(:, :, :), H_k(:, :, :)
  REAL(DP), ALLOCATABLE::eigval(:, :)
  COMPLEX(DP), ALLOCATABLE::eigvec(:, :, :)
  !
  CALL env_start()
  CALL read_input()
  CALL w90data%build_Hq()
  CALL w90data%build_bvec()
  ! CALL w90data%build_Aq()
  !
  CALL R_vec%build_R(w90data)
  CALL R_vec%fft_q2R(w90data, w90data%Hq, H_R)
  ! TODO: k parallelization
  ALLOCATE (H_k(Nw, Nw, t_kpt%nkpt))
  CALL R_vec%fft_R2k(H_R, H_k)
  !
  ! TODO: write_band from input
  ALLOCATE (eigval(Nw, t_kpt%nkpt))
  ALLOCATE (eigvec(Nw, Nw, t_kpt%nkpt))
  CALL write_band(H_k, eigval, eigvec)
  !
  CALL w90data%clear()
  CALL R_vec%clear()
  CALL env_end()
END PROGRAM main
