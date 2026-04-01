PROGRAM main
  USE env, ONLY: env_start, env_end, print_k_info
  USE io_input, ONLY: read_input
  USE kpoints, ONLY: t_kpt, t_iks
  USE itg_q, ONLY: make_q_data, clear_q_data
  USE itg_R, ONLY: make_R_data, clear_R_data
  USE itg_k, ONLY: allocate_k_data, make_k_data, &
                   clear_k_data, write_k_data
  IMPLICIT NONE
  !
  CALL env_start(__DATE__, __TIME__)
  CALL itg_warnings()
  CALL read_input()
  !
  CALL make_q_data()
  CALL make_R_data()
  CALL clear_q_data()
  !
  CALL allocate_k_data()
  DO t_iks = 1, t_kpt%nkpt
    CALL print_k_info(t_iks, t_kpt%nkpt)
    CALL make_k_data()
  END DO
  CALL print_k_info(t_iks, t_kpt%nkpt)

  CALL clear_R_data()
  CALL write_k_data()
  !
  CALL clear_k_data()
  CALL env_end()
END PROGRAM main
