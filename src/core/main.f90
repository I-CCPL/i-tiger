PROGRAM main
  USE env, ONLY: env_start, env_end
  USE io_input, ONLY: read_input
  USE itg_q, ONLY: make_q, clear_q_data
  USE itg_R, ONLY: make_R, clear_R_data
  USE itg_k, ONLY: make_k, clear_k_data
  IMPLICIT NONE
  !
  CALL env_start(__DATE__, __TIME__)
  CALL itg_warnings()
  CALL read_input()
  !
  CALL make_q()
  CALL make_R()
  CALL make_k()
  !
  CALL clear_q_data()
  CALL clear_R_data()
  CALL clear_k_data()
  CALL env_end()
END PROGRAM main
