PROGRAM main
  USE env, ONLY: env_start, env_end, print_k_info
  USE io_global, ONLY: stdout
  USE io_input, ONLY: read_input
  USE kpoints, ONLY: t_kpt, t_iks
  USE itg_q, ONLY: make_q, clear_q, bcast_q
  USE itg_R, ONLY: make_R, clear_R, bcast_R
  USE itg_k, ONLY: allocate_k, make_k, clear_k
  USE itg_f, ONLY: allocate_f, clear_f, make_f, write_f
  IMPLICIT NONE
  !
  CALL env_start(__DATE__, __TIME__)
  CALL itg_warnings()
  CALL read_input()
  !
  CALL make_q()
  CALL debug_q()
  ! CALL bcast_q()
  !
  CALL make_R()
  CALL debug_R()
  CALL clear_q()
  CALL bcast_R()
  !
  WRITE (stdout, '(2X, A)') 'Building data in k space...'
  CALL t_kpt%divide_k()
  CALL allocate_k()
  CALL allocate_f()
  DO t_iks = 1, t_kpt%nkpt
    CALL print_k_info(t_iks, t_kpt%nkpt)
    CALL make_k()
    CALL make_f()
    CALL debug_k()
  END DO
  CALL print_k_info(t_iks, t_kpt%nkpt)
  CALL write_sep_line()

  CALL clear_R()
  CALL write_f()
  !
  CALL clear_k()
  CALL clear_f()
  CALL env_end()
END PROGRAM main
