PROGRAM main
  USE env, ONLY: env_start, env_end
  USE io_input, ONLY: read_input
  IMPLICIT NONE
  !
  CALL env_start()
  CALL read_input()
  !
  CALL env_end()
END PROGRAM main
