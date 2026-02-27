PROGRAM main
  USE kinds, ONLY: DP
  USE env, ONLY: env_start, env_end
  USE io_input, ONLY: read_input
  USE system, ONLY: Nw
  USE w90, ONLY: w90data
  IMPLICIT NONE
  COMPLEX(DP), ALLOCATABLE::HR(:, :, :)
  !
  CALL env_start()
  CALL read_input()
  !
  CALL w90data%build_Hq()
  !
  CALL env_end()
END PROGRAM main
