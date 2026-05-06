MODULE f_base
  USE f_params
  IMPLICIT NONE
CONTAINS
  SUBROUTINE set_flags()
    USE itg_R, ONLY: set_R_flag
    USE itg_k, ONLY: set_k_flag
    USE itg_f, ONLY: set_f_flag
    CALL set_f_flag()
    CALL set_k_flag()
    CALL set_R_flag()
  END SUBROUTINE set_flags
END MODULE f_base
