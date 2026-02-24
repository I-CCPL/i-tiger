MODULE mp_base
  IMPLICIT NONE
#ifdef __MPI
  INCLUDE 'mpif.h'
#endif
  PRIVATE
  PUBLIC::mp_bcast, mp_sum
  !
  INTERFACE mp_bcast
    MODULE PROCEDURE mp_bcast_logical, mp_bcast_int, mp_bcast_real, mp_bcast_cmplx
  END INTERFACE mp_bcast
  !
  INTERFACE mp_sum
    MODULE PROCEDURE mp_sum_int, mp_sum_real, mp_sum_cmplx
  END INTERFACE mp_sum
  !
CONTAINS
  !
  SUBROUTINE mp_bcast_logical(msg)
    LOGICAL, INTENT(inout) :: msg(..)
    INTEGER::msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    CALL mp_ops_bcast_logical(msg, msg_size)
#endif
  END SUBROUTINE mp_bcast_logical
  !
  SUBROUTINE mp_bcast_int(msg)
    INTEGER, INTENT(INOUT) :: msg(..)
    INTEGER::msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    CALL mp_ops_bcast_int(msg, msg_size)
#endif
  END SUBROUTINE mp_bcast_int
  !
  SUBROUTINE mp_bcast_real(msg)
    USE kinds, ONLY: DP
    REAL(DP), INTENT(INOUT) :: msg(..)
    INTEGER::msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    CALL mp_ops_bcast_real(msg, msg_size)
#endif
  END SUBROUTINE mp_bcast_real
  !
  SUBROUTINE mp_bcast_cmplx(msg)
    USE kinds, ONLY: DP
    COMPLEX(DP), INTENT(INOUT) :: msg(..)
    INTEGER::msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    CALL mp_ops_bcast_cmplx(msg, msg_size)
#endif
  END SUBROUTINE mp_bcast_cmplx

  ! ================================================== !

  SUBROUTINE mp_sum_int(msg)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: msg(..)
    INTEGER::msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    CALL mp_ops_sum_int(msg, msg_size)
#endif
  END SUBROUTINE mp_sum_int
  !
  SUBROUTINE mp_sum_real(msg)
    USE kinds, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT) :: msg(..)
    INTEGER::msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    CALL mp_ops_sum_real(msg, msg_size)
#endif
  END SUBROUTINE mp_sum_real
  !
  SUBROUTINE mp_sum_cmplx(msg)
    USE kinds, ONLY: DP
    IMPLICIT NONE
    COMPLEX(DP), INTENT(INOUT) :: msg(..)
    INTEGER::msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    CALL mp_ops_sum_cmplx(msg, msg_size)
#endif
  END SUBROUTINE mp_sum_cmplx
END MODULE mp_base
