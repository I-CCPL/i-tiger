!------------------------------------------------------------------------------
!... IMPORTANT:
!...   This file uses Fortran 2018 assumed-rank dummy arguments,
!...   declared with `..` (ISO/IEC 1539-1:2018, 8.5.8.7, R825).
!...   Build it with a compiler and MPI interface
!...   that support assumed-rank dummy arguments.
!------------------------------------------------------------------------------

MODULE mp_base
  USE mp_global
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: mp_bcast, mp_sum, mp_min, mp_max

  INTERFACE mp_bcast
    MODULE PROCEDURE mp_bcast_logical, mp_bcast_int, &
      mp_bcast_real, mp_bcast_cmplx, &
      mp_bcast_char
  END INTERFACE mp_bcast

  INTERFACE mp_sum
    MODULE PROCEDURE mp_sum_int, mp_sum_real, mp_sum_cmplx
  END INTERFACE mp_sum

  INTERFACE mp_min
    MODULE PROCEDURE mp_min_int
  END INTERFACE mp_min

  INTERFACE mp_max
    MODULE PROCEDURE mp_max_int
  END INTERFACE mp_max
CONTAINS

  SUBROUTINE mp_bcast_logical(msg)
    LOGICAL, INTENT(INOUT), CONTIGUOUS :: msg(..)
    INTEGER :: msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    IF (msg_size <= 0) THEN
      CALL mp_abort(1, 'Invalid message size for mp_bcast.')
    END IF
    CALL MPI_BCAST(msg, msg_size, MPI_LOGICAL, mp_root, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI bcast logical failed.')
#endif
  END SUBROUTINE mp_bcast_logical

  SUBROUTINE mp_bcast_int(msg)
    INTEGER, INTENT(INOUT), CONTIGUOUS :: msg(..)
    INTEGER :: msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    IF (msg_size <= 0) THEN
      CALL mp_abort(1, 'Invalid message size for mp_bcast.')
    END IF
    CALL MPI_BCAST(msg, msg_size, MPI_INTEGER, mp_root, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI bcast integer failed.')
#endif
  END SUBROUTINE mp_bcast_int

  SUBROUTINE mp_bcast_real(msg)
    USE kinds, ONLY: DP
    REAL(DP), INTENT(INOUT), CONTIGUOUS :: msg(..)
    INTEGER :: msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    IF (msg_size <= 0) THEN
      CALL mp_abort(1, 'Invalid message size for mp_bcast.')
    END IF
    CALL MPI_BCAST(msg, msg_size, MPI_DOUBLE_PRECISION, mp_root, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI bcast real failed.')
#endif
  END SUBROUTINE mp_bcast_real

  SUBROUTINE mp_bcast_cmplx(msg)
    USE kinds, ONLY: DP
    COMPLEX(DP), INTENT(INOUT), CONTIGUOUS :: msg(..)
    INTEGER :: msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    IF (msg_size <= 0) THEN
      CALL mp_abort(1, 'Invalid message size for mp_bcast.')
    END IF
    CALL MPI_BCAST(msg, msg_size, MPI_DOUBLE_COMPLEX, mp_root, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI bcast complex failed.')
#endif
  END SUBROUTINE mp_bcast_cmplx

  SUBROUTINE mp_bcast_char(msg)
#ifdef __MPI
    USE mp_global, ONLY: mp_abort, mp_comm, mp_root, ierr, MPI_CHARACTER
#endif
    CHARACTER(LEN=*), INTENT(INOUT) :: msg
#ifdef __MPI
    IF (LEN(msg) <= 0) THEN
      CALL mp_abort(1, 'Invalid message size for mp_bcast.')
    END IF
    CALL MPI_BCAST(msg, LEN(msg), MPI_CHARACTER, mp_root, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI bcast character failed.')
#endif
  END SUBROUTINE mp_bcast_char
  ! ==================================================
  SUBROUTINE mp_sum_int(msg)
    INTEGER, INTENT(INOUT), CONTIGUOUS :: msg(..)
    INTEGER :: msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    IF (msg_size <= 0) THEN
      CALL mp_abort(1, 'Invalid message size for mp_sum.')
    END IF
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, msg, msg_size, MPI_INTEGER, MPI_SUM, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI sum integer failed.')
#endif
  END SUBROUTINE mp_sum_int

  SUBROUTINE mp_sum_real(msg)
    USE kinds, ONLY: DP
    REAL(DP), INTENT(INOUT), CONTIGUOUS :: msg(..)
    INTEGER :: msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    IF (msg_size <= 0) THEN
      CALL mp_abort(1, 'Invalid message size for mp_sum.')
    END IF
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, msg, msg_size, MPI_DOUBLE_PRECISION, MPI_SUM, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI sum real failed.')
#endif
  END SUBROUTINE mp_sum_real

  SUBROUTINE mp_sum_cmplx(msg)
    USE kinds, ONLY: DP
    COMPLEX(DP), INTENT(INOUT), CONTIGUOUS :: msg(..)
    INTEGER :: msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    IF (msg_size <= 0) THEN
      CALL mp_abort(1, 'Invalid message size for mp_sum.')
    END IF
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, msg, msg_size, MPI_DOUBLE_COMPLEX, MPI_SUM, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI sum complex failed.')
#endif
  END SUBROUTINE mp_sum_cmplx
  ! ==================================================
  SUBROUTINE mp_min_int(msg)
    INTEGER, INTENT(INOUT), CONTIGUOUS :: msg(..)
    INTEGER :: msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    IF (msg_size <= 0) THEN
      CALL mp_abort(1, 'Invalid message size for mp_min.')
    END IF

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, msg, msg_size, MPI_INTEGER, MPI_MIN, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI min integer failed.')
#endif
  END SUBROUTINE mp_min_int
  ! ==================================================
  SUBROUTINE mp_max_int(msg)
    INTEGER, INTENT(INOUT), CONTIGUOUS :: msg(..)
    INTEGER :: msg_size
#ifdef __MPI
    msg_size = SIZE(msg)
    IF (msg_size <= 0) THEN
      CALL mp_abort(1, 'Invalid message size for mp_max.')
    END IF

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, msg, msg_size, MPI_INTEGER, MPI_MAX, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI max integer failed.')
#endif
  END SUBROUTINE mp_max_int
END MODULE mp_base
