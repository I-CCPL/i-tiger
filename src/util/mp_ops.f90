
! Call: mp_global -> mp_blk -> mp_base
! File seperated for array flatten.
#define __MSG_MAXSIZE 100000

SUBROUTINE mp_ops_bcast_logical(msg, msg_size)
  USE mp_global, ONLY: mp_abort, mp_comm, mp_root, ierr, MPI_LOGICAL
  IMPLICIT NONE
  LOGICAL, INTENT(inout) :: msg(msg_size)
  INTEGER, INTENT(in) :: msg_size
#ifdef __MPI
  INTEGER::iblk, nblk, blk_size, istart
  IF (msg_size <= 0) THEN
    CALL mp_abort(1, 'Invalid message size for mp_bcast.')
  END IF
  !
  IF (msg_size <= __MSG_MAXSIZE) THEN
    CALL MPI_BCAST(msg, msg_size, MPI_LOGICAL, mp_root, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI bcast logical failed.')
  ELSE
    blk_size = __MSG_MAXSIZE
    nblk = msg_size/blk_size
    DO iblk = 1, nblk
      istart = (iblk - 1)*blk_size + 1
      CALL MPI_BCAST(msg(istart), blk_size, MPI_LOGICAL, mp_root, mp_comm, ierr)
      CALL mp_abort(ierr, 'MPI bcast logical failed.')
    END DO
    blk_size = MOD(msg_size, blk_size)
    IF (blk_size > 0) THEN
      istart = nblk*__MSG_MAXSIZE + 1
      CALL MPI_BCAST(msg(istart), blk_size, MPI_LOGICAL, mp_root, mp_comm, ierr)
      CALL mp_abort(ierr, 'MPI bcast logical failed.')
    END IF
  END IF
#endif
END SUBROUTINE mp_ops_bcast_logical
!
SUBROUTINE mp_ops_bcast_int(msg, msg_size)
  USE mp_global, ONLY: mp_abort, mp_comm, mp_root, ierr, MPI_INTEGER
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: msg(msg_size)
  INTEGER, INTENT(IN) :: msg_size
#ifdef __MPI
  INTEGER::iblk, nblk, blk_size, istart
  IF (msg_size <= 0) THEN
    CALL mp_abort(1, 'Invalid message size for mp_bcast.')
  END IF
  !
  IF (msg_size <= __MSG_MAXSIZE) THEN
    CALL MPI_BCAST(msg, msg_size, MPI_INTEGER, mp_root, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI bcast integer failed.')
  ELSE
    blk_size = __MSG_MAXSIZE
    nblk = msg_size/blk_size
    DO iblk = 1, nblk
      istart = (iblk - 1)*blk_size + 1
      CALL MPI_BCAST(msg(istart), blk_size, MPI_INTEGER, mp_root, mp_comm, ierr)
      CALL mp_abort(ierr, 'MPI bcast integer failed.')
    END DO
    blk_size = MOD(msg_size, blk_size)
    IF (blk_size > 0) THEN
      istart = nblk*__MSG_MAXSIZE + 1
      CALL MPI_BCAST(msg(istart), blk_size, MPI_INTEGER, mp_root, mp_comm, ierr)
      CALL mp_abort(ierr, 'MPI bcast integer failed.')
    END IF
  END IF
#endif
END SUBROUTINE mp_ops_bcast_int
!
SUBROUTINE mp_ops_bcast_real(msg, msg_size)
  USE kinds, ONLY: DP
  USE mp_global, ONLY: mp_abort, mp_comm, mp_root, ierr, MPI_DOUBLE_PRECISION
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: msg(msg_size)
  INTEGER, INTENT(IN) :: msg_size
#ifdef __MPI
  INTEGER::iblk, nblk, blk_size, istart
  IF (msg_size <= 0) THEN
    CALL mp_abort(1, 'Invalid message size for mp_bcast.')
  END IF
  !
  IF (msg_size <= __MSG_MAXSIZE) THEN
    CALL MPI_BCAST(msg, msg_size, MPI_DOUBLE_PRECISION, mp_root, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI bcast real failed.')
  ELSE
    blk_size = __MSG_MAXSIZE
    nblk = msg_size/blk_size
    DO iblk = 1, nblk
      istart = (iblk - 1)*blk_size + 1
      CALL MPI_BCAST(msg(istart), blk_size, MPI_DOUBLE_PRECISION, mp_root, mp_comm, ierr)
      CALL mp_abort(ierr, 'MPI bcast real failed.')
    END DO
    blk_size = MOD(msg_size, blk_size)
    IF (blk_size > 0) THEN
      istart = nblk*__MSG_MAXSIZE + 1
      CALL MPI_BCAST(msg(istart), blk_size, MPI_DOUBLE_PRECISION, mp_root, mp_comm, ierr)
      CALL mp_abort(ierr, 'MPI bcast real failed.')
    END IF
  END IF
#endif
END SUBROUTINE mp_ops_bcast_real
!
SUBROUTINE mp_ops_bcast_cmplx(msg, msg_size)
  USE kinds, ONLY: DP
  USE mp_global, ONLY: mp_abort, mp_comm, mp_root, ierr, MPI_DOUBLE_COMPLEX
  IMPLICIT NONE
  COMPLEX(DP), INTENT(INOUT) :: msg(msg_size)
  INTEGER, INTENT(IN) :: msg_size
#ifdef __MPI
  INTEGER::iblk, nblk, blk_size, istart
  IF (msg_size <= 0) THEN
    CALL mp_abort(1, 'Invalid message size for mp_bcast.')
  END IF
  !
  IF (msg_size <= __MSG_MAXSIZE) THEN
    CALL MPI_BCAST(msg, msg_size, MPI_DOUBLE_COMPLEX, mp_root, mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI bcast complex failed.')
  ELSE
    blk_size = __MSG_MAXSIZE
    nblk = msg_size/blk_size
    DO iblk = 1, nblk
      istart = (iblk - 1)*blk_size + 1
      CALL MPI_BCAST(msg(istart), blk_size, MPI_DOUBLE_COMPLEX, mp_root, mp_comm, ierr)
      CALL mp_abort(ierr, 'MPI bcast complex failed.')
    END DO
    blk_size = MOD(msg_size, blk_size)
    IF (blk_size > 0) THEN
      istart = nblk*__MSG_MAXSIZE + 1
      CALL MPI_BCAST(msg(istart), blk_size, MPI_DOUBLE_COMPLEX, mp_root, mp_comm, ierr)
      CALL mp_abort(ierr, 'MPI bcast complex failed.')
    END IF
  END IF
#endif
END SUBROUTINE mp_ops_bcast_cmplx

! ================================================== !

SUBROUTINE mp_ops_sum_int(msg, msg_size)
  USE mp_global, ONLY: mp_abort, mp_comm, ierr, MPI_INTEGER, MPI_SUM, MPI_IN_PLACE
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: msg(msg_size)
  INTEGER, INTENT(IN) :: msg_size
#ifdef __MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, msg, msg_size, MPI_INTEGER, MPI_SUM, mp_comm, ierr)
  CALL mp_abort(ierr, 'MPI sum integer failed.')
#endif
END SUBROUTINE mp_ops_sum_int
!
SUBROUTINE mp_ops_sum_real(msg, msg_size)
  USE kinds, ONLY: DP
  USE mp_global, ONLY: mp_abort, mp_comm, ierr, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_IN_PLACE
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: msg(msg_size)
  INTEGER, INTENT(IN) :: msg_size
#ifdef __MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, msg, msg_size, MPI_DOUBLE_PRECISION, MPI_SUM, mp_comm, ierr)
  CALL mp_abort(ierr, 'MPI sum real failed.')
#endif
END SUBROUTINE mp_ops_sum_real
!
SUBROUTINE mp_ops_sum_cmplx(msg, msg_size)
  USE kinds, ONLY: DP
  USE mp_global, ONLY: mp_abort, mp_comm, ierr, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_IN_PLACE
  IMPLICIT NONE
  COMPLEX(DP), INTENT(INOUT) :: msg(msg_size)
  INTEGER, INTENT(IN) :: msg_size
#ifdef __MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, msg, msg_size, MPI_DOUBLE_COMPLEX, MPI_SUM, mp_comm, ierr)
  CALL mp_abort(ierr, 'MPI sum complex failed.')
#endif
END SUBROUTINE mp_ops_sum_cmplx
