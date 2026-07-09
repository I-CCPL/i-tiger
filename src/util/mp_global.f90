MODULE mp_global
  IMPLICIT NONE
#ifdef __MPI
  INCLUDE 'mpif.h'
#endif
  INTEGER::mp_comm, mp_root, mp_rank, mp_size, ierr
  !
CONTAINS
  SUBROUTINE mp_start()
    LOGICAL::binit
#ifdef __MPI
    mp_comm = MPI_COMM_WORLD
    mp_root = 0
    CALL MPI_INITIALIZED(binit, ierr)
    IF (.NOT. binit) THEN
      CALL MPI_INIT(ierr)
      CALL mp_abort(ierr, 'MPI_Init failed.')
    END IF
    CALL MPI_COMM_RANK(mp_comm, mp_rank, ierr)
    CALL mp_abort(ierr, 'MPI_Comm_rank failed.')
    CALL MPI_COMM_SIZE(mp_comm, mp_size, ierr)
    CALL mp_abort(ierr, 'MPI_Comm_size failed.')
#else
    mp_comm = 0
    mp_root = 0
    mp_rank = 0
    mp_size = 1
#endif
  END SUBROUTINE mp_start
  !
  SUBROUTINE mp_end()
    LOGICAL::bfinal
#ifdef __MPI
    CALL MPI_FINALIZED(bfinal, ierr)
    IF (.NOT. bfinal) THEN
      CALL MPI_FINALIZE(ierr)
      CALL mp_abort(ierr, 'MPI_Finalize failed.')
    END IF
#endif
  END SUBROUTINE mp_end
  !
  SUBROUTINE mp_abort(code, msg)
    USE io_param, ONLY: stdout
#ifdef __INTEL_COMPILER
    USE ifcore, ONLY: tracebackqq
#endif
    INTEGER, INTENT(IN) :: code
    CHARACTER(len=*), INTENT(IN), OPTIONAL::msg
    IF (code == 0) RETURN
    WRITE (stdout, *)
    IF (PRESENT(msg)) THEN
      WRITE (stdout, '("Error: ",A)') TRIM(msg)
      WRITE (stdout, *)
    END IF
#ifdef __INTEL_COMPILER
    CALL tracebackqq(user_exit_code=1)
#elif defined __GFORTRAN
    CALL BACKTRACE()
#endif
#ifdef __MPI
    CALL MPI_ABORT(mp_comm, code, ierr)
#endif
    STOP code
  END SUBROUTINE mp_abort
  !
  SUBROUTINE mp_barrier()
#ifdef __MPI
    CALL MPI_BARRIER(mp_comm, ierr)
    CALL mp_abort(ierr, 'MPI_BARRIER failed.')
#endif
  END SUBROUTINE mp_barrier
  ! ================================================== !

END MODULE mp_global
