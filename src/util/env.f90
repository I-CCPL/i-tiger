MODULE env
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  CHARACTER(LEN=12)::itg_version = 'v0.0.11.3.2'
CONTAINS
  SUBROUTINE env_start(date, time)
    USE mp_global, ONLY: mp_start, mp_rank, mp_root, mp_size, mp_barrier
    USE io_global, ONLY: ionode, stdout
    CHARACTER(LEN=*), INTENT(IN)::date, time
    CHARACTER(LEN=10)::cdate, ctime
    !
    CALL mp_start()
    CALL start_clock('i-TIGER')
    ionode = (mp_rank == mp_root)
    IF (.NOT. ionode) THEN
      OPEN (unit=stdout, file='/dev/null', status='unknown')
    END IF
    !
    CALL current_date_time(cdate, ctime)
    WRITE (stdout, '(A)') 'Starting i-TIGER on '//TRIM(cdate)//' '//TRIM(ctime)
    WRITE (stdout, '(2X, 4A)') "= Compiled on ", date, ' ', time
    WRITE (stdout, '(2X,A,I0,A)') '= Running on ', mp_size, ' processors.'
    WRITE (stdout, *)
    CALL write_bold_line()
    WRITE (stdout, '(2X,A)') 'Incheon Tight-binding Induced Generalized Electronic Response'
    WRITE (stdout, '(2X,A,A)') 'i-TIGER ', TRIM(itg_version)
    CALL write_bold_line()
    WRITE (stdout, *)
    CALL mp_barrier()
  END SUBROUTINE env_start
  !
  SUBROUTINE print_k_info(ikpt, nkpt)
    USE kinds, ONLY: DP
    INTEGER, INTENT(IN)::ikpt, nkpt
    REAL(DP), EXTERNAL::get_clock
    IF (MOD(ikpt, 100) == 1 .OR. ikpt - 1 == nkpt) THEN
      WRITE (stdout, 2398) ikpt - 1, nkpt, get_clock('i-TIGER')
    END IF
2398 FORMAT(2X, "- k-point ", I0, " / ", I0, " (Elapsed: ", F0.1, "s)")
  END SUBROUTINE print_k_info
  !
  SUBROUTINE env_end()
    USE mp_global, ONLY: mp_end
    CHARACTER(LEN=10)::cdate, ctime
    CALL stop_clock('i-TIGER')
    CALL memory_report()
    CALL print_all_clocks()
    !
    CALL mp_end()
    CALL current_date_time(cdate, ctime)
    WRITE (stdout, *)
    WRITE (stdout, '(A)') 'Finished i-TIGER on '//TRIM(cdate)//' '//TRIM(ctime)
    FLUSH (stdout)
  END SUBROUTINE env_end
  !
  SUBROUTINE current_date_time(cdate, ctime)
    CHARACTER(LEN=10), INTENT(OUT)::cdate
    CHARACTER(LEN=10), INTENT(OUT)::ctime
    CHARACTER(LEN=1), PARAMETER::delim1 = '-', delim2 = ':'
    INTEGER::date_time(8)
    CALL DATE_AND_TIME(values=date_time)
    WRITE (cdate, '(I4.4,A1,I2.2,A1,I2.2)') date_time(1), delim1, date_time(2), delim1, date_time(3)
    WRITE (ctime, '(I2.2,A1,I2.2,A1,I2.2)') date_time(5), delim2, date_time(6), delim2, date_time(7)
  END SUBROUTINE current_date_time
  !
  SUBROUTINE memory_report()
    USE ISO_FORTRAN_ENV, ONLY: INT64
    USE io_global, ONLY: ionode, get_free_unit
    USE char_mod, ONLY: add_comma
    USE mp_global
    INTEGER::io_unit, ios
    LOGICAL::have_hwm
    INTEGER(INT64)::hwm_kb, hwm_max_kb
    CHARACTER(LEN=30)::msg
    CHARACTER(LEN=256)::line
    !... Get WmHWM from /proc/self/status
    io_unit = get_free_unit()
    OPEN (unit=io_unit, file='/proc/self/status', &
          status='old', action='read', iostat=ios)
    IF (ios /= 0) RETURN

    have_hwm = .FALSE.
    DO
      READ (io_unit, '(A)', iostat=ios) line
      IF (ios /= 0) EXIT
      IF (INDEX(line, 'VmHWM:') == 1) THEN
        READ (line(7:), *) hwm_kb
        have_hwm = .TRUE.
      END IF
    END DO
    CLOSE (io_unit)

    IF (have_hwm) THEN
#ifdef __MPI
      CALL MPI_REDUCE(hwm_kb, hwm_max_kb, 1, MPI_INTEGER8, MPI_MAX, mp_root, mp_comm, ierr)
#else
      hwm_max_kb = hwm_kb
#endif
      WRITE (msg, '(I0)') hwm_max_kb
      IF (ionode) CALL add_comma(msg)
      WRITE (stdout, '(2X,A)') 'Peak memory usage (VmHWM): '//TRIM(msg)//' KB'
      WRITE (stdout, *)
    END IF
  END SUBROUTINE memory_report
  !
  SUBROUTINE print_all_clocks()
    WRITE (stdout, '(2X,A)') 'TIMER REPORT: CPU_TIME and WALL_TIME'
    WRITE (stdout, '(2X,A)') 'utility functions'
    CALL print_clock('read_input')
    CALL print_clock('rotate_3d')
    CALL print_clock('fft_q2R')
    CALL print_clock('fft_R2k')

    WRITE (stdout, '(2X,A)') 'main functions'
    CALL print_clock('make_q')
    CALL print_clock('make_R')
    CALL print_clock('make_k')
    CALL print_clock('make_f')
    CALL print_clock('NLO_main')
    CALL print_clock('write_k_data')
    CALL print_clock('i-TIGER')
  END SUBROUTINE print_all_clocks
END MODULE
