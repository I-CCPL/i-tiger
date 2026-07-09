MODULE clocks
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER, PARAMETER::max_clocks = 32
  INTEGER::nclock = 0
  !
  CHARACTER(LEN=16)::clock_name(max_clocks)
  INTEGER::ncalled(max_clocks) = 0
  LOGICAL::is_on(max_clocks) = .FALSE.
  REAL(DP)::s_cpu(max_clocks), s_wall(max_clocks)
  !< start time for CPU and wall clock
  REAL(DP)::t_cpu(max_clocks), t_wall(max_clocks)
  !< total elapsed time for CPU and wall clock
END MODULE clocks

SUBROUTINE get_time(cpu_t, wall_t)
  USE kinds, ONLY: DP
#ifdef __MPI
  USE mp_global, ONLY: MPI_WTIME
#endif
  REAL(DP), INTENT(OUT)::cpu_t, wall_t
  INTEGER::cnt, rate
  !
  CALL CPU_TIME(cpu_t)
#ifdef __MPI
  wall_t = MPI_WTIME()
#else
  CALL SYSTEM_CLOCK(count=cnt, count_rate=rate)
  wall_t = REAL(cnt, DP)/REAL(rate, DP)
#endif
END SUBROUTINE get_time

SUBROUTINE start_clock(label)
  USE clocks
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)::label
  CHARACTER(LEN=16)::name
  INTEGER::i
  !
  name = TRIM(label)
  DO i = 1, nclock
    IF (clock_name(i) == name) THEN
      IF (is_on(i)) THEN
        WRITE (stdout, '("start_clock: ",'// &
               'A, "(#", I0, ") already running.")') TRIM(name), i
        RETURN
      ELSE
        EXIT
      END IF
    END IF
  END DO
  ! initialize new clock
  IF (i > nclock) THEN
    ! maximum reached
    IF (nclock == max_clocks) THEN
      WRITE (stdout, '(A)') "Maximum number of clocks reached."
      RETURN
    END IF
    i = nclock + 1
    nclock = i
    clock_name(i) = name
    ncalled(i) = 0
    t_cpu(i) = 0.0_DP
    t_wall(i) = 0.0_DP
  END IF
  !
  ncalled(i) = ncalled(i) + 1
  CALL get_time(s_cpu(i), s_wall(i))
  is_on(i) = .TRUE.
END SUBROUTINE start_clock

SUBROUTINE stop_clock(label)
  USE clocks
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)::label
  CHARACTER(LEN=16)::name
  INTEGER::i
  REAL(DP)::e_cpu, e_wall
  !
  name = TRIM(label)
  DO i = 1, nclock
    IF (clock_name(i) == name) THEN
      IF (is_on(i)) THEN
        CALL get_time(e_cpu, e_wall)
        t_cpu(i) = t_cpu(i) + e_cpu - s_cpu(i)
        t_wall(i) = t_wall(i) + e_wall - s_wall(i)
        is_on(i) = .FALSE.
      ELSE
        WRITE (stdout, '("stop_clock: ",'// &
               'A, "(#", I0, ") not running.")') TRIM(name), i
      END IF
      RETURN
    END IF
  END DO
  !
  WRITE (stdout, '("stop_clock: ",'// &
         'A, " not found.")') TRIM(name)
END SUBROUTINE stop_clock

FUNCTION get_clock(label) RESULT(wall_time)
  USE clocks
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)::label
  REAL(DP)::wall_time
  REAL(DP)::e_cpu, e_wall
  CHARACTER(LEN=16)::name
  INTEGER::i
  !
  name = TRIM(label)
  DO i = 1, nclock
    IF (clock_name(i) == name) THEN
      IF (is_on(i)) THEN
        CALL get_time(e_cpu, e_wall)
        wall_time = t_wall(i) + e_wall - s_wall(i)
      ELSE
        wall_time = t_wall(i)
      END IF
      RETURN
    END IF
  END DO
  !
  wall_time = 0.0_DP
END FUNCTION get_clock

SUBROUTINE print_clock(label)
  USE clocks
  USE mp_global, ONLY: mp_size
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)::label
  CHARACTER(LEN=16)::name
  INTEGER::i
  REAL(DP)::e_cpu, e_wall
  !
  name = TRIM(label)
  DO i = 1, nclock
    IF (clock_name(i) == name) THEN
      IF (is_on(i)) THEN
        CALL get_time(e_cpu, e_wall)
        t_cpu(i) = t_cpu(i) + e_cpu - s_cpu(i)
        t_wall(i) = t_wall(i) + e_wall - s_wall(i)
      END IF
      WRITE (stdout, 2352) name, t_cpu(i)*mp_size, t_wall(i), ncalled(i)
      RETURN
    END IF
  END DO
2352 FORMAT(2X, "- ", A8, ": ", F10.2, "s ", F10.2, "s ", "(", I10, " calls)")
END SUBROUTINE print_clock
