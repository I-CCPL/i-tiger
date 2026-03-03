MODULE io_input
  IMPLICIT NONE
  PRIVATE
  PUBLIC::read_input
  !
CONTAINS
  SUBROUTINE read_input()
    USE char_mod, ONLY: captital
    USE io_global, ONLY: stdout, prefix, write_sep_line
    USE wannier90, ONLY: w90data, read_w90
    CHARACTER(LEN=256)::line
    CHARACTER(LEN=80)::card
    LOGICAL::tend
    INTEGER::i
    !
    !... Read Namelists
    CALL read_control()
    w90data%prefix = TRIM(prefix)
    CALL read_w90()

    !... Read Cards
    DO
      CALL read_line(line, tend)
      IF (tend) EXIT

      DO i = 1, LEN_TRIM(line)
        line(i:i) = captital(line(i:i))
      END DO
      READ (line, *) card

      SELECT CASE (TRIM(card))
      CASE ("K_POINTS")
        CALL read_kpts(line)
      CASE DEFAULT
        WRITE (stdout, '(A)') 'Warning: card '//TRIM(card)//' ignored.'
      END SELECT
    END DO
    WRITE (stdout, '(2X, A)') 'Reading input completed.'
    CALL write_sep_line()
  END SUBROUTINE read_input
  !
  SUBROUTINE read_control()
    USE io_global, ONLY: ionode, stdin, stdout, &
                         prefix, debug
    NAMELIST /control/ prefix, debug
    !
    WRITE (stdout, '(2X, A)') 'Reading &CONTROL Namelist...'
    IF (ionode) THEN
      READ (stdin, nml=control)
      IF (debug) THEN
        WRITE (stdout, '(2X, A)') 'Debug mode is ON.'
        CALL read_debug()
      END IF
    END IF
  CONTAINS
    SUBROUTINE read_debug()
      USE wannier90, ONLY: chk_w90, Hq_band
      NAMELIST /debug/ chk_w90, Hq_band
      WRITE (stdout, '(2X, A)') 'Reading &DEBUG Namelist...'
      READ (stdin, nml=debug)
    END SUBROUTINE read_debug
  END SUBROUTINE read_control
  !
  SUBROUTINE read_line(line, tend)
    USE io_global, ONLY: ionode, stdin, stdout
    USE mp_base, ONLY: mp_bcast
    CHARACTER(LEN=*), INTENT(OUT)::line
    LOGICAL, INTENT(OUT)::tend
    tend = .FALSE.
    IF (ionode) THEN
6400  READ (stdin, '(A256)', END=6401) line
      IF (line == ' ' .OR. line(1:1) == '#' .OR. line(1:1) == '!') GOTO 6400
      GOTO 6410
6401  tend = .TRUE.
      WRITE (stdout, *)
      WRITE (stdout, '(2X, A)') 'End of input file'
6410  CONTINUE
    END IF
    CALL mp_bcast(tend)
    CALL mp_bcast(line)
  END SUBROUTINE read_line
  !
  SUBROUTINE read_kpts(line)
    USE kinds, ONLY: DP
    USE char_mod, ONLY: match
    USE io_global, ONLY: stdout
    USE kpoints, ONLY: kpts
    CHARACTER(LEN=256), INTENT(INOUT)::line
    LOGICAL::tend
    INTEGER::i
    !> mesh
    INTEGER::nk1, nk2, nk3, sk1, sk2, sk3
    !> path
    INTEGER::npath
    REAL(DP), ALLOCATABLE::skp(:, :)
    INTEGER, ALLOCATABLE::nkpps(:)
    !
    WRITE (stdout, '(2X, A)') 'Reading K_POINTS Cards...'
    IF (match('AUTOMATIC', line)) THEN
      CALL read_line(line, tend)
      IF (tend) GOTO 10
      READ (line, *) nk1, nk2, nk3, sk1, sk2, sk3
      CALL kpts%build_mesh(nk1, nk2, nk3, sk1, sk2, sk3)
      !
    ELSE IF (match('CRYSTAL_B', line)) THEN
      CALL read_line(line, tend)
      IF (tend) GOTO 10
      READ (line, *, END=10) npath
      ALLOCATE (skp(3, npath), nkpps(npath))

      DO i = 1, npath
        CALL read_line(line, tend)
        IF (tend) GOTO 10
        READ (line, *, END=10) skp(:, i), nkpps(i)
      END DO
      CALL kpts%build_path(npath, skp, nkpps)
    END IF

    RETURN
10  CALL errore(1, 'read_kpts', 'end of file')
  END SUBROUTINE read_kpts
END MODULE io_input
