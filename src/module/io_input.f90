MODULE io_input
  USE kinds, ONLY: DP
  USE io_global, ONLY: stdin, stdout, ionode
  USE mp_base, ONLY: mp_bcast
  IMPLICIT NONE
  PRIVATE::read_control, read_itg, read_line, read_kpts
  !
  LOGICAL::lBand = .FALSE.
  ! LOGICAL::lDOS = .FALSE.
  ! LOGICAL::lPDOS = .FALSE.
  ! INTEGER::dos_dE
  ! INTEGER::dos_Emin
  ! INTEGER::dos_Emax
  LOGICAL::lOAM = .FALSE.
  REAL(DP)::OAM_thr = 1D-8
  !
CONTAINS
  SUBROUTINE read_input()
    USE char_mod, ONLY: captital
    USE io_global, ONLY: stdout, prefix, write_sep_line
    USE wannier90, ONLY: w90data
    CHARACTER(LEN=256)::line
    CHARACTER(LEN=80)::card
    LOGICAL::tend
    INTEGER::i
    !
    !... Read Namelists
    CALL read_control()
    CALL read_itg()
    w90data%prefix = TRIM(prefix)
    CALL w90data%read_files()

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
    USE io_global, ONLY: prefix, debug
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
  SUBROUTINE read_itg()
    NAMELIST /itg/ lBand, lOAM, OAM_thr
    WRITE (stdout, '(2X, A)') 'Reading &ITG Namelist...'
    IF (ionode) THEN
      READ (stdin, nml=itg)
      IF (lOAM) THEN
        WRITE (stdout, '(2X, A, ES11.4)') '- OAM with threshold: ', OAM_thr
      END IF
    END IF
    CALL mp_bcast(lBand)
    CALL mp_bcast(lOAM)
    CALL mp_bcast(OAM_thr)
  END SUBROUTINE read_itg
  !
  SUBROUTINE read_line(line, tend)
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
    USE char_mod, ONLY: match
    USE kpoints, ONLY: t_kpt
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
      CALL mp_bcast(nk1)
      CALL mp_bcast(nk2)
      CALL mp_bcast(nk3)
      CALL mp_bcast(sk1)
      CALL mp_bcast(sk2)
      CALL mp_bcast(sk3)
      CALL t_kpt%build_mesh(nk1, nk2, nk3, sk1, sk2, sk3)
      !
      WRITE (stdout, '(2X, A, (1X, 3I0))') '- K-Mesh: ', nk1, nk2, nk3
      WRITE (stdout, '(2X, A, (1X, 3I0))') '- Shift:  ', sk1, sk2, sk3
      WRITE (stdout, '(2X,A, 1X, I0)') '- Total k-points: ', t_kpt%nktot
      !
    ELSE IF (match('CRYSTAL_B', line)) THEN
      CALL read_line(line, tend)
      IF (tend) GOTO 10
      READ (line, *, END=10) npath

      CALL mp_bcast(npath)
      ALLOCATE (skp(3, npath), nkpps(npath))
      DO i = 1, npath
        CALL read_line(line, tend)
        IF (tend) GOTO 10
        READ (line, *, END=10) skp(:, i), nkpps(i)
      END DO

      CALL mp_bcast(skp)
      CALL mp_bcast(nkpps)
      CALL t_kpt%build_path(npath, skp, nkpps)
      !
      WRITE (stdout, '(2X, A)') '- k-points along the path'
      WRITE (stdout, '(2X, A, I0)') '- Total k-points: ', t_kpt%nktot
      !
    END IF

    ! CALL t_kpt%divide_k()

    RETURN
10  CALL errore(1, 'read_kpts', 'end of file')
  END SUBROUTINE read_kpts
END MODULE io_input
