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
  LOGICAL::lBerry = .FALSE.
  REAL(DP)::dE_thr = 1D-8
  !< threshold for identifying degenerate states in eV
  REAL(DP)::dE_eta = 0.04
  !< broadening parameter for 1/dE
  REAL(DP)::E_fermi = 0.0_DP

  !... Shift current calculation parameters
  !... Frequency: eV units from input
  LOGICAL::lShift = .FALSE.
  REAL(DP)::shift_wmin = 0.0_DP
  REAL(DP)::shift_wmax = 0.0_DP
  REAL(DP)::shift_dw = 0.0_DP
  INTEGER::shift_nw
  REAL(DP)::shift_eta = 0.01_DP
  !
CONTAINS
  SUBROUTINE read_input()
    USE char_mod, ONLY: captital
    USE io_global, ONLY: stdout, prefix
    USE wannier90, ONLY: w90data
    CHARACTER(LEN=256)::line
    CHARACTER(LEN=80)::card
    LOGICAL::tend
    INTEGER::i
    !
    CALL start_clock('read_input')
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
    !
    WRITE (stdout, '(2X, A)') 'Reading input completed.'
    CALL stop_clock('read_input')
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
    USE system, ONLY: dim
    NAMELIST /itg/ lBand, lOAM, lBerry, dE_thr, dE_eta, E_fermi, &
      dim, lShift, shift_wmin, shift_wmax, shift_dw, shift_eta
    WRITE (stdout, '(2X, A)') 'Reading &ITG Namelist...'
    IF (ionode) THEN
      READ (stdin, nml=itg)
      WRITE (stdout, '(2X, A, ES11.4)') '- dE threshold: ', dE_thr
      WRITE (stdout, '(2X, A, ES11.4)') '- Fermi energy: ', E_fermi
      WRITE (stdout, '(2X, A, 1X, I0)') '- Dimension: ', dim
      IF (lShift) THEN
        shift_nw = CEILING((shift_wmax - shift_wmin)/shift_dw) + 1
        WRITE (stdout, '(2X, A, 2(1X, ES11.4))') '- Shift energy window (eV): ', shift_wmin, shift_wmax
        WRITE (stdout, '(2X, A, 1X, ES11.4)') '- Shift broadening (eV): ', shift_eta
        WRITE (stdout, '(2X, A, 1X, ES11.4)') '- Shift dw (eV): ', shift_dw
        WRITE (stdout, '(2X, A, 1X, I0)') '- Shift w points: ', shift_nw
      END IF
    END IF
    CALL mp_bcast(lBand)
    CALL mp_bcast(lOAM)
    CALL mp_bcast(lBerry)
    CALL mp_bcast(dE_thr)
    CALL mp_bcast(E_fermi)
    CALL mp_bcast(dim)
    IF (dim < 1 .OR. dim > 3) THEN
      CALL errore(1, 'read_itg', 'dimensionality must be 1, 2, or 3')
    END IF
    CALL mp_bcast(lShift)
    IF (lShift) THEN
      CALL mp_bcast(shift_wmin)
      CALL mp_bcast(shift_wmax)
      CALL mp_bcast(shift_eta)
      CALL mp_bcast(shift_dw)
      CALL mp_bcast(shift_nw)
      IF (shift_wmax <= shift_wmin) &
        CALL errore(1, 'read_itg', 'shift_wmax must be greater than shift_wmin')
      IF (shift_dw <= 0.0_DP) &
        CALL errore(1, 'read_itg', 'shift_dw must be positive')
      IF (shift_nw <= 0) &
        CALL errore(1, 'read_itg', 'shift_nw must be positive')
      IF (shift_eta <= 0.0_DP) &
        CALL errore(1, 'read_itg', 'shift_eta must be positive')
    END IF
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
    USE system, ONLY: red2cart_recip
    CHARACTER(LEN=256), INTENT(INOUT)::line
    LOGICAL::tend
    INTEGER::i
    !> mesh
    INTEGER::nk1, nk2, nk3, sk1, sk2, sk3
    !> path
    INTEGER::npath
    REAL(DP), ALLOCATABLE::skp(:, :)
    INTEGER, ALLOCATABLE::nkpps(:)
    REAL(DP)::dk_red(3), dk_cart(3), k_pos
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
      WRITE (stdout, '(2X, A, 3(1X, I0))') '- K-Mesh: ', nk1, nk2, nk3
      WRITE (stdout, '(2X, A, 3(1X, I0))') '- Shift:  ', sk1, sk2, sk3
      WRITE (stdout, '(2X,A, 1X, I0)') '- Total k-points: ', t_kpt%nktot
      !
    ELSE IF (match('CRYSTAL_B', line)) THEN
      IF (lShift) CALL errore(1, 'read_kpts', &
                              'shift current requires a uniform k-mesh.' &
                              //' Use AUTOMATIC keyword.')

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

      k_pos = 0.0_DP
      WRITE (stdout, 5413) 1, k_pos
      DO i = 2, npath
        dk_red = skp(:, i) - skp(:, i - 1)
        CALL red2cart_recip(dk_red, dk_cart)
        k_pos = k_pos + NORM2(dk_cart)
        WRITE (stdout, 5413) i, k_pos
      END DO
5413  FORMAT(2X, '- high sym. k pos(', I0, '): ', F10.4)
      !
    END IF

    RETURN
10  CALL errore(1, 'read_kpts', 'end of file')
  END SUBROUTINE read_kpts
END MODULE io_input
