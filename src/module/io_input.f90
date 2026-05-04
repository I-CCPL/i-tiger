MODULE io_input
  USE kinds, ONLY: DP, inf_DP
  USE io_global, ONLY: stdin, stdout, ionode
  USE mp_base, ONLY: mp_bcast
  IMPLICIT NONE
  PRIVATE::read_control, read_itg, read_line, read_kpts
  !
  !... debug
  LOGICAL::debug_q = .FALSE.
  LOGICAL::debug_R = .FALSE.
  LOGICAL::debug_k = .FALSE.
  !... itg
  CHARACTER(LEN=8)::FFT_conv = 'atomic'
  !< 'periodic' or 'atomic' or 'wannier' FFT convention for R2k
  INTEGER::convention = 0
  !< 0 for standard convention, consider degenerate R0s (build_ws)
  !< 1 for TB convention, consider degenerate R0s (build_ws)
  !< 2 for Wannier convention, consider only the nearest R (build_R)
  CHARACTER(LEN=20)::formula = 'kubo'
  !< 'kubo' or 'projection'

  LOGICAL::lBand = .FALSE.
  ! LOGICAL::lDOS = .FALSE.
  ! LOGICAL::lPDOS = .FALSE.
  ! INTEGER::dos_dE
  ! INTEGER::dos_Emin
  ! INTEGER::dos_Emax
  LOGICAL::lOAM = .FALSE.
  LOGICAL::lBerry = .FALSE.
  LOGICAL::lBCD = .FALSE.
  !< Berry Curvature Dipole
  REAL(DP)::dE_thr = 1D-8
  !< threshold for identifying degenerate states in eV
  REAL(DP)::dE_eta = 0.04
  !< broadening parameter for 1/dE
  REAL(DP)::E_fermi = 0.0_DP
  !< Now use for occupation number.

  REAL(DP)::Ef_min = inf_DP
  REAL(DP)::Ef_max = -inf_DP
  REAL(DP)::Ef_step = -1.0_DP
  INTEGER::Ef_nE
  !< Fermi energy range for BCD calculation.

  !... Nonlinear optics calculation parameters
  !... Frequency: eV units from input
  LOGICAL::lNLO = .FALSE.
  REAL(DP)::NLO_Emin = 0.0_DP
  REAL(DP)::NLO_Emax = 0.0_DP
  REAL(DP)::NLO_dE = 0.0_DP
  INTEGER::NLO_nE
  REAL(DP)::NLO_eta = 0.01_DP
  REAL(DP)::NLO_w_thr = 5.0_DP
  !< speeding up frequency integration
  !
  !... Wannier90 data
  LOGICAL::lreq_mmn
  !< whether mmn data is required based on requested calculations
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
    lreq_mmn = (lOAM .OR. lBerry .OR. lBCD .OR. lNLO)
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
    END IF
    CALL mp_bcast(debug)
    IF (debug) THEN
      WRITE (stdout, '(2X, A)') 'Debug mode is ON.'
      CALL read_debug()
    END IF
  CONTAINS
    SUBROUTINE read_debug()
      USE wannier90, ONLY: chk_w90, Hq_band
      NAMELIST /debug/ chk_w90, Hq_band, debug_q, debug_R, debug_k
      WRITE (stdout, '(2X, A)') 'Reading &DEBUG Namelist...'
      IF (ionode) THEN
        READ (stdin, nml=debug)
      END IF
      CALL mp_bcast(debug_k)
    END SUBROUTINE read_debug
  END SUBROUTINE read_control
  !
  SUBROUTINE read_itg()
    ! USE system, ONLY: dim
    NAMELIST /itg/ FFT_conv, convention, lBand, lOAM, lBerry, lBCD, &
      formula, dE_thr, dE_eta, E_fermi, Ef_min, Ef_max, Ef_step, & ! dim &
      lNLO, NLO_Emin, NLO_Emax, NLO_dE, NLO_eta, NLO_w_thr
    WRITE (stdout, '(2X, A)') 'Reading &ITG Namelist...'
    IF (ionode) THEN
      READ (stdin, nml=itg)
      WRITE (stdout, '(2X, A, ES11.4)') '- dE threshold: ', dE_thr
      WRITE (stdout, '(2X, A, ES11.4)') '- Fermi energy: ', E_fermi
      ! WRITE (stdout, '(2X, A, 1X, I0)') '- Dimension: ', dim
      IF (lNLO) THEN
        NLO_nE = CEILING((NLO_Emax - NLO_Emin)/NLO_dE) + 1
        WRITE (stdout, '(2X, A, 2(1X, ES11.4))') '- NLO energy window (eV): ', NLO_Emin, NLO_Emax
        WRITE (stdout, '(2X, A, 1X, ES11.4)') '- NLO broadening (eV): ', NLO_eta
        WRITE (stdout, '(2X, A, 1X, ES11.4)') '- NLO dE (eV): ', NLO_dE
        WRITE (stdout, '(2X, A, 1X, I0)') '- NLO E points: ', NLO_nE
        WRITE (stdout, '(2X, A, 1X, ES11.4)') '- NLO frequency threshold (eV): ', NLO_w_thr
      END IF
    END IF
    CALL mp_bcast(FFT_conv)
    IF (TRIM(FFT_conv) /= 'periodic' &
        .AND. TRIM(FFT_conv) /= 'atomic' &
        .AND. TRIM(FFT_conv) /= 'wannier') THEN
      CALL errore(1, 'read_itg', 'Unknown FFT convention: "'//TRIM(FFT_conv)//'"')
    END IF
    CALL mp_bcast(convention)
    !
    CALL mp_bcast(lBand)
    CALL mp_bcast(lOAM)
    CALL mp_bcast(lBerry)
    CALL mp_bcast(lBCD)
    CALL mp_bcast(formula)
    IF (TRIM(formula) /= 'kubo' &
        .AND. TRIM(formula) /= 'projection') THEN
      CALL errore(1, 'read_itg', 'Unknown formula: "'//TRIM(formula)//'"')
    END IF
    CALL mp_bcast(dE_thr)
    CALL mp_bcast(dE_eta)
    CALL mp_bcast(E_fermi)
    CALL mp_bcast(Ef_min)
    CALL mp_bcast(Ef_max)
    CALL mp_bcast(Ef_step)
    IF (lBCD) THEN
      IF (Ef_min > Ef_max) THEN
        CALL errore(1, 'read_itg', 'Ef_max >= Ef_min required.')
      END IF
      IF (Ef_step <= 0.0_DP) THEN
        CALL errore(1, 'read_itg', 'Ef_step > 0 required.')
      END IF
      Ef_nE = CEILING((Ef_max - Ef_min)/Ef_step) + 1
    END IF
    ! CALL mp_bcast(dim)
    ! IF (dim < 1 .OR. dim > 3) THEN
    !   CALL errore(1, 'read_itg', 'dimensionality must be 1, 2, or 3')
    ! END IF
    CALL mp_bcast(lNLO)
    IF (lNLO) THEN
      CALL mp_bcast(NLO_Emin)
      CALL mp_bcast(NLO_Emax)
      CALL mp_bcast(NLO_dE)
      CALL mp_bcast(NLO_eta)
      CALL mp_bcast(NLO_nE)
      CALL mp_bcast(NLO_w_thr)
      IF (NLO_Emax <= NLO_Emin) &
        CALL errore(1, 'read_itg', 'NLO_Emax must be greater than NLO_Emin')
      IF (NLO_dE <= 0.0_DP) &
        CALL errore(1, 'read_itg', 'NLO_dE must be positive')
      IF (NLO_nE <= 0) &
        CALL errore(1, 'read_itg', 'NLO_nE must be positive')
      IF (NLO_eta <= 0.0_DP) &
        CALL errore(1, 'read_itg', 'NLO_eta must be positive')
      IF (NLO_w_thr <= 0.0_DP) &
        CALL errore(1, 'read_itg', 'NLO_w_thr must be positive')
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
      IF (lNLO) CALL errore(1, 'read_kpts', &
                            'NLO requires a uniform k-mesh.' &
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
