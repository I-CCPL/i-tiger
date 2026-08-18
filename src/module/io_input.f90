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
    USE f_params, ONLY: FFT_conv, convention, formula
    NAMELIST /control/ prefix, debug, FFT_conv, convention, formula
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

    CALL mp_bcast(FFT_conv)
    CALL mp_bcast(convention)
    CALL mp_bcast(formula)
    IF (TRIM(FFT_conv) /= 'simple' &
        .AND. TRIM(FFT_conv) /= 'TB' &
        .AND. TRIM(FFT_conv) /= 'Test') THEN
      CALL errore(1, 'read_itg', 'Unknown FFT convention: "'//TRIM(FFT_conv)//'"')
    END IF
    IF (TRIM(formula) /= 'gauge' &
        .AND. TRIM(formula) /= 'projection') THEN
      CALL errore(1, 'read_itg', 'Unknown formula: "'//TRIM(formula)//'"')
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
    USE f_base, ONLY: set_flags
    USE f_params
    NAMELIST /itg/ lBand, lOAM, lOAM_g, &
      lBerry, lBerry_p, lBerry_g, &
      lBCD, lBCD_p, &
      dE_thr, dE_eta, E_fermi, Ef_min, Ef_max, Ef_step, & ! dim &
      lNLO, lNLO_g, lshift_g, lshift_E_g, shift_hw, ldielec_E_g, &
      lshift_k_g, ldielec_k_g, lshift_vec_g, shift_vec_b, &
      NLO_Emin, NLO_Emax, NLO_dE, NLO_eta, NLO_w_thr
    WRITE (stdout, '(2X, A)') 'Reading &ITG Namelist...'
    IF (ionode) READ (stdin, nml=itg)
    CALL mp_bcast(lBand)
    CALL mp_bcast(lOAM)
    CALL mp_bcast(lOAM_g)
    CALL mp_bcast(lBerry)
    CALL mp_bcast(lBerry_p)
    CALL mp_bcast(lBerry_g)
    CALL mp_bcast(lBCD)
    CALL mp_bcast(lBCD_p)
    CALL mp_bcast(lNLO)
    CALL mp_bcast(lNLO_g)
    CALL mp_bcast(lshift_g)
    CALL mp_bcast(lshift_E_g)
    CALL mp_bcast(ldielec_E_g)
    CALL mp_bcast(lshift_k_g)
    CALL mp_bcast(ldielec_k_g)
    CALL mp_bcast(lshift_vec_g)
    CALL mp_bcast(shift_vec_b)
    CALL set_flags()
    !
    WRITE (stdout, '(2X, A, ES11.4)') '- dE threshold: ', dE_thr
    WRITE (stdout, '(2X, A, ES11.4)') '- Fermi energy: ', E_fermi
    ! WRITE (stdout, '(2X, A, 1X, I0)') '- Dimension: ', dim
    IF (lNLO .OR. lshift_E_g .OR. ldielec_E_g .OR. &
        lshift_k_g .OR. ldielec_k_g .OR. lshift_vec_g) THEN
      NLO_nE = CEILING((NLO_Emax - NLO_Emin)/NLO_dE) + 1
      WRITE (stdout, '(2X, A, 2(1X, ES11.4))') '- NLO energy window (eV): ', NLO_Emin, NLO_Emax
      WRITE (stdout, '(2X, A, 1X, ES11.4)') '- NLO broadening (eV): ', NLO_eta
      WRITE (stdout, '(2X, A, 1X, ES11.4)') '- NLO dE (eV): ', NLO_dE
      WRITE (stdout, '(2X, A, 1X, I0)') '- NLO E points: ', NLO_nE
      WRITE (stdout, '(2X, A, 1X, ES11.4)') '- NLO frequency threshold (eV): ', NLO_w_thr
    END IF
    !
    CALL mp_bcast(dE_thr)
    CALL mp_bcast(dE_eta)
    CALL mp_bcast(E_fermi)
    CALL mp_bcast(Ef_min)
    CALL mp_bcast(Ef_max)
    CALL mp_bcast(Ef_step)
    IF (lBCD_p) THEN
      IF (Ef_min > Ef_max) THEN
        CALL errore(1, 'read_itg', 'Ef_max >= Ef_min required.')
      END IF
      IF (Ef_step <= 0.0_DP) THEN
        CALL errore(1, 'read_itg', 'Ef_step > 0 required.')
      END IF
      Ef_nE = NINT((Ef_max - Ef_min)/Ef_step) + 1
    END IF
    ! CALL mp_bcast(dim)
    ! IF (dim < 1 .OR. dim > 3) THEN
    !   CALL errore(1, 'read_itg', 'dimensionality must be 1, 2, or 3')
    ! END IF
    IF (lNLO .OR. lshift_E_g .OR. ldielec_E_g .OR. &
        lshift_k_g .OR. ldielec_k_g .OR. lshift_vec_g) THEN
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
    IF (lshift_E_g .OR. ldielec_E_g .OR. &
        lshift_k_g .OR. ldielec_k_g) THEN
      CALL mp_bcast(shift_hw)
      IF (shift_hw <= 0.0_DP) THEN
        CALL errore(1, 'read_itg', 'shift_hw must be positive')
      END IF
    END IF
    IF (lshift_vec_g) THEN
      IF (shift_vec_b < 1 .OR. shift_vec_b > 3) &
        CALL errore(1, 'read_itg', 'shift_vec_b must be 1, 2, or 3')
      WRITE (stdout, '(2X, A, 1X, I0)') '- Shift-vector polarization: ', shift_vec_b
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
    USE io_global, ONLY: ionode, get_free_unit
    USE char_mod, ONLY: match
    USE kpoints, ONLY: t_kpt
    USE system, ONLY: red2cart_recip
    USE f_params, ONLY: lNLO, lshift_E_g
    CHARACTER(LEN=256), INTENT(INOUT)::line
    LOGICAL::tend
    INTEGER::i
    INTEGER::iun_kpts
    CHARACTER(LEN=*), PARAMETER::fn_kpts = 'itg.kpts.dat'
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
      IF (ionode) THEN
        iun_kpts = get_free_unit()
        OPEN (unit=iun_kpts, file=fn_kpts, &
              status='replace', action='write')
        WRITE (iun_kpts, '(A)') '# k_idx, k_red'
        DO i = 1, t_kpt%nktot
          WRITE (iun_kpts, '(I8, 3(1X, ES11.4))') i, t_kpt%k_red(:, i)
        END DO
        CLOSE (iun_kpts)
      END IF
      !
    ELSE IF (match('CRYSTAL_B', line)) THEN
      IF (lNLO .OR. lshift_E_g) THEN
        CALL errore(1, 'read_kpts', &
                    'NLO requires a uniform k-mesh.' &
                    //' Use AUTOMATIC keyword.')
      END IF

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
      IF (ionode) THEN
        iun_kpts = get_free_unit()
        OPEN (unit=iun_kpts, file=fn_kpts, &
              status='replace', action='write')
        WRITE (iun_kpts, '(A)') '# k_idx, k_red, kpos'
        k_pos = 0.0_DP
        DO i = 1, t_kpt%nktot
          WRITE (iun_kpts, '(3(1X, ES11.4))') i, t_kpt%k_red(:, i), k_pos
          IF (i < t_kpt%nktot) &
            k_pos = k_pos + NORM2(t_kpt%k_cart(:, i + 1) - t_kpt%k_cart(:, i))
        END DO
        CLOSE (iun_kpts)
      END IF
      !
    END IF

    RETURN
10  CALL errore(1, 'read_kpts', 'end of file')
  END SUBROUTINE read_kpts
END MODULE io_input
