MODULE io_output
  USE kinds, ONLY: DP
  USE io_global, ONLY: ionode, stdout, get_free_unit
  USE system, ONLY: Nw
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  LOGICAL, PRIVATE::binit = .FALSE.
  REAL(DP), ALLOCATABLE::k_pos(:)

  PRIVATE::writing_info
CONTAINS
  SUBROUTINE io_output_init()
    USE kpoints, ONLY: t_kpt
    INTEGER::ikpt
    !
    IF (binit) RETURN
    IF (.NOT. ionode) RETURN
    ALLOCATE (k_pos(t_kpt%nktot))
    k_pos(1) = 0.0_DP
    DO ikpt = 2, t_kpt%nktot
      k_pos(ikpt) = k_pos(ikpt - 1) + NORM2(t_kpt%k_cart(:, ikpt) - t_kpt%k_cart(:, ikpt - 1))
    END DO

    binit = .TRUE.
  END SUBROUTINE io_output_init
  !
  SUBROUTINE writing_info(data_name, fname)
    CHARACTER(LEN=*), INTENT(IN) :: data_name, fname
    WRITE (stdout, '(2X, A, A)') '- Writing '//TRIM(data_name)//' data to "'//TRIM(fname)//'"...'
  END SUBROUTINE writing_info
  !
  SUBROUTINE write_band(fname, eigval)
    USE io_global, ONLY: stdout
    CHARACTER(LEN=*), INTENT(IN) :: fname
    REAL(DP), INTENT(IN)::eigval(Nw, t_kpt%nktot)
    INTEGER::io_unit, ikpt, iw
    IF (.NOT. ionode) RETURN
    !
    io_unit = get_free_unit()
    OPEN (unit=io_unit, file=fname)
    CALL writing_info('band structure', fname)
    WRITE (io_unit, '("#", A)') 'k_pos, eigenvalue (eV)'

    DO iw = 1, Nw
      DO ikpt = 1, t_kpt%nktot
        WRITE (io_unit, '(2F10.4)') k_pos(ikpt), eigval(iw, ikpt)
      END DO
      WRITE (io_unit, *) ! blank line
    END DO
    CLOSE (io_unit)
  END SUBROUTINE write_band
  !
  SUBROUTINE write_OAM(fname, L_k)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    REAL(DP), INTENT(IN) :: L_k(3, Nw, t_kpt%nktot)
    INTEGER :: io_unit, ikpt, iw
    IF (.NOT. ionode) RETURN
    !
    io_unit = get_free_unit()
    OPEN (unit=io_unit, file=fname)
    CALL writing_info('OAM', fname)
    WRITE (io_unit, '("#", A)') 'k_pos, OAM (hbar)'

    DO iw = 1, Nw
      DO ikpt = 1, t_kpt%nktot
        WRITE (io_unit, '(F10.4, 3(1X, ES11.4))') k_pos(ikpt), L_k(:, iw, ikpt)
      END DO
      WRITE (io_unit, *) ! blank line
    END DO
    CLOSE (io_unit)
  END SUBROUTINE write_OAM
  !
  SUBROUTINE write_Berry(fname, O_k)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    REAL(DP), INTENT(IN) :: O_k(3, t_kpt%nktot)
    INTEGER :: io_unit, ikpt
    IF (.NOT. ionode) RETURN
    !
    io_unit = get_free_unit()
    OPEN (unit=io_unit, file=fname)
    CALL writing_info('Berry curvature', fname)
    WRITE (io_unit, '("#", A)') 'k_pos, Berry  (Ang^2)'

    DO ikpt = 1, t_kpt%nktot
      WRITE (io_unit, '(F10.4, 3(1X, ES12.4E3))') k_pos(ikpt), O_k(:, ikpt)
    END DO
    WRITE (io_unit, *) ! blank line
    CLOSE (io_unit)
  END SUBROUTINE write_Berry
  !
  SUBROUTINE write_Berry_k(fname, berry_k)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    REAL(DP), INTENT(IN) :: berry_k(3, Nw, t_kpt%nktot)
    INTEGER :: io_unit, ikpt, iw
    IF (.NOT. ionode) RETURN
    !
    io_unit = get_free_unit()
    OPEN (unit=io_unit, file=fname)
    CALL writing_info('Berry curvature', fname)
    WRITE (io_unit, '("#", A)') 'k_pos, Berry  (arb.)'

    DO iw = 1, Nw
      DO ikpt = 1, t_kpt%nktot
        WRITE (io_unit, '(F10.4, 3(1X, ES12.4E3))') k_pos(ikpt), berry_k(:, iw, ikpt)
      END DO
      WRITE (io_unit, *) ! blank line
    END DO
    WRITE (io_unit, *) ! blank line
    CLOSE (io_unit)
  END SUBROUTINE write_Berry_k
  !
  SUBROUTINE write_shift(fname, hw, sigma_w)
    USE system, ONLY: dim
    CHARACTER(LEN=*), INTENT(IN) :: fname
    REAL(DP), INTENT(IN) :: hw(:)
    REAL(DP), INTENT(IN) :: sigma_w(3, 6, SIZE(hw))
    INTEGER :: io_unit, iom, ia
    CHARACTER(LEN=20) :: unit_str
    CHARACTER(LEN=256) :: fname_a
    CHARACTER(LEN=1), PARAMETER :: a_lab(3) = (/'x', 'y', 'z'/)
    IF (.NOT. ionode) RETURN
    !
    SELECT CASE (dim)
    CASE (3)
      unit_str = 'muA/V^2'
    CASE (2)
      unit_str = 'muA*Ang/V^2'
    CASE (1)
      unit_str = 'muA*Ang^2/V^2'
    END SELECT
    !
    DO ia = 1, 3
      fname_a = TRIM(fname)//'_'//a_lab(ia)//'.dat'
      io_unit = get_free_unit()
      OPEN (unit=io_unit, file=fname_a)
      CALL writing_info('shift current', fname_a)
      WRITE (io_unit, 0947) 'shift current units: '//TRIM(unit_str)
      WRITE (io_unit, 0947) 'current direction: '//a_lab(ia)
      WRITE (io_unit, 0947) 'hw (eV), xx, xy, yy, yz, zz, zx'

      DO iom = 1, SIZE(hw)
        WRITE (io_unit, '(F13.6, 6(1X, ES16.8E3))') hw(iom), sigma_w(ia, :, iom)
      END DO
      CLOSE (io_unit)
    END DO
0947 FORMAT("# ", A)
  END SUBROUTINE write_shift

END MODULE io_output
