MODULE w90_chk_dump
  USE kinds, ONLY: DP
  USE io_w90, ONLY: w90data
  IMPLICIT NONE
  PRIVATE
  PUBLIC::write_w90_chk_dump
  LOGICAL, PUBLIC::chk_dump = .FALSE.
  CHARACTER(LEN=256), PUBLIC::chk_dump_file = 'itg_input_check_dump.txt'
CONTAINS
  SUBROUTINE write_w90_chk_dump()
    USE io_global, ONLY: ionode, get_free_unit
    INTEGER :: iu, ios, excl_sum
    CHARACTER(LEN=256) :: fbase
    REAL(DP), ALLOCATABLE :: rv(:)
    COMPLEX(DP), ALLOCATABLE :: cv(:)
    INTEGER, ALLOCATABLE :: iv(:)
    LOGICAL, ALLOCATABLE :: lv(:)

    IF (.NOT. ionode) RETURN
    fbase = TRIM(chk_dump_file)
    IF (LEN_TRIM(fbase) == 0) fbase = 'itg_input_check_dump.txt'

    iu = get_free_unit()
    OPEN (unit=iu, file=TRIM(fbase), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'write_w90_chk_dump', 'Failed to open dump file: '//TRIM(fbase))

    IF (ALLOCATED(w90data%excl_bands)) THEN
      excl_sum = SUM(w90data%excl_bands)
    ELSE
      excl_sum = 0
    END IF

    WRITE (iu, '(A)') 'header='//TRIM(w90data%chk_header)
    WRITE (iu, '(A)') 'checkpoint='//TRIM(w90data%chk_checkpoint)
    WRITE (iu, '(A,I0)') 'nbnd=', w90data%nbnd
    WRITE (iu, '(A,I0)') 'nbnd_excl=', w90data%nbnd_excl
    WRITE (iu, '(A,I0)') 'nkpt=', w90data%nkpt
    WRITE (iu, '(A,I0)') 'nnb=', w90data%nnb
    WRITE (iu, '(A,I0)') 'Nw=', w90data%Nw
    WRITE (iu, '(A,3(1X,I0))') 'k_grid=', w90data%k_grid
    WRITE (iu, '(A,L1)') 'have_disentangled=', w90data%have_disentangled
    WRITE (iu, '(A,1X,ES24.16E3)') 'omega_invariant=', w90data%omega_invariant
    WRITE (iu, '(A,9(1X,ES24.16E3))') 'real_lattice=', w90data%real_lattice
    WRITE (iu, '(A,9(1X,ES24.16E3))') 'recip_lattice=', w90data%recip_lattice
    WRITE (iu, '(A,9(1X,ES24.16E3))') 'recip_lattice_inv=', w90data%recip_lattice_inv
    WRITE (iu, '(A,I0)') 'excl_bands_sum=', excl_sum
    WRITE (iu, '(A,1X,I0)') 'dims_excl_bands=', w90data%nbnd_excl
    WRITE (iu, '(A,2(1X,I0))') 'dims_k_cart=', 3, w90data%nkpt
    WRITE (iu, '(A,2(1X,I0))') 'dims_k_red=', 3, w90data%nkpt
    WRITE (iu, '(A,2(1X,I0))') 'dims_eigval=', w90data%nbnd, w90data%nkpt
    IF (w90data%have_disentangled) THEN
      WRITE (iu, '(A,2(1X,I0))') 'dims_lwindow=', w90data%nbnd, w90data%nkpt
      WRITE (iu, '(A,1X,I0)') 'dims_ndimwin=', w90data%nkpt
      WRITE (iu, '(A,3(1X,I0))') 'dims_u_matrix_opt=', w90data%nbnd, w90data%Nw, w90data%nkpt
    ELSE
      WRITE (iu, '(A,2(1X,I0))') 'dims_lwindow=', 0, 0
      WRITE (iu, '(A,1X,I0)') 'dims_ndimwin=', 0
      WRITE (iu, '(A,3(1X,I0))') 'dims_u_matrix_opt=', 0, 0, 0
    END IF
    WRITE (iu, '(A,3(1X,I0))') 'dims_u_matrix=', w90data%Nw, w90data%Nw, w90data%nkpt
    WRITE (iu, '(A,4(1X,I0))') 'dims_m_matrix=', w90data%Nw, w90data%Nw, w90data%nnb, w90data%nkpt
    WRITE (iu, '(A,2(1X,I0))') 'dims_wannier_centers=', 3, w90data%Nw
    WRITE (iu, '(A,1X,I0)') 'dims_wannier_spreads=', w90data%Nw
    CLOSE (iu)

    IF (ALLOCATED(w90data%excl_bands)) THEN
      ALLOCATE (iv(SIZE(w90data%excl_bands)))
      iv = RESHAPE(w90data%excl_bands, [SIZE(w90data%excl_bands)])
    ELSE
      ALLOCATE (iv(0))
    END IF
    CALL dump_i(TRIM(fbase)//'.excl_bands', iv)
    DEALLOCATE (iv)

    ALLOCATE (rv(SIZE(w90data%k_cart)))
    rv = RESHAPE(w90data%k_cart, [SIZE(w90data%k_cart)])
    CALL dump_r(TRIM(fbase)//'.k_cart', rv)
    DEALLOCATE (rv)

    ALLOCATE (rv(SIZE(w90data%k_red)))
    rv = RESHAPE(w90data%k_red, [SIZE(w90data%k_red)])
    CALL dump_r(TRIM(fbase)//'.k_red', rv)
    DEALLOCATE (rv)

    ALLOCATE (rv(SIZE(w90data%eigval)))
    rv = RESHAPE(w90data%eigval, [SIZE(w90data%eigval)])
    CALL dump_r(TRIM(fbase)//'.eigval', rv)
    DEALLOCATE (rv)

    IF (ALLOCATED(w90data%lwindow)) THEN
      ALLOCATE (lv(SIZE(w90data%lwindow)))
      lv = RESHAPE(w90data%lwindow, [SIZE(w90data%lwindow)])
    ELSE
      ALLOCATE (lv(0))
    END IF
    CALL dump_l(TRIM(fbase)//'.lwindow', lv)
    DEALLOCATE (lv)

    IF (ALLOCATED(w90data%ndimwin)) THEN
      ALLOCATE (iv(SIZE(w90data%ndimwin)))
      iv = RESHAPE(w90data%ndimwin, [SIZE(w90data%ndimwin)])
    ELSE
      ALLOCATE (iv(0))
    END IF
    CALL dump_i(TRIM(fbase)//'.ndimwin', iv)
    DEALLOCATE (iv)

    IF (ALLOCATED(w90data%u_matrix_opt)) THEN
      ALLOCATE (cv(SIZE(w90data%u_matrix_opt)))
      cv = RESHAPE(w90data%u_matrix_opt, [SIZE(w90data%u_matrix_opt)])
    ELSE
      ALLOCATE (cv(0))
    END IF
    CALL dump_c(TRIM(fbase)//'.u_matrix_opt', cv)
    DEALLOCATE (cv)

    ALLOCATE (cv(SIZE(w90data%u_matrix)))
    cv = RESHAPE(w90data%u_matrix, [SIZE(w90data%u_matrix)])
    CALL dump_c(TRIM(fbase)//'.u_matrix', cv)
    DEALLOCATE (cv)

    ALLOCATE (cv(SIZE(w90data%m_matrix)))
    cv = RESHAPE(w90data%m_matrix, [SIZE(w90data%m_matrix)])
    CALL dump_c(TRIM(fbase)//'.m_matrix', cv)
    DEALLOCATE (cv)

    ALLOCATE (rv(SIZE(w90data%wannier_centers)))
    rv = RESHAPE(w90data%wannier_centers, [SIZE(w90data%wannier_centers)])
    CALL dump_r(TRIM(fbase)//'.wannier_centers', rv)
    DEALLOCATE (rv)

    ALLOCATE (rv(SIZE(w90data%wannier_spreads)))
    rv = RESHAPE(w90data%wannier_spreads, [SIZE(w90data%wannier_spreads)])
    CALL dump_r(TRIM(fbase)//'.wannier_spreads', rv)
    DEALLOCATE (rv)
  END SUBROUTINE write_w90_chk_dump

  SUBROUTINE dump_r(fname, vec)
    USE io_global, ONLY: get_free_unit
    CHARACTER(LEN=*), INTENT(IN) :: fname
    REAL(DP), INTENT(IN) :: vec(:)
    INTEGER :: iu, ios, i
    iu = get_free_unit()
    OPEN (unit=iu, file=TRIM(fname), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'dump_r', 'Failed to open '//TRIM(fname))
    WRITE (iu, '(I0)') SIZE(vec)
    DO i = 1, SIZE(vec)
      WRITE (iu, '(ES24.16E3)') vec(i)
    END DO
    CLOSE (iu)
  END SUBROUTINE dump_r

  SUBROUTINE dump_c(fname, vec)
    USE io_global, ONLY: get_free_unit
    CHARACTER(LEN=*), INTENT(IN) :: fname
    COMPLEX(DP), INTENT(IN) :: vec(:)
    INTEGER :: iu, ios, i
    iu = get_free_unit()
    OPEN (unit=iu, file=TRIM(fname), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'dump_c', 'Failed to open '//TRIM(fname))
    WRITE (iu, '(I0)') SIZE(vec)
    DO i = 1, SIZE(vec)
      WRITE (iu, '(2(1X,ES24.16E3))') REAL(vec(i), DP), AIMAG(vec(i))
    END DO
    CLOSE (iu)
  END SUBROUTINE dump_c

  SUBROUTINE dump_i(fname, vec)
    USE io_global, ONLY: get_free_unit
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: vec(:)
    INTEGER :: iu, ios, i
    iu = get_free_unit()
    OPEN (unit=iu, file=TRIM(fname), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'dump_i', 'Failed to open '//TRIM(fname))
    WRITE (iu, '(I0)') SIZE(vec)
    DO i = 1, SIZE(vec)
      WRITE (iu, '(I0)') vec(i)
    END DO
    CLOSE (iu)
  END SUBROUTINE dump_i

  SUBROUTINE dump_l(fname, vec)
    USE io_global, ONLY: get_free_unit
    CHARACTER(LEN=*), INTENT(IN) :: fname
    LOGICAL, INTENT(IN) :: vec(:)
    INTEGER :: iu, ios, i, iv
    iu = get_free_unit()
    OPEN (unit=iu, file=TRIM(fname), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'dump_l', 'Failed to open '//TRIM(fname))
    WRITE (iu, '(I0)') SIZE(vec)
    DO i = 1, SIZE(vec)
      IF (vec(i)) THEN
        iv = 1
      ELSE
        iv = 0
      END IF
      WRITE (iu, '(I0)') iv
    END DO
    CLOSE (iu)
  END SUBROUTINE dump_l
END MODULE w90_chk_dump
