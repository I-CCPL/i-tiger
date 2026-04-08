SUBMODULE(wannier90) w90_chk_dum
CONTAINS
  MODULE SUBROUTINE clear_chk_dum(self)
    CLASS(chk_dum_type), INTENT(INOUT) :: self
    IF (ALLOCATED(self%excl_bands)) DEALLOCATE (self%excl_bands)
    IF (ALLOCATED(self%lwindow)) DEALLOCATE (self%lwindow)
    IF (ALLOCATED(self%u_matrix_opt)) DEALLOCATE (self%u_matrix_opt)
    IF (ALLOCATED(self%u_matrix)) DEALLOCATE (self%u_matrix)
    IF (ALLOCATED(self%m_matrix)) DEALLOCATE (self%m_matrix)
  END SUBROUTINE clear_chk_dum

  MODULE SUBROUTINE write_chk_dump(self, chk_dum)
    USE dump_vec_io, ONLY: dump_r, dump_c, dump_i, dump_l
    USE system, ONLY: Nw, real_lattice, recip_lattice
    TYPE(w90data_type), INTENT(IN) :: self
    TYPE(chk_dum_type), INTENT(IN) :: chk_dum
    INTEGER :: io_unit, ios, excl_sum
    CHARACTER(LEN=256) :: fbase
    REAL(DP), ALLOCATABLE :: rv(:)
    COMPLEX(DP), ALLOCATABLE :: cv(:)
    INTEGER, ALLOCATABLE :: iv(:)
    LOGICAL, ALLOCATABLE :: lv(:)

    fbase = TRIM(self%prefix)//'.chk_dump'
    IF (LEN_TRIM(fbase) == 0) fbase = 'itg_input_check_dump.txt'

    io_unit = get_free_unit()
    OPEN (unit=io_unit, file=TRIM(fbase), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'write_chk_dump', 'Failed to open dump file: '//TRIM(fbase))

    IF (SIZE(chk_dum%excl_bands) > 0) THEN
      excl_sum = SUM(chk_dum%excl_bands)
    ELSE
      excl_sum = 0
    END IF

    WRITE (io_unit, '(A)') 'header='//TRIM(chk_dum%header)
    WRITE (io_unit, '(A)') 'checkpoint='//TRIM(chk_dum%checkpoint)
    WRITE (io_unit, '(A,I0)') 'nbnd=', self%nbnd
    WRITE (io_unit, '(A,I0)') 'nbnd_excl=', chk_dum%nbnd_excl
    WRITE (io_unit, '(A,I0)') 'nkpt=', self%kpts%nkpt
    WRITE (io_unit, '(A,I0)') 'nnb=', self%nnb
    WRITE (io_unit, '(A,I0)') 'Nw=', Nw
    WRITE (io_unit, '(A,3(1X,I0))') 'k_grid=', self%k_grid
    WRITE (io_unit, '(A,L1)') 'have_disentangled=', chk_dum%have_disentangled
    WRITE (io_unit, '(A,1X,ES24.16E3)') 'omega_invariant=', chk_dum%omega_invariant
    WRITE (io_unit, '(A,9(1X,ES24.16E3))') 'real_lattice=', real_lattice
    WRITE (io_unit, '(A,9(1X,ES24.16E3))') 'recip_lattice=', recip_lattice
    WRITE (io_unit, '(A,I0)') 'excl_bands_sum=', excl_sum
    WRITE (io_unit, '(A,1X,I0)') 'dims_excl_bands=', SIZE(chk_dum%excl_bands)
    WRITE (io_unit, '(A,2(1X,I0))') 'dims_k_cart=', 3, self%kpts%nkpt
    WRITE (io_unit, '(A,2(1X,I0))') 'dims_k_red=', 3, self%kpts%nkpt
    WRITE (io_unit, '(A,2(1X,I0))') 'dims_eigval=', self%nbnd, self%kpts%nkpt
    IF (chk_dum%have_disentangled) THEN
      WRITE (io_unit, '(A,2(1X,I0))') 'dims_lwindow=', SIZE(chk_dum%lwindow, 1), SIZE(chk_dum%lwindow, 2)
      WRITE (io_unit, '(A,1X,I0)') 'dims_ndimwin=', SIZE(self%ndimwin)
      WRITE (io_unit, '(A,3(1X,I0))') 'dims_u_matrix_opt=', SIZE(chk_dum%u_matrix_opt, 1), SIZE(chk_dum%u_matrix_opt, 2), SIZE(chk_dum%u_matrix_opt, 3)
    ELSE
      WRITE (io_unit, '(A,2(1X,I0))') 'dims_lwindow=', 0, 0
      WRITE (io_unit, '(A,1X,I0)') 'dims_ndimwin=', 0
      WRITE (io_unit, '(A,3(1X,I0))') 'dims_u_matrix_opt=', 0, 0, 0
    END IF
   WRITE (io_unit, '(A,3(1X,I0))') 'dims_u_matrix=', SIZE(chk_dum%u_matrix, 1), SIZE(chk_dum%u_matrix, 2), SIZE(chk_dum%u_matrix, 3)
    WRITE (io_unit, '(A,4(1X,I0))') 'dims_m_matrix=', SIZE(chk_dum%m_matrix, 1), SIZE(chk_dum%m_matrix, 2), SIZE(chk_dum%m_matrix, 3), SIZE(chk_dum%m_matrix, 4)
    WRITE (io_unit, '(A,2(1X,I0))') 'dims_wannier_centers=', 3, Nw
    WRITE (io_unit, '(A,1X,I0)') 'dims_wannier_spreads=', Nw
    CLOSE (io_unit)

    ALLOCATE (iv(SIZE(chk_dum%excl_bands)))
    iv = RESHAPE(chk_dum%excl_bands, [SIZE(chk_dum%excl_bands)])
    CALL dump_i(TRIM(fbase)//'.excl_bands', iv)
    DEALLOCATE (iv)

    ALLOCATE (rv(SIZE(self%kpts%k_cart)))
    rv = RESHAPE(self%kpts%k_cart, [SIZE(self%kpts%k_cart)])
    CALL dump_r(TRIM(fbase)//'.k_cart', rv)
    DEALLOCATE (rv)

    ALLOCATE (rv(SIZE(self%kpts%k_red)))
    rv = RESHAPE(self%kpts%k_red, [SIZE(self%kpts%k_red)])
    CALL dump_r(TRIM(fbase)//'.k_red', rv)
    DEALLOCATE (rv)

    ALLOCATE (rv(SIZE(self%eigval)))
    rv = RESHAPE(self%eigval, [SIZE(self%eigval)])
    CALL dump_r(TRIM(fbase)//'.eigval', rv)
    DEALLOCATE (rv)

    ALLOCATE (lv(SIZE(chk_dum%lwindow)))
    lv = RESHAPE(chk_dum%lwindow, [SIZE(chk_dum%lwindow)])
    CALL dump_l(TRIM(fbase)//'.lwindow', lv)
    DEALLOCATE (lv)

    ALLOCATE (iv(SIZE(self%ndimwin)))
    iv = RESHAPE(self%ndimwin, [SIZE(self%ndimwin)])
    CALL dump_i(TRIM(fbase)//'.ndimwin', iv)
    DEALLOCATE (iv)

    ALLOCATE (cv(SIZE(chk_dum%u_matrix_opt)))
    cv = RESHAPE(chk_dum%u_matrix_opt, [SIZE(chk_dum%u_matrix_opt)])
    CALL dump_c(TRIM(fbase)//'.u_matrix_opt', cv)
    DEALLOCATE (cv)

    ALLOCATE (cv(SIZE(chk_dum%u_matrix)))
    cv = RESHAPE(chk_dum%u_matrix, [SIZE(chk_dum%u_matrix)])
    CALL dump_c(TRIM(fbase)//'.u_matrix', cv)
    DEALLOCATE (cv)

    ALLOCATE (cv(SIZE(chk_dum%m_matrix)))
    cv = RESHAPE(chk_dum%m_matrix, [SIZE(chk_dum%m_matrix)])
    CALL dump_c(TRIM(fbase)//'.m_matrix', cv)
    DEALLOCATE (cv)

    ALLOCATE (rv(SIZE(self%wannier_center_cart)))
    rv = RESHAPE(self%wannier_center_cart, [SIZE(self%wannier_center_cart)])
    CALL dump_r(TRIM(fbase)//'.wannier_centers', rv)
    DEALLOCATE (rv)

    ALLOCATE (rv(SIZE(self%wannier_spread)))
    rv = RESHAPE(self%wannier_spread, [SIZE(self%wannier_spread)])
    CALL dump_r(TRIM(fbase)//'.wannier_spreads', rv)
    DEALLOCATE (rv)
  END SUBROUTINE write_chk_dump

END SUBMODULE
