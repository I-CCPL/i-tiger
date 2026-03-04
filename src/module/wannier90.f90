MODULE wannier90
  USE kinds, ONLY: DP
  USE io_global, ONLY: ionode, stdout, get_free_unit
  IMPLICIT NONE
  PRIVATE
  PUBLIC::read_w90, w90data_type, w90data
  LOGICAL, PUBLIC::chk_w90 = .FALSE.
  LOGICAL, PUBLIC::Hq_band = .FALSE.
  !
  TYPE::w90data_type
    CHARACTER(LEN=256) :: prefix

    INTEGER::nbnd
    !< Number of bands. \
    !< nbnd==Nw for non-disentangled case \
    !< nbnd>=Nw for disentangled case.
    INTEGER :: nkpt
    !< Number of k-points
    INTEGER :: nnb
    !< Number of nearest neighbor k-points
    INTEGER :: k_grid(3)
    !< Number of k-points along each reciprocal lattice vector \
    !< (Monkhorst-Pack grid)
    REAL(DP), ALLOCATABLE :: k_cart(:, :)
    !< k-points in Cartesian coordinates (3, nkpt)
    REAL(DP), ALLOCATABLE :: k_red(:, :)
    !< k-points in reduced coordinates (3, nkpt)

    REAL(DP), ALLOCATABLE :: wannier_center_cart(:, :)
    !< Wannier center in Cartesian coordinates (3, Nw)
    REAL(DP), ALLOCATABLE :: wannier_spread(:)
    !< (Nw)

    COMPLEX(DP), ALLOCATABLE :: v_matrix(:, :, :)
    !< (nbnd, Nw, nkpt) for disentangled case
    REAL(DP), ALLOCATABLE::eigval(:, :)
    !< Eigenvalues (nbnd, nkpt)
    COMPLEX(DP), ALLOCATABLE :: eigvec(:, :, :)
    !< Unitary matrix corresponding eigenvector \
    !< (Nw, Nw, nkpt)
    COMPLEX(DP), ALLOCATABLE :: Hq(:, :, :)
    !< Hamiltonian in Wannier gauge (Nw, Nw, nkpt)

  CONTAINS
    PROCEDURE::read_eig => read_w90_eig
    PROCEDURE::read_chk => read_w90_chk
    PROCEDURE::clear => clear_w90_data
    PROCEDURE::build_Hq => build_w90_Hq
  END TYPE w90data_type
  TYPE(w90data_type)::w90data

  TYPE::chk_dum_type
    CHARACTER(LEN=33) :: header
    INTEGER :: nbnd_excl
    INTEGER, ALLOCATABLE :: excl_bands(:)
    CHARACTER(LEN=20) :: checkpoint
    LOGICAL :: have_disentangled
    REAL(DP) :: omega_invariant
    LOGICAL, ALLOCATABLE :: lwindow(:, :)
    INTEGER, ALLOCATABLE :: ndimwin(:)
    COMPLEX(DP), ALLOCATABLE :: u_matrix_opt(:, :, :)
    COMPLEX(DP), ALLOCATABLE::u_matrix(:, :, :)
    COMPLEX(DP), ALLOCATABLE::m_matrix(:, :, :, :)
  CONTAINS
    PROCEDURE::clear => clear_chk_dum
  END TYPE chk_dum_type
CONTAINS
  SUBROUTINE read_w90()
    TYPE(chk_dum_type) :: chk_dum
    !
    CALL w90data%clear()
    CALL w90data%read_chk(chk_dum)
    CALL w90data%read_eig()

    IF (ionode .AND. chk_w90) THEN
      CALL write_chk_dump(w90data, chk_dum)
    END IF

    CALL chk_dum%clear()
  END SUBROUTINE read_w90
  !
  SUBROUTINE read_w90_chk(self, chk_dum)
    !< Ref. wannier90/src/wannier90_readwrite.F90
    USE io_global, ONLY: check_file
    USE mp_base, ONLY: mp_bcast
    USE system, ONLY: Nw, cell_setup, real_lattice, recip_lattice, &
                      red2cart_recip
    CLASS(w90data_type), INTENT(INOUT) :: self
    TYPE(chk_dum_type), INTENT(OUT) :: chk_dum
    !
    INTEGER::io_unit, ios, ikpt
    INTEGER::ndw
    !
    WRITE (stdout, '(2X, A)') 'Reading .chk file...'
    CALL check_file(TRIM(self%prefix)//'.chk')
    IF (ionode) THEN
      io_unit = get_free_unit()
      OPEN (unit=io_unit, file=TRIM(self%prefix)//'.chk', form='unformatted', action='read', iostat=ios)
      CALL errore(ios, 'read_w90_chk', 'Failed to open '//TRIM(self%prefix)//'.chk')

      READ (io_unit) chk_dum%header
      WRITE (stdout, '(2X, A)') '- header: '//TRIM(chk_dum%header)
      READ (io_unit) self%nbnd
      WRITE (stdout, '(2X, A, I0)') '- nbnd: ', self%nbnd
      READ (io_unit) chk_dum%nbnd_excl
      WRITE (stdout, '(2X, A, I0)') '- nbnd_excl: ', chk_dum%nbnd_excl
      IF (chk_dum%nbnd_excl > 0) THEN
        ALLOCATE (chk_dum%excl_bands(chk_dum%nbnd_excl))
        READ (io_unit) chk_dum%excl_bands
      ELSE
        ALLOCATE (chk_dum%excl_bands(0))
        ! .chk always writes the excluded-band record; consume empty record to keep alignment.
        READ (io_unit)
      END IF

      READ (io_unit) real_lattice
      WRITE (stdout, '(2X, A, 3F12.6)') '- real_lattice:', real_lattice(:, 1)
      WRITE (stdout, '(2X, A, 3F12.6)') '               ', real_lattice(:, 2)
      WRITE (stdout, '(2X, A, 3F12.6)') '               ', real_lattice(:, 3)
      READ (io_unit) recip_lattice
      WRITE (stdout, '(2X, A, 3F12.6)') '- recip_lattice:', recip_lattice(:, 1)
      WRITE (stdout, '(2X, A, 3F12.6)') '                ', recip_lattice(:, 2)
      WRITE (stdout, '(2X, A, 3F12.6)') '                ', recip_lattice(:, 3)
      CALL cell_setup()

      READ (io_unit) self%nkpt
      WRITE (stdout, '(2X, A, I0)') '- nkpt: ', self%nkpt
      READ (io_unit) self%k_grid
      WRITE (stdout, '(2X, A, 3(1X,I0))') '- k_grid: ', self%k_grid
      ALLOCATE (self%k_cart(3, self%nkpt))
      ALLOCATE (self%k_red(3, self%nkpt))
      READ (io_unit) self%k_red
      CALL red2cart_recip(self%k_red, self%k_cart, self%nkpt)
      READ (io_unit) self%nnb
      WRITE (stdout, '(2X, A, I0)') '- nnb: ', self%nnb
      READ (io_unit) Nw
      WRITE (stdout, '(2X, A, I0)') '- Nw: ', Nw

      READ (io_unit) chk_dum%checkpoint
      WRITE (stdout, '(2X, A)') '- checkpoint: '//TRIM(chk_dum%checkpoint)
      READ (io_unit) chk_dum%have_disentangled
      WRITE (stdout, '(2X, A, L1)') '- have_disentangled: ', chk_dum%have_disentangled
      IF (chk_dum%have_disentangled) THEN
        READ (io_unit) chk_dum%omega_invariant
        WRITE (stdout, '(2X, A, F12.6)') '- omega_invariant: ', chk_dum%omega_invariant
        ALLOCATE (chk_dum%lwindow(self%nbnd, self%nkpt))
        ALLOCATE (chk_dum%ndimwin(self%nkpt))
        ALLOCATE (chk_dum%u_matrix_opt(self%nbnd, Nw, self%nkpt))
        READ (io_unit) chk_dum%lwindow
        READ (io_unit) chk_dum%ndimwin
        READ (io_unit) chk_dum%u_matrix_opt
      ELSE
        chk_dum%omega_invariant = 0.0_DP
        ALLOCATE (chk_dum%lwindow(0, 0))
        ALLOCATE (chk_dum%ndimwin(0))
        ALLOCATE (chk_dum%u_matrix_opt(0, 0, 0))
      END IF

      ALLOCATE (chk_dum%u_matrix(Nw, Nw, self%nkpt))
      ALLOCATE (chk_dum%m_matrix(Nw, Nw, self%nnb, self%nkpt))
      READ (io_unit) chk_dum%u_matrix
      READ (io_unit) chk_dum%m_matrix

      IF (chk_dum%have_disentangled) THEN
        ALLOCATE (self%v_matrix(self%nbnd, Nw, self%nkpt))
        DO ikpt = 1, self%nkpt
          ndw = chk_dum%ndimwin(ikpt)
          self%v_matrix(1:ndw, :, ikpt) = MATMUL(chk_dum%u_matrix_opt(1:ndw, :, ikpt), chk_dum%u_matrix(:, :, ikpt))
        END DO
      ELSE
        ! TODO: Check it is consistent with the non-disentangled case.
        IF (Nw /= self%nbnd) THEN
          CALL errore(1, 'read_w90_chk', 'Nw must equal nbnd for non-disentangled case')
        END IF
        !
        ALLOCATE (self%v_matrix(Nw, Nw, self%nkpt))
        DO ikpt = 1, self%nkpt
          self%v_matrix(:, :, ikpt) = chk_dum%u_matrix(:, :, ikpt)
        END DO
      END IF

      ALLOCATE (self%wannier_center_cart(3, Nw))
      ALLOCATE (self%wannier_spread(Nw))
      READ (io_unit) self%wannier_center_cart
      READ (io_unit) self%wannier_spread

      CLOSE (io_unit)
    END IF

    CALL mp_bcast(self%nbnd)
    CALL mp_bcast(real_lattice)
    CALL mp_bcast(recip_lattice)
    CALL mp_bcast(self%nkpt)
    CALL mp_bcast(self%k_grid)
    CALL mp_bcast(self%nnb)
    CALL mp_bcast(Nw)

    IF (.NOT. ionode) THEN
      CALL cell_setup()
      ALLOCATE (self%k_cart(3, self%nkpt))
      ALLOCATE (self%k_red(3, self%nkpt))
      ALLOCATE (self%v_matrix(self%nbnd, Nw, self%nkpt))
      ALLOCATE (self%wannier_center_cart(3, Nw))
      ALLOCATE (self%wannier_spread(Nw))
    END IF

    CALL mp_bcast(self%k_cart)
    CALL mp_bcast(self%k_red)
    CALL mp_bcast(self%v_matrix)
    CALL mp_bcast(self%wannier_center_cart)
    CALL mp_bcast(self%wannier_spread)
  END SUBROUTINE read_w90_chk
  !
  SUBROUTINE read_w90_eig(self)
    !< Ref. wannier90/src/readwrite.F90
    USE io_global, ONLY: check_file
    USE mp_base, ONLY: mp_bcast
    CLASS(w90data_type), INTENT(INOUT) :: self
    INTEGER::io_unit, ikpt, ibnd, jkpt, jbnd, ios
    CHARACTER(LEN=256)::msg
    !
    WRITE (stdout, '(2X, A)') 'Reading .eig file...'
    CALL check_file(TRIM(self%prefix)//'.eig')
    IF (ionode) THEN
      io_unit = get_free_unit()
      OPEN (unit=io_unit, file=TRIM(self%prefix)//'.eig', form='formatted', action='read', iostat=ios)
      CALL errore(ios, 'read_w90_eig', 'Failed to open '//TRIM(self%prefix)//'.eig')
      ALLOCATE (self%eigval(self%nbnd, self%nkpt))

      DO ikpt = 1, self%nkpt
        DO ibnd = 1, self%nbnd
          READ (io_unit, *, iostat=ios) jbnd, jkpt, self%eigval(ibnd, ikpt)
          IF (ios /= 0) THEN
            WRITE (msg, '(A,1X,I0,1X,A,1X,I0)') 'Failed reading .eig at ibnd=', ibnd, 'ikpt=', ikpt
            CALL errore(ios, 'read_w90_eig', TRIM(msg))
          END IF
          IF (jbnd /= ibnd .OR. jkpt /= ikpt) THEN
            WRITE (msg, '(A,1X,I0,A,I0,A,I0,A,I0,A)') &
              'Unexpected (.eig) indices. Expected (', ibnd, ',', ikpt, '), got (', jbnd, ',', jkpt, ').'
            CALL errore(1, 'read_w90_eig', TRIM(msg))
          END IF
        END DO
      END DO
      CLOSE (io_unit)
    ELSE
      ALLOCATE (self%eigval(self%nbnd, self%nkpt))
    END IF
    IF (self%nbnd > 0 .AND. self%nkpt > 0) CALL mp_bcast(self%eigval)
  END SUBROUTINE read_w90_eig

  SUBROUTINE write_chk_dump(self, chk_dum)
    USE dump_vec_io, ONLY: dump_r, dump_c, dump_i, dump_l
    USE system, ONLY: Nw, real_lattice, recip_lattice
    CLASS(w90data_type), INTENT(IN) :: self
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
    WRITE (io_unit, '(A,I0)') 'nkpt=', self%nkpt
    WRITE (io_unit, '(A,I0)') 'nnb=', self%nnb
    WRITE (io_unit, '(A,I0)') 'Nw=', Nw
    WRITE (io_unit, '(A,3(1X,I0))') 'k_grid=', self%k_grid
    WRITE (io_unit, '(A,L1)') 'have_disentangled=', chk_dum%have_disentangled
    WRITE (io_unit, '(A,1X,ES24.16E3)') 'omega_invariant=', chk_dum%omega_invariant
    WRITE (io_unit, '(A,9(1X,ES24.16E3))') 'real_lattice=', real_lattice
    WRITE (io_unit, '(A,9(1X,ES24.16E3))') 'recip_lattice=', recip_lattice
    WRITE (io_unit, '(A,I0)') 'excl_bands_sum=', excl_sum
    WRITE (io_unit, '(A,1X,I0)') 'dims_excl_bands=', SIZE(chk_dum%excl_bands)
    WRITE (io_unit, '(A,2(1X,I0))') 'dims_k_cart=', 3, self%nkpt
    WRITE (io_unit, '(A,2(1X,I0))') 'dims_k_red=', 3, self%nkpt
    WRITE (io_unit, '(A,2(1X,I0))') 'dims_eigval=', self%nbnd, self%nkpt
    IF (chk_dum%have_disentangled) THEN
      WRITE (io_unit, '(A,2(1X,I0))') 'dims_lwindow=', SIZE(chk_dum%lwindow, 1), SIZE(chk_dum%lwindow, 2)
      WRITE (io_unit, '(A,1X,I0)') 'dims_ndimwin=', SIZE(chk_dum%ndimwin)
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

    ALLOCATE (rv(SIZE(self%k_cart)))
    rv = RESHAPE(self%k_cart, [SIZE(self%k_cart)])
    CALL dump_r(TRIM(fbase)//'.k_cart', rv)
    DEALLOCATE (rv)

    ALLOCATE (rv(SIZE(self%k_red)))
    rv = RESHAPE(self%k_red, [SIZE(self%k_red)])
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

    ALLOCATE (iv(SIZE(chk_dum%ndimwin)))
    iv = RESHAPE(chk_dum%ndimwin, [SIZE(chk_dum%ndimwin)])
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

  SUBROUTINE clear_w90_data(self)
    CLASS(w90data_type), INTENT(INOUT) :: self
    self%nbnd = 0
    self%nkpt = 0
    self%k_grid = 0
    self%nnb = 0

    IF (ALLOCATED(self%eigval)) DEALLOCATE (self%eigval)
    IF (ALLOCATED(self%k_cart)) DEALLOCATE (self%k_cart)
    IF (ALLOCATED(self%k_red)) DEALLOCATE (self%k_red)
    IF (ALLOCATED(self%v_matrix)) DEALLOCATE (self%v_matrix)
    IF (ALLOCATED(self%wannier_center_cart)) DEALLOCATE (self%wannier_center_cart)
    IF (ALLOCATED(self%wannier_spread)) DEALLOCATE (self%wannier_spread)
    IF (ALLOCATED(self%eigvec)) DEALLOCATE (self%eigvec)
    IF (ALLOCATED(self%Hq)) DEALLOCATE (self%Hq)
  END SUBROUTINE clear_w90_data

  SUBROUTINE clear_chk_dum(chk_dum)
    CLASS(chk_dum_type), INTENT(INOUT) :: chk_dum
    IF (ALLOCATED(chk_dum%excl_bands)) DEALLOCATE (chk_dum%excl_bands)
    IF (ALLOCATED(chk_dum%lwindow)) DEALLOCATE (chk_dum%lwindow)
    IF (ALLOCATED(chk_dum%ndimwin)) DEALLOCATE (chk_dum%ndimwin)
    IF (ALLOCATED(chk_dum%u_matrix_opt)) DEALLOCATE (chk_dum%u_matrix_opt)
    IF (ALLOCATED(chk_dum%u_matrix)) DEALLOCATE (chk_dum%u_matrix)
    IF (ALLOCATED(chk_dum%m_matrix)) DEALLOCATE (chk_dum%m_matrix)
  END SUBROUTINE clear_chk_dum

  SUBROUTINE build_w90_Hq(self)
    !< Build Hamiltonian in q-space
    USE io_global, ONLY: write_sep_line
    USE lin_eig_H, ONLY: write_band
    USE mp_base, ONLY: mp_bcast
    USE system, ONLY: Nw
    CLASS(w90data_type), INTENT(INOUT) :: self
    INTEGER::ikpt, iw, jw, ibnd
    COMPLEX(DP)::hval
    REAL(DP)::herm_abs_max, h_abs_max, herm_rel
    CHARACTER(LEN=256)::msg
    !
    WRITE (stdout, '(2X, A)') 'Building H(q) in Wannier gauge...'
    ALLOCATE (self%Hq(Nw, Nw, self%nkpt))
    !
    IF (ionode) THEN
      self%Hq = CMPLX(0.0_DP, 0.0_DP, DP)

      DO ikpt = 1, self%nkpt
        DO iw = 1, Nw
          DO jw = 1, Nw
            hval = CMPLX(0.0_DP, 0.0_DP, DP)
            DO ibnd = 1, self%nbnd
              hval = hval + CONJG(self%v_matrix(ibnd, iw, ikpt))*self%eigval(ibnd, ikpt)* &
                     self%v_matrix(ibnd, jw, ikpt)
            END DO
            self%Hq(iw, jw, ikpt) = hval
          END DO
        END DO
      END DO
      !
      CALL check_hermiticity(Nw, self%nkpt, self%Hq, herm_abs_max, h_abs_max, herm_rel)
      WRITE (stdout, '(2X, A, 1X, ES12.4E3)') 'H(q) hermiticity |H-H^+|_max:', herm_abs_max
      WRITE (stdout, '(2X, A, 1X, ES12.4E3)') 'H(q) max element magnitude  :', h_abs_max
      WRITE (stdout, '(2X, A, 1X, ES12.4E3)') 'H(q) hermiticity relative   :', herm_rel
      IF (herm_rel > 1.0D-10) THEN
        WRITE (msg, '(A,1X,ES12.4E3)') 'H(q) hermiticity check failed. relative=', herm_rel
        CALL errore(1, 'build_w90_Hq', TRIM(msg))
      END IF
      IF (Hq_band) THEN
        ALLOCATE (self%eigvec(Nw, Nw, self%nkpt))
        CALL write_band(self%Hq, self%nkpt, self%eigval, self%eigvec)
      END IF
      CALL write_sep_line()
    END IF
    CALL mp_bcast(self%Hq)
  END SUBROUTINE build_w90_Hq
END MODULE wannier90
