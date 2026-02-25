MODULE io_w90
  USE kinds, ONLY: DP
  USE mat3x3_util, ONLY: inv3x3
  USE io_global, ONLY: ionode, stdout, get_free_unit
  IMPLICIT NONE
  PRIVATE
  PUBLIC::read_w90, w90data
  LOGICAL, PUBLIC::chk_w90 = .FALSE.
  LOGICAL, PUBLIC::Hq_band = .FALSE.
  !
  TYPE::w90data_type
    CHARACTER(LEN=256) :: prefix

    REAL(DP) :: real_lattice(3, 3)
    !< Real lattice
    REAL(DP) :: recip_lattice(3, 3)
    !< Reciprocal lattice
    REAL(DP) :: recip_lattice_inv(3, 3)
    !< Inverse reciprocal lattice

    INTEGER::nbnd
    !< Number of bands
    INTEGER :: Nw
    !< Number of Wannier functions
    INTEGER :: nkpt
    !< Number of k-points
    INTEGER :: nnb
    !< Number of nearest neighbor k-points
    INTEGER :: k_grid(3)
    !< Monkhorst-Pack grid
    REAL(DP), ALLOCATABLE :: k_cart(:, :)
    !< k-points in Cartesian coordinates (3, nkpt)
    REAL(DP), ALLOCATABLE :: k_red(:, :)
    !< k-points in reduced coordinates (3, nkpt)

    REAL(DP), ALLOCATABLE :: wannier_center_cart(:, :)
    !< (3, Nw)
    REAL(DP), ALLOCATABLE :: wannier_spread(:)
    !< (Nw)

    COMPLEX(DP), ALLOCATABLE :: v_matrix(:, :, :)
    !< (nbnd, Nw, nkpt) for disentangled case \
    !< (Nw, Nw, nkpt) for non-disentangled case.
    REAL(DP), ALLOCATABLE::eigval(:, :)
    !< Eigenvalues (nbnd, nkpt)
    COMPLEX(DP), ALLOCATABLE :: eigvec(:, :, :)
    !< Unitary matrix corresponding eigenvector \
    !< (Nw, Nw, nkpt)
    COMPLEX(DP), ALLOCATABLE :: Hq(:, :, :)
  CONTAINS
    PROCEDURE::read_eig => read_w90_eig
    PROCEDURE::read_chk => read_w90_chk
    PROCEDURE::clear_data => clear_w90_data
    PROCEDURE::build_Hq, write_band
  END TYPE w90data_type
  TYPE(w90data_type)::w90data
CONTAINS
  SUBROUTINE read_w90()
    INTEGER :: nbnd_excl
    CHARACTER(LEN=33) :: chk_header
    CHARACTER(LEN=20) :: chk_checkpoint
    LOGICAL :: have_disentangled
    REAL(DP) :: omega_invariant
    INTEGER, ALLOCATABLE :: excl_bands(:), ndimwin(:)
    LOGICAL, ALLOCATABLE :: lwindow(:, :)
    COMPLEX(DP), ALLOCATABLE :: u_matrix_opt(:, :, :), u_matrix(:, :, :), m_matrix(:, :, :, :)
    !
    CALL w90data%clear_data()
    CALL w90data%read_chk(nbnd_excl, excl_bands, chk_header, chk_checkpoint, have_disentangled, &
                          omega_invariant, lwindow, ndimwin, u_matrix_opt, u_matrix, m_matrix)
    CALL w90data%read_eig()

    IF (ionode .AND. chk_w90) THEN
      CALL write_chk_dump(w90data, nbnd_excl, excl_bands, chk_header, chk_checkpoint, have_disentangled, &
                          omega_invariant, lwindow, ndimwin, u_matrix_opt, u_matrix, m_matrix)
    END IF

    IF (ALLOCATED(excl_bands)) DEALLOCATE (excl_bands)
    IF (ALLOCATED(lwindow)) DEALLOCATE (lwindow)
    IF (ALLOCATED(ndimwin)) DEALLOCATE (ndimwin)
    IF (ALLOCATED(u_matrix_opt)) DEALLOCATE (u_matrix_opt)
    IF (ALLOCATED(u_matrix)) DEALLOCATE (u_matrix)
    IF (ALLOCATED(m_matrix)) DEALLOCATE (m_matrix)
  END SUBROUTINE read_w90
  !
  SUBROUTINE read_w90_chk(self, nbnd_excl, excl_bands, chk_header, chk_checkpoint, have_disentangled, &
                          omega_invariant, lwindow, ndimwin, u_matrix_opt, u_matrix, m_matrix)
    !< Ref. wannier90/src/wannier90_readwrite.F90
    USE io_global, ONLY: check_file
    USE mp_base, ONLY: mp_bcast
    CLASS(w90data_type), INTENT(INOUT) :: self
    INTEGER, INTENT(OUT) :: nbnd_excl
    CHARACTER(LEN=33), INTENT(OUT) :: chk_header
    CHARACTER(LEN=20), INTENT(OUT) :: chk_checkpoint
    LOGICAL, INTENT(OUT) :: have_disentangled
    REAL(DP), INTENT(OUT) :: omega_invariant
    INTEGER, ALLOCATABLE, INTENT(OUT) :: excl_bands(:), ndimwin(:)
    LOGICAL, ALLOCATABLE, INTENT(OUT) :: lwindow(:, :)
    COMPLEX(DP), ALLOCATABLE, INTENT(OUT) :: u_matrix_opt(:, :, :), u_matrix(:, :, :), m_matrix(:, :, :, :)
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

      READ (io_unit) chk_header
      WRITE (stdout, '(2X, A)') '- header: '//TRIM(chk_header)
      READ (io_unit) self%nbnd
      WRITE (stdout, '(2X, A, I0)') '- nbnd: ', self%nbnd
      READ (io_unit) nbnd_excl
      WRITE (stdout, '(2X, A, I0)') '- nbnd_excl: ', nbnd_excl
      IF (nbnd_excl > 0) THEN
        ALLOCATE (excl_bands(nbnd_excl))
        READ (io_unit) excl_bands
      ELSE
        ALLOCATE (excl_bands(0))
        ! .chk always writes the excluded-band record; consume empty record to keep alignment.
        READ (io_unit)
      END IF

      READ (io_unit) self%real_lattice
      WRITE (stdout, '(2X, A, 3F12.6)') '- real_lattice:', self%real_lattice(:, 1)
      WRITE (stdout, '(2X, A, 3F12.6)') '               ', self%real_lattice(:, 2)
      WRITE (stdout, '(2X, A, 3F12.6)') '               ', self%real_lattice(:, 3)
      READ (io_unit) self%recip_lattice
      self%recip_lattice_inv = inv3x3(self%recip_lattice)
      WRITE (stdout, '(2X, A, 3F12.6)') '- recip_lattice:', self%recip_lattice(:, 1)
      WRITE (stdout, '(2X, A, 3F12.6)') '                ', self%recip_lattice(:, 2)
      WRITE (stdout, '(2X, A, 3F12.6)') '                ', self%recip_lattice(:, 3)

      READ (io_unit) self%nkpt
      WRITE (stdout, '(2X, A, I0)') '- nkpt: ', self%nkpt
      READ (io_unit) self%k_grid
      WRITE (stdout, '(2X, A, 3(1X,I0))') '- k_grid: ', self%k_grid
      ALLOCATE (self%k_cart(3, self%nkpt))
      ALLOCATE (self%k_red(3, self%nkpt))
      READ (io_unit) self%k_cart
      self%k_red = MATMUL(self%recip_lattice_inv, self%k_cart)
      READ (io_unit) self%nnb
      WRITE (stdout, '(2X, A, I0)') '- nnb: ', self%nnb
      READ (io_unit) self%Nw
      WRITE (stdout, '(2X, A, I0)') '- Nw: ', self%Nw

      READ (io_unit) chk_checkpoint
      WRITE (stdout, '(2X, A)') '- checkpoint: '//TRIM(chk_checkpoint)
      READ (io_unit) have_disentangled
      WRITE (stdout, '(2X, A, L1)') '- have_disentangled: ', have_disentangled
      IF (have_disentangled) THEN
        READ (io_unit) omega_invariant
        WRITE (stdout, '(2X, A, F12.6)') '- omega_invariant: ', omega_invariant
        ALLOCATE (lwindow(self%nbnd, self%nkpt))
        ALLOCATE (ndimwin(self%nkpt))
        ALLOCATE (u_matrix_opt(self%nbnd, self%Nw, self%nkpt))
        READ (io_unit) lwindow
        READ (io_unit) ndimwin
        READ (io_unit) u_matrix_opt
      ELSE
        omega_invariant = 0.0_DP
        ALLOCATE (lwindow(0, 0))
        ALLOCATE (ndimwin(0))
        ALLOCATE (u_matrix_opt(0, 0, 0))
      END IF

      ALLOCATE (u_matrix(self%Nw, self%Nw, self%nkpt))
      ALLOCATE (m_matrix(self%Nw, self%Nw, self%nnb, self%nkpt))
      READ (io_unit) u_matrix
      READ (io_unit) m_matrix

      IF (have_disentangled) THEN
        ALLOCATE (self%v_matrix(self%nbnd, self%Nw, self%nkpt))
        DO ikpt = 1, self%nkpt
          ndw = ndimwin(ikpt)
          self%v_matrix(1:ndw, :, ikpt) = MATMUL(u_matrix_opt(1:ndw, :, ikpt), u_matrix(:, :, ikpt))
        END DO
      ELSE
        ! TODO: Check it is consistent with the non-disentangled case.
        IF (self%Nw /= self%nbnd) THEN
          CALL errore(1, 'read_w90_chk', 'Nw must equal nbnd for non-disentangled case')
        END IF
        !
        ALLOCATE (self%v_matrix(self%Nw, self%Nw, self%nkpt))
        DO ikpt = 1, self%nkpt
          self%v_matrix(:, :, ikpt) = u_matrix(:, :, ikpt)
        END DO
      END IF

      ALLOCATE (self%wannier_center_cart(3, self%Nw))
      ALLOCATE (self%wannier_spread(self%Nw))
      READ (io_unit) self%wannier_center_cart
      READ (io_unit) self%wannier_spread

      CLOSE (io_unit)
    END IF

    CALL mp_bcast(self%nbnd)
    CALL mp_bcast(self%real_lattice)
    CALL mp_bcast(self%recip_lattice)
    CALL mp_bcast(self%recip_lattice_inv)
    CALL mp_bcast(self%nkpt)
    CALL mp_bcast(self%k_grid)
    CALL mp_bcast(self%nnb)
    CALL mp_bcast(self%Nw)

    IF (.NOT. ionode) THEN
      ALLOCATE (self%k_cart(3, self%nkpt))
      ALLOCATE (self%k_red(3, self%nkpt))
      ALLOCATE (self%v_matrix(self%nbnd, self%Nw, self%nkpt))
      ALLOCATE (self%wannier_center_cart(3, self%Nw))
      ALLOCATE (self%wannier_spread(self%Nw))
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

  SUBROUTINE write_chk_dump(self, nbnd_excl, excl_bands, chk_header, chk_checkpoint, have_disentangled, &
                            omega_invariant, lwindow, ndimwin, u_matrix_opt, u_matrix, m_matrix)
    USE dump_vec_io, ONLY: dump_r, dump_c, dump_i, dump_l
    CLASS(w90data_type), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: nbnd_excl
    CHARACTER(LEN=*), INTENT(IN) :: chk_header, chk_checkpoint
    LOGICAL, INTENT(IN) :: have_disentangled
    REAL(DP), INTENT(IN) :: omega_invariant
    INTEGER, INTENT(IN) :: excl_bands(:), ndimwin(:)
    LOGICAL, INTENT(IN) :: lwindow(:, :)
    COMPLEX(DP), INTENT(IN) :: u_matrix_opt(:, :, :), u_matrix(:, :, :), m_matrix(:, :, :, :)
    INTEGER :: iu, ios, excl_sum
    CHARACTER(LEN=256) :: fbase
    REAL(DP), ALLOCATABLE :: rv(:)
    COMPLEX(DP), ALLOCATABLE :: cv(:)
    INTEGER, ALLOCATABLE :: iv(:)
    LOGICAL, ALLOCATABLE :: lv(:)

    fbase = TRIM(self%prefix)//'.chk_dump'
    IF (LEN_TRIM(fbase) == 0) fbase = 'itg_input_check_dump.txt'

    iu = get_free_unit()
    OPEN (unit=iu, file=TRIM(fbase), status='replace', action='write', iostat=ios)
    CALL errore(ios, 'write_chk_dump', 'Failed to open dump file: '//TRIM(fbase))

    IF (SIZE(excl_bands) > 0) THEN
      excl_sum = SUM(excl_bands)
    ELSE
      excl_sum = 0
    END IF

    WRITE (iu, '(A)') 'header='//TRIM(chk_header)
    WRITE (iu, '(A)') 'checkpoint='//TRIM(chk_checkpoint)
    WRITE (iu, '(A,I0)') 'nbnd=', self%nbnd
    WRITE (iu, '(A,I0)') 'nbnd_excl=', nbnd_excl
    WRITE (iu, '(A,I0)') 'nkpt=', self%nkpt
    WRITE (iu, '(A,I0)') 'nnb=', self%nnb
    WRITE (iu, '(A,I0)') 'Nw=', self%Nw
    WRITE (iu, '(A,3(1X,I0))') 'k_grid=', self%k_grid
    WRITE (iu, '(A,L1)') 'have_disentangled=', have_disentangled
    WRITE (iu, '(A,1X,ES24.16E3)') 'omega_invariant=', omega_invariant
    WRITE (iu, '(A,9(1X,ES24.16E3))') 'real_lattice=', self%real_lattice
    WRITE (iu, '(A,9(1X,ES24.16E3))') 'recip_lattice=', self%recip_lattice
    WRITE (iu, '(A,9(1X,ES24.16E3))') 'recip_lattice_inv=', self%recip_lattice_inv
    WRITE (iu, '(A,I0)') 'excl_bands_sum=', excl_sum
    WRITE (iu, '(A,1X,I0)') 'dims_excl_bands=', SIZE(excl_bands)
    WRITE (iu, '(A,2(1X,I0))') 'dims_k_cart=', 3, self%nkpt
    WRITE (iu, '(A,2(1X,I0))') 'dims_k_red=', 3, self%nkpt
    WRITE (iu, '(A,2(1X,I0))') 'dims_eigval=', self%nbnd, self%nkpt
    IF (have_disentangled) THEN
      WRITE (iu, '(A,2(1X,I0))') 'dims_lwindow=', SIZE(lwindow, 1), SIZE(lwindow, 2)
      WRITE (iu, '(A,1X,I0)') 'dims_ndimwin=', SIZE(ndimwin)
      WRITE (iu, '(A,3(1X,I0))') 'dims_u_matrix_opt=', SIZE(u_matrix_opt, 1), SIZE(u_matrix_opt, 2), SIZE(u_matrix_opt, 3)
    ELSE
      WRITE (iu, '(A,2(1X,I0))') 'dims_lwindow=', 0, 0
      WRITE (iu, '(A,1X,I0)') 'dims_ndimwin=', 0
      WRITE (iu, '(A,3(1X,I0))') 'dims_u_matrix_opt=', 0, 0, 0
    END IF
    WRITE (iu, '(A,3(1X,I0))') 'dims_u_matrix=', SIZE(u_matrix, 1), SIZE(u_matrix, 2), SIZE(u_matrix, 3)
    WRITE (iu, '(A,4(1X,I0))') 'dims_m_matrix=', SIZE(m_matrix, 1), SIZE(m_matrix, 2), SIZE(m_matrix, 3), SIZE(m_matrix, 4)
    WRITE (iu, '(A,2(1X,I0))') 'dims_wannier_centers=', 3, self%Nw
    WRITE (iu, '(A,1X,I0)') 'dims_wannier_spreads=', self%Nw
    CLOSE (iu)

    ALLOCATE (iv(SIZE(excl_bands)))
    iv = RESHAPE(excl_bands, [SIZE(excl_bands)])
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

    ALLOCATE (lv(SIZE(lwindow)))
    lv = RESHAPE(lwindow, [SIZE(lwindow)])
    CALL dump_l(TRIM(fbase)//'.lwindow', lv)
    DEALLOCATE (lv)

    ALLOCATE (iv(SIZE(ndimwin)))
    iv = RESHAPE(ndimwin, [SIZE(ndimwin)])
    CALL dump_i(TRIM(fbase)//'.ndimwin', iv)
    DEALLOCATE (iv)

    ALLOCATE (cv(SIZE(u_matrix_opt)))
    cv = RESHAPE(u_matrix_opt, [SIZE(u_matrix_opt)])
    CALL dump_c(TRIM(fbase)//'.u_matrix_opt', cv)
    DEALLOCATE (cv)

    ALLOCATE (cv(SIZE(u_matrix)))
    cv = RESHAPE(u_matrix, [SIZE(u_matrix)])
    CALL dump_c(TRIM(fbase)//'.u_matrix', cv)
    DEALLOCATE (cv)

    ALLOCATE (cv(SIZE(m_matrix)))
    cv = RESHAPE(m_matrix, [SIZE(m_matrix)])
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
    self%real_lattice = 0.0_DP
    self%recip_lattice = 0.0_DP
    self%recip_lattice_inv = 0.0_DP
    self%nkpt = 0
    self%k_grid = 0
    self%nnb = 0
    self%Nw = 0

    IF (ALLOCATED(self%eigval)) DEALLOCATE (self%eigval)
    IF (ALLOCATED(self%k_cart)) DEALLOCATE (self%k_cart)
    IF (ALLOCATED(self%k_red)) DEALLOCATE (self%k_red)
    IF (ALLOCATED(self%v_matrix)) DEALLOCATE (self%v_matrix)
    IF (ALLOCATED(self%wannier_center_cart)) DEALLOCATE (self%wannier_center_cart)
    IF (ALLOCATED(self%wannier_spread)) DEALLOCATE (self%wannier_spread)
    IF (ALLOCATED(self%eigvec)) DEALLOCATE (self%eigvec)
    IF (ALLOCATED(self%Hq)) DEALLOCATE (self%Hq)
  END SUBROUTINE clear_w90_data

  SUBROUTINE build_Hq(self)
    !< Build Hamiltonian in q-space
    USE io_global, ONLY: write_sep_line
    CLASS(w90data_type), INTENT(INOUT) :: self
    INTEGER::ikpt, iw, jw, ibnd
    COMPLEX(DP)::hval
    REAL(DP)::herm_abs_max, h_abs_max, herm_rel
    CHARACTER(LEN=256)::msg
    !
    IF (ionode) THEN
      WRITE (stdout, '(2X, A)') 'Building H(q) in Wannier gauge...'
      !
      ALLOCATE (self%Hq(self%Nw, self%Nw, self%nkpt))
      self%Hq = CMPLX(0.0_DP, 0.0_DP, DP)

      DO ikpt = 1, self%nkpt
        DO iw = 1, self%Nw
          DO jw = 1, self%Nw
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
      CALL check_hermiticity(self%Nw, self%nkpt, self%Hq, herm_abs_max, h_abs_max, herm_rel)
      WRITE (stdout, '(2X, A, 1X, ES12.4E3)') 'H(q) hermiticity |H-H^+|_max:', herm_abs_max
      WRITE (stdout, '(2X, A, 1X, ES12.4E3)') 'H(q) max element magnitude  :', h_abs_max
      WRITE (stdout, '(2X, A, 1X, ES12.4E3)') 'H(q) hermiticity relative   :', herm_rel
      IF (herm_rel > 1.0D-10) THEN
        WRITE (msg, '(A,1X,ES12.4E3)') 'H(q) hermiticity check failed. relative=', herm_rel
        CALL errore(1, 'build_hq', TRIM(msg))
      END IF
      IF (Hq_band) CALL self%write_band()
      CALL write_sep_line()
    END IF
  END SUBROUTINE build_Hq

  SUBROUTINE write_band(self)
    USE lin_eig_H, ONLY: eig_H
    CLASS(w90data_type), INTENT(INOUT) :: self
    INTEGER::ikpt, iw, ibnd, io_unit
    !
    WRITE (stdout, '(2X, A)') 'Diagonalizing H(q) to build bands...'

    IF (.NOT. ALLOCATED(self%Hq)) THEN
      CALL errore(1, 'write_band', 'H(q) is not built. Call build_Hq first.')
    END IF

    IF (.NOT. ALLOCATED(self%eigvec)) THEN
      ALLOCATE (self%eigvec(self%Nw, self%Nw, self%nkpt))
    END IF

    DO ikpt = 1, self%nkpt
      CALL eig_H(self%Nw, self%Hq(:, :, ikpt), &
                 self%eigval(:, ikpt), self%eigvec(:, :, ikpt))
    END DO

    io_unit = get_free_unit()
    OPEN (unit=io_unit, file=TRIM(self%prefix)//'.eigval')

    WRITE (io_unit, '("#", A)') 'ibnd, ikpt, eigval'
    DO iw = 1, self%Nw
      DO ikpt = 1, self%nkpt
        WRITE (io_unit, '(I6, I6, ES13.4E3)') iw, ikpt, self%eigval(iw, ikpt)
      END DO
    END DO
    CLOSE (io_unit)
  END SUBROUTINE write_band
END MODULE io_w90
