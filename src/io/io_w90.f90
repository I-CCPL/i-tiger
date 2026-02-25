MODULE io_w90
  USE kinds, ONLY: DP
  USE mat3x3_util, ONLY: inv3x3
  IMPLICIT NONE
  PRIVATE
  PUBLIC::read_w90, w90data
  !
  TYPE::w90data_type
    REAL(DP), ALLOCATABLE::eigval(:, :)
    !< Eigenvalues (nbnd, nkpt)
    ! ========================== !
    ! ===== .chk file data ===== !
    ! ========================== !
    INTEGER::nbnd
    !< Number of bands
    INTEGER::nbnd_excl
    !< Number of excluded bands
    INTEGER, ALLOCATABLE :: excl_bands(:)
    !< Excluded bands

    REAL(DP) :: real_lattice(3, 3)
    !< Real lattice
    REAL(DP) :: recip_lattice(3, 3)
    !< Reciprocal lattice
    REAL(DP) :: recip_lattice_inv(3, 3)
    !< Inverse reciprocal lattice
    CHARACTER(LEN=33) :: chk_header
    !< Header string in .chk
    CHARACTER(LEN=20) :: chk_checkpoint
    !< Checkpoint string in .chk

    INTEGER :: nkpt
    !< Number of k-points
    INTEGER :: k_grid(3)
    !< Monkhorst-Pack grid
    REAL(DP), ALLOCATABLE :: k_cart(:, :)
    !< k-points in Cartesian coordinates (3, nkpt)
    REAL(DP), ALLOCATABLE :: k_red(:, :)
    !< k-points in reduced coordinates (3, nkpt)
    INTEGER :: nnb
    !< Number of nearest neighbor k-points
    INTEGER :: Nw
    !< Number of Wannier functions

    ! ???? :: chkpt1 !< Position of checkpoint
    LOGICAL :: have_disentangled
    !< Whether disentanglement has been performed

    REAL(DP) :: omega_invariant
    !< Omega invariant (gauge-invariant part of the spread functional)
    LOGICAL, ALLOCATABLE :: lwindow(:, :)
    !< (nbnd, nkpt)
    INTEGER, ALLOCATABLE :: ndimwin(:)
    !< (nkpt)
    COMPLEX(DP), ALLOCATABLE :: u_matrix_opt(:, :, :)
    !< (nbnd, Nw, nkpt)

    COMPLEX(DP), ALLOCATABLE :: u_matrix(:, :, :)
    !< (Nw, Nw, nkpt)
    COMPLEX(DP), ALLOCATABLE :: m_matrix(:, :, :, :)
    !< (Nw, Nw, nnb, nkpt)
    COMPLEX(DP), ALLOCATABLE :: v_matrix(:, :, :)
    !< (nbnd, Nw, nkpt)

    REAL(DP), ALLOCATABLE :: wannier_centers(:, :)
    !< (3, Nw)
    REAL(DP), ALLOCATABLE :: wannier_spreads(:)
    !< (Nw)
  END TYPE w90data_type
  TYPE(w90data_type)::w90data
CONTAINS
  SUBROUTINE read_w90()
    CALL clear_w90data()
    CALL read_w90_chk()
    CALL read_w90_eig()
    ! CALL read_w90_unk() ! TODO
    ! CALL read_w90_mmn() ! TODO
  END SUBROUTINE read_w90
  !
  SUBROUTINE read_w90_chk()
    !< Ref. wannier90/src/wannier90_readwrite.F90
    USE io_global, ONLY: ionode, stdout, prefix, get_free_unit, check_file
    USE mp_base, ONLY: mp_bcast
    INTEGER::io_unit, ios
    !
    WRITE (stdout, '(2X, A)') 'Reading .chk file...'
    CALL check_file(TRIM(prefix)//'.chk')
    IF (ionode) THEN
      io_unit = get_free_unit()
      OPEN (unit=io_unit, file=TRIM(prefix)//'.chk', form='unformatted', action='read', iostat=ios)
      CALL errore(ios, 'read_w90_chk', 'Failed to open '//TRIM(prefix)//'.chk')

      READ (io_unit) w90data%chk_header
      WRITE (stdout, '(2X, A)') '- header: '//TRIM(w90data%chk_header)
      READ (io_unit) w90data%nbnd
      WRITE (stdout, '(2X, A, I0)') '- nbnd: ', w90data%nbnd
      READ (io_unit) w90data%nbnd_excl
      WRITE (stdout, '(2X, A, I0)') '- nbnd_excl: ', w90data%nbnd_excl
      IF (w90data%nbnd_excl > 0) THEN
        ALLOCATE (w90data%excl_bands(w90data%nbnd_excl))
        READ (io_unit) w90data%excl_bands
      ELSE
        ! .chk always writes the excluded-band record; consume empty record to keep alignment.
        READ (io_unit)
      END IF
      !
      READ (io_unit) w90data%real_lattice
      WRITE (stdout, '(2X, A, 3F12.6)') '- real_lattice:', w90data%real_lattice(:, 1)
      WRITE (stdout, '(2X, A, 3F12.6)') '               ', w90data%real_lattice(:, 2)
      WRITE (stdout, '(2X, A, 3F12.6)') '               ', w90data%real_lattice(:, 3)
      READ (io_unit) w90data%recip_lattice
      w90data%recip_lattice_inv = inv3x3(w90data%recip_lattice)
      WRITE (stdout, '(2X, A, 3F12.6)') '- recip_lattice:', w90data%recip_lattice(:, 1)
      WRITE (stdout, '(2X, A, 3F12.6)') '                ', w90data%recip_lattice(:, 2)
      WRITE (stdout, '(2X, A, 3F12.6)') '                ', w90data%recip_lattice(:, 3)
      !
      READ (io_unit) w90data%nkpt
      WRITE (stdout, '(2X, A, I0)') '- nkpt: ', w90data%nkpt
      READ (io_unit) w90data%k_grid
      WRITE (stdout, '(2X, A, 3(1X,I0))') '- k_grid: ', w90data%k_grid
      ALLOCATE (w90data%k_cart(3, w90data%nkpt))
      ALLOCATE (w90data%k_red(3, w90data%nkpt))
      READ (io_unit) w90data%k_cart
      w90data%k_red = MATMUL(w90data%recip_lattice_inv, w90data%k_cart)
      READ (io_unit) w90data%nnb
      WRITE (stdout, '(2X, A, I0)') '- nnb: ', w90data%nnb
      READ (io_unit) w90data%Nw
      WRITE (stdout, '(2X, A, I0)') '- Nw: ', w90data%Nw
      ! Checkpoint string.
      READ (io_unit) w90data%chk_checkpoint
      WRITE (stdout, '(2X, A)') '- checkpoint: '//TRIM(w90data%chk_checkpoint)
      !
      READ (io_unit) w90data%have_disentangled
      WRITE (stdout, '(2X, A, L1)') '- have_disentangled: ', w90data%have_disentangled
      IF (w90data%have_disentangled) THEN
        READ (io_unit) w90data%omega_invariant
        WRITE (stdout, '(2X, A, F12.6)') '- omega_invariant: ', w90data%omega_invariant
        ALLOCATE (w90data%lwindow(w90data%nbnd, w90data%nkpt))
        ALLOCATE (w90data%ndimwin(w90data%nkpt))
        ALLOCATE (w90data%u_matrix_opt(w90data%nbnd, w90data%Nw, w90data%nkpt))
        READ (io_unit) w90data%lwindow
        READ (io_unit) w90data%ndimwin
        READ (io_unit) w90data%u_matrix_opt
      END IF
      ALLOCATE (w90data%u_matrix(w90data%Nw, w90data%Nw, w90data%nkpt))
      ALLOCATE (w90data%m_matrix(w90data%Nw, w90data%Nw, w90data%nnb, w90data%nkpt))
      READ (io_unit) w90data%u_matrix
      READ (io_unit) w90data%m_matrix
      ALLOCATE (w90data%wannier_centers(3, w90data%Nw))
      ALLOCATE (w90data%wannier_spreads(w90data%Nw))
      READ (io_unit) w90data%wannier_centers
      READ (io_unit) w90data%wannier_spreads

      CLOSE (io_unit)
    END IF

    CALL mp_bcast(w90data%nbnd)
    CALL mp_bcast(w90data%nbnd_excl)
    CALL mp_bcast(w90data%real_lattice)
    CALL mp_bcast(w90data%recip_lattice)
    CALL mp_bcast(w90data%recip_lattice_inv)
    CALL mp_bcast(w90data%nkpt)
    CALL mp_bcast(w90data%k_grid)
    CALL mp_bcast(w90data%nnb)
    CALL mp_bcast(w90data%Nw)
    CALL mp_bcast(w90data%chk_header)
    CALL mp_bcast(w90data%chk_checkpoint)
    CALL mp_bcast(w90data%have_disentangled)
    CALL mp_bcast(w90data%omega_invariant)

    IF (.NOT. ionode) THEN
      IF (w90data%nbnd_excl > 0) ALLOCATE (w90data%excl_bands(w90data%nbnd_excl))
      IF (w90data%nkpt > 0) THEN
        ALLOCATE (w90data%k_cart(3, w90data%nkpt))
        ALLOCATE (w90data%k_red(3, w90data%nkpt))
      END IF
      IF (w90data%have_disentangled) THEN
        ALLOCATE (w90data%lwindow(w90data%nbnd, w90data%nkpt))
        ALLOCATE (w90data%ndimwin(w90data%nkpt))
        ALLOCATE (w90data%u_matrix_opt(w90data%nbnd, w90data%Nw, w90data%nkpt))
      END IF
      IF (w90data%nkpt > 0 .AND. w90data%Nw > 0) THEN
        ALLOCATE (w90data%u_matrix(w90data%Nw, w90data%Nw, w90data%nkpt))
        ALLOCATE (w90data%m_matrix(w90data%Nw, w90data%Nw, w90data%nnb, w90data%nkpt))
        ALLOCATE (w90data%wannier_centers(3, w90data%Nw))
        ALLOCATE (w90data%wannier_spreads(w90data%Nw))
      END IF
    END IF

    IF (w90data%nbnd_excl > 0) CALL mp_bcast(w90data%excl_bands)
    IF (w90data%nkpt > 0) THEN
      CALL mp_bcast(w90data%k_cart)
      CALL mp_bcast(w90data%k_red)
    END IF
    IF (w90data%have_disentangled) THEN
      CALL mp_bcast(w90data%lwindow)
      CALL mp_bcast(w90data%ndimwin)
      CALL mp_bcast(w90data%u_matrix_opt)
    END IF
    IF (w90data%nkpt > 0 .AND. w90data%Nw > 0) THEN
      CALL mp_bcast(w90data%u_matrix)
      IF (w90data%nnb > 0) CALL mp_bcast(w90data%m_matrix)
      CALL mp_bcast(w90data%wannier_centers)
      CALL mp_bcast(w90data%wannier_spreads)
    END IF
  END SUBROUTINE read_w90_chk
  !
  SUBROUTINE read_w90_eig()
    !< Ref. wannier90/src/readwrite.F90
    USE io_global, ONLY: ionode, stdout, prefix, get_free_unit, check_file
    USE mp_base, ONLY: mp_bcast
    INTEGER::io_unit, ikpt, ibnd, jkpt, jbnd, ios
    CHARACTER(LEN=256)::msg
    !
    WRITE (stdout, '(2X, A)') 'Reading .eig file...'
    CALL check_file(TRIM(prefix)//'.eig')
    IF (ionode) THEN
      io_unit = get_free_unit()
      OPEN (unit=io_unit, file=TRIM(prefix)//'.eig', form='formatted', action='read', iostat=ios)
      CALL errore(ios, 'read_w90_eig', 'Failed to open '//TRIM(prefix)//'.eig')
      ALLOCATE (w90data%eigval(w90data%nbnd, w90data%nkpt))

      DO ikpt = 1, w90data%nkpt
        DO ibnd = 1, w90data%nbnd
          READ (io_unit, *, iostat=ios) jbnd, jkpt, w90data%eigval(ibnd, ikpt)
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
      ALLOCATE (w90data%eigval(w90data%nbnd, w90data%nkpt))
    END IF
    IF (w90data%nbnd > 0 .AND. w90data%nkpt > 0) CALL mp_bcast(w90data%eigval)
  END SUBROUTINE read_w90_eig

  SUBROUTINE clear_w90data()
    w90data%nbnd = 0
    w90data%nbnd_excl = 0
    w90data%real_lattice = 0.0_DP
    w90data%recip_lattice = 0.0_DP
    w90data%recip_lattice_inv = 0.0_DP
    w90data%chk_header = ''
    w90data%chk_checkpoint = ''
    w90data%nkpt = 0
    w90data%k_grid = 0
    w90data%nnb = 0
    w90data%Nw = 0
    w90data%have_disentangled = .FALSE.
    w90data%omega_invariant = 0.0_DP

    IF (ALLOCATED(w90data%eigval)) DEALLOCATE (w90data%eigval)
    IF (ALLOCATED(w90data%excl_bands)) DEALLOCATE (w90data%excl_bands)
    IF (ALLOCATED(w90data%k_cart)) DEALLOCATE (w90data%k_cart)
    IF (ALLOCATED(w90data%k_red)) DEALLOCATE (w90data%k_red)
    IF (ALLOCATED(w90data%lwindow)) DEALLOCATE (w90data%lwindow)
    IF (ALLOCATED(w90data%ndimwin)) DEALLOCATE (w90data%ndimwin)
    IF (ALLOCATED(w90data%u_matrix_opt)) DEALLOCATE (w90data%u_matrix_opt)
    IF (ALLOCATED(w90data%u_matrix)) DEALLOCATE (w90data%u_matrix)
    IF (ALLOCATED(w90data%m_matrix)) DEALLOCATE (w90data%m_matrix)
    IF (ALLOCATED(w90data%v_matrix)) DEALLOCATE (w90data%v_matrix)
    IF (ALLOCATED(w90data%wannier_centers)) DEALLOCATE (w90data%wannier_centers)
    IF (ALLOCATED(w90data%wannier_spreads)) DEALLOCATE (w90data%wannier_spreads)
  END SUBROUTINE clear_w90data

END MODULE io_w90
