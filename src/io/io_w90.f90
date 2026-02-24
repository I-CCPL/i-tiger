MODULE io_w90
  USE kinds, ONLY: DP
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
    COMPLEX(DP), ALLOCATABLE :: m_matrix(:, :, :)
    !< (Nw, Nw, nkpt)
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
    CALL read_w90_chk()
    CALL read_w90_eig()
    ! CALL read_w90_unk() ! TODO
    ! CALL read_w90_mmn() ! TODO
  END SUBROUTINE read_w90
  !
  SUBROUTINE read_w90_chk()
    !< Ref. wannier90/src/wannier90_readwrite.F90
    USE io_global, ONLY: ionode, stdout, prefix, get_free_unit, check_file
    INTEGER::io_unit
    CHARACTER(LEN=33)::header
    CHARACTER(LEN=20)::checkpoint
    !
    WRITE (stdout, '(2X, A)') 'Reading .chk file...'
    CALL check_file(TRIM(prefix)//'.chk')
    IF (ionode) THEN
      io_unit = get_free_unit()
      OPEN (unit=io_unit, file=TRIM(prefix)//'.chk', form='unformatted', action='read')

      READ (io_unit) header
      WRITE (stdout, '(2X, A)') '- header: '//TRIM(header)
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
      WRITE (stdout, '(2X, A, 3F12.6)') '- recip_lattice:', w90data%recip_lattice(:, 1)
      WRITE (stdout, '(2X, A, 3F12.6)') '                ', w90data%recip_lattice(:, 2)
      WRITE (stdout, '(2X, A, 3F12.6)') '                ', w90data%recip_lattice(:, 3)
      !
      READ (io_unit) w90data%nkpt
      WRITE (stdout, '(2X, A, I0)') '- nkpt: ', w90data%nkpt
      READ (io_unit) w90data%k_grid
      WRITE (stdout, '(2X, A, 3(1X,I0))') '- k_grid: ', w90data%k_grid
      ALLOCATE (w90data%k_cart(3, w90data%nkpt))
      ! ALLOCATE (w90data%k_red(3, w90data%nkpt))
      READ (io_unit) w90data%k_cart
      READ (io_unit) w90data%nnb
      WRITE (stdout, '(2X, A, I0)') '- nnb: ', w90data%nnb
      READ (io_unit) w90data%Nw
      WRITE (stdout, '(2X, A, I0)') '- Nw: ', w90data%Nw
      ! Checkpoint string.
      READ (io_unit) checkpoint
      WRITE (stdout, '(2X, A)') '- checkpoint: '//TRIM(checkpoint)
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
      ALLOCATE (w90data%m_matrix(w90data%Nw, w90data%Nw, w90data%nkpt))
      READ (io_unit) w90data%u_matrix
      READ (io_unit) w90data%m_matrix
      ALLOCATE (w90data%wannier_centers(3, w90data%Nw))
      ALLOCATE (w90data%wannier_spreads(w90data%Nw))
      READ (io_unit) w90data%wannier_centers
      READ (io_unit) w90data%wannier_spreads

      CLOSE (io_unit)
    END IF
    !
    ! TODO: Broadcast w90data
  END SUBROUTINE read_w90_chk
  !
  SUBROUTINE read_w90_eig()
    !< Ref. wannier90/src/readwrite.F90
    USE io_global, ONLY: ionode, stdout, prefix, get_free_unit, check_file
    INTEGER::io_unit, ikpt, ibnd, jkpt, jbnd
    !
    WRITE (stdout, '(2X, A)') 'Reading .eig file...'
    CALL check_file(TRIM(prefix)//'.eig')
    IF (ionode) THEN
      io_unit = get_free_unit()
      OPEN (unit=io_unit, file=TRIM(prefix)//'.eig', form='formatted', action='read')
      ALLOCATE (w90data%eigval(w90data%nbnd, w90data%nkpt))

      DO ikpt = 1, w90data%nkpt
        DO ibnd = 1, w90data%nbnd
          READ (io_unit, '(2i, f)') jbnd, jkpt, w90data%eigval(ibnd, ikpt)
        END DO
      END DO
      CLOSE (io_unit)
    END IF
    !
    ! TODO: Broadcast w90data%eigval
  END SUBROUTINE read_w90_eig
END MODULE io_w90
