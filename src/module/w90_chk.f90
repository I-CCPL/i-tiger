SUBMODULE(wannier90) w90_chk
  IMPLICIT NONE
CONTAINS
  MODULE SUBROUTINE read_w90_chk(self, chk_dum)
    !< Ref. wannier90/src/wannier90_readwrite.F90
    USE kinds, ONLY: DP
    USE io_global, ONLY: check_file
    USE mp_base, ONLY: mp_bcast
    USE system, ONLY: Nw, cell_setup, real_lattice, recip_lattice, &
                      red2cart_recip
    USE wannier90, ONLY: w90data_type, chk_dum_type
    CLASS(w90data_type), INTENT(INOUT) :: self
    TYPE(chk_dum_type), INTENT(OUT) :: chk_dum
    !
    INTEGER::io_unit, ios, ikpt, i, j, m
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

      READ (io_unit) ((real_lattice(i, j), j=1, 3), i=1, 3)
      WRITE (stdout, '(2X, A, 3F12.6)') '- real_lattice:', real_lattice(:, 1)
      WRITE (stdout, '(2X, A, 3F12.6)') '               ', real_lattice(:, 2)
      WRITE (stdout, '(2X, A, 3F12.6)') '               ', real_lattice(:, 3)
      READ (io_unit) ((recip_lattice(i, j), j=1, 3), i=1, 3)
      WRITE (stdout, '(2X, A, 3F12.6)') '- recip_lattice:', recip_lattice(:, 1)
      WRITE (stdout, '(2X, A, 3F12.6)') '                ', recip_lattice(:, 2)
      WRITE (stdout, '(2X, A, 3F12.6)') '                ', recip_lattice(:, 3)
      CALL cell_setup()

      READ (io_unit) self%kpts%nkpt
      WRITE (stdout, '(2X, A, I0)') '- nkpt: ', self%kpts%nkpt
      READ (io_unit) self%k_grid
      WRITE (stdout, '(2X, A, 3(1X,I0))') '- k_grid: ', self%k_grid
      ALLOCATE (self%kpts%k_cart(3, self%kpts%nkpt))
      ALLOCATE (self%kpts%k_red(3, self%kpts%nkpt))
      READ (io_unit) self%kpts%k_red
      CALL red2cart_recip(self%kpts%k_red, self%kpts%k_cart, self%kpts%nkpt)
      READ (io_unit) self%nnb
      WRITE (stdout, '(2X, A, I0)') '- nnb: ', self%nnb
      READ (io_unit) Nw
      WRITE (stdout, '(2X, A, I0)') '- Nw: ', Nw

      READ (io_unit) chk_dum%checkpoint
      WRITE (stdout, '(2X, A)') '- checkpoint: '//TRIM(chk_dum%checkpoint)
      READ (io_unit) chk_dum%have_disentangled
      WRITE (stdout, '(2X, A, L1)') '- have_disentangled: ', chk_dum%have_disentangled

      ALLOCATE (self%win_min(self%kpts%nkpt))
      ALLOCATE (self%ndimwin(self%kpts%nkpt))
      IF (chk_dum%have_disentangled) THEN
        READ (io_unit) chk_dum%omega_invariant
        WRITE (stdout, '(2X, A, F12.6)') '- omega_invariant: ', chk_dum%omega_invariant
        ALLOCATE (chk_dum%lwindow(self%nbnd, self%kpts%nkpt))
        ALLOCATE (chk_dum%u_matrix_opt(self%nbnd, Nw, self%kpts%nkpt))
        READ (io_unit) chk_dum%lwindow
        READ (io_unit) self%ndimwin
        READ (io_unit) chk_dum%u_matrix_opt
        DO ikpt = 1, self%kpts%nkpt
          DO j = 1, self%nbnd
            IF (chk_dum%lwindow(j, ikpt)) THEN
              self%win_min(ikpt) = j
              EXIT
            END IF
          END DO
        END DO
      ELSE
        chk_dum%omega_invariant = 0.0_DP
        ALLOCATE (chk_dum%lwindow(0, 0))
        ALLOCATE (chk_dum%u_matrix_opt(0, 0, 0))
        self%win_min = 1
        self%ndimwin = Nw
      END IF

      ALLOCATE (chk_dum%u_matrix(Nw, Nw, self%kpts%nkpt))
      ALLOCATE (chk_dum%m_matrix(Nw, Nw, self%nnb, self%kpts%nkpt))
      READ (io_unit) chk_dum%u_matrix
      READ (io_unit) chk_dum%m_matrix

      ALLOCATE (self%v_matrix(self%nbnd, Nw, self%kpts%nkpt))
      IF (chk_dum%have_disentangled) THEN
        DO ikpt = 1, self%kpts%nkpt
          DO j = 1, Nw
            ndw = self%ndimwin(ikpt)
            DO i = 1, Nw
              DO m = 1, ndw
                self%v_matrix(m, j, ikpt) = self%v_matrix(m, j, ikpt) &
                                            + chk_dum%u_matrix_opt(m, i, ikpt)*chk_dum%u_matrix(i, j, ikpt)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO ikpt = 1, self%kpts%nkpt
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
    CALL mp_bcast(self%kpts%nkpt)
    self%kpts%nktot = self%kpts%nkpt
    self%kpts%wk = 1.0_DP/REAL(self%kpts%nkpt, DP)
    CALL mp_bcast(self%k_grid)
    CALL mp_bcast(self%nnb)
    CALL mp_bcast(Nw)

    IF (.NOT. ionode) THEN
      CALL cell_setup()
      ALLOCATE (self%kpts%k_cart(3, self%kpts%nkpt))
      ALLOCATE (self%kpts%k_red(3, self%kpts%nkpt))
      ALLOCATE (self%v_matrix(self%nbnd, Nw, self%kpts%nkpt))
      ALLOCATE (self%wannier_center_cart(3, Nw))
      ALLOCATE (self%wannier_spread(Nw))
    END IF

    CALL mp_bcast(self%kpts%k_cart)
    CALL mp_bcast(self%kpts%k_red)
    CALL mp_bcast(self%v_matrix)
    CALL mp_bcast(self%wannier_center_cart)
    CALL mp_bcast(self%wannier_spread)
  END SUBROUTINE read_w90_chk
  ! ==================================================
  MODULE SUBROUTINE build_w90_Hq(self)
    !< Build Hamiltonian in q-space
    USE lin_eig_H, ONLY: eig_H
    USE mp_base, ONLY: mp_bcast
    USE system, ONLY: Nw
    USE io_output, ONLY: write_band
    CLASS(w90data_type), INTENT(INOUT) :: self
    INTEGER::ikpt, iw, jw, ibnd, ndw, mw
    COMPLEX(DP)::hval
    REAL(DP)::herm_abs_max, h_abs_max, herm_rel
    COMPLEX(DP), ALLOCATABLE::Hq(:, :, :)
    !
    WRITE (stdout, '(2X, A)') 'Building H(q) in Wannier gauge...'
    ALLOCATE (Hq(Nw, Nw, self%kpts%nkpt))
    ALLOCATE (self%Hq(Nw, Nw, self%kpts%nkpt))
    !
    IF (ionode) THEN
      DO ikpt = 1, self%kpts%nkpt
        ndw = self%ndimwin(ikpt)
        mw = self%win_min(ikpt)
        DO jw = 1, Nw
          DO iw = 1, Nw
            Hq(iw, jw, ikpt) &
              = wannier_gauge_diag(self%eigval(mw:mw + ndw - 1, ikpt), &
                                   self%v_matrix(1:ndw, iw, ikpt), self%v_matrix(1:ndw, jw, ikpt))
          END DO
        END DO
        !
        !... Force Hemiticity
        DO jw = 1, Nw
          DO iw = 1, Nw
            self%Hq(iw, jw, ikpt) = (Hq(iw, jw, ikpt) + CONJG(Hq(jw, iw, ikpt)))/2
          END DO
        END DO
      END DO
    END IF
    CALL mp_bcast(self%Hq)
    !
    ! CALL check_hermiticity(self%kpts%nkpt, 1, Hq, 1.0D-10)
    IF (ionode .AND. Hq_band) THEN
      ALLOCATE (self%eigvec(Nw, Nw, self%kpts%nkpt))
      DO ikpt = 1, self%kpts%nkpt
        CALL eig_H(Nw, Hq(:, :, ikpt), self%eigval(:, ikpt), self%eigvec(:, :, ikpt))
      END DO
      CALL write_band('itg.Hq.dat', self%eigval)
    END IF
    CALL write_sep_line()
  END SUBROUTINE build_w90_Hq
END SUBMODULE w90_chk
