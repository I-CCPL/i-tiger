SUBMODULE(wannier90) w90_mmn
CONTAINS
  MODULE SUBROUTINE read_w90_mmn(self)
    USE mp_base, ONLY: mp_bcast
    USE io_global, ONLY: check_file
    CLASS(w90data_type), INTENT(INOUT) :: self
    INTEGER::io_unit, ios
    CHARACTER(LEN=512)::header
    INTEGER::nbnd, nkpt, nnb
    INTEGER::block_size, ikpt, inb, ibnd, jbnd
    INTEGER::k_from, k_to, g1, g2, g3
    REAL(DP)::re, im
    !
    WRITE (stdout, '(2X, A)') 'Reading .mmn file...'
    CALL check_file(TRIM(self%prefix)//'.mmn')
    IF (ionode) THEN
      io_unit = get_free_unit()
      OPEN (unit=io_unit, file=TRIM(self%prefix)//'.mmn', form='formatted', action='read', iostat=ios)
      CALL errore(ios, 'read_w90_mmn', 'Failed to open '//TRIM(self%prefix)//'.mmn')

      READ (io_unit, '(A)') header
      READ (io_unit, *) nbnd, nkpt, nnb
      IF (nbnd /= self%nbnd) THEN
        CALL errore(1, 'read_w90_mmn', 'Inconsistent number of bands in .mmn file')
      ELSE IF (nkpt /= self%kpts%nkpt) THEN
        CALL errore(1, 'read_w90_mmn', 'Inconsistent number of k-points in .mmn file')
      ELSE IF (nnb /= self%nnb) THEN
        CALL errore(1, 'read_w90_mmn', 'Inconsistent number of nearest neighbors in .mmn file')
      END IF
    END IF

    ALLOCATE (self%neighbour_k(self%nnb, self%kpts%nkpt))
    ALLOCATE (self%neighbour_g(3, self%nnb, self%kpts%nkpt))
    ALLOCATE (self%overlap(self%nbnd, self%nbnd, self%nnb, self%kpts%nkpt))
    IF (ionode) THEN
      DO ikpt = 1, nkpt
        DO inb = 1, nnb
          READ (io_unit, *) k_from, k_to, g1, g2, g3
          IF (k_from /= ikpt) THEN
            CALL errore(1, 'read_w90_mmn', 'Unexpected k-point index in .mmn file.')
          END IF
          self%neighbour_k(inb, ikpt) = k_to
          self%neighbour_g(:, inb, ikpt) = (/g1, g2, g3/)

          DO jbnd = 1, nbnd
            DO ibnd = 1, nbnd
              READ (io_unit, *) re, im
              self%overlap(ibnd, jbnd, inb, ikpt) = CMPLX(re, im, DP)
            END DO
          END DO
        END DO
      END DO
    END IF
    CALL mp_bcast(self%neighbour_k)
    CALL mp_bcast(self%neighbour_g)
    CALL mp_bcast(self%overlap)
  END SUBROUTINE read_w90_mmn

  MODULE SUBROUTINE build_w90_bvec(self)
    !< Build b vectors, the index map from (nnb, nkpt) to (nnb) and weight factor w_b
    USE kinds, ONLY: eq_real, eq_vec_real
    USE algo_unique, ONLY: qsort_perm
    USE system, ONLY: red2cart_recip
    CLASS(w90data_type), INTENT(inout) :: self
    INTEGER::ikpt, inb, jnb, iknb
    REAL(DP), ALLOCATABLE::bvec_red(:, :)
    REAL(DP)::bvec_cart(3)
    REAL(DP), ALLOCATABLE::bvec_length(:)
    INTEGER, ALLOCATABLE::bvec_sort_index(:)
    INTEGER::ishell, nshell, inb_shell
    INTEGER, ALLOCATABLE::nnb_shell(:)
    !> matrices for w_b calculation (Aw=q, A=USV^T)
    REAL(DP), ALLOCATABLE::A(:, :), U(:, :), S(:), VT(:, :), work(:)
    REAL(DP), ALLOCATABLE::w_shell(:)
    REAL(DP)::I(9)
    INTEGER::ldim, lwork, info, idx, jdx, kdx
    CHARACTER(LEN=256) :: msg

    !... Build bvec_red
    WRITE (stdout, '(2X, A)') 'Building b vectors...'
    ALLOCATE (bvec_red(3, self%nnb))
    ALLOCATE (bvec_length(self%nnb))
    ALLOCATE (bvec_sort_index(self%nnb))
    ! List of bvec is independent of k-point, so we can just use the first k-point
    ikpt = 1
    DO inb = 1, self%nnb
      iknb = self%neighbour_k(inb, ikpt)
      bvec_red(:, inb) = REAL(self%neighbour_g(:, inb, ikpt), DP) &
                         + self%kpts%k_red(:, iknb) - self%kpts%k_red(:, ikpt)
      CALL red2cart_recip(bvec_red(:, inb), bvec_cart)
      bvec_length(inb) = NORM2(bvec_cart)
    END DO
    ! Sort by vector length
    CALL qsort_perm(bvec_length, bvec_sort_index)
    ALLOCATE (self%bvec_red(3, self%nnb))
    DO inb = 1, self%nnb
      self%bvec_red(:, inb) = bvec_red(:, bvec_sort_index(inb))
    END DO

    !... Build the index map from (nnb, nkpt) to (nnb)
    WRITE (stdout, '(2X,A)') '- Building index map...'
    ALLOCATE (self%bvec_index(self%nnb, self%kpts%nkpt))
    self%bvec_index = 0
    DO ikpt = 1, self%kpts%nkpt
      DO inb = 1, self%nnb
        iknb = self%neighbour_k(inb, ikpt)
        bvec_red(:, 1) = REAL(self%neighbour_g(:, inb, ikpt), DP) &
                         + self%kpts%k_red(:, iknb) - self%kpts%k_red(:, ikpt)
        DO jnb = 1, self%nnb
          IF (eq_vec_real(bvec_red(:, 1), self%bvec_red(:, jnb), 1.0D-6)) THEN
            self%bvec_index(inb, ikpt) = jnb
            EXIT
          END IF
        END DO
      END DO
    END DO

    ! Check if all bvecs are mapped
    ! Try to reduce the tolerance if not mapped.
    DO ikpt = 1, self%kpts%nkpt
      DO inb = 1, self%nnb
        IF (self%bvec_index(inb, ikpt) == 0) THEN
          WRITE (msg, '(A, I4, A, I4)') 'Failed to find bvec index for inb=', inb, ', ikpt=', ikpt
          CALL errore(1, 'build_w90_bvec', TRIM(msg))
        END IF
      END DO
    END DO

    !... Group bvecs into shells
    ALLOCATE (nnb_shell(self%nnb))
    nshell = 1
    inb_shell = 1
    DO inb = 2, self%nnb
      IF (eq_real(bvec_length(bvec_sort_index(inb - 1)), &
                  bvec_length(bvec_sort_index(inb)), 1.0D-6)) THEN
        inb_shell = inb_shell + 1
      ELSE
        nnb_shell(nshell) = inb_shell
        inb_shell = 1
        nshell = nshell + 1
      END IF
    END DO
    nnb_shell(nshell) = inb_shell
    WRITE (stdout, '(2X,A, I0)') '- Number of bvec shells: ', nshell
    DEALLOCATE (bvec_length, bvec_sort_index)

    !... Build weight factor w_b
    ldim = MIN(nshell, 9)
    ALLOCATE (A(nshell, 9))
    ALLOCATE (U(ldim, nshell))
    ALLOCATE (S(ldim))
    ALLOCATE (VT(ldim, 9))
    ALLOCATE (work(1))

    inb = 1
    A(:, :) = 0.0_DP
    DO ishell = 1, nshell
      inb_shell = nnb_shell(ishell)
      DO jnb = inb, inb + inb_shell - 1
        CALL red2cart_recip(self%bvec_red(:, jnb), bvec_cart)
        kdx = 0
        DO idx = 1, 3
          DO jdx = 1, 3
            kdx = kdx + 1
            A(ishell, kdx) = A(ishell, kdx) + bvec_cart(idx)*bvec_cart(jdx)
          END DO
        END DO
      END DO
      inb = inb + inb_shell
    END DO

    CALL DGESVD('S', 'S', nshell, 9, A, nshell, S, U, ldim, VT, ldim, work, -1, info)
    IF (info /= 0) CALL errore(info, 'build_w90_bvec: dgesvd query failed')
    lwork = INT(work(1))
    DEALLOCATE (work)
    ALLOCATE (work(lwork))
    CALL DGESVD('S', 'S', nshell, 9, A, nshell, S, U, ldim, VT, ldim, work, lwork, info)
    IF (info /= 0) CALL errore(info, 'build_w90_bvec: dgesvd failed')

    ALLOCATE (w_shell(nshell))
    I = 0.0_DP
    I(1::4) = 1.0_DP
    w_shell(:) = 0.0_DP
    DO idx = 1, 9
      DO jdx = 1, ldim
        w_shell(:) = w_shell(:) + I(idx)*U(jdx, :)/S(jdx)*VT(jdx, idx)
      END DO
    END DO
    WRITE (stdout, '(2X, A)') '- Weight factors for bvec shells:'
    DO ishell = 1, nshell
      WRITE (stdout, '(4X, "=", I2,". ", F9.4)') ishell, w_shell(ishell)
    END DO

    ALLOCATE (self%wb(self%nnb))
    inb = 1
    DO ishell = 1, nshell
      inb_shell = nnb_shell(ishell)
      self%wb(inb:inb + inb_shell - 1) = w_shell(ishell)
      inb = inb + inb_shell
    END DO
    DEALLOCATE (A, U, S, VT, work, w_shell)
    CALL write_sep_line()
  END SUBROUTINE build_w90_bvec

  MODULE SUBROUTINE build_w90_Aq(self)
    !< Build A(q) in Wannier gauge
    USE constants, ONLY: zi
    USE mp_base, ONLY: mp_bcast
    USE system, ONLY: Nw, red2cart_recip
    CLASS(w90data_type), INTENT(INOUT) :: self
    INTEGER::ikpt, inb, jnb, iknb, ibnd, jbnd, iw, jw, ipol
    INTEGER::ndw1, ndw2, mw1, mw2
    REAL(DP)::b_cart(3)
    COMPLEX(DP)::M_W, A_qb(3, Nw, Nw)
    !< overlap matrix in Wannier gauge
    !
    WRITE (stdout, '(2X, A)') 'Building A(q)...'
    ALLOCATE (self%Aq(3, Nw, Nw, self%kpts%nkpt))
    IF (ionode) THEN
      DO ikpt = 1, self%kpts%nkpt
        ndw1 = self%ndimwin(ikpt)
        mw1 = self%win_min(ikpt)
        A_qb(:, :, :) = 0.0_DP
        DO inb = 1, self%nnb
          iknb = self%neighbour_k(inb, ikpt)
          jnb = self%bvec_index(inb, ikpt)
          ndw2 = self%ndimwin(iknb)
          mw2 = self%win_min(iknb)
          CALL red2cart_recip(self%bvec_red(:, jnb), b_cart)
          DO jw = 1, Nw
            DO iw = 1, Nw
              M_W = wannier_gauge(self%overlap(mw1:mw1 + ndw1 - 1, mw2:mw2 + ndw2 - 1, inb, ikpt), &
                                  self%v_matrix(1:ndw1, iw, ikpt), self%v_matrix(1:ndw2, jw, iknb))
              A_qb(:, iw, jw) = A_qb(:, iw, jw) + zi*self%wb(jnb)*M_W*b_cart(:)
            END DO
          END DO
        END DO

        !... Enforce Hermicity
        DO ipol = 1, 3
          self%Aq(ipol, :, :, ikpt) = 0.5_DP*(A_qb(ipol, :, :) + CONJG(TRANSPOSE(A_qb(ipol, :, :))))
        END DO
      END DO
    END IF

    CALL mp_bcast(self%Aq)
    CALL write_sep_line()
  END SUBROUTINE build_w90_Aq

END SUBMODULE w90_mmn
