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

  MODULE SUBROUTINE build_w90_Aq(self)
    !< Build A(q) in Wannier gauge
    USE constants, ONLY: zi
    USE system, ONLY: Nw, red2cart_recip
    CLASS(w90data_type), INTENT(INOUT) :: self
    INTEGER::ikpt, inb, iknb, ibnd, jbnd, iw, jw
    REAL(DP)::b_red(3), b_cart(3)
    COMPLEX(DP)::M_W, A_qb(3)
    !< overlap matrix in Wannier gauge
    !
    WRITE (stdout, '(2X, A)') 'Building A(q)...'
    CALL errore(1, 'build_w90_Aq', 'Not implemented yet')
    ALLOCATE (self%Aq(3, Nw, Nw, self%kpts%nkpt))
    self%Aq = CMPLX(0.0_DP, 0.0_DP, DP)
    DO inb = 1, self%nnb
      DO ikpt = 1, self%kpts%nkpt
        iknb = self%neighbour_k(inb, ikpt)
        b_red(:) = REAL(self%neighbour_g(:, inb, ikpt), DP) + self%kpts%k_red(:, iknb) - self%kpts%k_red(:, ikpt)
        CALL red2cart_recip(b_red, b_cart)
        DO iw = 1, Nw
          DO jw = 1, Nw
            M_W = wannier_gauge(self%nbnd, self%overlap(:, :, inb, ikpt), &
                                self%v_matrix(:, iw, ikpt), self%v_matrix(:, jw, iknb))

            A_qb(:) = zi*self%wb(inb)*M_W*b_cart(:)
            self%Aq(:, iw, jw, ikpt) = self%Aq(:, iw, jw, ikpt) + A_qb(:)
            WRITE (300, '(4I3,3(2X,SP,E11.4,E11.4,"j"))') inb, ikpt, iw, jw, A_qb(1), A_qb(2), A_qb(3)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE build_w90_Aq

END SUBMODULE w90_mmn
