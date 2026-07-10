MODULE R_vector
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  ! TODO: cell_expand from input
  INTEGER::cell_expand(3) = (/2, 2, 2/)
  !
  TYPE::R_vec_type
    INTEGER::nu_shift
    !< Number of unique Wannier center shifts
    REAL(DP), ALLOCATABLE :: shift_cart(:, :, :)
    REAL(DP), ALLOCATABLE :: shift_cart_u(:, :)
    !< Unique list of Wannier center shifts in Cartesian coordinates (3, nu_shift)
    INTEGER, ALLOCATABLE:: shift_map_inv(:, :)
    !< Map from Wannier function pair (iw, jw) to shift index

    INTEGER::nR0pt
    !< Assume nR0pt == nkpt
    INTEGER::R0_grid(3)
    REAL(DP), ALLOCATABLE::R0_red(:, :)
    !< R0_vec in reduced coordinates (3, nR0pt)
    REAL(DP), ALLOCATABLE::R0_cart(:, :)
    !< R0_vec in Cartesian coordinates (3, nR0pt)

    INTEGER::nRpt
    !< Number of R0+T
    REAL(DP), ALLOCATABLE::R_red(:, :)
    !< R_vec in reduced coordinates (3, nRpt)
    REAL(DP), ALLOCATABLE::R_cart(:, :)
    !< R_vec in Cartesian coordinates (3, nRpt)
    REAL(DP), ALLOCATABLE::w_R(:, :, :)
    !< weight for each R vector (Nw, Nw, nRpt)

    INTEGER::max_degen
    INTEGER, ALLOCATABLE::map_r02r(:, :, :)
    !< map from r0 to corresponding arbitrary r (nR0pt, iuw, idegen)
    INTEGER, ALLOCATABLE::map_r2r0(:, :, :)
    !< map from r to corresponding arbitrary r0 (2, iuw, nRpt)
  CONTAINS
    PROCEDURE::clear => clear_Rvec
    PROCEDURE::build_shift => build_shift_vecs
    PROCEDURE::bcast_shift => bcast_shift_vecs
    PROCEDURE::build_R => build_Rvecs
    PROCEDURE::bcast_R => bcast_Rvecs
    PROCEDURE::build_ws => build_ws
    PROCEDURE::bcast_ws => bcast_ws
  END TYPE R_vec_type
  PUBLIC::R_vec_type
  !
CONTAINS
  SUBROUTINE clear_Rvec(self)
    CLASS(R_vec_type), INTENT(INOUT)::self
    IF (ALLOCATED(self%shift_cart)) DEALLOCATE (self%shift_cart)
    IF (ALLOCATED(self%shift_cart_u)) DEALLOCATE (self%shift_cart_u)
    IF (ALLOCATED(self%shift_map_inv)) DEALLOCATE (self%shift_map_inv)
    IF (ALLOCATED(self%R0_red)) DEALLOCATE (self%R0_red)
    IF (ALLOCATED(self%R0_cart)) DEALLOCATE (self%R0_cart)
    IF (ALLOCATED(self%R_red)) DEALLOCATE (self%R_red)
    IF (ALLOCATED(self%R_cart)) DEALLOCATE (self%R_cart)
    IF (ALLOCATED(self%w_R)) DEALLOCATE (self%w_R)
    IF (ALLOCATED(self%map_r02r)) DEALLOCATE (self%map_r02r)
    IF (ALLOCATED(self%map_r2r0)) DEALLOCATE (self%map_r2r0)
  END SUBROUTINE clear_Rvec
  !
  SUBROUTINE build_shift_vecs(self, wannier_center_cart)
    !< Build Wannier center shift (r_m-r_n) in cartesian coordinates
    USE io_global, ONLY: stdout
    USE system, ONLY: Nw
    USE algo_unique, ONLY: unique_vec3_inv
    CLASS(R_vec_type), INTENT(INOUT) :: self
    REAL(DP), INTENT(IN) :: wannier_center_cart(3, Nw)
    INTEGER::shift_map_inv(Nw*Nw)
    INTEGER::iw, jw
    !
    IF (ALLOCATED(self%shift_cart_u)) RETURN
    WRITE (stdout, '(2X, A)') '- Building Wannier center shift vectors...'
    !
    ALLOCATE (self%shift_cart(3, Nw, Nw))
    DO jw = 1, Nw
      DO iw = 1, Nw
        self%shift_cart(:, iw, jw) = -wannier_center_cart(:, iw) + wannier_center_cart(:, jw)
      END DO
    END DO

    ! reduce shifts to unique ones
    CALL unique_vec3_inv(self%shift_cart(1, 1, 1), Nw*Nw, 1D-8, self%shift_cart_u, self%nu_shift, shift_map_inv)
    WRITE (stdout, '(2X, A, I0)') '- Number of unique shifts: ', self%nu_shift

    ALLOCATE (self%shift_map_inv(Nw, Nw))
    DO jw = 1, Nw
      DO iw = 1, Nw
        self%shift_map_inv(iw, jw) = shift_map_inv(iw + (jw - 1)*Nw)
      END DO
    END DO
  END SUBROUTINE build_shift_vecs
  SUBROUTINE bcast_shift_vecs(self)
    USE mp_base, ONLY: mp_bcast
    USE io_global, ONLY: ionode
    USE system, ONLY: Nw
    CLASS(R_vec_type), INTENT(inout) :: self
    CALL mp_bcast(self%nu_shift)
    IF (.NOT. ionode) THEN
      ALLOCATE (self%shift_cart(3, Nw, Nw))
      ALLOCATE (self%shift_cart_u(3, self%nu_shift))
      ALLOCATE (self%shift_map_inv(Nw, Nw))
    END IF
    CALL mp_bcast(self%shift_cart)
    CALL mp_bcast(self%shift_cart_u)
    CALL mp_bcast(self%shift_map_inv)
  END SUBROUTINE bcast_shift_vecs

  SUBROUTINE build_Rvecs(self, w90data)
    USE kinds, ONLY: eq_real
    USE constants, ONLY: vec_0
    USE io_global, ONLY: stdout
    USE wannier90, ONLY: w90data_type
    USE system, ONLY: Nw, red2cart_real, cart2red_real
    CLASS(R_vec_type), INTENT(INOUT)::self
    TYPE(w90data_type), INTENT(INOUT)::w90data
    INTEGER::iuw, iw, jw
    INTEGER::ir0pt, irpt, jrpt, nRpt, r0pt_x, r0pt_y, r0pt_z
    INTEGER::icell, ncell, cell_x, cell_y, cell_z, ndegen
    REAL(DP)::R_vector(3), R_dist
    REAL(DP), ALLOCATABLE::Tvec_red(:, :), Tvec_cart(:, :)
    REAL(DP), ALLOCATABLE::Rvec_cart(:, :), dist_min(:)
    REAL(DP), ALLOCATABLE::w_R(:, :, :)
    INTEGER::tmp_arr(3)
    INTEGER::cell_range(3)
    REAL(DP)::dr_nmR0(3) !< R0 + r_m - r_n
    LOGICAL, ALLOCATABLE::bRvec_selected(:)
    !< whether Rvec is selected for given R0
    CALL errore(1, 'build_Rvecs', 'Deprecated subroutine. Use build_ws instead.')

    self%nR0pt = w90data%kpts%nkpt

    CALL self%build_shift(w90data%wannier_center_cart)

    WRITE (stdout, '(2X, A)') '- Building R vectors for Fourier transform...'
    ! Build T vectors
    cell_range(:) = 2*cell_expand(:) + 1
    ncell = PRODUCT(cell_range)
    WRITE (stdout, '(2X, A)') '- Building T vectors for periodic images...'
    WRITE (stdout, '(2X, A, I0)') '- Number of considered T vectors: ', ncell
    ALLOCATE (Tvec_red(3, ncell))
    ALLOCATE (Tvec_cart(3, ncell))
    tmp_arr(:) = -cell_expand(:)
    DO icell = 1, ncell
      CALL grid_idx2xyz(cell_range, icell, cell_x, cell_y, cell_z, tmp_arr)
      Tvec_red(:, icell) = (/cell_x*w90data%k_grid(1), &
                             cell_y*w90data%k_grid(2), &
                             cell_z*w90data%k_grid(3)/)
      CALL red2cart_real(Tvec_red(:, icell), Tvec_cart(:, icell))
    END DO

    ! Find T such that minimize |R0+T+r_n-r_m| for given (R0, n, m)
    WRITE (stdout, '(2X, A)') '- Building R vectors for each shift and R0...'
    ALLOCATE (self%R0_red(3, self%nR0pt))
    ALLOCATE (self%R0_cart(3, self%nR0pt))
    ALLOCATE (dist_min(self%nu_shift))
    nRpt = self%nR0pt*ncell
    ALLOCATE (Rvec_cart(3, nRpt))
    ALLOCATE (w_R(Nw, Nw, nRpt))
    ALLOCATE (bRvec_selected(nRpt))
    self%R0_grid(:) = w90data%k_grid(:)
    bRvec_selected = .FALSE.

    DO ir0pt = 1, self%nR0pt
      ! Build R0 vector
      CALL grid_idx2xyz(self%R0_grid, ir0pt, r0pt_x, r0pt_y, r0pt_z, vec_0)
      self%R0_red(:, ir0pt) = (/REAL(r0pt_x, DP), REAL(r0pt_y, DP), REAL(r0pt_z, DP)/)
      CALL red2cart_real(self%R0_red(:, ir0pt), self%R0_cart(:, ir0pt))
      ! Build R vector
      DO icell = 1, ncell
        irpt = ir0pt + (icell - 1)*self%nR0pt
        Rvec_cart(:, irpt) = self%R0_cart(:, ir0pt) + Tvec_cart(:, icell)
      END DO
      DO iuw = 1, self%nu_shift
        ! minimize |R0+T+r_n-r_m| by searching T vectors
        dist_min(iuw) = 1.0D10
        DO icell = 1, ncell
          irpt = ir0pt + (icell - 1)*self%nR0pt
          R_vector = Rvec_cart(:, irpt) + self%shift_cart_u(:, iuw)
          R_dist = SUM(R_vector**2)
          IF (R_dist < dist_min(iuw)) THEN
            dist_min(iuw) = R_dist
          END IF
        END DO
      END DO

      DO jw = 1, Nw
        DO iw = 1, Nw
          ! R vectors for each T vector with distance close to dist_min
          ndegen = 0
          iuw = self%shift_map_inv(iw, jw)
          DO icell = 1, ncell
            irpt = ir0pt + (icell - 1)*self%nR0pt
            R_vector = Rvec_cart(:, irpt) + self%shift_cart_u(:, iuw)
            R_dist = SUM(R_vector**2)
            IF (eq_real(R_dist, dist_min(iuw), 1D-12)) THEN
              ndegen = ndegen + 1
              w_R(iw, jw, irpt) = 1.0D0
              bRvec_selected(irpt) = .TRUE.
            ELSE
              w_R(iw, jw, irpt) = 0.0D0
            END IF
          END DO
          ! weight for degenerate T vectors
          IF (ndegen > 0) THEN
            DO icell = 1, ncell
              irpt = ir0pt + (icell - 1)*self%nR0pt
              w_R(iw, jw, irpt) = w_R(iw, jw, irpt)/REAL(ndegen, DP)
            END DO
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE (Tvec_red)
    DEALLOCATE (Tvec_cart)

    WRITE (stdout, '(2X, A)') '- Selecting R vectors...'
    WRITE (stdout, '(2X, A, I0)') '- Number of R=R0+T vectors: ', nRpt
    ! unique R vectors are selected
    self%nRpt = COUNT(bRvec_selected)
    WRITE (stdout, '(2X, A, I0)') '- Number of selected R vectors: ', self%nRpt
    ALLOCATE (self%R_red(3, self%nRpt))
    ALLOCATE (self%R_cart(3, self%nRpt))
    ALLOCATE (self%w_R(Nw, Nw, self%nRpt))
    jrpt = 0
    DO irpt = 1, nRpt
      IF (bRvec_selected(irpt)) THEN
        jrpt = jrpt + 1
        self%R_cart(:, jrpt) = Rvec_cart(:, irpt)
        CALL cart2red_real(self%R_cart(:, jrpt), self%R_red(:, jrpt))
        self%w_R(:, :, jrpt) = w_R(:, :, irpt)
      END IF
    END DO
    IF (ALLOCATED(Rvec_cart)) DEALLOCATE (Rvec_cart)
    IF (ALLOCATED(w_R)) DEALLOCATE (w_R)
    IF (ALLOCATED(bRvec_selected)) DEALLOCATE (bRvec_selected)
  END SUBROUTINE build_Rvecs

  SUBROUTINE bcast_Rvecs(self)
    USE mp_base, ONLY: mp_bcast
    USE io_global, ONLY: ionode
    USE system, ONLY: Nw
    CLASS(R_vec_type), INTENT(inout) :: self
    CALL self%bcast_shift()
    !
    CALL mp_bcast(self%nR0pt)
    CALL mp_bcast(self%nRpt)
    IF (.NOT. ionode) THEN
      ALLOCATE (self%R0_red(3, self%nR0pt))
      ALLOCATE (self%R0_cart(3, self%nR0pt))
      ALLOCATE (self%R_red(3, self%nRpt))
      ALLOCATE (self%R_cart(3, self%nRpt))
      ALLOCATE (self%w_R(Nw, Nw, self%nRpt))
    END IF
    CALL mp_bcast(self%R0_grid)
    CALL mp_bcast(self%R0_red)
    CALL mp_bcast(self%R0_cart)
    CALL mp_bcast(self%R_red)
    CALL mp_bcast(self%R_cart)
    CALL mp_bcast(self%w_R)
  END SUBROUTINE bcast_Rvecs

  SUBROUTINE build_ws(self, w90data)
    USE io_global, ONLY: stdout
    USE system, ONLY: Nw, red2cart_real
    USE wannier90, ONLY: w90data_type
    CLASS(R_vec_type), INTENT(INOUT)::self
    TYPE(w90data_type), INTENT(INOUT)::w90data
    INTEGER::n1, n2, n3, icnt, i1, i2, i3, degen, i
    INTEGER::iw, jw, iuw, ir0pt
    REAL(DP)::ndiff_red(3), ndiff_cart(3)
    REAL(DP)::dist_min
    REAL(DP), ALLOCATABLE::dist(:)

    REAL(DP), ALLOCATABLE::w_R0(:), minRT(:, :, :), R_red(:, :)
    INTEGER::nRpt, ncell, irpt, idx
    INTEGER::cell_range(3)
    INTEGER, ALLOCATABLE::ndegen(:, :), degen_idx(:, :)
    LOGICAL, ALLOCATABLE::bRvec_selected(:)
    cell_range(:) = 2*cell_expand(:) + 1
    ncell = PRODUCT(cell_range)
    ALLOCATE (dist(ncell))

    CALL self%build_shift(w90data%wannier_center_cart)

    !... find T vectors to find equivalent R vectors for each R0
    self%nR0pt = 0
    DO n1 = -w90data%k_grid(1), w90data%k_grid(1)
      DO n2 = -w90data%k_grid(2), w90data%k_grid(2)
        DO n3 = -w90data%k_grid(3), w90data%k_grid(3)
          ! Loop over the 125 points R. R=0 corresponds to i1=i2=i3=0,
          ! or icnt=63
          icnt = 0
          DO i1 = -cell_expand(1), cell_expand(1)
            DO i2 = -cell_expand(2), cell_expand(2)
              DO i3 = -cell_expand(3), cell_expand(3)
                icnt = icnt + 1
                ! Calculate distance squared |r-R|^2
                ndiff_red(1) = n1 - i1*w90data%k_grid(1)
                ndiff_red(2) = n2 - i2*w90data%k_grid(2)
                ndiff_red(3) = n3 - i3*w90data%k_grid(3)
                CALL red2cart_real(ndiff_red, ndiff_cart)
                dist(icnt) = SUM(ndiff_cart**2)
              END DO
            END DO
          END DO
          dist_min = MINVAL(dist)
          IF (ABS(dist(ncell/2 + 1) - dist_min) < 1.E-7_DP) THEN
            self%nR0pt = self%nR0pt + 1
          END IF
        END DO
      END DO
    END DO

    WRITE (stdout, '(2X, A, I0)') '- Searching R0 vectors: ', self%nR0pt
    ALLOCATE (self%R0_red(3, self%nR0pt))
    ALLOCATE (self%R0_cart(3, self%nR0pt))
    ALLOCATE (w_R0(self%nR0pt))
    self%nR0pt = 0
    DO n1 = -w90data%k_grid(1), w90data%k_grid(1)
      DO n2 = -w90data%k_grid(2), w90data%k_grid(2)
        DO n3 = -w90data%k_grid(3), w90data%k_grid(3)
          ! Loop over the 125 points R. R=0 corresponds to i1=i2=i3=0,
          ! or icnt=63
          icnt = 0
          DO i1 = -cell_expand(1), cell_expand(1)
            DO i2 = -cell_expand(2), cell_expand(2)
              DO i3 = -cell_expand(3), cell_expand(3)
                icnt = icnt + 1
                ! Calculate distance squared |r-R|^2
                ndiff_red(1) = n1 - i1*w90data%k_grid(1)
                ndiff_red(2) = n2 - i2*w90data%k_grid(2)
                ndiff_red(3) = n3 - i3*w90data%k_grid(3)
                CALL red2cart_real(ndiff_red, ndiff_cart)
                dist(icnt) = SUM(ndiff_cart**2)
              END DO
            END DO
          END DO
          dist_min = MINVAL(dist)
          IF (ABS(dist(ncell/2 + 1) - dist_min) < 1.E-7_DP) THEN
            self%nR0pt = self%nR0pt + 1
            degen = 0
            DO i = 1, ncell
              IF (ABS(dist(i) - dist_min) < 1.E-7_DP) &
                degen = degen + 1
            END DO
            w_R0(self%nR0pt) = 1.0D0/REAL(degen, DP)
            self%R0_red(1, self%nR0pt) = n1
            self%R0_red(2, self%nR0pt) = n2
            self%R0_red(3, self%nR0pt) = n3
            CALL red2cart_real(self%R0_red(:, self%nR0pt), self%R0_cart(:, self%nR0pt))
          END IF
        END DO
      END DO
    END DO

    !... Search T vectors to minimize |R0+T+r_n-r_m| for each R0 and shift
    nRpt = self%nR0pt*ncell
    WRITE (stdout, '(2X, A, I0)') '- Searching R vectors: ', nRpt
    ALLOCATE (bRvec_selected(nRpt))
    ALLOCATE (R_red(3, nRpt))
    ALLOCATE (ndegen(self%nu_shift, self%nR0pt))
    ALLOCATE (degen_idx(self%nu_shift, nRpt))
    ALLOCATE (minRT(3, self%nu_shift, self%nR0pt))

    bRvec_selected = .FALSE.
    degen_idx = 0

    DO ir0pt = 1, self%nR0pt
      DO jw = 1, Nw
        DO iw = 1, Nw
          iuw = self%shift_map_inv(iw, jw)
          dist_min = SUM((self%R0_cart(:, ir0pt) + self%shift_cart_u(:, iuw))**2)
          minRT(:, iuw, ir0pt) = self%R0_red(:, ir0pt)
          ndegen(iuw, ir0pt) = 0
          ! evaluate minimum distance for R0+T+r_n-r_m
          DO i1 = -cell_expand(1), cell_expand(1)
            DO i2 = -cell_expand(2), cell_expand(2)
              DO i3 = -cell_expand(3), cell_expand(3)
                ndiff_red(1) = self%R0_red(1, ir0pt) + i1*w90data%k_grid(1)
                ndiff_red(2) = self%R0_red(2, ir0pt) + i2*w90data%k_grid(2)
                ndiff_red(3) = self%R0_red(3, ir0pt) + i3*w90data%k_grid(3)
                CALL red2cart_real(ndiff_red, ndiff_cart)
                IF (SUM((ndiff_cart + self%shift_cart_u(:, iuw))**2) < dist_min) THEN
                  dist_min = SUM((ndiff_cart + self%shift_cart_u(:, iuw))**2)
                END IF
              END DO
            END DO
          END DO
          ! select R0+T
          DO i1 = -cell_expand(1), cell_expand(1)
            DO i2 = -cell_expand(2), cell_expand(2)
              DO i3 = -cell_expand(3), cell_expand(3)
                ! Calculate distance squared |r-R|^2
                ndiff_red(1) = self%R0_red(1, ir0pt) + i1*w90data%k_grid(1)
                ndiff_red(2) = self%R0_red(2, ir0pt) + i2*w90data%k_grid(2)
                ndiff_red(3) = self%R0_red(3, ir0pt) + i3*w90data%k_grid(3)
                CALL red2cart_real(ndiff_red, ndiff_cart)
                IF (ABS(SQRT(SUM((ndiff_cart + self%shift_cart_u(:, iuw))**2)) &
                        - SQRT(dist_min)) < 1E-5) THEN
                  idx = (((i3 + cell_expand(3)) &
                          *cell_range(2) + (i2 + cell_expand(2))) &
                         *cell_range(1) + (i1 + cell_expand(1))) &
                        *self%nR0pt + ir0pt
                  bRvec_selected(idx) = .TRUE.
                  R_red(:, idx) = ndiff_red(:)
                  ndegen(iuw, ir0pt) = ndegen(iuw, ir0pt) + 1
                  degen_idx(iuw, idx) = ndegen(iuw, ir0pt)
                  !< idx is iuw's ndegen-th degenerate R vector for given R0
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    ! count number of R vectors
    self%nRpt = COUNT(bRvec_selected)
    self%max_degen = MAXVAL(ndegen)
    WRITE (stdout, '(2X, A, I0)') '- Number of selected R vectors: ', self%nRpt
    ALLOCATE (self%R_red(3, self%nRpt))
    ALLOCATE (self%R_cart(3, self%nRpt))
    ALLOCATE (self%w_R(Nw, Nw, self%nRpt))
    ALLOCATE (self%map_r02r(self%nR0pt, self%nu_shift, self%max_degen))
    ALLOCATE (self%map_r2r0(2, self%nu_shift, self%nRpt))
    irpt = 0
    DO idx = 1, nRpt
      IF (bRvec_selected(idx)) THEN
        irpt = irpt + 1
        self%R_red(:, irpt) = R_red(:, idx)
        CALL red2cart_real(self%R_red(:, irpt), self%R_cart(:, irpt))
        ir0pt = MOD(idx - 1, self%nR0pt) + 1
        DO jw = 1, Nw
          DO iw = 1, Nw
            iuw = self%shift_map_inv(iw, jw)
            i = degen_idx(iuw, idx)
            IF (i == 0) THEN
              self%w_R(iw, jw, irpt) = 0.0D0
            ELSE
              self%w_R(iw, jw, irpt) = 1.0D0/REAL(ndegen(iuw, ir0pt), DP)*w_R0(ir0pt)
              self%map_r02r(ir0pt, iuw, i) = irpt
            END IF
            self%map_r2r0(1, iuw, irpt) = ir0pt
            self%map_r2r0(2, iuw, irpt) = i
          END DO
        END DO
      END IF
    END DO

    WRITE (stdout, '(2X, A)') '- Finished building R vectors.'
  END SUBROUTINE build_ws

  SUBROUTINE bcast_ws(self)
    USE mp_base, ONLY: mp_bcast
    USE io_global, ONLY: ionode
    USE system, ONLY: Nw
    CLASS(R_vec_type), INTENT(inout) :: self
    CALL self%bcast_shift()
    !
    CALL mp_bcast(self%max_degen)
    CALL mp_bcast(self%nR0pt)
    CALL mp_bcast(self%nRpt)
    IF (.NOT. ionode) THEN
      ALLOCATE (self%R0_red(3, self%nR0pt))
      ALLOCATE (self%R0_cart(3, self%nR0pt))
      ALLOCATE (self%R_red(3, self%nRpt))
      ALLOCATE (self%R_cart(3, self%nRpt))
      ALLOCATE (self%w_R(Nw, Nw, self%nRpt))
      ALLOCATE (self%map_r02r(self%nR0pt, self%nu_shift, self%max_degen))
      ALLOCATE (self%map_r2r0(2, self%nu_shift, self%nRpt))
    END IF
    CALL mp_bcast(self%R0_red)
    CALL mp_bcast(self%R0_cart)
    CALL mp_bcast(self%R_red)
    CALL mp_bcast(self%R_cart)
    CALL mp_bcast(self%w_R)
    CALL mp_bcast(self%map_r02r)
    CALL mp_bcast(self%map_r2r0)
  END SUBROUTINE bcast_ws
END MODULE R_vector
