MODULE R_vector
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  !
  TYPE::R_vec_type
    INTEGER::nu_shift
    !< Number of unique Wannier center shifts
    REAL(DP), ALLOCATABLE :: shift_cart(:, :)
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
  CONTAINS
    PROCEDURE::clear => clear_Rvec
    PROCEDURE::build_shift => build_shift_vecs
    PROCEDURE::build_R => build_Rvecs
  END TYPE R_vec_type
  PUBLIC::R_vec_type
  !
CONTAINS
  SUBROUTINE clear_Rvec(self)
    CLASS(R_vec_type), INTENT(INOUT)::self
    IF (ALLOCATED(self%shift_cart)) DEALLOCATE (self%shift_cart)
    IF (ALLOCATED(self%shift_map_inv)) DEALLOCATE (self%shift_map_inv)
    IF (ALLOCATED(self%R0_red)) DEALLOCATE (self%R0_red)
    IF (ALLOCATED(self%R0_cart)) DEALLOCATE (self%R0_cart)
    IF (ALLOCATED(self%R_red)) DEALLOCATE (self%R_red)
    IF (ALLOCATED(self%R_cart)) DEALLOCATE (self%R_cart)
  END SUBROUTINE clear_Rvec
  !
  SUBROUTINE build_shift_vecs(self, wannier_center_cart)
    !< Build Wannier center shift (r_m-r_n) in cartesian coordinates
    USE io_global, ONLY: stdout
    USE system, ONLY: Nw
    USE algo_unique, ONLY: unique_vec3_inv
    CLASS(R_vec_type), INTENT(INOUT) :: self
    REAL(DP), INTENT(IN) :: wannier_center_cart(3, Nw)
    REAL(DP)::all_shift(3, Nw*Nw)
    INTEGER::shift_map_inv(Nw*Nw)
    INTEGER::iw, jw
    !
    IF (ALLOCATED(self%shift_cart)) RETURN
    WRITE (stdout, '(2X, A)') 'Building Wannier center shift vectors...'
    !
    DO jw = 1, Nw
      DO iw = 1, Nw
        all_shift(:, iw + (jw - 1)*Nw) = -wannier_center_cart(:, iw) + wannier_center_cart(:, jw)
      END DO
    END DO

    ! reduce shifts to unique ones
    CALL unique_vec3_inv(all_shift, 1D-8, self%shift_cart, self%nu_shift, shift_map_inv)
    WRITE (stdout, '(2X, A, I0)') '- Number of unique shifts: ', self%nu_shift

    ALLOCATE (self%shift_map_inv(Nw, Nw))
    DO jw = 1, Nw
      DO iw = 1, Nw
        self%shift_map_inv(iw, jw) = shift_map_inv(iw + (jw - 1)*Nw)
      END DO
    END DO
  END SUBROUTINE build_shift_vecs

  SUBROUTINE build_Rvecs(self, w90data)
    USE kinds, ONLY: eq_real
    USE constants, ONLY: vec_0
    USE io_global, ONLY: stdout, write_sep_line
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
    ! TODO: cell_expand from input
    INTEGER::cell_expand(3) = (/1, 1, 1/)
    INTEGER::cell_range(3)
    REAL(DP)::dr_nmR0(3) !< R0 + r_m - r_n
    LOGICAL, ALLOCATABLE::bRvec_selected(:)
    !< whether Rvec is selected for given R0

    self%nR0pt = w90data%kpts%nkpt

    CALL self%build_shift(w90data%wannier_center_cart)

    CALL write_sep_line()
    WRITE (stdout, '(2X, A)') 'Building R vectors for Fourier transform...'
    ! Build T vectors
    cell_range(:) = 2*cell_expand(:) + 1
    ncell = PRODUCT(cell_range)
    WRITE (stdout, '(2X, A)') '- Building T vectors for periodic images...'
    WRITE (stdout, '(2X, A, I0)') '- Number of considered T vectors: ', ncell
    ALLOCATE (Tvec_red(3, ncell))
    ALLOCATE (Tvec_cart(3, ncell))
    DO icell = 1, ncell
      CALL grid_idx2xyz(cell_range, icell, cell_x, cell_y, cell_z, -cell_expand)
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
          R_vector = Rvec_cart(:, irpt) + self%shift_cart(:, iuw)
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
            R_vector = Rvec_cart(:, irpt) + self%shift_cart(:, iuw)
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
    ! unique R vectors are selected
    self%nRpt = 0
    DO irpt = 1, nRpt
      IF (bRvec_selected(irpt)) THEN
        self%nRpt = self%nRpt + 1
      END IF
    END DO
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

    CALL write_sep_line()
  END SUBROUTINE build_Rvecs
END MODULE R_vector
