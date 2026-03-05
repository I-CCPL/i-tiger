MODULE R_vector
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  TYPE::R_vec_type
    INTEGER::nu_shift
    !< Number of unique Wannier center shifts
    REAL(DP), ALLOCATABLE :: shift_cart(:, :)
    !< Unique list of Wannier center shifts in Cartesian coordinates (3, nu_shift)
    INTEGER, ALLOCATABLE:: shift_map_inv(:, :)
    !< Map from Wannier function pair (iw, jw) to shift index

    INTEGER::nrpt
    !< Assume nrpt == nkpt
    INTEGER::R0_grid(3)
    REAL(DP), ALLOCATABLE::R0_red(:, :)
    !< R0_vec in reduced coordinates (3, nrpt)
    REAL(DP), ALLOCATABLE::R0_cart(:, :)
    !< R0_vec in Cartesian coordinates (3, nrpt)

    INTEGER::ncell
    !< Number of T
    REAL(DP), ALLOCATABLE::R_red(:, :, :, :)
    !< R_vec in reduced coordinates (3,ncell, nu_shift, nrpt)
    REAL(DP), ALLOCATABLE::R_cart(:, :, :, :)
    !< R_vec in reduced coordinates (3,ncell, nu_shift, nrpt)
    INTEGER, ALLOCATABLE::nRvec(:, :)
    !< Number of R_vec for each shift & R (nu_shift, nrpt)
  CONTAINS
    PROCEDURE::clear => clear_Rvec
    PROCEDURE::build_shift => build_shift_vecs
    PROCEDURE::build_R => build_Rvecs
    PROCEDURE::fft_q2R => Rvec_fft_q2R
    PROCEDURE::fft_R2k => Rvec_fft_R2k
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
    IF (ALLOCATED(self%nRvec)) DEALLOCATE (self%nRvec)
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
        all_shift(:, iw + (jw - 1)*Nw) = wannier_center_cart(:, iw) - wannier_center_cart(:, jw)
      END DO
    END DO

    ! reduce shifts to unique ones
    CALL unique_vec3_inv(all_shift, 1D-8, self%shift_cart, self%nu_shift, shift_map_inv)
    WRITE (stdout, '(2X, A, I0)') '- Number of unique shifts: ', self%nu_shift

    ALLOCATE (self%shift_map_inv(Nw, Nw))
    DO iw = 1, Nw
      DO jw = 1, Nw
        self%shift_map_inv(iw, jw) = shift_map_inv(iw + (jw - 1)*Nw)
      END DO
    END DO
  END SUBROUTINE build_shift_vecs

  SUBROUTINE build_Rvecs(self, w90data)
    USE kinds, ONLY: eq_real
    USE constants, ONLY: vec_0
    USE io_global, ONLY: stdout, write_sep_line
    USE wannier90, ONLY: w90data_type
    USE system, ONLY: Nw, red2cart_real
    CLASS(R_vec_type), INTENT(INOUT)::self
    TYPE(w90data_type), INTENT(INOUT)::w90data
    INTEGER::iuw, iRvec
    INTEGER::icell, cell_x, cell_y, cell_z
    INTEGER::irpt, rpt_x, rpt_y, rpt_z
    REAL(DP)::R_vec(3), R_dist, dist_min
    REAL(DP), ALLOCATABLE::Tvec_red(:, :), Tvec_cart(:, :)
    REAL(DP)::shift_red(3)
    ! TODO: cell_expand from input
    INTEGER::cell_expand(3) = (/1, 1, 1/)
    INTEGER::cell_range(3)
    REAL(DP)::dr_nmR0(3) !< R0 + r_m - r_n

    self%nrpt = w90data%kpts%nkpt

    CALL self%build_shift(w90data%wannier_center_cart)

    WRITE (stdout, '(2X, A)') 'Building R vectors for Fourier transform...'
    ! Build T vectors
    cell_range(:) = 2*cell_expand(:) + 1
    self%ncell = PRODUCT(cell_range)
    WRITE (stdout, '(2X, A)') '- Building T vectors for periodic images...'
    WRITE (stdout, '(2X, A, I0)') '- Number of considered T vectors: ', self%ncell
    ALLOCATE (Tvec_red(3, self%ncell))
    ALLOCATE (Tvec_cart(3, self%ncell))
    DO icell = 1, self%ncell
      CALL grid_idx2xyz(cell_range, icell, cell_x, cell_y, cell_z, -cell_expand)
      Tvec_red(:, icell) = (/cell_x*w90data%k_grid(1), &
                             cell_y*w90data%k_grid(2), &
                             cell_z*w90data%k_grid(3)/)
      CALL red2cart_real(Tvec_red(:, icell), Tvec_cart(:, icell))
    END DO

    WRITE (stdout, '(2X, A)') '- Building R vectors for each shift and R0...'
    ! Find T such that minimize |R0+T+r_n-r_m| for given (R0, n, m)
    ALLOCATE (self%R0_red(3, self%nrpt))
    ALLOCATE (self%R0_cart(3, self%nrpt))
    ALLOCATE (self%nRvec(self%nu_shift, self%nrpt))
    ALLOCATE (self%R_red(3, self%ncell, self%nu_shift, self%nrpt))
    ALLOCATE (self%R_cart(3, self%nrpt*self%ncell, self%nu_shift, self%nrpt))
    self%R0_grid(:) = w90data%k_grid(:)

    DO irpt = 1, self%nrpt
      ! Build R0 vector
      CALL grid_idx2xyz(self%R0_grid, irpt, rpt_x, rpt_y, rpt_z, vec_0)
      self%R0_red(:, irpt) = (/REAL(rpt_x, DP), REAL(rpt_y, DP), REAL(rpt_z, DP)/)
      CALL red2cart_real(self%R0_red(:, irpt), self%R0_cart(:, irpt))

      DO iuw = 1, self%nu_shift
        dr_nmR0(:) = self%R0_cart(:, irpt) + self%shift_cart(:, iuw)

        ! minimize |R0+T+r_n-r_m| by searching T vectors
        dist_min = 1.0D10
        DO icell = 1, self%ncell
          R_vec = Tvec_cart(:, icell) + dr_nmR0
          R_dist = SQRT(SUM(R_vec**2))
          IF (R_dist < dist_min) THEN
            dist_min = R_dist
          END IF
        END DO

        ! R vectors for each T vector with distance close to dist_min
        iRvec = 0
        DO icell = 1, self%ncell
          R_vec = Tvec_cart(:, icell) + dr_nmR0
          R_dist = SQRT(SUM(R_vec**2))
          IF (eq_real(R_dist, dist_min, 1.0D-6)) THEN
            iRvec = iRvec + 1
            self%R_red(:, iRvec, iuw, irpt) = Tvec_red(:, icell) + self%R0_red(:, irpt)
            self%R_cart(:, iRvec, iuw, irpt) = R_vec
          END IF
        END DO
        self%nRvec(iuw, irpt) = iRvec
      END DO
    END DO
    CALL write_sep_line()
  END SUBROUTINE build_Rvecs

  SUBROUTINE Rvec_fft_q2R(self, w90data, X_q, X_R)
    USE kinds, ONLY: DP
    USE io_global, ONLY: stdout
    USE constants, ONLY: tpi
    USE system, ONLY: Nw
    USE wannier90, ONLY: w90data_type
    CLASS(R_vec_type), INTENT(INOUT)::self
    TYPE(w90data_type), INTENT(IN)::w90data
    COMPLEX(DP), INTENT(IN)::X_q(Nw, Nw, w90data%kpts%nkpt)
    COMPLEX(DP), ALLOCATABLE, INTENT(OUT)::X_R(:, :, :)
    INTEGER::iw, jw, iuw, irpt, ikpt
    REAL(DP)::wk, phase
    COMPLEX(DP)::exp_phase

    WRITE (stdout, '(2X, A)') '- Performing Fourier transform from q to R space...'
    IF (.NOT. ALLOCATED(X_R)) ALLOCATE (X_R(Nw, Nw, self%nrpt))
    wk = 1/REAL(w90data%kpts%nkpt, DP)
    DO iw = 1, Nw
      DO jw = 1, Nw
        DO irpt = 1, self%nrpt
          X_R(iw, jw, irpt) = 0.0D0
          DO ikpt = 1, w90data%kpts%nkpt
            phase = -tpi*DOT_PRODUCT(w90data%kpts%k_red(:, ikpt), self%R0_red(:, irpt))
            exp_phase = CMPLX(COS(phase), SIN(phase))
            X_R(iw, jw, irpt) = X_R(iw, jw, irpt) + wk*exp_phase*X_q(iw, jw, ikpt)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE Rvec_fft_q2R

  SUBROUTINE Rvec_fft_R2k(self, k_list, X_R, X_k)
    USE kinds, ONLY: DP
    USE constants, ONLY: tpi
    USE io_global, ONLY: stdout
    USE system, ONLY: Nw
    USE kpoints, ONLY: kpoint_type
    CLASS(R_vec_type), INTENT(INOUT)::self
    TYPE(kpoint_type), INTENT(IN)::k_list
    COMPLEX(DP), INTENT(IN)::X_R(Nw, Nw, self%nrpt)
    COMPLEX(DP), INTENT(OUT)::X_k(Nw, Nw, k_list%nkpt)
    INTEGER::iw, jw, iuw, ikpt, irpt, irtpt
    REAL(DP)::phase, degen
    COMPLEX(DP)::exp_phase

    WRITE (stdout, '(2X, A)') '- Performing Fourier transform from R to k space...'
    DO iw = 1, Nw
      DO jw = 1, Nw
        iuw = self%shift_map_inv(iw, jw)
        DO ikpt = 1, k_list%nkpt
          X_k(iw, jw, ikpt) = 0.0D0
          DO irpt = 1, self%nrpt
            degen = 1/REAL(self%nRvec(iuw, irpt), DP)
            DO irtpt = 1, self%nRvec(iuw, irpt)
              phase = tpi*DOT_PRODUCT(k_list%k_red(:, ikpt), self%R_red(:, irtpt, iuw, irpt))
              exp_phase = CMPLX(COS(phase), SIN(phase))
              X_k(iw, jw, ikpt) = X_k(iw, jw, ikpt) + degen*exp_phase*X_R(iw, jw, irpt)
            END DO
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE Rvec_fft_R2k
END MODULE R_vector
