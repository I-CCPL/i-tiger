MODULE itg_R
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  USE wannier90, ONLY: w90data
  USE R_vector, ONLY: R_vec_type
  IMPLICIT NONE
  TYPE(R_vec_type)::R_vec

  COMPLEX(DP), ALLOCATABLE::H_R(:, :, :), dH_R(:, :, :, :)
  COMPLEX(DP), ALLOCATABLE::A_R(:, :, :, :), A_R_b(:, :, :, :)
CONTAINS
  SUBROUTINE make_R_data()
    USE kinds, ONLY: eq_vec_real
    USE constants, ONLY: zero
    USE fft_base, ONLY: fft_q2R
    USE debug_data, ONLY: write_matrix, write_diag_matrix
    USE der_base, ONLY: der_R
    INTEGER::inb, ikpt, irpt, jrpt, iw, jw
    REAL(DP)::tmp_vec(3)
    CALL start_clock('make_R_data')
    !
    CALL R_vec%build_R(w90data)
    !
    ALLOCATE (H_R(Nw, Nw, R_vec%nRpt))
    CALL fft_q2R(w90data, R_vec, w90data%Hq, H_R)
    ! CALL write_matrix('H_R.itg', H_R, R_vec%nRpt, 1)

    ALLOCATE (A_R(3, Nw, Nw, R_vec%nRpt))
    ALLOCATE (A_R_b(3, Nw, Nw, R_vec%nRpt))
    A_R = zero
    DO inb = 1, w90data%nnb
      CALL fft_q2R(w90data, R_vec, w90data%Aq(:, :, :, :, inb), A_R_b)

      ikpt = 0
      DO jrpt = 1, R_vec%nRpt
        tmp_vec = -R_vec%R_red(:, jrpt)
        DO irpt = 1, R_vec%nRpt
          IF (.NOT. eq_vec_real(R_Vec%R_red(:, irpt), &
                                tmp_vec, 1.0D-12)) THEN
            CYCLE
          END IF
          ikpt = ikpt + 1
          DO jw = 1, Nw
            DO iw = 1, Nw
              ! Enforce Hermiticity of A_R
              A_R(:, iw, jw, irpt) = A_R(:, iw, jw, irpt) + (A_R_b(:, iw, jw, irpt) + CONJG(A_R_b(:, jw, iw, jrpt)))/2
            END DO
          END DO
        END DO
      END DO
    END DO
    ! CALL write_matrix('A_R.itg', A_R, R_vec%nRpt, 1)

    ALLOCATE (dH_R(3, Nw, Nw, R_vec%nrpt))
    CALL der_R(R_vec, H_R, dH_R)
    ! CALL write_matrix('dH_R.itg', dH_R, R_vec%nrpt, 1)
    CALL stop_clock('make_R_data')
  END SUBROUTINE make_R_data
  !
  SUBROUTINE clear_R_data()
    CALL R_vec%clear()
    IF (ALLOCATED(H_R)) DEALLOCATE (H_R)
    IF (ALLOCATED(dH_R)) DEALLOCATE (dH_R)
    IF (ALLOCATED(A_R)) DEALLOCATE (A_R)
    IF (ALLOCATED(A_R_b)) DEALLOCATE (A_R_b)
  END SUBROUTINE clear_R_data
END MODULE itg_R
