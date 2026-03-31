MODULE itg_R
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  USE wannier90, ONLY: w90data
  USE R_vector, ONLY: R_vec_type
  IMPLICIT NONE
  TYPE(R_vec_type)::R_vec

  COMPLEX(DP), ALLOCATABLE::H_R(:, :, :), dH_R(:, :, :, :)
  COMPLEX(DP), ALLOCATABLE::A_R(:, :, :, :), A_R_b(:, :, :, :)
  COMPLEX(DP), ALLOCATABLE::dA_R(:, :, :, :, :)
CONTAINS
  SUBROUTINE make_R_data()
    USE constants, ONLY: zero
    USE fft_base, ONLY: fft_q2R
    USE der_base, ONLY: der_R
    USE io_input, ONLY: lShift
    USE wannier90, ONLY: lreq_mmn
    INTEGER::inb
    CALL start_clock('make_R_data')
    !
    CALL R_vec%build_R(w90data)
    !
    ALLOCATE (H_R(Nw, Nw, R_vec%nRpt))
    CALL fft_q2R(w90data, R_vec, w90data%Hq, H_R)

    IF (lreq_mmn) THEN
      ALLOCATE (A_R(3, Nw, Nw, R_vec%nRpt))
      ALLOCATE (A_R_b(3, Nw, Nw, R_vec%nRpt))
      A_R = zero
      DO inb = 1, w90data%nnb
        CALL fft_q2R(w90data, R_vec, w90data%Aq(:, :, :, :, inb), A_R_b)
        CALL enforce_Hemiticity_R(A_R_b, A_R)
      END DO

      ! IF (lShift) THEN
      !   ALLOCATE (dA_R(3, 3, Nw, Nw, R_vec%nRpt))
      !   CALL der_R(R_vec, A_R, dA_R)
      ! END IF

      ALLOCATE (dH_R(3, Nw, Nw, R_vec%nRpt))
      dH_R = zero
      CALL der_R(R_vec, H_R, dH_R)
    END IF
    CALL stop_clock('make_R_data')
  END SUBROUTINE make_R_data
  !
  SUBROUTINE enforce_Hemiticity_R(mat_in, mat_out)
    USE kinds, ONLY: DP, eq_vec_real
    USE constants, ONLY: zero
    COMPLEX(DP), INTENT(IN)::mat_in(:, :, :, :)
    COMPLEX(DP), INTENT(OUT)::mat_out(:, :, :, :)
    INTEGER::iw, jw, irpt, jrpt
    REAL(DP)::tmp_vec(3)
    !
    DO jrpt = 1, R_vec%nRpt
      tmp_vec = -R_vec%R_red(:, jrpt)
      DO irpt = 1, R_vec%nRpt
        IF (.NOT. eq_vec_real(R_Vec%R_red(:, irpt), &
                              tmp_vec, 1.0D-12)) THEN
          CYCLE
        END IF
        !
        DO jw = 1, Nw
          DO iw = 1, Nw
            ! Enforce Hermiticity. mat_out = (mat_in + mat_in^dagger)/2
            mat_out(:, iw, jw, irpt) = mat_out(:, iw, jw, irpt) + &
                                       (mat_in(:, iw, jw, irpt) + CONJG(mat_in(:, jw, iw, jrpt)))/2
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE enforce_Hemiticity_R
  !
  SUBROUTINE clear_R_data()
    CALL R_vec%clear()
    IF (ALLOCATED(H_R)) DEALLOCATE (H_R)
    IF (ALLOCATED(dH_R)) DEALLOCATE (dH_R)
    IF (ALLOCATED(A_R)) DEALLOCATE (A_R)
    IF (ALLOCATED(dA_R)) DEALLOCATE (dA_R)
    IF (ALLOCATED(A_R_b)) DEALLOCATE (A_R_b)
  END SUBROUTINE clear_R_data
END MODULE itg_R
