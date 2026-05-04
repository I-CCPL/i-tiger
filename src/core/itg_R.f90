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
  COMPLEX(DP), ALLOCATABLE::d2H_R(:, :, :, :, :)
  COMPLEX(DP), ALLOCATABLE::O_R(:, :, :, :), dO_R(:, :, :, :, :)
CONTAINS
  SUBROUTINE make_R_data()
    USE constants, ONLY: zero
    USE fft_base, ONLY: fft_q2R
    USE der_base, ONLY: der_R
    USE io_input, ONLY: lBCD, lNLO, lreq_mmn, convention
    INTEGER::iw, jw, irpt, a, b, c
    CALL start_clock('make_R_data')
    !
    SELECT CASE (convention)
    CASE (0)
      CALL R_vec%build_ws(w90data)
    CASE (1)
      CALL R_vec%build_ws(w90data)
    CASE (2)
      CALL R_vec%build_R(w90data)
    CASE default
      CALL errore(1, 'make_R_data', 'invalid convention')
    END SELECT
    !
    ALLOCATE (H_R(Nw, Nw, R_vec%nRpt))
    CALL fft_q2R(w90data, R_vec, w90data%Hq, H_R)!, dH_R, shift)

    IF (lreq_mmn) THEN
      ALLOCATE (A_R(Nw, Nw, R_vec%nRpt, 3))
      A_R = zero
      CALL fft_q2R(w90data, R_vec, w90data%Aq(:, :, :, :), A_R)
      ! CALL enforce_Hemiticity_R(A_R_b, A_R)

      ALLOCATE (dH_R(Nw, Nw, R_vec%nRpt, 3))
      dH_R = zero
      CALL der_R(R_vec, H_R, dH_R)

      ALLOCATE (dA_R(Nw, Nw, R_vec%nRpt, 3, 3))
      ALLOCATE (O_R(Nw, Nw, R_vec%nRpt, 3))
      CALL der_R(R_vec, A_R, dA_R)
      DO irpt = 1, R_vec%nRpt
        DO iw = 1, Nw
          DO jw = 1, Nw
            DO c = 1, 3
              a = MOD(c, 3) + 1
              b = MOD(a, 3) + 1
              O_R(iw, jw, irpt, c) = dA_R(iw, jw, irpt, b, a) - dA_R(iw, jw, irpt, a, b)
            END DO
          END DO
        END DO
      END DO

      IF (lBCD) THEN
        ALLOCATE (d2H_R(Nw, Nw, R_vec%nRpt, 3, 3))
        ! ALLOCATE (dA_R(Nw, Nw, R_vec%nRpt, 3, 3))
        ! ALLOCATE (O_R(Nw, Nw, R_vec%nRpt, 3))
        ALLOCATE (dO_R(Nw, Nw, R_vec%nRpt, 3, 3))
        CALL der_R(R_vec, dH_R, d2H_R)
        ! CALL der_R(R_vec, A_R, dA_R, shift)

        ! DO irpt = 1, R_vec%nRpt
        !   DO iw = 1, Nw
        !     DO jw = 1, Nw
        !       DO c = 1, 3
        !         a = MOD(c, 3) + 1
        !         b = MOD(a, 3) + 1
        !         O_R(iw, jw, irpt, c) = dA_R(iw, jw, irpt, b, a) - dA_R(iw, jw, irpt, a, b)
        !       END DO
        !     END DO
        !   END DO
        ! END DO
        CALL der_R(R_vec, O_R, dO_R)
      END IF

      IF (lNLO) THEN
        ALLOCATE (d2H_R(Nw, Nw, R_vec%nRpt, 3, 3))
        ALLOCATE (dA_R(Nw, Nw, R_vec%nRpt, 3, 3))
        CALL der_R(R_vec, dH_R, d2H_R)
        CALL der_R(R_vec, A_R, dA_R)
      END IF
    END IF
    CALL write_sep_line()
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
            mat_out(iw, jw, irpt, :) = mat_out(iw, jw, irpt, :) + &
                                       (mat_in(iw, jw, irpt, :) + CONJG(mat_in(jw, iw, jrpt, :)))/2
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
    IF (ALLOCATED(d2H_R)) DEALLOCATE (d2H_R)
    IF (ALLOCATED(A_R_b)) DEALLOCATE (A_R_b)
    IF (ALLOCATED(O_R)) DEALLOCATE (O_R)
    IF (ALLOCATED(dO_R)) DEALLOCATE (dO_R)
  END SUBROUTINE clear_R_data
END MODULE itg_R
