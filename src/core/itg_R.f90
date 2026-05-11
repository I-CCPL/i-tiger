MODULE itg_R
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  USE wannier90, ONLY: w90data
  USE R_vector, ONLY: R_vec_type
  IMPLICIT NONE
  TYPE(R_vec_type)::R_vec

  !> intermediate data in R space
  TYPE::R_data_type
    COMPLEX(DP), ALLOCATABLE::mH_R(:, :, :)
    !< Hamiltonian in R space (Nw, Nw, nRpt)
    LOGICAL::bA_R
    COMPLEX(DP), ALLOCATABLE::mA_R(:, :, :, :)
    !< Berry connection in R space (Nw, Nw, nRpt, 3)
  END TYPE R_data_type
  TYPE(R_data_type)::R_data
CONTAINS
  SUBROUTINE set_R_flag()
    USE f_params, ONLY: lOAM, lBerry, lBCD, lNLO
    R_data%bA_R = (lOAM .OR. lBerry .OR. lBCD .OR. lNLO)
  END SUBROUTINE set_R_flag
  !
  SUBROUTINE make_R()
    USE constants, ONLY: cmplx_0
    USE fft_base, ONLY: fft_q2R
    USE io_global, ONLY: ionode, stdout
    USE f_params, ONLY: lBCD, lNLO, convention
    INTEGER::a
    IF (.NOT. ionode) RETURN
    CALL start_clock('make_R')
    WRITE (stdout, '(2X, A)') 'Building data in R space...'
    !
    SELECT CASE (convention)
    CASE (0, 1)
      CALL R_vec%build_ws(w90data)
    CASE (2)
      CALL R_vec%build_R(w90data)
    CASE default
      CALL errore(1, 'make_R_data', 'invalid convention')
    END SELECT
    !
    ALLOCATE (R_data%mH_R(Nw, Nw, R_vec%nRpt))
    CALL fft_q2R(w90data, R_vec, w90data%Hq, R_data%mH_R)

    IF (R_data%bA_R) THEN
      ALLOCATE (R_data%mA_R(Nw, Nw, R_vec%nRpt, 3))
      DO a = 1, 3
        CALL fft_q2R(w90data, R_vec, &
                     w90data%Aq(:, :, :, a), R_data%mA_R(:, :, :, a))
      END DO
      ! CALL enforce_Hemiticity_R(R_data%mA_R_b, R_data%mA_R)
    END IF
    CALL write_sep_line()
    CALL stop_clock('make_R')
  END SUBROUTINE make_R
  !
  SUBROUTINE bcast_R()
    USE mp_base, ONLY: mp_bcast
    USE io_global, ONLY: ionode
    USE f_params, ONLY: convention
    !
    SELECT CASE (convention)
    CASE (0, 1)
      CALL R_vec%bcast_ws()
    CASE (2)
      CALL R_vec%bcast_R()
    END SELECT

    IF (.NOT. ionode) THEN
      ALLOCATE (R_data%mH_R(Nw, Nw, R_vec%nRpt))
      IF (R_data%bA_R) THEN
        ALLOCATE (R_data%mA_R(Nw, Nw, R_vec%nRpt, 3))
      END IF
    END IF
    !
    CALL mp_bcast(R_data%mH_R)
    IF (R_data%bA_R) THEN
      CALL mp_bcast(R_data%mA_R)
    END IF
  END SUBROUTINE bcast_R
  !
  SUBROUTINE enforce_Hemiticity_R(mat_in, mat_out)
    USE kinds, ONLY: DP, eq_vec_real
    USE constants, ONLY: cmplx_0
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
  SUBROUTINE clear_R()
    IF (ALLOCATED(R_data%mH_R)) DEALLOCATE (R_data%mH_R)
    IF (ALLOCATED(R_data%mA_R)) DEALLOCATE (R_data%mA_R)
    CALL R_vec%clear()
  END SUBROUTINE clear_R
END MODULE itg_R
