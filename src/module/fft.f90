MODULE fft_base
  USE kinds, ONLY: DP
  USE constants, ONLY: cmplx_0, tpi, cmplx_i
  USE io_global, ONLY: stdout
  USE f_params, ONLY: FFT_conv
  USE system, ONLY: Nw
  USE wannier90, ONLY: w90data_type
  USE R_vector, ONLY: R_vec_type
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  PRIVATE
  REAL(DP), ALLOCATABLE::shift_cart(:, :, :)
  REAL(DP), ALLOCATABLE::shift_red(:, :, :)
  PUBLIC::fft_init, fft_q2R, fft_R2k, fft_R2k_vec
CONTAINS
  SUBROUTINE fft_init(R_vec, A_R)
    USE system, ONLY: cart2red_real
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(INOUT) :: A_R(Nw, Nw, R_vec%nRpt, 3)
    REAL(DP)::center(3, Nw)
    INTEGER::irpt, iw, jw
    IF (TRIM(FFT_conv) /= 'TB') RETURN
    ALLOCATE (shift_cart(3, Nw, Nw))
    ALLOCATE (shift_red(3, Nw, Nw))
    DO irpt = 1, R_vec%nRpt
      IF (ALL(R_vec%R_red(:, irpt) == 0)) THEN
        DO iw = 1, Nw
          center(:, iw) = DBLE(A_R(iw, iw, irpt, :))
        END DO
        EXIT
      END IF
    END DO

    DO irpt = 1, R_vec%nRpt
      IF (ALL(R_vec%R_red(:, irpt) == 0)) THEN
        DO iw = 1, Nw
          A_R(iw, iw, irpt, :) = cmplx_0
        END DO
        EXIT
      END IF
    END DO

    DO iw = 1, Nw
      DO jw = 1, Nw
        shift_cart(:, iw, jw) = center(:, jw) - center(:, iw)
      END DO
    END DO
    CALL cart2red_real(shift_cart, shift_red, Nw, Nw)
  END SUBROUTINE fft_init
  !
  SUBROUTINE fft_q2R(w90data, R_vec, X_q, X_R)
    TYPE(w90data_type), INTENT(IN) :: w90data
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_q(:, :, :)
    COMPLEX(DP), INTENT(OUT) :: X_R(:, :, :)
    INTEGER::iw, jw, irpt, ikpt
    REAL(DP)::phase, Rvec(3)
    COMPLEX(DP)::exp_phase, fac

    ! IF (ANY(SHAPE(X_q) /= [Nw, Nw, w90data%kpts%nkpt])) THEN
    !   CALL errore(1, 'fft_q2R', 'X_q has wrong shape')
    ! END IF
    ! IF (ANY(SHAPE(X_R) /= [Nw, Nw, R_vec%nRpt])) THEN
    !   CALL errore(1, 'fft_q2R', 'X_R has wrong shape')
    ! END IF

    CALL start_clock('fft_q2R')
    X_R = cmplx_0
    DO irpt = 1, R_vec%nRpt
      DO ikpt = 1, w90data%kpts%nkpt
        phase = tpi*DOT_PRODUCT(w90data%kpts%k_red(:, ikpt), R_vec%R_red(:, irpt))
        exp_phase = EXP(-cmplx_i*phase)
        fac = exp_phase*w90data%kpts%wk
        DO jw = 1, Nw
          DO iw = 1, Nw
            X_R(iw, jw, irpt) = X_R(iw, jw, irpt) + X_q(iw, jw, ikpt)*fac
          END DO
        END DO
      END DO
    END DO
    CALL stop_clock('fft_q2R')
  END SUBROUTINE fft_q2R
  ! ================================================== !
  SUBROUTINE fft_R2k(R_vec, X_R, X_k, dX_k, d2X_k)
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(:, :, :)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: X_k(:, :)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_k(:, :, :)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: d2X_k(:, :, :, :)

    ! IF (ANY(SHAPE(X_R) /= [Nw, Nw, R_vec%nRpt])) THEN
    !   CALL errore(1, 'fft_R2k', 'X_R has wrong shape')
    ! END IF
    ! IF (PRESENT(X_k)) THEN
    !   IF (ANY(SHAPE(X_k) /= [Nw, Nw])) THEN
    !     CALL errore(1, 'fft_R2k', 'X_k has wrong shape')
    !   END IF
    ! END IF
    ! IF (PRESENT(dX_k)) THEN
    !   IF (ANY(SHAPE(dX_k) /= [Nw, Nw, 3])) THEN
    !     CALL errore(1, 'fft_R2k', 'dX_k has wrong shape')
    !   END IF
    ! END IF
    ! IF (PRESENT(d2X_k)) THEN
    !   IF (ANY(SHAPE(d2X_k) /= [Nw, Nw, 3, 3])) THEN
    !     CALL errore(1, 'fft_R2k', 'd2X_k has wrong shape')
    !   END IF
    ! END IF

    CALL start_clock('fft_R2k')
    SELECT CASE (TRIM(FFT_conv))
    CASE ('simple')
      CALL fft_R2k_simple(R_vec, X_R, X_k, dX_k, d2X_k)
    CASE ('TB')
      CALL fft_R2k_TB(R_vec, X_R, X_k, dX_k, d2X_k)
    CASE ('Test')
      CALL fft_R2k_Test(R_vec, X_R, X_k, dX_k, d2X_k)
    END SELECT
    CALL stop_clock('fft_R2k')
  END SUBROUTINE fft_R2k
  ! ================================================== !
  SUBROUTINE fft_R2k_body(R_vec, X_R, X_k, dX_k, d2X_k, &
                          iw, jw, irpt, R_cart)
    !< X_R FFT to X_k and k-derivatives dX_k, d2X_k
    USE kpoints, ONLY: t_iks
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(Nw, Nw, R_vec%nRpt)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: X_k(Nw, Nw)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: d2X_k(Nw, Nw, 3, 3)
    INTEGER, INTENT(IN)::iw, jw, irpt
    REAL(DP)::R_cart(3)
    INTEGER::a, b
    REAL(DP)::phase
    COMPLEX(DP)::exp_phase, fac
    phase = DOT_PRODUCT(t_kpt%k_cart(:, t_iks), R_cart(:))
    exp_phase = EXP(cmplx_i*phase)
    fac = exp_phase*R_vec%w_R(iw, jw, irpt)
    IF (PRESENT(X_k)) X_k(iw, jw) = X_k(iw, jw) + fac*X_R(iw, jw, irpt)
    IF (PRESENT(dX_k)) THEN
      dX_k(iw, jw, :) = dX_k(iw, jw, :) &
                        + cmplx_i*R_cart(:)*fac*X_R(iw, jw, irpt)
    END IF
    IF (PRESENT(d2X_k)) THEN
      DO a = 1, 3
        DO b = 1, 3
          d2X_k(iw, jw, a, b) = d2X_k(iw, jw, a, b) &
                                - R_cart(a)*R_cart(b)*fac*X_R(iw, jw, irpt)
        END DO
      END DO
    END IF
  END SUBROUTINE fft_R2k_body
  ! ================================================== !
  SUBROUTINE fft_R2k_simple(R_vec, X_R, X_k, dX_k, d2X_k)
    USE kpoints, ONLY: t_iks
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(Nw, Nw, R_vec%nRpt)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: X_k(Nw, Nw)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: d2X_k(Nw, Nw, 3, 3)
    INTEGER::iw, jw, irpt, a, b
    REAL(DP)::phase
    COMPLEX(DP)::exp_phase, fac
    !
    IF (PRESENT(X_k)) X_k(:, :) = cmplx_0
    IF (PRESENT(dX_k)) dX_k(:, :, :) = cmplx_0
    IF (PRESENT(d2X_k)) d2X_k(:, :, :, :) = cmplx_0
    DO irpt = 1, R_vec%nRpt
      DO jw = 1, Nw
        DO iw = 1, Nw
          CALL fft_R2k_body(R_vec, X_R, X_k, dX_k, d2X_k, &
                            iw, jw, irpt, R_vec%R_cart(:, irpt))
        END DO
      END DO
    END DO
  END SUBROUTINE fft_R2k_simple
  ! ================================================== !
  SUBROUTINE fft_R2k_TB(R_vec, X_R, X_k, dX_k, d2X_k)
    USE kpoints, ONLY: t_iks
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(Nw, Nw, R_vec%nRpt)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: X_k(Nw, Nw)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: d2X_k(Nw, Nw, 3, 3)
    INTEGER::iw, jw, irpt, a, b, c
    REAL(DP)::phase, R_cart(3)
    COMPLEX(DP)::exp_phase, fac
    !
    IF (PRESENT(X_k)) X_k(:, :) = cmplx_0
    IF (PRESENT(dX_k)) dX_k(:, :, :) = cmplx_0
    IF (PRESENT(d2X_k)) d2X_k(:, :, :, :) = cmplx_0
    DO irpt = 1, R_vec%nRpt
      DO jw = 1, Nw
        DO iw = 1, Nw
          R_cart = R_vec%R_cart(:, irpt) + shift_cart(:, iw, jw)
          CALL fft_R2k_body(R_vec, X_R, X_k, dX_k, d2X_k, &
                            iw, jw, irpt, R_cart)
        END DO
      END DO
    END DO
  END SUBROUTINE fft_R2k_TB
  ! ================================================== !
  SUBROUTINE fft_R2k_Test(R_vec, X_R, X_k, dX_k, d2X_k)
    USE kpoints, ONLY: t_iks
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(Nw, Nw, R_vec%nRpt)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: X_k(Nw, Nw)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: d2X_k(Nw, Nw, 3, 3)
    INTEGER::iw, jw, irpt, a, b, c
    REAL(DP)::phase, R_cart(3)
    COMPLEX(DP)::exp_phase, fac
    !
    IF (PRESENT(X_k)) X_k(:, :) = cmplx_0
    IF (PRESENT(dX_k)) dX_k(:, :, :) = cmplx_0
    IF (PRESENT(d2X_k)) d2X_k(:, :, :, :) = cmplx_0
    DO irpt = 1, R_vec%nRpt
      DO jw = 1, Nw
        DO iw = 1, Nw
          R_cart = R_vec%R_cart(:, irpt) + R_vec%shift_cart(:, iw, jw)
          CALL fft_R2k_body(R_vec, X_R, X_k, dX_k, d2X_k, &
                            iw, jw, irpt, R_cart)
        END DO
      END DO
    END DO
  END SUBROUTINE fft_R2k_Test
  ! ================================================== !
  SUBROUTINE fft_R2k_vec(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k)
    !< X_R FFT to X_k and k-derivatives dX_k, d2X_k, curl_X_k, curl_dX_k
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(Nw, Nw, R_vec%nRpt, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: X_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_k(Nw, Nw, 3, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: d2X_k(Nw, Nw, 3, 3, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: curl_X_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: curl_dX_k(Nw, Nw, 3, 3)
    CALL start_clock('fft_R2k')
    SELECT CASE (TRIM(FFT_conv))
    CASE ('simple')
      CALL fft_R2k_vec_simple(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k)
    CASE ('TB')
      CALL fft_R2k_vec_TB(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k)
    CASE ('Test')
      CALL fft_R2k_vec_Test(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k)
    END SELECT
    CALL stop_clock('fft_R2k')
  END SUBROUTINE fft_R2k_vec
  ! ================================================== !
  SUBROUTINE fft_R2k_vec_body(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k, &
                              iw, jw, irpt, R_cart)
    USE kpoints, ONLY: t_iks
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(Nw, Nw, R_vec%nRpt, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: X_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_k(Nw, Nw, 3, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: d2X_k(Nw, Nw, 3, 3, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: curl_X_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: curl_dX_k(Nw, Nw, 3, 3)
    INTEGER, INTENT(IN) :: iw, jw, irpt
    REAL(DP), INTENT(IN) :: R_cart(3)
    INTEGER::a, b, c, d
    REAL(DP)::phase
    COMPLEX(DP)::exp_phase, fac
    phase = DOT_PRODUCT(t_kpt%k_cart(:, t_iks), R_cart(:))
    exp_phase = EXP(cmplx_i*phase)
    fac = exp_phase*R_vec%w_R(iw, jw, irpt)
    IF (PRESENT(X_k)) THEN
      X_k(iw, jw, :) = X_k(iw, jw, :) + X_R(iw, jw, irpt, :)*fac
    END IF
    IF (PRESENT(dX_k)) THEN
      DO a = 1, 3
        dX_k(iw, jw, a, :) = dX_k(iw, jw, a, :) &
                             + cmplx_i*R_cart(a) &
                             *fac*X_R(iw, jw, irpt, :)
      END DO
    END IF
    IF (PRESENT(d2X_k)) THEN
      DO c = 1, 3
        DO b = 1, 3
          DO a = 1, 3
            d2X_k(iw, jw, c, b, a) = d2X_k(iw, jw, c, b, a) &
                                     - R_cart(a) &
                                     *R_cart(b) &
                                     *fac*X_R(iw, jw, irpt, c)
          END DO
        END DO
      END DO
    END IF
    IF (PRESENT(curl_X_k)) THEN
      DO c = 1, 3
        a = MOD(c, 3) + 1
        b = MOD(a, 3) + 1
        curl_X_k(iw, jw, c) = curl_X_k(iw, jw, c) &
                              + cmplx_i*fac*(R_cart(a)*X_R(iw, jw, irpt, b) &
                                             - R_cart(b)*X_R(iw, jw, irpt, a))
      END DO
    END IF
    IF (PRESENT(curl_dX_k)) THEN
      DO c = 1, 3
        a = MOD(c, 3) + 1
        b = MOD(a, 3) + 1
        DO d = 1, 3
          curl_dX_k(iw, jw, c, d) = curl_dX_k(iw, jw, c, d) &
                                    - R_cart(d)*fac*(R_cart(a)*X_R(iw, jw, irpt, b) &
                                                     - R_cart(b)*X_R(iw, jw, irpt, a))
        END DO
      END DO
    END IF
  END SUBROUTINE fft_R2k_vec_body
  ! ================================================== !
  SUBROUTINE fft_R2k_vec_simple(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k)
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(Nw, Nw, R_vec%nRpt, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: X_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_k(Nw, Nw, 3, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: d2X_k(Nw, Nw, 3, 3, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: curl_X_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: curl_dX_k(Nw, Nw, 3, 3)
    INTEGER::iw, jw, irpt
    !
    IF (PRESENT(X_k)) X_k(:, :, :) = cmplx_0
    IF (PRESENT(dX_k)) dX_k(:, :, :, :) = cmplx_0
    IF (PRESENT(d2X_k)) d2X_k(:, :, :, :, :) = cmplx_0
    IF (PRESENT(curl_X_k)) curl_X_k(:, :, :) = cmplx_0
    IF (PRESENT(curl_dX_k)) curl_dX_k(:, :, :, :) = cmplx_0
    DO irpt = 1, R_vec%nRpt
      DO jw = 1, Nw
        DO iw = 1, Nw
          CALL fft_R2k_vec_body(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k, &
                                iw, jw, irpt, R_vec%R_cart(:, irpt))
        END DO
      END DO
    END DO
  END SUBROUTINE fft_R2k_vec_simple
  ! ================================================== !
  SUBROUTINE fft_R2k_vec_TB(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k)
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(Nw, Nw, R_vec%nRpt, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: X_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_k(Nw, Nw, 3, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: d2X_k(Nw, Nw, 3, 3, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: curl_X_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: curl_dX_k(Nw, Nw, 3, 3)
    INTEGER::iw, jw, irpt
    REAL(DP)::Rvec(3)
    !
    IF (PRESENT(X_k)) X_k(:, :, :) = cmplx_0
    IF (PRESENT(dX_k)) dX_k(:, :, :, :) = cmplx_0
    IF (PRESENT(d2X_k)) d2X_k(:, :, :, :, :) = cmplx_0
    IF (PRESENT(curl_X_k)) curl_X_k(:, :, :) = cmplx_0
    IF (PRESENT(curl_dX_k)) curl_dX_k(:, :, :, :) = cmplx_0

    DO irpt = 1, R_vec%nRpt
      DO jw = 1, Nw
        DO iw = 1, Nw
          Rvec = R_vec%R_cart(:, irpt) + shift_cart(:, iw, jw)
          CALL fft_R2k_vec_body(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k, &
                                iw, jw, irpt, Rvec)
        END DO
      END DO
    END DO
  END SUBROUTINE fft_R2k_vec_TB
  ! ================================================== !
  SUBROUTINE fft_R2k_vec_Test(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k)
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(Nw, Nw, R_vec%nRpt, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: X_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_k(Nw, Nw, 3, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: d2X_k(Nw, Nw, 3, 3, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: curl_X_k(Nw, Nw, 3)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: curl_dX_k(Nw, Nw, 3, 3)
    INTEGER::iw, jw, irpt
    REAL(DP)::R_cart(3)
    !
    IF (PRESENT(X_k)) X_k(:, :, :) = cmplx_0
    IF (PRESENT(dX_k)) dX_k(:, :, :, :) = cmplx_0
    IF (PRESENT(d2X_k)) d2X_k(:, :, :, :, :) = cmplx_0
    IF (PRESENT(curl_X_k)) curl_X_k(:, :, :) = cmplx_0
    IF (PRESENT(curl_dX_k)) curl_dX_k(:, :, :, :) = cmplx_0
    DO irpt = 1, R_vec%nRpt
      DO jw = 1, Nw
        DO iw = 1, Nw
          R_cart = R_vec%R_cart(:, irpt) + R_vec%shift_cart(:, iw, jw)
          CALL fft_R2k_vec_body(R_vec, X_R, X_k, dX_k, d2X_k, curl_X_k, curl_dX_k, &
                                iw, jw, irpt, R_cart)
        END DO
      END DO
    END DO
  END SUBROUTINE fft_R2k_vec_Test
  ! ================================================== !
END MODULE fft_base
