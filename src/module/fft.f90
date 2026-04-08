MODULE fft_base
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, tpi
  USE io_global, ONLY: stdout
  USE system, ONLY: Nw
  USE wannier90, ONLY: w90data_type
  USE R_vector, ONLY: R_vec_type
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fft_q2R(w90data, R_vec, X_q, X_R, dX_R, shift)
    TYPE(w90data_type), INTENT(IN) :: w90data
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_q(..)
    !< (ldX, Nw, Nw, Nkpt)
    COMPLEX(DP), INTENT(OUT) :: X_R(..)
    !< (ldX, Nw, Nw, nRpt)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_R(..)
    !< (3, ldX, Nw, Nw, nRpt)
    REAL(DP), OPTIONAL, INTENT(IN) :: shift(3, Nw, Nw)
    INTEGER::ldX, ldY
    ldX = SIZE(X_q)/Nw/Nw/w90data%kpts%nkpt
    ldY = SIZE(X_R)/Nw/Nw/R_vec%nRpt
    IF (ldX /= ldY) THEN
      CALL errore(1, "fft_q2R", "invalid size")
    END IF
    CALL fft_q2R_4d(w90data, R_vec, ldX, X_q, X_R, dX_R, shift)
  END SUBROUTINE fft_q2R
  !
  SUBROUTINE fft_R2k(R_vec, X_R, X_k, shift, AA)
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(..)
    !< (ldX, Nw, Nw, nRpt)
    COMPLEX(DP), INTENT(OUT) :: X_k(..)
    !< (ldX, Nw, Nw)
    REAL(DP), INTENT(IN) :: shift(3, Nw, Nw)
    LOGICAL::AA
    INTEGER::ldX, ldY
    ldY = SIZE(X_R)/Nw/Nw/R_vec%nRpt
    ldX = SIZE(X_k)/Nw/Nw
    IF (ldX /= ldY) THEN
      CALL errore(1, "fft_R2k", "invalid size")
    END IF
    CALL fft_R2k_4d(R_vec, ldX, X_R, X_k, shift, AA)
  END SUBROUTINE fft_R2k
END MODULE fft_base

SUBROUTINE fft_q2R_4d(w90data, R_vec, ldX, X_q, X_R, dX_R, shift)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, tpi, zi
  USE io_global, ONLY: stdout
  USE system, ONLY: Nw
  USE wannier90, ONLY: w90data_type
  USE R_vector, ONLY: R_vec_type
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  TYPE(w90data_type), INTENT(IN) :: w90data
  TYPE(R_vec_type), INTENT(IN) :: R_vec
  INTEGER, INTENT(IN)::ldX
  COMPLEX(DP), INTENT(IN) :: X_q(ldX, Nw, Nw, w90data%kpts%nkpt)
  COMPLEX(DP), INTENT(OUT) :: X_R(ldX, Nw, Nw, R_vec%nRpt)
  COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_R(3, ldX, Nw, Nw, R_vec%nRpt)
  REAL(DP), OPTIONAL, INTENT(IN) :: shift(3, Nw, Nw)
  INTEGER::idX, iw, jw, irpt, ikpt
  REAL(DP)::phase, Rvec(3)
  COMPLEX(DP)::exp_phase, fac
  !
  CALL start_clock('fft_q2R')
  DO irpt = 1, R_vec%nRpt
    X_R(:, :, :, irpt) = zero
    IF (PRESENT(dX_R)) THEN
      dX_R(:, :, :, :, irpt) = zero
    END IF
    DO ikpt = 1, w90data%kpts%nkpt
      phase = tpi*DOT_PRODUCT(w90data%kpts%k_red(:, ikpt), R_vec%R_red(:, irpt))
      exp_phase = CMPLX(COS(phase), -SIN(phase), KIND=DP)
      DO jw = 1, Nw
        DO iw = 1, Nw
          fac = exp_phase*w90data%kpts%wk

          DO idX = 1, ldX
            X_R(idX, iw, jw, irpt) = X_R(idX, iw, jw, irpt) + X_q(idX, iw, jw, ikpt)*fac
            IF (PRESENT(dX_R)) THEN
              Rvec(:) = R_vec%R_cart(:, irpt) + shift(:, iw, jw)
              dX_R(1:3, idX, iw, jw, irpt) = dX_R(1:3, idX, iw, jw, irpt) &
                                             + zi*Rvec(1:3)*X_q(idX, iw, jw, ikpt)*fac
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
  CALL stop_clock('fft_q2R')
END SUBROUTINE fft_q2R_4d

SUBROUTINE fft_R2k_4d(R_vec, ldX, X_R, X_k, shift, AA)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, tpi
  USE io_global, ONLY: stdout
  USE system, ONLY: Nw, cart2red_real
  USE R_vector, ONLY: R_vec_type
  USE kpoints, ONLY: t_kpt, t_iks
  IMPLICIT NONE
  TYPE(R_vec_type), INTENT(IN) :: R_vec
  INTEGER, INTENT(IN)::ldX
  COMPLEX(DP), INTENT(IN) :: X_R(ldX, Nw, Nw, R_vec%nRpt)
  COMPLEX(DP), INTENT(OUT) :: X_k(ldX, Nw, Nw)
  REAL(DP), INTENT(IN) :: shift(3, Nw, Nw)
  LOGICAL, INTENT(IN) :: AA
  INTEGER::iw, jw, irpt, iuw
  REAL(DP)::phase, shift_red(3)
  COMPLEX(DP)::exp_phase, fac
  !
  CALL start_clock('fft_R2k')
  ! WRITE (stdout, '(2X, A)') '- Performing Fourier transform from R to k space...'
  X_k(:, :, :) = zero
  DO irpt = 1, R_vec%nRpt
    DO jw = 1, Nw
      DO iw = 1, Nw
        IF (AA .AND. ALL(R_vec%R_red(:, irpt) == 0) .AND. iw == jw) THEN
          CYCLE
        END IF
        iuw = R_vec%shift_map_inv(iw, jw)
        CALL cart2red_real(shift(:, iw, jw), shift_red)
        phase = tpi*DOT_PRODUCT(t_kpt%k_red(:, t_iks), R_vec%R_red(:, irpt) + shift_red)
        exp_phase = CMPLX(COS(phase), SIN(phase), KIND=DP)

        fac = exp_phase*R_vec%w_R(iw, jw, irpt)
        X_k(:, iw, jw) = X_k(:, iw, jw) + X_R(:, iw, jw, irpt)*fac
      END DO
    END DO
  END DO
  CALL stop_clock('fft_R2k')
END SUBROUTINE fft_R2k_4d
