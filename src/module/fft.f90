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
  SUBROUTINE fft_q2R(w90data, R_vec, X_q, X_R)
    TYPE(w90data_type), INTENT(IN) :: w90data
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_q(..)
    !< (ldX, Nw, Nw, Nkpt)
    COMPLEX(DP), INTENT(OUT) :: X_R(..)
    !< (ldX, Nw, Nw, nRpt)
    INTEGER::ldX, ldY
    ldX = SIZE(X_q)/Nw/Nw/w90data%kpts%nkpt
    ldY = SIZE(X_R)/Nw/Nw/R_vec%nRpt
    IF (ldX /= ldY) THEN
      CALL errore(1, "fft_q2R", "invalid size")
    END IF
    CALL fft_q2R_4d(w90data, R_vec, ldX, X_q, X_R)
  END SUBROUTINE fft_q2R
  !
  SUBROUTINE fft_R2k(R_vec, X_R, X_k)
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(..)
    !< (ldX, Nw, Nw, nRpt)
    COMPLEX(DP), INTENT(OUT) :: X_k(..)
    !< (ldX, Nw, Nw, Nkpt)
    INTEGER::ldX, ldY
    ldY = SIZE(X_R)/Nw/Nw/R_vec%nRpt
    ldX = SIZE(X_k)/Nw/Nw
    IF (ldX /= ldY) THEN
      CALL errore(1, "fft_R2k", "invalid size")
    END IF
    CALL fft_R2k_4d(R_vec, ldX, X_R, X_k)
  END SUBROUTINE fft_R2k
END MODULE fft_base

SUBROUTINE fft_q2R_4d(w90data, R_vec, ldX, X_q, X_R)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, tpi
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
  INTEGER::iw, jw, irpt, ikpt
  REAL(DP)::phase
  COMPLEX(DP)::exp_phase
  !
  ! WRITE (stdout, '(2X, A)') '- Performing Fourier transform from q to R space...'
  DO irpt = 1, R_vec%nRpt
    X_R(:, :, :, irpt) = zero
    DO ikpt = 1, w90data%kpts%nkpt
      phase = -tpi*DOT_PRODUCT(w90data%kpts%k_red(:, ikpt), R_vec%R_red(:, irpt))
      exp_phase = CMPLX(COS(phase), SIN(phase), KIND=DP)
      DO jw = 1, Nw
        DO iw = 1, Nw
          X_R(:, iw, jw, irpt) = X_R(:, iw, jw, irpt) &
                                 + X_q(:, iw, jw, ikpt)*exp_phase*w90data%kpts%wk*R_vec%w_R(iw, jw, irpt)
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE fft_q2R_4d

SUBROUTINE fft_R2k_4d(R_vec, ldX, X_R, X_k)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, tpi
  USE io_global, ONLY: stdout
  USE system, ONLY: Nw
  USE R_vector, ONLY: R_vec_type
  USE kpoints, ONLY: t_kpt, t_iks
  IMPLICIT NONE
  TYPE(R_vec_type), INTENT(IN) :: R_vec
  INTEGER, INTENT(IN)::ldX
  COMPLEX(DP), INTENT(IN) :: X_R(ldX, Nw, Nw, R_vec%nRpt)
  COMPLEX(DP), INTENT(OUT) :: X_k(ldX, Nw, Nw)
  INTEGER::iw, jw, irpt
  REAL(DP)::phase
  COMPLEX(DP)::exp_phase
  !
  ! WRITE (stdout, '(2X, A)') '- Performing Fourier transform from R to k space...'
  X_k(:, :, :) = zero
  DO irpt = 1, R_vec%nRpt
    phase = tpi*DOT_PRODUCT(t_kpt%k_red(:, t_iks), R_vec%R_red(:, irpt))
    exp_phase = CMPLX(COS(phase), SIN(phase), KIND=DP)
    DO jw = 1, Nw
      DO iw = 1, Nw
        X_k(:, iw, jw) = X_k(:, iw, jw) + X_R(:, iw, jw, irpt)*exp_phase
      END DO
    END DO
  END DO
END SUBROUTINE fft_R2k_4d
