MODULE fft_base
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, tpi
  USE io_global, ONLY: stdout
  USE io_input, ONLY: convention
  USE system, ONLY: Nw
  USE wannier90, ONLY: w90data_type
  USE R_vector, ONLY: R_vec_type
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  REAL(DP), ALLOCATABLE::shift_cart(:, :, :)
  REAL(DP), ALLOCATABLE::shift_red(:, :, :)
CONTAINS
  SUBROUTINE fft_init(R_vec, A_R)
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    REAL(DP), INTENT(IN) :: A_R(3, Nw, Nw, R_vec%nRpt)
    REAL(DP)::center(3, Nw)
    INTEGER::irpt, iw, jw
    ALLOCATE (shift_cart(3, Nw, Nw))
    ALLOCATE (shift_red(3, Nw, Nw))
    DO irpt = 1, R_vec%nRpt
      IF (ALL(R_vec%R_red(:, irpt) == 0)) THEN
        DO iw = 1, Nw
          center(:, iw) = DBLE(A_R(:, iw, iw, irpt))
        END DO
        EXIT
      END IF
    END DO

    DO iw = 1, Nw
      DO jw = 1, Nw
        shift_cart(:, iw, jw) = center(:, jw) - center(:, iw)
      END DO
    END DO

  END SUBROUTINE fft_init
  !
  SUBROUTINE fft_q2R(w90data, R_vec, X_q, X_R, dX_R)
    TYPE(w90data_type), INTENT(IN) :: w90data
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_q(..)
    !< (Nw, Nw, Nkpt, ldX)
    COMPLEX(DP), INTENT(OUT) :: X_R(..)
    !< (ldX, Nw, Nw, nRpt)
    COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_R(..)
    !< (3, ldX, Nw, Nw, nRpt)
    INTEGER::ldX, ldY
    ldX = SIZE(X_q)/Nw/Nw/w90data%kpts%nkpt
    ldY = SIZE(X_R)/Nw/Nw/R_vec%nRpt
    IF (ldX /= ldY) THEN
      CALL errore(1, "fft_q2R", "invalid size")
    END IF

    CALL start_clock('fft_q2R')
    SELECT CASE (convention)
    CASE (0)
      CALL fft_q2R_4d_0(w90data, R_vec, ldX, X_q, X_R, dX_R)
    CASE (1)
      CALL fft_q2R_4d_1(w90data, R_vec, ldX, X_q, X_R, dX_R, shift_cart)
    CASE (2)
      CALL fft_q2R_4d_2(w90data, R_vec, ldX, X_q, X_R)
    CASE default
      CALL errore(1, 'fft_q2R', 'invalid convention')
    END SELECT
    CALL stop_clock('fft_q2R')
  END SUBROUTINE fft_q2R
  !
  SUBROUTINE fft_R2k(R_vec, X_R, X_k, AA)
    TYPE(R_vec_type), INTENT(IN) :: R_vec
    COMPLEX(DP), INTENT(IN) :: X_R(..)
    !< (ldX, Nw, Nw, nRpt)
    COMPLEX(DP), INTENT(OUT) :: X_k(..)
    !< (ldX, Nw, Nw)
    LOGICAL, INTENT(IN)::AA
    INTEGER::ldX, ldY
    ldY = SIZE(X_R)/Nw/Nw/R_vec%nRpt
    ldX = SIZE(X_k)/Nw/Nw
    IF (ldX /= ldY) THEN
      CALL errore(1, "fft_R2k", "invalid size")
    END IF

    CALL start_clock('fft_R2k')
    SELECT CASE (convention)
    CASE (0)
      CALL fft_R2k_4d_0(R_vec, ldX, X_R, X_k)
    CASE (1)
      CALL fft_R2k_4d_1(R_vec, ldX, X_R, X_k, shift_red, AA)
    CASE (2)
      CALL fft_R2k_4d_2(R_vec, ldX, X_R, X_k)
    CASE default
      CALL errore(1, 'fft_R2k', 'invalid convention')
    END SELECT
    CALL stop_clock('fft_R2k')
  END SUBROUTINE fft_R2k
END MODULE fft_base

SUBROUTINE fft_q2R_4d_0(w90data, R_vec, ldX, X_q, X_R, dX_R)
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
  COMPLEX(DP), INTENT(IN) :: X_q(Nw, Nw, w90data%kpts%nkpt, ldX)
  COMPLEX(DP), INTENT(OUT) :: X_R(ldX, Nw, Nw, R_vec%nRpt)
  COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_R(3, ldX, Nw, Nw, R_vec%nRpt)
  INTEGER::idX, iw, jw, irpt, ikpt
  REAL(DP)::phase, Rvec(3)
  COMPLEX(DP)::exp_phase, fac
  !
  DO irpt = 1, R_vec%nRpt
    X_R(:, :, :, irpt) = zero
    IF (PRESENT(dX_R)) THEN
      dX_R(:, :, :, :, irpt) = zero
    END IF
    DO ikpt = 1, w90data%kpts%nkpt
      phase = tpi*DOT_PRODUCT(w90data%kpts%k_red(:, ikpt), R_vec%R_red(:, irpt))
      exp_phase = CMPLX(COS(phase), -SIN(phase), KIND=DP)
      fac = exp_phase*w90data%kpts%wk
      DO jw = 1, Nw
        DO iw = 1, Nw

          DO idX = 1, ldX
            X_R(idX, iw, jw, irpt) = X_R(idX, iw, jw, irpt) + X_q(iw, jw, ikpt, idX)*fac
            IF (PRESENT(dX_R)) THEN
              dX_R(1:3, idX, iw, jw, irpt) = dX_R(1:3, idX, iw, jw, irpt) &
                                             + zi*R_vec%R_cart(1:3, irpt)*X_q(iw, jw, ikpt, idX)*fac
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE fft_q2R_4d_0

SUBROUTINE fft_q2R_4d_1(w90data, R_vec, ldX, X_q, X_R, dX_R, shift)
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
  COMPLEX(DP), INTENT(IN) :: X_q(Nw, Nw, w90data%kpts%nkpt, ldX)
  COMPLEX(DP), INTENT(OUT) :: X_R(ldX, Nw, Nw, R_vec%nRpt)
  COMPLEX(DP), OPTIONAL, INTENT(OUT) :: dX_R(3, ldX, Nw, Nw, R_vec%nRpt)
  REAL(DP), OPTIONAL, INTENT(IN) :: shift(3, Nw, Nw)
  INTEGER::idX, iw, jw, irpt, ikpt
  REAL(DP)::phase, Rvec(3)
  COMPLEX(DP)::exp_phase, fac
  !
  DO irpt = 1, R_vec%nRpt
    X_R(:, :, :, irpt) = zero
    IF (PRESENT(dX_R)) THEN
      dX_R(:, :, :, :, irpt) = zero
    END IF
    DO ikpt = 1, w90data%kpts%nkpt
      phase = tpi*DOT_PRODUCT(w90data%kpts%k_red(:, ikpt), R_vec%R_red(:, irpt))
      exp_phase = CMPLX(COS(phase), -SIN(phase), KIND=DP)
      fac = exp_phase*w90data%kpts%wk
      DO jw = 1, Nw
        DO iw = 1, Nw

          DO idX = 1, ldX
            X_R(idX, iw, jw, irpt) = X_R(idX, iw, jw, irpt) + X_q(iw, jw, ikpt, idX)*fac
            IF (PRESENT(dX_R)) THEN
              Rvec(:) = R_vec%R_cart(:, irpt) + shift(:, iw, jw)
              dX_R(1:3, idX, iw, jw, irpt) = dX_R(1:3, idX, iw, jw, irpt) &
                                             + zi*Rvec(1:3)*X_q(iw, jw, ikpt, idX)*fac
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE fft_q2R_4d_1

SUBROUTINE fft_R2k_4d_0(R_vec, ldX, X_R, X_k)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, tpi, zi
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
  COMPLEX(DP)::exp_phase, fac
  !
  X_k(:, :, :) = zero
  DO irpt = 1, R_vec%nRpt
    DO jw = 1, Nw
      DO iw = 1, Nw
        phase = tpi*DOT_PRODUCT(t_kpt%k_red(:, t_iks), R_vec%R_red(:, irpt))
        exp_phase = CMPLX(COS(phase), SIN(phase), KIND=DP)
        fac = exp_phase*R_vec%w_R(iw, jw, irpt)
        X_k(:, iw, jw) = X_k(:, iw, jw) + X_R(:, iw, jw, irpt)*fac
      END DO
    END DO
  END DO
END SUBROUTINE fft_R2k_4d_0

SUBROUTINE fft_R2k_4d_1(R_vec, ldX, X_R, X_k, shift_red, AA)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, tpi, zi
  USE io_global, ONLY: stdout
  USE system, ONLY: Nw, cart2red_real
  USE R_vector, ONLY: R_vec_type
  USE kpoints, ONLY: t_kpt, t_iks
  IMPLICIT NONE
  TYPE(R_vec_type), INTENT(IN) :: R_vec
  INTEGER, INTENT(IN)::ldX
  COMPLEX(DP), INTENT(IN) :: X_R(ldX, Nw, Nw, R_vec%nRpt)
  COMPLEX(DP), INTENT(OUT) :: X_k(ldX, Nw, Nw)
  REAL(DP), INTENT(IN) :: shift_red(3, Nw, Nw)
  LOGICAL, INTENT(IN) :: AA
  INTEGER::iw, jw, irpt
  REAL(DP)::phase
  LOGICAL::R_zero
  !
  X_k(:, :, :) = zero
  DO irpt = 1, R_vec%nRpt
    R_zero = AA .AND. ALL(R_vec%R_red(:, irpt) == 0)
    DO jw = 1, Nw
      DO iw = 1, Nw
        IF (R_zero .AND. iw == jw) THEN
          CYCLE
        END IF
        phase = tpi*DOT_PRODUCT(t_kpt%k_red(:, t_iks), R_vec%R_red(:, irpt) + shift_red(:, iw, jw))
        X_k(:, iw, jw) = X_k(:, iw, jw) &
                         + X_R(:, iw, jw, irpt)*EXP(zi*phase)*R_vec%w_R(iw, jw, irpt)
      END DO
    END DO
  END DO
END SUBROUTINE fft_R2k_4d_1

SUBROUTINE fft_q2R_4d_2(w90data, R_vec, ldX, X_q, X_R)
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
  COMPLEX(DP), INTENT(IN) :: X_q(Nw, Nw, w90data%kpts%nkpt, ldX)
  COMPLEX(DP), INTENT(OUT) :: X_R(ldX, Nw, Nw, R_vec%nRpt)
  INTEGER::iw, jw, irpt, ikpt
  REAL(DP)::phase
  COMPLEX(DP)::exp_phase, fac
  !
  DO irpt = 1, R_vec%nRpt
    X_R(:, :, :, irpt) = zero
    DO ikpt = 1, w90data%kpts%nkpt
      phase = -tpi*DOT_PRODUCT(w90data%kpts%k_red(:, ikpt), R_vec%R_red(:, irpt))
      exp_phase = CMPLX(COS(phase), SIN(phase), KIND=DP)
      fac = exp_phase*w90data%kpts%wk
      DO jw = 1, Nw
        DO iw = 1, Nw
          X_R(:, iw, jw, irpt) = X_R(:, iw, jw, irpt) &
                                 + X_q(iw, jw, ikpt, :)*fac*R_vec%w_R(iw, jw, irpt)
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE fft_q2R_4d_2

SUBROUTINE fft_R2k_4d_2(R_vec, ldX, X_R, X_k)
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
  INTEGER::iw, jw, irpt, iuw
  REAL(DP)::phase, shift_red(3)
  COMPLEX(DP)::exp_phase
  !
  X_k(:, :, :) = zero
  DO irpt = 1, R_vec%nRpt
    DO jw = 1, Nw
      DO iw = 1, Nw
        iuw = R_vec%shift_map_inv(iw, jw)
        CALL cart2red_real(R_vec%shift_cart(:, iuw), shift_red)
        phase = tpi*DOT_PRODUCT(t_kpt%k_red(:, t_iks), R_vec%R_red(:, irpt) + shift_red)
        exp_phase = CMPLX(COS(phase), SIN(phase), KIND=DP)

        X_k(:, iw, jw) = X_k(:, iw, jw) + X_R(:, iw, jw, irpt)*exp_phase
      END DO
    END DO
  END DO
END SUBROUTINE fft_R2k_4d_2
