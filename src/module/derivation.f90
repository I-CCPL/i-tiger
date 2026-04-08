MODULE der_base
  USE kinds, ONLY: DP
  USE system, ONLY: Nw
  IMPLICIT NONE
CONTAINS
  SUBROUTINE der_q(w90data, X_k, dX_k)
    USE wannier90, ONLY: w90data_type
    TYPE(w90data_type), INTENT(IN)::w90data
    COMPLEX(DP), INTENT(IN)::X_k(..)
    !< (ldX, Nw, Nw, Nkpt)
    COMPLEX(DP), INTENT(OUT)::dX_k(..)
    !< (3, ldX, Nw, Nw, Nkpt)
    INTEGER::ldX
    ldX = SIZE(X_k)
    IF (ldX*3 /= SIZE(dX_k)) THEN
      CALL errore(1, 'der_q', 'invalid size')
    END IF
    ldX = ldX/Nw/NW/w90data%kpts%nkpt
    CALL derivation_q_4D(w90data, ldX, X_k, dX_k)
  END SUBROUTINE der_q
  !
  SUBROUTINE der_R(R_vec, X_R, dX_R, shift_cart)
    USE R_vector, ONLY: R_vec_type
    TYPE(R_vec_type), INTENT(IN)::R_vec
    COMPLEX(DP), INTENT(IN)::X_R(..)
    !< (ldX, Nw, Nw, nRpt)
    COMPLEX(DP), INTENT(OUT)::dX_R(..)
    !< (3, ldX, Nw, Nw, NRpt)
    REAL(DP), INTENT(IN) :: shift_cart(3, Nw, Nw)
    INTEGER::ldX
    ldX = SIZE(X_R)
    IF (ldX*3 /= SIZE(dX_R)) THEN
      CALL errore(1, 'der_R', 'invalid size')
    END IF
    ldX = ldX/Nw/Nw/R_vec%nRpt
    CALL derivation_R_4D(R_vec, ldX, X_R, dX_R, shift_cart)
  END SUBROUTINE der_R
END MODULE der_base

SUBROUTINE derivation_q_4D(w90data, ldX, X_k, dX_k)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero
  USE wannier90, ONLY: w90data_type
  USE system, ONLY: Nw, red2cart_recip
  IMPLICIT NONE
  TYPE(w90data_type), INTENT(IN)::w90data
  INTEGER, INTENT(IN)::ldX
  COMPLEX(DP), INTENT(IN)::X_k(ldX, Nw, Nw, w90data%kpts%nkpt)
  COMPLEX(DP), INTENT(OUT)::dX_k(3, ldX, Nw, Nw, w90data%kpts%nkpt)
  COMPLEX(DP):: dX_bk(3)
  REAL(DP)::bvec_cart(3)
  INTEGER::inb, jnb, ikpt, jkpt, iw, jw, idim

  DO ikpt = 1, w90data%kpts%nkpt
    DO jw = 1, Nw
      DO iw = 1, Nw
        DO idim = 1, ldX
          dX_bk(:) = zero
          DO inb = 1, w90data%nnb
            jnb = w90data%bvec_index(inb, ikpt)
            jkpt = w90data%neighbour_k(inb, ikpt)
            CALL red2cart_recip(w90data%bvec_red(:, jnb), bvec_cart)
            dX_bk(:) = dX_bk(:) + X_k(idim, iw, jw, jkpt)*w90data%wb(jnb)*bvec_cart(:)
          END DO
          dX_k(:, idim, iw, jw, ikpt) = dX_bk(:)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE derivation_q_4D

SUBROUTINE derivation_R_4D(R_vec, ldX, X_R, dX_R, shift_cart)
  USE kinds, ONLY: DP
  USE constants, ONLY: zero, zi
  USE system, ONLY: Nw
  USE R_vector, ONLY: R_vec_type
  IMPLICIT NONE
  TYPE(R_vec_type), INTENT(IN)::R_vec
  INTEGER, INTENT(IN)::ldX
  COMPLEX(DP), INTENT(IN)::X_R(ldX, Nw, Nw, R_vec%nRpt)
  COMPLEX(DP), INTENT(OUT)::dX_R(3, ldX, Nw, Nw, R_vec%nRpt)
  REAL(DP), INTENT(IN) :: shift_cart(3, Nw, Nw)
  INTEGER::iRpt, iw, jw, iuw, idx
  REAL(DP)::Rvec(3)
  !
  DO iRpt = 1, R_vec%nRpt
    DO jw = 1, Nw
      DO iw = 1, Nw
        Rvec = shift_cart(:, iw, jw) + R_vec%R_cart(:, iRpt)
        DO idx = 1, ldX
          dX_R(:, idx, iw, jw, iRpt) = zi*X_R(idx, iw, jw, iRpt)*Rvec(:)
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE derivation_R_4D
