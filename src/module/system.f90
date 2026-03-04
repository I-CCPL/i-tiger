MODULE system
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER::Nw
  LOGICAL, PRIVATE::bInit = .FALSE.
  !< flag for initialized.
  REAL(DP) :: real_lattice(3, 3)
  !< Real lattice (a)
  REAL(DP) :: real_lattice_inv(3, 3)
  !< Inverse of real lattice
  REAL(DP) :: recip_lattice(3, 3)
  !< Reciprocal lattice (a*b=2pi)
  REAL(DP) :: recip_lattice_inv(3, 3)
  !< Inverse reciprocal lattice
  !
  INTERFACE red2cart_real
    MODULE PROCEDURE red2cart_real_1D, red2cart_real_2D
  END INTERFACE red2cart_real
  PRIVATE::red2cart_real_1D, red2cart_real_2D
  INTERFACE cart2red_real
    MODULE PROCEDURE cart2red_real_1D, cart2red_real_2D
  END INTERFACE cart2red_real
  PRIVATE::cart2red_real_1D, cart2red_real_2D
  !
  INTERFACE red2cart_recip
    MODULE PROCEDURE red2cart_recip_1D, red2cart_recip_2D
  END INTERFACE red2cart_recip
  PRIVATE::red2cart_recip_1D, red2cart_recip_2D
  INTERFACE cart2red_recip
    MODULE PROCEDURE cart2red_recip_1D, cart2red_recip_2D
  END INTERFACE cart2red_recip
  PRIVATE::cart2red_recip_1D, cart2red_recip_2D
CONTAINS
  SUBROUTINE cell_setup()
    USE lin_mat3x3, ONLY: inv3x3
    bInit = .TRUE.
    real_lattice_inv = inv3x3(real_lattice)
    recip_lattice_inv = inv3x3(recip_lattice)
  END SUBROUTINE cell_setup
  !
  SUBROUTINE red2cart_real_1D(A_red, A_cart)
    !< Convert reduced coordinates to cartesian coordinates in real space
    REAL(DP), INTENT(IN) :: A_red(3)
    REAL(DP), INTENT(OUT) :: A_cart(3)
    IF (.NOT. bInit) CALL errore(1, 'red2cart', 'call before init')
    A_cart = MATMUL(real_lattice, A_red)
  END SUBROUTINE red2cart_real_1D
  SUBROUTINE red2cart_real_2D(A_red, A_cart, ndim)
    !< Convert reduced coordinates to cartesian coordinates in real space
    REAL(DP), INTENT(IN) :: A_red(3, ndim)
    REAL(DP), INTENT(OUT) :: A_cart(3, ndim)
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: idim
    DO idim = 1, ndim
      CALL red2cart_real_1D(A_red(:, idim), A_cart(:, idim))
    END DO
  END SUBROUTINE red2cart_real_2D
  !
  SUBROUTINE cart2red_real_1D(A_cart, A_red)
    !< Convert cartesian coordinates to reduced coordinates in real space
    REAL(DP), INTENT(IN) :: A_cart(3)
    REAL(DP), INTENT(OUT) :: A_red(3)
    IF (.NOT. bInit) CALL errore(1, 'cart2red', 'call before init')
    A_red = MATMUL(real_lattice_inv, A_cart)
  END SUBROUTINE cart2red_real_1D
  SUBROUTINE cart2red_real_2D(A_cart, A_red, ndim)
    !< Convert cartesian coordinates to reduced coordinates in real space
    REAL(DP), INTENT(IN) :: A_cart(3, ndim)
    REAL(DP), INTENT(OUT) :: A_red(3, ndim)
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: idim
    DO idim = 1, ndim
      CALL cart2red_real_1D(A_cart(:, idim), A_red(:, idim))
    END DO
  END SUBROUTINE cart2red_real_2D
  !
  SUBROUTINE red2cart_recip_1D(A_red, A_cart)
    !< Convert reduced coordinates to cartesian coordinates in reciprocal space
    REAL(DP), INTENT(IN) :: A_red(3)
    REAL(DP), INTENT(OUT) :: A_cart(3)
    IF (.NOT. bInit) CALL errore(1, 'red2cart_recip', 'call before init')
    A_cart = MATMUL(recip_lattice, A_red)
  END SUBROUTINE red2cart_recip_1D
  SUBROUTINE red2cart_recip_2D(A_red, A_cart, ndim)
    !< Convert reduced coordinates to cartesian coordinates in reciprocal space
    REAL(DP), INTENT(IN) :: A_red(3, ndim)
    REAL(DP), INTENT(OUT) :: A_cart(3, ndim)
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: idim
    DO idim = 1, ndim
      CALL red2cart_recip_1D(A_red(:, idim), A_cart(:, idim))
    END DO
  END SUBROUTINE red2cart_recip_2D
  !
  SUBROUTINE cart2red_recip_1D(A_cart, A_red)
    !< Convert cartesian coordinates to reduced coordinates in reciprocal space
    REAL(DP), INTENT(IN) :: A_cart(3)
    REAL(DP), INTENT(OUT) :: A_red(3)
    IF (.NOT. bInit) CALL errore(1, 'cart2red_recip', 'call before init')
    A_red = MATMUL(recip_lattice_inv, A_cart)
  END SUBROUTINE cart2red_recip_1D
  SUBROUTINE cart2red_recip_2D(A_cart, A_red, ndim)
    !< Convert cartesian coordinates to reduced coordinates in reciprocal space
    REAL(DP), INTENT(IN) :: A_cart(3, ndim)
    REAL(DP), INTENT(OUT) :: A_red(3, ndim)
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: idim
    DO idim = 1, ndim
      CALL cart2red_recip_1D(A_cart(:, idim), A_red(:, idim))
    END DO
  END SUBROUTINE cart2red_recip_2D
END MODULE system
