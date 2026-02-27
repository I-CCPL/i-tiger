MODULE cell
  USE kinds, ONLY: DP
  IMPLICIT NONE
  REAL(DP) :: real_lattice(3, 3)
  !< Real lattice
  REAL(DP) :: real_lattice_inv(3, 3)
  !< Inverse of real lattice
  REAL(DP) :: recip_lattice(3, 3)
  !< Reciprocal lattice
  REAL(DP) :: recip_lattice_inv(3, 3)
  !< Inverse reciprocal lattice
  !
  INTERFACE red2cart
    MODULE PROCEDURE red2cart_1D, red2cart_2D
  END INTERFACE red2cart
  PRIVATE::red2cart_1D, red2cart_2D
  INTERFACE cart2red
    MODULE PROCEDURE cart2red_1D, cart2red_2D
  END INTERFACE cart2red
  PRIVATE::cart2red_1D, cart2red_2D
CONTAINS
  SUBROUTINE cell_setup()
    USE lin_mat3x3, ONLY: inv3x3
    real_lattice_inv = inv3x3(real_lattice)
    recip_lattice_inv = inv3x3(recip_lattice)
  END SUBROUTINE cell_setup
  !
  SUBROUTINE red2cart_1D(A_red, A_cart)
    !< Convert reduced coordinates to cartesian coordinates
    REAL(DP), INTENT(IN) :: A_red(3)
    REAL(DP), INTENT(OUT) :: A_cart(3)
    A_cart = MATMUL(real_lattice, A_red)
  END SUBROUTINE red2cart_1D
  SUBROUTINE red2cart_2D(A_red, A_cart, ndim)
    !< Convert reduced coordinates to cartesian coordinates
    REAL(DP), INTENT(IN) :: A_red(3, ndim)
    REAL(DP), INTENT(OUT) :: A_cart(3, ndim)
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: idim
    DO idim = 1, ndim
      CALL red2cart_1D(A_red(:, idim), A_cart(:, idim))
    END DO
  END SUBROUTINE red2cart_2D
  !
  SUBROUTINE cart2red_1D(A_cart, A_red)
    !< Convert cartesian coordinates to reduced coordinates
    REAL(DP), INTENT(IN) :: A_cart(3)
    REAL(DP), INTENT(OUT) :: A_red(3)
    A_red = MATMUL(real_lattice_inv, A_cart)
  END SUBROUTINE cart2red_1D
  SUBROUTINE cart2red_2D(A_cart, A_red, ndim)
    !< Convert cartesian coordinates to reduced coordinates
    REAL(DP), INTENT(IN) :: A_cart(3, ndim)
    REAL(DP), INTENT(OUT) :: A_red(3, ndim)
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: idim
    DO idim = 1, ndim
      CALL cart2red_1D(A_cart(:, idim), A_red(:, idim))
    END DO
  END SUBROUTINE cart2red_2D
END MODULE cell
