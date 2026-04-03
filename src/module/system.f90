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
  REAL(DP)::V_cell_3D
  !< Volume of unit cell

  !... NOT PROPERLY IMPLEMENTED YET
  ! REAL(DP)::V_cell_nD
  !< Volume of unit cell for nD system (n=1, 2, 3)
  ! INTEGER::dim = 3
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
    USE lin_mat, ONLY: inv3x3
    USE lin_vec, ONLY: dot3, cross3
    REAL(DP)::tmp_vec(3)
    bInit = .TRUE.
    real_lattice_inv = inv3x3(real_lattice)
    recip_lattice_inv = inv3x3(recip_lattice)

    tmp_vec = cross3(real_lattice(:, 1), real_lattice(:, 2))
    V_cell_3D = ABS(dot3(tmp_vec, real_lattice(:, 3)))
    ! SELECT CASE (dim)
    ! CASE (3)
    !   V_cell_nD = V_cell_3D
    ! CASE (2)
    !   ! Assume align in the xy-plane
    !   V_cell_nD = NORM2(tmp_vec)
    ! CASE (1)
    !   ! Assume align in the z-axis
    !   V_cell_nD = NORM2(real_lattice(:, 3))
    ! END SELECT
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
