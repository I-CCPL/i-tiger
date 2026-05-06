MODULE wannier90
  USE kinds, ONLY: DP
  USE io_global, ONLY: ionode, stdout, get_free_unit
  USE kpoints, ONLY: kpoint_type
  IMPLICIT NONE
  PRIVATE
  LOGICAL, PUBLIC::chk_w90 = .FALSE.
  LOGICAL, PUBLIC::Hq_band = .FALSE.
  PUBLIC::chk_dum_type, w90data_type, w90data

  ! ==================================================
  ! chk_dum
  ! ==================================================
  TYPE::chk_dum_type
    CHARACTER(LEN=33) :: header
    INTEGER :: nbnd_excl
    INTEGER, ALLOCATABLE :: excl_bands(:)
    CHARACTER(LEN=20) :: checkpoint
    LOGICAL :: have_disentangled
    REAL(DP) :: omega_invariant
    LOGICAL, ALLOCATABLE :: lwindow(:, :)
    COMPLEX(DP), ALLOCATABLE :: u_matrix_opt(:, :, :)
    COMPLEX(DP), ALLOCATABLE::u_matrix(:, :, :)
    COMPLEX(DP), ALLOCATABLE::m_matrix(:, :, :, :)
  CONTAINS
    PROCEDURE::clear => clear_chk_dum
  END TYPE chk_dum_type
  INTERFACE
    MODULE SUBROUTINE clear_chk_dum(self)
      CLASS(chk_dum_type), INTENT(INOUT) :: self
    END SUBROUTINE clear_chk_dum
  END INTERFACE

  ! ==================================================
  ! w90data
  ! ==================================================
  TYPE::w90data_type
    !... General data
    CHARACTER(LEN=256) :: prefix

    INTEGER::nbnd
    !< Number of bands. \
    !< nbnd==Nw for non-disentangled case \
    !< nbnd>=Nw for disentangled case.
    INTEGER :: nnb
    !< Number of nearest neighbor k-points
    INTEGER :: k_grid(3)
    !< Number of k-points along each reciprocal lattice vector \
    !< (Monkhorst-Pack grid)
    TYPE(kpoint_type)::kpts
    REAL(DP), ALLOCATABLE :: wannier_center_cart(:, :)
    !< Wannier center in Cartesian coordinates (3, Nw)
    REAL(DP), ALLOCATABLE :: wannier_spread(:)
    !< (Nw)

    !... chk data
    COMPLEX(DP), ALLOCATABLE :: v_matrix(:, :, :)
    !< (nbnd, Nw, nkpt) for disentangled case
    INTEGER, ALLOCATABLE :: ndimwin(:)
    !< (nkpt) Number of bands in the outer window
    INTEGER, ALLOCATABLE::win_min(:)
    !< (nkpt) Minimum band index of the outer window
    REAL(DP), ALLOCATABLE::eigval(:, :)
    !< Eigenvalues (nbnd, nkpt)
    COMPLEX(DP), ALLOCATABLE :: eigvec(:, :, :)
    !< Unitary matrix corresponding eigenvector \
    !< (Nw, Nw, nkpt)
    COMPLEX(DP), ALLOCATABLE :: Hq(:, :, :)
    !< Hamiltonian in Wannier gauge (Nw, Nw, nkpt)

    !... mmn data
    INTEGER, ALLOCATABLE :: neighbour_k(:, :)
    !< (nnb, nkpt)
    INTEGER, ALLOCATABLE :: neighbour_g(:, :, :)
    !< (3, nnb, nkpt)
    COMPLEX(DP), ALLOCATABLE :: overlap(:, :, :, :)
    !< (nbnd, nbnd, nnb, nkpt)
    INTEGER, ALLOCATABLE::bvec_index(:, :)
    !< index map from (nnb, nkpt) to (nnb)
    REAL(DP), ALLOCATABLE :: bvec_red(:, :)
    !< b vectors in reduced coordinates (3, nnb)
    REAL(DP), ALLOCATABLE::wb(:)
    !< weight of b vector (nnb)
    COMPLEX(DP), ALLOCATABLE :: Aq(:, :, :, :)
    !< (Nw, Nw, nkpt, 3)
  CONTAINS
    PROCEDURE::clear => clear_w90_data
    PROCEDURE::read_files => read_w90_files
    PROCEDURE::read_chk => read_w90_chk
    PROCEDURE::build_Hq => build_w90_Hq
    PROCEDURE::bcast_Hq => bcast_w90_Hq
    PROCEDURE::read_eig => read_w90_eig

    PROCEDURE::read_mmn => read_w90_mmn
    PROCEDURE::build_bvec => build_w90_bvec
    PROCEDURE::bcast_bvec => bcast_w90_bvec
    PROCEDURE::build_Aq => build_w90_Aq
    PROCEDURE::bcast_Aq => bcast_w90_Aq
  END TYPE w90data_type
  TYPE(w90data_type)::w90data

  INTERFACE
    MODULE SUBROUTINE clear_w90_data(self)
      CLASS(w90data_type), INTENT(INOUT)::self
    END SUBROUTINE clear_w90_data

    MODULE SUBROUTINE read_w90_files(self)
      CLASS(w90data_type), INTENT(INOUT) :: self
    END SUBROUTINE read_w90_files
    ! ----------------------------------------
    MODULE SUBROUTINE read_w90_chk(self, chk_dum)
      CLASS(w90data_type), INTENT(INOUT)::self
      TYPE(chk_dum_type), INTENT(OUT)::chk_dum
    END SUBROUTINE read_w90_chk

    MODULE SUBROUTINE build_w90_Hq(self)
      CLASS(w90data_type), INTENT(INOUT) :: self
    END SUBROUTINE build_w90_Hq

    MODULE SUBROUTINE bcast_w90_Hq(self)
      CLASS(w90data_type), INTENT(INOUT) :: self
    END SUBROUTINE bcast_w90_Hq

    MODULE SUBROUTINE read_w90_eig(self)
      CLASS(w90data_type), INTENT(INOUT) :: self
    END SUBROUTINE read_w90_eig
    ! ----------------------------------------
    MODULE SUBROUTINE read_w90_mmn(self)
      CLASS(w90data_type), INTENT(INOUT) :: self
    END SUBROUTINE read_w90_mmn

    MODULE SUBROUTINE build_w90_bvec(self)
      CLASS(w90data_type), INTENT(inout) :: self
    END SUBROUTINE build_w90_bvec
    MODULE SUBROUTINE bcast_w90_bvec(self)
      CLASS(w90data_type), INTENT(inout) :: self
    END SUBROUTINE bcast_w90_bvec

    MODULE SUBROUTINE build_w90_Aq(self)
      CLASS(w90data_type), INTENT(INOUT) :: self
    END SUBROUTINE build_w90_Aq
    MODULE SUBROUTINE bcast_w90_Aq(self)
      CLASS(w90data_type), INTENT(INOUT) :: self
    END SUBROUTINE bcast_w90_Aq
  END INTERFACE

  ! ==================================================
  ! General w90
  ! ==================================================
  INTERFACE

    MODULE FUNCTION wannier_gauge_diag(mat_H, v1, v2) RESULT(retval)
      REAL(DP), INTENT(IN) :: mat_H(:)
      COMPLEX(DP), INTENT(IN) :: v1(:), v2(:)
      COMPLEX(DP):: retval
    END FUNCTION wannier_gauge_diag

    MODULE FUNCTION wannier_gauge(mat_H, v1, v2) RESULT(retval)
      COMPLEX(DP), INTENT(IN) :: mat_H(:, :)
      COMPLEX(DP), INTENT(IN) :: v1(:), v2(:)
      COMPLEX(DP) :: retval
    END FUNCTION wannier_gauge

    MODULE SUBROUTINE write_chk_dump(self, chk_dum)
      TYPE(w90data_type), INTENT(IN) :: self
      TYPE(chk_dum_type), INTENT(IN) :: chk_dum
    END SUBROUTINE write_chk_dump
  END INTERFACE
END MODULE wannier90
