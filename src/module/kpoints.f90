MODULE kpoints
  USE kinds, ONLY: DP
  IMPLICIT NONE
  TYPE::kpoint_type
    INTEGER::nkpt
    !< Number of k-points in this node
    INTEGER::nktot
    !< Total number of k-points
    REAL(DP)::wk
    !< Weight of the k-point (1/nkpt)
    REAL(DP), ALLOCATABLE:: k_cart(:, :)
    !< k-points in Cartesian coordinates (3, nkpt)
    REAL(DP), ALLOCATABLE:: k_red(:, :)
    !< k-points in reduced coordinates (3, nkpt)
  CONTAINS
    PROCEDURE::build_path => build_kpath
    PROCEDURE::build_mesh => build_kmesh
  END TYPE kpoint_type
CONTAINS
  SUBROUTINE build_kpath(self, npath, skp, nkpps)
    USE cell, ONLY: red2cart
    CLASS(kpoint_type), INTENT(INOUT)::self
    INTEGER, INTENT(IN)::npath
    REAL(DP), INTENT(IN)::skp(3, npath)
    !< Edges for each path segment (3, npath)
    INTEGER, INTENT(IN)::nkpps(npath)
    !< Number of k-points for each path segment (npath)
    INTEGER::ikpt, ipath, ikpp, nkpp
    INTEGER::nkpt

    nkpt = SUM(nkpps(1:npath - 1)) + 1
    !< Total number of k-points along the path (including the last point)
    ALLOCATE (self%k_cart(3, nkpt))
    ALLOCATE (self%k_red(3, nkpt))

    ikpt = 1
    DO ipath = 1, npath - 1
      nkpp = nkpps(ipath)
      DO ikpp = 1, nkpp
        self%k_red(:, ikpt) = skp(1:3, ipath) + &
                              (skp(1:3, ipath + 1) - skp(1:3, ipath))*REAL(ikpp - 1, DP)/REAL(nkpp, DP)
        ikpt = ikpt + 1
      END DO
    END DO
    self%k_red(:, ikpt) = skp(:, npath)
    CALL red2cart(self%k_red, self%k_cart, nkpt)
    self%nkpt = nkpt
    self%nktot = nkpt
    self%wk = 1.0_DP/REAL(nkpt, DP)
  END SUBROUTINE build_kpath
  !
  SUBROUTINE build_kmesh(self, nk1, nk2, nk3)
    USE cell, ONLY: red2cart
    CLASS(kpoint_type), INTENT(INOUT)::self
    INTEGER, INTENT(IN)::nk1, nk2, nk3
    REAL(DP)::kx, ky, kz
    INTEGER::ik1, ik2, ik3, ikpt
    INTEGER::nktot

    nktot = nk1*nk2*nk3
    ALLOCATE (self%k_cart(3, nktot))
    ALLOCATE (self%k_red(3, nktot))
    ikpt = 0
    DO ik1 = 1, nk1
      kx = REAL(ik1 - 1, DP)/REAL(nk1 - 1, DP)
      IF (2*ik1 > nk1 + 1) kx = kx - 1.0_DP
      DO ik2 = 1, nk2
        ky = REAL(ik2 - 1, DP)/REAL(nk2 - 1, DP)
        IF (2*ik2 > nk2 + 1) ky = ky - 1.0_DP
        DO ik3 = 1, nk3
          ikpt = ikpt + 1
          kz = REAL(ik3 - 1, DP)/REAL(nk3 - 1, DP)
          IF (2*ik3 > nk3 + 1) kz = kz - 1.0_DP
          self%k_red(:, ikpt) = (/kx, ky, kz/)
        END DO
      END DO
    END DO
    CALL red2cart(self%k_red, self%k_cart, nktot)
    self%nkpt = nktot
    self%nktot = nktot
    self%wk = 1.0_DP/REAL(nktot, DP)
  END SUBROUTINE build_kmesh
END MODULE kpoints
