SUBROUTINE grid_idx2xyz(grid, idx, x, y, z, shift)
  !< idx starts from 1 \
  !< x,y,z start from shift
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: grid(3)
  INTEGER, INTENT(IN) :: idx
  INTEGER, INTENT(OUT) :: x, y, z
  INTEGER, INTENT(IN) :: shift(3)
  INTEGER::idx0, xy
  idx0 = idx - 1
  z = MOD(idx0, grid(1)) + shift(1)
  xy = idx0/grid(1)
  y = MOD(xy, grid(2)) + shift(2)
  x = xy/grid(2) + shift(3)
END SUBROUTINE grid_idx2xyz
