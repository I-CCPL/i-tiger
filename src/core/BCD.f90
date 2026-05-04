!... Berry curvature dipole (BCD) calculations
SUBROUTINE compute_BCD_sea()
  !< Computes the BCD contribution from the Fermi-sea (occupied) states.
  IMPLICIT NONE
  !
  CALL errore(1, 'compute_BCD_sea', 'Not implemented yet')
END SUBROUTINE compute_BCD_sea

SUBROUTINE compute_BCD_surf()
  !< Computes the BCD contribution from the Fermi-surface (boundary) states.
  IMPLICIT NONE
  CALL errore(1, 'compute_BCD_surf', 'Not implemented yet')
END SUBROUTINE compute_BCD_surf
