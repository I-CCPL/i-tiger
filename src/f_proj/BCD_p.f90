!... Berry curvature dipole (BCD) calculations
SUBROUTINE compute_BCD_sea_p()
  !< Computes the BCD contribution from the Fermi-sea (occupied) states.
  IMPLICIT NONE
  !
  CALL errore(1, 'compute_BCD_sea_p', 'Not implemented yet')
END SUBROUTINE compute_BCD_sea_p

SUBROUTINE compute_BCD_surf_p(O_bar, A_bar, dH_bar, BCD_surf)
  !< Computes the BCD contribution from the Fermi-surface (boundary) states.
  USE kinds, ONLY: DP
  USE f_params, ONLY: Ef_min, Ef_max, Ef_step, Ef_nE, E_fermi, dE_thr, dE_eta
  USE kpoints, ONLY: t_kpt, t_iks
  USE system, ONLY: Nw
  USE delta_func, ONLY: w1gauss
  USE kpoints, ONLY: t_kpt
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN)::O_bar(Nw, Nw, 3)
  COMPLEX(DP), INTENT(IN)::A_bar(Nw, Nw, 3)
  COMPLEX(DP), INTENT(IN)::dH_bar(Nw, Nw, 3)
  REAL(DP), INTENT(INOUT)::BCD_surf(3, 3, Ef_nE)
  INTEGER::ief, j, n
  REAL(DP)::delta
  REAL(DP)::imf_k(3, 3, Ef_nE)
  REAL(DP)::curv_nk(3)
  LOGICAL::bdone
  DO n = 1, Nw
    ! IF (n > 1) THEN
    !   IF (t_kpt%eigval(n, t_iks) - t_kpt%eigval(n - 1, t_iks) <= dE_thr) CYCLE
    ! END IF
    ! IF (n < Nw) THEN
    !   IF (t_kpt%eigval(n + 1, t_iks) - t_kpt%eigval(n, t_iks) <= dE_thr) CYCLE
    ! END IF
    bdone = .FALSE.
    DO ief = 1, Ef_nE
      E_fermi = Ef_min + (ief - 1)*Ef_step
      IF (ABS(t_kpt%eigval(n, t_iks) - E_fermi) > 5*dE_eta) CYCLE

      IF (.NOT. bdone) THEN
        CALL get_berry_nk_p(n, t_kpt%eigval(:, t_iks), A_bar, dH_bar, curv_nk(:))
        bdone = .TRUE.
      END IF
      delta = w1gauss(t_kpt%eigval(n, t_iks) - E_fermi, dE_eta, 0)*t_kpt%wk
      ! dH_bar [eV*Ang], curv_nk [Ang^2], delta [1/eV] => BCD_surf [Ang^3]
      DO j = 1, 3
        BCD_surf(:, j, ief) = BCD_surf(:, j, ief) + DBLE(dH_bar(n, n, :))*curv_nk(j)*delta
      END DO
    END DO
  END DO
END SUBROUTINE compute_BCD_surf_p
