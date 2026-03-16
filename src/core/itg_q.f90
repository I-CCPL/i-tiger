MODULE itg_q
  USE wannier90, ONLY: w90data
  IMPLICIT NONE
CONTAINS
  SUBROUTINE make_q_data()
    USE debug_data, ONLY: write_matrix
    CALL start_clock('make_q_data')
    !
    CALL w90data%build_Hq()
    ! CALL write_matrix('H_q.itg', w90data%Hq, w90data%kpts%nkpt, 1)
    CALL w90data%build_bvec()
    CALL w90data%build_Aq()
    !
    CALL stop_clock('make_q_data')
  END SUBROUTINE make_q_data
  !
  SUBROUTINE clear_q_data()
    CALL w90data%clear()
  END SUBROUTINE clear_q_data
END MODULE itg_q
