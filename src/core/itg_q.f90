MODULE itg_q
  USE wannier90, ONLY: w90data
  IMPLICIT NONE
CONTAINS
  SUBROUTINE make_q()
    USE io_global, ONLY: ionode, stdout
    USE debug_data, ONLY: write_matrix
    USE itg_R, ONLY: R_data
    IF (.NOT. ionode) RETURN
    CALL start_clock('make_q')
    WRITE (stdout, '(2X, A)') 'Building data in q space...'
    !
    CALL w90data%build_Hq()
    ! CALL write_matrix('H_q.itg', w90data%Hq, w90data%kpts%nkpt, 1)
    IF (R_data%bA_R) THEN
      CALL w90data%build_bvec()
      CALL w90data%build_Aq()
    END IF
    !
    CALL write_sep_line()
    CALL stop_clock('make_q')
  END SUBROUTINE make_q
  !
  SUBROUTINE bcast_q()
    !> Broadcasting q-space data to all nodes is unnecessary.
    ! CALL w90data%bcast_Hq()
    ! CALL w90data%bcast_bvec()
    ! CALL w90data%bcast_Aq()
  END SUBROUTINE bcast_q
  !
  SUBROUTINE clear_q()
    CALL w90data%clear()
  END SUBROUTINE clear_q
END MODULE itg_q
