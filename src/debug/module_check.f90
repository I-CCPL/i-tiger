MODULE data_check
  USE kinds, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  TYPE, PUBLIC::check_data
    INTEGER::l_idx1 
    !< Number of loof index
    INTEGER::l_idx2
    INTEGER::l_idx3
    INTEGER::l_idx4
  END TYPE check_data

  TYPE(check_data), PUBLIC::d_chk
 
END MODULE
