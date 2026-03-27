MODULE lin_vec
  USE kinds, ONLY: DP
  IMPLICIT NONE
CONTAINS
  PURE FUNCTION dot3(vec1, vec2) RESULT(retval)
    REAL(DP), INTENT(IN) :: vec1(3), vec2(3)
    REAL(DP) :: retval
    retval = SUM(vec1*vec2)
  END FUNCTION dot3
  !
  PURE FUNCTION cross3(vec1, vec2) RESULT(retval)
    REAL(DP), INTENT(IN) :: vec1(3), vec2(3)
    REAL(DP) :: retval(3)
    retval(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    retval(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    retval(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  END FUNCTION cross3
END MODULE lin_vec
