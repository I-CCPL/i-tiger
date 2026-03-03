MODULE char_mod
  IMPLICIT NONE
CONTAINS
  FUNCTION captital(c_in) RESULT(c_out)
    CHARACTER(LEN=1), INTENT(IN) :: c_in
    CHARACTER(LEN=1):: c_out
    CHARACTER(LEN=26)::lower = 'abcdefghijklmnopqrstuvwxyz', &
                        upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    INTEGER::i
    DO i = 1, 26
      IF (c_in == lower(i:i)) THEN
        c_out = upper(i:i)
        RETURN
      END IF
    END DO
    c_out = c_in
  END FUNCTION captital
  !
  FUNCTION match(c1, c2) RESULT(bMatch)
    !< .TRUE. if c2 contains c1, .FALSE. otherwise.
    CHARACTER(LEN=*), INTENT(IN)::c1, c2
    LOGICAL::bMatch
    INTEGER::len1, len2, i
    len1 = LEN_TRIM(c1)
    len2 = LEN_TRIM(c2)

    DO i = 1, (len2 - len1 + 1)
      IF (c1(1:len1) == c2(i:(i + len1 - 1))) THEN
        bMatch = .TRUE.
        RETURN
      END IF
    END DO
    bMatch = .FALSE.
  END FUNCTION match
END MODULE char_mod
