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
  !
  SUBROUTINE add_comma(str)
    CHARACTER(LEN=*), INTENT(INOUT)::str
    CHARACTER(LEN=30)::tmp
    INTEGER::len1, len2, i, j, q, r
    len1 = LEN_TRIM(str)
    i = 1
    IF (str(1:1) == '-') THEN
      tmp(1:1) = '-'
      i = 2
    END IF

    len2 = 0
    DO WHILE (i <= len1)
      IF (str(i:i) < '0' .OR. str(i:i) > '9') THEN
        EXIT
      END IF
      len2 = len2 + 1
      len1 = len1 - 1
    END DO

    j = i
    q = (len2 - 1)/3
    r = MOD(len2 - 1, 3)
    IF (q /= 0) THEN
      tmp(j:j + r) = str(i:i + r)
      j = j + r + 1
      i = i + r + 1
    END IF
    DO WHILE (q > 0)
      q = q - 1
      tmp(j:j) = ','
      tmp(j + 1:j + 4) = str(i:i + 3)
      j = j + 4
      i = i + 3
    END DO
    str = tmp
  END SUBROUTINE add_comma
END MODULE char_mod
