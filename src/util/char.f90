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
    CHARACTER(LEN=LEN(str))::tmp
    INTEGER::first, in_pos, out_pos, len_num, n_digits

    tmp = ' '
    first = 1
    IF (str(1:1) == '-') THEN
      tmp(1:1) = '-'
      first = 2
    END IF

    len_num = 0
    DO in_pos = first, LEN_TRIM(str)
      IF (str(in_pos:in_pos) < '0' .OR. str(in_pos:in_pos) > '9') EXIT
      len_num = len_num + 1
    END DO

    IF (len_num == 0) RETURN

    out_pos = first + len_num + (len_num - 1)/3 - 1
    in_pos = first + len_num - 1
    n_digits = 0

    DO WHILE (in_pos >= first)
      tmp(out_pos:out_pos) = str(in_pos:in_pos)
      out_pos = out_pos - 1
      in_pos = in_pos - 1
      n_digits = n_digits + 1

      IF (MOD(n_digits, 3) == 0 .AND. in_pos >= first) THEN
        tmp(out_pos:out_pos) = ','
        out_pos = out_pos - 1
      END IF
    END DO

    str = ' '
    str(1:LEN_TRIM(tmp)) = tmp(1:LEN_TRIM(tmp))
  END SUBROUTINE add_comma
END MODULE char_mod
