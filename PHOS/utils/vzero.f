      SUBROUTINE VZERO (A,N)
      DIMENSION A(*)
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 A(I)= 0.
      RETURN
      END
