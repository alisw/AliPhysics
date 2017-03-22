*CMZ :          23/08/93  13.30.11  by  Jonathan Butterworth
*-- Author :
      DOUBLE PRECISION FUNCTION DVDOT(N, A, B)
      INTEGER N
      DOUBLE PRECISION A(N), B(N)
      INTEGER I
      DOUBLE PRECISION SUM
      DVDOT=0.0D+0
      IF(N.LE.0) RETURN
      SUM=0.0D+0
      DO 10 I=1,N
      SUM=SUM+A(I)*B(I)
 10   CONTINUE
      DVDOT=SUM
      RETURN
      END
