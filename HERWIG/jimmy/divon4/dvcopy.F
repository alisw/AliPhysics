*CMZ :          23/08/93  13.30.10  by  Jonathan Butterworth
*-- Author :
      SUBROUTINE DVCOPY(N,A,B)
      INTEGER N
      DOUBLE PRECISION A(N), B(N)
      INTEGER I
      DO 10 I=1,N
      B(I)=A(I)
 10   CONTINUE
      RETURN
      END
