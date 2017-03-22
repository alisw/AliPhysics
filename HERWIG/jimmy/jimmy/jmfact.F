*CMZ :          31/08/93  17.51.25  by  Jonathan Butterworth
*-- Author :
      DOUBLE PRECISION FUNCTION JMFACT( N )
C --------------------------------------------------------------------------
C Factorial of N, given as double precision to allow for very large numbers.
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N, I
      DOUBLE PRECISION TOTAL

      TOTAL = 1.D0
      DO I = 2, N
        TOTAL = TOTAL*DBLE(I)
      ENDDO
      JMFACT = TOTAL
      RETURN
      END
