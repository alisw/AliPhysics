*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:35  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE FIGI2(NM,N,T,D,E,Z,IERR)
      INTEGER I,J,N,NM,IERR
      REAL T(NM,3),D(N),E(N),Z(NM,N)
      REAL H
      IERR = 0
      DO 100 I = 1, N
         DO 50 J = 1, N
   50    Z(I,J) = 0.0
         IF (I .EQ. 1) GO TO 70
         H = T(I,1) * T(I-1,3)
         IF (H) 900, 60, 80
   60    IF (T(I,1) .NE. 0.0 .OR. T(I-1,3) .NE. 0.0) GO TO 1000
         E(I) = 0.0
   70    Z(I,I) = 1.0
         GO TO 90
   80    E(I) = SQRT(H)
         Z(I,I) = Z(I-1,I-1) * E(I) / T(I-1,3)
   90    D(I) = T(I,2)
  100 CONTINUE
      GO TO 1001
  900 IERR = N + I
      GO TO 1001
 1000 IERR = 2 * N + I
 1001 RETURN
      END
