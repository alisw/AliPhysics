*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:34  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE FIGI(NM,N,T,D,E,E2,IERR)
      INTEGER I,N,NM,IERR
      REAL T(NM,3),D(N),E(N),E2(N)
      IERR = 0
      DO 100 I = 1, N
         IF (I .EQ. 1) GO TO 90
         E2(I) = T(I,1) * T(I-1,3)
         IF (E2(I)) 1000, 60, 80
   60    IF (T(I,1) .EQ. 0.0 .AND. T(I-1,3) .EQ. 0.0) GO TO 80
         IERR = -(3 * N + I)
   80    E(I) = SQRT(E2(I))
   90    D(I) = T(I,2)
  100 CONTINUE
      GO TO 1001
 1000 IERR = N + I
 1001 RETURN
      END
