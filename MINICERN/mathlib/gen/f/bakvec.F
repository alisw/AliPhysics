*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:32  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE BAKVEC(NM,N,T,E,M,Z,IERR)
      INTEGER I,J,M,N,NM,IERR
      REAL T(NM,3),E(N),Z(NM,M)
      IERR = 0
      E(1) = 1.0
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
         IF (E(I) .NE. 0.0) GO TO 80
         IF (T(I,1) .NE. 0.0 .OR. T(I-1,3) .NE. 0.0) GO TO 1000
         E(I) = 1.0
         GO TO 100
   80    E(I) = E(I-1) * E(I) / T(I-1,3)
  100 CONTINUE
      DO 120 J = 1, M
         DO 120 I = 2, N
         Z(I,J) = Z(I,J) * E(I)
  120 CONTINUE
      GO TO 1001
 1000 IERR = 2 * N + I
 1001 RETURN
      END
