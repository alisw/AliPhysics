*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:32  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL SCALE(N),Z(NM,M)
      REAL S
      IF (IGH .EQ. LOW) GO TO 120
      DO 110 I = LOW, IGH
         S = SCALE(I)
         DO 100 J = 1, M
  100    Z(I,J) = Z(I,J) * S
  110 CONTINUE
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
         DO 130 J = 1, M
            S = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = S
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
