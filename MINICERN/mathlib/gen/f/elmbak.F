*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:34  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ELMBAK(NM,LOW,IGH,A,INT,M,Z)
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL A(NM,IGH),Z(NM,M)
      REAL X
      INTEGER INT(IGH)
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         MP1 = MP + 1
         DO 110 I = MP1, IGH
            X = A(I,MP-1)
            IF (X .EQ. 0.0) GO TO 110
            DO 100 J = 1, M
  100       Z(I,J) = Z(I,J) + X * Z(MP,J)
  110    CONTINUE
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
         DO 130 J = 1, M
            X = Z(I,J)
            Z(I,J) = Z(MP,J)
            Z(MP,J) = X
  130    CONTINUE
  140 CONTINUE
  200 RETURN
      END
