*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:33  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE COMBAK(NM,LOW,IGH,AR,AI,INT,M,ZR,ZI)
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL AR(NM,IGH),AI(NM,IGH),ZR(NM,M),ZI(NM,M)
      REAL XR,XI
      INTEGER INT(IGH)
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         MP1 = MP + 1
         DO 110 I = MP1, IGH
            XR = AR(I,MP-1)
            XI = AI(I,MP-1)
            IF (XR .EQ. 0.0 .AND. XI .EQ. 0.0) GO TO 110
            DO 100 J = 1, M
               ZR(I,J) = ZR(I,J) + XR * ZR(MP,J) - XI * ZI(MP,J)
               ZI(I,J) = ZI(I,J) + XR * ZI(MP,J) + XI * ZR(MP,J)
  100       CONTINUE
  110    CONTINUE
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
         DO 130 J = 1, M
            XR = ZR(I,J)
            ZR(I,J) = ZR(MP,J)
            ZR(MP,J) = XR
            XI = ZI(I,J)
            ZI(I,J) = ZI(MP,J)
            ZI(MP,J) = XI
  130    CONTINUE
  140 CONTINUE
  200 RETURN
      END
