*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:37  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ORTBAK(NM,LOW,IGH,A,ORT,M,Z)
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL A(NM,IGH),ORT(IGH),Z(NM,M)
      REAL G,H
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         H = A(MP,MP-1) * ORT(MP)
         IF (H .EQ. 0.0) GO TO 140
         MP1 = MP + 1
         DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
         DO 130 J = 1, M
            G = 0.0
            DO 110 I = MP, IGH
  110       G = G + ORT(I) * Z(I,J)
            G = G / H
            DO 120 I = MP, IGH
  120       Z(I,J) = Z(I,J) + G * ORT(I)
  130    CONTINUE
  140 CONTINUE
  200 RETURN
      END
