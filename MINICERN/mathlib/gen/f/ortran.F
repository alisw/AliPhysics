*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:37  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ORTRAN(NM,N,LOW,IGH,A,ORT,Z)
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL A(NM,IGH),ORT(IGH),Z(NM,N)
      REAL G,H
      DO 80 I = 1, N
         DO 60 J = 1, N
   60    Z(I,J) = 0.0
         Z(I,I) = 1.0
   80 CONTINUE
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
      DO 140 MM = 1, KL
         MP = IGH - MM
         H = A(MP,MP-1) * ORT(MP)
         IF (H .EQ. 0.0) GO TO 140
         MP1 = MP + 1
         DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
         DO 130 J = MP, IGH
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
