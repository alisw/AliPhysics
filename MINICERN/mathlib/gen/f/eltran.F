*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:34  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ELTRAN(NM,N,LOW,IGH,A,INT,Z)
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL A(NM,IGH),Z(NM,N)
      INTEGER INT(IGH)
      DO 80 I = 1, N
         DO 60 J = 1, N
   60    Z(I,J) = 0.0
         Z(I,I) = 1.0
   80 CONTINUE
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
      DO 140 MM = 1, KL
         MP = IGH - MM
         MP1 = MP + 1
         DO 100 I = MP1, IGH
  100    Z(I,MP) = A(I,MP-1)
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
         DO 130 J = MP, IGH
            Z(MP,J) = Z(I,J)
            Z(I,J) = 0.0
  130    CONTINUE
         Z(I,MP) = 1.0
  140 CONTINUE
  200 RETURN
      END
