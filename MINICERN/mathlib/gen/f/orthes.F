*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:37  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ORTHES(NM,N,LOW,IGH,A,ORT)
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      REAL A(NM,N),ORT(IGH)
      REAL F,G,H,SCALE
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
      DO 180 M = KP1, LA
         H = 0.0
         ORT(M) = 0.0
         SCALE = 0.0
         DO 90 I = M, IGH
   90    SCALE = SCALE + ABS(A(I,M-1))
         IF (SCALE .EQ. 0.0) GO TO 180
         MP = M + IGH
         DO 100 II = M, IGH
            I = MP - II
            ORT(I) = A(I,M-1) / SCALE
            H = H + ORT(I) * ORT(I)
  100    CONTINUE
         G = -SIGN(SQRT(H),ORT(M))
         H = H - ORT(M) * G
         ORT(M) = ORT(M) - G
         DO 130 J = M, N
            F = 0.0
            DO 110 II = M, IGH
               I = MP - II
               F = F + ORT(I) * A(I,J)
  110       CONTINUE
            F = F / H
            DO 120 I = M, IGH
  120       A(I,J) = A(I,J) - F * ORT(I)
  130    CONTINUE
         DO 160 I = 1, IGH
            F = 0.0
            DO 140 JJ = M, IGH
               J = MP - JJ
               F = F + ORT(J) * A(I,J)
  140       CONTINUE
            F = F / H
            DO 150 J = M, IGH
  150       A(I,J) = A(I,J) - F * ORT(J)
  160    CONTINUE
         ORT(M) = SCALE * ORT(M)
         A(M,M-1) = SCALE * G
  180 CONTINUE
  200 RETURN
      END
