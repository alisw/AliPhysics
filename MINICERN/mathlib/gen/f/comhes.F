*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:33  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE COMHES(NM,N,LOW,IGH,AR,AI,INT)
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL AR(NM,N),AI(NM,N)
      REAL XR,XI,YR,YI
      INTEGER INT(IGH)
      COMPLEX X,Y
      REAL T1(2),T2(2)
      EQUIVALENCE (X,T1(1),XR),(T1(2),XI),(Y,T2(1),YR),(T2(2),YI)
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
      DO 180 M = KP1, LA
         MM1 = M - 1
         XR = 0.0
         XI = 0.0
         I = M
         DO 100 J = M, IGH
            IF (ABS(AR(J,MM1)) + ABS(AI(J,MM1))
     X         .LE. ABS(XR) + ABS(XI)) GO TO 100
            XR = AR(J,MM1)
            XI = AI(J,MM1)
            I = J
  100    CONTINUE
         INT(M) = I
         IF (I .EQ. M) GO TO 130
         DO 110 J = MM1, N
            YR = AR(I,J)
            AR(I,J) = AR(M,J)
            AR(M,J) = YR
            YI = AI(I,J)
            AI(I,J) = AI(M,J)
            AI(M,J) = YI
  110    CONTINUE
         DO 120 J = 1, IGH
            YR = AR(J,I)
            AR(J,I) = AR(J,M)
            AR(J,M) = YR
            YI = AI(J,I)
            AI(J,I) = AI(J,M)
            AI(J,M) = YI
  120    CONTINUE
  130    IF (XR .EQ. 0.0 .AND. XI .EQ. 0.0) GO TO 180
         MP1 = M + 1
         DO 160 I = MP1, IGH
            YR = AR(I,MM1)
            YI = AI(I,MM1)
            IF (YR .EQ. 0.0 .AND. YI .EQ. 0.0) GO TO 160
            Y = Y / X
            AR(I,MM1) = YR
            AI(I,MM1) = YI
            DO 140 J = M, N
               AR(I,J) = AR(I,J) - YR * AR(M,J) + YI * AI(M,J)
               AI(I,J) = AI(I,J) - YR * AI(M,J) - YI * AR(M,J)
  140       CONTINUE
            DO 150 J = 1, IGH
               AR(J,M) = AR(J,M) + YR * AR(J,I) - YI * AI(J,I)
               AI(J,M) = AI(J,M) + YR * AI(J,I) + YI * AR(J,I)
  150       CONTINUE
  160    CONTINUE
  180 CONTINUE
  200 RETURN
      END
