*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:34  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL A(NM,N)
      REAL X,Y
      INTEGER INT(IGH)
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
      DO 180 M = KP1, LA
         MM1 = M - 1
         X = 0.0
         I = M
         DO 100 J = M, IGH
            IF (ABS(A(J,MM1)) .LE. ABS(X)) GO TO 100
            X = A(J,MM1)
            I = J
  100    CONTINUE
         INT(M) = I
         IF (I .EQ. M) GO TO 130
         DO 110 J = MM1, N
            Y = A(I,J)
            A(I,J) = A(M,J)
            A(M,J) = Y
  110    CONTINUE
         DO 120 J = 1, IGH
            Y = A(J,I)
            A(J,I) = A(J,M)
            A(J,M) = Y
  120    CONTINUE
  130    IF (X .EQ. 0.0) GO TO 180
         MP1 = M + 1
         DO 160 I = MP1, IGH
            Y = A(I,MM1)
            IF (Y .EQ. 0.0) GO TO 160
            Y = Y / X
            A(I,MM1) = Y
            DO 140 J = M, N
  140       A(I,J) = A(I,J) - Y * A(M,J)
            DO 150 J = 1, IGH
  150       A(J,M) = A(J,M) + Y * A(J,I)
  160    CONTINUE
  180 CONTINUE
  200 RETURN
      END
