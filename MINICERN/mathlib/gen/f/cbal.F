*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:33  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE CBAL(NM,N,AR,AI,LOW,IGH,SCALE)
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL AR(NM,N),AI(NM,N),SCALE(N)
      REAL C,F,G,R,S,B2,RADIX
      LOGICAL NOCONV
      RADIX = 2.
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
      DO 30 I = 1, L
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   30 CONTINUE
      DO 40 I = K, N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   40 CONTINUE
   50 GO TO (80,130), IEXC
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (AR(J,I) .NE. 0.0 .OR. AI(J,I) .NE. 0.0) GO TO 120
  110    CONTINUE
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
      GO TO 140
  130 K = K + 1
  140 DO 170 J = K, L
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (AR(I,J) .NE. 0.0 .OR. AI(I,J) .NE. 0.0) GO TO 170
  150    CONTINUE
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
      DO 180 I = K, L
  180 SCALE(I) = 1.0
  190 NOCONV = .FALSE.
      DO 270 I = K, L
         C = 0.0
         R = 0.0
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + ABS(AR(J,I)) + ABS(AI(J,I))
            R = R + ABS(AR(I,J)) + ABS(AI(I,J))
  200    CONTINUE
         G = R / RADIX
         F = 1.0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
  240    IF ((C + R) / F .GE. 0.95 * S) GO TO 270
         G = 1.0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
         DO 250 J = K, N
            AR(I,J) = AR(I,J) * G
            AI(I,J) = AI(I,J) * G
  250    CONTINUE
         DO 260 J = 1, L
            AR(J,I) = AR(J,I) * F
            AI(J,I) = AI(J,I) * F
  260    CONTINUE
  270 CONTINUE
      IF (NOCONV) GO TO 190
  280 LOW = K
      IGH = L
      RETURN
      END
