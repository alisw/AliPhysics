*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:32  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL A(NM,N),SCALE(N)
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
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   30 CONTINUE
      DO 40 I = K, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
   40 CONTINUE
   50 GO TO (80,130), IEXC
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (A(J,I) .NE. 0.0) GO TO 120
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
            IF (A(I,J) .NE. 0.0) GO TO 170
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
            C = C + ABS(A(J,I))
            R = R + ABS(A(I,J))
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
  250    A(I,J) = A(I,J) * G
         DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F
  270 CONTINUE
      IF (NOCONV) GO TO 190
  280 LOW = K
      IGH = L
      RETURN
      END
