*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:38  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),E2(N)
      REAL F,G,H,SCALE
      DO 100 I = 1, N
  100 D(I) = A(I,I)
      DO  300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0
         SCALE = 0.0
         IF (L .LT. 1) GO TO 130
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(A(I,K))
         IF (SCALE .NE. 0.0) GO TO 140
  130    E(I) = 0.0
         E2(I) = 0.0
         GO TO 290
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            H = H + A(I,K) * A(I,K)
  150    CONTINUE
         E2(I) = SCALE * SCALE * H
         F = A(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         A(I,L) = F - G
         IF (L .EQ. 1) GO TO 270
         F = 0.0
         DO 240 J = 1, L
            G = 0.0
            DO 180 K = 1, J
  180       G = G + A(J,K) * A(I,K)
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
  200       G = G + A(K,J) * A(I,K)
  220       E(J) = G / H
            F = F + E(J) * A(I,J)
  240    CONTINUE
         H = F / (H + H)
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - H * F
            E(J) = G
            DO 260 K = 1, J
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
  260    CONTINUE
  270    DO 280 K = 1, L
  280    A(I,K) = SCALE * A(I,K)
  290    H = D(I)
         D(I) = A(I,I)
         A(I,I) = H
  300 CONTINUE
      RETURN
      END
