*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:38  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),Z(NM,N)
      REAL F,G,H,HH,SCALE
      DO 100 I = 1, N
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
      IF (N .EQ. 1) GO TO 320
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0
         SCALE = 0.0
         IF (L .LT. 2) GO TO 130
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
         IF (SCALE .NE. 0.0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / (SCALE * H)
            G = 0.0
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
         HH = F / (H + H)
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
         DO 280 K = 1, L
  280    Z(I,K) = SCALE * Z(I,K)
  290    D(I) = H
  300 CONTINUE
  320 D(1) = 0.0
      E(1) = 0.0
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0) GO TO 380
         DO 360 J = 1, L
            G = 0.0
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0
         IF (L .LT. 1) GO TO 500
         DO 400 J = 1, L
            Z(I,J) = 0.0
            Z(J,I) = 0.0
  400    CONTINUE
  500 CONTINUE
      RETURN
      END
