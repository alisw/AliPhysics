*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:37  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE IMTQL2(NM,N,D,E,Z,IERR)
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
      REAL D(N),E(N),Z(NM,N)
      REAL B,C,F,G,P,R,S,MACHEP
#if defined(CERNLIB_CDC)
      MACHEP=2.**(-47)
#endif
#if !defined(CERNLIB_CDC)
      MACHEP=2.**(-23)
#endif
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
  100 E(I-1) = E(I)
      E(N) = 0.0
      DO 240 L = 1, N
         J = 0
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            IF (ABS(E(M)) .LE. MACHEP * (ABS(D(M)) + ABS(D(M+1))))
     X         GO TO 120
  110    CONTINUE
  120    P = D(L)
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
         G = (D(L+1) - P) / (2.0 * E(L))
         R = SQRT(G*G+1.0)
         G = D(M) - P + E(L) / (G + SIGN(R,G))
         S = 1.0
         C = 1.0
         P = 0.0
         MML = M - L
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF (ABS(F) .LT. ABS(G)) GO TO 150
            C = G / F
            R = SQRT(C*C+1.0)
            E(I+1) = F * R
            S = 1.0 / R
            C = C * S
            GO TO 160
  150       S = F / G
            R = SQRT(S*S+1.0)
            E(I+1) = G * R
            C = 1.0 / R
            S = S * C
  160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
            DO 180 K = 1, N
               F = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * F
               Z(K,I) = C * Z(K,I) - S * F
  180       CONTINUE
  200    CONTINUE
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0
         GO TO 105
  240 CONTINUE
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
  300 CONTINUE
      GO TO 1001
 1000 IERR = L
 1001 RETURN
      END
