*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:36  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE IMTQL1(N,D,E,IERR)
      INTEGER I,J,L,M,N,II,MML,IERR
      REAL D(N),E(N)
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
      DO 290 L = 1, N
         J = 0
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            IF (ABS(E(M)) .LE. MACHEP * (ABS(D(M)) + ABS(D(M+1))))
     X         GO TO 120
  110    CONTINUE
  120    P = D(L)
         IF (M .EQ. L) GO TO 215
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
  200    CONTINUE
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0
         GO TO 105
  215    IF (L .EQ. 1) GO TO 250
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
  250    I = 1
  270    D(I) = P
  290 CONTINUE
      GO TO 1001
 1000 IERR = L
 1001 RETURN
      END
