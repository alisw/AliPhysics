*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:37  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE TQL1(N,D,E,IERR)
      INTEGER I,J,L,M,N,II,MML,IERR
      REAL D(N),E(N)
      REAL B,C,F,G,H,P,R,S,MACHEP
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
      F = 0.0
      B = 0.0
      E(N) = 0.0
      DO 290 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + ABS(E(L)))
         IF (B .LT. H) B = H
         DO 110 M = L, N
            IF (ABS(E(M)) .LE. B) GO TO 120
  110    CONTINUE
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
         P = (D(L+1) - D(L)) / (2.0 * E(L))
         R = SQRT(P*P+1.0)
         H = D(L) - E(L) / (P + SIGN(R,P))
         DO 140 I = L, N
  140    D(I) = D(I) - H
         F = F + H
         P = D(M)
         C = 1.0
         S = 0.0
         MML = M - L
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0)
            E(I+1) = S * E(I) * R
            S = 1.0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
  200    CONTINUE
         E(L) = S * P
         D(L) = C * P
         IF (ABS(E(L)) .GT. B) GO TO 130
  210    P = D(L) + F
         IF (L .EQ. 1) GO TO 250
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
