*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:37  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RATQR(N,EPS1,D,E,E2,M,W,IND,BD,TYPE,IDEF,IERR)
      INTEGER I,J,K,M,N,II,JJ,K1,IDEF,IERR,JDEF
      REAL D(N),E(N),E2(N),W(N),BD(N)
      REAL F,P,Q,R,S,EP,QP,ERR,TOT,EPS1,DELTA,MACHEP
      INTEGER IND(M)
      LOGICAL TYPE
#if defined(CERNLIB_CDC)
      MACHEP=2.**(-47)
#endif
#if !defined(CERNLIB_CDC)
      MACHEP=2.**(-23)
#endif
      IERR = 0
      JDEF = IDEF
      DO 20 I = 1, N
   20 W(I) = D(I)
      IF (TYPE) GO TO 40
      J = 1
      GO TO 400
   40 ERR = 0.0
      S = 0.0
      TOT = W(1)
      Q = 0.0
      J = 0
      DO 100 I = 1, N
         P = Q
         IF (I .EQ. 1) GO TO 60
         IF (P .GT. MACHEP * (ABS(D(I)) + ABS(D(I-1)))) GO TO 80
   60    E2(I) = 0.0
         J = J + 1
   80    BD(I) = E2(I)
         IND(I) = J
         Q = 0.0
         IF (I .NE. N) Q = ABS(E(I+1))
         TOT = MIN(W(I)-P-Q,TOT)
  100 CONTINUE
      IF (JDEF .EQ. 1 .AND. TOT .LT. 0.0) GO TO 140
      DO 110 I = 1, N
  110 W(I) = W(I) - TOT
      GO TO 160
  140 TOT = 0.0
  160 DO 360 K = 1, M
  180    TOT = TOT + S
         DELTA = W(N) - S
         I = N
         F = ABS(MACHEP*TOT)
         IF (EPS1 .LT. F) EPS1 = F
         IF (DELTA .GT. EPS1) GO TO 190
         IF (DELTA .LT. (-EPS1)) GO TO 1000
         GO TO 300
  190    IF (K .EQ. N) GO TO 210
         K1 = K + 1
         DO 200 J = K1, N
            IF (BD(J) .LE. (MACHEP*(W(J)+W(J-1))) ** 2) BD(J) = 0.0
  200    CONTINUE
  210    F = BD(N) / DELTA
         QP = DELTA + F
         P = 1.0
         IF (K .EQ. N) GO TO 260
         K1 = N - K
         DO 240 II = 1, K1
            I = N - II
            Q = W(I) - S - F
            R = Q / QP
            P = P * R + 1.0
            EP = F * R
            W(I+1) = QP + EP
            DELTA = Q - EP
            IF (DELTA .GT. EPS1) GO TO 220
            IF (DELTA .LT. (-EPS1)) GO TO 1000
            GO TO 300
  220       F = BD(I) / Q
            QP = DELTA + F
            BD(I+1) = QP * EP
  240    CONTINUE
  260    W(K) = QP
         S = QP / P
         IF (TOT + S .GT. TOT) GO TO 180
         IERR = 5 * N + K
         S = 0.0
         DELTA = QP
         DO 280 J = K, N
            IF (W(J) .GT. DELTA) GO TO 280
            I = J
            DELTA = W(J)
  280    CONTINUE
  300    IF (I .LT. N) BD(I+1) = BD(I) * F / QP
         II = IND(I)
         IF (I .EQ. K) GO TO 340
         K1 = I - K
         DO 320 JJ = 1, K1
            J = I - JJ
            W(J+1) = W(J) - S
            BD(J+1) = BD(J)
            IND(J+1) = IND(J)
  320    CONTINUE
  340    W(K) = TOT
         ERR = ERR + ABS(DELTA)
         BD(K) = ERR
         IND(K) = II
  360 CONTINUE
      IF (TYPE) GO TO 1001
      F = BD(1)
      E2(1) = 2.0
      BD(1) = F
      J = 2
  400 DO 500 I = 1, N
  500 W(I) = -W(I)
      JDEF = -JDEF
      GO TO (40,1001), J
 1000 IERR = 6 * N + 1
 1001 RETURN
      END
