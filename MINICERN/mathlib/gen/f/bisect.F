*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:33  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE BISECT(N,EPS1,D,E,E2,LB,UB,MM,M,W,IND,IERR,RV4,RV5)
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,MM,M1,M2,TAG,IERR,ISTURM
      REAL D(N),E(N),E2(N),W(MM),RV4(N),RV5(N)
      REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP
      INTEGER IND(MM)
#if defined(CERNLIB_CDC)
      MACHEP=2.**(-47)
#endif
#if !defined(CERNLIB_CDC)
      MACHEP=2.**(-23)
#endif
      IERR = 0
      TAG = 0
      T1 = LB
      T2 = UB
      DO 40 I = 1, N
         IF (I .EQ. 1) GO TO 20
         IF (ABS(E(I)) .GT. MACHEP * (ABS(D(I)) + ABS(D(I-1))))
     X      GO TO 40
   20    E2(I) = 0.0
   40 CONTINUE
      P = 1
      Q = N
      X1 = UB
      ISTURM = 1
      GO TO 320
   60 M = S
      X1 = LB
      ISTURM = 2
      GO TO 320
   80 M = M - S
      IF (M .GT. MM) GO TO 980
      Q = 0
      R = 0
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0
      DO 120 Q = P, N
         X1 = U
         U = 0.0
         V = 0.0
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = MIN(D(Q)-(X1+U),XU)
         X0 = MAX(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0) GO TO 140
  120 CONTINUE
  140 X1 = MAX(ABS(XU),ABS(X0)) * MACHEP
      IF (EPS1 .LE. 0.0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * (Q-P+1)
      LB = MAX(T1,XU-X1)
      UB = MIN(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
      X0 = UB
      ISTURM = 5
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
      K = M2
  250    XU = LB
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
  300    X1 = (XU + X0) * 0.5
         IF ((X0 - XU) .LE. (2.0 * MACHEP *
     X      (ABS(XU) + ABS(X0)) + ABS(EPS1))) GO TO 420
  320    S = P - 1
         U = 1.0
         DO 340 I = P, Q
            IF (U .NE. 0.0) GO TO 325
            V = ABS(E(I)) / MACHEP
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0) S = S + 1
  340    CONTINUE
         GO TO (60,80,200,220,360), ISTURM
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
  980 IERR = 3 * N + 1
 1001 LB = T1
      UB = T2
      RETURN
      END
