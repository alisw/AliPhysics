*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:38  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE TSTURM(NM,N,EPS1,D,E,E2,LB,UB,MM,M,W,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV5,RV6)
      INTEGER I,J,K,M,N,P,Q,R,S,II,IP,JJ,MM,M1,M2,NM,ITS,
     X        IERR,GROUP,ISTURM
      REAL D(N),E(N),E2(N),W(MM),Z(NM,MM),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV5(N),RV6(N)
      REAL U,V,LB,T1,T2,UB,UK,XU,X0,X1,EPS1,EPS2,EPS3,EPS4,
     X       NORM,MACHEP
#if defined(CERNLIB_CDC)
      MACHEP=2.**(-47)
#endif
#if !defined(CERNLIB_CDC)
      MACHEP=2.**(-23)
#endif
      IERR = 0
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
      R = R + 1
      DO 160 I = 1, N
  160 Z(I,R) = 0.0
      W(R) = D(P)
      Z(P,R) = 1.0
      GO TO 940
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
      NORM = ABS(D(P))
      IP = P + 1
      DO 500 I = IP, Q
  500 NORM = NORM + ABS(D(I)) + ABS(E(I))
      EPS2 = 1.0E-3 * NORM
      EPS3 = MACHEP * NORM
      UK = Q-P+1
      EPS4 = UK * EPS3
      UK = EPS4 / SQRT(UK)
      GROUP = 0
      S = P
      DO 920 K = M1, M2
         R = R + 1
         ITS = 1
         W(R) = RV5(K)
         X1 = RV5(K)
         IF (K .EQ. M1) GO TO 520
         IF (X1 - X0 .GE. EPS2) GROUP = -1
         GROUP = GROUP + 1
         IF (X1 .LE. X0) X1 = X0 + EPS3
  520    V = 0.0
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (ABS(E(I)) .LT. ABS(U)) GO TO 540
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
         IF (U .EQ. 0.0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0
         RV3(Q) = 0.0
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
         IF (GROUP .EQ. 0) GO TO 700
         DO 680 JJ = 1, GROUP
            J = R - GROUP - 1 + JJ
            XU = 0.0
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
  680    CONTINUE
  700    NORM = 0.0
         DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
         IF (NORM .GE. 1.0) GO TO 840
         IF (ITS .EQ. 5) GO TO 960
         IF (NORM .NE. 0.0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
  780    DO 820 I = IP, Q
            U = RV6(I)
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
         ITS = ITS + 1
         GO TO 600
  840    U = 0.0
         DO 860 I = P, Q
  860    U = U + RV6(I)**2
         XU = 1.0 / SQRT(U)
         DO 880 I = 1, N
  880    Z(I,R) = 0.0
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
         X0 = X1
  920 CONTINUE
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
  960 IERR = 4 * N + R
      GO TO 1001
  980 IERR = 3 * N + 1
 1001 LB = T1
      UB = T2
      RETURN
      END
