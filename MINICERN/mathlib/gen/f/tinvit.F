*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:37  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV6)
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      REAL D(N),E(N),E2(N),W(M),Z(NM,M),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      REAL U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,MACHEP
      INTEGER IND(M)
#if defined(CERNLIB_CDC)
      MACHEP=2.**(-47)
#endif
#if !defined(CERNLIB_CDC)
      MACHEP=2.**(-23)
#endif
      IERR = 0
      TAG = 0
      ORDER = 1.0 - E2(1)
      Q = 0
  100 P = Q + 1
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. 0.0) GO TO 140
  120 CONTINUE
  140 TAG = TAG + 1
      S = 0
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
         XU = 1.0
         IF (P .NE. Q) GO TO 490
         RV6(P) = 1.0
         GO TO 870
  490    NORM = ABS(D(P))
         IP = P + 1
         DO 500 I = IP, Q
  500    NORM = NORM + ABS(D(I)) + ABS(E(I))
         EPS2 = 1.0E-3 * NORM
         EPS3 = MACHEP * NORM
         UK = Q-P+1
         EPS4 = UK * EPS3
         UK = EPS4 / SQRT(UK)
         S = P
  505    GROUP = 0
         GO TO 520
  510    IF (ABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0) X1 = X0 + ORDER * EPS3
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
         J = R
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
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
         IF (ITS .EQ. 5) GO TO 830
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
  830    IERR = -R
         XU = 0.0
         GO TO 870
  840    U = 0.0
         DO 860 I = P, Q
  860    U = U + RV6(I)**2
         XU = 1.0 / SQRT(U)
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
         X0 = X1
  920 CONTINUE
      IF (Q .LT. N) GO TO 100
      RETURN
      END
