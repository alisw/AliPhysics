*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:33  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE CINVIT(NM,N,AR,AI,WR,WI,SELECT,MM,M,ZR,ZI,
     X                  IERR,RM1,RM2,RV1,RV2)
      INTEGER I,J,K,M,N,S,II,MM,MP,NM,UK,IP1,ITS,KM1,IERR
      REAL AR(NM,N),AI(NM,N),WR(N),WI(N),ZR(NM,MM),ZI(NM,MM),
     X       RM1(N,N),RM2(N,N),RV1(N),RV2(N)
      REAL X,Y,EPS3,NORM,NORMV,GROWTO,ILAMBD,MACHEP,RLAMBD,UKROOT
      LOGICAL SELECT(N)
      COMPLEX Z3
      REAL T3(2)
      EQUIVALENCE (Z3,T3(1))
#if defined(CERNLIB_CDC)
      MACHEP=2.**(-47)
#endif
#if !defined(CERNLIB_CDC)
      MACHEP=2.**(-23)
#endif
      IERR = 0
      UK = 0
      S = 1
      DO 980 K = 1, N
         IF (.NOT. SELECT(K)) GO TO 980
         IF (S .GT. MM) GO TO 1000
         IF (UK .GE. K) GO TO 200
         DO 120 UK = K, N
            IF (UK .EQ. N) GO TO 140
            IF (AR(UK+1,UK) .EQ. 0.0 .AND. AI(UK+1,UK) .EQ. 0.0)
     X         GO TO 140
  120    CONTINUE
  140    NORM = 0.0
         MP = 1
         DO 180 I = 1, UK
            X = 0.0
            DO 160 J = MP, UK
  160       X = X + ABS(CMPLX(AR(I,J),AI(I,J)))
            IF (X .GT. NORM) NORM = X
            MP = I
  180    CONTINUE
         IF (NORM .EQ. 0.0) NORM = 1.0
         EPS3 = MACHEP * NORM
         UKROOT = SQRT(REAL(UK))
         GROWTO = 1.0E-1 / UKROOT
  200    RLAMBD = WR(K)
         ILAMBD = WI(K)
         IF (K .EQ. 1) GO TO 280
         KM1 = K - 1
         GO TO 240
  220    RLAMBD = RLAMBD + EPS3
  240    DO 260 II = 1, KM1
            I = K - II
            IF (SELECT(I) .AND. ABS(WR(I)-RLAMBD) .LT. EPS3 .AND.
     X         ABS(WI(I)-ILAMBD) .LT. EPS3) GO TO 220
  260    CONTINUE
         WR(K) = RLAMBD
  280    MP = 1
         DO 320 I = 1, UK
            DO 300 J = MP, UK
               RM1(I,J) = AR(I,J)
               RM2(I,J) = AI(I,J)
  300       CONTINUE
            RM1(I,I) = RM1(I,I) - RLAMBD
            RM2(I,I) = RM2(I,I) - ILAMBD
            MP = I
            RV1(I) = EPS3
  320    CONTINUE
         IF (UK .EQ. 1) GO TO 420
         DO 400 I = 2, UK
            MP = I - 1
            IF (ABS(CMPLX(RM1(I,MP),RM2(I,MP))) .LE.
     X          ABS(CMPLX(RM1(MP,MP),RM2(MP,MP)))) GO TO 360
            DO 340 J = MP, UK
               Y = RM1(I,J)
               RM1(I,J) = RM1(MP,J)
               RM1(MP,J) = Y
               Y = RM2(I,J)
               RM2(I,J) = RM2(MP,J)
               RM2(MP,J) = Y
  340       CONTINUE
  360       IF (RM1(MP,MP) .EQ. 0.0 .AND. RM2(MP,MP) .EQ. 0.0)
     X         RM1(MP,MP) = EPS3
            Z3 = CMPLX(RM1(I,MP),RM2(I,MP)) /
     X             CMPLX(RM1(MP,MP),RM2(MP,MP))
            X = T3(1)
            Y = T3(2)
            IF (X .EQ. 0.0 .AND. Y .EQ. 0.0) GO TO 400
            DO 380 J = I, UK
               RM1(I,J) = RM1(I,J) - X * RM1(MP,J) + Y * RM2(MP,J)
               RM2(I,J) = RM2(I,J) - X * RM2(MP,J) - Y * RM1(MP,J)
  380       CONTINUE
  400    CONTINUE
  420    IF (RM1(UK,UK) .EQ. 0.0 .AND. RM2(UK,UK) .EQ. 0.0)
     X      RM1(UK,UK) = EPS3
         ITS = 0
  660    DO 720 II = 1, UK
            I = UK + 1 - II
            X = RV1(I)
            Y = 0.0
            IF (I .EQ. UK) GO TO 700
            IP1 = I + 1
            DO 680 J = IP1, UK
               X = X - RM1(I,J) * RV1(J) + RM2(I,J) * RV2(J)
               Y = Y - RM1(I,J) * RV2(J) - RM2(I,J) * RV1(J)
  680       CONTINUE
  700       Z3 = CMPLX(X,Y) / CMPLX(RM1(I,I),RM2(I,I))
            RV1(I) = T3(1)
            RV2(I) = T3(2)
  720    CONTINUE
         ITS = ITS + 1
         NORM = 0.0
         NORMV = 0.0
         DO 780 I = 1, UK
            X = ABS(CMPLX(RV1(I),RV2(I)))
            IF (NORMV .GE. X) GO TO 760
            NORMV = X
            J = I
  760       NORM = NORM + X
  780    CONTINUE
         IF (NORM .LT. GROWTO) GO TO 840
         X = RV1(J)
         Y = RV2(J)
         DO 820 I = 1, UK
            Z3 = CMPLX(RV1(I),RV2(I)) / CMPLX(X,Y)
            ZR(I,S) = T3(1)
            ZI(I,S) = T3(2)
  820    CONTINUE
         IF (UK .EQ. N) GO TO 940
         J = UK + 1
         GO TO 900
  840    IF (ITS .GE. UK) GO TO 880
         X = UKROOT
         Y = EPS3 / (X + 1.0)
         RV1(1) = EPS3
         DO 860 I = 2, UK
  860    RV1(I) = Y
         J = UK - ITS + 1
         RV1(J) = RV1(J) - EPS3 * X
         GO TO 660
  880    J = 1
         IERR = -K
  900    DO 920 I = J, N
            ZR(I,S) = 0.0
            ZI(I,S) = 0.0
  920    CONTINUE
  940    S = S + 1
  980 CONTINUE
      GO TO 1001
 1000 IF (IERR .NE. 0) IERR = IERR - N
      IF (IERR .EQ. 0) IERR = -(2 * N + 1)
 1001 M = S - 1
      RETURN
      END
