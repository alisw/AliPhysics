*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:37  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE INVIT(NM,N,A,WR,WI,SELECT,MM,M,Z,IERR,RM3,RV1,RV2)
      INTEGER I,J,K,M,N,S,II,IP,MM,MP,NM,UK,IP1,ITS,KM1,IERR
      REAL A(NM,N),WR(N),WI(N),Z(NM,MM),RM3(N,*),RV1(N),RV2(N)
      REAL W,X,Y,EPS3,NORM,NORMV,GROWTO,ILAMBD,MACHEP,RLAMBD,UKROOT
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
      IP = 0
      DO 980 K = 1, N
         IF (WI(K) .EQ. 0.0 .OR. IP .LT. 0) GO TO 100
         IP = 1
         IF (SELECT(K) .AND. SELECT(K+1)) SELECT(K+1) = .FALSE.
  100    IF (.NOT. SELECT(K)) GO TO 960
         IF (WI(K) .NE. 0.0) S = S + 1
         IF (S .GT. MM) GO TO 1000
         IF (UK .GE. K) GO TO 200
         DO 120 UK = K, N
            IF (UK .EQ. N) GO TO 140
            IF (A(UK+1,UK) .EQ. 0.0) GO TO 140
  120    CONTINUE
  140    NORM = 0.0
         MP = 1
         DO 180 I = 1, UK
            X = 0.0
            DO 160 J = MP, UK
  160       X = X + ABS(A(I,J))
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
         IP1 = K + IP
         WR(IP1) = RLAMBD
  280    MP = 1
         DO 320 I = 1, UK
            DO 300 J = MP, UK
  300       RM3(J,I) = A(I,J)
            RM3(I,I) = RM3(I,I) - RLAMBD
            MP = I
            RV1(I) = EPS3
  320    CONTINUE
         ITS = 0
         IF (ILAMBD .NE. 0.0) GO TO 520
         IF (UK .EQ. 1) GO TO 420
         DO 400 I = 2, UK
            MP = I - 1
            IF (ABS(RM3(MP,I)) .LE. ABS(RM3(MP,MP))) GO TO 360
            DO 340 J = MP, UK
               Y = RM3(J,I)
               RM3(J,I) = RM3(J,MP)
               RM3(J,MP) = Y
  340       CONTINUE
  360       IF (RM3(MP,MP) .EQ. 0.0) RM3(MP,MP) = EPS3
            X = RM3(MP,I) / RM3(MP,MP)
            IF (X .EQ. 0.0) GO TO 400
            DO 380 J = I, UK
  380       RM3(J,I) = RM3(J,I) - X * RM3(J,MP)
  400    CONTINUE
  420    IF (RM3(UK,UK) .EQ. 0.0) RM3(UK,UK) = EPS3
  440    DO 500 II = 1, UK
            I = UK + 1 - II
            Y = RV1(I)
            IF (I .EQ. UK) GO TO 480
            IP1 = I + 1
            DO 460 J = IP1, UK
  460       Y = Y - RM3(J,I) * RV1(J)
  480       RV1(I) = Y / RM3(I,I)
  500    CONTINUE
         GO TO 740
  520    RM3(1,3) = -ILAMBD
         DO 540 I = 2, UK
  540    RM3(1,I+2) = 0.0
         DO 640 I = 2, UK
            MP = I - 1
            W = RM3(MP,I)
            X = RM3(MP,MP) * RM3(MP,MP) + RM3(MP,I+1) * RM3(MP,I+1)
            IF (W * W .LE. X) GO TO 580
            X = RM3(MP,MP) / W
            Y = RM3(MP,I+1) / W
            RM3(MP,MP) = W
            RM3(MP,I+1) = 0.0
            DO 560 J = I, UK
               W = RM3(J,I)
               RM3(J,I) = RM3(J,MP) - X * W
               RM3(J,MP) = W
               RM3(I,J+2) = RM3(MP,J+2) - Y * W
               RM3(MP,J+2) = 0.0
  560       CONTINUE
            RM3(MP,I+2) = -ILAMBD
            RM3(I,I) = RM3(I,I) - Y * ILAMBD
            RM3(I,I+2) = RM3(I,I+2) + X * ILAMBD
            GO TO 640
  580       IF (X .NE. 0.0) GO TO 600
            RM3(MP,MP) = EPS3
            RM3(MP,I+1) = 0.0
            X = EPS3 * EPS3
  600       W = W / X
            X = RM3(MP,MP) * W
            Y = -RM3(MP,I+1) * W
            DO 620 J = I, UK
               RM3(J,I) = RM3(J,I) - X * RM3(J,MP) + Y * RM3(MP,J+2)
               RM3(I,J+2) = -X * RM3(MP,J+2) - Y * RM3(J,MP)
  620       CONTINUE
            RM3(I,I+2) = RM3(I,I+2) - ILAMBD
  640    CONTINUE
         IF (RM3(UK,UK) .EQ. 0.0 .AND.
     X      RM3(UK,UK+2) .EQ. 0.0) RM3(UK,UK) = EPS3
  660    DO 720 II = 1, UK
            I = UK + 1 - II
            X = RV1(I)
            Y = 0.0
            IF (I .EQ. UK) GO TO 700
            IP1 = I + 1
            DO 680 J = IP1, UK
               X = X - RM3(J,I) * RV1(J) + RM3(I,J+2) * RV2(J)
               Y = Y - RM3(J,I) * RV2(J) - RM3(I,J+2) * RV1(J)
  680       CONTINUE
  700       Z3 = CMPLX(X,Y) / CMPLX(RM3(I,I),RM3(I,I+2))
            RV1(I) = T3(1)
            RV2(I) = T3(2)
  720    CONTINUE
  740    ITS = ITS + 1
         NORM = 0.0
         NORMV = 0.0
         DO 780 I = 1, UK
            IF (ILAMBD .EQ. 0.0) X = ABS(RV1(I))
            IF (ILAMBD .NE. 0.0) X = ABS(CMPLX(RV1(I),RV2(I)))
            IF (NORMV .GE. X) GO TO 760
            NORMV = X
            J = I
  760       NORM = NORM + X
  780    CONTINUE
         IF (NORM .LT. GROWTO) GO TO 840
         X = RV1(J)
         IF (ILAMBD .EQ. 0.0) X = 1.0 / X
         IF (ILAMBD .NE. 0.0) Y = RV2(J)
         DO 820 I = 1, UK
            IF (ILAMBD .NE. 0.0) GO TO 800
            Z(I,S) = RV1(I) * X
            GO TO 820
  800       Z3 = CMPLX(RV1(I),RV2(I)) / CMPLX(X,Y)
            Z(I,S-1) = T3(1)
            Z(I,S) = T3(2)
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
         IF (ILAMBD .EQ. 0.0) GO TO 440
         GO TO 660
  880    J = 1
         IERR = -K
  900    DO 920 I = J, N
            Z(I,S) = 0.0
            IF (ILAMBD .NE. 0.0) Z(I,S-1) = 0.0
  920    CONTINUE
  940    S = S + 1
  960    IF (IP .EQ. (-1)) IP = 0
         IF (IP .EQ. 1) IP = -1
  980 CONTINUE
      GO TO 1001
 1000 IF (IERR .NE. 0) IERR = IERR - N
      IF (IERR .EQ. 0) IERR = -(2 * N + 1)
 1001 M = S - 1 - ABS(IP)
      RETURN
      END
