*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:33  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE COMLR2(NM,N,LOW,IGH,INT,HR,HI,WR,WI,ZR,ZI,IERR)
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NM,NN,IGH,IM1,IP1,
     X        ITS,LOW,MP1,ENM1,IEND,IERR
      REAL HR(NM,N),HI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N)
      REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,MACHEP
      INTEGER INT(IGH)
      COMPLEX X,Y,Z
      REAL T1(2),T2(2),T3(2)
      EQUIVALENCE (X,T1(1),XR),(T1(2),XI),(Y,T2(1),YR),(T2(2),YI),
     X            (Z,T3(1),ZZR),(T3(2),ZZI)
#if defined(CERNLIB_CDC)
      MACHEP=2.**(-47)
#endif
#if !defined(CERNLIB_CDC)
      MACHEP=2.**(-23)
#endif
      IERR = 0
      DO 100 I = 1,N
         DO 100 J = 1, N
            ZR(I,J) = 0.0
            ZI(I,J) = 0.0
            IF (I .EQ. J) ZR(I,J) = 1.0
  100 CONTINUE
      IEND = IGH - LOW - 1
      IF (IEND .LE. 0) GO TO 180
      DO 160 II = 1, IEND
         I = IGH - II
         IP1 = I + 1
         DO 120 K = IP1, IGH
            ZR(K,I) = HR(K,I-1)
            ZI(K,I) = HI(K,I-1)
  120    CONTINUE
         J = INT(I)
         IF (I .EQ. J) GO TO 160
         DO 140 K = I, IGH
            ZR(I,K) = ZR(J,K)
            ZI(I,K) = ZI(J,K)
            ZR(J,K) = 0.0
            ZI(J,K) = 0.0
  140    CONTINUE
         ZR(J,I) = 1.0
  160 CONTINUE
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
      EN = IGH
      TR = 0.0
      TI = 0.0
  220 IF (EN .LT. LOW) GO TO 680
      ITS = 0
      ENM1 = EN - 1
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         IF (ABS(HR(L,L-1)) + ABS(HI(L,L-1)) .LE.
     X      MACHEP * (ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1))
     X             + ABS(HR(L,L)) +ABS(HI(L,L)))) GO TO 300
  260 CONTINUE
  300 IF (L .EQ. EN) GO TO 660
      IF (ITS .EQ. 30) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1) - HI(ENM1,EN) * HI(EN,ENM1)
      XI = HR(ENM1,EN) * HI(EN,ENM1) + HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0 .AND. XI .EQ. 0.0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0
      YI = (HI(ENM1,ENM1) - SI) / 2.0
      Z = SQRT(CMPLX(YR**2-YI**2+XR,2.0*YR*YI+XI))
      IF (YR * ZZR + YI * ZZI .LT. 0.0) Z = -Z
      X = X / (Y + Z)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = ABS(HI(EN,ENM1)) + ABS(HI(ENM1,EN-2))
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      XR = ABS(HR(ENM1,ENM1)) + ABS(HI(ENM1,ENM1))
      YR = ABS(HR(EN,ENM1)) + ABS(HI(EN,ENM1))
      ZZR = ABS(HR(EN,EN)) + ABS(HI(EN,EN))
      DO 380 MM = L, ENM1
         M = ENM1 + L - MM
         IF (M .EQ. L) GO TO 420
         YI = YR
         YR = ABS(HR(M,M-1)) + ABS(HI(M,M-1))
         XI = ZZR
         ZZR = XR
         XR = ABS(HR(M-1,M-1)) + ABS(HI(M-1,M-1))
         IF (YR .LE. MACHEP * ZZR / YI * (ZZR + XR + XI)) GO TO 420
  380 CONTINUE
  420 MP1 = M + 1
      DO 520 I = MP1, EN
         IM1 = I - 1
         XR = HR(IM1,IM1)
         XI = HI(IM1,IM1)
         YR = HR(I,IM1)
         YI = HI(I,IM1)
         IF (ABS(XR) + ABS(XI) .GE. ABS(YR) + ABS(YI)) GO TO 460
         DO 440 J = IM1, N
            ZZR = HR(IM1,J)
            HR(IM1,J) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(IM1,J)
            HI(IM1,J) = HI(I,J)
            HI(I,J) = ZZI
  440    CONTINUE
         Z = X / Y
         WR(I) = 1.0
         GO TO 480
  460    Z = Y / X
         WR(I) = -1.0
  480    HR(I,IM1) = ZZR
         HI(I,IM1) = ZZI
         DO 500 J = I, N
            HR(I,J) = HR(I,J) - ZZR * HR(IM1,J) + ZZI * HI(IM1,J)
            HI(I,J) = HI(I,J) - ZZR * HI(IM1,J) - ZZI * HR(IM1,J)
  500   CONTINUE
  520 CONTINUE
      DO 640 J = MP1, EN
         XR = HR(J,J-1)
         XI = HI(J,J-1)
         HR(J,J-1) = 0.0
         HI(J,J-1) = 0.0
         IF (WR(J) .LE. 0.0) GO TO 580
         DO 540 I = 1, J
            ZZR = HR(I,J-1)
            HR(I,J-1) = HR(I,J)
            HR(I,J) = ZZR
            ZZI = HI(I,J-1)
            HI(I,J-1) = HI(I,J)
            HI(I,J) = ZZI
  540    CONTINUE
         DO 560 I = LOW, IGH
            ZZR = ZR(I,J-1)
            ZR(I,J-1) = ZR(I,J)
            ZR(I,J) = ZZR
            ZZI = ZI(I,J-1)
            ZI(I,J-1) = ZI(I,J)
            ZI(I,J) = ZZI
  560    CONTINUE
  580    DO 600 I = 1, J
            HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)
            HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)
  600    CONTINUE
         DO 620 I = LOW, IGH
            ZR(I,J-1) = ZR(I,J-1) + XR * ZR(I,J) - XI * ZI(I,J)
            ZI(I,J-1) = ZI(I,J-1) + XR * ZI(I,J) + XI * ZR(I,J)
  620    CONTINUE
  640 CONTINUE
      GO TO 240
  660 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GO TO 220
  680 IF (N .EQ. 1) GO TO 1001
      NORM = 0.0
      DO 720 I = 1, N
         DO 720 J = I, N
            NORM = NORM + ABS(HR(I,J)) + ABS(HI(I,J))
  720 CONTINUE
      DO 800 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         ENM1 = EN - 1
         DO 780 II = 1, ENM1
            I = EN - II
            ZZR = HR(I,EN)
            ZZI = HI(I,EN)
            IF (I .EQ. ENM1) GO TO 760
            IP1 = I + 1
            DO 740 J = IP1, ENM1
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
  760       YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR .EQ. 0.0 .AND. YI .EQ. 0.0) YR = MACHEP * NORM
            Z = Z / Y
            HR(I,EN) = T3(1)
            HI(I,EN) = T3(2)
  780    CONTINUE
  800 CONTINUE
      ENM1 = N - 1
      DO  840 I = 1, ENM1
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO  840
         IP1 = I + 1
         DO 820 J = IP1, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  820    CONTINUE
  840 CONTINUE
      DO 880 JJ = LOW, ENM1
         J = N + LOW - JJ
         M = MIN(J-1,IGH)
         DO 880 I = LOW, IGH
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
            DO 860 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  880 CONTINUE
      GO TO 1001
 1000 IERR = EN
 1001 RETURN
      END
