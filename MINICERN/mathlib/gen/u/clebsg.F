*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:46  mclareni
* Mathlib gen
*
*
*FCA :          Fri Mar 26 17:27:50 CET 1999 by  Federico Carminati
*               removed noarg
#include "gen/pilot.h"
#if (defined(CERNLIB_IBM))&&(defined(CERNLIB_F4))
      FUNCTION CLEBSG(/A/,/B/,/C/,/XX/,/YY/,/ZZ/,/GG/,/HH/,/PP/)
#endif
#if !defined(CERNLIB_IBM)||!defined(CERNLIB_F4)
      FUNCTION CLEBSG( A , B , C , XX , YY , ZZ , GG , HH , PP )
#endif
C
C         ADAPTED FROM HARWELL LIBRARY BY T. LINDELOF AND F. JAMES
C         08/01/74 LAST UPDATE OF HARWELL LIBRARY
C
C     WIGN3J- WIGNER 3-J SYMBOL
C     CLEBSG- CLEBSCH-GORDAN COEFFICIENT
C     WIGN6J- WIGNER 6-J SYMBOL
C     RACAHC- RACAH COEFFICIENT
C     JAHNUF- U-FUNCTION (JAHN)
C     WIGN9J- WIGNER 9-J SYMBOL
C
#if !defined(CERNLIB_F4)
      REAL JAHNUF
#endif
      DIMENSION H(101),J(101)
      DIMENSION AY(4),IAY(4)
      COMMON/FGERCM/IERR,IERCT
      DATA JJJ/0/
      INTPTF(Q)=Q+Q+SIGN(.10,Q)
      IPARF(I)=4*(I/4)-I+1
C*IA  Q
      DATA Q/0/
      IERR= 0
      KEY= 2
*      CALL NOARG(NARG)
      NARG=9
      IF (NARG.EQ.3) KEY=4
      GO TO 1
C
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY WIGN3J
#endif
#if (defined(CERNLIB_IBM))&&(defined(CERNLIB_F4))
      ENTRY WIGN3J(/A/,/B/,/C/,/XX/,/YY/,/ZZ/,/GG/,/HH/,/PP/)
#endif
#if (!defined(CERNLIB_CDC)||!defined(CERNLIB_F4))&&(!defined(CERNLIB_IBM)||!defined(CERNLIB_F4))
      ENTRY WIGN3J(A,B,C,XX,YY,ZZ,GG,HH,PP)
#endif
      IERR=0
*      CALL NOARG(NARG)
      NARG=9
      IF(NARG.EQ.3) GOTO 2
      KEY=1
      GOTO 1
    2 KEY=3
      GOTO 1
C
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY RACAHC
#endif
#if (defined(CERNLIB_IBM))&&(defined(CERNLIB_F4))
      ENTRY RACAHC(/A/,/B/,/C/,/XX/,/YY/,/ZZ/,/GG/,/HH/,/PP/)
#endif
#if (!defined(CERNLIB_CDC)||!defined(CERNLIB_F4))&&(!defined(CERNLIB_IBM)||!defined(CERNLIB_F4))
      ENTRY RACAHC(A,B,C,XX,YY,ZZ,GG,HH,PP)
#endif
      KEY=12
      IERR=0
      GOTO 100
C
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY WIGN6J
#endif
#if (defined(CERNLIB_IBM))&&(defined(CERNLIB_F4))
      ENTRY WIGN6J(/A/,/B/,/C/,/XX/,/YY/,/ZZ/,/GG/,/HH/,/PP/)
#endif
#if (!defined(CERNLIB_CDC)||!defined(CERNLIB_F4))&&(!defined(CERNLIB_IBM)||!defined(CERNLIB_F4))
      ENTRY WIGN6J(A,B,C,XX,YY,ZZ,GG,HH,PP)
#endif
      KEY=11
      IERR=0
      GOTO 100
C
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY JAHNUF
#endif
#if (defined(CERNLIB_IBM))&&(defined(CERNLIB_F4))
      ENTRY JAHNUF(/A/,/B/,/C/,/XX/,/YY/,/ZZ/,/GG/,/HH/,/PP/)
#endif
#if (!defined(CERNLIB_CDC)||!defined(CERNLIB_F4))&&(!defined(CERNLIB_IBM)||!defined(CERNLIB_F4))
      ENTRY JAHNUF(A,B,C,XX,YY,ZZ,GG,HH,PP)
#endif
      KEY=13
      IERR=0
      GOTO 100
C
    1 K1=INTPTF(A)
      K2=INTPTF(B)
      K3=INTPTF(C)
      IF(KEY.GE.3) GOTO 100
      K4=INTPTF(XX)
      K5=INTPTF(YY)
      K6= INTPTF(ZZ)
      IF (KEY.EQ.1) K6= -K6
C
  100 IF(JJJ.NE.0) GOTO 500
      JJJ=1
      IERCT=0
      H(1)=1.0
      J(1)=0
      X=0.
      DO 400 I=2,101
      X=X+1.0
      H(I)=H(I-1)*X
      J(I)=J(I-1)
  200 IF(H(I).LT.10.0) GOTO 400
      H(I)=0.01*H(I)
      J(I)=J(I)+2
      GOTO 200
  400 CONTINUE
C
  500 IF(KEY.LT.-5) GOTO 750
      IF(KEY.GE.3) GOTO 320
      IF((K4+K5-K6).NE.0) GOTO 710
      M1=K1+K2-K3
      M2=K2+K3-K1
      M3=K3+K1-K2
      M4=K1+K4
      M5=K1-K4
      M6=K2+K5
      M7=K2-K5
      M8=K3+K6
      M9=K3-K6
      M10=K1+K2+K3+2
C
      IF(M1.LT.0) GOTO 710
      IF(M2.LT.0) GOTO 710
      IF(M3.LT.0) GOTO 710
      IF(M4.LT.0) GOTO 710
      IF(M5.LT.0) GOTO 710
      IF(M6.LT.0) GOTO 710
      IF(M7.LT.0) GOTO 710
      IF(M8.LT.0) GOTO 710
      IF(M9.LT.0) GOTO 710
      IF((M4-(M4/2)-(M4/2)).NE.0) GOTO 710
      IF((M6-(M6/2)-(M6/2)).NE.0) GOTO 710
      IF((M8-(M8/2)-(M8/2)).NE.0) GOTO 710
      IF((M10-(M10/2)-(M10/2)).NE.0) GOTO 710
C
      Y=K3+1
      M1=M1/2+1
      M2=M2/2+1
      M3=M3/2+1
      M4=M4/2+1
      M5=M5/2+1
      M6=M6/2+1
      M7=M7/2+1
      M8=M8/2+1
      M9=M9/2+1
      M10=M10/2+1
C
      Y= SQRT(Y*H(M1)*H(M2)*H(M3)*H(M4)*H(M5)*
     X H(M6)*H(M7)*H(M8)*H(M9)/H(M10))
      IY=(J(M1)+J(M2)+J(M3)+J(M4)+J(M5)+
     X J(M6)+J(M7)+J(M8)+J(M9)-J(M10))/2
C
      N4=M1
      IF(N4.GT.M5)N4=M5
      IF(N4.GT.M6)N4=M6
      N4=N4-1
      M2=K2-K3-K4
      M3=K1+K5-K3
      N5=0
      IF(N5.LT.M2) N5=M2
      IF(N5.LT.M3) N5=M3
      N5PAR=IPARF(N5)
      N5=N5/2
      Z=0.0
      GOTO 610
C
  700 MM1=M1-N5
      MM2=M5-N5
      MM3=M6-N5
      MM4=N5-(M2/2)+1
      MM5=N5-(M3/2)+1
C
      X=1./(H(MM1)*H(MM2)*H(MM3)*H(MM4)*H(MM5)*H(N5+1))
      IX=-J(MM1)-J(MM2)-J(MM3)-J(MM4)-J(MM5)-J(N5+1)
C
  800 IF(IX+IY)900,210,110
  900 X=0.1*X
      IX=IX+1
      GOTO 800
  110 X=10.0  *X
      IX=IX-1
      GOTO 800
C
  210 IF(N5PAR.LT.0) X=-X
      Z=Z+X
C*UL  510 N5PAR=-N5PAR
      N5PAR=-N5PAR
      N5=N5+1
C
  610 IF(N5-N4)700,700,810
C
 710  CLEBSH=0.0
      IERR=1
      IERCT=IERCT+1
      GOTO 220
C
 810  CLEBSH=Z*Y
C*UL  910 GOTO(120,220),KEY
      GOTO(120,220),KEY
C
  220 CLEBSG=CLEBSH
      RETURN
C
  120 JS=K1-K2+K6
      IF(JS.LT.0) JS=-JS
      JSPAR=IPARF(JS)
      CLEBSG=JSPAR*CLEBSH/ SQRT(K3+1.0  )
      RETURN
C
  320 IF(KEY.GE.10) GOTO 130
      KEY=KEY-2
      IF((K1-(K1/2)-(K1/2)).NE.0) GOTO 420
      IF((K2-(K2/2)-(K2/2)).NE.0) GOTO 420
      IF((K3-(K3/2)-(K3/2)).NE.0) GOTO 420
      IJ=K1+K2+K3
      IJPAR=IPARF(IJ)
      IF(IJPAR.LE.0) GOTO 420
      M1=IJ-K1-K1
      M2=IJ-K2-K2
      M3=IJ-K3-K3
      M4=IJ+2
      IF(M1.LT.0) GOTO 420
      IF(M2.LT.0) GOTO 420
      IF(M3.LT.0) GOTO 420
      M1=M1/2+1
      M2=M2/2+1
      M3=M3/2+1
      M4=IJ/2+2
      Y= SQRT(H(M1)*H(M2)*H(M3)/H(M4))
      IY=(J(M1)+J(M2)+J(M3)-J(M4))/2
      IJ=IJ/2
      IJPAR=IPARF(IJ)
      IJ=IJ/2+1
      M1=M1/2+1
      M2=M2/2+1
      M3=M3/2+1
      Z=H(IJ)/(H(M1)*H(M2)*H(M3))
      IZ=J(IJ)-J(M1)-J(M2)-J(M3)
      IZ=IZ+IY
      CLEBSH=IJPAR*Y*Z*10.0  **IZ
      GOTO(220,720),KEY
C
  720 JQ=K2-K1
      IF(JQ.LT.0) JQ=-JQ
      IJPAR=IPARF(JQ)
      CLEBSG=CLEBSH*IJPAR* SQRT(K3+1.0  )
      RETURN
C
  420 CLEBSH=0.0
      IERR=1
      IERCT=IERCT+1
      GOTO(220,720),KEY
C
  130 IF(KEY.EQ.11) GOTO 450
      IF(KEY.GT.19) GOTO 750
C*UL  550 K1=INTPTF(A)
      K1=INTPTF(A)
      K2=INTPTF(B)
      K3=INTPTF(YY)
      K4=INTPTF(XX)
      K5=INTPTF(C)
      K6=INTPTF(ZZ)
C
  750 KA=K1
      KB=K2
      KC=K3
      KEYTRI=1
      GOTO 630
C
  230 KA=K4
      KB=K5
      KEYTRI=2
      GOTO 630
C
  330 KB=K2
      KC=K6
      KEYTRI=3
      GOTO 630
C
  430 KA=K1
      KB=K5
      KEYTRI=4
      GOTO 630
C
  530 Y=AY(1)*AY(2)*AY(3)*AY(4)
      IYY=IAY(1)+IAY(2)+IAY(3)+IAY(4)
      M1=(K1+K2+K4+K5)/2+2
      M2=(K1+K2-K3)/2+1
      M3=(K4+K5-K3)/2+1
      M4=(K1+K5-K6)/2+1
      M5=(K2+K4-K6)/2+1
      M6=K1+K4-K3-K6
      M7=K2+K5-K3-K6
C
      N4=M1
      IF(N4.GT.M2) N4=M2
      IF(N4.GT.M3) N4=M3
      IF(N4.GT.M4) N4=M4
      IF(N4.GT.M5) N4=M5
      N4=N4-1
      N5=0
      IF(N5.LT.M6) N5=M6
      IF(N5.LT.M7) N5=M7
      N5PAR=IPARF(N5)
      N5=N5/2
      M6=M6/2-1
      M7=M7/2-1
      Z=0.0
      GOTO 730
C
  140 X=H(M1-N5)/(H(N5+1)*H(M2-N5)*H(M3-N5)*H(M4-N5)
     X *H(M5-N5)*H(N5-M6)*H(N5-M7))
      IX=J(M1-N5)-J(N5+1)-J(M2-N5)-J(M3-N5)-J(M4-N5)
     X -J(M5-N5)-J(N5-M6)-J(N5-M7)
  240 IF(IX+IYY)340,440,540
  340 X=0.1*X
      IX=IX+1
      GOTO 240
  540 X=10.0  *X
      IX=IX-1
      GOTO 240
  440 IF(N5PAR.LT.0) X=-X
      Z=Z+X
      N5PAR=-N5PAR
      N5=N5+1
C
  730 IF(N5.LE.N4) GOTO 140
C
      RACAH=Z*Y
  840 IF(KEY.LT.-5) GOTO 160
      KEY=KEY-10
      GOTO(150,250,350),KEY
C
  830 RACAH=0.0
      IERR=1
      IERCT=IERCT+1
      GOTO 840
C
  150 IJPAR=IPARF(K1+K2+K4+K5)
      IF(IJPAR.LT.0) RACAH=-RACAH
  250 CLEBSG=RACAH
      RETURN
C
  350 FACTOR= SQRT((K3+1.0  )*(K6+1))
      CLEBSG=FACTOR*RACAH
      RETURN
  450 K1=INTPTF(A)
      K2=INTPTF(B)
      K3=INTPTF(C)
      K4=INTPTF(XX)
      K5=INTPTF(YY)
      K6=INTPTF(ZZ)
      GOTO 750
C
C     TRIANGLE FUNCTION
C
  630 MA=KA+KB-KC
      MB=KA-KB+KC
      MC=-KA+KB+KC
      MD=KA+KB+KC+2
      IF(MA.LT.0) GOTO 830
      IF(MB.LT.0) GOTO 830
      IF(MC.LT.0) GOTO 830
      IF((MD-(MD/2)-(MD/2)).NE.0) GOTO 830
      MA=MA/2+1
      MB=MB/2+1
      MC=MC/2+1
      MD=MD/2+1
      AY(KEYTRI)= SQRT(H(MA)*H(MB)*H(MC)/H(MD))
      IAY(KEYTRI)=(J(MA)+J(MB)+J(MC)-J(MD))/2
      GOTO(230,330,430,530),KEYTRI
C
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY WIGN9J
#endif
#if (defined(CERNLIB_IBM))&&(defined(CERNLIB_F4))
      ENTRY WIGN9J(/A/,/B/,/C/,/XX/,/YY/,/ZZ/,/GG/,/HH/,/PP/)
#endif
#if (!defined(CERNLIB_CDC)||!defined(CERNLIB_F4))&&(!defined(CERNLIB_IBM)||!defined(CERNLIB_F4))
      ENTRY WIGN9J(A,B,C,XX,YY,ZZ,GG,HH,PP)
#endif
C
      KEY=-10
      IERR=0
C
      KK1=INTPTF(A)
      KK2=INTPTF(B)
      KK3=INTPTF(C)
      KK4=INTPTF(XX)
      KK5=INTPTF(YY)
      KK6=INTPTF(ZZ)
      KK7=INTPTF(GG)
      KK8=INTPTF(HH)
      KK9=INTPTF(PP)
C
      KUP=KK1+KK9
      M1=KK4+KK8
      M2=KK2+KK6
      IF(KUP.GT.M1) KUP=M1
      IF(KUP.GT.M2) KUP=M2
C
      K=KK1-KK9
      IF(K.LT.0) K=-K
      M1=KK4-KK8
      IF(M1.LT.0) M1=-M1
      M2=KK2-KK6
      IF(M2.LT.0) M2=-M2
      IF(K.LT.M1) K=M1
      IF(K.LT.M2) K=M2
C
      ANINE=0.0
C
  660 IF(K.GT.KUP) GOTO 260
      K1=KK1
      K2=KK4
      K3=KK7
      K4=KK8
      K5=KK9
      K6=K
      KEYRAC=1
      GOTO 100
C
  160 GOTO(360,460,560),KEYRAC
C
  360 RA=RACAH
      K1=KK2
      K2=KK8
      K3=KK5
      K4=KK4
      K5=KK6
      KEYRAC=2
      GOTO 750
C
  460 RB=RACAH
      K1=KK9
      K2=KK6
      K3=KK3
      K4=KK2
      K5=KK1
      KEYRAC=3
      GOTO 750
C
  560 ANINE=ANINE+RA*RB*RACAH*(K+1)
      K=K+2
      GOTO 660
C
  260 CLEBSG=ANINE
      RETURN
      END
