*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:32  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      FUNCTION GPRELM(Z,T,CUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       To calculate DE/DX in GeV*barn/atom for direct           *
C.    *       pair production by muons.                                *
C.    *                                                                *
C.    *    ==>Called by : GPRSGA                                       *
C.    *       Author    L.Urban  *********                             *
C.    *                                                                *
C.    ******************************************************************
C.
      PARAMETER (AKSI=1.19,BETA=1.24,DEL=0.13,VE=0.008)
      PARAMETER (ECMIN=2.044E-3,CONMX=0.130653)
C                               CONMX=0.75*SQRT(2.7182...)*EMMU
#include "geant321/gconsp.inc"
      DIMENSION C(100),D(28),C1(60),C2(40)
      EQUIVALENCE (C(1),C1(1)),(C(61),C2(1))
      DATA C1/0.199979E-04,-0.247381E-04, 0.901012E-05,-0.624746E-06
     +      , 0.306301E-08, 0.568134E-09,-0.154371E-04, 0.192663E-04
     +      ,-0.638613E-05, 0.325587E-06, 0.798057E-08,-0.623709E-09
     +      , 0.334831E-04,-0.305301E-04, 0.514764E-05, 0.273273E-07
     +      ,-0.455952E-07, 0.200990E-08,-0.617418E-05,-0.119758E-04
     +      , 0.505842E-05,-0.680982E-06, 0.335650E-07,-0.465426E-09
     +      ,-0.134652E-04, 0.268825E-05, 0.722810E-07,-0.648440E-07
     +      , 0.532560E-08,-0.122207E-09,-0.354308E-05, 0.125249E-05
     +      ,-0.182348E-06, 0.125659E-07,-0.390005E-09, 0.423919E-11
     +      , 0.427113E-05,-0.570105E-05, 0.156413E-05,-0.247880E-07
     +      ,-0.534990E-08, 0.172881E-09,-0.309866E-06, 0.654607E-06
     +      ,-0.885876E-07,-0.148160E-07, 0.814881E-09, 0.144373E-10
     +      ,-0.594089E-07, 0.643470E-08,-0.204298E-07, 0.346177E-08
     +      , 0.451711E-10,-0.115400E-10, 0.849857E-08,-0.485366E-08
     +      , 0.247710E-08,-0.224799E-09,-0.125246E-10, 0.116449E-11/
      DATA C2/0.406000E-08,-0.164080E-07, 0.109450E-07,-0.201483E-08
     +      , 0.889711E-10, 0.583750E-08, 0.523552E-08,-0.713086E-08
     +      , 0.151684E-08,-0.703329E-10, 0.863004E-07,-0.107717E-06
     +      , 0.344991E-07,-0.382381E-08, 0.128222E-09,-0.659685E-07
     +      , 0.283383E-07,-0.263676E-08,-0.978698E-10, 0.109724E-10
     +      , 0.209150E-09, 0.880502E-09,-0.235454E-09, 0.120280E-10
     +      ,-0.120162E-13,-0.814192E-08, 0.402793E-08, 0.952674E-09
     +      ,-0.377756E-09, 0.198358E-10, 0.254519E-08,-0.187283E-08
     +      , 0.168056E-09, 0.293299E-10,-0.221041E-11,-0.233497E-09
     +      , 0.197097E-09,-0.321342E-10, 0.329587E-12, 0.691294E-13/
      DATA D/-0.790941E-07, 0.840429E-07, 0.738033E-08,-0.407627E-07
     +      ,-0.116392E-07, 0.358398E-08, 0.141013E-07, 0.814070E-08
     +      ,-0.423749E-08, 0.173692E-08,-0.174475E-08,-0.185530E-08
     +      , 0.737768E-09, 0.160269E-09,-0.598384E-09, 0.926201E-10
     +      , 0.133740E-09, 0.957004E-11,-0.850466E-10, 0.214267E-10
     +      , 0.967441E-10,-0.181431E-11,-0.324662E-11,-0.579085E-12
     +      ,-0.237465E-11, 0.108081E-10,-0.902849E-11,-0.408422E-11/
C.
C.    ------------------------------------------------------------------
C.
      GPRELM=0.
      IF(CUT.LT.ECMIN) GOTO 999
C
      Z3=Z**0.333333
      E=T+EMMU
      ECMAX=E-CONMX*Z3
      ECMA5=ECMAX*0.2
      X=LOG(E/EMMU)
      IF(CUT.GE.ECMAX) GOTO 160
C
      IF(CUT.LE.ECMA5) THEN
        CCUT=CUT
      ELSE
        CCUT=ECMA5
      ENDIF
      Y=LOG(CCUT/(VE*E))
C
      S=0.
      YY=1.
      DO 30 I=1,2
         XX=1.
         DO 20 J=1,6
            K=6*I+J-6
            S=S+C(K)*XX*YY
            XX=XX*X
  20     CONTINUE
         YY=YY*Y
  30  CONTINUE
      DO 50 I=3,6
         XX=1.
         DO 40 J=1,6
            K=6*I+J-6
            IF(Y.LE.0.) THEN
               S=S+C(K)*XX*YY
            ELSE
               S=S+C(K+24)*XX*YY
            ENDIF
            XX=XX*X
  40     CONTINUE
         YY=YY*Y
  50  CONTINUE
      SS=0.
      YY=1.
      DO 70 I=1,2
         XX=1.
         DO 60 J=1,5
            K=5*I+J+55
            SS=SS+C(K)*XX*YY
            XX=XX*X
  60     CONTINUE
         YY=YY*Y
  70  CONTINUE
      DO 90 I=3,5
         XX=1.
         DO 80 J=1,5
            K=5*I+J+55
            IF(Y.LE.0.) THEN
               SS=SS+C(K)*XX*YY
            ELSE
               SS=SS+C(K+15)*XX*YY
            ENDIF
            XX=XX*X
  80     CONTINUE
         YY=YY*Y
  90  CONTINUE
C
      S=S+Z*SS
      IF(S.LE.0.) GOTO 999
C
      FAC=E*((CCUT-ECMIN)/E)**BETA
      FAC=Z*(Z+AKSI*(1.+DEL*LOG(Z)))*FAC
      GPRELM=FAC*S
      IF(CUT.LE.ECMA5) GOTO 999
      GPREL5=GPRELM
C
 160  Y=Z3
      S=D(1)
      K=1
      XX=1.
      YPX=Y/X
      DO 180 I=1,6
         II=I+1
         XX=XX*X
         XXYY=XX
         DO 170 J=1,II
            K=K+1
            S=S+D(K)*XXYY
            XXYY=XXYY*YPX
 170     CONTINUE
 180  CONTINUE
      IF(S.LE.0.) GOTO 999
C
      FAC=Z*(Z+1.)*ECMAX
      IF(FAC.LE.0.) GOTO 999
      GPRELM=FAC*S
C
      IF((CUT.LE.ECMA5).OR.(CUT.GE.ECMAX)) GOTO 999
      DEL1=(GPRELM-GPREL5)/(ECMAX-ECMA5)
      GPRELM=GPREL5+DEL1*(CUT-ECMA5)
C
999   RETURN
      END
