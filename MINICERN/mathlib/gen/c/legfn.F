*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:59  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE LEGFN (V,Z,P,Q,NC,NF)
#include "legbl.inc"
     1,V,Z,P,Q
     2,ZZS,CISP,CISM,PT,VVP,ZZZ
     3,CSQRTK
      DATA AC/0.0000001/
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      DATA QINF/37770000000000000000B/
#endif
#if (defined(CERNLIB_CDC))&&(!defined(CERNLIB_F4))
      DATA QINF/O"37770000000000000000"/
#endif
#if defined(CERNLIB_IBM)
      DATA QINF/.1E70/
#endif
#if (!defined(CERNLIB_CDC))&&(!defined(CERNLIB_IBM))
      DATA QINF / .1E36/
#endif
      CALL C311BD
      ACC=AC**2
      IF(REAL(V)+0.5) 16,17,17
   16 VV=-V-1.0
      GOTO 18
   17 VV=V
C---  It seems that NIC311 expects a REAL argument
C     NVV=NIC311(VV)
   18 NVV=NIC311(REAL(VV))
C*UL 30      VVP=VV*PI
      VVP=VV*PI
      GOTO (21,22,23,24,25),NVV
   21 SV=1.0
      GOTO 26
   23 SV=-1.0
   26 CV=0.0
      GOTO 28
   22 CV=-1.0
      GOTO 27
   24 CV=1.0
   27 SV=0.0
      GOTO 28
   25 SV=SIN(REAL(VVP))
      CV=COS(REAL(VVP))
   28 EIP=EXP(AIMAG(VVP))
      EIM=EXP(-AIMAG(VVP))
      SHV=0.5*(EIP-EIM)
      CHV=0.5*(EIP+EIM)
      SVV=SV*CHV+U*CV*SHV
      CVV=CV*CHV-U*SV*SHV
      CISP=CVV+U*SVV
      CISM=CVV-U*SVV
      IF(REAL(Z))9,10,11
    9 ZZ=-Z
      N23=3
      NFF=-NF
      GOTO 12
   10 N24=4
      ZZ=Z
      IF(AIMAG(ZZ))13,77,13
   11 ZZ=Z
      N23=2
      NFF=NF
   12 N24=2
   13 ZZS=ZZ**2
      Z1=(1.0-ZZ)/2.0
      IF(AIMAG(ZZ))7,6,7
    6 IF(REAL(Z1))4,70,5
    4 NFRIG=SIGN(N23,NFF)
      GOTO 8
    5 NFRIG=NFF
      GOTO 8
    7 NFRIG=SIGN(REAL(N24),AIMAG(ZZ))
    8 Z2=1.0/ZZS
      SRZ=CSQRTK(ZZS-1.0,NFRIG,1)
      ZZ1=(ZZ+SRZ)/(2.0*SRZ)
      ZZ2=(-ZZ+SRZ)/(2.0*SRZ)
      IF(REAL(Z))1,2,2
    1 SRZ=-SRZ
    2 ZZZ=ZZ
      ZZ=Z
      VR=REAL(VV)**2
      VI=AIMAG(VV)**2
      R1=ABS(Z1)
      R2=ABS(Z2)
      R3=ABS(ZZ1)
      R4=ABS(ZZ2)
      IF(R1-R2)61,61,62
   61 RR=MAX(R3,R4)/0.045-19.5
      IF(RR)66,66,65
   65 VRI=VR+VI
      IF(VRI)165,67,165
  165 IF(RR**2-VRI+VI**2/(2.0*VRI))66,67,67
   66 CALL LEGV
      IF(NCVG)166,97,166
  166 ZZ=ZZZ
      CALL LEG1
      GOTO 90
   67 ZZ=ZZZ
      CALL LEG1
      IF(NCVG)167,90,167
  167 ZZ=Z
      CALL LEGV
      GOTO 97
   62 IF(VR+VI-16.0*R2**2)63,64,64
   64 CALL LEGV
      IF(NCVG)164,97,164
  164 CALL LEGZ
      GOTO 97
   63 CALL LEGZ
      IF(NCVG)163,97,163
  163 CALL LEGV
      GOTO 97
   90 IF(REAL(Z)) 93,97,97
   93 PT=-2.0/PI*SVV*QQ
      IF(NFRIG) 94,95,96
   94 PP=CISP*PP+PT
      QQ=-CISM*QQ
      GOTO 97
   95 QQ=-CVV*QQ-PI/2.0*SVV*PP
      PP=CVV*PP+PT
      GOTO 97
   96 PP=CISM*PP+PT
      QQ=-CISP*QQ
      GOTO 97
   97 IF(REAL(V)+0.5) 91,92,92
   91 IF(ABS(SVV).NE.0.) GOTO 98
      QQ=QINF
      GOTO 92
   98 QQ=(QQ*SVV-PI*CVV*PP)/SVV
   92 P=PP
      Q=QQ
      NC=NCVG
      RETURN
   70 QQ=QINF
C
      IF(REAL(Z))173,173,74
  173 IF(AIMAG(V))71,73,71
   73 GOTO (71,72,71,74,71),NVV
   71 PP=QQ
      GOTO 82
   72 PP=-1.0
      GOTO 82
   74 PP=1.0
      GOTO 82
   77 NFRIG=NF
      CALL LEGOR
   82 NCVG=0
      GOTO 92
      END
