#include "isajet/pilot.h"
      SUBROUTINE SIGSSZ
C
C          Calculate d(sigma)/d(pt**2)d(y1)d(y2) for supersymmetric
C          zino or wino plus squark or gluino in MSSM using cross
C          sections from Baer, Karatas, and Tata, PR D42, 2259.
C          Also include wino and zino pairs.
C
C          SIGMA    = cross section summed over types allowed by
C                     JETTYPE cards.
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 +I1
C          JETTYP -> IDENT mapping:
C          GLSS, UPSSL, UBSSL, ..., UPSSR, UBSSR, ...,
C          W1SS+, W1SS-, WS22+, W2SS-, Z1SS, Z2SS, Z3SS, Z4SS
C
C          Extra factor of 1/2 needed for nonidentical final jets.
C          Y=-log(tan(theta/2)) gives jacobean P1*P2/E1*E2
C
C          Called from SIGSSY and so does not reinitialize /JETSIG/.
C
C          Ver 7.23: Add test setting SIG=0 for Z_i pairs if 
C          ABS(ZZ)>0.999 and SIG<0.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/const.inc"
#include "isajet/jetpar.inc"
#include "isajet/jetsig.inc"
#include "isajet/primar.inc"
#include "isajet/q1q2.inc"
#include "isajet/qcdpar.inc"
#include "isajet/sspar.inc"
#include "isajet/sssm.inc"
#include "isajet/sstype.inc"
#include "isajet/wcon.inc"
C
      REAL X(2)
      EQUIVALENCE (X(1),X1)
      COMPLEX AQZ(2,4),BQZ(2,4),AQW(2,2),WIJ
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT)
      INTEGER JS2JT(25),IW2JS(4),IW2IM(4),IZ2JS(4),IS2UD(25)
      SAVE JS2JT,IW2JS,IW2IM,IZ2JS,IS2UD
      INTEGER IDQSS(25),IDZSS(4),IDWSS(4)
      SAVE IDQSS,IDZSS,IDWSS
      INTEGER ITHZ(4),ITHW(2)
      REAL AMWISS(2)
      REAL XZIWJ(4,2),YZIWJ(4,2)
      REAL SIG,SIG0,CON,AMQIQ,S,T,U,AMWIW,FAC,AM22,AM12,TT,GP,G,
     $E1,E2,AMG,YM,XM,GS,THX,THY,AMZIZ,AMSQK
      INTEGER IX,JQ,IQ,IQ1,IQ2,JW,IW,JTYPW,IH,JTYPZ,IZ,ITHG,IWM
      COMPLEX ZONE,ZI
      SAVE ZONE,ZI
      REAL QFCN,STRUC,PSIFCN,AMASS
      REAL CON11,CON22,CON12,AMQIQ1,AMQIQ2
      INTEGER IX1,IX2
      REAL CS2THW,TNTHW,CTTHW,AL(2),BE(2),ESQ,XWI(2),YWI(2)
      REAL X12,Y12,SN12,AMWIW1,AMWIW2,EQ1,ZZ,XMGG,XMZZ
      REAL XMGZ,XMUU,XMGU,XMZU,XMDD,XMGD,XMZD,DEL,RSH,SR2
      REAL SIGUT,SIGTU,EHAT,PHAT,EBM,TPP,AMWI,AMQ,PROPW
      REAL SIGUT1,SIGUT2,SIGUT3,SGUT12,SGUT13,SGUT23
      REAL SIGTU1,SIGTU2,SIGTU3,SGTU12,SGTU13,SGTU23
      REAL AMSQL,AMSQR,KK,AMZIZ1,AMZIZ2
      REAL SIGLL,SIGRR,SIGZZ,SIGLZ,SIGRZ,SSGT,SSGST,PROPZ,SSXLAM
      INTEGER IZ1,JTYPZ1,IZ2,JTYPZ2
      INTEGER IW1,JW1,JTYPW1,IDW1,IW2,JW2,JTYPW2,IDW2,IFLQ,IUD(13)
C
C          IDENT codes from /SSTYPE/. (Fortran 77 allows - signs in
C          parameter statements but not data statements.)
      INTEGER MSUPL,MSDNL,MSSTL,MSCHL,MSBT1,MSTP1,
     $MSUPR,MSDNR,MSSTR,MSCHR,MSBT2,MSTP2,MSW1,MSW2
      PARAMETER (MSUPL=-ISUPL)
      PARAMETER (MSDNL=-ISDNL)
      PARAMETER (MSSTL=-ISSTL)
      PARAMETER (MSCHL=-ISCHL)
      PARAMETER (MSBT1=-ISBT1)
      PARAMETER (MSTP1=-ISTP1)
      PARAMETER (MSUPR=-ISUPR)
      PARAMETER (MSDNR=-ISDNR)
      PARAMETER (MSSTR=-ISSTR)
      PARAMETER (MSCHR=-ISCHR)
      PARAMETER (MSBT2=-ISBT2)
      PARAMETER (MSTP2=-ISTP2)
      PARAMETER (MSW1=-ISW1)
      PARAMETER (MSW2=-ISW2)
      DATA IDQSS/0,
     $ISUPL,MSUPL,ISDNL,MSDNL,ISSTL,MSSTL,ISCHL,MSCHL,ISBT1,MSBT1,
     $ISTP1,MSTP1,
     $ISUPR,MSUPR,ISDNR,MSDNR,ISSTR,MSSTR,ISCHR,MSCHR,ISBT2,MSBT2,
     $ISTP2,MSTP2/
      DATA IDZSS/ISZ1,ISZ2,ISZ3,ISZ4/
      DATA IDWSS/ISW1,MSW1,ISW2,MSW2/
      DATA IUD/0,1,-1,2,-2,2,-2,1,-1,2,-2,1,-1/
C
C          JS2JT: Susy jettype -> normal jettype
      DATA JS2JT/1,
     $2,3,4,5,6,7,8,9,10,11,12,13,2,3,4,5,6,7,8,9,10,11,12,13/
C          IW2JS: Wino index -> susy jettype
      DATA IW2JS/26,27,28,29/
C          IW2IM: Wino index -> match code
      DATA IW2IM/2,3,2,3/
C          IZ2JS: Zino index -> susy jettype
      DATA IZ2JS/30,31,32,33/
C          IS2UD: Susy jettype -> u/d code
      DATA IS2UD/0,1,1,2,2,2,2,1,1,2,2,1,1,1,1,2,2,2,2,1,1,2,2,1,1/
C
      DATA ZONE,ZI/(1.,0.),(0.,1.)/
C
C          Functions
      QFCN(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
      PSIFCN(AM12,AM22,TT)=((S+TT-AM12)/(2*S)
     $-AM12*(AM22-TT)/(AM12-TT)**2
     $+(TT*(AM22-AM12)+AM22*(S-AM22+AM12))/(S*(AM12-TT)))
C
C          Constants from Baer, Barger, Karatas, and Tata,
C          PR D36, 96, using results from SSMIX
C
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
C     GS=SQRT(4.*PI*ALFA3)
      XM=1./TAN(GAMMAL)
      YM=1./TAN(GAMMAR)
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
      AMG=AMASS(ISGL)
      ITHG=+1
C          Signed masses
      AMWISS(1)=AMW1SS
      AMWISS(2)=AMW2SS
C          Zi couplings
      DO 100 IZ=1,4
        ITHZ(IZ)=0
        IF(AMZISS(IZ).LT.0) ITHZ(IZ)=1
        AQZ(1,IZ)=ZI**(ITHZ(IZ)-1)*(-ZONE)**(ITHZ(IZ)+1)
     $  *(+G/SQRT2*ZMIXSS(3,IZ)+GP/(3*SQRT2)*ZMIXSS(4,IZ))
        AQZ(2,IZ)=ZI**(ITHZ(IZ)-1)*(-ZONE)**(ITHZ(IZ)+1)
     $  *(-G/SQRT2*ZMIXSS(3,IZ)+GP/(3*SQRT2)*ZMIXSS(4,IZ))
        BQZ(1,IZ)=+(4./3.)*ZI**(ITHZ(IZ)-1)*GP/SQRT2*ZMIXSS(4,IZ)
        BQZ(2,IZ)=-(2./3.)*ZI**(ITHZ(IZ)-1)*GP/SQRT2*ZMIXSS(4,IZ)
100   CONTINUE
C          Wi couplings
      ITHW(1)=0
      IF(AMW1SS.LT.0.) ITHW(1)=1
      AQW(1,1)=ZI*G*SIN(GAMMAL)
      AQW(2,1)=ZI*G*(-ZONE)**ITHW(1)*SIN(GAMMAR)
      ITHW(2)=0
      IF(AMW2SS.LT.0.) ITHW(2)=1
      AQW(1,2)=ZI*G*THX*COS(GAMMAL)
      AQW(2,2)=ZI*G*(-ZONE)**ITHW(2)*THY*COS(GAMMAR)
C          Quark couplings to Z
      CS2THW=1.-SN2THW
      TNTHW=SQRT(SN2THW/CS2THW)
      CTTHW=1./TNTHW
      AL(1)=CTTHW/4.-5*TNTHW/12.
      AL(2)=TNTHW/12.-CTTHW/4.
      BE(1)=-(CTTHW+TNTHW)/4.
      BE(2)=-BE(1)
      ESQ=4*PI*ALFAEM
C           Chargino couplings to Z
      XWI(1)=1.-(COS(GAMMAL)**2+COS(GAMMAR)**2)/4./CS2THW
      XWI(2)=1.-(SIN(GAMMAL)**2+SIN(GAMMAR)**2)/4./CS2THW
      YWI(1)=(COS(GAMMAR)**2-COS(GAMMAL)**2)/4./CS2THW
      YWI(2)=(SIN(GAMMAR)**2-SIN(GAMMAL)**2)/4./CS2THW
      X12=.5*(THX*SIN(GAMMAL)*COS(GAMMAL)-
     $    THY*SIN(GAMMAR)*COS(GAMMAR))
      Y12=.5*(THX*SIN(GAMMAL)*COS(GAMMAL)+
     $    THY*SIN(GAMMAR)*COS(GAMMAR))
      SN12=-1.*SIGN(1.,AMW1SS)*SIGN(1.,AMW2SS)
C
C         qk qb --> ziss glss
C
      DO 200 IZ=1,4
        AMZIZ=ABS(AMZISS(IZ))
        JTYPZ=IZ2JS(IZ)
C          Jet 1 = ziss, jet 2 = glss
        IF(.NOT.(GOQ(JTYPZ,1).AND.GOQ(1,2))) GO TO 220
        CALL TWOKIN(0.,0.,AMZIZ,AMG)
        IF(X1.GE.1..OR.X2.GE.1.) GO TO 220
        GS=SQRT(4*PI*ALFQSQ)
        E1=SQRT(P(1)**2+AMZIZ**2)
        E2=SQRT(P(2)**2+AMG**2)
        FAC=1./(16.*PI*S**2)
        FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
C          Sum over initial quarks (no top quarks)
        DO 210 IQ=2,11
          IQ1=IQ
          IQ2=MATCH(IQ1,4)
          AMQIQ=AMASS(IDQSS(IQ))
          SIG0=(AMZIZ**2-T)*(AMG**2-T)/(AMQIQ**2-T)**2
     $    +(AMZIZ**2-U)*(AMG**2-U)/(AMQIQ**2-U)**2
     $    -2*(-1)**(ITHZ(IZ)+ITHG)*AMG*AMZIZ*S
     $    /((AMQIQ**2-T)*(AMQIQ**2-U))
          SIG0=SIG0*2*GS**2/9
          CON=AQZ(IS2UD(IQ),IZ)*CONJG(AQZ(IS2UD(IQ),IZ))
     $    +BQZ(IS2UD(IQ),IZ)*CONJG(BQZ(IS2UD(IQ),IZ))
          SIG=FAC*CON*SIG0*QFCN(IQ1,1)*QFCN(IQ2,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,IQ1,IQ2,JTYPZ,1)
210     CONTINUE
C          Jet 1 = glss, jet 2 = ziss
220     IF(.NOT.(GOQ(1,1).AND.GOQ(JTYPZ,2))) GO TO 200
        CALL TWOKIN(0.,0.,AMG,AMZIZ)
        IF(X1.GE.1..OR.X2.GE.1.) GO TO 200
        GS=SQRT(4*PI*ALFQSQ)
        E1=SQRT(P(1)**2+AMG**2)
        E2=SQRT(P(2)**2+AMZIZ**2)
        FAC=1./(16.*PI*S**2)
        FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
        DO 230 IQ=2,11
          IQ1=IQ
          IQ2=MATCH(IQ1,4)
          AMQIQ=AMASS(IDQSS(IQ))
          SIG0=(AMZIZ**2-T)*(AMG**2-T)/(AMQIQ**2-T)**2
     $    +(AMZIZ**2-U)*(AMG**2-U)/(AMQIQ**2-U)**2
     $    -2*(-1)**(ITHZ(IZ)+ITHG)*AMG*AMZIZ*S
     $    /((AMQIQ**2-T)*(AMQIQ**2-U))
          SIG0=SIG0*2*GS**2/9
          CON=AQZ(IS2UD(IQ),IZ)*CONJG(AQZ(IS2UD(IQ),IZ))
     $    +BQZ(IS2UD(IQ),IZ)*CONJG(BQZ(IS2UD(IQ),IZ))
          SIG=FAC*CON*SIG0*QFCN(IQ1,1)*QFCN(IQ2,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,IQ1,IQ2,1,JTYPZ)
230     CONTINUE
200   CONTINUE
C
C          qk gl -> ziss qkss 
C
      DO 300 IZ=1,4
        AMZIZ=ABS(AMZISS(IZ))
        JTYPZ=IZ2JS(IZ)
        DO 310 IQ=2,25
          JQ=JS2JT(IQ)
          IF(IABS(JQ).GE.12) GO TO 310
          AMQIQ=AMASS(IDQSS(IQ))
C          Jet 1 = ziss, jet 2 = qkss
          IF(.NOT.(GOQ(JTYPZ,1).AND.GOQ(IQ,2))) GO TO 320
          CALL TWOKIN(0.,0.,AMZIZ,AMQIQ)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 320
          GS=SQRT(4*PI*ALFQSQ)
          E1=SQRT(P(1)**2+AMZIZ**2)
          E2=SQRT(P(2)**2+AMQIQ**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          IX=IS2UD(IQ)
C          Use AQZ for left squarks, BQZ for right
          IF(IQ.LE.13) THEN
            CON=AQZ(IX,IZ)*CONJG(AQZ(IX,IZ))
          ELSE
            CON=BQZ(IX,IZ)*CONJG(BQZ(IX,IZ))
          ENDIF
          SIG=GS**2/6*FAC*CON*PSIFCN(AMQIQ**2,AMZIZ**2,T)
     $    *QFCN(JQ,1)*QFCN(1,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ,1,JTYPZ,IQ)
          SIG=GS**2/6*FAC*CON*PSIFCN(AMQIQ**2,AMZIZ**2,U)
     $    *QFCN(1,1)*QFCN(JQ,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,1,JQ,JTYPZ,IQ)
C          Jet 1 = qkss, jet 2 = ziss
320       IF(.NOT.(GOQ(IQ,1).AND.GOQ(JTYPZ,2))) GO TO 310
          CALL TWOKIN(0.,0.,AMQIQ,AMZIZ)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 310
          GS=SQRT(4*PI*ALFQSQ)
          E1=SQRT(P(1)**2+AMQIQ**2)
          E2=SQRT(P(2)**2+AMZIZ**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          IX=IS2UD(IQ)
C          Use AQZ for left squarks, BQZ for right
          IF(IQ.LE.13) THEN
            CON=AQZ(IX,IZ)*CONJG(AQZ(IX,IZ))
          ELSE
            CON=BQZ(IX,IZ)*CONJG(BQZ(IX,IZ))
          ENDIF
          SIG=GS**2/6*CON*FAC*PSIFCN(AMQIQ**2,AMZIZ**2,U)
     $    *QFCN(JQ,1)*QFCN(1,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ,1,IQ,JTYPZ)
          SIG=GS**2/6*CON*FAC*PSIFCN(AMQIQ**2,AMZIZ**2,T)
     $    *QFCN(1,1)*QFCN(JQ,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,1,JQ,IQ,JTYPZ)
310     CONTINUE
300   CONTINUE
C
C          qk gl -> wiss qkss 
C
      DO 400 IW=1,4
        JW=(IW+1)/2
        AMWIW=ABS(AMWISS(JW))
        JTYPW=IW2JS(IW)
        IWM=IW2IM(IW)
C          Left squarks only - 
        DO 410 IQ=2,11
          AMQIQ=AMASS(IDQSS(IQ))
C          JQ is the matching incoming quark
          JQ=JS2JT(IQ)
          JQ=MATCH(JQ,4)
          JQ=MATCH(JQ,IWM)
          IF(JQ.EQ.0.OR.JQ.GE.12) GO TO 410
C          Jet 1 = wiss, jet 2 = qkss
          IF(.NOT.(GOQ(JTYPW,1).AND.GOQ(IQ,2))) GO TO 420
          CALL TWOKIN(0.,0.,AMWIW,AMQIQ)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 420
          GS=SQRT(4*PI*ALFQSQ)
          E1=SQRT(P(1)**2+AMWIW**2)
          E2=SQRT(P(2)**2+AMQIQ**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          IX=IS2UD(JQ)
          CON=AQW(IX,JW)*CONJG(AQW(IX,JW))
          SIG=GS**2/6*FAC*CON*PSIFCN(AMQIQ**2,AMWIW**2,T)
     $    *QFCN(JQ,1)*QFCN(1,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ,1,JTYPW,IQ)
          SIG=GS**2/6*FAC*CON*PSIFCN(AMQIQ**2,AMWIW**2,U)
     $    *QFCN(1,1)*QFCN(JQ,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,1,JQ,JTYPW,IQ)
C          Jet 1 = qkss, jet 2 = wiss
420       IF(.NOT.(GOQ(IQ,1).AND.GOQ(JTYPW,2))) GO TO 410
          CALL TWOKIN(0.,0.,AMQIQ,AMWIW)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 410
          GS=SQRT(4*PI*ALFQSQ)
          E1=SQRT(P(1)**2+AMQIQ**2)
          E2=SQRT(P(2)**2+AMWIW**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          IX=IS2UD(JQ)
          CON=AQW(IX,JW)*CONJG(AQW(IX,JW))
          SIG=GS**2/6*FAC*CON*PSIFCN(AMQIQ**2,AMWIW**2,U)
     $    *QFCN(JQ,1)*QFCN(1,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,JQ,1,IQ,JTYPW)
          SIG=GS**2/6*FAC*CON*PSIFCN(AMQIQ**2,AMWIW**2,T)
     $    *QFCN(1,1)*QFCN(JQ,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,1,JQ,IQ,JTYPW)
410     CONTINUE
400   CONTINUE
C
C          qk qb -> wiss glss 
C
      DO 500 IW=1,4
        JW=(IW+1)/2
        AMWIW=ABS(AMWISS(JW))
        JTYPW=IW2JS(IW)
        IWM=IW2IM(IW)
C          Jet 1 = wiss, jet 2 = glss
        IF(.NOT.(GOQ(JTYPW,1).AND.GOQ(1,2))) GO TO 520
        CALL TWOKIN(0.,0.,AMWIW,AMG)
        IF(X1.GE.1..OR.X2.GE.1.) GO TO 520
        GS=SQRT(4*PI*ALFQSQ)
        E1=SQRT(P(1)**2+AMWIW**2)
        E2=SQRT(P(2)**2+AMG**2)
        FAC=1./(16.*PI*S**2)
        FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
C          Loop over quarks (no top quarks)
        DO 510 IQ=2,11
          IQ1=IQ
          IQ2=MATCH(IQ1,IWM)
          IF(IQ2.EQ.0.OR.IQ2.GE.12) GO TO 510
          AMQIQ1=AMASS(IDQSS(IQ1))
          IX1=IS2UD(IQ1)
          AMQIQ2=AMASS(IDQSS(IQ2))
          IX2=IS2UD(IQ2)
          CON11=AQW(IX1,JW)*CONJG(AQW(IX1,JW))
          CON22=AQW(IX2,JW)*CONJG(AQW(IX2,JW))
          CON12=2*(-1)**ITHG*REAL(AQW(IX1,JW)*AQW(IX2,JW))
          SIG=CON11*(AMWIW**2-T)*(AMG**2-T)/(AMQIQ2**2-T)**2
     $    +CON22*(AMWIW**2-U)*(AMG**2-U)/(AMQIQ1**2-U)**2
     $    +CON12*AMG*AMWIW*S/((AMQIQ2**2-T)*(AMQIQ1**2-U))
          SIG=2*GS**2/9*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,IQ1,IQ2,JTYPW,1)
C          No interchange needed here
510     CONTINUE
C          Jet 1 = glss, jet 2 = wiss
520     IF(.NOT.(GOQ(1,1).AND.GOQ(JTYPW,2))) GO TO 500
        CALL TWOKIN(0.,0.,AMG,AMWIW)
        IF(X1.GE.1..OR.X2.GE.1.) GO TO 500
        GS=SQRT(4*PI*ALFQSQ)
        E1=SQRT(P(1)**2+AMG**2)
        E2=SQRT(P(2)**2+AMWIW**2)
        FAC=1./(16.*PI*S**2)
        FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
C          Loop over quarks (no top quarks)
        DO 530 IQ=2,11
          IQ1=IQ
          IQ2=MATCH(IQ1,IWM)
          IF(IQ2.EQ.0.OR.IQ2.GE.12) GO TO 530
          AMQIQ1=AMASS(IDQSS(IQ1))
          IX1=IS2UD(IQ1)
          AMQIQ2=AMASS(IDQSS(IQ2))
          IX2=IS2UD(IQ2)
          CON11=AQW(IX1,JW)*CONJG(AQW(IX1,JW))
          CON22=AQW(IX2,JW)*CONJG(AQW(IX2,JW))
          CON12=2*(-1)**ITHG*REAL(AQW(IX1,JW)*AQW(IX2,JW))
          SIG=CON11*(AMWIW**2-U)*(AMG**2-U)/(AMQIQ2**2-U)**2
     $    +CON22*(AMWIW**2-T)*(AMG**2-T)/(AMQIQ1**2-T)**2
     $    +CON12*AMG*AMWIW*S/((AMQIQ2**2-U)*(AMQIQ1**2-T))
          SIG=2*GS**2/9*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
          SIG=.5*SIG
          CALL SIGFIL(SIG,IQ1,IQ2,1,JTYPW)
C       NO INTERCHANGE NEEDED HERE
530     CONTINUE
500   CONTINUE
C
C          Gaugino pair production. The W,Z poles are assumed
C          to be outside the physical region.
C          Constants from SSWZBF:
C
      SR2=SQRT(2.)
      DO 601 IZ=1,4
        XZIWJ(IZ,1)=.5*(SIGN(1.,AMWISS(1))*SIGN(1.,AMZISS(IZ))
     $  *(COS(GAMMAR)*ZMIXSS(1,IZ)/SR2+SIN(GAMMAR)*ZMIXSS(3,IZ))
     $  -COS(GAMMAL)*ZMIXSS(2,IZ)/SR2+SIN(GAMMAL)*ZMIXSS(3,IZ))
        YZIWJ(IZ,1)=.5*(-SIGN(1.,AMWISS(1))*SIGN(1.,AMZISS(IZ))
     $  *(COS(GAMMAR)*ZMIXSS(1,IZ)/SR2+SIN(GAMMAR)*ZMIXSS(3,IZ))
     $  -COS(GAMMAL)*ZMIXSS(2,IZ)/SR2+SIN(GAMMAL)*ZMIXSS(3,IZ))
        XZIWJ(IZ,2)=.5*(SIGN(1.,AMWISS(2))*SIGN(1.,AMZISS(IZ))*THY
     $  *(-SIN(GAMMAR)*ZMIXSS(1,IZ)/SR2+COS(GAMMAR)*ZMIXSS(3,IZ))
     $  +THX*(SIN(GAMMAL)*ZMIXSS(2,IZ)/SR2+COS(GAMMAL)*ZMIXSS(3,IZ)))
        YZIWJ(IZ,2)=.5*(-SIGN(1.,AMWISS(2))*SIGN(1.,AMZISS(IZ))
     $  *THY*(-SIN(GAMMAR)*ZMIXSS(1,IZ)/SR2+COS(GAMMAR)*ZMIXSS(3,IZ))
     $  +THX*(SIN(GAMMAL)*ZMIXSS(2,IZ)/SR2+COS(GAMMAL)*ZMIXSS(3,IZ)))
601   CONTINUE
C
C          Zino + wino: W* and squark graphs included
C
      DO 610 IW=1,4
        JW=(IW+1)/2
        AMWIW=ABS(AMWISS(JW))
        JTYPW=IW2JS(IW)
        IWM=IW2IM(IW)
        DO 620 IZ=1,4
          AMZIZ=ABS(AMZISS(IZ))
          JTYPZ=IZ2JS(IZ)
          AMQ=AMASS(IDQSS(2))
C          Jet 1 = wiss, jet 2 = zjss
          IF(.NOT.(GOQ(JTYPW,1).AND.GOQ(JTYPZ,2))) GO TO 630
          CALL TWOKIN(0.,0.,AMWIW,AMZIZ)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 630
          E1=SQRT(P(1)**2+AMWIW**2)
          E2=SQRT(P(2)**2+AMZIZ**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
C          Loop over quarks (no top quarks)
          SIGUT1=(XZIWJ(IZ,JW)**2+YZIWJ(IZ,JW)**2)
     $    *((AMWIW**2-U)*(AMZIZ**2-U)+(AMWIW**2-T)*(AMZIZ**2-T))/4.
     $    +2*XZIWJ(IZ,JW)*YZIWJ(IZ,JW)
     $    *((AMWIW**2-U)*(AMZIZ**2-U)-(AMWIW**2-T)*(AMZIZ**2-T))/4.
     $    +AMWIW*AMZIZ*(XZIWJ(IZ,JW)**2-YZIWJ(IZ,JW)**2)*S/2.
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGUT1=2*G**4/3./PROPW*SIGUT1
          SIGUT2=(AQZ(2,IZ)*CONJG(AQZ(2,IZ)))*
     $    (AQW(1,JW)*CONJG(AQW(1,JW)))
     $    *(AMWIW**2-U)*(AMZIZ**2-U)/4./3./(U-AMQ**2)**2
          SIGUT3=(AQZ(1,IZ)*CONJG(AQZ(1,IZ)))*
     $    (AQW(2,JW)*CONJG(AQW(2,JW)))
     $    *(AMWIW**2-T)*(AMZIZ**2-T)/4./3./(T-AMQ**2)**2
          SGUT12=-G**2*SR2*(S-AMW**2)/PROPW/(U-AMQ**2)/12.*
     $    REAL(CONJG(AQZ(2,IZ))*AQW(1,JW)*(-ZI)**(ITHZ(IZ)))*
     $    (8*(XZIWJ(IZ,JW)+YZIWJ(IZ,JW))*(AMZIZ**2-U)*(AMWIW**2-U)/4.
     $    +4*(XZIWJ(IZ,JW)-YZIWJ(IZ,JW))*AMWIW*AMZIZ*S/2.)
          SGUT13=G**2*SR2*(S-AMW**2)/PROPW/(T-AMQ**2)/12.*
     $    REAL(CONJG(AQW(2,JW))*AQZ(1,IZ)*(-ZI)**(ITHZ(IZ)))*
     $    (8*(XZIWJ(IZ,JW)-YZIWJ(IZ,JW))*(AMZIZ**2-T)*(AMWIW**2-T)/4.
     $    +4*(XZIWJ(IZ,JW)+YZIWJ(IZ,JW))*AMWIW*AMZIZ*S/2.)
          SGUT23=-4*AMWIW*AMZIZ*S/2./(U-AMQ**2)/(T-AMQ**2)/12.*
     $    REAL(AQZ(1,IZ)*AQZ(2,IZ)*CONJG(AQW(1,JW)*AQW(2,JW)))
          SIGUT=SIGUT1+SIGUT2+SIGUT3+SGUT12+SGUT13+SGUT23
C
          SIGTU1=(XZIWJ(IZ,JW)**2+YZIWJ(IZ,JW)**2)
     $    *((AMWIW**2-T)*(AMZIZ**2-T)+(AMWIW**2-U)*(AMZIZ**2-U))/4.
     $    +2*XZIWJ(IZ,JW)*YZIWJ(IZ,JW)
     $    *((AMWIW**2-T)*(AMZIZ**2-T)-(AMWIW**2-U)*(AMZIZ**2-U))/4.
     $    +AMWIW*AMZIZ*(XZIWJ(IZ,JW)**2-YZIWJ(IZ,JW)**2)*S/2.
          SIGTU1=2*G**4/3./PROPW*SIGTU1
          SIGTU2=(AQZ(2,IZ)*CONJG(AQZ(2,IZ)))*
     $    (AQW(1,JW)*CONJG(AQW(1,JW)))
     $    *(AMWIW**2-T)*(AMZIZ**2-T)/4./3./(T-AMQ**2)**2
          SIGTU3=(AQZ(1,IZ)*CONJG(AQZ(1,IZ)))*
     $    (AQW(2,JW)*CONJG(AQW(2,JW)))
     $    *(AMWIW**2-U)*(AMZIZ**2-U)/4./3./(U-AMQ**2)**2
          SGTU12=-G**2*SR2*(S-AMW**2)/PROPW/(T-AMQ**2)/12.*
     $    REAL(CONJG(AQZ(2,IZ))*AQW(1,JW)*(-ZI)**(ITHZ(IZ)))*
     $    (8*(XZIWJ(IZ,JW)+YZIWJ(IZ,JW))*(AMZIZ**2-T)*(AMWIW**2-T)/4.
     $    +4*(XZIWJ(IZ,JW)-YZIWJ(IZ,JW))*AMWIW*AMZIZ*S/2.)
          SGTU13=G**2*SR2*(S-AMW**2)/PROPW/(U-AMQ**2)/12.*
     $    REAL(CONJG(AQW(2,JW))*AQZ(1,IZ)*(-ZI)**(ITHZ(IZ)))*
     $    (8*(XZIWJ(IZ,JW)-YZIWJ(IZ,JW))*(AMZIZ**2-U)*(AMWIW**2-U)/4.
     $    +4*(XZIWJ(IZ,JW)+YZIWJ(IZ,JW))*AMWIW*AMZIZ*S/2.)
          SGTU23=-4*AMWIW*AMZIZ*S/2./(T-AMQ**2)/(U-AMQ**2)/12.*
     $    REAL(AQZ(1,IZ)*AQZ(2,IZ)*CONJG(AQW(1,JW)*AQW(2,JW)))
          SIGTU=SIGTU1+SIGTU2+SIGTU3+SGTU12+SGTU13+SGTU23
          IF (IWM.EQ.2) THEN
            SIG=.5*SIGUT*FAC*QFCN(5,1)*QFCN(2,2)
            CALL SIGFIL(SIG,5,2,JTYPW,JTYPZ)
            SIG=.5*SIGUT*FAC*QFCN(7,1)*QFCN(8,2)
            CALL SIGFIL(SIG,7,8,JTYPW,JTYPZ)
            SIG=.5*SIGTU*FAC*QFCN(2,1)*QFCN(5,2)
            CALL SIGFIL(SIG,2,5,JTYPW,JTYPZ)
            SIG=.5*SIGTU*FAC*QFCN(8,1)*QFCN(7,2)
            CALL SIGFIL(SIG,8,7,JTYPW,JTYPZ)
          ELSE
            SIG=.5*SIGTU*FAC*QFCN(4,1)*QFCN(3,2)
            CALL SIGFIL(SIG,4,3,JTYPW,JTYPZ)
            SIG=.5*SIGTU*FAC*QFCN(6,1)*QFCN(9,2)
            CALL SIGFIL(SIG,6,9,JTYPW,JTYPZ)
            SIG=.5*SIGUT*FAC*QFCN(3,1)*QFCN(4,2)
            CALL SIGFIL(SIG,3,4,JTYPW,JTYPZ)
            SIG=.5*SIGUT*FAC*QFCN(9,1)*QFCN(6,2)
            CALL SIGFIL(SIG,9,6,JTYPW,JTYPZ)
          END IF
C          Jet 1 = zjss, jet 2 = wiss
630       IF(.NOT.(GOQ(JTYPZ,1).AND.GOQ(JTYPW,2))) GO TO 620
          CALL TWOKIN(0.,0.,AMZIZ,AMWIW)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 610
          E1=SQRT(P(1)**2+AMZIZ**2)
          E2=SQRT(P(2)**2+AMWIW**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
C          Loop over quarks (no top quarks)
          SIGUT1=(XZIWJ(IZ,JW)**2+YZIWJ(IZ,JW)**2)
     $    *((AMWIW**2-U)*(AMZIZ**2-U)+(AMWIW**2-T)*(AMZIZ**2-T))/4.
     $    +2*XZIWJ(IZ,JW)*YZIWJ(IZ,JW)
     $    *((AMWIW**2-U)*(AMZIZ**2-U)-(AMWIW**2-T)*(AMZIZ**2-T))/4.
     $    +AMWIW*AMZIZ*(XZIWJ(IZ,JW)**2-YZIWJ(IZ,JW)**2)*S/2.
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGUT1=2*G**4/3./PROPW*SIGUT1
          SIGUT2=(AQZ(2,IZ)*CONJG(AQZ(2,IZ)))*
     $    (AQW(1,JW)*CONJG(AQW(1,JW)))
     $    *(AMWIW**2-U)*(AMZIZ**2-U)/4./3./(U-AMQ**2)**2
          SIGUT3=(AQZ(1,IZ)*CONJG(AQZ(1,IZ)))*
     $    (AQW(2,JW)*CONJG(AQW(2,JW)))
     $    *(AMWIW**2-T)*(AMZIZ**2-T)/4./3./(T-AMQ**2)**2
          SGUT12=-G**2*SR2*(S-AMW**2)/PROPW/(U-AMQ**2)/12.*
     $    REAL(CONJG(AQZ(2,IZ))*AQW(1,JW)*(-ZI)**(ITHZ(IZ)))*
     $    (8*(XZIWJ(IZ,JW)+YZIWJ(IZ,JW))*(AMZIZ**2-U)*(AMWIW**2-U)/4.
     $    +4*(XZIWJ(IZ,JW)-YZIWJ(IZ,JW))*AMWIW*AMZIZ*S/2.)
          SGUT13=G**2*SR2*(S-AMW**2)/PROPW/(T-AMQ**2)/12.*
     $    REAL(CONJG(AQW(2,JW))*AQZ(1,IZ)*(-ZI)**(ITHZ(IZ)))*
     $    (8*(XZIWJ(IZ,JW)-YZIWJ(IZ,JW))*(AMZIZ**2-T)*(AMWIW**2-T)/4.
     $    +4*(XZIWJ(IZ,JW)+YZIWJ(IZ,JW))*AMWIW*AMZIZ*S/2.)
          SGUT23=-4*AMWIW*AMZIZ*S/2./(U-AMQ**2)/(T-AMQ**2)/12.*
     $    REAL(AQZ(1,IZ)*AQZ(2,IZ)*CONJG(AQW(1,JW)*AQW(2,JW)))
          SIGUT=SIGUT1+SIGUT2+SIGUT3+SGUT12+SGUT13+SGUT23
C
          SIGTU1=(XZIWJ(IZ,JW)**2+YZIWJ(IZ,JW)**2)
     $    *((AMWIW**2-T)*(AMZIZ**2-T)+(AMWIW**2-U)*(AMZIZ**2-U))/4.
     $    +2*XZIWJ(IZ,JW)*YZIWJ(IZ,JW)
     $    *((AMWIW**2-T)*(AMZIZ**2-T)-(AMWIW**2-U)*(AMZIZ**2-U))/4.
     $    +AMWIW*AMZIZ*(XZIWJ(IZ,JW)**2-YZIWJ(IZ,JW)**2)*S/2.
          SIGTU1=2*G**4/3./PROPW*SIGTU1
          SIGTU2=(AQZ(2,IZ)*CONJG(AQZ(2,IZ)))*
     $    (AQW(1,JW)*CONJG(AQW(1,JW)))
     $    *(AMWIW**2-T)*(AMZIZ**2-T)/4./3./(T-AMQ**2)**2
          SIGTU3=(AQZ(1,IZ)*CONJG(AQZ(1,IZ)))*
     $    (AQW(2,JW)*CONJG(AQW(2,JW)))
     $    *(AMWIW**2-U)*(AMZIZ**2-U)/4./3./(U-AMQ**2)**2
          SGTU12=-G**2*SR2*(S-AMW**2)/PROPW/(T-AMQ**2)/12.*
     $    REAL(CONJG(AQZ(2,IZ))*AQW(1,JW)*(-ZI)**(ITHZ(IZ)))*
     $    (8*(XZIWJ(IZ,JW)+YZIWJ(IZ,JW))*(AMZIZ**2-T)*(AMWIW**2-T)/4.
     $    +4*(XZIWJ(IZ,JW)-YZIWJ(IZ,JW))*AMWIW*AMZIZ*S/2.)
          SGTU13=G**2*SR2*(S-AMW**2)/PROPW/(U-AMQ**2)/12.*
     $    REAL(CONJG(AQW(2,JW))*AQZ(1,IZ)*(-ZI)**(ITHZ(IZ)))*
     $    (8*(XZIWJ(IZ,JW)-YZIWJ(IZ,JW))*(AMZIZ**2-U)*(AMWIW**2-U)/4.
     $    +4*(XZIWJ(IZ,JW)+YZIWJ(IZ,JW))*AMWIW*AMZIZ*S/2.)
          SGTU23=-4*AMWIW*AMZIZ*S/2./(T-AMQ**2)/(U-AMQ**2)/12.*
     $    REAL(AQZ(1,IZ)*AQZ(2,IZ)*CONJG(AQW(1,JW)*AQW(2,JW)))
          SIGTU=SIGTU1+SIGTU2+SIGTU3+SGTU12+SGTU13+SGTU23
          IF (IWM.EQ.2) THEN
            SIG=.5*SIGTU*FAC*QFCN(5,1)*QFCN(2,2)
            CALL SIGFIL(SIG,5,2,JTYPZ,JTYPW)
            SIG=.5*SIGTU*FAC*QFCN(7,1)*QFCN(8,2)
            CALL SIGFIL(SIG,7,8,JTYPZ,JTYPW)
            SIG=.5*SIGUT*FAC*QFCN(2,1)*QFCN(5,2)
            CALL SIGFIL(SIG,2,5,JTYPZ,JTYPW)
            SIG=.5*SIGUT*FAC*QFCN(8,1)*QFCN(7,2)
            CALL SIGFIL(SIG,8,7,JTYPZ,JTYPW)
          ELSE
            SIG=.5*SIGUT*FAC*QFCN(4,1)*QFCN(3,2)
            CALL SIGFIL(SIG,4,3,JTYPZ,JTYPW)
            SIG=.5*SIGUT*FAC*QFCN(6,1)*QFCN(9,2)
            CALL SIGFIL(SIG,6,9,JTYPZ,JTYPW)
            SIG=.5*SIGTU*FAC*QFCN(3,1)*QFCN(4,2)
            CALL SIGFIL(SIG,3,4,JTYPZ,JTYPW)
            SIG=.5*SIGTU*FAC*QFCN(9,1)*QFCN(6,2)
            CALL SIGFIL(SIG,9,6,JTYPZ,JTYPW)
          END IF
620     CONTINUE
610   CONTINUE
C
C          Chargino pair production
C          added squark exchange contribution 7/11/97
C
      DO 700 IW1=1,4
        JW1=(IW1+1)/2
        AMWIW1=ABS(AMWISS(JW1))
        JTYPW1=IW2JS(IW1)
        IDW1=IDWSS(IW1)
        DO 710 IW2=1,4
          JW2=(IW2+1)/2
          AMWIW2=ABS(AMWISS(JW2))
          JTYPW2=IW2JS(IW2)
          IDW2=IDWSS(IW2)
          IF (.NOT.(GOQ(JTYPW1,1).AND.GOQ(JTYPW2,2))) GO TO 710
          CALL TWOKIN(0.,0.,AMWIW1,AMWIW2)
          IF (X1.GE.1..OR.X2.GE.1.) GO TO 710
          E1=SQRT(P(1)**2+AMWIW1**2)
          E2=SQRT(P(2)**2+AMWIW2**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          DO 720 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IF (IFLQ.EQ.1) THEN
              EQ1=2./3.
            ELSE
              EQ1=-1./3.
            END IF
            IQ2=MATCH(IQ1,4)
            IF (IQ1.EQ.2.OR.IQ1.EQ.3) AMSQK=AMDLSS
            IF (IQ1.EQ.4.OR.IQ1.EQ.5) AMSQK=AMULSS
            IF (IQ1.EQ.6.OR.IQ1.EQ.7) AMSQK=AMCLSS
            IF (IQ1.EQ.8.OR.IQ1.EQ.9) AMSQK=AMSLSS
            IF (IQ1.EQ.10.OR.IQ1.EQ.11) AMSQK=AMB1SS
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 720
            IF (IDW1.EQ.-IDW2) THEN
C          Convert ISAJET t_hat to particle-particle t_hat
              IF (IUD(IQ1)*IDW1.GT.0) THEN
                TPP=U
              ELSE
                TPP=T
              END IF
              ZZ=(2*TPP-2*AMWIW1**2+S)/SQRT(S*S-4*S*AMWIW1**2)
              EHAT=SQRT(S)/2.
              PHAT=SQRT(EHAT**2-AMWIW1**2)
              XMGG=16.*ESQ*ESQ*(EHAT**2*(1.+ZZ**2)+
     $        AMWIW1**2*(1.-ZZ**2))/S*EQ1**2
              XMZZ=16*ESQ*ESQ*CTTHW**2*S/((S-AMZ**2)**2+
     $        (GAMZ*AMZ)**2)*((XWI(JW1)**2+YWI(JW1)**2)*
     $        (AL(IFLQ)**2+BE(IFLQ)**2)*
     $        (EHAT**2*(1.+ZZ**2)+AMWIW1**2*(1.-ZZ**2))-2.*
     $        YWI(JW1)**2*(AL(IFLQ)**2+
     $        BE(IFLQ)**2)*AMWIW1**2-8*XWI(JW1)*YWI(JW1)*
     $        AL(IFLQ)*BE(IFLQ)*EHAT*PHAT*ZZ)
              XMGZ=(-EQ1)*(-32.)*ESQ*ESQ*CTTHW*(S-AMZ**2)/
     $        ((S-AMZ**2)**2+(GAMZ*AMZ)**2)*
     $        (AL(IFLQ)*XWI(JW1)*(EHAT**2*
     $        (1.+ZZ**2)+AMWIW1**2*(1.-ZZ**2))-2*
     $        BE(IFLQ)*YWI(JW1)*EHAT*PHAT*ZZ)
              XMUU=ESQ*ESQ*SIN(GAMMAR)**4*S*(EHAT-PHAT*ZZ)**2/
     $         SN2THW**2/(EHAT**2+PHAT**2-2*EHAT*PHAT*ZZ+
     $         AMSQK**2)**2
              XMGU=EQ1*4*ESQ*ESQ*SIN(GAMMAR)**2*
     $         ((EHAT-PHAT*ZZ)**2+AMWIW1**2)/SN2THW/
     $         (EHAT**2+PHAT**2-2*EHAT*PHAT*ZZ+AMSQK**2)
              XMZU=4*ESQ*ESQ*CTTHW*SIN(GAMMAR)**2*(S-AMZ**2)
     $         *(AL(IFLQ)-BE(IFLQ))*S/SN2THW/((S-AMZ**2)**2+
     $         (GAMZ*AMZ)**2)*((XWI(JW1)-YWI(JW1))*
     $         ((EHAT-PHAT*ZZ)**2+AMWIW1**2)+2*YWI(JW1)*
     $         AMWIW1**2)/(EHAT**2+PHAT**2-2*EHAT*PHAT*ZZ+
     $         AMSQK**2)
              XMDD=ESQ*ESQ*SIN(GAMMAL)**4*S*(EHAT+PHAT*ZZ)**2/
     $         SN2THW**2/(EHAT**2+PHAT**2+2*EHAT*PHAT*ZZ+
     $         AMSQK**2)**2
              XMGD=-4*EQ1*ESQ*ESQ*SIN(GAMMAL)**2*
     $         ((EHAT+PHAT*ZZ)**2+AMWIW1**2)/SN2THW/
     $         (EHAT**2+PHAT**2+2*EHAT*PHAT*ZZ+AMSQK**2)
              XMZD=-4*ESQ*ESQ*CTTHW*SIN(GAMMAL)**2*(S-AMZ**2)
     $         *(AL(IFLQ)-BE(IFLQ))*S/SN2THW/((S-AMZ**2)**2+
     $         (GAMZ*AMZ)**2)*((XWI(JW1)+YWI(JW1))*
     $         ((EHAT+PHAT*ZZ)**2+AMWIW1**2)-2*YWI(JW1)*
     $         AMWIW1**2)/(EHAT**2+PHAT**2+2*EHAT*PHAT*ZZ+
     $         AMSQK**2)
              IF (IFLQ.EQ.1) THEN
               SIG=(XMGG+XMZZ+XMGZ+XMDD+XMGD+XMZD)/12.
              ELSE
               SIG=(XMGG+XMZZ+XMGZ+XMUU+XMGU+XMZU)/12.
              END IF
              SIG=SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
              SIG=.5*SIG
C              IF(SIG.LT.0.AND.ABS(ZZ).GT.0.999) SIG=0
              CALL SIGFIL(SIG,IQ1,IQ2,JTYPW1,JTYPW2)
            ELSEIF (IDW1*IDW2.LT.0) THEN
              PHAT=SQRT(S*S+AMWIW1**4+AMWIW2**4-2*S*AMWIW1**2
     $        -2*S*AMWIW2**2-2*AMWIW1**2*AMWIW2**2)/2./SQRT(S)
              IF (IUD(IQ1)*IDW1.GT.0) THEN
                TPP=U
              ELSE
                TPP=T
              END IF
              IF (IDW1.LT.0) THEN
                AMWI=AMWIW1
              ELSE
                AMWI=AMWIW2
              END IF
              EHAT=SQRT(PHAT**2+AMWI**2)
              EBM=SQRT(S)/2.
              ZZ=(TPP-AMWI**2+SQRT(S)*EHAT)/SQRT(S)/PHAT
              DEL=(AMW2SS**2-AMW1SS**2)/4./EBM
              XMZZ=4*(CTTHW+TNTHW)**2/((S-AMZ**2)**2+
     $        (GAMZ*AMZ)**2)*((X12**2+Y12**2)*
     $        (AL(IFLQ)**2+BE(IFLQ)**2)*
     $        (EBM**2+PHAT**2*ZZ**2-DEL**2-SN12*AMWIW1*AMWIW2)+
     $        2*X12**2*SN12*(AL(IFLQ)**2+ BE(IFLQ)**2)*AMWIW1*
     $        AMWIW2-8*X12*Y12*AL(IFLQ)*BE(IFLQ)*EBM*PHAT*ZZ)
              XMUU=SIN(GAMMAR)**2*COS(GAMMAR)**2*((EBM-PHAT*ZZ)
     $         **2-DEL**2)/SN2THW**2/(2*EBM*(EBM-DEL)-2*EBM*PHAT*
     $         ZZ+AMSQK**2-AMW1SS**2)**2
              XMZU=-2*THY*(CTTHW+TNTHW)*SIN(GAMMAR)*COS(GAMMAR)*
     $         (S-AMZ**2)*(AL(IFLQ)-BE(IFLQ))/SN2THW/((S-AMZ**2)
     $         **2+(GAMZ*AMZ)**2)*((X12-Y12)*((EBM-PHAT*ZZ)**2-
     $         DEL**2-SN12*AMWIW1*AMWIW2)+2*X12*SN12*AMWIW1*
     $         AMWIW2)/(2*EBM*(EBM-DEL)-2*EBM*PHAT*ZZ+AMSQK**2
     $         -AMW1SS**2)
              XMDD=SIN(GAMMAL)**2*COS(GAMMAL)**2*((EBM+PHAT*ZZ)
     $         **2-DEL**2)/SN2THW**2/(2*EBM*(EBM-DEL)+2*EBM*PHAT*
     $         ZZ+AMSQK**2-AMW1SS**2)**2
              XMZD=-2*THX*(CTTHW+TNTHW)*SIN(GAMMAL)*COS(GAMMAL)*
     $         (S-AMZ**2)*(AL(IFLQ)-BE(IFLQ))/SN2THW/((S-AMZ**2)
     $         **2+(GAMZ*AMZ)**2)*((X12+Y12)*((EBM+PHAT*ZZ)**2-
     $         DEL**2+SN12*AMWIW1*AMWIW2)-2*Y12*SN12*AMWIW1*
     $         AMWIW2)/(2*EBM*(EBM-DEL)+2*EBM*PHAT*ZZ+AMSQK**2
     $         -AMW1SS**2)
              IF (IFLQ.EQ.1) THEN
               SIG=ESQ*ESQ*(XMZZ+XMDD+XMZD)*S/12.
              ELSE
               SIG=ESQ*ESQ*(XMZZ+XMUU+XMZU)*S/12.
              END IF
              SIG=SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
              SIG=.5*SIG
              CALL SIGFIL(SIG,IQ1,IQ2,JTYPW1,JTYPW2)
            END IF
720       CONTINUE
710     CONTINUE
700   CONTINUE
C
C         qk qb --> ziss zjss
C
      DO 800 IZ1=1,4
        AMZIZ1=ABS(AMZISS(IZ1))
        JTYPZ1=IZ2JS(IZ1)
        DO 810 IZ2=1,4
          AMZIZ2=ABS(AMZISS(IZ2))
          JTYPZ2=IZ2JS(IZ2)
          IF(.NOT.(GOQ(JTYPZ1,1).AND.GOQ(JTYPZ2,2))) GO TO 810
          CALL TWOKIN(0.,0.,AMZIZ1,AMZIZ2)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 810
          E1=SQRT(P(1)**2+AMZIZ1**2)
          E2=SQRT(P(2)**2+AMZIZ2**2)
          FAC=1./(16.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          WIJ=SQRT(G**2+GP**2)*ZI**(ITHZ(IZ2))*(-ZI)**(ITHZ(IZ1))*
     $    (ZMIXSS(1,IZ1)*ZMIXSS(1,IZ2)-ZMIXSS(2,IZ1)*
     $    ZMIXSS(2,IZ2))/4.
          RSH=SQRT(S)
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          KK=SQRT(S*S+(AMZIZ1**2-AMZIZ2**2)**2-2*S*
     $    (AMZIZ1**2+AMZIZ2**2))/2./RSH
C          Sum over initial quarks (no top quarks)
          DO 820 IQ=2,11
            IQ1=IQ
            IQ2=MATCH(IQ1,4)
            AMSQL=AMASS(IDQSS(IQ))
            AMSQR=AMASS(IDQSS(IQ+12))
            PHAT=SQRT(SSXLAM(S,AMZIZ1**2,AMZIZ2**2))/2./RSH
            EHAT=SQRT(PHAT**2+AMZIZ1**2)
            ZZ=(T-AMZIZ1**2+RSH*EHAT)/RSH/PHAT
            IF (IUD(IQ).LT.0) ZZ=-ZZ
            IFLQ=IS2UD(IQ)
            SIGLL=AQZ(IFLQ,IZ1)*CONJG(AQZ(IFLQ,IZ1))*AQZ(IFLQ,IZ2)*
     $      CONJG(AQZ(IFLQ,IZ2))*SSGT(S,AMSQL,ZZ,IZ1,IZ2)
            SIGRR=BQZ(IFLQ,IZ1)*CONJG(BQZ(IFLQ,IZ1))*BQZ(IFLQ,IZ2)*
     $      CONJG(BQZ(IFLQ,IZ2))*SSGT(S,AMSQR,ZZ,IZ1,IZ2)
            SIGZZ=4*ESQ*WIJ*CONJG(WIJ)*(AL(IFLQ)**2+BE(IFLQ)**2)*
     $      (S*S-(AMZIZ1**2-AMZIZ2**2)**2+4*(-1.)**(ITHZ(IZ1)+
     $      ITHZ(IZ2)+1)*S*AMZIZ1*AMZIZ2+4*S*KK*KK*ZZ*ZZ)/PROPZ
            SIGLZ=-SQRT(ESQ)*(AL(IFLQ)-BE(IFLQ))*(S-AMZ**2)/2./
     $      PROPZ*(REAL(WIJ*CONJG(AQZ(IFLQ,IZ1))*AQZ(IFLQ,IZ2))*
     $      SSGST(S,AMSQL,ZZ,IZ1,IZ2)+(-1.)**(ITHZ(IZ1)+ITHZ(IZ2))*
     $      REAL(WIJ*AQZ(IFLQ,IZ1)*CONJG(AQZ(IFLQ,IZ2)))*
     $      SSGST(S,AMSQL,-ZZ,IZ1,IZ2))
            SIGRZ=-SQRT(ESQ)*(-1.)**(ITHZ(IZ1)+ITHZ(IZ2)+1)*
     $      (AL(IFLQ)+BE(IFLQ))*(S-AMZ**2)/2./
     $      PROPZ*(REAL(WIJ*CONJG(BQZ(IFLQ,IZ1))*BQZ(IFLQ,IZ2))*
     $      SSGST(S,AMSQR,ZZ,IZ1,IZ2)+(-1.)**(ITHZ(IZ1)+ITHZ(IZ2))*
     $      REAL(WIJ*BQZ(IFLQ,IZ1)*CONJG(BQZ(IFLQ,IZ2)))*
     $      SSGST(S,AMSQR,-ZZ,IZ1,IZ2))
            SIG=KK*(SIGLL+SIGRR+SIGZZ+SIGLZ+SIGRZ)/3./PHAT
C          Below factor of 2 for id particles and jettyp switch
            SIG=SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)/2.
            IF(SIG.LT.0.AND.ABS(ZZ).GT.0.999) SIG=0
            CALL SIGFIL(SIG,IQ1,IQ2,JTYPZ1,JTYPZ2)
820       CONTINUE
810     CONTINUE
800   CONTINUE
      RETURN
      END
