#include "isajet/pilot.h"
      SUBROUTINE SIGWW
C
C          Calculate D(SIGMA)/D(PT**2)D(Y1)D(Y2) for QK+QB-->W+W
C          summed over W types allowed on JETTYPE cards and
C          including branching ratio implied by WMODE cards.
C
C          SIGMA    = cross section summed over quark types allowed by
C                     JETTYPE card.
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4.
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 + I1
C                     using JETTYPE code.
C
C          Cross sections from Brown and Mikaelian,
C          Phys Rev D19, 922, D20, 1164.
C          Include extra factor of 1/2 for double counting.
C
C          Double precision needed for 32-bit machines.
C
C          Ver. 6.22: Modified to used W + GM decay distributions from
C                     Cortes, Hagiwara, and Herzog, NP B278, 26 (1986)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/qcdpar.inc"
#include "isajet/jetpar.inc"
#include "isajet/primar.inc"
#include "isajet/q1q2.inc"
#include "isajet/jetsig.inc"
#include "isajet/const.inc"
#include "isajet/qsave.inc"
#include "isajet/wcon.inc"
#include "isajet/wwpar.inc"
C
      DIMENSION X(2),LISTW(4),QSGN(6)
      EQUIVALENCE (X(1),X1)
      EQUIVALENCE (S,SWW),(T,TWW),(U,UWW)
#if defined(CERNLIB_SINGLE)
      REAL S,T,U,TX,UX,TT,UU
     $,WWA,WWI,WWE,WZA,WZI,WZE,TERM
     $,GA,GI,GE,GJ,GZ
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION S,T,U,TX,UX,TT,UU
     $,WWA,WWI,WWE,WZA,WZI,WZE,TERM
     $,GA,GI,GE,GJ,GZ
#endif
      REAL WM2S,ZM2S,X,STRUC,FJAC,SGN,QSGN,SIG,FACTOR,EQ3(12)
      INTEGER I,IH,IQ,IW1,IW2,JW,JZ,IW,IQ1,IQ2,JG,LISTW,IFOUR
      INTEGER IFLI,IFLJ
      LOGICAL LQK1
C
      DATA LISTW/10,80,-80,90/
      DATA QSGN/1.,-1.,-1.,1.,-1.,1./
      DATA EQ3/2.,-1.,-1.,2.,-1.,2.,0.,-3.,0.,-3.,0.,-3./
C
C          Functions for W+W-
      WWA(S,T,U)=(U*T/WM2**2-1.)*(.25-WM2/S+3.*(WM2/S)**2)+S/WM2-4.
      WWI(S,T,U)=(U*T/WM2**2-1.)*(.25-.5*WM2/S-WM2**2/(S*T))
     $+S/WM2-2.+2.*WM2/T
      WWE(S,T,U)=(U*T/WM2**2-1.)*(.25+(WM2/T)**2)+S/WM2
C          Functions for W+-Z0
      WZA(S,T,U)=(U*T/(WM2*ZM2)-1.)*(.25-(WM2+ZM2)/(2.*S)
     $+((WM2+ZM2)**2+8.*WM2*ZM2)/(4.*S**2))
     $+(WM2+ZM2)/(WM2*ZM2)*(.5*S-WM2-ZM2+(WM2-ZM2)**2/(2.*S))
      WZI(S,T,U)=.25*(U*T/(WM2*ZM2)-1.)*(1.-(WM2+ZM2)/S
     $-4.*WM2*ZM2/(S*T))
     $+(WM2+ZM2)/(2.*WM2*ZM2)*(S-WM2-ZM2+2.*WM2*ZM2/T)
      WZE(S,T,U)=.25*(U*T/(WM2*ZM2)-1.)+.5*S*(WM2+ZM2)/(WM2*ZM2)
C
C          Initialize
      DO 10 I=1,MXSIGS
10    SIGS(I)=0.
      SIGMA=0.
      NSIGS=0
C
C          Convention is that even for double precision single
C          precision mass is exact.
      WM2=WMASS(2)
      WM2=WM2**2
      ZM2=WMASS(4)
      ZM2=ZM2**2
C          Also need single precision mass**2.
      WM2S=WM2
      ZM2S=ZM2
C
C          W+ W- pairs
C
      IF(.NOT.((GOQ(2,1).AND.GOQ(3,2)).OR.(GOQ(3,1).AND.GOQ(2,2))))
     $GO TO 200
      CALL WWKIN(WMASS(2),WMASS(2))
      IF(X1.GE.1..OR.X2.GE.1.) GO TO 200
      DO 110 IH=1,2
      DO 110 IQ=2,9
110   QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
      FJAC=S/SCM*UNITS
      FJAC=FJAC*PI*ALFA**2/(3.*S**2)
      FJAC=FJAC*P(1)*P(2)/SQRT((P(1)**2+WM2S)*(P(2)**2+WM2S))
      FJAC=.5*FJAC
C          Sum over jet1 = W+ and jet2 = W+.
C          Swap t and u in latter case.
      DO 120 IW1=2,3
      IW2=5-IW1
      IF(.NOT.(GOQ(IW1,1).AND.GOQ(IW2,2))) GO TO 120
      IF(IW1.EQ.3) GO TO 121
      TX=T
      UX=U
      GO TO 122
121   TX=U
      UX=T
C
C          Sum over quarks, swapping t and u for negative charge.
122   DO 130 IQ=1,4
      GA=2.*(AQDP(IQ,1)+EZDP*AQDP(IQ,4)*S/(S-ZM2))**2
     $+2.*(EZDP*BQDP(IQ,4)*S/(S-ZM2))**2
      GI=8.*(AQDP(IQ,1)+EZDP*(AQDP(IQ,4)+BQDP(IQ,4))*S/(S-ZM2))
     $*(AQDP(IQ,2))**2
      GE=16.*(AQDP(IQ,2))**4
      SGN=QSGN(IQ)
      IF(SGN.LT.0.) GO TO 131
      TT=TX
      UU=UX
      GO TO 132
131   TT=UX
      UU=TX
132   SIG=QSAVE(2*IQ,1)*QSAVE(2*IQ+1,2)*FJAC*TBRWW(IW1,1)*TBRWW(IW2,2)
     $*(GA*WWA(S,TT,UU)-SGN*GI*WWI(S,TT,UU)+GE*WWE(S,TT,UU))
      CALL SIGFIL(SIG,2*IQ,2*IQ+1,IW1,IW2)
      SIG=QSAVE(2*IQ+1,1)*QSAVE(2*IQ,2)*FJAC*TBRWW(IW1,1)*TBRWW(IW2,2)
     $*(GA*WWA(S,UU,TT)-SGN*GI*WWI(S,UU,TT)+GE*WWE(S,UU,TT))
      CALL SIGFIL(SIG,2*IQ+1,2*IQ,IW1,IW2)
130   CONTINUE
120   CONTINUE
C
C          Z0 Z0 pairs
C
200   IF(.NOT.(GOQ(4,1).AND.GOQ(4,2))) GO TO 300
      CALL WWKIN(WMASS(4),WMASS(4))
      IF(X1.GE.1..OR.X2.GE.1.) RETURN
      DO 210 IH=1,2
      DO 210 IQ=2,9
210   QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
C          Jacobean -- including factor of 1/2 for identical particles.
      FJAC=.5*S/SCM*UNITS
      FJAC=FJAC*PI*ALFA**2/(3.*S**2)
      FJAC=FJAC*P(1)*P(2)/SQRT((P(1)**2+ZM2S)*(P(2)**2+ZM2S))
      DO 220 IQ=1,4
      GZ=2.*(AQDP(IQ,4)**4+BQDP(IQ,4)**4
     $+6.*AQDP(IQ,4)**2*BQDP(IQ,4)**2)
      FACTOR=(T/U+U/T+4.*ZM2*S/(T*U)-ZM2**2*(1./T**2+1./U**2))
      FACTOR=FACTOR*FJAC*GZ*TBRWW(4,1)*TBRWW(4,2)
      SIG=FACTOR*QSAVE(2*IQ,1)*QSAVE(2*IQ+1,2)
      CALL SIGFIL(SIG,2*IQ,2*IQ+1,4,4)
      SIG=FACTOR*QSAVE(2*IQ+1,1)*QSAVE(2*IQ,2)
      CALL SIGFIL(SIG,2*IQ+1,2*IQ,4,4)
220   CONTINUE
C
C          W+- Z0 pairs
C
C          JW and JZ are W+- and Z0 jet numbers.
300   DO 310 JW=1,2
        JZ=3-JW
        IF(.NOT.((GOQ(2,JW).OR.GOQ(3,JW)).AND.GOQ(4,JZ))) GO TO 310
C
C          Must swap t and u if JW=2.
        IF(JW.EQ.1) THEN
          CALL WWKIN(WMASS(2),WMASS(4))
          TX=T
          UX=U
          FJAC=S/SCM*UNITS
          FJAC=FJAC*PI*ALFA**2/(3.*S**2)
          FJAC=.5*FJAC
          FJAC=FJAC*P(1)*P(2)/SQRT((P(1)**2+WM2S)*(P(2)**2+ZM2S))
        ELSE
          CALL WWKIN(WMASS(4),WMASS(2))
          TX=U
          UX=T
          FJAC=S/SCM*UNITS
          FJAC=FJAC*PI*ALFA**2/(3.*S**2)
          FJAC=.5*FJAC
          FJAC=FJAC*P(1)*P(2)/SQRT((P(1)**2+ZM2S)*(P(2)**2+WM2S))
        ENDIF
        IF(X1.GE.1..OR.X2.GE.1.) GO TO 310
        DO 320 IH=1,2
        DO 320 IQ=1,9
320     QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
C
C          Sum over W+ and W-
        DO 340 IW=2,3
          IF(IW.EQ.2) THEN
            SGN=+1
          ELSE
            SGN=+1
          ENDIF
C
C          Sum over quarks, swapping t and u as needed.
          DO 350 IQ1=2,9
            IQ2=MATCH(IQ1,IW)
            IF(IQ2.EQ.0) GO TO 350
            IQ=IQ1/2
            IF(2*IQ.EQ.IQ1) THEN
              LQK1=.TRUE.
            ELSE
              LQK1=.FALSE.
            ENDIF
            IF((LQK1.AND.IW.EQ.3).OR.(.NOT.LQK1.AND.IW.EQ.2)) THEN
              TT=TX
              UU=UX
              IFLI=IQ1/2
              IFLJ=IQ2/2
            ELSE
              TT=UX
              UU=TX
              IFLI=IQ2/2
              IFLJ=IQ1/2
            ENDIF
C
            GA=AQDP(IQ,IW)*EZDP*S/(S-WM2)
            GI=AQDP(IQ,IW)*(AQDP(IFLI,4)+BQDP(IFLI,4))
            GJ=AQDP(IQ,IW)*(AQDP(IFLJ,4)+BQDP(IFLJ,4))
            TERM=GA**2*WZA(S,TT,UU)
            TERM=TERM+2.*GA*SGN*(-GJ*WZI(S,TT,UU)+GI*WZI(S,UU,TT))
            TERM=TERM+(GI-GJ)**2*WZE(S,TT,UU)
            TERM=TERM+GI**2*(UU*TT-WM2*ZM2)/UU**2
     $      +2.*GI*GJ*S*(WM2+ZM2)/(TT*UU)+GJ**2*(UU*TT-WM2*ZM2)/TT**2
            TERM=TERM*4.*FJAC*QSAVE(IQ1,1)*QSAVE(IQ2,2)
            TERM=TERM*TBRWW(IW,JW)*TBRWW(4,JZ)
            SIG=TERM
            IF(JW.EQ.1) THEN
              CALL SIGFIL(SIG,IQ1,IQ2,IW,4)
            ELSE
              CALL SIGFIL(SIG,IQ1,IQ2,4,IW)
            ENDIF
350       CONTINUE
340     CONTINUE
310   CONTINUE
C
C          W+- GM pairs.
C
400   DO 410 JW=1,2
        JG=3-JW
        IF(.NOT.((GOQ(2,JW).OR.GOQ(3,JW)).AND.GOQ(1,JG))) GO TO 410
C
C          Must swap t and u if JW=2.
        IF(JW.EQ.1) THEN
          CALL WWKIN(WMASS(2),0.)
          TX=T
          UX=U
          FJAC=S/SCM*UNITS
          FJAC=FJAC*PI*ALFA**2/S**2
          FJAC=.5*FJAC
          FJAC=FJAC*P(1)/SQRT(P(1)**2+WM2S)
        ELSE
          CALL WWKIN(0.,WMASS(2))
          TX=U
          UX=T
          FJAC=S/SCM*UNITS
          FJAC=FJAC*PI*ALFA**2/S**2
          FJAC=.5*FJAC
          FJAC=FJAC*P(2)/SQRT(P(2)**2+WM2S)
        ENDIF
C
        IF(X1.GE.1..OR.X2.GE.1.) GO TO 410
        DO 420 IH=1,2
        DO 420 IQ=1,9
420     QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
C
C          Sum over W+ and W-
        DO 440 IW=2,3
C
C          Sum over quarks, swapping t and u as needed.
          DO 450 IQ1=2,9
            IQ2=MATCH(IQ1,IW)
            IF(IQ2.EQ.0) GO TO 450
            IQ=IQ1/2
            IF(2*IQ.EQ.IQ1) THEN
              LQK1=.TRUE.
            ELSE
              LQK1=.FALSE.
            ENDIF
            IF((LQK1.AND.IW.EQ.3).OR.(.NOT.LQK1.AND.IW.EQ.2)) THEN
              TT=TX
              UU=UX
            ELSE
              TT=UX
              UU=TX
            ENDIF
C
            SIG=TBRWW(IW,JW)/(6.*SIN2W)*(-1./3.+UU/(TT+UU))**2
     $      *(UU**2+TT**2+2.*S*WM2)/(TT*UU)
            SIG=SIG*FJAC*QSAVE(IQ1,1)*QSAVE(IQ2,2)
            IF(JW.EQ.1) CALL SIGFIL(SIG,IQ1,IQ2,IW,1)
            IF(JW.EQ.2) CALL SIGFIL(SIG,IQ1,IQ2,1,IW)
450       CONTINUE
440     CONTINUE
410   CONTINUE
C
C          Z0 GM pairs
C
      IF (GOQ(4,1).AND.GOQ(1,2)) THEN
      CALL WWKIN(WMASS(4),0.)
      IF(X1.GE.1..OR.X2.GE.1.) GO TO 500
      DO 510 IH=1,2
      DO 510 IQ=2,9
510   QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
      FJAC=S/SCM*P(1)/SQRT(P(1)**2+ZM2S)*UNITS
      FJAC=FJAC*PI*ALFA**2/(3.*S**2)
      DO 520 IQ=1,4
      GZ=AQDP(IQ,4)**2+(AQDP(IQ,4)-BQDP(IQ,4))**2
      FACTOR=(S**2+ZM2**2)/2./T/U+1.
      FACTOR=(EQ3(IQ)/3.)**2*FACTOR*FJAC*GZ*TBRWW(4,1)
      SIG=FACTOR*QSAVE(2*IQ,1)*QSAVE(2*IQ+1,2)
      CALL SIGFIL(SIG,2*IQ,2*IQ+1,4,1)
      SIG=FACTOR*QSAVE(2*IQ+1,1)*QSAVE(2*IQ,2)
      CALL SIGFIL(SIG,2*IQ+1,2*IQ,4,1)
520   CONTINUE
500   CONTINUE
      END IF
C
      IF (GOQ(1,1).AND.GOQ(4,2)) THEN
      CALL WWKIN(0.,WMASS(4))
      IF(X1.GE.1..OR.X2.GE.1.) GO TO 600
      DO 610 IH=1,2
      DO 610 IQ=2,9
610   QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
      FJAC=S/SCM*P(2)/SQRT(P(2)**2+ZM2S)*UNITS
      FJAC=FJAC*PI*ALFA**2/(3.*S**2)
      DO 620 IQ=1,4
      GZ=AQDP(IQ,4)**2+(AQDP(IQ,4)-BQDP(IQ,4))**2
      FACTOR=(S**2+ZM2**2)/2./T/U+1.
      FACTOR=(EQ3(IQ)/3.)**2*FACTOR*FJAC*GZ*TBRWW(4,2)
      SIG=FACTOR*QSAVE(2*IQ,1)*QSAVE(2*IQ+1,2)
      CALL SIGFIL(SIG,2*IQ,2*IQ+1,1,4)
      SIG=FACTOR*QSAVE(2*IQ+1,1)*QSAVE(2*IQ,2)
      CALL SIGFIL(SIG,2*IQ+1,2*IQ,1,4)
620   CONTINUE
600   CONTINUE
      END IF
C
      RETURN
      END
