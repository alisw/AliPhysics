#include "isajet/pilot.h"
      SUBROUTINE SIGWH
C
C          Calculate d(sigma)/d(pt**2)d(y1)d(y2) for WH and ZH
C          associated production.
C
C          SIGMA    = cross section summed over types allowed by
C                     JETTYPE cards.
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 +I1
C
C          Extra factor of 1/2 needed for nonidentical final jets.
C          Y=-log(tan(theta/2)) gives jacobean P1*P2/E1*E2
C
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
#include "isajet/wcon.inc"
#include "isajet/hcon.inc"
#include "isajet/xmssm.inc"
C
      REAL X(2)
      EQUIVALENCE (X(1),X1)
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT)
      REAL SIG,S,T,U,FAC,AMW,AMZ,AMW2,AMZ2,E1,E2
      REAL QFCN,STRUC,SIGHW
      REAL PROPZ,PROPW,GV(2),GA(2),AMH,AMH2,GAMW,GAMZ
      INTEGER IS2UD(25),IQ,IH,I,IQ1,IQ2,IFLQ
      SAVE IS2UD
C
C          IS2UD: Susy jettype -> u/d code
      DATA IS2UD/0,1,1,2,2,2,2,1,1,2,2,1,1,1,1,2,2,2,2,1,1,2,2,1,1/

C          Functions
      QFCN(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
C
      IF (GOMSSM) THEN
        CALL SIGWHS
        RETURN
      END IF
C          Initialize
      DO 10 I=1,MXSIGS
10    SIGS(I)=0.
      SIGMA=0.
      NSIGS=0
C
      AMW=WMASS(2)
      AMW2=AMW**2
      AMZ=WMASS(4)
      AMZ2=AMZ**2
      AMH=HMASS
      AMH2=AMH**2
      GAMW=WGAM(2)
      GAMZ=WGAM(4)
      GV(1)=.25-2*SIN2W/3.
      GV(2)=-.25+SIN2W/3.
      GA(1)=-.25
      GA(2)=.25
C
C          WH production via W-*
C
      IF (GOQ(28,1).AND.GOQ(30,2)) THEN
          CALL TWOKIN(0.,0.,AMW,AMH)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 100
          E1=SQRT(P(1)**2+AMW**2)
          E2=SQRT(P(2)**2+AMH**2)
          FAC=1./(12.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGHW=GF**2*AMW**8*(S/AMW2+(1.-T/AMW2)*(1.-U/AMW2))/
     $          PROPW*TBRWW(3,1)
          SIG=.5*SIGHW*FAC*QFCN(3,1)*QFCN(4,2)
          CALL SIGFIL(SIG,3,4,28,30)
          SIG=.5*SIGHW*FAC*QFCN(4,1)*QFCN(3,2)
          CALL SIGFIL(SIG,4,3,28,30)
          SIG=.5*SIGHW*FAC*QFCN(9,1)*QFCN(6,2)
          CALL SIGFIL(SIG,9,6,28,30)
          SIG=.5*SIGHW*FAC*QFCN(6,1)*QFCN(9,2)
          CALL SIGFIL(SIG,6,9,28,30)
100       CONTINUE
      END IF
C
      IF (GOQ(30,1).AND.GOQ(28,2)) THEN
          CALL TWOKIN(0.,0.,AMH,AMW)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 110
          E1=SQRT(P(1)**2+AMH**2)
          E2=SQRT(P(2)**2+AMW**2)
          FAC=1./(12.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGHW=GF**2*AMW**8*(S/AMW2+(1.-T/AMW2)*(1.-U/AMW2))/
     $          PROPW*TBRWW(3,2)
          SIG=.5*SIGHW*FAC*QFCN(3,1)*QFCN(4,2)
          CALL SIGFIL(SIG,3,4,30,28)
          SIG=.5*SIGHW*FAC*QFCN(4,1)*QFCN(3,2)
          CALL SIGFIL(SIG,4,3,30,28)
          SIG=.5*SIGHW*FAC*QFCN(9,1)*QFCN(6,2)
          CALL SIGFIL(SIG,9,6,30,28)
          SIG=.5*SIGHW*FAC*QFCN(6,1)*QFCN(9,2)
          CALL SIGFIL(SIG,6,9,30,28)
110       CONTINUE
      END IF
C
C
C          WH production via W+*
C
      IF (GOQ(27,1).AND.GOQ(30,2)) THEN
          CALL TWOKIN(0.,0.,AMW,AMH)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 120
          E1=SQRT(P(1)**2+AMW**2)
          E2=SQRT(P(2)**2+AMH**2)
          FAC=1./(12.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGHW=GF**2*AMW**8*(S/AMW2+(1.-T/AMW2)*(1.-U/AMW2))/
     $          PROPW*TBRWW(2,1)
          SIG=.5*SIGHW*FAC*QFCN(2,1)*QFCN(5,2)
          CALL SIGFIL(SIG,2,5,27,30)
          SIG=.5*SIGHW*FAC*QFCN(5,1)*QFCN(2,2)
          CALL SIGFIL(SIG,5,2,27,30)
          SIG=.5*SIGHW*FAC*QFCN(8,1)*QFCN(7,2)
          CALL SIGFIL(SIG,8,7,27,30)
          SIG=.5*SIGHW*FAC*QFCN(7,1)*QFCN(8,2)
          CALL SIGFIL(SIG,7,8,27,30)
120       CONTINUE
      END IF
C
      IF (GOQ(30,1).AND.GOQ(27,2)) THEN
          CALL TWOKIN(0.,0.,AMH,AMW)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 130
          E1=SQRT(P(1)**2+AMH**2)
          E2=SQRT(P(2)**2+AMW**2)
          FAC=1./(12.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGHW=GF**2*AMW**8*(S/AMW2+(1.-T/AMW2)*(1.-U/AMW2))/
     $          PROPW*TBRWW(2,2)
          SIG=.5*SIGHW*FAC*QFCN(2,1)*QFCN(5,2)
          CALL SIGFIL(SIG,2,5,30,27)
          SIG=.5*SIGHW*FAC*QFCN(5,1)*QFCN(2,2)
          CALL SIGFIL(SIG,5,2,30,27)
          SIG=.5*SIGHW*FAC*QFCN(8,1)*QFCN(7,2)
          CALL SIGFIL(SIG,8,7,30,27)
          SIG=.5*SIGHW*FAC*QFCN(7,1)*QFCN(8,2)
          CALL SIGFIL(SIG,7,8,30,27)
130       CONTINUE
      END IF
C
C          ZH production via Z*
C          
      IF (GOQ(29,1).AND.GOQ(30,2)) THEN
          CALL TWOKIN(0.,0.,AMZ,AMH)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 200
          E1=SQRT(P(1)**2+AMZ2)
          E2=SQRT(P(2)**2+AMH2)
          FAC=1./(3.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 210 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 210
            SIG=GF**2*AMZ**8*(GV(IFLQ)**2+GA(IFLQ)**2)*
     $     (S/AMZ2+(1.-T/AMZ2)*(1.-U/AMZ2))/PROPZ*TBRWW(4,1)
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,29,30)
210         CONTINUE
200       CONTINUE
      END IF
C          HZ production via Z*
C          
      IF (GOQ(30,1).AND.GOQ(29,2)) THEN
          CALL TWOKIN(0.,0.,AMH,AMZ)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 220
          E1=SQRT(P(1)**2+AMH2)
          E2=SQRT(P(2)**2+AMZ2)
          FAC=1./(3.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 230 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 230
            SIG=GF**2*AMZ**8*(GV(IFLQ)**2+GA(IFLQ)**2)*
     $     (S/AMZ2+(1.-T/AMZ2)*(1.-U/AMZ2))/PROPZ*TBRWW(4,2)
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,30,29)
230         CONTINUE
220       CONTINUE
      END IF
      RETURN
      END
