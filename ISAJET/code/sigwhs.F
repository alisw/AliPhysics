#include "isajet/pilot.h"
      SUBROUTINE SIGWHS
C
C          Calculate d(sigma)/d(pt**2)d(y1)d(y2) for 
C          Wh, WH, Zh, ZH, hA, HA and H+H- production in SUSY
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
#include "isajet/sspar.inc"
C
      REAL X(2)
      EQUIVALENCE (X(1),X1)
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT)
      REAL SIG,S,T,U,FAC,AMW,AMZ,AMW2,AMZ2,E1,E2,EQ1
      REAL QFCN,STRUC,SIGHW,SCFAC,BETA,SINW,COS2W
      REAL PROPZ,PROPW,GV(2),GA(2),AMH,GAMW,GAMZ
      INTEGER IS2UD(25),IQ,IH,I,IQ1,IQ2,IFLQ
      SAVE IS2UD
C
C          IS2UD: Susy jettype -> u/d code
      DATA IS2UD/0,1,1,2,2,2,2,1,1,2,2,1,1,1,1,2,2,2,2,1,1,2,2,1,1/

C          Functions
      QFCN(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
C
C          Initialize
      DO 10 I=1,MXSIGS
10    SIGS(I)=0.
      SIGMA=0.
      NSIGS=0
C
      BETA=ATAN(1./RV2V1)
      AMW=WMASS(2)
      AMW2=AMW**2
      AMZ=WMASS(4)
      AMZ2=AMZ**2
      GAMW=WGAM(2)
      GAMZ=WGAM(4)
      GV(1)=.25-2*SIN2W/3.
      GV(2)=-.25+SIN2W/3.
      GA(1)=-.25
      GA(2)=.25
      SINW=SQRT(SIN2W)
      THW=ASIN(SINW)
      COS2W=COS(2*THW)
      DO IH=81,82
      IF (IH.EQ.81) THEN
        SCFAC=SIN(ALFAH+BETA)**2
        AMH=AMHL
      ELSE
        SCFAC=COS(ALFAH+BETA)**2
        AMH=AMHH
      END IF
C
C          Wh, WH production via W-*
C
      IF (GOQ(79,1).AND.GOQ(IH,2)) THEN
          CALL TWOKIN(0.,0.,AMW,AMH)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 100
          E1=SQRT(P(1)**2+AMW**2)
          E2=SQRT(P(2)**2+AMH**2)
          FAC=1./(12.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGHW=GF**2*AMW**8*(S/AMW2+(1.-T/AMW2)*(1.-U/AMW2))/
     $          PROPW*TBRWW(3,1)*SCFAC
          SIG=.5*SIGHW*FAC*QFCN(3,1)*QFCN(4,2)
          CALL SIGFIL(SIG,3,4,79,IH)
          SIG=.5*SIGHW*FAC*QFCN(4,1)*QFCN(3,2)
          CALL SIGFIL(SIG,4,3,79,IH)
          SIG=.5*SIGHW*FAC*QFCN(9,1)*QFCN(6,2)
          CALL SIGFIL(SIG,9,6,79,IH)
          SIG=.5*SIGHW*FAC*QFCN(6,1)*QFCN(9,2)
          CALL SIGFIL(SIG,6,9,79,IH)
100       CONTINUE
      END IF
C
      IF (GOQ(IH,1).AND.GOQ(79,2)) THEN
          CALL TWOKIN(0.,0.,AMH,AMW)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 110
          E1=SQRT(P(1)**2+AMH**2)
          E2=SQRT(P(2)**2+AMW**2)
          FAC=1./(12.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGHW=GF**2*AMW**8*(S/AMW2+(1.-T/AMW2)*(1.-U/AMW2))/
     $          PROPW*TBRWW(3,2)*SCFAC
          SIG=.5*SIGHW*FAC*QFCN(3,1)*QFCN(4,2)
          CALL SIGFIL(SIG,3,4,IH,79)
          SIG=.5*SIGHW*FAC*QFCN(4,1)*QFCN(3,2)
          CALL SIGFIL(SIG,4,3,IH,79)
          SIG=.5*SIGHW*FAC*QFCN(9,1)*QFCN(6,2)
          CALL SIGFIL(SIG,9,6,IH,79)
          SIG=.5*SIGHW*FAC*QFCN(6,1)*QFCN(9,2)
          CALL SIGFIL(SIG,6,9,IH,79)
110       CONTINUE
      END IF
C
C
C          Wh, WH production via W+*
C
      IF (GOQ(78,1).AND.GOQ(IH,2)) THEN
          CALL TWOKIN(0.,0.,AMW,AMH)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 120
          E1=SQRT(P(1)**2+AMW**2)
          E2=SQRT(P(2)**2+AMH**2)
          FAC=1./(12.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGHW=GF**2*AMW**8*(S/AMW2+(1.-T/AMW2)*(1.-U/AMW2))/
     $          PROPW*TBRWW(2,1)*SCFAC
          SIG=.5*SIGHW*FAC*QFCN(2,1)*QFCN(5,2)
          CALL SIGFIL(SIG,2,5,78,IH)
          SIG=.5*SIGHW*FAC*QFCN(5,1)*QFCN(2,2)
          CALL SIGFIL(SIG,5,2,78,IH)
          SIG=.5*SIGHW*FAC*QFCN(8,1)*QFCN(7,2)
          CALL SIGFIL(SIG,8,7,78,IH)
          SIG=.5*SIGHW*FAC*QFCN(7,1)*QFCN(8,2)
          CALL SIGFIL(SIG,7,8,78,IH)
120       CONTINUE
      END IF
C
      IF (GOQ(IH,1).AND.GOQ(78,2)) THEN
          CALL TWOKIN(0.,0.,AMH,AMW)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 130
          E1=SQRT(P(1)**2+AMH**2)
          E2=SQRT(P(2)**2+AMW**2)
          FAC=1./(12.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPW=(S-AMW**2)**2+AMW**2*GAMW**2
          SIGHW=GF**2*AMW**8*(S/AMW2+(1.-T/AMW2)*(1.-U/AMW2))/
     $          PROPW*TBRWW(2,2)*SCFAC
          SIG=.5*SIGHW*FAC*QFCN(2,1)*QFCN(5,2)
          CALL SIGFIL(SIG,2,5,IH,78)
          SIG=.5*SIGHW*FAC*QFCN(5,1)*QFCN(2,2)
          CALL SIGFIL(SIG,5,2,IH,78)
          SIG=.5*SIGHW*FAC*QFCN(8,1)*QFCN(7,2)
          CALL SIGFIL(SIG,8,7,IH,78)
          SIG=.5*SIGHW*FAC*QFCN(7,1)*QFCN(8,2)
          CALL SIGFIL(SIG,7,8,IH,78)
130       CONTINUE
      END IF
C
C          Zh, ZH production via Z*
C          
      IF (GOQ(80,1).AND.GOQ(IH,2)) THEN
          CALL TWOKIN(0.,0.,AMZ,AMH)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 200
          E1=SQRT(P(1)**2+AMZ2)
          E2=SQRT(P(2)**2+AMH**2)
          FAC=1./(3.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 210 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 210
            SIG=GF**2*AMZ**8*(GV(IFLQ)**2+GA(IFLQ)**2)*
     $     (S/AMZ2+(1.-T/AMZ2)*(1.-U/AMZ2))/PROPZ*TBRWW(4,1)*SCFAC
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,80,IH)
210         CONTINUE
200       CONTINUE
      END IF
C          hZ, HZ production via Z*
C          
      IF (GOQ(IH,1).AND.GOQ(80,2)) THEN
          CALL TWOKIN(0.,0.,AMH,AMZ)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 220
          E1=SQRT(P(1)**2+AMH**2)
          E2=SQRT(P(2)**2+AMZ2)
          FAC=1./(3.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 230 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 230
            SIG=GF**2*AMZ**8*(GV(IFLQ)**2+GA(IFLQ)**2)*
     $     (S/AMZ2+(1.-T/AMZ2)*(1.-U/AMZ2))/PROPZ*TBRWW(4,2)*SCFAC
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,IH,80)
230         CONTINUE
220       CONTINUE
      END IF
C
C     Next, do Ah and AH production
C
      IF (GOQ(83,1).AND.GOQ(IH,2)) THEN
          CALL TWOKIN(0.,0.,AMHA,AMH)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 240
          E1=SQRT(P(1)**2+AMHA**2)
          E2=SQRT(P(2)**2+AMH**2)
          FAC=1./(12.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 250 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 250
            SIG=GF**2*AMZ**4*(GV(IFLQ)**2+GA(IFLQ)**2)*
     $          ((AMHA**2+U-T-AMH**2)*(AMHA**2+T-U-AMH**2)-
     $          S*(2*AMHA**2+2*AMH**2-S))/PROPZ*SCFAC
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,83,IH)
250         CONTINUE
240       CONTINUE
      END IF
      IF (GOQ(IH,1).AND.GOQ(83,2)) THEN
          CALL TWOKIN(0.,0.,AMH,AMHA)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 260
          E1=SQRT(P(1)**2+AMH**2)
          E2=SQRT(P(2)**2+AMHA**2)
          FAC=1./(12.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 270 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 270
            SIG=GF**2*AMZ**4*(GV(IFLQ)**2+GA(IFLQ)**2)*
     $          ((AMHA**2+U-T-AMH**2)*(AMHA**2+T-U-AMH**2)-
     $          S*(2*AMHA**2+2*AMH**2-S))/PROPZ*SCFAC
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,IH,83)
270         CONTINUE
260       CONTINUE
      END IF
      END DO
C
C     Next, do H+H- production
C
      IF (GOQ(84,1).AND.GOQ(85,2)) THEN
          CALL TWOKIN(0.,0.,AMHC,AMHC)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 300
          E1=SQRT(P(1)**2+AMHC**2)
          E2=SQRT(P(2)**2+AMHC**2)
          FAC=1./(96.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 310 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IFLQ.EQ.1) THEN
              EQ1=2./3.
            ELSE
              EQ1=-1./3.
            END IF
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 310
            SIG=((4*PI*ALFA)**2*EQ1**2/S/S+32*PI*ALFA*EQ1*GF*AMZ**2*
     $          COS2W*GV(IFLQ)*(S-AMZ**2)/S/PROPZ/SQRT2+8*GF**2*
     $          AMZ**4*COS2W**2*(GV(IFLQ)**2+GA(IFLQ)**2)/PROPZ)*
     $          ((U-T)*(T-U)-S*(4*AMHC**2-S))
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,84,85)
310         CONTINUE
300       CONTINUE
      END IF
      IF (GOQ(85,1).AND.GOQ(84,2)) THEN
          CALL TWOKIN(0.,0.,AMHC,AMHC)
          IF(X1.GE.1..OR.X2.GE.1.) GO TO 320
          E1=SQRT(P(1)**2+AMHC**2)
          E2=SQRT(P(2)**2+AMHC**2)
          FAC=1./(96.*PI*S**2)
          FAC=FAC*S/SCM*(P(1)*P(2)/(E1*E2))*UNITS
          PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
          DO 330 IQ1=2,11
            IFLQ=IS2UD(IQ1)
            IQ2=MATCH(IQ1,4)
            IF (IFLQ.EQ.1) THEN
              EQ1=2./3.
            ELSE
              EQ1=-1./3.
            END IF
            IF (IQ2.EQ.0.OR.IQ2.GE.12) GO TO 330
            SIG=((4*PI*ALFA)**2*EQ1**2/S/S+32*PI*ALFA*EQ1*GF*AMZ**2*
     $          COS2W*GV(IFLQ)*(S-AMZ**2)/S/PROPZ/SQRT2+8*GF**2*
     $          AMZ**4*COS2W**2*(GV(IFLQ)**2+GA(IFLQ)**2)/PROPZ)*
     $          ((U-T)*(T-U)-S*(4*AMHC**2-S))
            SIG=.5*SIG*FAC*QFCN(IQ1,1)*QFCN(IQ2,2)
            CALL SIGFIL(SIG,IQ1,IQ2,85,84)
330         CONTINUE
320       CONTINUE
      END IF
      RETURN
      END
