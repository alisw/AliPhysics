#include "isajet/pilot.h"
      SUBROUTINE SIGTC
C
C          Compute the integrated technirho cross section
C          d(sigma)/d(qmw**2)d(yw) = d(sigma)/d(qmw**2)*f(x1)*f(x2)/scm
C          including W-technirho mixing from EHLQ 6.22 and 6.23 and 
C          elastic resonance in longitudinal WW fusion.
C
C          Use WTYPE for control with
C          WTYPE = 2     3     4
C                  rho+  rho-  rho0
C
C          SIGMA    = cross section summed over allowed types.
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4.
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 + I1
C                     using JETTYPE code.
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
#include "isajet/qsave.inc"
#include "isajet/wcon.inc"
#include "isajet/const.inc"
#include "isajet/jetlim.inc"
#include "isajet/hcon.inc"
#include "isajet/tcpar.inc"
C
      REAL AMQCUR(6),WTHELI(4),FINT(9),X(2)
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT),(X(1),X1)
      INTEGER MATCHT(4,4)
      REAL ACOSH,Z,ATANH,AMASS,QMW2,QMZ,EHAT,ANEFF,Q2SAVE,YHAT,EY,AMW,
     $AMZ,STRUC,STRUCW,WM,ZM,PWWCM,SIG0,S,T,U,FACINV,RATZ,Q1L,Q1R,SIG1,
     $SIG,QZW
      INTEGER I,IH,IQ,IW,IQ1,IQ2,IQ3,IQ4,IRHO,LISTW(4)
C
      DATA AMQCUR/.005,.009,.175,1.25,4.50,30./
      DATA LISTW/10,80,-80,90/
      DATA MATCHT/0,0,0,0, 0,29,0,27, 0,0,29,28, 0,28,27,0/
C
C          Functions
      ACOSH(Z)=ALOG(Z+SQRT(Z**2-1.))
      ATANH(Z)=.5*ALOG((1.+Z)/(1.-Z))
C
C          Kinematics (identical to Drell-Yan)
C
      AMQCUR(6)=AMASS(6)
      QMW2=QMW**2
      QTMW=SQRT(QMW2+QTW**2)
      Q0W=QTMW*COSH(YW)
      QZW=QTMW*SINH(YW)
      QW=SQRT(QZW**2+QTW**2)
      IF(QW.NE.0.) THEN
        CTHW=QZW/QW
        STHW=QTW/QW
        IF(ABS(CTHW).LT.1.) THEN
          THW=ACOS(CTHW)
        ELSE
          CTHW=0.
          STHW=1.
          THW=.5*PI
        ENDIF
      ELSE
        CTHW=0.
        STHW=1.
        THW=.5*PI
      ENDIF
      EHAT=QMW
      SHAT=QMW**2
      QSQ=SHAT
      ANEFF=4.+QSQ/(QSQ+AMASS(5)**2)+QSQ/(QSQ+AMASS(6)**2)
      ALFQSQ=12.*PI/((33.-ANEFF)*ALOG(QSQ/ALAM2))
      Q2SAVE=QSQ
      YHAT=YW
      EY=EXP(YHAT)
      X1=EHAT/ECM*EY
      X2=EHAT/(ECM*EY)
C
C          Initialize
C
      SIGMA=0.
      NSIGS=0
      DO 100 I=1,MXSIGS
100   SIGS(I)=0
C
      IF(X1.GE.1..OR.X2.GE.1.) RETURN
      AMW=WMASS(2)
      AMZ=WMASS(4)
C
C          Compute structure functions
C
      DO 110 IH=1,2
        DO 120 IQ=1,13
120     QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
        DO 130 IQ=14,26
130     QSAVE(IQ,IH)=0.
        DO 140 IW=2,4
          AMW=AMASS(LISTW(IW))
          IF(QMW.GT.2.*AMW) THEN
            QSAVE(25+IW,IH)=STRUCW(X(IH),IW,IDIN(IH))/X(IH)
          ELSE
            QSAVE(25+IW,IH)=0.
          ENDIF
140     CONTINUE
110   CONTINUE
C
C          qk + qb --> technirho0
C
      IF(.NOT.((GOQ(27,1).AND.GOQ(28,2)).OR.(GOQ(28,1).AND.GOQ(27,2))))
     $GO TO 300
      WM=WMASS(2)
      ZM=WMASS(4)
      IF(QMW.LE.2.*AMW) GO TO 300
      PWWCM=.5*SQRT(QMW**2-4.*WM**2)
      SIG0=PI*ALFA**2/(72.*SIN2W*S)*(2.*PWWCM/QMW)**3*X1*X2*UNITS
      SIG0=SIG0*TCMRHO**2/((S-TCMRHO**2)**2+TCMRHO**2*TCGRHO**2)
C          Initial state sum
      DO 210 IQ1=2,13
        IQ2=MATCH(IQ1,4)
        IF(IQ2.EQ.0) GO TO 210
        FACINV=2.*SQRT(SIN2W*(1.-SIN2W))
        RATZ=S/(S-ZM**2)
        Q1L=AQ(IQ1/2,4)*FACINV
        Q1R=BQ(IQ1/2,4)*FACINV
        SIG1=.25*SIG0*(1.-RATZ*Q1L/(Q1R*(1.-SIN2W))
     $  +RATZ**2*(Q1L**2+Q1R**2)/(4.*(1-SIN2W)**2))
     $  *QSAVE(IQ1,1)*QSAVE(IQ2,2)
C          Final state sum
        DO 220 IQ3=27,28
          IQ4=MATCHT(IQ3-25,4)
          IF(GOQ(IQ3,1).AND.GOQ(IQ4,2)) THEN
            SIG=SIG1*TBRWW(IQ3-25,1)*TBRWW(IQ4-25,2)
            CALL SIGFIL(SIG,IQ1,IQ2,IQ3,IQ4)
          ENDIF
220     CONTINUE
210   CONTINUE
C
C          W+ + W- -> technirho0 -> W+ + W-
C
      SIG0=12*PI/PWWCM**2*TCGRHO**2*X1*X2*UNITS
     $/((S-TCMRHO**2)**2+TCMRHO**2*TCGRHO**2)
C          Initial state sum
      DO 230 IQ1=27,28
        IQ2=MATCHT(IQ1-25,4)
        SIG1=.25*SIG0*QSAVE(IQ1,1)*QSAVE(IQ2,2)
C          Final state sum
        DO 240 IQ3=27,28
          IQ4=MATCHT(IQ3-25,4)
          IF(GOQ(IQ3,1).AND.GOQ(IQ4,2)) THEN
            SIG=SIG1*TBRWW(IQ3-25,1)*TBRWW(IQ4-25,2)
            CALL SIGFIL(SIG,IQ1,IQ2,IQ3,IQ4)
          ENDIF
240     CONTINUE
230   CONTINUE
C
C          q + qbar -> technirho+-
C
300   IF(.NOT.((GOQ(27,1).AND.GOQ(29,2)).OR.(GOQ(28,1).AND.GOQ(29,2))
     $.OR.(GOQ(29,1).AND.GOQ(27,2)).OR.(GOQ(29,1).AND.GOQ(28,2))))
     $GO TO 400
      WM=WMASS(2)
      ZM=WMASS(4)
      IF(QMW.LE.WM+ZM) GO TO 400
      PWWCM=SQRT((S-WM**2-ZM**2)**2-4.*WM**2*ZM**2)/(2.*QMW)
      SIG0=PI*ALFA**2/(144.*SIN2W)*S/(S-WM**2)**2*(2.*PWWCM/QMW)**3
     $*X1*X2*UNITS
      SIG0=SIG0*TCMRHO**2/((S-TCMRHO**2)**2+TCMRHO**2*TCGRHO**2)
      DO 310 IRHO=2,3
C          Initial state sum
        DO 320 IQ1=2,13
          IQ2=MATCH(IQ1,IRHO)
          IF(IQ2.EQ.0) GO TO 320
          SIG1=.25*SIG0*QSAVE(IQ1,1)*QSAVE(IQ2,2)
C          Final state sum
          DO 330 IQ3=27,28
            IQ4=MATCHT(IQ3-25,IRHO)
            IF(IQ4.EQ.0) GO TO 330
            IF(GOQ(IQ3,1).AND.GOQ(IQ4,2)) THEN
              SIG=SIG1*TBRWW(IQ3-25,1)*TBRWW(IQ4-25,2)
              CALL SIGFIL(SIG,IQ1,IQ2,IQ3,IQ4)
            ENDIF
330       CONTINUE
320     CONTINUE
310   CONTINUE
C
C          W+- + Z0 -> technirho+- -> W+- + Z0
C
      SIG0=12*PI/PWWCM**2*TCGRHO**2*X1*X2*UNITS
     $/((S-TCMRHO**2)**2+TCMRHO**2*TCGRHO**2)
      DO 340 IRHO=2,3
C          Initial state sum
        DO 350 IQ1=27,29
          IQ2=MATCHT(IQ1-25,IRHO)
          IF(IQ2.EQ.0) GO TO 350
          SIG1=.25*SIG0*QSAVE(IQ1,1)*QSAVE(IQ2,2)
C          Final state sum
          DO 360 IQ3=27,29
            IQ4=MATCHT(IQ3-25,IRHO)
            IF(IQ4.EQ.0) GO TO 360
            IF(GOQ(IQ3,1).AND.GOQ(IQ4,2)) THEN
              SIG=SIG1*TBRWW(IQ3-25,1)*TBRWW(IQ4-25,2)
              CALL SIGFIL(SIG,IQ1,IQ2,IQ3,IQ4)
            ENDIF
360       CONTINUE
350     CONTINUE
340   CONTINUE
C
400   RETURN
      END
