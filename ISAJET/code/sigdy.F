#include "isajet/pilot.h"
      SUBROUTINE SIGDY
C
C          Compute the Drell-Yan and Drell-Yan plus jet cross sections
C          d(sigma)/d(qmw**2)d(qtw**2)d(yw)d(yj)
C
C          SIGMA    = cross section summed over quark types allowed by
C                     JETTYPE3 and WTYPE cards.
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4.
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 + I1
C                     using JETTYPE code.
C
C          QT cutoff for W+JET taken from Parisi and Petronzio,
C          Nucl Phys B154, 427
C          qk + gl --> qk + w suppressed at low QTW by extra factor
C          of qtw**2/(qtw**2+qt2cut(qmw))
C
C          Ver 7.17: include top mass for gb --> Wt and gt --> Zt 
C          with no extra qt suppression factor. Note we do NOT include
C          gt --> Wb; while this process makes sense for qt >> m_t,
C          it has a pole in the physical region at low qt from the 
C          on-shell decay t --> Wb. We let Q**2 --> Q**2 + m_t**2 
C          in the scale for the parton distributions.
C
C          Ver 7.32: Rewrite AJLWT for gb --> Wt, etc., in terms of 
C          scaled variables, and restore SWT**5 later to avoid 
C          floating errors on VMS.
C
C          Ver 7.41: Recalculate COUT for each mass(!).
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
#include "isajet/nodcay.inc"
C
      REAL X(2)
      REAL Z,S,T,U,QMW2,QZW,EHAT,Q2SAVE,YHAT,EY,P3Z,P1,P2,AMASS,ANEFF,
     $SIG0,DENOM,QT2CUT,SIGT,SIGU,FAC,PROP,FACTOR,SIG,AMT,AMT2,SWT,
     $P1WT,P2WT,X1WT,X2WT,TWT,UWT,Q2,QFCN,STRUC,XX,ACOSH,ATANH,P2M,P1M
      REAL AMI2,AMF2,EFWT
      REAL AJLWT,AJLZT1,AJLZT2,A2,A2B2,QQ,TM2
      INTEGER I,IQ,IH,IQ1,IFL,IQ2,IW
      INTEGER NZERO(4)
      REAL AMFAC(13)
      INTEGER NUTYP(25)
      INTEGER IFL1,IFL2
      REAL TERM
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT),(X(1),X1)
C
      DATA NZERO/11,9,9,11/
      DATA AMFAC/11*0.,2*1./
      DATA NUTYP/13*0,1,1,0,0,1,1,0,0,1,1,0,0/
C
C          Functions
      ACOSH(Z)=ALOG(Z+SQRT(Z**2-1.))
      ATANH(Z)=.5*ALOG((1.+Z)/(1.-Z))
      PROP(I)=1./((QMW2-WMASS(I)**2)**2+(WMASS(I)*WGAM(I))**2)
C          Qt cutoff function
      QT2CUT(QMW)=CUTOFF*QMW**CUTPOW
C          Parton distributions
      QFCN(XX,IQ,IH)=STRUC(XX,QSQ+AMT2,IQ,IDIN(IH))/XX
C          Integrated matrix elements JLint from FORM
      AJLWT(S,T,QQ,TM2)=
     $ - 32*QQ**3*S*T + 32*QQ**3*S*TM2 + 32*QQ**2*S**2*T 
     $ + 32*QQ**2*S*T**2 - 16*QQ**2*S*T*TM2 - 16*QQ**2*S*TM2**2 
     $ - 16*QQ*S**3*T + 16*QQ*S**3*TM2 - 16*QQ*S**2*T*TM2
     $ - 16*QQ*S*T**3 + 32*QQ*S*T**2*TM2 - 16*QQ*S*T*TM2**2 
     $ - 8*S**3*T*TM2 + 8*S**3*TM2**2 - 16*S**2*T**2*TM2 
     $ + 16*S**2*T*TM2**2 - 16*S**2*TM2**3 - 8*S*T**3*TM2 
     $ + 8*S*T**2*TM2**2 - 8*S*T*TM2**3 + 8*S*TM2**4
C
      AJLZT1(S,T,QQ,TM2)=
     $ + A2 * ( - 96*QQ**2*S*T*TM2 + 96*QQ**2*S*TM2**2
     $ + 96*QQ**2*T*TM2**2 - 96*QQ**2*TM2**3 + 96*QQ*S**2*T*TM2
     $ + 96*QQ*S*T**2*TM2 - 192*QQ*S*T*TM2**2 - 96*QQ*S*TM2**3
     $ - 96*QQ*T*TM2**3 + 192*QQ*TM2**4 + 16*S**3*T*TM2
     $ - 16*S**3*TM2**2 + 32*S**2*T**2*TM2 - 112*S**2*T*TM2**2
     $ + 80*S**2*TM2**3 + 16*S*T**3*TM2 - 112*S*T**2*TM2**2
     $ + 224*S*T*TM2**3 - 128*S*TM2**4 - 16*T**3*TM2**2
     $ + 80*T**2*TM2**3 - 128*T*TM2**4 + 64*TM2**5 )
      AJLZT2(S,T,QQ,TM2)=
     $ + A2B2 * ( - 16*QQ**3*S*T + 16*QQ**3*S*TM2 + 16*QQ**3*T*TM2
     $ - 16*QQ**3*TM2**2 + 16*QQ**2*S**2*T + 16*QQ**2*S*T**2
     $ + 32*QQ**2*S*T*TM2 - 80*QQ**2*S*TM2**2 - 80*QQ**2*T*TM2**2
     $ + 96*QQ**2*TM2**3 - 8*QQ*S**3*T + 8*QQ*S**3*TM2
     $ - 40*QQ*S**2*T*TM2 - 24*QQ*S**2*TM2**2 - 8*QQ*S*T**3
     $ - 40*QQ*S*T**2*TM2 + 80*QQ*S*T*TM2**2 + 96*QQ*S*TM2**3
     $ + 8*QQ*T**3*TM2 - 24*QQ*T**2*TM2**2 + 96*QQ*T*TM2**3
     $ - 144*QQ*TM2**4 - 16*S**3*T*TM2 + 16*S**3*TM2**2
     $ - 32*S**2*T**2*TM2 + 112*S**2*T*TM2**2 - 80*S**2*TM2**3
     $ - 16*S*T**3*TM2 + 112*S*T**2*TM2**2 - 224*S*T*TM2**3
     $ + 128*S*TM2**4 + 16*T**3*TM2**2 - 80*T**2*TM2**3
     $ + 128*T*TM2**4 - 64*TM2**5 )
C
C          Kinematics
C
      QMW2=QMW**2
      QTMW=SQRT(QMW2+QTW**2)
      Q0W=QTMW*COSH(YW)
      QZW=QTMW*SINH(YW)
      QW=SQRT(QZW**2+QTW**2)
C          Protect against errors
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
C
      IF(STDDY) THEN
C          Kinematics for standard Drell-Yan
        EHAT=QMW
        SHAT=QMW**2
        QSQ=SHAT
        Q2SAVE=QSQ
        YHAT=YW
        EY=EXP(YHAT)
        X1=EHAT/ECM*EY
        X2=EHAT/(ECM*EY)
      ELSE
C          Kinematics for Drell-Yan plus jet
        P3Z=P(3)*CTH(3)
        SHAT=QMW2+2.*Q0W*P(3)-2.*QZW*P3Z+2.*PT(3)**2
        P1=.5*(P(3)+P3Z+Q0W+QZW)
        P2=.5*(P(3)-P3Z+Q0W-QZW)
        X1=P1/HALFE
        X2=P2/HALFE
        THAT=-2.*P1*(P(3)-P3Z)
        UHAT=-2.*P2*(P(3)+P3Z)
        QSQ=QTW**2
        QSQ=AMAX1(QSQ,4.)
        ANEFF=4.+QSQ/(QSQ+AMASS(5)**2)+QSQ/(QSQ+AMASS(6)**2)
        ALFQSQ=12.*PI/((33.-2.*ANEFF)*ALOG(QSQ/ALAM2))
        Q2SAVE=QSQ
        QSQ=SHAT
      ENDIF
C
C          Initialize
C
      SIGMA=0.
      NSIGS=0
      DO 100 I=1,MXSIGS
        SIGS(I)=0.
100   CONTINUE
      IF(X1.GE.1..OR.X2.GE.1.) RETURN
C
C          Compute structure functions
C
      DO 110 IH=1,2
        DO 120 IQ=1,11
          QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
120     CONTINUE
        QSAVE(12,IH)=0
        QSAVE(13,IH)=0
110   CONTINUE
      QSQ=Q2SAVE
C
C          Recompute COUT for this mass
C
      DO 130 IW=1,4
        COUT(IW)=0.
        IF(.NOT.GODY(IW)) GO TO 130
        DO 140 IQ1=2,25
          IQ2=MATCH(IQ1,IW)
          IF(IQ2.EQ.0) GO TO 140
          IF(.NOT.(GOQ(IQ1,1).AND.GOQ(IQ2,2))) GO TO 140
          IF(NUTYP(IQ1)*NUTYP(IQ2).EQ.1.AND.NONUNU) GO TO 140
          IFL1=IQ1/2
          IFL2=IQ2/2
          IF(AMASS(IFL1)+AMASS(IFL2).GE.QMW) GO TO 140
          TERM=.5*(AQ(IFL1,IW)**2+BQ(IFL1,IW)**2)
          IF(IQ1.LE.13) TERM=3.*TERM
          COUT(IW)=COUT(IW)+TERM
140     CONTINUE
130   CONTINUE
C
      IF(STDDY) GO TO 400
C
C          Compute cross section for types allowed by WTYPE and
C          JETTYPE cards.
C
C          qk + gl --> qk + W
C
      SIG0=ALFA**2*ALFQSQ*QMW2/(9.*SCM*S)*UNITS
      DENOM=S**2*EXP(.5*ALOG(QTW**4+QT2CUT(QMW)**2))
      SIGT=SIG0*(S**2+U**2+2.*T*QMW2)*(-T)/DENOM
      SIGU=SIG0*(S**2+T**2+2.*U*QMW2)*(-U)/DENOM
      DO 200 IW=1,4
        IF(.NOT.GODY(IW)) GO TO 200
        FAC=COUT(IW)*PROP(IW)
        DO 210 IQ=2,NZERO(IW)
          IF(.NOT.GOQ(IQ,3)) GO TO 210
          IQ1=MATCH(IQ,4)
          IQ1=MATCH(IQ1,IW)
          IF(IQ1.EQ.0) GO TO 210
          IFL=IQ/2
          FACTOR=FAC*(AQ(IFL,IW)**2+BQ(IFL,IW)**2)
     $    *QTW**2/(QTW**2+QT2CUT(QMW))
          SIG=FACTOR*SIGT*QSAVE(IQ1,1)*QSAVE(1,2)
          CALL SIGFIL(SIG,IQ1,1,IW,IQ)
          SIG=FACTOR*SIGU*QSAVE(IQ1,2)*QSAVE(1,1)
          CALL SIGFIL(SIG,1,IQ1,IW,IQ)
210     CONTINUE
200   CONTINUE
C
C          bt,tp + gl -> bt,tp + W,Z
C
      AMT=AMASS(6)
      AMT2=AMT**2
      Q2=QMW2
      DO 220 IW=2,4
        IF(.NOT.GODY(IW)) GO TO 220
        DO 230 IQ=NZERO(IW)+1,13
          IF(.NOT.GOQ(IQ,3)) GO TO 230
          IQ1=MATCH(IQ,4)
          IQ1=MATCH(IQ1,IW)
          IF(IQ1.EQ.0) GO TO 230
          IF(IQ1.GE.12.AND.IW.NE.4) GO TO 230
C          Assign zero or top masses for initial/final quarks
          AMF2=AMT2*AMFAC(IQ)
          AMI2=AMT2*AMFAC(IQ1)
          EFWT=SQRT(P(3)**2+AMF2)
          SWT=QMW2+AMF2+2.*Q0W*EFWT-2.*QZW*P3Z+2.*PT(3)**2
C
C          qk + gl initial state
C          Do kinematics using p(small) = 0 for gluon
C
          P1WT=EFWT+P3Z+Q0W+QZW
          P1M=AMI2/P1WT
          P2WT=EFWT-P3Z+Q0W-QZW-P1M
          X1WT=.5*P1WT/HALFE
          X2WT=.5*P2WT/HALFE
          TWT=-P1WT*(EFWT-P3Z)-P1M*(EFWT+P3Z)+AMI2+AMF2
          UWT=-P2WT*(EFWT+P3Z)+AMF2
          IF(X1WT.LT.0.OR.X1WT.GT.1.OR.X2WT.LT.0.OR.X2WT.GT.1)
     $    GO TO 240
C          Cross sections
          IF(IW.EQ.2.OR.IW.EQ.3) THEN
            SIG0=ALFA**2*ALFQSQ/(144*SCM*SWT)*UNITS
            SIG0=SIG0*(AQ(5,IW)**2+BQ(5,IW)**2)*COUT(IW)*PROP(IW)
            SIGU=SIG0*AJLWT(SWT/SWT,UWT/SWT,Q2/SWT,AMT2/SWT)*SWT*
     $      (SWT/(SWT-AMI2))**2*(SWT/(UWT-AMF2))**2
            SIG=SIGU*QFCN(X1WT,IQ1,1)*QFCN(X2WT,1,2)
            CALL SIGFIL(SIG,IQ1,1,IW,IQ)
          ELSEIF(IW.EQ.4) THEN
            SIG0=ALFA**2*ALFQSQ/(144*SCM*SWT)*UNITS
            SIG0=SIG0*COUT(IW)*PROP(IW)
            A2=AQ(6,IW)**2
            A2B2=AQ(6,IW)**2+BQ(6,IW)**2
            SIGU=SIG0*(AJLZT1(SWT/SWT,UWT/SWT,Q2/SWT,AMT2/SWT)+
     $      AJLZT2(SWT/SWT,UWT/SWT,Q2/SWT,AMT2/SWT))*SWT*
     $      (SWT/(SWT-AMI2))**2*(SWT/(UWT-AMF2))**2
            SIG=SIGU*QFCN(X1WT,IQ1,1)*QFCN(X2WT,1,2)
            CALL SIGFIL(SIG,IQ1,1,IW,IQ)
          ENDIF
C
C          gl + qk initial state
C          Do kinematics  using p(small) = 0 for gluon
C
240       P2WT=EFWT-P3Z+Q0W-QZW
          P2M=AMI2/P2WT
          P1WT=EFWT+P3Z+Q0W+QZW-P2M
          X1WT=.5*P1WT/HALFE
          X2WT=.5*P2WT/HALFE
          TWT=-P1WT*(EFWT-P3Z)+AMF2
          UWT=-P2WT*(EFWT+P3Z)-P2M*(EFWT-P3Z)+AMI2+AMF2
          IF(X1WT.LT.0.OR.X1WT.GT.1.OR.X2WT.LT.0.OR.X2WT.GT.1)
     $    GO TO 230
C          Cross sections
          IF(IW.EQ.2.OR.IW.EQ.3) THEN
            SIG0=ALFA**2*ALFQSQ/(144*SCM*SWT)*UNITS
            SIG0=SIG0*(AQ(5,IW)**2+BQ(5,IW)**2)*COUT(IW)*PROP(IW)
            SIGT=SIG0*AJLWT(SWT/SWT,TWT/SWT,Q2/SWT,AMT2/SWT)*SWT*
     $      (SWT/(SWT-AMI2))**2*(SWT/(TWT-AMF2)**2)
            SIG=SIGT*QFCN(X1WT,1,1)*QFCN(X2WT,IQ1,2)
            CALL SIGFIL(SIG,1,IQ1,IW,IQ)
          ELSEIF(IW.EQ.4) THEN
            SIG0=ALFA**2*ALFQSQ/(144*SCM*SWT)*UNITS
            SIG0=SIG0*COUT(IW)*PROP(IW)
            A2=AQ(6,IW)**2
            A2B2=AQ(6,IW)**2+BQ(6,IW)**2
            SIGU=SIG0*(AJLZT1(SWT/SWT,TWT/SWT,Q2/SWT,AMT2/SWT)+
     $      AJLZT2(SWT/SWT,TWT/SWT,Q2/SWT,AMT2/SWT))*SWT*
     $      (SWT/(SWT-AMI2))**2*(SWT/(TWT-AMF2))**2
            SIG=SIGU*QFCN(X1WT,1,1)*QFCN(X2WT,IQ1,2)
            CALL SIGFIL(SIG,1,IQ1,IW,IQ)
          ENDIF
230     CONTINUE
220   CONTINUE
C
C          qk + qb --> gl + W
C
      IF(.NOT.GOQ(1,3)) RETURN
      SIG0=8.*ALFA**2*ALFQSQ*QMW2/(27.*SCM*S)*UNITS
      DENOM=S*EXP(.5*ALOG(QTW**4+QT2CUT(QMW)**2))
      SIG0=SIG0*(T**2+U**2+2.*S*QMW2)/DENOM
      DO 300 IW=1,4
        IF(.NOT.GODY(IW)) GO TO 300
        FAC=COUT(IW)*PROP(IW)
        DO 310 IQ1=2,11
          IQ2=MATCH(IQ1,IW)
          IF(IQ2.EQ.0) GO TO 310
          IFL=IQ1/2
          SIG=FAC*SIG0*(AQ(IFL,IW)**2+BQ(IFL,IW)**2)
     $    *QSAVE(IQ1,1)*QSAVE(IQ2,2)
          CALL SIGFIL(SIG,IQ1,IQ2,IW,1)
310     CONTINUE
300   CONTINUE
      RETURN
C
C          Standard Drell-Yan for QT=0.
C
400   CONTINUE
      SIG0=4.*PI*ALFA**2*QMW2/(9.*SCM)*UNITS
      DO 410 IW=1,4
        IF(.NOT.GODY(IW)) GO TO 410
        FAC=COUT(IW)*PROP(IW)
        DO 420 IQ1=2,13
          IQ2=MATCH(IQ1,IW)
          IF(IQ2.EQ.0) GO TO 420
          IFL=IQ1/2
          SIG=FAC*SIG0*(AQ(IFL,IW)**2+BQ(IFL,IW)**2)
     $    *QSAVE(IQ1,1)*QSAVE(IQ2,2)
          CALL SIGFIL(SIG,IQ1,IQ2,IW,0)
420     CONTINUE
410   CONTINUE
C
      RETURN
      END
