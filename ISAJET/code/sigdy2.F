#include "isajet/pilot.h"
      SUBROUTINE SIGDY2
C
C          Compute the lepton-lepton-jet cross or quark-antiquark-jet
C          cross section
C          d(sigma)/d(qmw**2)d(qtw**2)d(yw)d(yj)d(omega*)
C          for the specified W and jet types
C
C          Also fix the incoming partons to be the selected types.
C
C          QT cutoff from Parisi and Petronzio, Nucl Phys B154, 427
C          qk+gl-->qk+w suppressed at low QTW
C
C          Ver 6.40: Fix underflow in standard Drell-Yan
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
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/jetpar.inc"
#include "isajet/primar.inc"
#include "isajet/q1q2.inc"
#include "isajet/jetsig.inc"
#include "isajet/wsig.inc"
#include "isajet/qsave.inc"
#include "isajet/wcon.inc"
#include "isajet/const.inc"
C
      REAL PROP,AJLWT,FCDIR,FCINT,QT2CUT,QFCN,AMASS,STRUC
      REAL AJLZT1,AJLZT2,AJLZT3,AJLZT4
      REAL X1WT,X2WT,P1WT,P2WT,SWT,TWT,UWT,QZW,P3Z,
     $S,T,U,S1,T1,U1,TDIR,UDIR,TINT,UINT,COUPL,DENOM,
     $AMT,AMT2,TERM1,TERM2,SIG0,QMW2,Q2,XX,S1WT,T1WT,U1WT,P1M,P2M,
     $AMI2,AMF2,EFWT,A2,A2B2,AB,AL2BL2,ALBL,QQ,TM2
      REAL AMFAC(13)
      INTEGER I,JF,IFLQ,JQK,IQ1,IQ2,IFL1,IQ,IFLL,IH,IQ3
      INTEGER NZERO(4)
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT)
      EQUIVALENCE (S1,SHAT1),(T1,THAT1),(U1,UHAT1)
C
      DATA NZERO/13,9,9,11/
      DATA AMFAC/11*0.,2*1./

C          Functions.  FCDIR and FCINT are direct and interference
C          terms for virtual Compton cross section.
C
      PROP(I)=1./((QMW2-WMASS(I)**2)**2+(WMASS(I)*WGAM(I))**2)
      FCDIR(S,T,U,S1,T1,U1)=COUPL*(S*(2.*S1**2-2.*S1*U1-2.*T1*U1)
     1+T*(-2.*S1*T1-4.*S1*U1-2.*T1*U1)+U*(2.*U1**2-2.*S1*T1-2.*S1*U1)
     2+QMW2*(S**2+U**2+2.*T*QMW2))*(-T)/DENOM
      FCINT(S,T,U,S1,T1,U1)=-COUPL*(S1*(S*T-T*U+T*QMW2+QMW2**2)
     1+T1*(-S**2+U**2+2.*S*QMW2-2.*U*QMW2)+U1*(S*T-T*U-T*QMW2-QMW2**2))
     2*(-T)/DENOM
C          QT cutoff function
      QT2CUT(QMW)=CUTOFF*QMW**CUTPOW
C          Parton distributions for top processes
      QFCN(XX,IQ,IH)=STRUC(XX,QSQ+AMT2,IQ,IDIN(IH))/XX          
C          Matrix elements JL/128 from FORM
      AJLWT(S,T,T1,U1,QQ,TM2)=
     $ + 2*QQ**3*S*T - 2*QQ**3*S*TM2 - 2*QQ**2*S**2*TM2
     $ - 2*QQ**2*S*T**2 + 4*QQ**2*S*T*T1 + 2*QQ**2*S*T*U1
     $ - 4*QQ**2*S*T1*TM2 - 2*QQ**2*S*U1*TM2 + 2*QQ**2*S*TM2**2
     $ + QQ*S**2*T*TM2 - 4*QQ*S**2*T1*TM2 + QQ*S**2*TM2**2
     $ + QQ*S*T**3 - 2*QQ*S*T**2*T1 - QQ*S*T**2*TM2 + 2*QQ*S*T*T1**2
     $ + 2*QQ*S*T*T1*U1 - 2*QQ*S*T*T1*TM2 + QQ*S*T*U1**2
     $ - 3*QQ*S*T*U1*TM2 + QQ*S*T*TM2**2 - 2*QQ*S*T1**2*TM2
     $ - 2*QQ*S*T1*U1*TM2 + 4*QQ*S*T1*TM2**2 - QQ*S*U1**2*TM2
     $ + 3*QQ*S*U1*TM2**2 - QQ*S*TM2**3 + S**2*T*T1*TM2
     $ - S**2*T*U1*TM2 - 2*S**2*T1**2*TM2 + S**2*T1*TM2**2
     $ + S**2*U1*TM2**2 + S*T**2*T1*TM2 - 2*S*T*T1**2*TM2
     $ - 2*S*T*T1*U1*TM2 - S*T*U1**2*TM2 + S*T*U1*TM2**2
     $ + 2*S*T1**2*TM2**2 + 2*S*T1*U1*TM2**2 - S*T1*TM2**3
     $ + S*U1**2*TM2**2 - S*U1*TM2**3
C
      AJLZT1(S,T,T1,U1,QQ,TM2)=
     $ + A2*AL2BL2 * ( 8*QQ**2*S*T*TM2 - 8*QQ**2*S*TM2**2
     $ - 8*QQ**2*T*TM2**2 + 8*QQ**2*TM2**3 - 8*QQ*S**2*T*TM2
     $ - 8*QQ*S*T**2*TM2 + 16*QQ*S*T*TM2**2 + 8*QQ*S*TM2**3
     $ + 8*QQ*T*TM2**3 - 16*QQ*TM2**4 + 8*S**2*T*U1*TM2
     $ - 8*S**2*U1*TM2**2 + 8*S*T**2*U1*TM2 + 8*S*T*U1**2*TM2
     $ - 32*S*T*U1*TM2**2 - 8*S*U1**2*TM2**2 + 24*S*U1*TM2**3
     $ - 8*T**2*U1*TM2**2 - 8*T*U1**2*TM2**2 + 24*T*U1*TM2**3
     $ + 8*U1**2*TM2**3 - 16*U1*TM2**4 )/8.
      AJLZT2(S,T,T1,U1,QQ,TM2)=
     $ + A2B2*AL2BL2 * ( 2*QQ**3*S*T - 2*QQ**3*S*TM2
     $ - 2*QQ**3*T*TM2 + 2*QQ**3*TM2**2 - 2*QQ**2*S**2*TM2
     $ - 2*QQ**2*S*T**2 + 4*QQ**2*S*T*T1 + 2*QQ**2*S*T*U1
     $ - 8*QQ**2*S*T*TM2 - 4*QQ**2*S*T1*TM2 - 2*QQ**2*S*U1*TM2
     $ + 14*QQ**2*S*TM2**2 - 4*QQ**2*T*T1*TM2 - 2*QQ**2*T*U1*TM2
     $ + 12*QQ**2*T*TM2**2 + 4*QQ**2*T1*TM2**2 + 2*QQ**2*U1*TM2**2
     $ - 14*QQ**2*TM2**3 + QQ*S**3*T - QQ*S**3*TM2 + 2*QQ*S**2*T*T1
     $ + 2*QQ*S**2*T*U1 - QQ*S**2*T*TM2 - 6*QQ*S**2*T1*TM2
     $ - 2*QQ*S**2*U1*TM2 + 9*QQ*S**2*TM2**2 + QQ*S*T**3
     $ - 2*QQ*S*T**2*T1 + 3*QQ*S*T**2*TM2 + 4*QQ*S*T*T1**2
     $ + 4*QQ*S*T*T1*U1 - 16*QQ*S*T*T1*TM2 + 2*QQ*S*T*U1**2
     $ - 12*QQ*S*T*U1*TM2 + 10*QQ*S*T*TM2**2 - 4*QQ*S*T1**2*TM2
     $ - 4*QQ*S*T1*U1*TM2 + 26*QQ*S*T1*TM2**2 - 2*QQ*S*U1**2*TM2
     $ + 12*QQ*S*U1*TM2**2 - 30*QQ*S*TM2**3 - QQ*T**3*TM2
     $ - 2*QQ*T**2*T1*TM2 - 4*QQ*T**2*U1*TM2 + 5*QQ*T**2*TM2**2
     $ - 4*QQ*T*T1**2*TM2 - 4*QQ*T*T1*U1*TM2 + 22*QQ*T*T1*TM2**2
     $ - 2*QQ*T*U1**2*TM2 + 18*QQ*T*U1*TM2**2 )/8.
      AJLZT3(S,T,T1,U1,QQ,TM2)=
     $ + A2B2*AL2BL2 * ( - 26*QQ*T*TM2**3 + 4*QQ*T1**2*TM2**2
     $ + 4*QQ*T1*U1*TM2**2 - 24*QQ*T1*TM2**3 + 2*QQ*U1**2*TM2**2
     $ - 14*QQ*U1*TM2**3 + 30*QQ*TM2**4 - 8*S**2*T*U1*TM2
     $ - 4*S**2*T1**2*TM2 + 8*S**2*T1*TM2**2 + 8*S**2*U1*TM2**2
     $ - 4*S**2*TM2**3 - 8*S*T**2*U1*TM2 - 8*S*T*T1**2*TM2
     $ - 8*S*T*T1*U1*TM2 + 16*S*T*T1*TM2**2 - 8*S*T*U1**2*TM2
     $ + 40*S*T*U1*TM2**2 - 8*S*T*TM2**3 + 16*S*T1**2*TM2**2
     $ + 8*S*T1*U1*TM2**2 - 32*S*T1*TM2**3 + 8*S*U1**2*TM2**2
     $ - 32*S*U1*TM2**3 + 16*S*TM2**4 - 4*T**2*T1**2*TM2
     $ - 8*T**2*T1*U1*TM2 + 8*T**2*T1*TM2**2 - 4*T**2*U1**2*TM2
     $ + 16*T**2*U1*TM2**2 - 4*T**2*TM2**3 + 16*T*T1**2*TM2**2
     $ + 24*T*T1*U1*TM2**2 - 32*T*T1*TM2**3 + 16*T*U1**2*TM2**2
     $ - 48*T*U1*TM2**3 + 16*T*TM2**4 - 16*T1**2*TM2**3
     $ - 16*T1*U1*TM2**3 + 32*T1*TM2**4 - 12*U1**2*TM2**3
     $ + 32*U1*TM2**4 - 16*TM2**5 )/8.
      AJLZT4(S,T,T1,U1,QQ,TM2)=
     $ + AB*ALBL * ( 8*QQ**3*S*T - 8*QQ**3*S*TM2 - 8*QQ**3*T*TM2
     $ + 8*QQ**3*TM2**2 - 8*QQ**2*S**2*TM2 - 8*QQ**2*S*T**2
     $ + 16*QQ**2*S*T*T1 + 8*QQ**2*S*T*U1 - 16*QQ**2*S*T*TM2
     $ - 16*QQ**2*S*T1*TM2 - 8*QQ**2*S*U1*TM2 + 40*QQ**2*S*TM2**2
     $ - 16*QQ**2*T*T1*TM2 - 8*QQ**2*T*U1*TM2 + 32*QQ**2*T*TM2**2
     $ + 16*QQ**2*T1*TM2**2 + 8*QQ**2*U1*TM2**2 - 40*QQ**2*TM2**3
     $ - 4*QQ*S**3*T + 4*QQ*S**3*TM2 - 8*QQ*S**2*T*T1
     $ - 8*QQ*S**2*T*U1 + 20*QQ*S**2*T*TM2 - 8*QQ*S**2*T1*TM2
     $ + 8*QQ*S**2*U1*TM2 - 4*QQ*S**2*TM2**2 + 4*QQ*S*T**3
     $ - 8*QQ*S*T**2*T1 - 4*QQ*S*T**2*TM2 + 40*QQ*S*T1*TM2**2
     $ - 32*QQ*S*TM2**3 - 4*QQ*T**3*TM2 - 8*QQ*T**2*T1*TM2
     $ - 16*QQ*T**2*U1*TM2 + 20*QQ*T**2*TM2**2 + 40*QQ*T*T1*TM2**2
     $ + 40*QQ*T*U1*TM2**2 - 48*QQ*T*TM2**3 - 48*QQ*T1*TM2**3
     $ - 24*QQ*U1*TM2**3 + 48*QQ*TM2**4 )/8.
C
C          Find whether JETTYP(1) or JETTYP(2) is particle
C
      JF=1
      IF(2*(JETTYP(1)/2).NE.JETTYP(1)) JF=2
C
C          Kinematics
C
      QMW2=QMW**2
      QZW=QTMW*SINH(YW)
      Q0W=QTMW*COSH(YW)
      QW=SQRT(QZW**2+QTW**2)

      T1=-X2*ECM*PT(JF)*EXP(YJ(JF))
      U1=-X1*ECM*PT(JF)*EXP(-YJ(JF))
      S1=-T1-U1-QMW2
      SIGLLQ=0.
      IF(STDDY) GO TO 400
C
C          qk + qb --> gl + w
C
      IF(JETTYP(3).EQ.1) THEN
        IFLL=JETTYP(1)/2
        COUPL=-ALFA**2*ALFQSQ*PROP(JWTYP)/(9.*PI*SCM*S)
        DENOM=S**2*EXP(.5*ALOG(QTW**4+QT2CUT(QMW)**2))
        TDIR=FCDIR(T,S,U,T1,S1,U1)*(AQ(IFLL,JWTYP)**2+BQ(IFLL,JWTYP)**2)
        UDIR=FCDIR(U,S,T,U1,S1,T1)*(AQ(IFLL,JWTYP)**2+BQ(IFLL,JWTYP)**2)
        TINT=FCINT(T,S,U,T1,S1,U1)*2.*AQ(IFLL,JWTYP)*BQ(IFLL,JWTYP)
        UINT=FCINT(U,S,T,U1,S1,T1)*2.*AQ(IFLL,JWTYP)*BQ(IFLL,JWTYP)
        IQ1=INITYP(1)
        IQ2=INITYP(2)
        IFL1=IQ1/2
        IF(2*IFL1.EQ.IQ1) THEN
          TERM1=TDIR*(AQ(IFL1,JWTYP)**2+BQ(IFL1,JWTYP)**2)
     $    *QSAVE(IQ1,1)*QSAVE(IQ2,2)
          TERM2=TINT*2.*AQ(IFL1,JWTYP)*BQ(IFL1,JWTYP)
     $    *QSAVE(IQ1,1)*QSAVE(IQ2,2)
          SIGLLQ=SIGLLQ+TERM1+TERM2
        ELSE
          TERM1=UDIR*(AQ(IFL1,JWTYP)**2+BQ(IFL1,JWTYP)**2)
     $    *QSAVE(IQ1,1)*QSAVE(IQ2,2)
          TERM2=UINT*2.*AQ(IFL1,JWTYP)*BQ(IFL1,JWTYP)
     $    *QSAVE(IQ1,1)*QSAVE(IQ2,2)
          SIGLLQ=SIGLLQ+TERM1+TERM2
        ENDIF
        SIGLLQ=SIGLLQ*UNITS
        IF(JETTYP(1).LE.13) SIGLLQ=3.*SIGLLQ
        RETURN
C
C          qk + gl --> qk + w
C
      ELSEIF(JETTYP(3).LE.NZERO(JWTYP)) THEN
        JQK=MATCH(JETTYP(3),4)
        JQK=MATCH(JQK,JWTYP)
        IF(JQK.EQ.0) RETURN
        COUPL=ALFA**2*ALFQSQ*PROP(JWTYP)/(24.*PI*SCM*S)
        DENOM=S**2*EXP(.5*ALOG(QTW**4+QT2CUT(QMW)**2))
        IFLQ=JQK/2
        IFLL=JETTYP(1)/2
        IF(INITYP(2).EQ.1) THEN
          TDIR=FCDIR(S,T,U,S1,T1,U1)*QSAVE(JQK,1)*QSAVE(1,2)
        ELSE
          TDIR=FCDIR(S,U,T,S1,U1,T1)*QSAVE(JQK,2)*QSAVE(1,1)
        ENDIF
        TDIR=TDIR*(AQ(IFLQ,JWTYP)**2+BQ(IFLQ,JWTYP)**2)
     $  *(AQ(IFLL,JWTYP)**2+BQ(IFLL,JWTYP)**2)
        IF(INITYP(2).EQ.1) THEN
          TINT=FCINT(S,T,U,S1,T1,U1)*QSAVE(JQK,1)*QSAVE(1,2)
        ELSE
          TINT=FCINT(S,U,T,S1,U1,T1)*QSAVE(JQK,2)*QSAVE(1,1)
        ENDIF
        TINT=TINT*4.*AQ(IFLQ,JWTYP)*BQ(IFLQ,JWTYP)*AQ(IFLL,JWTYP)
     $  *BQ(IFLL,JWTYP)
        SIGLLQ=TDIR+TINT
        SIGLLQ=SIGLLQ*UNITS
        IF(JETTYP(1).LE.13) SIGLLQ=3.*SIGLLQ
        SIGLLQ=SIGLLQ*QTW**2/(QTW**2+QT2CUT(QMW))
        RETURN
C
C          bt,tp + gl --> bt,tp + W,Z
C
      ELSEIF(JETTYP(3).GE.NZERO(JWTYP)+1) THEN
        IQ3=JETTYP(3)
        JQK=MATCH(IQ3,4)
        JQK=MATCH(JQK,JWTYP)
        IF(JQK.EQ.0) RETURN
        AMT=AMASS(6)
        AMT2=AMT**2
        Q2=QMW2
        AMF2=AMFAC(IQ3)*AMT2
        AMI2=AMFAC(JQK)*AMT2
        EFWT=SQRT(P(3)**2+AMF2)
        P3Z=P(3)*CTH(3)
        SWT=QMW2+AMF2+2.*Q0W*EFWT-2.*QZW*P3Z+2.*PT(3)**2
C          Kinematics
        IF(INITYP(2).EQ.1) THEN
          P1WT=EFWT+P3Z+Q0W+QZW
          P1M=AMI2/P1WT
          P2WT=EFWT-P3Z+Q0W-QZW-P1M
          X1WT=.5*P1WT/HALFE
          X2WT=.5*P2WT/HALFE
          TWT=-P1WT*(EFWT-P3Z)-P1M*(P(3)+P3Z)+AMI2+AMF2
          UWT=-P2WT*(EFWT+P3Z)+AMF2
          T1WT=-X2WT*ECM*PT(JF)*EXP(YJ(JF))
          U1WT=-X1WT*ECM*PT(JF)*EXP(-YJ(JF))-P1M*PT(JF)*EXP(YJ(JF))
          S1WT=-T1WT-U1WT-QMW2+AMI2+AMF2
        ELSE
          P2WT=EFWT-P3Z+Q0W-QZW
          P2M=AMI2/P2WT
          P1WT=EFWT+P3Z+Q0W+QZW-P2M
          X1WT=.5*P1WT/HALFE
          X2WT=.5*P2WT/HALFE
          TWT=-P1WT*(EFWT-P3Z)+AMF2
          UWT=-P2WT*(EFWT+P3Z)-P2M*(EFWT-P3Z)+AMI2+AMF2
          T1WT=-X2WT*ECM*PT(JF)*EXP(YJ(JF))-P2M*PT(JF)*EXP(-YJ(JF))
          U1WT=-X1WT*ECM*PT(JF)*EXP(-YJ(JF))
          S1WT=-T1WT-U1WT-QMW2+AMI2+AMF2
        ENDIF
C          Cross section
        SIG0=-ALFA**2*ALFQSQ/(12*PI*SCM*SWT)*PROP(JWTYP)*UNITS
        IF(JETTYP(1).LE.13) SIG0=3*SIG0
        IF(JWTYP.EQ.2.OR.JWTYP.EQ.3) THEN
          SIG0=SIG0*(AQ(6,JWTYP)**2+BQ(6,JWTYP)**2)**2
          IF(INITYP(2).EQ.1.AND.(IQ3.EQ.12.OR.IQ3.EQ.13)) THEN
            SIGLLQ=AJLWT(SWT/SWT,TWT/SWT,T1WT/SWT,U1WT/SWT,Q2/SWT,
     $      AMT2/SWT)
            SIGLLQ=SIGLLQ*SWT*(SWT/(SWT-AMI2))**2*(SWT/(TWT-AMF2))**2
            SIGLLQ=SIGLLQ*SIG0*QFCN(X1WT,JQK,1)*QFCN(X2WT,1,2)
          ELSEIF(INITYP(1).EQ.1.AND.(IQ3.EQ.12.OR.IQ3.EQ.13)) THEN
            SIGLLQ=AJLWT(SWT/SWT,UWT/SWT,U1WT/SWT,T1WT/SWT,Q2/SWT,
     $      AMT2/SWT)
            SIGLLQ=SIGLLQ*SWT*(SWT/(SWT-AMI2))**2*(SWT/(TWT-AMF2))**2
            SIGLLQ=SIGLLQ*SIG0*QFCN(X1WT,JQK,2)*QFCN(X2WT,1,1)
          ENDIF
        ELSEIF(JWTYP.EQ.4) THEN
          A2=AQ(6,JWTYP)**2
          A2B2=AQ(6,JWTYP)**2+BQ(6,JWTYP)**2
          AB=AQ(6,JWTYP)*BQ(6,JWTYP)
          AL2BL2=AQ(JETTYP(1)/2,JWTYP)**2+BQ(JETTYP(1)/2,JWTYP)**2
          ALBL=AQ(JETTYP(1)/2,JWTYP)*BQ(JETTYP(1)/2,JWTYP)
          IF(INITYP(2).EQ.1) THEN
            SIGLLQ=AJLZT1(SWT/SWT,TWT/SWT,T1WT/SWT,U1WT/SWT,
     $      Q2/SWT,AMT2/SWT)
            SIGLLQ=SIGLLQ+AJLZT2(SWT/SWT,TWT/SWT,T1WT/SWT,U1WT/SWT,
     $      Q2/SWT,AMT2/SWT)
            SIGLLQ=SIGLLQ+AJLZT3(SWT/SWT,TWT/SWT,T1WT/SWT,U1WT/SWT,
     $      Q2/SWT,AMT2/SWT)
            SIGLLQ=SIGLLQ+AJLZT4(SWT/SWT,TWT/SWT,T1WT/SWT,U1WT/SWT,
     $      Q2/SWT,AMT2/SWT)
            SIGLLQ=SIGLLQ*SWT*(SWT/(SWT-AMI2))**2*(SWT/(TWT-AMF2))**2
            SIGLLQ=SIGLLQ*SIG0*QFCN(X1WT,JQK,1)*QFCN(X2WT,1,2)
          ELSEIF(INITYP(1).EQ.1) THEN
            SIGLLQ=AJLZT1(SWT/SWT,UWT/SWT,U1WT/SWT,T1WT/SWT,
     $      Q2/SWT,AMT2/SWT)
            SIGLLQ=SIGLLQ+AJLZT2(SWT/SWT,UWT/SWT,U1WT/SWT,T1WT/SWT,
     $      Q2/SWT,AMT2/SWT)
            SIGLLQ=SIGLLQ+AJLZT3(SWT/SWT,UWT/SWT,U1WT/SWT,T1WT/SWT,
     $      Q2/SWT,AMT2/SWT)
            SIGLLQ=SIGLLQ+AJLZT4(SWT/SWT,UWT/SWT,U1WT/SWT,T1WT/SWT,
     $      Q2/SWT,AMT2/SWT)
            SIGLLQ=SIGLLQ*SWT*(SWT/(SWT-AMI2))**2*(SWT/(UWT-AMF2))**2
            SIGLLQ=SIGLLQ*SIG0*QFCN(X1WT,JQK,2)*QFCN(X2WT,1,1)
          ENDIF
        ENDIF
      ENDIF
      RETURN
C
C          Standard Drell-Yan with QT=0.
C
400   CONTINUE
      IFLL=JETTYP(1)/2
      COUPL=ALFA**2*PROP(JWTYP)*UNITS
      TDIR=COUPL*(AQ(IFLL,JWTYP)**2+BQ(IFLL,JWTYP)**2)
     $*((U1**2+T1**2)/(6.*SCM*QMW2))
      TINT=COUPL*2.*AQ(IFLL,JWTYP)*BQ(IFLL,JWTYP)
     $*((U1**2-T1**2)/(6.*SCM*QMW2))
      IQ1=INITYP(1)
      IQ2=INITYP(2)
      IFL1=IQ1/2
      TERM1=TDIR*(AQ(IFL1,JWTYP)**2+BQ(IFL1,JWTYP)**2)
     $*QSAVE(IQ1,1)*QSAVE(IQ2,2)
      TERM2=-TINT*2.*AQ(IFL1,JWTYP)*BQ(IFL1,JWTYP)
     $*QSAVE(IQ1,1)*QSAVE(IQ2,2)
      IF(2*IFL1.EQ.IQ1) SIGLLQ=SIGLLQ+TERM1+TERM2
      IF(2*IFL1.NE.IQ1) SIGLLQ=SIGLLQ+TERM1-TERM2
      IF(JETTYP(1).LE.13) SIGLLQ=3.*SIGLLQ
      RETURN
      END
