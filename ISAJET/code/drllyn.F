#include "isajet/pilot.h"
      SUBROUTINE DRLLYN
C
C          Generate QMW (and QTW) for DRELLYAN or HIGGS event using
C          integrated cross section. Then generate decay -- for HIGGS,
C          the mode must be chosen using the integrated cross sections
C          because of interference with W+W->W+W scattering.
C
C          Note that NOGOOD calls the cross section.
C
C          Ver. 6.40: Add technicolor resonances. Use logs for QDEN,
C          PTDEN, WTFAC, etc. Also scale QMW generation by QMAX.
C
C          Ver. 7.01: Correct QDEN to correspond to correct fit form:
C          SIGMA = ANOMR(K)*(QMAX**2/Q**2)**QPOW
C          See QFUNC.
C
C          Ver. 7.14: Add SUSY Higgs
C          Ver. 7.15: Fix bug with THETAW limits by adding epsilon to
C          allowed range. Check for possible invalid Higgs decays.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/jetsig.inc"
#include "isajet/totals.inc"
#include "isajet/q1q2.inc"
#include "isajet/partcl.inc"
#include "isajet/pjets.inc"
#include "isajet/pinits.inc"
#include "isajet/wcon.inc"
#include "isajet/primar.inc"
#include "isajet/dylim.inc"
#include "isajet/const.inc"
#include "isajet/jetpar.inc"
#include "isajet/jetlim.inc"
#include "isajet/wgen.inc"
#include "isajet/dypar.inc"
#include "isajet/keys.inc"
#include "isajet/hcon.inc"
#include "isajet/isloop.inc"
#include "isajet/idrun.inc"
#include "isajet/xmssm.inc"
#include "isajet/listss.inc"
C
      DIMENSION X(2)
      EQUIVALENCE (X(1),X1)
      DIMENSION PREST(5),PL(5),EL(3),EML(3),EMSQL(3)
      DIMENSION WTFAC(3)
      LOGICAL NOGOOD
      LOGICAL YGENJ
      DIMENSION BRANCH(29),LISTJ(29),LISTW(5)
      REAL ACOSH,XXX,ASINH,CHOOSE,RANF,SUM,WTFAC,PTDEN,QDEN,ETA,QPW,
     $S12,BRANCH,SUMBR,BRMODE,AMASS,BRINV,TRY,EMSQL,EL,PL12,PREST,
     $COSTHL,THL,PHL,PTL,SGN,PL,BP,PLPL,PLMN,AMINI,AMFIN,PINI,PFIN,
     $ QPL,QMN,AM1SQ,AM2SQ,ROOT,P1PL,P1MN,P2PL,P2MN,X,EML
      INTEGER NTRY,K,IQ1,IQ2,IFL1,IFL2,LISTJ,IQ,NTRY2,IFL,LISTW,I
      REAL ZZSTAR
      INTEGER IZSTAR,JVIR,N0J
C
      DATA LISTJ/
     $9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,
     $10,80,-80,90/
      DATA LISTW/10,80,-80,90,92/
      ACOSH(XXX)=ALOG(XXX+SQRT(XXX**2-1.))
      ASINH(XXX)=ALOG(XXX+SQRT(XXX**2+1.))
C
C          Entry
C
      NPTCL=0
      NTRY=0
200   CONTINUE
      SIGMA=0.
      WT=1.
    1 CONTINUE
      NTRY=NTRY+1
      IF(NTRY.GT.NTRIES) GO TO 999
      SUMWT=SUMWT+SIGMA*WT/(NEVOLV*NFRGMN)
      NKINPT=NKINPT+1
      SIGMA=0.
      WT=1.
      DO 2 K=1,3
    2 SIGSL(K)=0
C            Choose interval for cross section calculation
      CHOOSE=RANF()
      SUM=0.
      DO 3 K=NKL,NKH
        SUM=SUM+QSELWT(K)
        IF(CHOOSE.LE.SUM) GO TO 30
3     CONTINUE
30    KSEL=K
C          Generate QTW in selected region
      IF(.NOT.FIXQT) THEN
        ETA=(PTGN(1,K)+PTGN(2,K)*RANF())**PTGN(3,K)
        PTSEL(K)=SQRT(ETA-RNU2(K))
        PTDEN=ALOG(ETA)*PTPOW(K)
        WTFAC(1)=ALOG(ABS(PTGN(2,K)))+ALOG(ABS(PTGN(3,K)))
     1  +ALOG(ABS(PTSEL(K)**2+RNU2(K)))*((PTGN(3,K)-1.)/PTGN(3,K))
        PT(3)=PTSEL(K)
      ELSE
        PTDEN=0.
        WTFAC(1)=-1000.
      ENDIF
C          Generate QMW
      IF(.NOT.FIXQM) THEN
        IF(.NOT.K.EQ.2) THEN
          QSEL(K)=QMAX**2*(QGEN(1,K)+QGEN(2,K)*RANF())**QGEN(3,K)
          QDEN=ALOG(QSEL(K)/QMAX**2)*QPOW(K)
          WTFAC(2)=ALOG(ABS(QGEN(2,K)))+ALOG(ABS(QGEN(3,K)))
     1    +ALOG(QSEL(K)/QMAX**2)*((QGEN(3,K)-1.)/QGEN(3,K))
     2    +ALOG(QMAX**2)
          QSEL(K)=SQRT(QSEL(K))
          QMW=QSEL(K)
        ELSE
          ETA=QGEN(3,K)*TAN(QGEN(1,K)+QGEN(2,K)*RANF())
          QSEL(K)=SQRT(ETA+EMSQ)
          WTFAC(2)=ALOG(QGEN(2,K))+ALOG(QGEN(3,K))
     1    +ALOG((ETA/QGEN(3,K))**2+1.)
          QMW=QSEL(K)
          QDEN=ALOG((QMW**2-EMSQ)**2+EMGAM**2)
        ENDIF
      ELSE
        QDEN=0.
        WTFAC(2)=-1000.
      ENDIF
      SIGSL(K)=EXP(ANORM(K)-PTDEN-QDEN)
C
      IF(STDDY) THEN
        WT=EXP(WTFAC(2)-ALOG(QSELWT(K)))
      ELSE
        WT=EXP(WTFAC(1)+WTFAC(2)-ALOG(QSELWT(K)))
      ENDIF      
      QTW=PT(3)
      YW=YWMIN+(YWMAX-YWMIN)*RANF()
      WT=WT*(YWMAX-YWMIN)
      PHIW=PHWMIN+(PHWMAX-PHWMIN)*RANF()
      PHI(3)=AMOD(PHIW+PI,2.*PI)
      QPW=SQRT(QTW**2+QMW**2)*SINH(YW)
      QW=SQRT(QTW**2+QPW**2)
      THW=QPW/QW
      IF(ABS(THW).GT.1.) THW=SIGN(1.,THW)
      THW=ACOS(THW)
      IF(THW.LT.THWMIN-1.E-6.OR.THW.GT.THWMAX+1.E-6) GOTO 1
      XW=QPW/HALFE
      IF(XW.LT.XWMIN.OR.XW.GT.XWMAX) GOTO 1
      IF(.NOT.STDDY) THEN
        IF(.NOT.YGENJ(3)) GOTO 1
        P(3)=PT(3)/STH(3)
        XJ(3)=P(3)*CTH(3)/HALFE
        IF(XJ(3).LT.XJMIN(3).OR.XJ(3).GT.XJMAX(3)) GOTO 1
      ENDIF
C
C          Check integrated cross section
C
      IF(NOGOOD(2)) GO TO 1
      SUMWT=SUMWT+SIGMA*WT/(NEVOLV*NFRGMN)
      NWGEN=NWGEN+1
      S12=QMW**2
C
C          No decay for KKG:
C          For compatibility reasons, the jet is still the 3rd one.
C          Jets 1 and 2 (W decay products) are voided; no decay step.
C
      IF(KEYS(11)) THEN
        DO 50 I=1,2
	  P(I)=0.
	  PT(I)=0.
	  CTH(I)=0.
	  PHI(I)=0.
	  EMSQL(I)=0.
	  IDJETS(I)=0
50      CONTINUE
	GOTO 350
      ENDIF
C
C          Select W decay mode
C          QMW dependence neglected in branching ratios
C          BRANCH is cum. br. with heavy modes subtracted.
C
      IF(KEYS(3)) THEN
        BRANCH(1)=0.
        SUMBR=0.
        DO 105 IQ1=2,25
          IQ2=MATCH(IQ1,JWTYP)
          IF(IQ2.EQ.0) THEN
            BRMODE=0.
          ELSE
            BRMODE=WCBR(IQ1,JWTYP)-WCBR(IQ1-1,JWTYP)
            IFL1=LISTJ(IQ1)
            IFL2=LISTJ(IQ2)
            IF(S12.LE.(AMASS(IFL1)+AMASS(IFL2))**2) BRMODE=0.
          ENDIF
          BRANCH(IQ1)=BRANCH(IQ1-1)+BRMODE
          SUMBR=SUMBR+BRMODE
105     CONTINUE
        BRINV=1./SUMBR
C
        TRY=RANF()
        DO 110 IQ=1,25
          IF(TRY.LT.BRANCH(IQ)*BRINV.AND.MATCH(IQ,JWTYP).NE.0) THEN
            JETTYP(1)=IQ
            JETTYP(2)=MATCH(IQ,JWTYP)
            GO TO 120
          ENDIF
110     CONTINUE
      ENDIF
C
120   IF(GOMSSM) THEN
        IFL1=LISTSS(JETTYP(1))
        IFL2=LISTSS(JETTYP(2))
      ELSE
        IFL1=LISTJ(JETTYP(1))
        IFL2=LISTJ(JETTYP(2))
      ENDIF
C
C          Select masses of decay products. These are just normal masses
C          except for Z Z* decay of Higgs, where one is virtual.
C
      EML(1)=AMASS(IFL1)
      EML(2)=AMASS(IFL2)
      IF(KEYS(7).AND.EML(1)+EML(2).GT.QMW) THEN
C          WW* or ZZ* decay - generate/check W* or Z* mass
        IF((IABS(IFL1).EQ.80.AND.IABS(IFL2).EQ.80)
     $  .OR.(IFL1.EQ.90.AND.IFL2.EQ.90)) THEN
          IZSTAR=3-2*RANF()
          IF(GOMSSM) THEN
            JVIR=JETTYP(IZSTAR)-76
          ELSE
            JVIR=JETTYP(IZSTAR)-25
          ENDIF
          EML(IZSTAR)=ZZSTAR(QMW,JVIR)
          IF(EML(IZSTAR).LT.ZSTARS(JVIR,IZSTAR)) GO TO 200
C          Other decay - invalid for this QMW
        ELSE
          GO TO 200
        ENDIF
      ENDIF
C
C          Generate W decay in its rest frame and compare with SIGDY2.
C          First set up momenta of decay products:
C
      EMSQL(1)=EML(1)**2
      EMSQL(2)=EML(2)**2
      EL(1)=(S12+EMSQL(1)-EMSQL(2))/(2.*QMW)
      EL(2)=(S12+EMSQL(2)-EMSQL(1))/(2.*QMW)
      PL12=SQRT((S12-(EML(1)+EML(2))**2)*(S12-(EML(1)-EML(2))**2))
     $/(2.*QMW)
C          W momentum
      PREST(1)=QTW*COS(PHIW)
      PREST(2)=QTW*SIN(PHIW)
      PREST(3)=QPW
      PREST(4)=SQRT(QW**2+QMW**2)
      PREST(5)=QMW
      NTRY2=0
C          Generate next W decay
20    CONTINUE
      NTRY2=NTRY2+1
      IF(NTRY2.GT.NTRIES) GO TO 999
      COSTHL=2.*RANF()-1.
      THL=ACOS(COSTHL)
      PHL=2.*PI*RANF()
      PTL=PL12*SIN(THL)
C
      DO 300 I=1,2
        SGN=3-2*I
        PL(1)=SGN*PTL*COS(PHL)
        PL(2)=SGN*PTL*SIN(PHL)
        PL(3)=SGN*PL12*COSTHL
        PL(4)=EL(I)
        PL(5)=EML(I)
C          Boost with W momentum
        BP=0.
        DO 310 K=1,3
310     BP=BP+PL(K)*PREST(K)
        BP=BP/PREST(5)
        DO 320 K=1,3
320     PL(K)=PL(K)+PREST(K)*PL(4)/PREST(5)
     $  +PREST(K)*BP/(PREST(4)+PREST(5))
        PL(4)=PL(4)*PREST(4)/PREST(5)+BP
C          Fill common blocks
        PT(I)=SQRT(PL(1)**2+PL(2)**2)
        P(I)=SQRT(PT(I)**2+PL(3)**2)
        IF(PT(I).GT.0.) THEN
          PHI(I)=ATAN2(PL(2),PL(1))
        ELSE
          PHI(I)=(I-1)*PI
        ENDIF
        IF(PHI(I).LT.0.) PHI(I)=PHI(I)+2.*PI
        CTH(I)=PL(3)/P(I)
        STH(I)=PT(I)/P(I)
        TH(I)=ACOS(CTH(I))
        XJ(I)=PL(3)/HALFE
        IF(CTH(I).GT.0.) THEN
          PLPL=PL(4)+PL(3)
          PLMN=(PT(I)**2+EMSQL(I))/PLPL
        ELSE
          PLMN=PL(4)-PL(3)
          PLPL=(PT(I)**2+EMSQL(I))/PLMN
        ENDIF
        YJ(I)=.5*ALOG(PLPL/PLMN)
300   CONTINUE
C
C          Test cross section
C          Extra kinematics for W+W->W+W
C
      IF(KEYS(7).OR.KEYS(9)) THEN
        SHAT=S12
        IF(GOMSSM) THEN
          AMINI=AMASS(LISTSS(INITYP(1)))
        ELSE
          AMINI=AMASS(LISTJ(INITYP(1)))
        ENDIF
        AMFIN=EML(1)
        PINI=.5*SQRT(S12-4.*AMINI**2)
        PFIN=PL12
        THAT=AMINI**2+AMFIN**2-.5*S12+2.*PINI*PFIN*COSTHL
        UHAT=AMINI**2+AMFIN**2-.5*S12-2.*PINI*PFIN*COSTHL
      ENDIF
C
C          Check W decay
C
      IF(NOGOOD(3)) GO TO 20
C
C          Check W decay with kinematic limits
C
      IF(NOGOOD(4)) GO TO 200
350   NKEEP=NKEEP+1
C
C            Set PBEAM
C
      PBEAM(1)=(1.-X1)*HALFE
      PBEAM(2)=(1.-X2)*HALFE
      IF(NJET.LT.3) GO TO 502
      IFL=LISTJ(JETTYP(3))
      EMSQL(3)=AMASS(IFL)**2
502   CONTINUE
C
C          Set PJETS
C
      IF(KEYS(11)) THEN
        N0J=3
      ELSE
        N0J=1
      ENDIF
      DO 501 I=N0J,NJET
        PJETS(3,I)=P(I)*CTH(I)
        PJETS(1,I)=PT(I)*COS(PHI(I))
        PJETS(2,I)=PT(I)*SIN(PHI(I))
        PJETS(4,I)=SQRT(P(I)**2+EMSQL(I))
        PJETS(5,I)=SQRT(EMSQL(I))
        IF(KEYS(7).AND.GOMSSM) THEN
          IDJETS(I)=LISTSS(JETTYP(I))
        ELSE
          IDJETS(I)=LISTJ(JETTYP(I))
        ENDIF
501   CONTINUE
C          No technicolor IDENT's defined, so...
      IF(KEYS(3)) THEN
        IDENTW=LISTW(JWTYP)
      ELSEIF(KEYS(7).AND..NOT.GOMSSM) THEN
        IDENTW=81
      ELSEIF(KEYS(7).AND.GOMSSM) THEN
        IDENTW=IHTYPE
      ELSEIF(KEYS(11)) THEN
        IDENTW=92
      ELSE
        IDENTW=0
      ENDIF
C          W momentum in /PJETS/
      IF(KEYS(11)) THEN
        QWJET(1)=QTW*COS(PHIW)
        QWJET(2)=QTW*SIN(PHIW)
        QWJET(3)=QPW
        QWJET(4)=SQRT(QW**2+QMW**2)
        QWJET(5)=QMW
      ELSE
        DO 503 K=1,4
503     QWJET(K)=PJETS(K,1)+PJETS(K,2)
        QWJET(5)=QMW
      ENDIF
C
C          Set PINITS
      DO 504 I=1,2
        IF(KEYS(7).AND.GOMSSM) THEN
          IDINIT(I)=LISTSS(INITYP(I))
        ELSE
          IDINIT(I)=LISTJ(INITYP(I))
        ENDIF
        PINITS(5,I)=AMASS(IDINIT(I))
        PINITS(1,I)=0.
        PINITS(2,I)=0.
504   CONTINUE
C          Calculate total momentum
      QPL=QWJET(4)+QWJET(3)
      QMN=QWJET(4)-QWJET(3)
      IF(NJET.EQ.3) THEN
        QPL=QPL+PJETS(4,3)+PJETS(3,3)
        QMN=QMN+PJETS(4,3)-PJETS(3,3)
      ENDIF
C          and solve initial kinematics
      AM1SQ=PINITS(5,1)**2
      AM2SQ=PINITS(5,2)**2
      ROOT=SQRT((QPL*QMN-AM1SQ-AM2SQ)**2-4.*AM1SQ*AM2SQ)
      P1PL=(QPL*QMN+AM1SQ-AM2SQ+ROOT)/(2.*QMN)
      P1MN=AM1SQ/P1PL
      P2MN=(QPL*QMN+AM2SQ-AM1SQ+ROOT)/(2.*QPL)
      P2PL=AM2SQ/P2MN
      PINITS(3,1)=.5*(P1PL-P1MN)
      PINITS(4,1)=.5*(P1PL+P1MN)
      PINITS(3,2)=.5*(P2PL-P2MN)
      PINITS(4,2)=.5*(P2PL+P2MN)
      RETURN
C
999   CALL PRTEVT(0)
      WRITE(ITLIS,9999) NTRIES
9999  FORMAT(//' IT IS TAKING MORE THAN',I5,' TRIES TO GENERATE AN',
     C' EVENT. CHECK LIMITS OR INCREASE NTRIES')
      STOP 99
      END
