#include "isajet/pilot.h"
      SUBROUTINE TWOJET
C
C          Driving routine to generate initial parameters for jets,
C          assuming zero initial transverse momentum, ie PT(1)=PT(2).
C
C          Parameters are PT,YJ,PHI with P,YJ,XJ as dependent variables,
C          where YJ=RAPIDITY, XJ=Feynman X.
C          All parameters are stored in COMMON/JETPAR/.
C          Cross section is called from NOGOOD.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/idrun.inc"
#include "isajet/itapes.inc"
#include "isajet/keys.inc"
#include "isajet/mbpar.inc"
#include "isajet/pjets.inc"
#include "isajet/pinits.inc"
#include "isajet/jetlim.inc"
#include "isajet/ptpar.inc"
#include "isajet/jetpar.inc"
#include "isajet/primar.inc"
#include "isajet/partcl.inc"
#include "isajet/const.inc"
#include "isajet/jetsig.inc"
#include "isajet/totals.inc"
#include "isajet/isloop.inc"
#include "isajet/sstype.inc"
#include "isajet/xmssm.inc"
C
      REAL ACOSH,XXX,WTFCN,PPP,RANF,SIGN,SGN,AMQ1,AMASS,AMQ2
      REAL PPLUS,PMINUS,PSUM3,PSUM4,PPL,PMN,SQ1,SQ2,ROOT,P1PL,P1MN
      REAL P2PL,P2MN,AMI1,AMI2
      INTEGER NREJ,I,II,IS,IFL1,IFL2
      REAL X(2)
      EQUIVALENCE (X(1),X1)
      LOGICAL NOGOOD
      LOGICAL YGENJ
      INTEGER LISTJ(17),LISTW(4),LISTSS(85),LISTSM(30)
C
C          SUSY IDENT codes from /SSTYPE/. (Fortran 77 allows - signs
C          in parameter statements but not data statements.)
      INTEGER MSUPL,MSDNL,MSSTL,MSCHL,MSBT1,MSTP1,
     $MSUPR,MSDNR,MSSTR,MSCHR,MSBT2,MSTP2,MSW1,MSW2,
     $MSNEL,MSEL,MSNML,MSMUL,MSNTL,MSTAU1,MSER,MSMUR,MSTAU2
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
      PARAMETER (MSNEL=-ISNEL)
      PARAMETER (MSEL=-ISEL)
      PARAMETER (MSNML=-ISNML)
      PARAMETER (MSMUL=-ISMUL)
      PARAMETER (MSNTL=-ISNTL)
      PARAMETER (MSTAU1=-ISTAU1)
      PARAMETER (MSER=-ISER)
      PARAMETER (MSMUR=-ISMUR)
      PARAMETER (MSTAU2=-ISTAU2)
C
      DATA LISTSS/ISGL,
     $ISUPL,MSUPL,ISDNL,MSDNL,ISSTL,MSSTL,ISCHL,MSCHL,ISBT1,MSBT1,
     $ISTP1,MSTP1,
     $ISUPR,MSUPR,ISDNR,MSDNR,ISSTR,MSSTR,ISCHR,MSCHR,ISBT2,MSBT2,
     $ISTP2,MSTP2,
     $ISW1,MSW1,ISW2,MSW2,ISZ1,ISZ2,ISZ3,ISZ4,
     $ISNEL,MSNEL,ISEL,MSEL,ISNML,MSNML,ISMUL,MSMUL,ISNTL,MSNTL,
     $ISTAU1,MSTAU1,ISER,MSER,ISMUR,MSMUR,ISTAU2,MSTAU2,
     $9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,11,-11,12,-12,13,-13,
     $14,-14,15,-15,16,-16,10,80,-80,90,82,83,84,86,-86/
      DATA LISTSM/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,11,-11,12,-12,13,-13,
     $14,-14,15,-15,16,-16,10,80,-80,90,81/
      DATA LISTJ/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,7,-7,8,-8/
      DATA LISTW/10,80,-80,90/
C          Inverse hyperbolic cosine function
      ACOSH(XXX)=ALOG(XXX+SQRT(XXX**2-1.))
      WTFCN(PPP)=2.*PPP*PTGEN2*PTGEN3*PPP**((PTGEN3-1.)/PTGEN3)
C
C          Initialize
C
      NPTCL=0
      PHI(1)=PHIMIN(1)+(PHIMAX(1)-PHIMIN(1))*RANF()
      PHI(2)=AMOD(PHI(1)+PI,2.*PI)
      NREJ=-1
      SIGMA=0.
      WT=1.
      IF(.NOT.FIXPT(2)) GOTO 101
      FIXPT(1)=.TRUE.
      PT(1)=PT(2)
  101 CONTINUE
      IF(FIXPT(1)) GOTO 400
      DO 110 I=1,2
      IF(FIXP(I)) GOTO 200
      IF(FIXXJ(I)) GOTO 300
  110 CONTINUE
C
C          Genetate PT and YJ with no variables fixed
C
  111 NREJ=NREJ+1
      IF(NREJ.GT.NTRIES) GO TO 910
      SUMWT=SUMWT+SIGMA*WT/(NEVOLV*NFRGMN)
      NKINPT=NKINPT+1
      SIGMA=0.
      WT=1.
C            Generate PT with a power law distribution
      PT(1)=(PTGEN1+PTGEN2*RANF())**PTGEN3
      PT(2)=PT(1)
      SIGMAX=PTFUN1*PT(1)**PTFUN2
C          GENERATE FLAT IN YJ, CALCULATE CORRESPONDING TH
      DO 115 I=1,2
      IF(FIXYJ(I)) GOTO 115
      IF(.NOT.YGENJ(I)) GOTO 111
  115 CONTINUE
      DO 116 I=1,2
      P(I)=PT(I)/STH(I)
      IF(P(I).LT.PMIN(I).OR.P(I).GT.PMAX(I)) GOTO 111
      XJ(I)=P(I)*CTH(I)/HALFE
      IF(XJ(I).LT.XJMIN(I).OR.XJ(I).GT.XJMAX(I)) GOTO 111
  116 CONTINUE
      WT=WT*WTFCN(PT(1))
      IF(NOGOOD(1)) GOTO 111
      SUMWT=SUMWT+SIGMA*WT/(NEVOLV*NFRGMN)
      NKEEP=NKEEP+1
      GO TO 500
C
C          Generate PT and YJ fixing P
C
  200 CONTINUE
      II=3-I
  211 NREJ=NREJ+1
      IF(NREJ.GT.NTRIES) GO TO 910
      NKINPT=NKINPT+1
      WT=0.
      IF(FIXYJ(I)) GOTO 212
C          Generate PT with a power law distribution
      PT(1)=(PTGEN1+PTGEN2*RANF())**PTGEN3
      SIGMAX=PTFUN1*PT(1)**PTFUN2
      PT(2)=PT(1)
C          Given PT, TH is fixed except for a sign
      STH(I)=PT(I)/P(I)
      SIGN=1.0
      IF(RANF().GT.0.5) SIGN=-1.0
      CTH(I)=SIGN*SQRT(1.-STH(I)**2)
      TH(I)=ATAN2(STH(I),CTH(I))
      YJ(I)=-ALOG(TAN(TH(I)/2.))
      IF(YJ(I).LT.YJMIN(I).OR.YJ(I).GT.YJMAX(I)) GOTO 211
      GOTO 213
  212 PT(1)=P(I)*STH(I)
  213 CONTINUE
      XJ(I)=P(I)*CTH(I)/HALFE
      IF(XJ(I).LT.XJMIN(I).OR.XJ(I).GT.XJMAX(I)) GOTO 211
      IF(FIXP(II)) GOTO 220
      IF(FIXXJ(II)) GOTO 230
      IF(FIXYJ(II)) GOTO 215
      IF(.NOT.YGENJ(II)) GOTO 211
  215 CONTINUE
      P(II)=PT(II)/STH(II)
      IF(P(II).LT.PMIN(II).OR.P(II).GT.PMAX(II)) GOTO 211
      XJ(II)=P(II)*CTH(II)/HALFE
      IF(XJ(II).LT.XJMIN(II).OR.XJ(II).GT.XJMAX(II)) GOTO 211
      GOTO 250
220   STH(II)=PT(II)/P(II)
      SGN=1.0
      IF(RANF().GT.0.5) SGN=-1.0
      CTH(II)=SGN*SQRT(1.-STH(II)**2)
      TH(II)=ATAN2(STH(II),CTH(II))
      YJ(II)=-ALOG(TAN(TH(II)/2.))
      IF(YJ(II).LT.YJMIN(II).OR.YJ(II).GT.YJMAX(II)) GOTO 211
      XJ(II)=P(II)*CTH(II)/HALFE
      IF(XJ(II).LT.XJMIN(II).OR.XJ(II).GT.XJMAX(II)) GOTO 211
      GOTO 250
  230 TH(II)=ATAN2(PT(II),XJ(II)*HALFE)
      YJ(II)=-ALOG(TAN(TH(II)/2.))
      IF(YJ(II).LT.YJMIN(II).OR.YJ(II).GT.YJMAX(II)) GOTO 211
      CTH(II)=COS(TH(II))
      STH(II)=SIN(TH(II))
  250 CONTINUE
      IF(NOGOOD(1)) GOTO 211
      NKEEP=NKEEP+1
      GO TO 500
C
C          Generate PT and YJ at fixed XJ
C
  300 CONTINUE
      II=3-I
  311 NREJ=NREJ+1
      IF(NREJ.GT.NTRIES) GO TO 910
      NKINPT=NKINPT+1
      WT=0.
C          Generate PT with a power law distribution
      PT(1)=(PTGEN1+PTGEN2*RANF())**PTGEN3
      SIGMAX=PTFUN1*PT(1)**PTFUN2
      PT(2)=PT(1)
      TH(I)=ATAN2(PT(I),XJ(I)*HALFE)
      YJ(I)=-ALOG(TAN(TH(I)/2.))
      IF(YJ(I).LT.YJMIN(I).OR.YJ(I).GT.YJMAX(I)) GOTO 311
      CTH(I)=COS(TH(I))
      STH(I)=SIN(TH(I))
      P(I)=PT(I)/STH(I)
      IF(FIXYJ(II)) GOTO 315
      IF(FIXP(II)) GOTO 314
      YJ(II)=YJMIN(II)+(YJMAX(II)-YJMIN(II))*RANF()
      TH(II)=2.*ATAN(EXP(-YJ(II)))
      CTH(II)=COS(TH(II))
      STH(II)=SIN(TH(II))
      GOTO 315
  314 CONTINUE
      STH(II)=PT(II)/P(II)
      CTH(II)=SQRT(1.-STH(II)**2)
      IF(RANF().GT.0.5) CTH(II)=-CTH(II)
      TH(II)=ATAN2(STH(II),CTH(II))
      YJ(II)=-ALOG(TAN(TH(II)/2.))
  315 CONTINUE
      P(II)=PT(II)/STH(II)
      XJ(II)=P(II)*CTH(II)/HALFE
      IF(XJ(II).LT.XJMIN(II).OR.XJ(II).GT.XJMAX(II)) GOTO 311
      IF(NOGOOD(1)) GOTO 311
      NKEEP=NKEEP+1
      GO TO 500
C
C          Generate YJ at fixed PT
C
  400 CONTINUE
      PT(2)=PT(1)
  411 NREJ=NREJ+1
      IF(NREJ.GT.NTRIES) GO TO 910
      NKINPT=NKINPT+1
      WT=0.
      DO 415 I=1,2
      IF(FIXYJ(I)) GOTO 415
      IF(FIXP(I)) GOTO 413
      IF(.NOT.YGENJ(I)) GO TO 411
      GOTO 414
  413 CONTINUE
      IS=1
      IF(RANF().GT.0.5) IS=2
      CTH(I)=CTHS(IS,I)
      TH(I)=THS(IS,I)
      YJ(I)=YJS(IS,I)
  414 CONTINUE
      P(I)=PT(I)/STH(I)
      XJ(I)=P(I)*CTH(I)/HALFE
  415 CONTINUE
      IF(NOGOOD(1)) GOTO 411
      NKEEP=NKEEP+1
C
C          Reset /JETPAR/
C
  500 CONTINUE
      IF(KEYS(1)) THEN
        IFL1=LISTJ(JETTYP(1))
        IFL2=LISTJ(JETTYP(2))
        AMQ1=AMASS(IFL1)
        AMQ2=AMASS(IFL2)
        AMI1=AMASS(LISTJ(INITYP(1)))
        AMI2=AMASS(LISTJ(INITYP(2)))
        CALL TWOKIN(AMI1,AMI2,AMQ1,AMQ2)
      ELSEIF(KEYS(5).OR.(KEYS(10).AND.GOMSSM)) THEN
        IFL1=LISTSS(JETTYP(1))
        IFL2=LISTSS(JETTYP(2))
        AMQ1=AMASS(IFL1)
        AMQ2=AMASS(IFL2)
        CALL TWOKIN(0.,0.,AMQ1,AMQ2)
      ELSEIF(KEYS(6)) THEN
        IFL1=LISTW(JETTYP(1))
        IFL2=LISTW(JETTYP(2))
        AMQ1=AMASS(IFL1)
        AMQ2=AMASS(IFL2)
        CALL TWOKIN(0.,0.,AMQ1,AMQ2)
      ELSEIF(KEYS(8)) THEN
        IF(JETTYP(1).LE.13) THEN
          IFL1=LISTJ(JETTYP(1))
        ELSE
          IFL1=10
        ENDIF
        IF(JETTYP(2).LE.13) THEN
          IFL2=LISTJ(JETTYP(2))
        ELSE
          IFL2=10
        ENDIF
        AMQ1=AMASS(IFL1)
        AMQ2=AMASS(IFL2)
        CALL TWOKIN(0.,0.,AMQ1,AMQ2)
      ELSEIF(KEYS(10).AND.(.NOT.GOMSSM)) THEN
        IFL1=LISTSM(JETTYP(1))
        IFL2=LISTSM(JETTYP(2))
        AMQ1=AMASS(IFL1)
        AMQ2=AMASS(IFL2)
        CALL TWOKIN(0.,0.,AMQ1,AMQ2)
      ENDIF
C
C            Set PBEAM and PJETS
C
      PBEAM(1)=(1.-X1)*HALFE
      PBEAM(2)=(1.-X2)*HALFE
      DO 501 I=1,2
        PJETS(3,I)=P(I)*CTH(I)
        PJETS(1,I)=PT(I)*COS(PHI(I))
        PJETS(2,I)=PT(I)*SIN(PHI(I))
        IF(KEYS(1)) THEN
          IDJETS(I)=LISTJ(JETTYP(I))
        ELSEIF(KEYS(5).OR.(KEYS(10).AND.GOMSSM)) THEN
          IDJETS(I)=LISTSS(JETTYP(I))
        ELSEIF(KEYS(6)) THEN
          IDJETS(I)=LISTW(JETTYP(I))
        ELSEIF(KEYS(8)) THEN
          IDJETS(1)=IFL1
          IDJETS(2)=IFL2
        ELSEIF(KEYS(10)) THEN
          IDJETS(I)=LISTSM(JETTYP(I))
        ENDIF
        PJETS(5,I)=AMASS(IDJETS(I))
        PJETS(4,I)=SQRT(P(I)**2+PJETS(5,I)**2)
  501 CONTINUE
C
C          Set PINITS
C
      DO 600 I=1,2
      IDINIT(I)=LISTJ(INITYP(I))
      PINITS(5,I)=AMASS(IDINIT(I))
      PPLUS=X(I)*ECM
      PMINUS=PINITS(5,I)**2/PPLUS
      PINITS(4,I)=.5*(PPLUS+PMINUS)
      PINITS(3,I)=.5*(PPLUS-PMINUS)*(3-2*I)
      PINITS(2,I)=0.
      PINITS(1,I)=0.
600   CONTINUE
C          Calculate PINITS exactly.
      PSUM3=PJETS(3,1)+PJETS(3,2)
      PSUM4=PJETS(4,1)+PJETS(4,2)
      IF(PSUM3.GT.0.) THEN
        PPL=PSUM4+PSUM3
        PMN=SHAT/PPL
      ELSE
        PMN=PSUM4-PSUM3
        PPL=SHAT/PMN
      ENDIF
      SQ1=PINITS(5,1)**2
      SQ2=PINITS(5,2)**2
      ROOT=SQRT((PPL*PMN-SQ1-SQ2)**2-4.*SQ1*SQ2)
      P1PL=(PPL*PMN+SQ1-SQ2+ROOT)/(2.*PMN)
      P1MN=SQ1/P1PL
      P2MN=(PPL*PMN+SQ2-SQ1+ROOT)/(2.*PPL)
      P2PL=SQ2/P2MN
      PINITS(4,1)=.5*(P1PL+P1MN)
      PINITS(3,1)=.5*(P1PL-P1MN)
      PINITS(4,2)=.5*(P2PL+P2MN)
      PINITS(3,2)=.5*(P2PL-P2MN)
      RETURN
C
C          Error
C
910   CALL PRTEVT(0)
      WRITE(ITLIS,1000) NREJ
 1000 FORMAT(//' IT IS TAKING MORE THAN',I5,' TRIES TO GENERATE AN',
     $' EVENT. CHECK LIMITS OR INCREASE NTRIES.')
      STOP 99
      END
