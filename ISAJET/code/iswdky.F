#include "isajet/pilot.h"
      SUBROUTINE ISWDKY
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-       decay W's and Z's as done in ISAJET
C-
C-   Created   6-MAY-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/const.inc"
#include "isajet/frame.inc"
#include "isajet/jetpar.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/pjets.inc"
#include "isajet/partcl.inc"
#include "isajet/primar.inc"
#include "isajet/wcon.inc"
      REAL X(2)    
      EQUIVALENCE (X(1),X1) 
      REAL PREST(5),PL(5),EL(3),EML(3),EMSQL(3)    
      REAL WTFAC(3)    
      REAL BRANCH(29)
      INTEGER LISTJ(29),LISTW(4)   
      REAL RANF,SUM,PTDEN,QDEN,ETA,
     $S12,SUMBR,BRMODE,AMASS,BRINV,TRY,PL12,
     $COSTHL,THL,PHL,PTL,SGN,BP,PLPL,PLMN,AMINI,AMFIN,PINI,PFIN, 
     $ QPL,QMN,AM1SQ,AM2SQ,ROOT,P1PL,P1MN,P2PL,P2MN
      INTEGER NADD,K,IQ1,IQ2,IFL1,IFL2,IQ,IFL,I   
      REAL EY
      REAL QWPL,QWMN
C   
      DATA LISTJ/   
     $9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,  
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,    
     $10,80,-80,90/ 
      DATA LISTW/10,80,-80,90/  
C----------------------------------------------------------------------
C   
C          Entry    
C   
      NPTCL=0   
C   
C          Kinematics. Note that YW is the true rapidity and QW is
C          the true 3-momentum. See DRLLYN.
C   
      QMW=QWJET(5)
      QTW=SQRT(QWJET(1)**2+QWJET(2)**2)
      QW=SQRT(QWJET(1)**2+QWJET(2)**2+QWJET(3)**2)
      IF(QTW.NE.0) THEN
        PHIW=ATAN2(QWJET(2),QWJET(1))
        IF(PHIW.LT.0) PHIW=PHIW+2*PI
      ELSE
        PHIW=0
      ENDIF
      QWPL=QWJET(4)+QWJET(3)
      QWMN=QWJET(4)-QWJET(3)
      IF(QWPL.GT.0..AND.QWMN.GT.0.) THEN
        YW=0.5*ALOG(QWPL/QWMN)
      ELSE
        YW=999.*SIGN(1.,QWJET(3))
      ENDIF
      IF(QW.NE.0.) THEN
        THW=ACOS(QWJET(3)/QW)
      ELSE
        THW=0.
      ENDIF
C   
C          Select W decay mode  
C          QMW dependence neglected in branching ratios 
C          BRANCH is cum. br. with heavy modes subtracted.  
C   
      S12=QMW**2    
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
105   CONTINUE    
      BRINV=1./SUMBR  
C   
      TRY=RANF()  
      DO 110 IQ=1,25  
        IF(TRY.LT.BRANCH(IQ)*BRINV.AND.MATCH(IQ,JWTYP).NE.0) THEN 
          JETTYP(1)=IQ    
          JETTYP(2)=MATCH(IQ,JWTYP)   
          GO TO 120   
        ENDIF 
110   CONTINUE    
C   
120   IFL1=LISTJ(JETTYP(1)) 
      IFL2=LISTJ(JETTYP(2)) 
C   
C          Select masses of decay products. 
C   
      EML(1)=AMASS(IFL1)    
      EML(2)=AMASS(IFL2)    
C   
C          Generate W decay in its rest frame 
C          First set up momenta of decay products:  
C   
      EMSQL(1)=EML(1)**2    
      EMSQL(2)=EML(2)**2    
      EL(1)=(S12+EMSQL(1)-EMSQL(2))/(2.*QMW)    
      EL(2)=(S12+EMSQL(2)-EMSQL(1))/(2.*QMW)    
      PL12=SQRT((S12-(EML(1)+EML(2))**2)*(S12-(EML(1)-EML(2))**2))  
     $/(2.*QMW) 
C          W momentum   
      DO 140 K=1,5
140   PREST(K)=QWJET(K)
C          Generate next W decay    
20    CONTINUE  
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
C          Set PJETS    
C   
      DO 501 I=1,2
        PJETS(3,I)=P(I)*CTH(I)  
        PJETS(1,I)=PT(I)*COS(PHI(I))    
        PJETS(2,I)=PT(I)*SIN(PHI(I))    
        PJETS(4,I)=SQRT(P(I)**2+EMSQL(I))   
        PJETS(5,I)=SQRT(EMSQL(I))   
        IDJETS(I)=LISTJ(JETTYP(I))  
501   CONTINUE  
  999 RETURN
      END
