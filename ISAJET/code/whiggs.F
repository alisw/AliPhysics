#include "isajet/pilot.h"
      SUBROUTINE WHIGGS
C
C          Finish generation of whiggs events started bY TWOJET.
C          Select W decay modes as allowed by WMODE1, WMODE2.
C          Generate W decay angles and put vectors in PPAIR.
C
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
#include "isajet/pjets.inc"
#include "isajet/pinits.inc"
#include "isajet/keys.inc"
#include "isajet/hcon.inc"
#include "isajet/wwpar.inc"
#include "isajet/xmssm.inc"
C
      DIMENSION X(2),LIST(25),P1WCM(4),P2WCM(4),P1LAB(4),P2LAB(4)
     1,PBOOST(4)
      EQUIVALENCE (X(1),X1)
      DIMENSION JWWTYP(2)
      REAL GVQ(2),GAQ(2),GVL(2),GAL(2)
      REAL X,RND,RANF,CBRWW,AMASS,AM0,AM1,AM2,
     $E1CM,E2CM,P12CM,CTHCM,STHCM,PHICM,CPHICM,SPHICM,P1WCM,P2WCM,
     $PBOOST,P1LAB,P2LAB,ZHSIG,ZHMAX
      INTEGER JWWTYP,JET,JWT,JQ,IQ1,IQ2,LIST,NREJ,NJ0,K
      REAL BRANCH(12),SUMBR,BETAWH,GAMWH,PZWHCM,CTHD,WHSIG
      INTEGER IDABS,IDABS1,IDABSJ,IDIABS
C
      DATA LIST/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16/
C
      GVQ(1)=.25-2*SIN2W/3.
      GVQ(2)=-.25+SIN2W/3.
      GAQ(1)=-.25
      GAQ(2)=.25
      GVL(1)=.25
      GVL(2)=-.25+SIN2W
      GAL(1)=-.25
      GAL(2)=.25
      NPAIR=0
      IF(KEYS(10).AND..NOT.GOMSSM) THEN
        JWWTYP(1)=JETTYP(1)-25
        JWWTYP(2)=JETTYP(2)-25
      ELSEIF(KEYS(10).AND.GOMSSM) THEN
        JWWTYP(1)=JETTYP(1)-76
        JWWTYP(2)=JETTYP(2)-76
      ENDIF
C
C          Select W decay modes and put in /JETSET/. First particle
C          is always the fermion.
 
      DO 200 JET=1,2
        IDABS=IABS(IDJETS(JET))
        IF(IDABS.NE.80.AND.IDABS.NE.90) GO TO 200
        RND=RANF()
        JWT=JWWTYP(JET)
C          Must only consider allowed decays for this mass
        SUMBR=0.
        DO 201 JQ=1,12
          IQ1=2*JQ
          IQ2=MATCH(IQ1,JWT)
          IF(IQ2.EQ.0) THEN
            BRANCH(JQ)=0.
            GO TO 201
          ENDIF
          AM1=AMASS(LIST(IQ1))
          AM2=AMASS(LIST(IQ2))
          IF(AM1+AM2.LT.PJETS(5,JET)) THEN
            BRANCH(JQ)=RBRWW(JQ,JWT,JET)
            SUMBR=SUMBR+BRANCH(JQ)
          ELSE
            BRANCH(JQ)=0.
          ENDIF
201     CONTINUE
        IF(SUMBR.LE.0.) GO TO 998
        DO 202 JQ=1,12
202     BRANCH(JQ)=BRANCH(JQ)/SUMBR
C
        CBRWW=0.
        DO 210 JQ=1,12
          CBRWW=CBRWW+BRANCH(JQ)
          IF(RND.GT.CBRWW) GO TO 210
          IQ1=2*JQ
          IQ2=MATCH(IQ1,JWT)
          IDPAIR(NPAIR+1)=LIST(IQ1)
          IDPAIR(NPAIR+2)=LIST(IQ2)
          PPAIR(5,NPAIR+1)=AMASS(LIST(IQ1))
          PPAIR(5,NPAIR+2)=AMASS(LIST(IQ2))
          JPAIR(NPAIR+1)=JET
          JPAIR(NPAIR+2)=JET
          NPAIR=NPAIR+2
          JQWW(JET)=JQ
          GO TO 200
210     CONTINUE
200   CONTINUE
C
C          Generate decay uniformly in angle and put in PPAIR.
C          Will check cross section later.
C
      NREJ=0
300   NJ0=2
      DO 310 JET=1,2
        IDABS=IABS(IDJETS(JET))
        IF(IDABS.NE.80.AND.IDABS.NE.90) GO TO 310
C          Construct W com momenta.
        IDABSJ=IDABS
        AM0=PJETS(5,JET)
        AM1=PPAIR(5,NJ0-1)
        AM2=PPAIR(5,NJ0)
        E1CM=(AM0**2+AM1**2-AM2**2)/(2.*AM0)
        E2CM=(AM0**2+AM2**2-AM1**2)/(2.*AM0)
        P12CM=(AM0**2-AM1**2-AM2**2)**2-4.*(AM1*AM2)**2
        P12CM=SQRT(P12CM)/(2.*AM0)
        CTHCM=2.*RANF()-1.
        STHCM=SQRT(1.-CTHCM**2)
        PHICM=2.*PI*RANF()
        CPHICM=COS(PHICM)
        SPHICM=SIN(PHICM)
        P1WCM(1)=P12CM*STHCM*CPHICM
        P2WCM(1)=-P1WCM(1)
        P1WCM(2)=P12CM*STHCM*SPHICM
        P2WCM(2)=-P1WCM(2)
        P1WCM(3)=P12CM*CTHCM
        P2WCM(3)=-P1WCM(3)
        P1WCM(4)=E1CM
        P2WCM(4)=E2CM
C          Boost to lab frame.
        DO 320 K=1,3
320     PBOOST(K)=-PJETS(K,JET)
        PBOOST(4)=PJETS(4,JET)
        CALL LBOOST(PBOOST,1,P1WCM,P1LAB)
        CALL LBOOST(PBOOST,1,P2WCM,P2LAB)
        DO 330 K=1,4
          PPAIR(K,NJ0-1)=P1LAB(K)
          PPAIR(K,NJ0)=P2LAB(K)
330     CONTINUE
      NJ0=NJ0+2
310   CONTINUE
C
C          Impose simple (1+-cos(theta))**2 decay distribution for WH
C          Must use P1 in WH CoM frame
      IF (IDABSJ.NE.80.AND.IDABSJ.NE.90) GO TO 400
      BETAWH=(PJETS(3,1)+PJETS(3,2))/(PJETS(4,1)+PJETS(4,2))
      GAMWH=1./SQRT(1.-BETAWH**2)
      PZWHCM=GAMWH*(P1LAB(3)-BETAWH*P1LAB(4))
      CTHD=PZWHCM/SQRT(P1LAB(1)**2+P1LAB(2)**2+PZWHCM**2)
      IF (IDINIT(1).LT.0) CTHD=-CTHD      
      IDIABS=IABS(IDINIT(1))
      IDABS1=IABS(IDPAIR(1))
      IF (IDABSJ.EQ.80) THEN
        WHSIG=(1.+CTHD)**2
        IF(WHSIG.GT.4*RANF()) GO TO 400
      END IF
      IF (IDABSJ.EQ.90) THEN
        IF (IDIABS.EQ.1.OR.IDIABS.EQ.4) THEN
          IF (IDABS1.EQ.1.OR.IDABS1.EQ.4) THEN
            ZHSIG=(GVQ(1)**2+GAQ(1)**2)**2*(1.+CTHD**2)
     $             +8*GVQ(1)*GAQ(1)*GVQ(1)*GAQ(1)*CTHD
            ZHMAX=2*(GVQ(1)**2+GAQ(1)**2)**2
     $             +8*GVQ(1)*GAQ(1)*GVQ(1)*GAQ(1)
          ELSEIF (IDABS1.EQ.2.OR.IDABS1.EQ.3.OR.IDABS1.EQ.5) THEN
            ZHSIG=(GVQ(1)**2+GAQ(1))*(GVQ(2)**2+GAQ(2)**2)*(1.+CTHD**2)
     $             +8*GVQ(1)*GAQ(1)*GVQ(2)*GAQ(2)*CTHD
            ZHMAX=(GVQ(1)**2+GAQ(1))*(GVQ(2)**2+GAQ(2)**2)*2
     $             +8*GVQ(1)*GAQ(1)*GVQ(2)*GAQ(2)
          ELSEIF (IDABS1.EQ.11.OR.IDABS1.EQ.13.OR.IDABS1.EQ.15) THEN
            ZHSIG=(GVQ(1)**2+GAQ(1))*(GVL(1)**2+GAL(1)**2)*(1.+CTHD**2)
     $             +8*GVQ(1)*GAQ(1)*GVL(1)*GAL(1)*CTHD
            ZHMAX=(GVQ(1)**2+GAQ(1))*(GVL(1)**2+GAL(1)**2)*2
     $             +8*GVQ(1)*GAQ(1)*GVL(1)*GAL(1)
          ELSEIF (IDABS1.EQ.12.OR.IDABS1.EQ.14.OR.IDABS1.EQ.16) THEN
            ZHSIG=(GVQ(1)**2+GAQ(1))*(GVL(2)**2+GAL(2)**2)*(1.+CTHD**2)
     $             +8*GVQ(1)*GAQ(1)*GVL(2)*GAL(2)*CTHD
            ZHMAX=(GVQ(1)**2+GAQ(1))*(GVL(2)**2+GAL(2)**2)*2
     $             +8*GVQ(1)*GAQ(1)*GVL(2)*GAL(2)
          END IF
        ELSE IF (IDIABS.EQ.2.OR.IDIABS.EQ.3.OR.IDIABS.EQ.5) THEN
          IF (IDABS1.EQ.1.OR.IDABS1.EQ.4) THEN
            ZHSIG=(GVQ(2)**2+GAQ(2)**2)**2*(1.+CTHD**2)
     $             +8*GVQ(2)*GAQ(2)*GVQ(1)*GAQ(1)*CTHD
            ZHMAX=(GVQ(2)**2+GAQ(2)**2)**2*2
     $             +8*GVQ(2)*GAQ(2)*GVQ(1)*GAQ(1)
          ELSEIF (IDABS1.EQ.2.OR.IDABS1.EQ.3.OR.IDABS1.EQ.5) THEN
            ZHSIG=(GVQ(2)**2+GAQ(2))*(GVQ(2)**2+GAQ(2)**2)*(1.+CTHD**2)
     $             +8*GVQ(2)*GAQ(2)*GVQ(2)*GAQ(2)*CTHD
            ZHMAX=(GVQ(2)**2+GAQ(2))*(GVQ(2)**2+GAQ(2)**2)*2
     $             +8*GVQ(2)*GAQ(2)*GVQ(2)*GAQ(2)
          ELSEIF (IDABS1.EQ.11.OR.IDABS1.EQ.13.OR.IDABS1.EQ.15) THEN
            ZHSIG=(GVQ(2)**2+GAQ(2))*(GVL(1)**2+GAL(1)**2)*(1.+CTHD**2)
     $             +8*GVQ(2)*GAQ(2)*GVL(1)*GAL(1)*CTHD
            ZHMAX=(GVQ(2)**2+GAQ(2))*(GVL(1)**2+GAL(1)**2)*2
     $             +8*GVQ(2)*GAQ(2)*GVL(1)*GAL(1)
          ELSEIF (IDABS1.EQ.12.OR.IDABS1.EQ.14.OR.IDABS1.EQ.16) THEN
            ZHSIG=(GVQ(2)**2+GAQ(2))*(GVL(2)**2+GAL(2)**2)*(1.+CTHD**2)
     $             +8*GVQ(2)*GAQ(2)*GVL(2)*GAL(2)*CTHD
            ZHMAX=(GVQ(2)**2+GAQ(2))*(GVL(2)**2+GAL(2)**2)*2
     $             +8*GVQ(2)*GAQ(2)*GVL(2)*GAL(2)
          END IF
        END IF      
        IF(ZHSIG.GT.RANF()*ZHMAX) GO TO 400
      END IF
      NREJ=NREJ+1
      IF(NREJ.LT.NTRIES) GO TO 300
      GO TO 999
C
C          Good event
C
400   CONTINUE
      RETURN
C
999   CALL PRTEVT(0)
      WRITE(ITLIS,9991) NREJ
9991  FORMAT(//' IT IS TAKING MORE THAN',I5,' TRIES TO GENERATE ',
     1'A GOOD WHIGGS EVENT.'/' CHECK LIMITS OR INCREASE NTRIES.')
      STOP 99
998   CALL PRTEVT(0)
      WRITE(ITLIS,9981) JET
9981  FORMAT(//' ERROR IN WHIGGS ... NO DECAY POSSIBLE FOR JET',I3)
      STOP 99
      END
