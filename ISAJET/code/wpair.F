#include "isajet/pilot.h"
      SUBROUTINE WPAIR
C
C          Finish generation of wpair events started bY TWOJET.
C          Select W decay modes as allowed by WMODE1, WMODE2.
C          Generate W decay angles and put vectors in PPAIR.
C
C          Also generate massless decay vectors PZERO for matrix
C          element -- double precision for 32-bit machines.
C
C          Ver 6.26: Check kinematics for W -> ff decay, since Z0 from
C                    Higgs decay can be virtual.
C          Ver. 6.30: Added check in loop 201.
C          Ver. 7.14: Add MSSM Higgs hooks
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
#include "isajet/wwsig.inc"
#include "isajet/wwpar.inc"
#include "isajet/const.inc"
#include "isajet/qsave.inc"
#include "isajet/wcon.inc"
#include "isajet/pjets.inc"
#include "isajet/pinits.inc"
#include "isajet/keys.inc"
#include "isajet/wsig.inc"
#include "isajet/hcon.inc"
#include "isajet/xmssm.inc"
C
      DIMENSION X(2),LIST(25),P1WCM(4),P2WCM(4),P1LAB(4),P2LAB(4)
     $,P1CM0(4),P2CM0(4),P1LAB0(4),P2LAB0(4)
     1,PBOOST(4)
      EQUIVALENCE (X(1),X1)
      DIMENSION PWW(5,2)
      EQUIVALENCE (PWW(1,1),P3WW(1))
      DIMENSION JWWTYP(2),THWFF(2),PHIWFF(2)
#if defined(CERNLIB_SINGLE)
      REAL P1CM0,P2CM0,DPHI,DCTH,DSTH,DAM0,PWW,BOOST
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION P1CM0,P2CM0,DPHI,DCTH,DSTH,DAM0,PWW,BOOST
#endif
      REAL AMWW1,AMWW2,X,STRUC,STRUCW,RND,RANF,CBRWW,AMASS,AM0,AM1,AM2,
     $E1CM,E2CM,P12CM,CTHCM,STHCM,PHICM,CPHICM,SPHICM,P1WCM,P2WCM,
     $PBOOST,P1LAB,P2LAB,AFX,SGWWMX,P1LAB0,P2LAB0,THWFF,PHIWFF
      INTEGER IH,IQ,JWWTYP,JET,JWT,JQ,IQ1,IQ2,LIST,NREJ,NJ0,K
      REAL BRANCH(12),SUMBR
      INTEGER IDABS,IDABS1,IDABS2
C
      DATA LIST/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16/
C
C          Initialize for given W type.
      AMWW1=PJETS(5,1)
      AMWW2=PJETS(5,2)
      CALL WWKIN(AMWW1,AMWW2)
      NPAIR=0
C
C          Calculate and save structure functions.
      DO 120 IH=1,2
      DO 121 IQ=1,13
121   QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
      DO 122 IQ=14,25
122   QSAVE(IQ,IH)=0.
      IF(KEYS(7).OR.KEYS(9)) THEN
        DO 123 IQ=26,29
123     QSAVE(IQ,IH)=STRUCW(X(IH),IQ-25,IDIN(IH))/X(IH)
      ENDIF
120   CONTINUE
C          JWWTYP points to W types 1,2,3,4
      IF(KEYS(6)) THEN
        JWWTYP(1)=JETTYP(1)
        JWWTYP(2)=JETTYP(2)
      ELSEIF((KEYS(7).AND..NOT.GOMSSM).OR.KEYS(9)) THEN
        JWWTYP(1)=JETTYP(1)-25
        JWWTYP(2)=JETTYP(2)-25
      ELSEIF(KEYS(7).AND.GOMSSM) THEN
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
C          Also construct zero mass vectors at same angle
#if defined(CERNLIB_SINGLE)
C          Single precision.
        P1CM0(1)=.5*AM0*STHCM*CPHICM
        P2CM0(1)=-P1CM0(1)
        P1CM0(2)=.5*AM0*STHCM*SPHICM
        P2CM0(2)=-P1CM0(2)
        P1CM0(3)=.5*AM0*CTHCM
        P2CM0(3)=-P1CM0(3)
        P1CM0(4)=.5*AM0
        P2CM0(4)=P1CM0(4)
#endif
#if defined(CERNLIB_DOUBLE)
C          Double precision.
        DAM0=AM0
        DCTH=CTHCM
        DSTH=DSQRT(1.D0-DCTH**2)
        DPHI=PHICM
        P1CM0(1)=.5*AM0*DSTH*DCOS(DPHI)
        P2CM0(1)=-P1CM0(1)
        P1CM0(2)=.5*AM0*DSTH*DSIN(DPHI)
        P2CM0(2)=-P1CM0(2)
        P1CM0(3)=.5*AM0*DCTH
        P2CM0(3)=-P1CM0(3)
        P1CM0(4)=.5*AM0
        P2CM0(4)=P1CM0(4)
#endif
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
C          Boost zero mass vectors -- double precision for 32 bits.
        PZERO(4,NJ0-1)=(P1CM0(4)*PWW(4,JET)+P1CM0(1)*PWW(1,JET)
     $  +P1CM0(2)*PWW(2,JET)+P1CM0(3)*PWW(3,JET))/PWW(5,JET)
        BOOST=(P1CM0(4)+PZERO(4,NJ0-1))/(PWW(4,JET)+PWW(5,JET))
        DO 340 K=1,3
340     PZERO(K,NJ0-1)=P1CM0(K)+BOOST*PWW(K,JET)
        PZERO(4,NJ0)=(P2CM0(4)*PWW(4,JET)+P2CM0(1)*PWW(1,JET)
     $  +P2CM0(2)*PWW(2,JET)+P2CM0(3)*PWW(3,JET))/PWW(5,JET)
        BOOST=(P2CM0(4)+PZERO(4,NJ0))/(PWW(4,JET)+PWW(5,JET))
        DO 350 K=1,3
350     PZERO(K,NJ0)=P2CM0(K)+BOOST*PWW(K,JET)
        NJ0=NJ0+2
310   CONTINUE
C
C          Calculate cross section SIGWW2 containing TBRWW*RBRWW.
C          Compare with WW cross section containing TBRWW. Ratio
C          must be bounded by 3/(4*PI) for each W.
C
      AFX=3./(2.*PI)
      IF(KEYS(6)) THEN
        CALL SIGWW2
        SGWWMX=SIGEVT
        IF(IDJETS(1).NE.10) SGWWMX=SGWWMX*RBRWW(JQWW(1),JWWTYP(1),1)*AFX
        IF(IDJETS(2).NE.10) SGWWMX=SGWWMX*RBRWW(JQWW(2),JWWTYP(2),2)*AFX
      ELSEIF(KEYS(7)) THEN
C          Note that except for WW -> WW processes, SIGH3 just computes
C          the decay angular distribution, so it can be used for both 
C          for SM and SUSY HL0/HH0 decays; HA0 -> WW is forbidden.
C          For Z + HL0 decays, we just return, ie use phase space.
        IDABS1=IABS(IDJETS(1))
        IDABS2=IABS(IDJETS(2))
        IF(.NOT.(IDABS1.EQ.10.OR.IDABS1.EQ.80.OR.IDABS1.EQ.90)) RETURN
        IF(.NOT.(IDABS2.EQ.10.OR.IDABS2.EQ.80.OR.IDABS2.EQ.90)) RETURN
        CALL SIGH3
        SGWWMX=SIGLLQ*AFX**2
      ELSEIF(KEYS(9)) THEN
        CALL SIGTC3
        SGWWMX=SIGLLQ*AFX**2
      ENDIF
      IF(WWSIG.GT.SGWWMX*RANF()) GO TO 400
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
     1'A GOOD WPAIR EVENT.'/' CHECK LIMITS OR INCREASE NTRIES.')
      STOP 99
998   CALL PRTEVT(0)
      WRITE(ITLIS,9981) JET
9981  FORMAT(//' ERROR IN WPAIR ... NO DECAY POSSIBLE FOR JET',I3)
      STOP 99
      END
