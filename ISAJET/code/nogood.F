#include "isajet/pilot.h"
      LOGICAL FUNCTION NOGOOD(KK)
C
C          Insure proper distribution and check kinematics.
C          Select jet types.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/keys.inc"
#include "isajet/wcon.inc"
#include "isajet/const.inc"
#include "isajet/wsig.inc"
#include "isajet/wgen.inc"
#include "isajet/dylim.inc"
#include "isajet/jetlim.inc"
#include "isajet/jetpar.inc"
#include "isajet/jetsig.inc"
#include "isajet/ptpar.inc"
#include "isajet/hcon.inc"
#include "isajet/xmssm.inc"
C
      REAL RANF,SIGINV,SUM,TRY,BRANCH
      INTEGER KK,I,II,K,IFL
C
      NOGOOD=.TRUE.
      GO TO (1,2,3,4,5,6),KK
C
C          TWOJET, SUPERSYM, WPAIR or PHOTON events
C
1     CONTINUE
      IF(KEYS(1)) THEN
        CALL SIGQCD
      ELSEIF(KEYS(5)) THEN
        CALL SIGSSY
      ELSEIF(KEYS(6)) THEN
        CALL SIGWW
      ELSEIF(KEYS(8)) THEN
        CALL SIGGAM
      ELSEIF(KEYS(10)) THEN
        CALL SIGWH
      ENDIF
      IF(SIGMA.LE.0) RETURN
      IF(SIGMAX*RANF().GT.SIGMA) RETURN
      NOGOOD=.FALSE.
      SIGINV=1./SIGMA
      SUM=0.
      TRY=RANF()
      DO 100 I=1,NSIGS
        SUM=SUM+SIGS(I)*SIGINV
        IF(SUM.LT.TRY) GO TO 100
C          Find reaction
        ISIGS=I
        SIGEVT=SIGS(ISIGS)
        II=INOUT(I)
        DO 110 K=1,2
        INITYP(K)=MOD(II,IOPAK)
110     II=II/IOPAK
        DO 120 K=1,2
        JETTYP(K)=MOD(II,IOPAK)
120     II=II/IOPAK
        RETURN
100   CONTINUE
      RETURN
C
C          DRELLYAN events--test of SIGDY
C
2     CONTINUE
      IF(KEYS(3)) THEN
        CALL SIGDY
      ELSEIF(KEYS(7).AND..NOT.GOMSSM) THEN
        CALL SIGH
      ELSEIF(KEYS(7).AND.GOMSSM) THEN
        CALL SIGHSS
      ELSEIF(KEYS(9)) THEN
        CALL SIGTC
      ELSEIF(KEYS(11)) THEN
        CALL SIGKKG
      ENDIF
      IF(SIGMA.LE.0.) RETURN
      IF(SIGSL(KSEL)*RANF().GT.SIGMA) RETURN
      NOGOOD=.FALSE.
      SIGINV=1./SIGMA
      SUM=0.
      TRY=RANF()
C          Find reaction.
      DO 200 I=1,NSIGS
        SUM=SUM+SIGS(I)*SIGINV
        IF(SUM.LT.TRY) GO TO 200
        ISIGS=I
        SIGEVT=SIGS(ISIGS)
        GO TO 210
200   CONTINUE
C          Unpack INOUT to find JETTYP and INITYP
210   IF(KEYS(3).OR.KEYS(11)) THEN
        II=INOUT(I)
        DO 220 K=1,2
        INITYP(K)=MOD(II,IOPAK)
220     II=II/IOPAK
        JWTYP=MOD(II,IOPAK)
        II=II/IOPAK
        JETTYP(3)=MOD(II,IOPAK)
      ELSEIF(KEYS(7).OR.KEYS(9)) THEN
        II=INOUT(ISIGS)
        DO 230 I=1,2
        INITYP(I)=MOD(II,IOPAK)
230     II=II/IOPAK
        DO 240 I=1,2
        JETTYP(I)=MOD(II,IOPAK)
240     II=II/IOPAK
      ENDIF
      RETURN
C
C          DRELLYAN events--test of SIGDY2
C
3     CONTINUE
      IF(KEYS(3)) THEN
        CALL SIGDY2
        IFL=JETTYP(1)/2
        BRANCH=(AQ(IFL,JWTYP)**2+BQ(IFL,JWTYP)**2)/COUT(JWTYP)
      ELSEIF(KEYS(7).AND..NOT.GOMSSM) THEN
        CALL SIGH2
        BRANCH=1.
      ELSEIF(KEYS(7).AND.GOMSSM) THEN
        SIGLLQ=SIGMA/(4*PI)
        NOGOOD=.FALSE.
        RETURN
      ELSEIF(KEYS(9)) THEN
        CALL SIGTC2
        BRANCH=1.
      ENDIF
      IF(SIGLLQ.GT.SIGS(ISIGS)*BRANCH*3.*RANF()/(4.*PI))
     1NOGOOD=.FALSE.
      RETURN
C
C          DRELLYAN events--test of kinematics
C
4     CONTINUE
      DO 400 I=1,2
        IF(P(I).LT.PMIN(I).OR.P(I).GT.PMAX(I)) GO TO 410
        IF(PT(I).LT.PTMIN(I).OR.PT(I).GT.PTMAX(I)) GO TO 410
        IF(YJ(I).LT.YJMIN(I).OR.YJ(I).GT.YJMAX(I)) GO TO 410
        IF(PHI(I).LT.PHIMIN(I).OR.PHI(I).GT.PHIMAX(I)) GO TO 410
400   CONTINUE
      NOGOOD=.FALSE.
410   RETURN
C
5     CONTINUE
6     CONTINUE
      RETURN
C
      END
