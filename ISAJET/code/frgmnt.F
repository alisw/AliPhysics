#include "isajet/pilot.h"
      SUBROUTINE FRGMNT
C
C          Control jet fragmentation.  Boost to frames defined in
C          EVOLVE and call JETGEN.
C
C          EVOLVE initializes /PJSET/ as follows--
C          1      - 2            = PINITS (except for E+E-)
C          N0W    - N0W          = QWJET  (for DRELLYAN, NJET=3)
C          N0JETS - N0JETS+NJET  = PJETS
C          N0PAIR - N0PAIR+NPAIR = PPAIR  (for WPAIR)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/pjets.inc"
#include "isajet/pinits.inc"
#include "isajet/partcl.inc"
#include "isajet/const.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/keys.inc"
#include "isajet/q1q2.inc"
#include "isajet/frame.inc"
C
      REAL PSUM(5),PALLJ(5),P12(5),PIN(5,2),PWREST(5),PADD(5)
      REAL POLD(5),PNEW(5)
      REAL PINPL,PINMN,BP,PT2AVE,PTADD,RANF,PHIADD,PALLPL,PALLMN
      REAL PALLX,PALLY
      INTEGER K,J,JJET,NZERO,IB,NPTCL1,NPTCL2,IFAIL,JET,NPJET1,NPLV1
      INTEGER NPJET3,IP,NP1,NP2,NFIRST,IP1,IFR,NLJ
      DOUBLE PRECISION DSUM(5),DPASS(5)
C
C          Initialize
      DO 100 K=1,5
100   DSUM(K)=0.
      NLJ=NJET
      IF(KEYS(3)) NLJ=NJET+1
      DO 101 J=1,NLJ
        JJET=N0JETS+J-1
        IF(JJET.EQ.N0W) GOTO 101
        CALL DBLVEC(PJSET(1,JJET),DPASS)
        DO 102 K=1,4
102     DSUM(K)=DSUM(K)+DPASS(K)
101   CONTINUE
      DSUM(5)=DSQRT(DSUM(4)**2-DSUM(1)**2-DSUM(2)**2-DSUM(3)**2)
      DO 103 K=1,5
103   PALLJ(K)=DSUM(K)
C
      NZERO=NPTCL+1
C
C          Fragment partons from initial state shower
C
      IF(.NOT.KEYS(2)) THEN
        DO 110 J=1,NJSET
          IF(JDCAY(J).EQ.JPACK*J+J) THEN
            IB=JORIG(J)/JPACK-10
            DO 120 K=1,5
120         PIN(K,IB)=PJSET(K,J)
          ENDIF
110     CONTINUE
C
        CALL FRGJET(11)
        CALL FRGJET(12)
C
        NPTCL1=NPTCL+1
        NPTCL2=NPTCL1+1
        IF(NPTCL1.GT.MXPTCL) GO TO 9999
        PINPL=.5*(PIN(4,1)+PIN(3,1)+PIN(4,2)+PIN(3,2))
        PINMN=.5*(PIN(4,1)-PIN(3,1)+PIN(4,2)-PIN(3,2))
        PPTCL(1,NPTCL1)=0.
        PPTCL(2,NPTCL1)=0.
        PPTCL(3,NPTCL1)=HALFE-PINPL
        PPTCL(4,NPTCL1)=HALFE-PINPL
        PPTCL(5,NPTCL1)=0.
        PPTCL(1,NPTCL2)=0.
        PPTCL(2,NPTCL2)=0.
        PPTCL(3,NPTCL2)=-(HALFE-PINMN)
        PPTCL(4,NPTCL2)=HALFE-PINMN
        PPTCL(5,NPTCL2)=0.
        DO 130 K=1,4
130     PSUM(K)=-PALLJ(K)
        PSUM(4)=PSUM(4)+ECM
        PSUM(5)=PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2
        IF(PSUM(5).GE.0.) THEN
          PSUM(5)=SQRT(PSUM(5))
          CALL RESCAL(NZERO,NPTCL2,PSUM,IFAIL)
        ENDIF
C
        DO 140 K=1,4
140     PBEAMS(K)=PPTCL(K,NPTCL1)+PPTCL(K,NPTCL2)
        PBEAMS(5)=SQRT(PBEAMS(4)**2-PBEAMS(1)**2-PBEAMS(2)**2
     $  -PBEAMS(3)**2)
      ENDIF
C
C          Boost partons from final jets with -FRAME
C
200   DO 210 J=1,NJSET
        JET=JORIG(J)/JPACK
        IF ( JET.EQ.0 ) THEN
          IFR=1
        ELSE
          IF(JET.GT.10) GO TO 210
          IF(KEYS(6)) THEN
            IF(IDJETS(JET).EQ.10) GO TO 210
          ENDIF
          IFR=IFRAME(JET)
        ENDIF
C
C          Do this boost in double precision for 32-bit machines
        CALL DBOOST(-1,FRAME(1,IFR),PJSET(1,J))
210   CONTINUE
C
C          Fragment partons from final jets
C
      NPJET1=NPTCL+1
      DO 220 K=1,4
220   PSUM(K)=0
C
C          Conserve mass of 1+2 for DRELLYAN (automatic for WPAIR)
C
      IF(KEYS(3)) THEN
        CALL FRGJET(1)
        CALL FRGJET(2)
        IF(STDDY) THEN
          DO 230 K=1,4
            PSUM(K)=PJSET(K,3)+PJSET(K,4)
  230     CONTINUE
        ELSE
          DO 240 K=1,4
            PSUM(K)=PJSET(K,N0W+1)+PJSET(K,N0W+2)
  240     CONTINUE
        ENDIF
        PSUM(5)=SQRT(PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2)
        NPLV1=NPTCL
        CALL RESCAL(NPJET1,NPLV1,PSUM,IFAIL)
C          EXTRADIM has only jet3 + graviton
      ELSEIF(KEYS(11)) THEN
        CALL FRGJET(3)
        CALL FRGJET(0)
        NPLV1=NPTCL
        DO 241 K=1,4
          PSUM(K)=PJSET(K,3)+PJSET(K,4)
241     CONTINUE
        PSUM(5)=SQRT(PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2)
        CALL RESCAL (NPJET1,NPLV1,PSUM,IFAIL)
      ELSE
C          All other processes
        DO 242 J=1,NJET
          JJET=N0JETS+J-1
          CALL FRGJET(J)
          DO 243 K=1,4
  243     PSUM(K)=PSUM(K)+PJSET(K,JJET)
  242   CONTINUE
        PSUM(5)=SQRT(PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2)
        NPLV1=NPTCL
        CALL RESCAL(NPJET1,NPLV1,PSUM,IFAIL)
      ENDIF
C
C           Add extra jets for DRELLYAN
      IF(KEYS(3).AND..NOT.STDDY) THEN
        NPJET3=NPTCL+1
        DO 245 J=3,NJET
245     CALL FRGJET(J)
        NPTCL1=NPTCL+1
        IF(NPTCL1.GT.MXPTCL) GO TO 9999
        DO 250 K=1,4
          PPTCL(K,NPTCL1)=PJSET(K,N0W)
250     PSUM(K)=PJSET(K,N0W)
        PPTCL(5,NPTCL1)=PJSET(5,N0W)
        DO 246 J=3,NJET
          JJET=N0JETS+J-3
          DO 246 K=1,4
            PSUM(K)=PSUM(K)+PJSET(K,JJET)
  246   CONTINUE    
        PSUM(5)=SQRT(PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2)
        CALL RESCAL(NPJET3,NPTCL1,PSUM,IFAIL)
        DO 260 K=1,5
260     PWREST(K)=PPTCL(K,NPTCL1)
      ENDIF
C
C          Boost partons back to PP COM
C
      DO 300 J=1,NJSET
        JET=JORIG(J)/JPACK
        IF ( JET.EQ.0 ) THEN
          IFR=1
        ELSE
          IF(JET.GT.10) GO TO 300
          IF(KEYS(6)) THEN
            IF(IDJETS(JET).EQ.10) GO TO 300
          ENDIF
          IFR=IFRAME(JET)
        ENDIF
        BP=0.
        DO 305 K=1,3
305     BP=BP+FRAME(K,IFR)*PJSET(K,J)
        BP=BP/FRAME(5,IFR)
        DO 310 K=1,3
310     PJSET(K,J)=PJSET(K,J)+FRAME(K,IFR)*PJSET(4,J)/FRAME(5,IFR)
     $  +FRAME(K,IFR)*BP/(FRAME(4,IFR)+FRAME(5,IFR))
        PJSET(4,J)=FRAME(4,IFR)*PJSET(4,J)/FRAME(5,IFR)+BP
300   CONTINUE
C
C          Reset FRAME to boost hadrons to PP COM
C
      IF(KEYS(1).OR.KEYS(2).OR.(KEYS(3).AND.NJET.EQ.2).OR.KEYS(5)
     $.OR.(KEYS(7).AND.NPAIR.EQ.0).OR.KEYS(8)) THEN
        DO 410 K=1,5
          FRAME(K,1)=PALLJ(K)
410     CONTINUE
      ELSEIF(KEYS(3).AND.NJET.GT.2) THEN
        DO 420 K=1,5
420     FRAME(K,1)=PALLJ(K)
        BP=0.
        DO 430 K=1,3
430     BP=BP+FRAME(K,1)*PWREST(K)
        BP=BP/FRAME(5,1)
        DO 440 K=1,3
          FRAME(K,2)=PWREST(K)+FRAME(K,1)*PWREST(4)/FRAME(5,1)
     $    +FRAME(K,1)*BP/(FRAME(4,1)+FRAME(5,1))
440     CONTINUE
        FRAME(4,2)=FRAME(4,1)*PWREST(4)/FRAME(5,1)+BP
      ENDIF
C
C          Boost hadrons back to PP COM
C
      DO 500 IP=NZERO,NPTCL
        JET=IABS(IORIG(IP))/IPACK
        IF(JET.GT.10) GO TO 500
        IF(KEYS(6)) THEN
          IF(IDJETS(JET).EQ.10) GO TO 500
        ENDIF
        IF(JET.EQ.0) THEN
          IFR=1
        ELSE
          IFR=IFRAME(JET)
        ENDIF
        BP=0.
        DO 510 K=1,3
510     BP=BP+FRAME(K,IFR)*PPTCL(K,IP)
        BP=BP/FRAME(5,IFR)
        DO 520 K=1,3
520     PPTCL(K,IP)=PPTCL(K,IP)+FRAME(K,IFR)*PPTCL(4,IP)/FRAME(5,IFR)
     $  +FRAME(K,IFR)*BP/(FRAME(4,IFR)+FRAME(5,IFR))
        PPTCL(4,IP)=FRAME(4,IFR)*PPTCL(4,IP)/FRAME(5,IFR)+BP
500   CONTINUE
C
C          Add intrinsic PT
C
      IF(.NOT.KEYS(2)) THEN
        PT2AVE=.1*SQRT(QSQ)
        PTADD=SQRT(-PT2AVE*ALOG(RANF()))
        PHIADD=2.*PI*RANF()
        PADD(1)=2.*PTADD*COS(PHIADD)
        PADD(2)=2.*PTADD*SIN(PHIADD)
C          Must use large and small components carefully to calculate
C          mass on 32-bit machines.
        PALLPL=0.
        PALLMN=0.
        PALLX=0.
        PALLY=0.
        DO 525 IP=NZERO,NPTCL
          PALLX=PALLX+PPTCL(1,IP)
          PALLY=PALLY+PPTCL(2,IP)
          IF(PPTCL(3,IP).GT.0.) THEN
            PALLPL=PALLPL+(PPTCL(4,IP)+PPTCL(3,IP))
            PALLMN=PALLMN+(PPTCL(1,IP)**2+PPTCL(2,IP)**2+PPTCL(5,IP)**2)
     $      /(PPTCL(4,IP)+PPTCL(3,IP))
          ELSE
            PALLMN=PALLMN+(PPTCL(4,IP)-PPTCL(3,IP))
            PALLPL=PALLPL+(PPTCL(1,IP)**2+PPTCL(2,IP)**2+PPTCL(5,IP)**2)
     $      /(PPTCL(4,IP)-PPTCL(3,IP))
          ENDIF
525     CONTINUE
        POLD(1)=PALLX
        POLD(2)=PALLY
        POLD(3)=.5*(PALLPL-PALLMN)
        POLD(4)=.5*(PALLPL+PALLMN)
        POLD(5)=SQRT(PALLPL*PALLMN-PALLX**2-PALLY**2)
        PNEW(1)=PADD(1)+POLD(1)
        PNEW(2)=PADD(2)+POLD(2)
        PNEW(3)=POLD(3)
        PNEW(4)=SQRT(PNEW(1)**2+PNEW(2)**2+PNEW(3)**2+POLD(5)**2)
        PNEW(5)=POLD(5)
C
        DO 530 IP=NZERO,NPTCL
          BP=0.
          DO 531 K=1,3
531       BP=BP+POLD(K)*PPTCL(K,IP)
          BP=BP/POLD(5)
          DO 532 K=1,3
532       PPTCL(K,IP)=PPTCL(K,IP)-POLD(K)*PPTCL(4,IP)/POLD(5)
     $    +POLD(K)*BP/(POLD(4)+POLD(5))
          PPTCL(4,IP)=PPTCL(4,IP)*POLD(4)/POLD(5)-BP
C
          BP=0.
          DO 533 K=1,3
533       BP=BP+PNEW(K)*PPTCL(K,IP)
          BP=BP/PNEW(5)
          DO 534 K=1,3
534       PPTCL(K,IP)=PPTCL(K,IP)+PNEW(K)*PPTCL(4,IP)/PNEW(5)
     $    +PNEW(K)*BP/(PNEW(4)+PNEW(5))
          PPTCL(4,IP)=PPTCL(4,IP)*PNEW(4)/PNEW(5)+BP
530     CONTINUE
C
C            Add opposite PT to beam jets
        DO 541 K=1,4
541     PBEAMS(K)=-PNEW(K)
        PBEAMS(4)=PBEAMS(4)+ECM
        PBEAMS(5)=PBEAMS(4)**2-PBEAMS(1)**2-PBEAMS(2)**2 -PBEAMS(3)**2
        IF ( PBEAMS(5).GT.0 ) THEN
          PBEAMS(5)=SQRT(PBEAMS(5))
        ELSE
          PBEAMS(4)=SQRT(PBEAMS(4)**2-PBEAMS(5)+4.)
          PBEAMS(5)=2.
        ENDIF
      ENDIF
C
C          Decay hadrons
C
      NP1=NZERO
600   NP2=NPTCL
      DO 610 IP=NP1,NP2
        NFIRST=NPTCL+1
        JET=IABS(IORIG(IP))/IPACK
        CALL DECAY(IP)
        DO 620 IP1=NFIRST,NPTCL
620     IORIG(IP1)=ISIGN(IABS(IORIG(IP1))+IPACK*JET,IORIG(IP1))
610   CONTINUE
      NP1=NP2+1
      IF(NP1.LE.NPTCL) GO TO 600
      RETURN
C
C          Error
C
9999  CALL PRTEVT(0)
      WRITE(ITLIS,9998) NPTCL
9998  FORMAT(//' ERROR IN FRGMNT ... NPTCL > ',I6)
      RETURN
      END
