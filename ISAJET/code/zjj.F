#include "isajet/pilot.h"
      SUBROUTINE ZJJ
C-----------------------------------------------------------------------
C
C          Use MadGraph/Helas to generate Z + 2 jets after setup by
C          ZJJ0 using cross section routines from MadGraph:
C          ZJJ1:   q1 q1b -> Z q2 q2b, q1 != q2
C          ZJJ2:   g  g   -> Z q2 q2b
C          ZJJ3:   q1 q1b -> Z g  g
C          ZJJ4:   q1 q1b -> Z q1 q1b
C          ZJJ5:   q1 q2  -> Z q1 q2
C          ZJJ6:   q1 q1  -> Z q1 q1
C          ZJJ7:   g  q   -> Z g  q
C
C          Note: The Z is always jet1, but the other two jets are 
C          symmetrized so a symmetry factor of 1/2 is needed for every
C          subprocess. This is included by MadGraph for identical
C          particles!
C
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/const.inc"
#include "isajet/q1q2.inc"
#include "isajet/itapes.inc"
#include "isajet/jetlim.inc"
#include "isajet/jetpar.inc"
#include "isajet/jetset.inc"
#include "isajet/partcl.inc"
#include "isajet/pinits.inc"
#include "isajet/pjets.inc"
#include "isajet/primar.inc"
#include "isajet/sstype.inc"
#include "isajet/totals.inc"
#include "isajet/mgkin.inc"
#include "isajet/mgcoms.inc"
#include "isajet/mgsigs.inc"
C
      INTEGER IMAD(6)
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3)
      EQUIVALENCE (P1(0),PJETS8(0,1))
      EQUIVALENCE (P2(0),PJETS8(0,2))
      EQUIVALENCE (P3(0),PJETS8(0,3))
      EQUIVALENCE (P4(0),PJETS8(0,4))
      EQUIVALENCE (P5(0),PJETS8(0,5))
      REAL QFCN,XX,QQ,RND,RANF,SIG,FJAC,STRUC,ALQCD
      REAL*8 SZJJ1,WT8,TERM,SUM,SZJJ2,SZJJ3,SZJJ4,SIG8,SIGI8
      REAL*8 SZJJ5,SZJJ6,SZJJ7
      INTEGER IQ,IH,ISIG8,IFL1,IFL2,IM1,IM2,IQ1,IQ2,NTRY,I,II,K,IWT8
C
C          Map Jettype/2 to MadGraph
      DATA IMAD/3,4,8,7,12,11/
C
C          Parton distributions
      QFCN(XX,IQ,IH)=STRUC(XX,QQ,IQ,IDIN(IH))/XX
C
C          Begin
C
      NTRY=0
      NJSET=0
      NPTCL=0
C
C          Select process
C
      RND=RANF()
      ISIG8=0
      SIG=0
      DO 10 I=1,NSIG8
        SIG=SIG+WTSUM8(I)/NWT8(I)
10    CONTINUE
      SUM=0
      DO 20 I=1,NSIG8
        II=ISORT8(NSIG8+1-I)
        SUM=SUM+WTSUM8(II)/NWT8(II)
        IF(SUM.GE.RND*SIG) THEN
          ISIG8=II
          GO TO 100
        ENDIF
20    CONTINUE
      WRITE(ITLIS,*) 'ERROR IN ZJJ: NO MODE FOUND'
      STOP99
C
100   CONTINUE
      SIG8=0
      FJAC=UNITS/SCM
      NTRY=NTRY+1
      IF(NTRY.GT.NTRIES) THEN
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NTRY = ',NTRY
        WRITE(ITLIS,*) 'PROCESS WAS ',(IDENT8(K,ISIG8),K=1,5)
        SIGI8=WTSUM8(ISIG8)/NWT8(ISIG8)
        WRITE(ITLIS,*) 'PROCESS SIGMA/MAX = ',SIGI8,WTMAX8(ISIG8)
        WRITE(ITLIS,*) 'CHECK YOUR LIMITS OR INCREASE NTRIES'
        STOP99
      ENDIF
C
C          Cases 1,4: q1 q1b -> z q2 q2b
C
      IF(IFUNC8(ISIG8).EQ.1.OR.IFUNC8(ISIG8).EQ.4) THEN
        AMJET8(3)=ZMASS
        IFL1=IABS(IDENT8(1,ISIG8))
        IM1=IMAD(IFL1)
        IQ1=2*IFL1
        IQ2=IQ1+1
        AMJET8(1)=FMASS(IM1)
        AMJET8(2)=FMASS(IM1)
        IFL2=IABS(IDENT8(4,ISIG8))
        IM2=IMAD(IFL2)
        AMJET8(4)=FMASS(IM2)
        AMJET8(5)=FMASS(IM2)
        DO 210 I=1,NTRIES
          IWT8=I
          CALL MULJET(WT8)
          IF(WT8.GT.0) GO TO 220
210     CONTINUE
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NO PHASE SPACE POINT IN ',
     $  NTRIES,' TRIES'
        STOP99
220     NWTTOT=NWTTOT+IWT8-1
        NWT8(ISIG8)=NWT8(ISIG8)+IWT8-1
        X1=(P1(0)+P1(3))/ECM
        X2=(P2(0)-P2(3))/ECM
        QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $  P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
C
C          Subcases
C
        IF(IDENT8(1,ISIG8).GT.0.AND.IDENT8(4,ISIG8).GT.0) THEN
          IF(IFUNC8(ISIG8).EQ.1) THEN
            TERM=SZJJ1(P1,P2,P3,P4,P5,IM1,IM2)
          ELSE
            TERM=SZJJ4(P1,P2,P3,P4,P5,IM1)
          ENDIF
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(1,ISIG8).GT.0.AND.IDENT8(4,ISIG8).LT.0) THEN
          IF(IFUNC8(ISIG8).EQ.1) THEN
            TERM=SZJJ1(P1,P2,P3,P5,P4,IM1,IM2)
          ELSE
            TERM=SZJJ4(P1,P2,P3,P5,P4,IM1)
          ENDIF
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(1,ISIG8).LT.0.AND.IDENT8(4,ISIG8).GT.0) THEN
          IF(IFUNC8(ISIG8).EQ.1) THEN
            TERM=SZJJ1(P1,P2,P3,P4,P5,IM1,IM2)
          ELSE
            TERM=SZJJ4(P1,P2,P3,P4,P5,IM1)
          ENDIF
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ2,1)*QFCN(X2,IQ1,2)
        ELSEIF(IDENT8(1,ISIG8).LT.0.AND.IDENT8(4,ISIG8).LT.0) THEN
          IF(IFUNC8(ISIG8).EQ.1) THEN
            TERM=SZJJ1(P1,P2,P3,P5,P4,IM1,IM2)
          ELSE
            TERM=SZJJ4(P1,P2,P3,P5,P4,IM1)
          ENDIF
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ2,1)*QFCN(X2,IQ1,2)
        ELSE
          WRITE(ITLIS,*) 'ERROR IN ZJJ...INVALID FLAVOR FOR ZJJ1'
          STOP99
        ENDIF
        SIG8=0.5*TERM
        GO TO 900
      ENDIF
C
C          Case 2: g g -> z q2 q2b
C
      IF(IFUNC8(ISIG8).EQ.2) THEN
        AMJET8(3)=ZMASS
        IFL1=IABS(IDENT8(1,ISIG8))
        AMJET8(1)=0
        AMJET8(2)=0
        IFL2=IABS(IDENT8(4,ISIG8))
        IM2=IMAD(IFL2)
        AMJET8(4)=FMASS(IM2)
        AMJET8(5)=FMASS(IM2)
        DO 310 I=1,NTRIES
          IWT8=I
          CALL MULJET(WT8)
          IF(WT8.GT.0) GO TO 320
310     CONTINUE
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NO PHASE SPACE POINT IN ',
     $  NTRIES,' TRIES'
        STOP99
320     NWTTOT=NWTTOT+IWT8-1
        NWT8(ISIG8)=NWT8(ISIG8)+IWT8-1
        X1=(P1(0)+P1(3))/ECM
        X2=(P2(0)-P2(3))/ECM
        QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $  P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
C
C          Subcases
C
        IF(IDENT8(4,ISIG8).GT.0) THEN
          TERM=SZJJ2(P1,P2,P3,P4,P5,IM2)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,1,1)*QFCN(X2,1,2)
        ELSEIF(IDENT8(4,ISIG8).LT.0) THEN
          TERM=SZJJ2(P1,P2,P3,P5,P4,IM2)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,1,1)*QFCN(X2,1,2)
        ELSE
          WRITE(ITLIS,*) 'ERROR IN ZJJ...INVALID FLAVOR FOR ZJJ2'
          STOP99
        ENDIF
        SIG8=0.5*TERM
        GO TO 900
      ENDIF
C
C          Case 3: q1 q1b -> z g g
C
      IF(IFUNC8(ISIG8).EQ.3) THEN
        AMJET8(3)=ZMASS
        IFL1=IABS(IDENT8(1,ISIG8))
        IQ1=2*IFL1
        IQ2=IQ1+1
        IM1=IMAD(IFL1)
        AMJET8(1)=FMASS(IM1)
        AMJET8(2)=FMASS(IM1)
        IFL2=9
        AMJET8(4)=0
        AMJET8(5)=0
        DO 410 I=1,NTRIES
          IWT8=I
          CALL MULJET(WT8)
          IF(WT8.GT.0) GO TO 420
410     CONTINUE
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NO PHASE SPACE POINT IN ',
     $  NTRIES,' TRIES'
        STOP99
420     NWTTOT=NWTTOT+IWT8-1
        NWT8(ISIG8)=NWT8(ISIG8)+IWT8-1
        X1=(P1(0)+P1(3))/ECM
        X2=(P2(0)-P2(3))/ECM
        QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $  P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
C
C          Subcases
C
        IF(IDENT8(1,ISIG8).GT.0) THEN
          TERM=SZJJ3(P1,P2,P3,P4,P5,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(1,ISIG8).LT.0) THEN
          TERM=SZJJ3(P2,P1,P3,P4,P5,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ2,1)*QFCN(X2,IQ1,2)
        ELSE
          WRITE(ITLIS,*) 'ERROR IN ZJJ...INVALID FLAVOR FOR ZJJ3'
          STOP99
        ENDIF
        SIG8=TERM
        GO TO 900
      ENDIF
C
C          Cases 5,6: q1 q2 -> z q1 q2
C
      IF(IFUNC8(ISIG8).EQ.5.OR.IFUNC8(ISIG8).EQ.6) THEN
        IFL1=IABS(IDENT8(1,ISIG8))
        IM1=IMAD(IFL1)
        IFL2=IABS(IDENT8(2,ISIG8))
        IM2=IMAD(IFL2)
        IQ1=2*IFL1
        IQ2=2*IFL2
        IF(IDENT8(1,ISIG8).LT.0) IQ1=IQ1+1
        IF(IDENT8(2,ISIG8).LT.0) IQ2=IQ2+1
        AMJET8(1)=FMASS(IM1)
        AMJET8(2)=FMASS(IM2)
        AMJET8(3)=ZMASS
        AMJET8(4)=FMASS(IM1)
        AMJET8(5)=FMASS(IM2)
        DO 510 I=1,NTRIES
          IWT8=I
          CALL MULJET(WT8)
          IF(WT8.GT.0) GO TO 520
510     CONTINUE
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NO PHASE SPACE POINT IN ',
     $  NTRIES,' TRIES'
        STOP99
520     NWTTOT=NWTTOT+IWT8-1
        NWT8(ISIG8)=NWT8(ISIG8)+IWT8-1
        X1=(P1(0)+P1(3))/ECM
        X2=(P2(0)-P2(3))/ECM
        QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $  P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
C
C          Subcases
C
        IF(IDENT8(1,ISIG8).EQ.IDENT8(4,ISIG8)) THEN
          IF(IFUNC8(ISIG8).EQ.5) THEN
            TERM=SZJJ5(P1,P2,P3,P4,P5,IM1,IM2)
          ELSE
            TERM=SZJJ6(P1,P2,P3,P4,P5,IM1)
          ENDIF
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(1,ISIG8).EQ.IDENT8(5,ISIG8)) THEN
          TERM=SZJJ5(P1,P2,P3,P5,P4,IM1,IM2)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSE
          WRITE(ITLIS,*) 'ERROR IN ZJJ...INVALID FLAVOR FOR ZJJ1'
          STOP99
        ENDIF
        SIG8=TERM
        IF(IFL1.NE.IFL2) SIG8=0.5*SIG8
        GO TO 900
      ENDIF
C
C          Case 7: g q -> z g q
C
      IF(IFUNC8(ISIG8).EQ.7) THEN
        IF(IDENT8(1,ISIG8).EQ.9) THEN
          IFL1=IABS(IDENT8(2,ISIG8))
          IM1=IMAD(IFL1)
          AMJET8(1)=0
          AMJET8(2)=FMASS(IM1)
          IQ1=1
          IQ2=2*IFL1
          IF(IDENT8(2,ISIG8).LT.0) IQ2=IQ2+1
        ELSE
          IFL1=IABS(IDENT8(1,ISIG8))
          IM1=IMAD(IFL1)
          AMJET8(1)=FMASS(IM1)
          AMJET8(2)=0
          IQ2=1
          IQ1=2*IFL1
          IF(IDENT8(1,ISIG8).LT.0) IQ1=IQ1+1
        ENDIF
        AMJET8(3)=ZMASS
        IF(IDENT8(4,ISIG8).EQ.9) THEN
          AMJET8(4)=0
          AMJET8(5)=FMASS(IM1)
        ELSE
          AMJET8(4)=FMASS(IM1)
          AMJET8(5)=0
        ENDIF
        DO 610 I=1,NTRIES
          IWT8=I
          CALL MULJET(WT8)
          IF(WT8.GT.0) GO TO 620
610     CONTINUE
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NO PHASE SPACE POINT IN ',
     $  NTRIES,' TRIES'
        STOP99
620     NWTTOT=NWTTOT+IWT8-1
        NWT8(ISIG8)=NWT8(ISIG8)+IWT8-1
        X1=(P1(0)+P1(3))/ECM
        X2=(P2(0)-P2(3))/ECM
        QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $  P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
C
C          Subcases
C
        IF(IDENT8(1,ISIG8).EQ.9.AND.IDENT8(4,ISIG8).EQ.9) THEN
          TERM=SZJJ7(P1,P2,P3,P4,P5,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(2,ISIG8).EQ.9.AND.IDENT8(4,ISIG8).EQ.9) THEN
          TERM=SZJJ7(P2,P1,P3,P4,P5,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(1,ISIG8).EQ.9.AND.IDENT8(5,ISIG8).EQ.9) THEN
          TERM=SZJJ7(P1,P2,P3,P5,P4,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(2,ISIG8).EQ.9.AND.IDENT8(5,ISIG8).EQ.9) THEN
          TERM=SZJJ7(P2,P1,P3,P5,P4,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSE
          WRITE(ITLIS,*) 'ERROR IN ZJJ...INVALID FLAVOR FOR ZJJ1'
          STOP99
        ENDIF
        SIG8=0.5*TERM
        GO TO 900
      ENDIF
C
C          Increment totals and test
C
900   WTTOT8=WTTOT8+SIG8
      NWTTOT=NWTTOT+1
      WTSUM8(ISIG8)=WTSUM8(ISIG8)+SIG8
      WTMAX8(ISIG8)=MAX(WTMAX8(ISIG8),SIG8)
      NWT8(ISIG8)=NWT8(ISIG8)+1
      IF(SIG8.LT.RANF()*WTMAX8(ISIG8)) GO TO 100
C
C          Good event
C
      DO 910 I=1,3
        DO 911 K=1,3
          PJETS(K,I)=PJETS8(K,I+2)
911     CONTINUE
        PJETS(4,I)=PJETS8(0,I+2)
        PJETS(5,I)=AMJET8(I+2)
        IDJETS(I)=IDENT8(I+2,ISIG8)
910   CONTINUE
      DO 920 I=1,2
        DO 921 K=1,3
          PINITS(K,I)=PJETS8(K,I)
921     CONTINUE
        PINITS(4,I)=PJETS8(0,I)
        PINITS(5,I)=AMJET8(I)
        IDINIT(I)=IDENT8(I,ISIG8)
920   CONTINUE
C
      QSQ=QQ
      SHAT=(P1(0)+P2(0))**2-(P1(3)+P2(3))**2
      PBEAM(1)=(1.-X1)*HALFE
      PBEAM(2)=(1.-X2)*HALFE
C
C          Set /TOTALS/
C
      NKINPT=NWTTOT
      SUMWT=WTTOT8
C
      RETURN
      END
