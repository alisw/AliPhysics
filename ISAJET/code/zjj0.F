#include "isajet/pilot.h"
      SUBROUTINE ZJJ0
C-----------------------------------------------------------------------
C
C          Initialize MadGraph/Helas to generate Z + 2 jets.
C          Cross section routines from MadGraph:
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
#include "isajet/primar.inc"
#include "isajet/sstype.inc"
#include "isajet/mgkin.inc"
#include "isajet/mgcoms.inc"
#include "isajet/mgsigs.inc"
C
      INTEGER IMAD(6)
      INTEGER IQ1,IQ2,IQ4,IQ5,IFL1,IFL2,IM1,IM2,I,NEV,KK,II,J
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3)
      EQUIVALENCE (P1(0),PJETS8(0,1))
      EQUIVALENCE (P2(0),PJETS8(0,2))
      EQUIVALENCE (P3(0),PJETS8(0,3))
      EQUIVALENCE (P4(0),PJETS8(0,4))
      EQUIVALENCE (P5(0),PJETS8(0,5))
      REAL*8 WT,TERM,FJAC,SZJJ1,SZJJ2,SZJJ3,SZJJ4,SUM
      REAL*8 SZJJ5,SZJJ6,SZJJ7
      REAL X1,X2,XX,QQ,TMP(MXSIG8),QFCN,STRUC,ALQCD
      INTEGER IQ,IH,NPT,NNN,ISUM
C
C          Map Jettype/2 to MadGraph
      DATA IMAD/3,4,8,7,12,11/
C
C          Parton distributions
      QFCN(XX,IQ,IH)=STRUC(XX,QQ,IQ,IDIN(IH))/XX
C
C          Begin
C
      NPT=MAX(NSIGMA,100)
      WRITE(ITLIS,1000) NTRIES,NPT
1000  FORMAT(//' INITIALIZING CROSS SECTIONS WITH',I6,
     $' TRIES FOR',I6,' POINTS EACH:')
      FJAC=UNITS/SCM
      WTTOT8=0
      NSIG8=0
C
C          Cases 1,4: q1 q1b -> z q2 q2b
C
      AMJET8(3)=ZMASS
      DO 100 IFL1=1,5
        IM1=IMAD(IFL1)
        IQ1=2*IFL1
        IQ2=IQ1+1
        AMJET8(1)=FMASS(IM1)
        AMJET8(2)=FMASS(IM1)
        DO 110 IFL2=1,6
          IM2=IMAD(IFL2)
          IQ4=2*IFL2
          IQ5=IQ4+1
          AMJET8(4)=FMASS(IM2)
          AMJET8(5)=FMASS(IM2)
C
C          Subcase 1a: 3=z, 4=q, 5=qb
C
          IF(GOQ(IQ4,2).AND.GOQ(IQ5,3)) THEN
            IF(NSIG8+2.GT.MXSIG8) GO TO 999
            DO 120 I=1,2
              WTSUM8(NSIG8+I)=0
              WTMAX8(NSIG8+I)=0
              NWT8(NSIG8+I)=0
              IDENT8(1,NSIG8+I)=(3-2*I)*IFL1
              IDENT8(2,NSIG8+I)=-(3-2*I)*IFL1
              IDENT8(3,NSIG8+I)=IDZ
              IDENT8(4,NSIG8+I)=IFL2
              IDENT8(5,NSIG8+I)=-IFL2
              IF(IFL1.EQ.IFL2) THEN
                IFUNC8(NSIG8+I)=4
              ELSE
                IFUNC8(NSIG8+I)=1
              ENDIF
              NNN=0
              DO 125 NEV=1,NTRIES
                IF(NNN.GT.NPT) GO TO 120
                CALL MULJET(WT)
                NWT8(NSIG8+I)=NWT8(NSIG8+I)+1
                NWTTOT=NWTTOT+1
                IF(WT.GT.0) THEN
                  NNN=NNN+1
                  X1=(P1(0)+P1(3))/ECM
                  X2=(P2(0)-P2(3))/ECM
                  QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $            P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
                  IF(I.EQ.1) THEN
                    IF(IFL2.EQ.IFL1) THEN
                      TERM=SZJJ4(P1,P2,P3,P4,P5,IM1)
                    ELSE
                      TERM=SZJJ1(P1,P2,P3,P4,P5,IM1,IM2)
                    ENDIF
                    TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
                    TERM=TERM*WT*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
                  ELSE
                    IF(IFL2.EQ.IFL1) THEN
                      TERM=SZJJ4(P2,P1,P3,P4,P5,IM1)
                    ELSE
                      TERM=SZJJ1(P2,P1,P3,P4,P5,IM1,IM2)
                    ENDIF
                    TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
                    TERM=TERM*WT*FJAC*QFCN(X1,IQ2,1)*QFCN(X2,IQ1,2)
                  ENDIF
                  TERM=0.5*TERM
                  WTTOT8=WTTOT8+TERM
                  WTSUM8(NSIG8+I)=WTSUM8(NSIG8+I)+TERM
                  WTMAX8(NSIG8+I)=MAX(WTMAX8(NSIG8+I),TERM)
                ENDIF
125           CONTINUE
              WRITE(ITLIS,*) ' ZJJ0 WARNING: INSUFFICIENT TRIES FOR ',
     $        (IDENT8(KK,NSIG8+I),KK=1,5)
120         CONTINUE
            NSIG8=NSIG8+2
          ENDIF
C
C          Subcase 1b: 3=z, 4=qb, 5=q
C
          IF(GOQ(IQ5,2).AND.GOQ(IQ4,3)) THEN
            IF(NSIG8+2.GT.MXSIG8) GO TO 999
            DO 130 I=1,2
              WTSUM8(NSIG8+I)=0
              WTMAX8(NSIG8+I)=0
              NWT8(NSIG8+I)=0
              IDENT8(1,NSIG8+I)=(3-2*I)*IFL1
              IDENT8(2,NSIG8+I)=-(3-2*I)*IFL1
              IDENT8(3,NSIG8+I)=IDZ
              IDENT8(4,NSIG8+I)=-IFL2
              IDENT8(5,NSIG8+I)=IFL2
              IF(IFL1.EQ.IFL2) THEN
                IFUNC8(NSIG8+I)=4
              ELSE
                IFUNC8(NSIG8+I)=1
              ENDIF
              NNN=0
              DO 135 NEV=1,NTRIES
                IF(NNN.GT.NPT) GO TO 130
                CALL MULJET(WT)
                NWT8(NSIG8+I)=NWT8(NSIG8+I)+1
                NWTTOT=NWTTOT+1
                IF(WT.GT.0) THEN
                  NNN=NNN+1
                  X1=(P1(0)+P1(3))/ECM
                  X2=(P2(0)-P2(3))/ECM
                  QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $            P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
                  IF(I.EQ.1) THEN
                    IF(IFL1.EQ.IFL2) THEN
                      TERM=SZJJ4(P1,P2,P3,P5,P4,IM1)
                    ELSE
                      TERM=SZJJ1(P1,P2,P3,P5,P4,IM1,IM2)
                    ENDIF
                    TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
                    TERM=TERM*WT*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
                  ELSE
                    IF(IFL1.EQ.IFL2) THEN
                      TERM=SZJJ4(P2,P1,P3,P5,P4,IM1)
                    ELSE
                      TERM=SZJJ1(P2,P1,P3,P5,P4,IM1,IM2)
                    ENDIF
                    TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
                    TERM=TERM*WT*FJAC*QFCN(X1,IQ2,1)*QFCN(X2,IQ1,2)
                  ENDIF
                  TERM=0.5*TERM
                  WTTOT8=WTTOT8+TERM
                  WTSUM8(NSIG8+I)=WTSUM8(NSIG8+I)+TERM
                  WTMAX8(NSIG8+I)=MAX(WTMAX8(NSIG8+I),TERM)
                ENDIF
135           CONTINUE
              WRITE(ITLIS,*) ' ZJJ0 WARNING: INSUFFICIENT TRIES FOR ',
     $        (IDENT8(KK,NSIG8+I),KK=1,5)
130         CONTINUE
            NSIG8=NSIG8+2
          ENDIF
110     CONTINUE
100   CONTINUE
C
C          Case 2: g g -> z q2 q2b
C
      AMJET8(3)=ZMASS
      IFL1=9
      AMJET8(1)=0
      AMJET8(2)=0
      DO 210 IFL2=1,6
        IM2=IMAD(IFL2)
        IQ4=2*IFL2
        IQ5=IQ4+1
        AMJET8(4)=FMASS(IM2)
        AMJET8(5)=FMASS(IM2)
C
C          Subcase 2a: 3=z, 4=q, 5=qb
C
        IF(GOQ(IQ4,2).AND.GOQ(IQ5,3)) THEN
          IF(NSIG8+1.GT.MXSIG8) GO TO 999
          WTSUM8(NSIG8+1)=0
          WTMAX8(NSIG8+1)=0
          NWT8(NSIG8+1)=0
          IDENT8(1,NSIG8+1)=IDGL
          IDENT8(2,NSIG8+1)=IDGL
          IDENT8(3,NSIG8+1)=IDZ
          IDENT8(4,NSIG8+1)=IFL2
          IDENT8(5,NSIG8+1)=-IFL2
          IFUNC8(NSIG8+1)=2
          NNN=0
          DO 225 NEV=1,NTRIES
            IF(NNN.GT.NPT) GO TO 220
            CALL MULJET(WT)
            NWT8(NSIG8+1)=NWT8(NSIG8+1)+1
            NWTTOT=NWTTOT+1
            IF(WT.GT.0) THEN
              NNN=NNN+1
              X1=(P1(0)+P1(3))/ECM
              X2=(P2(0)-P2(3))/ECM
              QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $        P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
              TERM=SZJJ2(P1,P2,P3,P4,P5,IM2)
              TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
              TERM=TERM*WT*FJAC*QFCN(X1,1,1)*QFCN(X2,1,2)
              TERM=0.5*TERM
              WTTOT8=WTTOT8+TERM
              WTSUM8(NSIG8+1)=WTSUM8(NSIG8+1)+TERM
              WTMAX8(NSIG8+1)=MAX(WTMAX8(NSIG8+1),TERM)
            ENDIF
225       CONTINUE
          WRITE(ITLIS,*) ' ZJJ0 WARNING: INSUFFICIENT TRIES FOR ',
     $    (IDENT8(KK,NSIG8+1),KK=1,5)
220       CONTINUE
          NSIG8=NSIG8+1
        ENDIF
C
C          Subcase 2b: 3=z, 4=qb, 5=q
C
        IF(GOQ(IQ5,2).AND.GOQ(IQ4,3)) THEN
          IF(NSIG8+1.GT.MXSIG8) GO TO 999
          WTSUM8(NSIG8+1)=0
          WTMAX8(NSIG8+1)=0
          NWT8(NSIG8+1)=0
          IDENT8(1,NSIG8+1)=IDGL
          IDENT8(2,NSIG8+1)=IDGL
          IDENT8(3,NSIG8+1)=IDZ
          IDENT8(4,NSIG8+1)=-IFL2
          IDENT8(5,NSIG8+1)=IFL2
          IFUNC8(NSIG8+1)=2
          NNN=0
          DO 235 NEV=1,NTRIES
            IF(NNN.GT.NPT) GO TO 230
            CALL MULJET(WT)
            NWT8(NSIG8+1)=NWT8(NSIG8+1)+1
            NWTTOT=NWTTOT+1
            IF(WT.GT.0) THEN
              NNN=NNN+1
              X1=(P1(0)+P1(3))/ECM
              X2=(P2(0)-P2(3))/ECM
              QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $        P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
              TERM=SZJJ2(P1,P2,P3,P5,P4,IM2)
              TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
              TERM=TERM*WT*FJAC*QFCN(X1,1,1)*QFCN(X2,1,2)
              TERM=0.5*TERM
              WTTOT8=WTTOT8+TERM
              WTSUM8(NSIG8+1)=WTSUM8(NSIG8+1)+TERM
              WTMAX8(NSIG8+1)=MAX(WTMAX8(NSIG8+1),TERM)
            ENDIF
235       CONTINUE
          WRITE(ITLIS,*) ' ZJJ0 WARNING: INSUFFICIENT TRIES FOR ',
     $    (IDENT8(KK,NSIG8+1),KK=1,5)
230       CONTINUE
          NSIG8=NSIG8+1
        ENDIF
210   CONTINUE
C
C          Case 3: q1 q1b -> z g g
C
      AMJET8(3)=ZMASS
      AMJET8(4)=0
      AMJET8(5)=0
      DO 310 IFL1=1,5
        IM1=IMAD(IFL1)
        IQ1=2*IFL1
        IQ2=IQ1+1
        AMJET8(1)=FMASS(IM1)
        AMJET8(2)=FMASS(IM1)
C
        IF(GOQ(1,2).AND.GOQ(1,3)) THEN
          IF(NSIG8+2.GT.MXSIG8) GO TO 999
          DO 320 I=1,2
            WTSUM8(NSIG8+I)=0
            WTMAX8(NSIG8+I)=0
            NWT8(NSIG8+I)=0
            IDENT8(1,NSIG8+I)=(3-2*I)*IFL1
            IDENT8(2,NSIG8+I)=-(3-2*I)*IFL1
            IDENT8(3,NSIG8+I)=IDZ
            IDENT8(4,NSIG8+I)=IDGL
            IDENT8(5,NSIG8+I)=IDGL
            IFUNC8(NSIG8+I)=3
            NNN=0
            DO 325 NEV=1,NTRIES
              IF(NNN.GT.NPT) GO TO 320
              CALL MULJET(WT)
              NWT8(NSIG8+I)=NWT8(NSIG8+I)+1
              NWTTOT=NWTTOT+1
              IF(WT.GT.0) THEN
                NNN=NNN+1
                X1=(P1(0)+P1(3))/ECM
                X2=(P2(0)-P2(3))/ECM
                QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $          P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
                IF(I.EQ.1) THEN
                  TERM=SZJJ3(P1,P2,P3,P4,P5,IM1)
                  TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
                  TERM=TERM*WT*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
                ELSE
                  TERM=SZJJ3(P2,P1,P3,P4,P5,IM1)
                  TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
                  TERM=TERM*WT*FJAC*QFCN(X1,IQ2,1)*QFCN(X2,IQ1,2)
                ENDIF
                WTTOT8=WTTOT8+TERM
                WTSUM8(NSIG8+I)=WTSUM8(NSIG8+I)+TERM
                WTMAX8(NSIG8+I)=MAX(WTMAX8(NSIG8+I),TERM)
              ENDIF
325         CONTINUE
            WRITE(ITLIS,*) ' ZJJ0 WARNING: INSUFFICIENT TRIES FOR ',
     $      (IDENT8(KK,NSIG8+I),KK=1,5)
320       CONTINUE
          NSIG8=NSIG8+2
        ENDIF
310   CONTINUE
C
C          Cases 5,6: q1 q2 -> Z q1 q2, q1 != q2
C          Since we integrate over the Z decay, we can use the same
C          cross sections for quarks (I=1) and antiquarks (I=2).
C
      DO 400 IFL1=1,5
        IM1=IMAD(IFL1)
        DO 410 IFL2=1,5
          IM2=IMAD(IFL2)
          AMJET8(1)=FMASS(IM1)
          AMJET8(2)=FMASS(IM2)
          AMJET8(3)=ZMASS
          AMJET8(4)=FMASS(IM1)
          AMJET8(5)=FMASS(IM2)
          DO 420 I=1,2
            IQ1=2*IFL1+I-1
            IQ2=2*IFL2+I-1
            IF(GOQ(IQ1,1).AND.GOQ(IQ2,2)) THEN
              WTSUM8(NSIG8+1)=0
              WTMAX8(NSIG8+1)=0
              NWT8(NSIG8+1)=0
              IDENT8(1,NSIG8+1)=(3-2*I)*IFL1
              IDENT8(2,NSIG8+1)=(3-2*I)*IFL2
              IDENT8(3,NSIG8+1)=IDZ
              IDENT8(4,NSIG8+1)=(3-2*I)*IFL1
              IDENT8(5,NSIG8+1)=(3-2*I)*IFL2
              IF(IFL1.EQ.IFL2) THEN
                IFUNC8(NSIG8+1)=6
              ELSE
                IFUNC8(NSIG8+1)=5
              ENDIF
              NNN=0
              DO 425 NEV=1,NTRIES
                IF(NNN.GT.NPT) THEN
                  NSIG8=NSIG8+1
                  GO TO 420
                ENDIF
                CALL MULJET(WT)
                NWT8(NSIG8+1)=NWT8(NSIG8+1)+1
                NWTTOT=NWTTOT+1
                IF(WT.GT.0) THEN
                  NNN=NNN+1
                  X1=(P1(0)+P1(3))/ECM
                  X2=(P2(0)-P2(3))/ECM
                  QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $            P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
                  IF(IFL1.EQ.IFL2) THEN
                    TERM=SZJJ6(P1,P2,P3,P4,P5,IM1)
                    TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
                    TERM=TERM*WT*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
                  ELSE
                    TERM=SZJJ5(P2,P1,P3,P4,P5,IM1,IM2)
                    TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
                    TERM=TERM*WT*FJAC*QFCN(X1,IQ2,1)*QFCN(X2,IQ1,2)
                    TERM=0.5*TERM
                  ENDIF
                  WTTOT8=WTTOT8+TERM
                  WTSUM8(NSIG8+1)=WTSUM8(NSIG8+1)+TERM
                  WTMAX8(NSIG8+1)=MAX(WTMAX8(NSIG8+1),TERM)
                ENDIF
425           CONTINUE
              WRITE(ITLIS,*) ' ZJJ0 WARNING: INSUFFICIENT TRIES FOR ',
     $        (IDENT8(KK,NSIG8+1),KK=1,5)
              NSIG8=NSIG8+1
            ENDIF
420       CONTINUE
410     CONTINUE
400   CONTINUE
C
C          Case 7: g q -> z g q
C          Since we integrate over the Z decay, we can use the same
C          cross sections for quarks (I=1) and antiquarks (I=2).
C
      DO 500 IFL2=1,5
        IM2=IMAD(IFL2)
        DO 510 I=1,2
          IQ5=2*IFL2+I-1
C
C          Subcase 7a: 3=z, 4=g, 5=q (J=1,2 for initial states)
C
          IF(GOQ(1,2).AND.GOQ(IQ5,3)) THEN
            IF(NSIG8+2.GT.MXSIG8) GO TO 999
            AMJET8(3)=ZMASS
            AMJET8(4)=0
            AMJET8(5)=FMASS(IM2)
            DO 520 J=1,2
              WTSUM8(NSIG8+J)=0
              WTMAX8(NSIG8+J)=0
              NWT8(NSIG8+J)=0
              IF(J.EQ.1) THEN
                IDENT8(1,NSIG8+J)=IDGL
                IDENT8(2,NSIG8+J)=(3-2*I)*IFL2
                IQ1=1
                IQ2=IQ5
                AMJET8(1)=0
                AMJET8(2)=FMASS(IM2)
              ELSE
                IDENT8(2,NSIG8+J)=IDGL
                IDENT8(1,NSIG8+J)=(3-2*I)*IFL2
                IQ1=IQ5
                IQ2=1
                AMJET8(2)=0
                AMJET8(1)=FMASS(IM2)
              ENDIF
              IDENT8(3,NSIG8+J)=IDZ
              IDENT8(4,NSIG8+J)=IDGL
              IDENT8(5,NSIG8+J)=(3-2*I)*IFL2
              IFUNC8(NSIG8+J)=7
              NNN=0
              DO 525 NEV=1,NTRIES
                IF(NNN.GT.NPT) GO TO 520
                CALL MULJET(WT)
                NWT8(NSIG8+J)=NWT8(NSIG8+J)+1
                NWTTOT=NWTTOT+1
                IF(WT.GT.0) THEN
                  NNN=NNN+1
                  X1=(P1(0)+P1(3))/ECM
                  X2=(P2(0)-P2(3))/ECM
                  QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $            P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
                  IF(J.EQ.1) THEN
                    TERM=SZJJ7(P1,P2,P3,P4,P5,IM2)
                  ELSE
                    TERM=SZJJ7(P2,P1,P3,P4,P5,IM2)
                  ENDIF
                  TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
                  TERM=TERM*WT*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
                  TERM=0.5*TERM
                  WTTOT8=WTTOT8+TERM
                  WTSUM8(NSIG8+J)=WTSUM8(NSIG8+J)+TERM
                  WTMAX8(NSIG8+J)=MAX(WTMAX8(NSIG8+J),TERM)
                ENDIF
525           CONTINUE
              WRITE(ITLIS,*) ' ZJJ0 WARNING: INSUFFICIENT TRIES FOR ',
     $        (IDENT8(KK,NSIG8+1),KK=1,5)
520         CONTINUE
            NSIG8=NSIG8+2
          ENDIF
C
C          Subcase 7b: 3=z, 4=q, 5=g
C
          IF(GOQ(IQ5,2).AND.GOQ(1,3)) THEN
            IF(NSIG8+2.GT.MXSIG8) GO TO 999
            AMJET8(3)=ZMASS
            AMJET8(4)=FMASS(IM2)
            AMJET8(5)=0
            DO 530 J=1,2
              WTSUM8(NSIG8+J)=0
              WTMAX8(NSIG8+J)=0
              NWT8(NSIG8+J)=0
              IF(J.EQ.1) THEN
                IDENT8(1,NSIG8+J)=IDGL
                IDENT8(2,NSIG8+J)=(3-2*I)*IFL2
                IQ1=1
                IQ2=IQ5
                AMJET8(1)=0
                AMJET8(2)=FMASS(IM2)
              ELSE
                IDENT8(2,NSIG8+J)=IDGL
                IDENT8(1,NSIG8+J)=(3-2*I)*IFL2
                IQ1=IQ5
                IQ2=1
                AMJET8(2)=0
                AMJET8(1)=FMASS(IM2)
              ENDIF
              IDENT8(3,NSIG8+J)=IDZ
              IDENT8(4,NSIG8+J)=(3-2*I)*IFL2
              IDENT8(5,NSIG8+J)=IDGL
              IFUNC8(NSIG8+J)=7
              NNN=0
              DO 535 NEV=1,NTRIES
                IF(NNN.GT.NPT) GO TO 530
                CALL MULJET(WT)
                NWT8(NSIG8+J)=NWT8(NSIG8+J)+1
                NWTTOT=NWTTOT+1
                IF(WT.GT.0) THEN
                  NNN=NNN+1
                  X1=(P1(0)+P1(3))/ECM
                  X2=(P2(0)-P2(3))/ECM
                  QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $            P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
                  IF(J.EQ.1) THEN
                    TERM=SZJJ2(P1,P2,P3,P5,P4,IM2)
                  ELSE
                    TERM=SZJJ2(P2,P1,P3,P5,P4,IM2)
                  ENDIF
                  TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
                  TERM=TERM*WT*FJAC*QFCN(X1,1,1)*QFCN(X2,1,2)
                  TERM=0.5*TERM
                  WTTOT8=WTTOT8+TERM
                  WTSUM8(NSIG8+J)=WTSUM8(NSIG8+J)+TERM
                  WTMAX8(NSIG8+J)=MAX(WTMAX8(NSIG8+J),TERM)
                ENDIF
535           CONTINUE
              WRITE(ITLIS,*) ' ZJJ0 WARNING: INSUFFICIENT TRIES FOR ',
     $        (IDENT8(KK,NSIG8+1),KK=1,5)
530         CONTINUE
            NSIG8=NSIG8+2
          ENDIF
510     CONTINUE
500   CONTINUE
C
C          Sort using initial cross sections
C
      SUM=0
      ISUM=0
      DO 991 I=1,NSIG8
        ISORT8(I)=NSIG8+1-I
        TMP(I)=WTSUM8(I)/NWT8(I)
        SUM=SUM+WTSUM8(I)
        ISUM=ISUM+NWT8(I)
991   CONTINUE
      IF(NSIG8.GT.1) CALL SORTTF(TMP,ISORT8,NSIG8)
      WRITE(ITLIS,*)
      WRITE(ITLIS,9001)
9001  FORMAT(6X,'INITIAL MULTIJET CROSS SECTIONS'/
     $6X,'PROCESS',18X,'SIGMA',10X,'MAX(SIGMA)')
      DO 992 I=1,NSIG8
        II=ISORT8(I)
        WRITE(ITLIS,9002) (IDENT8(KK,II),KK=1,5),TMP(II),WTMAX8(II)
9002    FORMAT(2X,5I5,2E15.5)
992   CONTINUE
      RETURN
C
C          Errors
C
999   WRITE(ITLIS,*) 'ERROR IN ZJJ0, NSIG8 = ',NSIG8
      STOP 99
      END
