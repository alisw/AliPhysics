*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:21  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE FUMILI (S,M,N1,N2,N3,EPS,AKAPPA,ALAMBD,IT,MC)
C-----ENTRY FOR CHISQ MINIMISATION
#include "d510pl.inc"
#include "d510si.inc"
#include "d510ui.inc"
#include "d510uo.inc"
C-----10.*MAXIMUM RELATIVE PRECISION ON CDC 6000
      DATA RP/1.E-14/
      INDFLG(3)=0
 1    IF (IT.GE.0) WRITE(6,84)
#if defined(CERNLIB_IBM)||defined(CERNLIB_CDC)
      CALL VFILL(AMX,100,1.0E75)
#endif
#if (!defined(CERNLIB_CDC))&&(!defined(CERNLIB_IBM))
      CALL VFILL(AMX,100,1.0E37)
#endif
      CALL VCOPYN(AMX,AMN,100)
      NN2=0
      N=M
      FIXFLG=0.
      ENDFLG=0.
      INDFLG(2)=0
      IFIX1=0.
      FI=0.
      NN3=0
      DO 2 I=1,N
      R(I)=0.
      IF (EPS.GT.0.) SIGMA(I)=0.
 2    PL(I)=PL0(I)
C-----START NEW ITERATION
 3    NN1=1
      T1=1.
C-----REPEAT ITERATION WITH SMALLER STEP
 4    S=0.
      N0=0
      DO 7 I=1,N
      G(I)=0.
      IF (PL0(I)) 7,7,5
 5    N0=N0+1
      IF (PL(I)) 7,7,6
 6    PL0(I)=PL(I)
 7    CONTINUE
      NN0=N0*(N0+1)/2
      IF (NN0.LT.1) GO TO 9
      DO 8 I=1,NN0
      Z(I)=0.
 8    CONTINUE
 9    NA=M
      INDFLG(1)=0
C-----CALCULATE OBJECTIVE FUNCTION
      CALL SGZ (M,S)
      SP=RP*ABS(S)
      IF (NN0.LT.1) GO TO 11
      DO 10 I=1,NN0
      Z0(I)=Z(I)
 10   CONTINUE
 11   IF (NN3) 19,19,12
 12   IF (NN1-N1) 13,13,19
 13   T=2.*(S-OLDS-GT)
      IF (INDFLG(1)) 16,14,16
 14   IF (ABS(S-OLDS).LE.SP.AND.-GT.LE.SP) GO TO 19
      IF (0.59*T+GT) 19,15,15
 15   T=-GT/T
      IF (T-0.25) 16,17,17
 16   T=0.25
 17   GT=GT*T
      T1=T1*T
      NN2=0
      DO 18 I=1,N
      IF (PL(I).LE.0.) GO TO 18
      A(I)=A(I)-DA(I)
      PL(I)=PL(I)*T
      DA(I)=DA(I)*T
      A(I)=A(I)+DA(I)
 18   CONTINUE
      NN1=NN1+1
      GO TO 4
C-----REMOVE CONTRIBUTION OF FIXED PARAMETERS FROM Z
   19 IF(INDFLG(1).EQ.0) GO TO 20
      ENDFLG=-4.
      GO TO 85
 20   K1=1
      K2=1
      I1=1
      DO 30 I=1,N
      IF (PL0(I)) 30,30,21
 21   IF (PL(I).EQ.0.) PL(I)=PL0(I)
      IF (PL(I)) 23,23,24
 22   PL(I)=0.
 23   K1=K1+I1
      GO TO 29
 24   IF (A(I).GE.AMX(I).AND.G(I).LT.0.) GO TO 22
      IF (A(I).LE.AMN(I).AND.G(I).GT.0.) GO TO 22
      DO 28 J=1,I
      IF (PL0(J)) 28,28,25
 25   IF (PL(J)) 27,27,26
 26   Z(K2)=Z0(K1)
      K2=K2+1
 27   K1=K1+1
 28   CONTINUE
 29   I1=I1+1
 30   CONTINUE
C-----INVERT Z
      I1=1
      L=I1
      DO 32 I=1,N
      IF (PL(I)) 32,32,31
 31   R(I)=Z(L)
      I1=I1+1
      L=L+I1
 32   CONTINUE
      N0=I1-1
      CALL MCONV (N0)
      IF (INDFLG(1)) 33,34,33
 33   INDFLG(1)=0
      INDFLG(2)=1
      GO TO 49
 34   CONTINUE
C-----CALCULATE THEORETICAL STEP TO MINIMUM
      I1=1
      DO 41 I=1,N
      DA(I)=0.
      IF (PL(I)) 41,41,35
 35   L1=1
      DO 40 L=1,N
      IF (PL(L)) 40,40,36
 36   IF (I1-L1) 37,37,38
 37   K=L1*(L1-1)/2+I1
      GO TO 39
 38   K=I1*(I1-1)/2+L1
 39   DA(I)=DA(I)-G(L)*Z(K)
      L1=L1+1
 40   CONTINUE
      I1=I1+1
 41   CONTINUE
C-----CHECK FOR PARAMETERS ON BOUNDARY
      AFIX=0.
      IFIX=0
      I1=1
      L=I1
      DO 47 I=1,N
      IF (PL(I)) 47,47,42
 42   SIGI=SQRT(ABS(Z(L)))
      R(I)=R(I)*Z(L)
      IF (EPS) 44,44,43
 43   SIGMA(I)=SIGI
 44   IF ((A(I).LT.AMX(I).OR.DA(I).LE.0.).AND.(A(I).GT.AMN(I).OR.DA(I).G
     1E.0.)) GO TO 46
      AKAP=ABS(DA(I)/SIGI)
      IF (AKAP-AFIX) 46,46,45
 45   AFIX=AKAP
      IFIX=I
      IFIX1=I
 46   I1=I1+1
      L=L+I1
 47   CONTINUE
      IF (IFIX) 48,50,48
 48   PL(IFIX)=-1.
 49   FIXFLG=FIXFLG+1.
      FI=0.
C-----REPEAT CALCULATION OF THEORETICAL STEP AFTER FIXING EACH PARAMETER
      GO TO 19
C-----CALCULATE STEP CORRECTION FACTOR
 50   ALAMBD=1.
      AKAPPA=0.
      IMAX=0
      DO 60 I=1,N
      IF (PL(I)) 60,60,51
 51   BM=AMX(I)-A(I)
      ABI=A(I)+PL(I)
      ABM=AMX(I)
      IF (DA(I)) 52,52,53
 52   BM=A(I)-AMN(I)
      ABI=A(I)-PL(I)
      ABM=AMN(I)
 53   BI=PL(I)
      IF (BI-BM) 55,55,54
 54   BI=BM
      ABI=ABM
 55   IF (ABS(DA(I))-BI) 58,58,56
 56   AL=ABS(BI/DA(I))
      IF (ALAMBD-AL) 58,58,57
 57   IMAX=I
      AIMAX=ABI
      ALAMBD=AL
 58   AKAP=ABS(DA(I)/SIGMA(I))
      IF (AKAP-AKAPPA) 60,60,59
 59   AKAPPA=AKAP
 60   CONTINUE
C-----CALCULATE NEW CORRECTED STEP
      GT=0.
      AMB=1.E18
      IF (ALAMBD) 62,62,61
 61   AMB=0.25/ALAMBD
 62   CONTINUE
      DO 67 I=1,N
      IF (PL(I)) 67,67,63
 63   IF (NN2-N2) 66,66,64
 64   IF (ABS(DA(I)/PL(I))-AMB) 66,65,65
 65   PL(I)=4.*PL(I)
      T1=4.
 66   DA(I)=DA(I)*ALAMBD
      GT=GT+DA(I)*G(I)
 67   CONTINUE
C-----CHECK IF MINIMUM ATTAINED AND SET EXIT MODE
      IF (-GT.GT.SP.OR.T1.GE.1..OR.ALAMBD.GE.1.) GO TO 68
      ENDFLG=-1.
 68   IF (ENDFLG) 85,69,69
 69   IF (AKAPPA-ABS(EPS)) 70,75,75
 70   IF (FIXFLG) 72,71,72
 71   ENDFLG=1.
      GO TO 85
 72   IF (ENDFLG) 85,77,73
 73   IF (IFIX1) 85,85,76
 74   IF (FI-FIXFLG) 76,76,77
 75   IF (FIXFLG) 74,76,74
 76   FI=FI+1.
      ENDFLG=0.
 85   IF(ENDFLG.EQ.0..AND.NN3.GE.N3) ENDFLG=-3.
      IF(ENDFLG.GT.0..AND.INDFLG(2).GT.0) ENDFLG=-2.
      CALL MONITO (S,M,NN3,IT,EPS,GT,AKAPPA,ALAMBD)
      IF (ENDFLG) 83,79,83
C-----CHECK IF FIXING ON BOUND IS CORRECT
 77   ENDFLG=1.
      FIXFLG=0.
      IFIX1=0
      DO 78 I=1,M
 78   PL(I)=PL0(I)
      INDFLG(2)=0
      GO TO 19
C-----NEXT ITERATION
 79   ENDFLG=0.
      DO 80 I=1,N
      A(I)=A(I)+DA(I)
 80   CONTINUE
      IF (IMAX) 82,82,81
 81   A(IMAX)=AIMAX
 82   OLDS=S
      NN2=NN2+1
      NN3=NN3+1
      GO TO 3
 83   MC=ENDFLG
      RETURN
C-----ENTRY FOR MAXIMUM LIKLEHOOD
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY LIKELM
#endif
#if !defined(CERNLIB_CDC)||!defined(CERNLIB_F4)
      ENTRY LIKELM(S,M,N1,N2,N3,EPS,AKAPPA,ALAMBD,IT,MC)
#endif
      INDFLG(3)=1
      GO TO 1
C
 84   FORMAT('1',43X,'FUNCTION MINIMISATION BY SUBROUTINE FUMILI/LIKE',
     +'LM'/'0',55X,'IN THE FOLLOWING PRINT-OUT'/
     +     '0',27X,'S = VALUE OF OBJECTIVE FUNCTION,',
     + 'EC = EXPECTED CHANGE IN S DURING NEXT ITERATION'/
     +   '0',34X,'KAPPA = ESTIMATED DISTANCE TO MINIMUM,  LAMBDA =',
     + 'STEP LENGTH MODIFIER'///)
      END
