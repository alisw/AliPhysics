*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE D501L2(K,N,X,NX,MODE,EPS,MAXIT,IPRT,M,A,AL,AU,
     1                  PHI1,DPHI,IAFR,MFR,STD,P,LAMU,DSCAL,B,W1,W2,
     2                  AM,COV,SUB,NERROR)

#include "gen/imp64.inc"
#include "gen/def64.inc"
     +    JP2,LAMBDA,LAMU,LK,MY
      LOGICAL LFN,LID,LRP,LPR
      DIMENSION X(*),A(*),AL(*),AU(*),DPHI(*),STD(*),P(*),LAMU(*)
      DIMENSION DSCAL(*),B(*),W1(*),W2(*),AM(M,*),COV(M,*),IAFR(*)
      PARAMETER (Z0 = 0, Z1 = 1, HALF = Z1/2, R3 = Z1/3, R10 = Z1/10)
      PARAMETER (SIG1 = R10, SIG2 = 11*R10, COEF = R10**3, STEP = Z1)
      PARAMETER (RHO1 = R10**4, RHO2 = Z1/4, RHO3 = 3*Z1/4)
      EXTERNAL SUB

************************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   THIS ROUTINE IS ONE OF THE MAIN ROUTINES OF THE  LEAMAX  PACKAGE
*   FOR SOLVING THE PROBLEM OF MAXIMUM LIKELIHOOD ESTIMATION.
*   BOUNDS ON THE VARIABLES MAY BE SET.
*
************************************************************************

************************************************************************
*   COMPUTE AN APPROXIMATION  EPS0  TO THE RELATIVE MACHINE PRECISION
************************************************************************

      EPS0=Z1
   10 EPS0=EPS0/10
      IF (Z1+EPS0 .NE. Z1) GO TO 10
      EPS0=10*EPS0

************************************************************************
*   CHECK THE VALUES OF INPUT PARAMETERS
************************************************************************

      NERROR=0

      CALL D501P1(K,N,M,X,NX,Y,SY,MODE,EPS0,EPS,MAXIT,IPRT,M,A,AL,AU,
     +            NERROR,'DMAXLK')

      IF(NERROR .NE. 0) RETURN

************************************************************************
*   SET INITIAL VALUES
************************************************************************

      EPS1=10*EPS0

      LFN=.FALSE.
      LID=.FALSE.
      LRP=.FALSE.
      LPR=IPRT .NE. 0

      ITER=0

      CALL DVSET(M,Z1,DSCAL(1),DSCAL(2))
      CALL DVSET(M,Z0,LAMU(1),LAMU(2))
      CALL DVSET(M,Z0,STD(1),STD(2))

************************************************************************
*   COMPUTE INITIAL VALUE  PHI1  OF OBJECTIVE FUNCTION
************************************************************************

      PHI1=0
      IX=1
      DO 20 I=1,N
      CALL SUB(K,X(IX),M,A,F0,ZZ,0,NERROR)
      IF(F0 .LE. 0  .OR.  NERROR .NE. 0) THEN
       NERROR=3
       RETURN
      ENDIF
      PHI1=PHI1-LOG(F0)
   20 IX=IX+NX

************************************************************************
*   COMPUTE J(TRANS) X J, B, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N2(K,N,M,A,AL,AU,X,NX,W1,B,DPHI,DSCAL,LAMU,AM,COV,
     +            IAFR,MFR,SUB,EPS0,EPS1,MODE,NERROR)
      IF(NERROR .NE. 0) RETURN

************************************************************************
*   IF  MFR = 0  MINIMUM IN A CORNER; STOP ITERATION
************************************************************************

      IF(MFR .EQ. 0) GO TO 190

************************************************************************
*   ITERATION BEGINS
************************************************************************

      DELTA=0
      LAMBDA=0

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=SQRT(DVMPY(MFR,B(1),B(2),B(1),B(2)))

      IF(LPR) CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     +                    MODE,'DMAXLK')

      DA=0
      DO 40 I=1,MFR
   40 DA=DA+(DSCAL(I)*A(IAFR(I)))**2
      DA=SQRT(DA)

************************************************************************
*   ITERATION WITH GAUSS-NEWTON STEP
************************************************************************

   50 LAMBDA=0

************************************************************************
*   SOLVE NORMAL EQUATIONS
************************************************************************

      CALL DVSCL(MFR,-Z1,B(1),B(2),W1(1),W1(2))
      CALL DSINV(MFR,AM,M,NERROR)

      IF(NERROR .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF
      CALL DMMPY(MFR,MFR,AM(1,1),AM(1,2),AM(2,1),W1(1),W1(2),P(1),P(2))


************************************************************************
*   COMPUTE THE L2-NORM OF THE SCALED VECTOR  P
************************************************************************

      DP=0
      DO 60 I=1,MFR
   60 DP=DP+(DSCAL(I)*P(I))**2
      DP=SQRT(DP)

************************************************************************
*   COMPUTE THE STEP SIZE  ALFA
************************************************************************

      ALFA=1
      DO 70 I=1,MFR
      IF(P(I) .NE. 0) THEN
       IF(P(I) .GT. 0) THEN
        ALFA1=AU(IAFR(I))
       ELSE
        ALFA1=AL(IAFR(I))
       ENDIF
       ALFA1=(ALFA1-A(IAFR(I)))/P(I)
       IF(ALFA1 .EQ. 0) THEN
        P(I)=0
       ELSE
        ALFA=MIN(ALFA,ALFA1)
       ENDIF
      ENDIF
   70 CONTINUE

************************************************************************
*   COMPUTE INITIAL DELTA IF NECESSARY
************************************************************************

      IF(.NOT.LID) THEN
       DELTA=STEP*MAX(DA,DP/SIG2)
       LID=.TRUE.
      ENDIF
      IF(DELTA .LE. EPS*DA) GO TO 190

************************************************************************
*   CONTINUATION WITH GAUSS-NEWTON OR SWITCHING TO LEVENBERG-MARQUARDT?
************************************************************************

      IF(DP .GT. SIG2*DELTA) THEN

************************************************************************
*   DO THE LEVENBERG - MARQUARDT STEP, (HEBDEN'S METHOD).
*   - COMPUTE THE LM - PARAMETER LAMBDA
*   - COMPUTE THE CORRESPONDING STEP P, ITS DP AND ALFA.
************************************************************************

       CALL DVSCL(MFR,-Z1,B(1),B(2),STD(1),STD(2))
       UK=0
       DO 80 I=1,MFR
   80  UK=UK+(DPHI(IAFR(I))/DSCAL(I))**2
       UK=SQRT(UK)/DELTA

************************************************************************
*   COMPUTE INITIAL LAMBDA
************************************************************************

       LAMBDA=COEF*UK

       LK=0
       ITERA=0

   90  ITERA=ITERA+1
       IF(ITERA .GE. 50) GO TO 190

************************************************************************
*   RESET LAMBDA IF NECESSARY
************************************************************************

       IF(LK .GE. LAMBDA .OR. LAMBDA .GE. UK)
     +    LAMBDA=MAX(COEF*UK,SQRT(LK*UK))

************************************************************************
*   COMPUTE NEW P FOR NEW LAMBDA
************************************************************************

       CALL DMCPY(MFR,MFR,COV(1,1),COV(1,2),COV(2,1),
     +                    AM(1,1),AM(1,2),AM(2,1))
       DO 100 I=1,MFR
  100  AM(I,I)=AM(I,I)+LAMBDA*DSCAL(I)**2

************************************************************************
*   SOLVE NORMAL EQUATIONS
************************************************************************

       CALL DSINV(MFR,AM,M,NERROR)
       IF(NERROR .NE. 0) THEN
        NERROR=4
        RETURN
       ENDIF
       CALL DMMPY(MFR,MFR,AM(1,1),AM(1,2),AM(2,1),STD(1),STD(2),
     +                    P(1),P(2))

************************************************************************
*   COMPUTE THE L2-NORM OF THE SCALED VECTOR  P
************************************************************************

       DP=0
       DO 110 I=1,MFR
  110  DP=DP+(DSCAL(I)*P(I))**2
       DP=SQRT(DP)

************************************************************************
*   STOP ITERATION IN THE CASE OF NORM EQUAL TO ZERO
************************************************************************

       IF (DP .LE. 0) GO TO 190

       IF(SIG1*DELTA .GT. DP .OR. DP .GT. SIG2*DELTA) THEN

************************************************************************
*   CONTINUE ITERATION FOR LAMBDA
************************************************************************

        P1=DP-DELTA
        DO 120 I=1,MFR
  120   W1(I)=DSCAL(I)**2*P(I)
        P1P=-DMBIL(MFR,W1(1),W1(2),AM(1,1),AM(1,2),AM(2,1),
     +                 W1(1),W1(2))/DP

************************************************************************
*   UPDATE LK, UK, LAMBDA
************************************************************************

        IF(P1 .LT. 0) UK=LAMBDA
        LK=MAX(LK,LAMBDA-P1/P1P)
        IF(LK .GE. UK) UK=2*LK
        LAMBDA=LAMBDA-(DP/DELTA)*(P1/P1P)
        GO TO 90
       ENDIF
      ENDIF

************************************************************************
*   END OF LEVENBERG - MARQUARDT STEP
************************************************************************

      ALFA=1
      DO 130 I=1,MFR
      IF(P(I) .NE. 0) THEN
       IF(P(I) .GT. 0) THEN
        ALFA1=AU(IAFR(I))
       ELSE
        ALFA1=AL(IAFR(I))
       ENDIF
       ALFA1=(ALFA1-A(IAFR(I)))/P(I)
       IF(ALFA1 .EQ. 0) THEN
        P(I)=0
       ELSE
        ALFA=MIN(ALFA,ALFA1)
       ENDIF
      ENDIF
  130 CONTINUE

************************************************************************
*   COMPUTE  A + ALPHA * P
************************************************************************
      CALL DVCPY(M,A(1),A(2),W1(1),W1(2))
      DO 140 I=1,MFR
  140 W1(IAFR(I))=A(IAFR(I))+ALFA*P(I)

************************************************************************
*   COMPUTE VALUE  PHI2  OF OBJECTIVE FUNCTION
************************************************************************

      PHI2=0
      IX=1
      DO 150 I=1,N
      CALL SUB(K,X(IX),M,W1,F0,ZZ,0,NERROR)
      IF(F0 .LE. 0  .OR.  NERROR .NE. 0) THEN
       NERROR=3
       RETURN
      ENDIF
      PHI2=PHI2-LOG(F0)
  150 IX=IX+NX

      PHMAXI=1
      AAU=MAX(ABS(PHI1),ABS(PHI2))
      IF(AAU .GT. 0) PHMAXI=1/AAU

      CALL DVSCL(MFR,PHMAXI,P(1),P(2),W2(1),W2(2))
      JP2=DMBIL(MFR,W2(1),W2(2),COV(1,1),COV(1,2),COV(2,1),W2(1),W2(2))

************************************************************************
*   COMPUTE THE APPROXIMATION MEASURE  RHO  AND THE UPDATING FACTOR  MY
*   FOR  DELTA
************************************************************************

       IF(PHI1 .LE. PHI2) THEN
        RHO=0
        MY=R10
       ELSE
        S2=2*(PHI1-PHI2)*PHMAXI**2
        S3=LAMBDA*(DP*PHMAXI)**2
        RHO=S2/(JP2+2*S3)
        MY=-JP2-S3
        S2=S2+2*MY
        IF(S2 .EQ. 0) THEN
         MY=R10
        ELSE
         MY=MIN(MAX(MY/S2,R10),HALF)
        ENDIF
       ENDIF

************************************************************************
*   END OF COMPUTATTION OF RHO AND MY
************************************************************************

************************************************************************
*   IF RHO .LE. RHO1, REDUCE DELTA BY FACTOR MY AND MAKE NEW LEVENBERG-
*   MARQUARDT STEP, OTHERWISE ACCEPT P
************************************************************************

      IF(RHO .LE. RHO1) THEN
       DELTA=MY*DELTA
       DA=0
       DO 160 I=1,MFR
  160  DA=DA+(DSCAL(I)*A(IAFR(I)))**2
       DA=SQRT(DA)
       GO TO 50
      ENDIF
      CALL DVCPY(M,W1(1),W1(2),A(1),A(2))
      DA=0
      DO 170 I=1,MFR
  170 DA=DA+(DSCAL(I)*A(IAFR(I)))**2
      DA=SQRT(DA)
      MFROLD=MFR

************************************************************************
*   COMPUTE J(TRANS) X J, B, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N2(K,N,M,A,AL,AU,X,NX,W1,B,DPHI,DSCAL,LAMU,AM,COV,
     +            IAFR,MFR,SUB,EPS0,EPS1,MODE,NERROR)
      IF(NERROR .NE. 0) RETURN

************************************************************************
*   IF  MFR = 0  MINIMUM IN A CORNER; STOP ITERATION
************************************************************************

      IF(MFR .EQ. 0) THEN
       ITER=ITER+1
       GO TO 190
      ENDIF

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=SQRT(DVMPY(MFR,B(1),B(2),B(1),B(2)))

************************************************************************
*   TERMINATION CRITERION
************************************************************************

      IF (     PHI2      .LE. PHI1
     1   .AND. PHI1-PHI2 .LE. EPS*(1+ABS(PHI2))
     2   .AND. DP        .LE. SQRT(EPS)*(1+DA)
     3   .AND. DPHINO    .LE. EPS**R3*(1+ABS(PHI2)))    LFN=.TRUE.

      ITER=ITER+1
      PHI1=PHI2

      IF(.NOT.LFN) THEN
       IF(ITER .GE. MAXIT) THEN
        IF(LPR) CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     +                      MODE,'DMAXLK')
        NERROR=2
        GO TO 190
       ENDIF

       IF(LPR) THEN
          IF(MOD(ITER,IPRT) .EQ. 0)
     1       CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     2                   MODE,'DMAXLK')
       ENDIF

************************************************************************
*   UPDATE DELTA AND GO BACK TO GAUSS-NEWTON STEP
************************************************************************

       IF(MFROLD .NE. MFR) LID=.FALSE.
       IF(RHO .LE. RHO2) THEN
        DELTA=MY*DELTA
       ELSE IF(RHO .GE. RHO3 .OR. LAMBDA .EQ. 0) THEN
        DELTA=2*DP
       ENDIF

       GO TO 50

      ENDIF

************************************************************************
*   END OF ITERATION
************************************************************************

  190 LFN=.TRUE.

************************************************************************
*   COMPUTE J(TRANS) X J, B, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N2(K,N,M,A,AL,AU,X,NX,W1,B,DPHI,DSCAL,LAMU,AM,COV,
     +            IAFR,MFR,SUB,EPS0,EPS1,MODE,MERROR)
      IF(MERROR .NE. 0) THEN
       NERROR=MERROR
       RETURN
      ENDIF

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=SQRT(DVMPY(MFR,B(1),B(2),B(1),B(2)))

************************************************************************
*   COMPUTE THE VALUE  PHI1  OF THE OBJECTIVE FUNCTION
************************************************************************

      PHI1=0
      IX=1
      DO 200 I=1,N
       CALL SUB(K,X(IX),M,A,F0,ZZ,0,MERROR)
       IF(F0 .LE. 0  .OR.  MERROR .NE. 0) THEN
        NERROR=3
        RETURN
       ENDIF
       PHI1=PHI1-LOG(F0)
  200 IX=IX+NX


************************************************************************
*   PRINT LAST ITERATION RESULTS
************************************************************************

      IF(LPR) CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     +                    MODE,'DMAXLK')

      RETURN

      END



