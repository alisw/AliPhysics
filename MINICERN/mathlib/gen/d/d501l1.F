*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE D501L1(VERS,SUB,K,N,X,NX,Y,SY,MODE,EPS,MAXIT,
     1                  IPRT,M,A,AL,AU,PHI1,DPHI,IAFR,MFR,
     2                  COV,NC,STD,P,LAMU,DSCAL,W1,W2,W3,TAU,
     3                  COPYF,COPYDF,R2,R1,F,DF,JPVT,NERROR)

#include "gen/imp64.inc"
#include "gen/def64.inc"
     +   JP2,LAMBDA,LAMU,LK,MY
      CHARACTER VERS*6
      LOGICAL LFN,LID,LRP,LPR
      DIMENSION X(*),Y(*),SY(*),A(*),AL(*),AU(*),DPHI(*),F(*),DF(N,*)
      DIMENSION STD(*),P(*),LAMU(*),DSCAL(*),W1(*),W2(*),W3(*),COV(NC,*)
      DIMENSION IAFR(*),TAU(*),JPVT(*),COPYF(*),COPYDF(N+M,*)
      DIMENSION R1(M,*),R2(M,*)
      DIMENSION W64(64)

      PARAMETER (Z0 = 0, Z1 = 1, HALF = Z1/2, R3 = Z1/3, R10 = Z1/10)
      PARAMETER (SIG1 = R10, SIG2 = 11*R10, COEF = R10**3, STEP = Z1)
      PARAMETER (RHO1 = R10**4, RHO2 = Z1/4, RHO3 = 3*Z1/4)

      EXTERNAL SUB

************************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   THIS ROUTINE IS ONE OF THE MAIN ROUTINE OF THE LEAMAX PACKAGE.
*   IT SOLVES TWO DIFFERENT PROBLEMS DEPENDING ON THE VALUE
*   OF THE PARAMETER VERS.
*   ( VERS = DSUMSQ : GENERAL NONLINEAR LEAST SQUARES PROBLEM
*     VERS = DFUNFT : LEAST SQUARES DATA FITTING PROBLEM      )
*
*   IN ALL CASES BOUNDS ON THE VARIABLES MAY BE SET.
*
************************************************************************

************************************************************************
*   COMPUTE AN APPROXIMATION  EPS0  TO THE RELATIVE MACHINE PRECISION
************************************************************************

      EPS0=Z1
    5 EPS0=EPS0/10
      IF (Z1+EPS0 .NE. Z1) GO TO 5
      EPS0=10*EPS0

************************************************************************
*   CHECK THE VAUES OF INPUT PARAMETERS
************************************************************************

      NERROR=0

      CALL D501P1(K,N,NC,X,NX,Y,SY,MODE,EPS0,EPS,MAXIT,IPRT,M,A,AL,AU,
     +            NERROR,VERS)

      IF (NERROR .NE. 0) RETURN

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

       CALL D501SF(VERS,SUB,0,M,A,N,F,DF,K,NX,X,Y,SY,W2,NERROR)
       IF (NERROR .NE. 0) RETURN
       PHI1=HALF*DVMPY(N,F(1),F(2),F(1),F(2))

************************************************************************
*   COMPUTE  F, DF, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N1(K,N,M,A,AL,AU,X,NX,Y,SY,W2,DPHI,DSCAL,LAMU,F,DF,IAFR,
     +            MFR,SUB,EPS0,EPS1,MODE,VERS,NERROR)
      IF(NERROR .NE. 0) RETURN

************************************************************************
*   IF MFR = 0 MINIMUM IN A CORNER; STOP ITERATION
************************************************************************

       IF(MFR .EQ. 0) GO TO 230

************************************************************************
*   ITERATION BEGINS
************************************************************************

      DELTA=0
      LAMBDA=0

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=0
      DO 10 I=1,MFR
   10 DPHINO=DPHINO+DPHI(IAFR(I))**2
      DPHINO=SQRT(DPHINO)

      IF(LPR) CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     +                    MODE,VERS)

      DA=0
      DO 20 I=1,MFR
   20 DA=DA+(DSCAL(I)*A(IAFR(I)))**2
      DA=SQRT(DA)

************************************************************************
*   ITERATION WITH GAUSS-NEWTON STEP
************************************************************************

   30 LAMBDA=0

************************************************************************
*   COPY  F  AND  DF
*   COMPUTE THE QR FACTORIZATION WITH COLUMN PIVOTING OF  DF
*   AND SOLVE THE LINEAR LEAST SQUARES PROBLEM USING
*   LAPACK ROUTINES  DGEQPF , DORMQR , DTRTRS
************************************************************************

      CALL DVSCL(N,-Z1,F(1),F(2),COPYF(1),COPYF(2))
      CALL DMCPY(N,MFR,DF(1,1),DF(1,2),DF(2,1),COPYDF(1,1),
     1           COPYDF(1,2),COPYDF(2,1))
C**** KSK 25.07.95
      DO 31 NN = 1,MFR
   31 JPVT(NN)=0
C     CALL DVSET(MFR,Z0,JPVT(1),JPVT(2))
C**** KSK 25.07.95

      CALL DGEQPF(N,MFR,COPYDF,N+M,JPVT,TAU,W3,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      CALL DORMQR('L','T',N,1,MFR,COPYDF,N+M,TAU,COPYF,N,W64,64,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      CALL DTRTRS('U','N','N',MFR,1,COPYDF,N+M,COPYF,N,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      DO 40 I=1,MFR
   40 P(JPVT(I))=COPYF(I)

************************************************************************
*   COMPUTE THE MATRIX  R2
************************************************************************

      DO  50 I=1,MFR
      DO  50 J=1,MFR
      IF (I .GT. J) THEN
       R1(I,J)=0
      ELSE
       R1(I,J)=COPYDF(I,J)
      ENDIF
   50 CONTINUE

      DO  60 I=1,MFR
      DO  60 J=1,MFR
   60 R2(I,J)=R1(I,JPVT(J))

      CALL DMMLT(MFR,MFR,MFR,R2(1,1),R2(2,1),R2(1,2),R2(1,1),R2(1,2),
     +           R2(2,1),COV(1,1),COV(1,2),COV(2,1),W2)

************************************************************************
*   COMPUTE THE L2-NORM OF THE SCALED VECTOR  P
************************************************************************

      DP=0
      DO 70 I=1,MFR
   70 DP=DP+(DSCAL(I)*P(I))**2
      DP=SQRT(DP)

************************************************************************
*   COMPUTE THE STEP SIZE  ALFA
************************************************************************

      ALFA=1
      DO 80 I=1,MFR
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
   80 CONTINUE

************************************************************************
*   COMPUTE INITIAL DELTA IF NECESSARY
************************************************************************

      IF(.NOT.LID) THEN
       DELTA=STEP*MAX(DA,DP/SIG2)
       LID=.TRUE.
      ENDIF
      IF(DELTA .LE. EPS*DA) GO TO 230

************************************************************************
*   CONTINUATION WITH GAUSS-NEWTON OR SWITCHING TO LEVENBERG-MARQUARDT?
************************************************************************

      IF(DP .GT. SIG2*DELTA) THEN

***********************************************************************
*   DO THE LEVENBERG - MARQUARDT STEP, (HEBDEN'S METHOD).
*   - COMPUTE THE LM - PARAMETER LAMBDA
*   - COMPUTE THE CORRESPONDING STEP P, ITS DP AND ALFA.
***********************************************************************

       DO 90 I=1,MFR
   90  STD(I)=-DPHI(IAFR(I))

       UK=0
       DO 100 I=1,MFR
  100  UK=UK+(DPHI(IAFR(I))/DSCAL(I))**2
       UK=SQRT(UK)/DELTA

************************************************************************
*   COMPUTE INITIAL LAMBDA
************************************************************************

       LAMBDA=COEF*UK

       LK=0
       ITERA=0

  110  ITERA=ITERA+1
       IF(ITERA .GE. 50) GO TO 230

************************************************************************
*   RESET LAMBDA IF NECESSARY
************************************************************************

       IF(LK .GE. LAMBDA .OR. LAMBDA .GE. UK)
     +    LAMBDA=MAX(COEF*UK,SQRT(LK*UK))

************************************************************************
*   COMPUTE NEW P FOR NEW LAMBDA
************************************************************************

************************************************************************
*   COPY F AND DF, AND EXTEND  DF  BY  SQRT(LAMBDA) * DIAG(DSCAL(I))
*   COMPUTE THE QR FACTORIZATION WITH COLUMN PIVOTING OF  EXTENDED  DF
*   AND SOLVE THE LINEAR LEAST SQUARES PROBLEM USING  LAPACK  ROUTINES
*   DGEQPF , DORMQR , DTRTRS
************************************************************************

      CALL DVSET(N+MFR,Z0,COPYF(1),COPYF(2))
      CALL DVSCL(N,-Z1,F(1),F(2),COPYF(1),COPYF(2))
      CALL DMSET(N+MFR,MFR,Z0,COPYDF(1,1),COPYDF(1,2),COPYDF(2,1))
      CALL DMCPY(N,MFR,DF(1,1),DF(1,2),DF(2,1),
     +           COPYDF(1,1),COPYDF(1,2),COPYDF(2,1))
      DO 120 I=1,MFR
  120 COPYDF(N+I,I)=SQRT(LAMBDA)*DSCAL(I)
C**** KSK 25.07.95
      DO 121 NN = 1,MFR
  121 JPVT(NN)=0
C     CALL DVSET(MFR,Z0,JPVT(1),JPVT(2))
C**** KSK 25.07.95

      CALL DGEQPF(N+MFR,MFR,COPYDF,N+M,JPVT,TAU,W3,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      CALL DORMQR('L','T',N+MFR,1,MFR,COPYDF,N+M,TAU,COPYF,
     +            N+MFR,W64,64,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      CALL DTRTRS('U','N','N',MFR,1,COPYDF,N+M,COPYF,N+MFR,INFO)
      IF (INFO .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      DO 130 I=1,MFR
  130 P(JPVT(I))=COPYF(I)


************************************************************************
*   STOP ITERATION?
************************************************************************

       DP=0
       DO 140 I=1,MFR
  140  DP=DP+(DSCAL(I)*P(I))**2
       DP=SQRT(DP)

       IF(SIG1*DELTA .GT. DP .OR. DP .GT. SIG2*DELTA) THEN

************************************************************************
*   CONTINUE ITERATION FOR LAMBDA
************************************************************************

        IF(DP .LE. 0) GO TO 230
        P1=DP-DELTA
        DO 150 I=1,MFR
  150   W1(I)=DSCAL(I)**2*P(I)

************************************************************************
*   COMPUTE THE MATRIX  R1
************************************************************************

      DO 160 I=1,MFR
      DO 160 J=1,MFR
      IF (I .GT. J) THEN
       R1(I,J)=0
      ELSE
       R1(I,J)=COPYDF(I,J)
      ENDIF
  160 CONTINUE

      DO 170 I=1,MFR
      DO 170 J=1,MFR
  170 R2(I,J)=R1(I,JPVT(J))

      CALL DMMLT(MFR,MFR,MFR,R2(1,1),R2(2,1),R2(1,2),R2(1,1),R2(1,2),
     +           R2(2,1),R1(1,1),R1(1,2),R1(2,1),W2)

      CALL DSINV(MFR,R1,M,NERROR)
      IF (NERROR .NE. 0) THEN
       NERROR=4
       RETURN
      ENDIF

      P1P=-DMBIL(MFR,W1(1),W1(2),R1(1,1),R1(1,2),R1(2,1),
     +               W1(1),W1(2))/DP

************************************************************************
*   UPDATE LK, UK, LAMBDA
************************************************************************

        IF(P1 .LT. 0) UK=LAMBDA
        LK=MAX(LK,LAMBDA-P1/P1P)
        IF(LK .GE. UK) UK=2*LK
        LAMBDA=LAMBDA-(DP/DELTA)*(P1/P1P)
        GO TO 110
       ENDIF
      ENDIF

************************************************************************
*   END OF LEVENBERG - MARQUARDT STEP
************************************************************************

      ALFA=1
      DO 180 I=1,MFR
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
  180 CONTINUE

************************************************************************
*   COMPUTE   A + ALPHA * P
************************************************************************

      CALL DVCPY(M,A(1),A(2),W1(1),W1(2))
      DO 190 I=1,MFR
  190 W1(IAFR(I))=A(IAFR(I))+ALFA*P(I)

************************************************************************
*   COMPUTE VALUE  PHI2  OF THE OBJECTIVE FUNCTION
************************************************************************

      CALL D501SF(VERS,SUB,0,M,W1,N,F,DF,K,NX,X,Y,SY,W2,NERROR)
      IF (NERROR .NE. 0) RETURN

      PHI2=HALF*DVMPY(N,F(1),F(2),F(1),F(2))

      PHMAXI=1
      IF(PHI1 .GT. 0) PHMAXI=1/SQRT(PHI1)
      CALL DVSCL(MFR,PHMAXI,P(1),P(2),W2(1),W2(2))
      JP2=DMBIL(MFR,W2(1),W2(2),COV(1,1),COV(1,2),COV(2,1),W2(1),W2(2))

************************************************************************
*   COMPUTE THE APPROXIMATION MEASURE  RHO  AND THE UPDATING FACTOR  MY
*   FOR DELTA
************************************************************************

       IF(PHI1 .LE. 0) THEN
        RHO=1
        MY=HALF
       ELSE
        S2=LAMBDA*DP**2/PHI1
        S3=1-PHI2/PHI1
        S4=HALF*JP2+S2
        IF(S4 .EQ. 0) THEN
         RHO=1
         MY=R10
        ELSE
         RHO=0
         IF(S3 .GT. 0) RHO=S3/S4
         MY=-HALF*(JP2+S2)
         S2=2*MY+S3
         IF(S2 .EQ. 0) THEN
          MY=R10
         ELSEIF(S3 .EQ. 0) THEN
          MY=HALF
         ELSE
          MY=MIN(MAX(MY/S2,R10),HALF)
         ENDIF
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
       DO 200 I=1,MFR
  200  DA=DA+(DSCAL(I)*A(IAFR(I)))**2
       DA=SQRT(DA)
       GO TO 30
      ENDIF
      CALL DVCPY(M,W1(1),W1(2),A(1),A(2))
      DA=0
      DO 210 I=1,MFR
  210 DA=DA+(DSCAL(I)*A(IAFR(I)))**2
      DA=SQRT(DA)
      MFROLD=MFR

************************************************************************
*   COMPUTE  F, DF, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N1(K,N,M,A,AL,AU,X,NX,Y,SY,W2,DPHI,DSCAL,LAMU,F,DF,IAFR,
     +            MFR,SUB,EPS0,EPS1,MODE,VERS,NERROR)
      IF(NERROR .NE. 0) RETURN

************************************************************************
*   IF MFR = 0  MINIMUM IN A CORNER; STOP ITERATION
************************************************************************

      IF(MFR .EQ. 0) THEN
       ITER=ITER+1
       GO TO 230
      ENDIF

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=0
      DO 220 I=1,MFR
  220 DPHINO=DPHINO+DPHI(IAFR(I))**2
      DPHINO=SQRT(DPHINO)

************************************************************************
*   TERMINATION CRITERION
************************************************************************

      IF (     PHI2      .LE. PHI1
     1   .AND. PHI1-PHI2 .LE. EPS*(1+ABS(PHI2))
     2   .AND. DP        .LE. SQRT(EPS)*(1+DA)
     3   .AND. DPHINO    .LE. EPS**R3*(1+ABS(PHI2))) LFN=.TRUE.

      ITER=ITER+1
      PHI1=PHI2

      IF(.NOT.LFN) THEN
       CALL D501SF(VERS,SUB,0,M,A,N,F,DF,K,NX,X,Y,SY,W2,NERROR)
       IF (NERROR .NE. 0) RETURN
       PHI1=HALF*DVMPY(N,F(1),F(2),F(1),F(2))

       IF(LPR) THEN
          IF((MOD(ITER,IPRT) .EQ. 0  .OR.  ITER .GE. MAXIT))
     1       CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     2                   MODE,VERS)
       ENDIF

       IF(ITER .GE. MAXIT) THEN
        NERROR=2
        GO TO 230
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
       GO TO 30
      ENDIF

************************************************************************
*   END OF ITERATION
************************************************************************

  230 LFN=.TRUE.

************************************************************************
*   COMPUTE  F, DF, DPHI, DSCAL, LAMU, MFR, IAFR
************************************************************************

      CALL D501N1(K,N,M,A,AL,AU,X,NX,Y,SY,W2,DPHI,DSCAL,LAMU,F,DF,IAFR,
     +            MFR,SUB,EPS0,EPS1,MODE,VERS,MERROR)
      IF(MERROR .NE. 0) THEN
       NERROR=MERROR
       RETURN
      ENDIF

************************************************************************
*   COMPUTE THE L2-NORM OF THE PROJECTED GRADIENT
************************************************************************

      DPHINO=0
      DO 240 I=1,MFR
  240 DPHINO=DPHINO+DPHI(IAFR(I))**2
      DPHINO=SQRT(DPHINO)

************************************************************************
*   COMPUTE THE VALUE PHI1 OF THE OBJECTIVE FUNCTION
************************************************************************

      PHI1=HALF*DVMPY(N,F(1),F(2),F(1),F(2))

************************************************************************
*   COMPUTE THE COVARIANCE MATRIX  COV  AND THE STANDARD DEVIATION  STD
*   FOR THE FREE VARIABLES
************************************************************************

      CALL DVSET(M,Z0,STD(1),STD(2))

      IF(MFR .GT. 0) THEN

       CALL DSINV(MFR,COV,NC,MERROR)
       IF(MERROR .NE. 0) THEN
        NERROR=4
        RETURN
       ENDIF

       S=2*PHI1
       IF(N .NE. MFR) S=S/(N-MFR)
       CALL DMSCL(MFR,MFR,S,COV(1,1),COV(1,2),COV(2,1),
     +                      COV(1,1),COV(1,2),COV(2,1))
       DO 250 I=1,MFR
  250  STD(IAFR(I))=SQRT(COV(I,I))
      ENDIF

************************************************************************
*   PRINT LAST ITERATION RESULTS
************************************************************************

      IF(LPR) CALL D501P2(LRP,M,A,DPHI,STD,LAMU,PHI1,DPHINO,ITER,LFN,
     +                    MODE,VERS)

      RETURN

      END



