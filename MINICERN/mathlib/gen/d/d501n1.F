*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE D501N1(K,N,M,A,AL,AU,X,NX,Y,SY,WORK,DPHI,DSCAL,LAMU,
     +                  F,DF,IAFR,MFR,SUB,EPS0,EPS,MODE,VERS,NERROR)

*************************************************************************
*   LEAMAX, VERSION: 15.03.1993
*************************************************************************
*
*   THIS ROUTINE COMPUTES FUNCTION VALUES, DERIVATIVES, THE GRADIENT,
*   AND THE SCALING PARAMETERS. IT ALSO DETERMINES THE ACTIVE SET OF
*   CONSTRAINTS AND THE LAGRANGE MULTIPLIER.
*
************************************************************************

#include "gen/imp64.inc"
#include "gen/def64.inc"
     +    LAMU
      CHARACTER VERS*6
      DIMENSION A(*),AL(*),AU(*),X(*),Y(*),SY(*),WORK(*),DPHI(*)
      DIMENSION DSCAL(*),LAMU(*),F(*),DF(N,*),IAFR(*)
      EXTERNAL SUB

      PARAMETER (Z0 = 0)

************************************************************************
*   COMPUTE INITIAL VALUES
************************************************************************

      HREL=SQRT(EPS0)
      HABS=10*EPS0

      NERROR=0

************************************************************************
*   COMPUTE FUNCTION VALUES AND DERIVATIVES (IF MODE NOTEQUAL ZERO)
************************************************************************

      CALL D501SF(VERS,SUB,MODE,M,A,N,F,DF,K,NX,X,Y,SY,WORK(N+1),NERROR)
      IF(NERROR .NE. 0) RETURN

      IF(MODE .EQ. 0) THEN

************************************************************************
*    APPROXIMATE DERIVATIVES
************************************************************************

       DO 10 J=1,M
       H =ABS(A(J))*HREL+HABS
       IF (A(J)+H .GT. AU(J)) H=-H
       A(J)=A(J)+H
       CALL D501SF
     +      (VERS,SUB,MODE,M,A,N,WORK,DF,K,NX,X,Y,SY,WORK(N+1),NERROR)
       IF(NERROR .NE. 0) RETURN
       A(J)=A(J)-H
       CALL DVSUB(N,WORK(1),WORK(2),F(1),F(2),DF(1,J),DF(2,J))
   10  CALL DVSCL(N,1/H,DF(1,J),DF(2,J),DF(1,J),DF(2,J))
      ENDIF

************************************************************************
*   COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION
************************************************************************

      CALL DMMPY(M,N,DF(1,1),DF(2,1),DF(1,2),F(1),F(2),DPHI(1),DPHI(2))

************************************************************************
*   DETERMINE THE DIAGONAL MATRIX   DSCAL   FOR SCALING THE PROBLEM
************************************************************************

      DO 30 I=1,M
      AI=0
      DO 20 J=1,N
   20 AI=AI+DF(J,I)**2
   30 DSCAL(I)=MAX(DSCAL(I),SQRT(AI))

************************************************************************
*     DETERMINE FREE VARIABLES AND STORE THEIR INDECES IN IAFR
*     DETERMINE LAGRANGE-MULTIPLIER   LAMU
************************************************************************

      GR=0
      DO 40 I=1,MFR
   40 GR=GR+(DSCAL(I)*A(IAFR(I)))**2
      GR=HREL*SQRT(GR)

      CALL DVSET(M,Z0,LAMU(1),LAMU(2))

      MFR=0

      DO 50 I=1,M
      IF(AU(I)-AL(I) .LT. EPS*(ABS(AU(I))+ABS(AL(I)))+2*HABS) THEN
        A(I)=AU(I)
        LAMU(I)=DPHI(I)
      ELSE
       IF(A(I) .GE. AU(I)-(EPS * ABS(AU(I)) + HABS )) THEN
        A(I)=AU(I)
        IF(DPHI(I) .GT. -GR) THEN
         MFR=MFR+1
         IAFR(MFR)=I
        ELSE
         LAMU(I)=DPHI(I)
        ENDIF
       ELSE IF(A(I) .LE. AL(I)+(EPS * ABS(AL(I)) + HABS )) THEN
        A(I)=AL(I)
        IF(DPHI(I) .LT. GR) THEN
         MFR=MFR+1
         IAFR(MFR)=I
        ELSE
         LAMU(I)=DPHI(I)
        ENDIF
       ELSE
        MFR=MFR+1
        IAFR(MFR)=I
       ENDIF
      ENDIF

   50 CONTINUE

************************************************************************
*   DELETE ROWS OF  DSCAL  AND COLUMNS  OF  DF
*   WHICH BELONG TO NON-FREE VARIABLES
************************************************************************

       DO 60 I=1,MFR
       DSCAL(I)=DSCAL(IAFR(I))
       DO 60 L=1,N
   60  DF(L,I)=DF(L,IAFR(I))

      RETURN
      END




