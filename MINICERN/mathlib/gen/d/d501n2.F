*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE D501N2(K,N,M,A,AL,AU,X,NX,WORK,B,DPHI,DSCAL,LAMU,
     1                  AM,COV,IAFR,MFR,SUB,EPS0,EPS,MODE,NERROR)

#include "gen/imp64.inc"
#include "gen/def64.inc"
     +   LAMU
      DIMENSION A(*),AL(*),AU(*),X(*),WORK(*),B(*),DPHI(*),DSCAL(*)
      DIMENSION LAMU(*),AM(M,*),COV(M,*),IAFR(*)
      PARAMETER (Z0 = 0)

*************************************************************************
*   LEAMAX, VERSION: 15.03.1993
*************************************************************************
*
*   THIS ROUTINE COMPUTES THE GRADIENT, THE JACOBIAN, AND IT SETS UP
*   THE MATRIX FOR THE NORMAL EQUATIONS. IT ALSO DETERMINES THE ACTIVE
*   SET OF CONSTRAINTS AND THE LAGRANGE-MULTIPLIER.
*
************************************************************************

************************************************************************
*   SET INITIAL VALUES
************************************************************************

      HREL=SQRT(EPS0)
      HABS=10*EPS0

************************************************************************
*   COMPUTE THE GRADIENT   B  OF THE OBJECTIVE FUNCTION
*   COMPUTE AN APPROXIMATION  AM  OF THE SECOND DERIVATIVE (THE HESSIAN)
*   OF THE OBJECTIVE FUNCTION
************************************************************************

      NERROR=0
      CALL DVSET(M,Z0,B(1),B(2))
      CALL DMSET(M,M,Z0,AM(1,1),AM(1,2),AM(2,1))
      IX=1

      DO 30 I=1,N

      CALL SUB(K,X(IX),M,A,F0,WORK,MODE,NERROR)
      IF(NERROR .NE. 0  .OR.  F0 .LE. 0) THEN
       NERROR=3
       RETURN
      ENDIF

      IF(MODE .EQ. 0) THEN

************************************************************************
*   APPROXIMATE DERIVATIVES
************************************************************************

       DO 10 J=1,M
        H =ABS(A(J))*HREL+HABS
        IF (A(J)+H .GT. AU(J)) H =-H
        A(J)=A(J)+H
        CALL SUB(K,X(IX),M,A,FH,WORK,MODE,NERROR)
        IF(NERROR .NE. 0) THEN
         NERROR=3
         RETURN
        ENDIF
        A(J)=A(J)-H
   10   WORK(J)=(FH-F0)/H
       ENDIF

       CALL DVSCL(M,1/F0,WORK(1),WORK(2),WORK(1),WORK(2))
       CALL DVSUB(M,B(1),B(2),WORK(1),WORK(2),B(1),B(2))

       DO 20 L=1,M
       DO 20 J=L,M
   20  AM(L,J)=AM(L,J)+WORK(L)*WORK(J)

   30  IX=IX+NX

       CALL DMUTL(M,AM(1,1),AM(1,2),AM(2,1))

************************************************************************
*   COPY THE GRADIENT OF THE OBJECTIVE FUNCTION TO  DPHI
************************************************************************

      CALL DVCPY(M,B(1),B(2),DPHI(1),DPHI(2))

************************************************************************
*   DETERMINE THE DIAGONAL MATRIX  DSCAL  FOR SCALING THE PROBLEM
************************************************************************

      DO 40 I=1,M
   40 DSCAL(I)=MAX(DSCAL(I),SQRT(AM(I,I)))

************************************************************************
*   DETERMINE FREE VARIABLES AND STORE THEIR INDICES IN IAFR
*   DETERMINE LAGRANGE MULTIPLIER  LAMU
************************************************************************

      GR=0
      DO 50 I=1,MFR
   50 GR=GR+(DSCAL(I)*A(IAFR(I)))**2
      GR=HREL*SQRT(GR)
      CALL DVSET(M,Z0,LAMU(1),LAMU(2))

      MFR=0

      DO 60 I=1,M
      IF(AU(I)-AL(I) .LT. EPS*(ABS(AU(I))+ABS(AL(I)))+2*HABS) THEN
        A(I)=AU(I)
        LAMU(I)=DPHI(I)
      ELSE
       IF(A(I) .GE. AU(I)-(EPS * ABS(AU(I)) + HABS)) THEN
        A(I)=AU(I)
        IF(DPHI(I) .GT. -GR) THEN
         MFR=MFR+1
         IAFR(MFR)=I
        ELSE
         LAMU(I)=DPHI(I)
        ENDIF
       ELSE IF(A(I) .LE. AL(I)+(EPS * ABS(AL(I)) + HABS)) THEN
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

   60 CONTINUE

***********************************************************************
*   DELETE ROWS AND COLUMNS OF  AM  AND  B  WHICH BELONG TO NON-FREE
*   VARIABLES
************************************************************************

      IF(MFR .EQ. 0 .OR. MFR .EQ. M) THEN
       MFC=M
      ELSE
       MFC=MFR
       DO 70 I =1,MFR
       B(I)=B(IAFR(I))
       DSCAL(I)=DSCAL(IAFR(I))
       DO 70 L = 1,M
   70  AM(L,I)=AM(L,IAFR(I))
       DO 80 I=1,MFR
       DO 80 L=1,M
   80  AM(I,L)=AM(IAFR(I),L)
      ENDIF

      CALL DMCPY(MFC,MFC,AM(1,1),AM(1,2),AM(2,1),
     +                   COV(1,1),COV(1,2),COV(2,1))
      RETURN

      END



