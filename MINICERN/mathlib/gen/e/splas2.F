*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:27  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE SPLAS2(N,NC,NCX,NCY,NX,NY,MX,MY,KX,KY,NDIMC,NDIMZ,
     +                  XI,YI,ZI,KNOT,TX,TY,A,S,VT,W1,W2,LW,C,NERR)

#include "gen/imp64.inc"
      DIMENSION XI(*),YI(*),ZI(NDIMZ,*),TX(*),TY(*)
      DIMENSION A(N,*),S(*),VT(NC,*),W1(*),W2(*),C(NDIMC,*)

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   THE SUBROUTINE  SPLAS2  IS USED BY  DSPAP2  FOR COMPUTING THE
*   COEFFICIENTS
*          C(I,J)   (I=1,...,NCX , J=1,...,NCY)
*   OF A TWO-DIMENSIONAL POLYNOMIAL APPROXIMATION SPLINE  Z = S(X,Y)  IN
*   REPRESENTATION OF NORMALIZED TWO-DIMENSIONAL B-SPLINES  B(I,J)(X,Y).
*
************************************************************************

      PARAMETER (Z0 = 0 , Z1 = 1 , Z2 = 2 , Z10 = 10 , HALF = Z1/Z2)
*
*   COMPUTE AN APPROXIMATION EPS0 TO THE RELATIVE MACHINE PRECISION
*
      EPS0=Z1
   10 EPS0=EPS0/Z10
      IF(Z1+EPS0 .NE. Z1) GO TO 10
      EPS0=Z10*EPS0
*
*   COMPUTE KNOTS BY MEANS OF GIVEN DATA POINTS (IF KNOT = 1 OR 2)
*
      IF (KNOT .EQ. 1) THEN
       CALL DSPKN2(KX,KY,MX,MY,XI(1),XI(NX),YI(1),YI(NY),TX,TY,NERR)
      ELSEIF(KNOT .EQ. 2) THEN
       DO 20 I=1,KX+1
       TX(I)=XI(1)
   20  TX(NCX+I)=XI(NX)
       DO 30 I=KX+2,NCX
   30  TX(I)=HALF*(XI(NX*(I-KX-2)/NCX+1)+XI(NX*I/NCX))
       DO 40 J=1,KY+1
       TY(J)=YI(1)
   40  TY(NCY+J)=YI(NY)
       DO 50 J=KY+2,NCY
   50  TY(J)=HALF*(YI(NY*(J-KY-2)/NCY+1)+YI(NY*J/NCY))
      ENDIF
*
*   COMPUTE MATRIX  A  AND SOLVE LINEAR LEAST SQUARES PROBLEM USING  SVD
*
      DO 60  I=1,NX
      X=XI(I)
      DO 60  J=1,NY
      Y=YI(J)
      DO 60 IC=1,NCX
      DO 60 JC=1,NCY
      A((I-1)*NY+J,(IC-1)*NCY+JC)=
     +            DSPNB2(KX,KY,MX,MY,IC,JC,0,0,X,Y,TX,TY,NERR)
   60 CONTINUE
      CALL DGESVD('O','A',N,NC,A,N,S,W1,1,VT,NC,W2,LW,INFO)
      DO 70 I=1,NX
      DO 70 J=1,NY
   70 W1((I-1)*NY+J)=ZI(I,J)
      CALL DMMPY(NC,N,A(1,1),A(2,1),A(1,2),W1(1),W1(2),W2(1),W2(2))
      DO 80 J=1,NC
      IF(S(J) .GT. EPS0*S(1)) THEN
       W1(J)=W2(J)/S(J)
      ELSE
       W1(J)=Z0
      ENDIF
   80 CONTINUE
      CALL DMMPY(NC,NC,VT(1,1),VT(2,1),VT(1,2),W1(1),W1(2),W2(1),W2(2))
      DO 90 I=1,NCX
      DO 90 J=1,NCY
   90 C(I,J)=W2((I-1)*NCY+J)
      NERR=0

      RETURN
      END



