*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE SPLAS1(N,NC,M,K,XI,YI,KNOT,T,A,S,VT,W,LW,C,NERR)

#include "gen/imp64.inc"
      DIMENSION XI(*),YI(*),T(*),A(N,*),S(*),VT(NC,*),W(*),C(*)

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   THE SUBROUTINE  SPLAS1  IS USED BY  DSPAP1  FOR COMPUTING THE
*   COEFFICIENTS  C(1),...,C(NC)  OF A POLYNOMIAL APPROXIMATION SPLINE
*   S(X)  IN B-SPLINE REPRESENTATION
*
************************************************************************

      PARAMETER (Z0 = 0 , Z1 = 1 , Z2 = 2 , Z10 = 10 , HALF = Z1/Z2)
*
*   COMPUTE AN APPROXIMATION EPS0 TO THE RELATIVE MACHINE PRECISION
*
      EPS0=Z1
   10 EPS0=EPS0/Z10
      IF (Z1+EPS0 .NE. Z1) GO TO 10
      EPS0=Z10*EPS0
*
*   COMPUTE KNOTS BY MEANS OF GIVEN DATA POINTS (IF KNOT = 1 OR 2)
*
      IF (KNOT .EQ. 1) THEN
       CALL DSPKN1(K,M,XI(1),XI(N),T,NERR)
      ELSEIF (KNOT .EQ. 2) THEN
       DO 20 I=1,K+1
       T(I)=XI(1)
   20  T(NC+I)=XI(N)
       DO 30 I=K+2,NC
   30  T(I)=HALF*(XI(N*(I-K-2)/NC+1)+XI(N*I/NC))
      ENDIF
*
*   COMPUTE MATRIX  A  AND SOLVE LINEAR LEAST SQUARES PROBLEM USING  SVD
*
      DO 40 I=1,N
      DO 40 J=1,NC
   40 A(I,J)=DSPNB1(K,M,J,0,XI(I),T,NERR)
      CALL DGESVD('O','A',N,NC,A,N,S,W,1,VT,NC,W,LW,INFO)
      CALL DMMPY(NC,N,A(1,1),A(2,1),A(1,2),YI(1),YI(2),W(1),W(2))
      DO 50 J=1,NC
      IF (S(J) .GT. EPS0*S(1)) THEN
       W(J)=W(J)/S(J)
      ELSE
       W(J)=Z0
      ENDIF
   50 CONTINUE
      CALL DMMPY(NC,NC,VT(1,1),VT(2,1),VT(1,2),W(1),W(2),C(1),C(2))
      NERR=0

      RETURN
      END



