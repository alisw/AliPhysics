*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:25  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DSPIN1(K,N,XI,YI,KNOT,T,C,W,IW,NERR)

#include "gen/imp64.inc"
      DIMENSION XI(*),YI(*),T(*),W(*),IW(*),C(*)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DSPIN1')

************************************************************************
*   NORBAS, VERSION: 10.02.1993
************************************************************************
*
*   DSPIN1 COMPUTES THE COEFFICIENTS  C(1),...,C(N)  OF A POLYNOMIAL
*   INTERPOLATION SPLINE  S(X)  IN B-SPLINE REPRESENTATION
*
*          S(X) =  SUMME(I=1,...,N)  C(I) * B(I,K)(X)
*
*   TO A USER SUPPLIED DATA SET
*
*          (XI(J),YI(J)) ,   J = 1,2,...,N  ,  N >= K+1
*
*   OF A FUNCTION  Y=F(X) , I.E.
*
*           S(XI(J)) = YI(J) , J = 1,2,...,N .
*
*   THE FUNCTIONS  B(I,K)(X)  ARE NORMALIZED B-SPLINES OF DEGREE  K
*   ( K <= N-1  AND  0<= K <= 25)  WITH INDEX  I  (1 <= I <= N) OVER A
*   SET OF SPLINE-KNOTS
*              T(1),T(2), ... ,T(M)    ( M = N+K+1 )
*   (KNOTS IN ASCENDING ORDER, WITH MULTIPLICITIES NOT GREATER
*   THAN  K+1).
*   FOR FURTHER DETAILS TO THE ONE-DIMENSIONAL NORMALIZED B-SPLINES SEE
*   THE COMMENTS TO  DSPNB1.
*
*   PARAMETERS:
*
*   N     (INTEGER) NUMBER OF INTERPOLATION POINTS .
*   K     (INTEGER) DEGREE (= ORDER - 1) OF B-SPLINES.
*   XI    (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER N .
*         XI MUST CONTAIN THE INTERPOLATION POINTS IN ASCENDING ORDER,
*         ON ENTRY.
*   YI    (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER N CONTAINING
*         THE FUNCTION VALUES YI(J), J=1,...,N,  ON ENTRY.
*   KNOT  (INTEGER) PARAMETER FOR STEERING THE CHOICE OF KNOTS.
*         ON ENTRY:
*         = 1 : KNOTS ARE COMPUTED BY  DSPIN1  IN THE FOLLOWING WAY:
*               T(J) = XI(1) ,                J = 1,...,K+1
*               T(J) = XI(1)+(J-K-1)*(XI(N)-XI(1)) ,
*                                             J = K+2,...,N
*               T(N+J) = XI(N) ,              J = 1,...,K+1
*         = 2 : KNOTS ARE COMPUTED BY  DSPIN1  IN THE FOLLOWING WAY:
*               T(J) = XI(1) ,                J = 1,...,K+1
*               T(J) = (XI(J-K-1)+XI(J))/2 ,  J = K+2,...,N
*               T(N+J) = XI(N) ,              J = 1,...,K+1
*         OTHERWISE KNOTS ARE USER SUPPLIED. RECOMMENDED CHOICE :
*               T(1) <= ... <= T(K+1) <= XI(1)
*               XI(1) < T(K+2) < ... < T(N) < XI(N)
*               XI(N) <= T(N+1) <= ... <= T(N+K+1)
*   T     (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER M .
*         IF THE INPUT VALUE OF THE PARAMETER KNOT IS  1 OR 2  THE
*         KNOTS ARE COMPUTED BY  DSPIN1  AND THEY ARE GIVEN IN THE
*         ARRAY T, ON EXIT.
*         IN THE OTHER CASES THE ARRAY  T  MUST CONTAIN THE USER
*         SUPPLIED KNOTS, ON ENTRY.
*   W     (DOUBLE PECISION) WORKING ARRAY OF AT LEAST ORDER
*         (3*K+1)*N .
*   IW    (INTEGER) WORKING ARRAY OF AT LEAST ORDER  N .
*   C     (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER  N . ON EXIT
*         C(1),...,C(N)  CONTAIN THE COEFFICIENTS OF THE B-SPLINE
*         REPRESENTATION OF  S(X).
*   NERR  (INTEGER) ERROR INDICATOR. ON EXIT:
*         = 0: NO ERROR DETECTED
*         = 1: AT LEAST ONE OF THE CONSTANTS  K , N  IS ILLEGAL
*         = 2: THE LAPACK ROUTINES  DGBTRF , DGBTRS  COULD NOT SOLVE
*              THE LINEAR SYSTEM OF EQUATIONS WITH BAND-MATRIX FOR
*              COMPUTING C(1),...,C(N) . IT INDICATES THAT A SOLUTION
*              OF THE INTERPOLATION PROBLEM DOES NOT EXIST.
*              (ESPECIALLY, THE EXISTENCE OF A SOLUTION DEPENDS ON THE
*              SET OF KNOTS!)
*
*   ERROR MESSAGES:
*
*   IF ONE OF THE FOLLOWING RELATION IS SATISFIED BY THE CHOSEN INPUT-
*   PARAMETERS THE PROGRAM RETURNS, AND AN ERROR MESSAGE IS PRINTED:
*     K  < 0      OR    K > 25    OR
*     N < K+1 .
*
************************************************************************

      PARAMETER (Z1 = 1 , Z2 = 2 , HALF = Z1/Z2)

      NERR=1
      IF(K .LT. 0 .OR. K .GT. 25) THEN
       WRITE(ERRTXT,101) 'K',K
       CALL MTLPRT(NAME,'E210.1',ERRTXT)
      ELSEIF(N .LT. K+1) THEN
       WRITE(ERRTXT,101) 'N',N
       CALL MTLPRT(NAME,'E210.4',ERRTXT)
      ELSE

       M=N+K+1
       L=3*K+1
*
*   COMPUTE KNOTS FROM INTERPOLATION POINTS (IF KNOT = 1 OR 2)
*
       IF (KNOT .EQ. 1) THEN
        CALL DSPKN1(K,M,XI(1),XI(N),T,NERR)
       ELSEIF (KNOT .EQ. 2) THEN
        DO 10 I=1,K+1
        T(I)=XI(1)
   10   T(N+I)=XI(N)
        DO 20 I=K+2,N
   20   T(I)=HALF*(XI(I-K-1)+XI(I))
       ENDIF
*
*   COMPUTE BAND-MATRIX  W
*
       DO 50 I=K+1,3*K+1
       DO 50 J=1,N
       IJ=I+J-2*K-1
       IF (1 .LE. IJ .AND. IJ .LE. N)
     +    W((J-1)*L+I)=DSPNB1(K,M,J,0,XI(IJ),T,NERR)
   50 CONTINUE
*
*   SOLVE SYSTEM OF EQUATIONS FOR COMPUTING  C
*
       DO 60 J=1,N
   60  IW(J)=J
       NERR=2
       CALL DGBTRF(N,N,K,K,W,L,IW,INFO)
       IF(INFO .NE. 0) RETURN
       CALL DVCPY(N,YI(1),YI(2),C(1),C(2))
       CALL DGBTRS('N',N,K,K,1,W,L,IW,C,N,INFO)
       IF(INFO .NE. 0) RETURN
       NERR=0
      ENDIF
      RETURN

  101 FORMAT(1X,A5,' =',I6,'   NOT IN RANGE')
      END



