*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:25  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DSPAP1(K,M,N,XI,YI,KNOT,T,C,W,NW,NERR)

#include "gen/imp64.inc"
      DIMENSION XI(*),YI(*),T(*),W(*),C(*)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DSPAP1')

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   DSPAP1 COMPUTES THE COEFFICIENTS  C(1),...,C(NC)  OF A POLYNOMIAL
*   APPROXIMATION SPLINE  S(X)  IN B-SPLINE REPRESENTATION
*
*          S(X) =  SUMME(I=1,...,NC)  C(I) * B(I,K)(X)  ,  NC=M-K-1
*
*   TO A USER SUPPLIED DATA SET
*
*          (XI(J),YI(J)) ,   J = 1,2,...,N  ,  N >= M-K-1 >= K+1
*
*   OF A FUNCTION  Y=F(X) , I.E.
*
*           S(XI(J)) = YI(J) , J = 1,2,...,N .
*
*   THE FUNCTIONS  B(I,K)(X)  ARE NORMALIZED B-SPLINES OF DEGREE  K
*   (0<= K <= 25) WITH INDEX I (1 <= I <= N) OVER A SET OF SPLINE-KNOTS
*              T(1),T(2), ... ,T(M)    ( M <= N+K+1 )
*   (KNOTS IN ASCENDING ORDER, WITH MULTIPLICITIES NOT GREATER
*   THAN  K+1).
*   FOR FURTHER DETAILS TO THE ONE-DIMENSIONAL NORMALIZED B-SPLINES SEE
*   THE COMMENTS TO  DSPNB1.
*
*   PARAMETERS:
*
*   N     (INTEGER) NUMBER OF APPROXIMATION POINTS .
*   M     (INTEGER) NUMBER OF KNOTS.
*   K     (INTEGER) DEGREE OF B-SPLINES.
*   XI    (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER N .
*         XI MUST CONTAIN THE APPROXIMATION POINTS IN ASCENDING ORDER,
*         ON ENTRY.
*   YI    (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER N CONTAINING
*         THE FUNCTION VALUES YI(J), J=1,...,N,  ON ENTRY.
*   KNOT  (INTEGER) PARAMETER FOR STEERING THE CHOICE OF KNOTS.
*         ON ENTRY:
*         = 1 : KNOTS ARE COMPUTED BY  DSPAP1  IN THE FOLLOWING WAY:
*               T(J) = XI(1) ,                J = 1,...,K+1
*               T(J) = XI(1)+(J-K-1)*(XI(N)-XI(1))/(NC-K) ,
*                                             J = K+2,...,NC
*               T(NC+J) = XI(N) ,             J = 1,...,K+1
*         = 2 : KNOTS ARE COMPUTED BY  DSPAP1  IN THE FOLLOWING WAY:
*               T(J) = XI(1) ,                J = 1,...,K+1
*               T(J) = (XI(J-K-1)+XI(J))/2 ,  J = K+2,...,N
*               T(N+J) = XI(N) ,              J = 1,...,K+1
*         OTHERWISE KNOTS ARE USER SUPPLIED. RECOMMENDED CHOICE :
*               T(1) <= ... <= T(K+1) <= XI(1)
*               XI(1) < T(K+2) < ... < T(NC) < XI(N)
*               XI(N) <= T(NC+1) <= ... <= T(NC+K+1)
*   T     (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER M .
*         IF THE INPUT VALUE OF THE PARAMETER KNOT IS  1 OR 2  THE
*         KNOTS ARE COMPUTED BY  DSPAP1  AND THEY ARE GIVEN IN THE
*         ARRAY T, ON EXIT.
*         IN THE OTHER CASES THE ARRAY  T  MUST CONTAIN THE USER
*         SUPPLIED KNOTS, ON ENTRY.
*   W     (DOUBLE PECISION) WORKING ARRAY OF AT LEAST ORDER  NW.
*   NW    (INTEGER) ORDER OF WORKING ARRAY  W .
*         NW >= N*(NC+5)+NC*(NC+1) ,  WITH   NC=M-K-1 .
*         FOR GOOD PERFORMANCE,  NW  SHOULD GENERALLY BE LARGER.
*   C     (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER  N . ON EXIT
*         C(1),...,C(NC)  CONTAIN THE COEFFICIENTS OF THE B-SPLINE
*         REPRESENTATION OF  S(X).
*   NERR  (INTEGER) ERROR INDICATOR. ON EXIT:
*         = 0: NO ERROR DETECTED
*         = 1: AT LEAST ONE OF THE CONSTANTS  K , M , N  IS ILLEGAL
*
*   ERROR MESSAGES:
*
*   IF ONE OF THE FOLLOWING RELATION IS SATISFIED BY THE CHOSEN INPUT-
*   PARAMETERS THE PROGRAM RETURNS, AND AN ERROR MESSAGE IS PRINTED:
*     K  < 0      OR    K > 25    OR
*     M  < 2*K+2  OR
*     N  < M-K-1 .
*
************************************************************************

      NERR=1
      IF(K .LT. 0 .OR. K .GT. 25) THEN
       WRITE(ERRTXT,101) 'K',K
       CALL MTLPRT(NAME,'E210.1',ERRTXT)
      ELSEIF(M .LT. 2*K+2) THEN
       WRITE(ERRTXT,101) 'M',M
       CALL MTLPRT(NAME,'E210.2',ERRTXT)
      ELSEIF(N .LT. M-K-1) THEN
       WRITE(ERRTXT,101) 'N',N
       CALL MTLPRT(NAME,'E210.4',ERRTXT)
      ELSE

       NC=M-K-1
       M1=1
       M2=M1+N*NC
       M3=M2+NC
       M4=M3+NC*NC
       LW=NW-M4+1

       CALL SPLAS1(N,NC,M,K,XI,YI,KNOT,T,W,W(M2),W(M3),W(M4),LW,C,NERR)

      ENDIF
      RETURN

  101 FORMAT(1X,A5,' =',I6,'   NOT IN RANGE')
      END



