*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DSPVD1(F,K,M,T,C,NERR)

#include "gen/imp64.inc"
      DIMENSION T(*),C(*)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DSPVD1')

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   DSPVD1 COMPUTES THE COEFFICIENTS  C(1),...,C(M-K-1)  OF A POLYNOMIAL
*   VARIATION DIMINISHING SPLINE APPROXIMATION  S(X)  IN B-SPLINE
*   REPRESENTATION
*
*          S(X) =  SUMME(I=1,...,M-K-1)  C(I) * B(I,K)(X)
*
*   TO A USER SUPPLIED FUNCTION  Y=F(X).
*
*   THE FUNCTIONS  B(I,K)(X)  ARE NORMALIZED B-SPLINES OF DEGREE  K
*   (0<= K <= 25)  WITH INDEX  I  (1 <= I <= M-K-1) OVER A SET OF
*   SPLINE-KNOTS
*              T(1),T(2), ... ,T(M)    ( M >= 2*K+2 )
*   (KNOTS IN ASCENDING ORDER, WITH MULTIPLICITIES NOT GREATER
*   THAN  K+1).
*   THE FUNCTION VALUE OF THE NORMALIZED B-SPLINE  B(I,K)(X)  IS
*   IDENTICALLY ZERO OUTSIDE THE INTERVAL  T(I) <= X < T(I+K+1).
*   THE NORMALIZATION OF  N(X) = B(I,K)(X)  IS SUCH THAT THE INTGRAL OF
*   N(X)  OVER THE WHOLE X-RANGE EQUALS
*                  ( T(I+K+1) - T(I) ) / (K+1)  .
*
*   PARAMETERS:
*
*   F     (DOUBLE PRECISION) USER SUPPLIED FUNCTION  F(X)  FOR WHICH
*         THE CORRESPONDING SPLINE APPROXIMATION HAS TO BE COMPUTED.
*   K     (INTEGER) DEGREE (= ORDER - 1) OF B-SPLINES.
*   M     (INTEGER) NUMBER OF KNOTS FOR THE B-SPLINES.
*   T     (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER M CONTAINING THE
*         KNOTS, ON ENTRY.
*   C     (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER  M-K-1. ON EXIT
*         C(1),...,C(M-K-1)  CONTAIN THE COEFFICIENTS OF THE B-SPLINE
*         REPRESENTATION OF  S(X).
*   NERR  (INTEGER) ERROR INDICATOR. ON EXIT:
*         = 0: NO ERROR DETECTED
*         = 1: AT LEAST ONE OF THE CONSTANTS  K , M  IS ILLEGAL
*
*   ERROR MESSAGES:
*
*   IF ONE OF THE FOLLOWING RELATION IS SATISFIED BY THE CHOSEN INPUT-
*   PARAMETERS THE PROGRAM RETURNS, AND AN ERROR MESSAGE IS PRINTED:
*     K < 1      OR    K > 25    OR
*     M < 2*K+2  .
*
************************************************************************

      PARAMETER (Z1 = 1)

      NERR=1
      IF(K .LT. 1 .OR. K .GT. 25) THEN
       WRITE(ERRTXT,101) 'K',K
       CALL MTLPRT(NAME,'E210.1',ERRTXT)
      ELSEIF(M .LT. 2*K+2) THEN
       WRITE(ERRTXT,101) 'M',M
       CALL MTLPRT(NAME,'E210.2',ERRTXT)
      ELSE
       NERR=0
       R=Z1/K
       DO 10 I = 1,M-K-1
   10  C(I)=F(R*DVSUM(K,T(I+1),T(I+2)))
      ENDIF

      RETURN

  101 FORMAT(1X,A5,' =',I6,'   NOT IN RANGE')
      END



