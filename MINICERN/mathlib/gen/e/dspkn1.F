*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:25  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DSPKN1(K,M,A,B,T,NERR)

#include "gen/imp64.inc"
      DIMENSION T(*)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DSPKN1')

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   DSPKN1 COMPUTES M ( M >= 2*K+2 ) KNOTS IN THE INTERVAL  A <= X <= B
*   AS KNOTS FOR NORMALIZED B-SPLINES OF DEGREE  K  ( 0 <= K <= 25 ).
*   THE FIRST (K+1) KNOTS ARE EQUAL TO A, THE LAST (K+1) KNOTS ARE EQUAL
*   TO B, AND THE OTHER KNOTS ARE EQUIDISTANT.
*
*   PARAMETERS:
*
*   K    (INTEGER) DEGREE OF B-SPLINES
*   M    (INTEGER) NUMBER OF KNOTS IN THE INTERVAL  A <= X <= B
*   A    (DOUBLE PRECISION) LEFT  ENDPOINT OF INTERVAL
*   B    (DOUBLE PRECISION) RIGHT ENDPOINT OF INTERVAL
*   T    (DOUBLE RECISION) ARRAY OF AT LEAST ORDER M CONTAINING THE
*         KNOTS, ON EXIT.
*   NERR (INTEGER) ERROR INDICATOR. ON EXIT:
*        = 0: NO ERROR DETECTED
*        = 1: AT LEAST ONE OF THE CONSTANTS K , M IS ILLEGAL
*
*   USAGE:
*
*       THE USER HAS TO PROVIDE K, M, A, AND B.
*       ON RETURN THE ARRAY T CONTAINS THE DESIRED KNOTS IN ASCENDING
*       ORDER.
*
*   ERROR MESSAGES:
*
*   IF ONE OF THE FOLLOWING RELATION IS SATISFIED BY THE CHOSEN INPUT-
*   PARAMETERS THE PROGRAM RETURNS, AND AN ERROR MESSAGE IS PRINTED:
*     K < 0      OR    K > 25    OR
*     M < 2*K+2  OR
*     B <= A .
*
*   ********************************************************************

      NERR=1
      IF(K .LT. 0 .OR. K .GT. 25) THEN
       WRITE(ERRTXT,101) 'K',K
       CALL MTLPRT(NAME,'E210.1',ERRTXT)
      ELSEIF(M .LT. 2*K+2) THEN
       WRITE(ERRTXT,101) 'M',M
       CALL MTLPRT(NAME,'E210.2',ERRTXT)
      ELSEIF(B .LE. A) THEN
       WRITE(ERRTXT,103) 'A',A,'B',B
       CALL MTLPRT(NAME,'E210.7',ERRTXT)
      ELSE
       NERR=0
       CALL DVSET(K+1,A,T(1),T(2))
       CALL DVSET(K+1,B,T(M-K),T(M-K+1))
       D=(B-A)/(M-2*K-1)
       DO 10 I = K+2,M-K-1
   10  T(I)=A+D*(I-K-1)
      ENDIF
      RETURN

  101 FORMAT(1X,A5,' =',I6,'   NOT IN RANGE')
  103 FORMAT(1X,A5,' =',1PD15.8,A7,' =',1PD15.8,'   INCONSISTENT')
      END


