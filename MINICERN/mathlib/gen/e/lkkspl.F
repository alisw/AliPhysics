*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION LKKSPL(X,T,M)

#include "gen/imp64.inc"
      DIMENSION T(*)

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   FOR GIVEN  X  AND A NON-DECREASING SEQUENCE  T(1) , T(2) ,..., T(M)
*   THE FUNCTION  LKKSPL  COMPUTES AND RETURNS THE UNIQUELY DETERMINED
*   INDEX  KK  SUCH THAT THE RELATION
*              T(KK-1) <=  X  < T(KK)
*   HOLDS.
*   IT RETURNS  KK = 1 ,    IF   X <  T(1) ,   AND
*               KK = M+1 ,  IF   X >= T(M) .
*
*   (THE SEARCH ALGORITHM IS OF REGULA FALSI - TYPE).
*
************************************************************************

      IF(X .GE. T(M)) THEN
       LKKSPL=M+1
      ELSEIF(X .LT. T(1)) THEN
       LKKSPL=1
      ELSE
       M1=1
       M2=M
       T1=T(M1)
       T2=T(M2)
   10  M3=M1+(M2-M1)*(X-T1)/(T2-T1)
       M4=M3+1
       T3=T(M3)
       T4=T(M4)
       IF((T3 .GT. X .OR. X .GE. T4) .AND.
     +    (T3 .NE. X .OR. X .NE. T4)) THEN
        IF(X .GE. T3) THEN
         M1=M4
         T1=T4
        ELSE
         M2=M3
         T2=T3
        ENDIF
        GO TO 10
       ENDIF
   20  LKKSPL=M4
       IF(X .GE. T(M4)) THEN
        M4=M4+1
        GO TO 20
       ENDIF
      ENDIF

      RETURN
      END



