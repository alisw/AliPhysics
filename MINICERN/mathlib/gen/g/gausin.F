*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:42  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION GAUSIN(P)

C     Computes a "Normal Deviate"
C     Based on G.W. Hill & A.W. Davis, Algorithm 442 Normal Deviate
C     Collected Algorithms from CACM

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'GAUSIN')
      PARAMETER (C = 2.50662 827)

      IF(P .LE. 0 .OR. P .GE. 1) THEN
       H=0
       WRITE(ERRTXT,101) P
       CALL MTLPRT(NAME,'G105.1',ERRTXT)
      ELSEIF(P .EQ. 0.5) THEN
       H=0
      ELSE
       X=P
       IF(P .GT. 0.5) X=1-P
       X=SQRT(-2*LOG(X))
       X=X-((7.45551*X+450.636)*X+1271.059)/
     1     (((X+110.4212)*X+750.365)*X+500.756)
       IF(P .LT. 0.5) X=-X
       H=C*(P-FREQ(X))*EXP(0.5*X**2)+X
      ENDIF
      GAUSIN=H
      RETURN
  101 FORMAT('ARGUMENT P =',1P,E15.5,' NOT IN RANGE')
      END
#endif
