*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:56  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RNPSSN(AMU,N,IERR)
      EXTERNAL RANLUX
 
      SAVE EMU,AMU0,AMAX
 
      DATA AMU0 /-12345.67/, AMAX /88/
      PARAMETER (AMXA = 88)
 
      IERR=0
      IF(AMU .LE. 0) THEN
       IERR=1
       J=0
      ELSEIF(AMU .GT. AMAX) THEN
       CALL RNORMX(R,1,RANLUX)
       J=R*SQRT(AMU)+AMU+0.5
      ELSE
       IF(AMU .NE. AMU0) THEN
        AMU0=AMU
        EMU=EXP(-AMU)
       ENDIF
       P=1
       J=-1
    1  J=J+1
       CALL RANLUX(R,1)
       P=P*R
       IF(P .GT. EMU) GO TO 1
      ENDIF
      N=J
      RETURN
 
      ENTRY RNPSET(AMX)
      AMAX=MIN(AMX,AMXA)
      WRITE(6,'(/7X,''+++++ CERN V136 RNPSSN :  SWITCH TO '',
     1 ''NORMAL APPROXIMATION FOR      AMU > '',F7.2/)') AMAX
      RETURN
      END
