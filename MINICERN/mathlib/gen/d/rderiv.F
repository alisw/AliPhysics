*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RDERIV(F,X,DELTA,DFDX,RERR)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'RDERIV')
      ENTRY DERIV(F,X,DELTA,DFDX,RERR)
      CALL MTLPRT(NAME,'C401.0','Not available on this machine')
      RETURN
      END
#endif
