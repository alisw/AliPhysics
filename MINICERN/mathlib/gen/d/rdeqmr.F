*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:17  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RDEQMR(N,XA,XZ,Y,H0,EPS,SUB,W)
C
      CHARACTER NAME*6
      NAME = 'RDEQMR'
      CALL MTLPRT(NAME,'C202.0','Not available on this machine')
      RETURN

      ENTRY DEQMR(N,XA,XZ,Y,H0,EPS,SUB,W)
      NAME = 'DEQMR'
      CALL MTLPRT(NAME,'C202.0','Not available on this machine')
      RETURN
      END
#endif
