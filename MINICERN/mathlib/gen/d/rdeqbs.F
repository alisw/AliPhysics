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
      SUBROUTINE RDEQBS(N,XA,XZ,Y,H0,EPS,SUB,W)
C
      CHARACTER NAME*6
      NAME = 'RDEQBS'
      CALL MTLPRT(NAME,'C201.0','Not available on this machine')
      RETURN

      ENTRY DEQBS(N,XA,XZ,Y,H0,EPS,SUB,W)
      NAME = 'DEQBS'
      CALL MTLPRT(NAME,'C201.0','Not available on this machine')
      END
#endif
