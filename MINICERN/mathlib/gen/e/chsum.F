*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:29  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION CHSUM(MODE,C,N,X)
      DIMENSION C(0:N)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'CHSUM')
      CALL MTLPRT(NAME,'E407.0','Not available on this machine')
      CHSUM=0
      RETURN
      END
#endif
