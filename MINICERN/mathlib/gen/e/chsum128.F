*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:29  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DCHSUM(MODE,C,N,X)
#include "gen/imp128.inc"
      DIMENSION C(0:N)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DCHSUM')
      CALL MTLPRT(NAME,'E407.0','Not available on this machine')
      DCHSUM=0
      RETURN
      END
#endif
