*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:53  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DFREQ(X)
C
#include "gen/imp128.inc"
      CHARACTER*(*) NAME
      PARAMETER(NAME='DFREQ')
C
      DFREQ=0
      CALL MTLPRT(NAME,'C301.0','Not available on this machine')
      RETURN
      END
#endif
