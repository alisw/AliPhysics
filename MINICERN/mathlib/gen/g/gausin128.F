*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:42  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DGAUSN(P)
#include "gen/imp128.inc"
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DGAUSN')
      DGAUSN=0
      CALL MTLPRT(NAME,'G105.0','Not available on this machine')
      DGAUSN=0
      RETURN
      END
#endif
