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
      FUNCTION DERF(X)
C
#include "gen/imp128.inc"
      CHARACTER*(*) NAME1,NAME2
      PARAMETER(NAME1='DERF')
      PARAMETER(NAME2='DERFC')
      DERF=0
      CALL MTLPRT(NAME,'C300.0','Not available on this machine')
C
      ENTRY DERFC(X)
      DERFC=0
C
      CALL MTLPRT(NAME,'C300.0','Not available on this machine')
      RETURN
      END
#endif
