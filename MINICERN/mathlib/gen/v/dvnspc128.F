*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DVNSPC(R,RHO,D)
C
#include "gen/imp128.inc"
      DVNSPC=0
      CALL MTLPRT('DVNSPC','V700.0','Not available on this machine')
      RETURN
      END
#endif
