*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:12  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DTHETA( K, X, Q )
C
#include "gen/imp128.inc"
C
      DTHETA=0
      CALL MTLPRT('DTHETA','C349.0','Not available on this machine')
      RETURN
      END
#endif
