*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:03  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DCLAUS(X)
      CHARACTER*6 NAME
      NAME = 'DCLAUS'
      CALL MTLPRT(NAME,'C326',
     +'not available on this machine - see documentation')
      DCLAUS = 0.D0
      RETURN
      END
#endif
