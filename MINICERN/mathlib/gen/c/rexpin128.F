*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:06  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DEXPIN(X)
#include "gen/imp128.inc"
C
      CHARACTER*(*) NAME1, NAME2
      PARAMETER( NAME1='DEXPIN')
      PARAMETER( NAME2='DEXPIE')
      DEXPIN=0
      CALL MTLPRT(NAME1,'C337.0','Not available on this machine')
      RETURN

      ENTRY DEXPIE(X)
      DEXPIE=0
      CALL MTLPRT(NAME2,'C337.0','Not available on this machine')
      END
#endif
