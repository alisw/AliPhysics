*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:02  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DFRSIN(X)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME1, NAME2
      PARAMETER( NAME1='DFRSIN')
      PARAMETER( NAME2='DFRCOS')
      DFRSIN=0
      CALL MTLPRT(NAME1,'C322.0','Not available on this machine')
      RETURN

      ENTRY DFRCOS(X)
      DFRCOS=0
      CALL MTLPRT(NAME2,'C322.0','Not available on this machine')
      RETURN
      END
#endif
