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
      FUNCTION DATANI(X)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME
      PARAMETER( NAME ='DATANI')
      DATANI=0
      CALL MTLPRT(NAME1,'C324.0','Not available on this machine')
      RETURN
      END
#endif
