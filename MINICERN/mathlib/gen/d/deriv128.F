*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:24  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DDERIV(F,X,DELTA,DFDX,RERR)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME
      PARAMETER( NAME='DDERIV')
      CALL MTLPRT(NAME,'D401.0','Not available on this machine')
      RETURN
      END
#endif
