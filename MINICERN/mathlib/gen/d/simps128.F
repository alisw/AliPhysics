*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DSIMPS(F,A,B,N2)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME
      PARAMETER( NAME='DSIMPS')
      DSIMPS=0
      CALL MTLPRT(NAME,'D101.0','Not available on this machine')
      RETURN
      END
#endif
