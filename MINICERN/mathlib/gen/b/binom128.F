*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:49  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DBINOM(X)
C
#include "gen/imp128.inc"
      CHARACTER*(*) NAME
      PARAMETER( NAME='DBINOM')
      DBINOM=0
      CALL MTLPRT(NAME,'B100.0','Not available on this machine')
      RETURN
      END
#endif
