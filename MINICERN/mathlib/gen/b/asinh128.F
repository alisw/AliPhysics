*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:50  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DASINH(X)
C
#include "gen/imp128.inc"
      CHARACTER*(*) NAME
      PARAMETER( NAME='DASINH')
      DASINH=0
      CALL MTLPRT(NAME,'B102.0','Not available on this machine')
      RETURN
      END
#endif
