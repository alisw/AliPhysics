*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:07  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DDAWSN(X)
C
#include "gen/imp128.inc"
C
      CHARACTER*(*) NAME
      PARAMETER( NAME='DDAWSN')
      DDAWSN=0
      CALL MTLPRT(NAME,'C339.0','Not available on this machine')
      RETURN
      END
#endif
