*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:08  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DSTRH0(X)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME1,NAME2
      PARAMETER( NAME1='DSTRH0')
      PARAMETER( NAME2='DSTRH1')
      DSTRH0=0
      CALL MTLPRT(NAME,'C342.0','Not available on this machine')
      RETURN
      ENTRY DSTRH1(X)
      DSTRH1=0
      CALL MTLPRT(NAME,'C342.0','Not available on this machine')
      RETURN
      END
#endif
