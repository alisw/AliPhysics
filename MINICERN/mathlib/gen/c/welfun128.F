*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:01  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE WELFUN(W,AK2,SN,CN,DN)
C
#include "gen/imp128.inc"
      CHARACTER*(*) NAME
      PARAMETER(NAME='WELFUN')
C
      CALL MTLPRT(NAME,'C320.0','Not available on this machine')
      RETURN
      END
#endif
