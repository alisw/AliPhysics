*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DDJMNB(AJ,AM,AN,BETA)
#include "gen/imp128.inc"
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DDJMNB')
      DDJMNB=0
      CALL MTLPRT(NAME,'U501.0','Not available on this machine')
      RETURN
      END
#endif
