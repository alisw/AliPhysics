*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:49  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RDJMNB(AJ,AM,AN,BETA)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'RDJMNB')
      ENTRY DJMNB(AJ,AM,AN,BETA)
      CALL MTLPRT(NAME,'U501.0','Not available on this machine')
      DJMNB=0
      RETURN
      END
#endif
