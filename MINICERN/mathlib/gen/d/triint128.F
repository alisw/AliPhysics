*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:14  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DTRINT(F,NSD,NPT,EPS,X1,Y1,X2,Y2,X3,Y3)
C
#include "gen/imp128.inc"
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DTRINT')
      DTRINT=0
      CALL MTLPRT(NAME,'D105.0','Not available on this machine')
      RETURN
      END
#endif
