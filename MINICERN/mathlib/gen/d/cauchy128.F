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
      FUNCTION DCAUCH(F,A,B,S,EPS)
C
#include "gen/imp128.inc"
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DCAUCH')
      DCAUCH=0
      CALL MTLPRT(NAME,'D104.0','Not available on this machine')
      RETURN
      END
#endif
