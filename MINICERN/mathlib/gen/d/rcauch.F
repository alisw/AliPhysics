*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:14  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RCAUCH(F,A,B,S,EPS)
C
      CHARACTER NAME*6
      NAME = 'RCAUCH'
      CALL MTLPRT(NAME,'D104.0','Not available on this machine')
      RCAUCH=0
      RETURN

      ENTRY CAUCHY(F,A,B,S,EPS)
      NAME = 'CAUCHY'
      CALL MTLPRT(NAME,'D104.0','Not available on this machine')
      CAUCHY=0
      END
#endif
