*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:52  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RRTEQ3(R,S,T,X,D)
#if defined(CERNLIB_DOUBLE)
C
      CHARACTER NAME*6
      NAME = 'RRTEQ3'
      CALL MTLPRT(NAME,'C207.0','Not available on this machine')
      RETURN
#endif
      ENTRY RTEQ3(R,S,T,X,D)
#if defined(CERNLIB_DOUBLE)
      NAME = 'RTEQ3'
      CALL MTLPRT(NAME,'C207.0','Not available on this machine')
#endif
#if !defined(CERNLIB_DOUBLE)
      CALL RRTEQ3(R,S,T,X,D)
#endif
      END
