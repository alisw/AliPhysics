*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:51  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RZEROX(A0,B0,EPS,MAXF,F,MODE)
C
      CHARACTER NAME*6
      NAME = 'RZEROX'
      CALL MTLPRT(NAME,'C200.0','Not available on this machine')
      RETURN
      ENTRY ZEROX(A0,B0,EPS,MAXF,F,MODE)
      NAME = 'ZEROX'
      CALL MTLPRT(NAME,'C200.0','Not available on this machine')
      END
#endif
