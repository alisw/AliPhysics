*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:52  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RRTEQ4(A,B,C,D,Z,DC,MT)
#if defined(CERNLIB_DOUBLE)
C
      CHARACTER NAME*6
      NAME = 'RRTEQ4'
      CALL MTLPRT(NAME,'C208.0','Not available on this machine')
      RETURN
#endif
      ENTRY RTEQ4(A,B,C,D,Z,DC,MT)
#if defined(CERNLIB_DOUBLE)
      NAME = 'RTEQ4'
      CALL MTLPRT(NAME,'C208.0','Not available on this machine')
      RETURN
#endif
#if !defined(CERNLIB_DOUBLE)
      CALL RRTEQ4(A,B,C,D,Z,DC,MT)
#endif
      END
