*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:52  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RSNLEQ(N,X,F,FTOL,XTOL,MAXF,IPRT,INFO,SUB,W)
#if defined(CERNLIB_DOUBLE)
C
      CHARACTER NAME*6
      NAME = 'RSNLEQ'
      CALL MTLPRT(NAME,'C201.0','Not available on this machine')
      RETURN
#endif
      ENTRY SNLEQ(N,X,F,FTOL,XTOL,MAXF,IPRT,INFO,SUB,W)
#if defined(CERNLIB_DOUBLE)
      NAME = 'SNLEQ'
      CALL MTLPRT(NAME,'C201.0','Not available on this machine')
#endif
#if !defined(CERNLIB_DOUBLE)
      CALL RSNLEQ(N,X,F,FTOL,XTOL,MAXF,IPRT,INFO,SUB,W)
#endif
      END
