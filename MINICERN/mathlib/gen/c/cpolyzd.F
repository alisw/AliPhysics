*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:53  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE WPOLYZ(C,N,MAXIT,Z0,H)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'WPOLYZ')
#endif
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE CPOLYZ(C,N,MAXIT,Z0,H)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'CPOLYZ')
#endif
      CALL MTLPRT(NAME,'C209.0','Not available on this machine')
      RETURN
      END
