*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:17  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RRKNYS(N,H,X,Y,YP,SUB,W)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'RRKNYS')
      ENTRY RKNYS(N,H,X,Y,YP,SUB,W)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DRKNYS(N,H,X,Y,YP,SUB,W)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DRKNYS')
#endif
      CALL MTLPRT(NAME,'D203.0','Not available on this machine')
      RETURN
      END
