*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:27  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RCSPLN(N,X,M,Y,NDIM,MODE,A,B,C,D)
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'RCSPLN/RCSPNT')

      ENTRY RCSPNT(N,X,M,Y,NDIM,MODE,A,B,C,D)
      CALL MTLPRT(NAME,'E211.0','Not available on this machine')
      RETURN
      END
#endif
