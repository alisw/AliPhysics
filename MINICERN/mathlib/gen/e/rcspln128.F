*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:27  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DCSPLN(N,X,M,Y,NDIM,MODE,A,B,C,D)
#include "gen/imp128.inc"
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DCSPLN/DCSPNT')

      ENTRY DCSPNT(N,X,M,Y,NDIM,MODE,A,B,C,D)
      CALL MTLPRT(NAME,'E211.0','Not available on this machine')
      RETURN
      END
#endif
