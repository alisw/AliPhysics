*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:40  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DLHOIN(A,MA,M,N,MAXV,V,NV,NVEC,EPS,IOUT,W,IW)
#include "gen/imp128.inc"
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DLHOIN')
      CALL MTLPRT(NAME,'F500.0','Not available on this machine')
      RETURN
      END
#endif
