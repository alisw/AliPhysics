*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:51  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DSRTNT(K,N,A,B,C,U1,V1,RES,LRL)
C
#include "gen/imp128.inc"
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DSRTNT')
      CALL MTLPRT(NAME,'B300.0','Not available on this machine')
      RETURN
      END
#endif
