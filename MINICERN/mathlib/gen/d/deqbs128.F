*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:17  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DDEQBS(N,XA,XZ,Y,H0,EPS,SUB,W)
#include "gen/imp128.inc"
      CHARACTER*(*) NAME
      PARAMETER( NAME='DDEQBS')
      CALL MTLPRT(NAME,'D201.0','Not available on this machine')
      RETURN
      END
#endif
