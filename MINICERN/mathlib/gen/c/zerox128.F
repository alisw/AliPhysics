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
      SUBROUTINE DZEROX(A0,B0,EPS,MAXF,F,MODE)
C
#include "gen/imp128.inc"
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DZEROX')
      CALL MTLPRT(NAME,'C200.0','Not available on this machine')
      RETURN
      END
#endif
