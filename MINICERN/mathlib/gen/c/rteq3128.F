*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:52  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DRTEQ3(R,S,T,X,D)
C
#include "gen/imp128.inc"
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DRTEQ3')
      CALL MTLPRT(NAME,'C207.0','Not available on this machine')
      RETURN
      END
#endif
