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
      SUBROUTINE DRTEQ4(A,B,C,D,Z,DC,MT)
C
#include "gen/imp128.inc"
C
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DRTEQ4')
      CALL MTLPRT(NAME,'C208.0','Not available on this machine')
      RETURN
      END
#endif
