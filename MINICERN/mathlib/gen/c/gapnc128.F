*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:05  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DGAPNC(A,X)
#include "gen/imp128.inc"
C
      CHARACTER*(*) NAME1, NAME2
      PARAMETER( NAME1='DGAPNC')
      PARAMETER( NAME2='DGAGNC')
      DGAPNC=0
      CALL MTLPRT(NAME1,'C334.0','Not available on this machine')
      RETURN

      ENTRY DGAGNC(A,X)
      DGAGNC=0
      CALL MTLPRT(NAME2,'C334.0','Not available on this machine')
      END
#endif
