*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DSIMPS(F,A,B,N2)
#include "gen/imp64.inc"

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DSIMPS')
#include "simpscod.inc"

      DSIMPS=H
      RETURN
  101 FORMAT('NON-POSITIVE OR EVEN NUMBER OF FUNCTION VALUES =',I6)
      END
#endif
