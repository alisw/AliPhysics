*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
C This will be SIMPS,IF=DOUBLE and SIMPS64,IF=-DOUBLE
      FUNCTION SIMPS(F,A,B,N2)

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'SIMPS')
#include "simpscod.inc"
      SIMPS=H
      RETURN
  101 FORMAT('NON-POSITIVE OR EVEN NUMBER OF FUNCTION VALUES =',I6)
      END
