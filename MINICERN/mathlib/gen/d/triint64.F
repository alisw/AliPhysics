*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:14  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DTRINT(F,NSD,NPT,EPS,X1,Y1,X2,Y2,X3,Y3)
C     INTEGRATION OVER A TRIANGLE USING A 7-, 25- OR 64-POINT FORMULA,
C     WITH OR WITHOUT SUBDIVISION OF THE TRIANGLE.

#include "gen/imp64.inc"
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DTRINT')
#include "triintcod.inc"

   99 DTRINT=H
      RETURN
  101 FORMAT('INCORRECT NUMBER OF POINTS =',I5)
      END
#endif
