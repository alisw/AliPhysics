*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:15  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
C This will be CGAUSS,IF=DOUBLE and CGAUSS64,IF=-DOUBLE.
      FUNCTION CGAUSS(F,A,B,EPS)
      COMPLEX CGAUSS,F,A,B,AA,BB,U,C1,C2,S8,S16,H
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'CGAUSS')
#include "cgausscod.inc"
   99 CGAUSS=H
      RETURN
      END
