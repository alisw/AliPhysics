*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:16  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION WGAUSS(F,A,B,EPS)
#include "gen/imp64.inc"
#include "gen/defc64.inc"
     +     WGAUSS,F,A,B,AA,BB,U,C1,C2,S8,S16,H
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'WGAUSS')
#include "cgausscod.inc"

   99 WGAUSS=H
      RETURN
      END
#endif
