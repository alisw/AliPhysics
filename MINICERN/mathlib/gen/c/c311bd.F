*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE C311BD
#include "legbl.inc"
      DATA IBD/0/
      IF(IBD.NE.0) RETURN
      IBD=1
      PISR=1.7724 53850 90552
      PI=3.141592653589793
      U=(0.0,1.0)
      RETURN
      END
