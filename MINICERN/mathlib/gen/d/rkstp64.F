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
      SUBROUTINE DRKSTP(N,H,X,Y,SUB,W)
#include "gen/imp64.inc"
#include "rkstpcod.inc"
#endif
