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
      SUBROUTINE D107D1(MODE,F,A,B,N,X,W)
#include "gen/imp64.inc"
      CHARACTER*6 NAME(2)
      DATA NAME /'DGQUAD','DGSET'/
#include "d107rd.inc"

#endif
