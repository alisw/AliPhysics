*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:20  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DMINFC(F,A,B,EPS,DELTA,X,Y,LLM)
#include "gen/imp64.inc"
      LOGICAL LLM,LLT,LGE
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT

      PARAMETER (NAME = 'DMINFC')
#include "minfccod.inc"
#endif
