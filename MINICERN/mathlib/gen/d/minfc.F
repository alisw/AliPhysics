*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:20  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
C This will be RMINFC,IF=DOUBLE and MINFC64 (RMINFC), IF=-DOUBLE
      SUBROUTINE RMINFC(F,A,B,EPS,DELTA,X,Y,LLM)
      LOGICAL LLM,LLT,LGE
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT

      PARAMETER (NAME = 'RMINFC')
#include "minfccod.inc"
