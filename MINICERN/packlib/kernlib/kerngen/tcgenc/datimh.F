*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMMPW)
#include "mpwgs/datimh.F"
#elif defined(CERNLIB_QMSUN)
#include "sungs/datimh.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/datimh.F"
#else
      SUBROUTINE DATIMH (ND,NT)
C
C CERN PROGLIB# Z007    DATIMH  DUMMY   .VERSION KERNFOR  4.03  821008
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE

      DIMENSION ND(9), NT(9)
      DIMENSION M(8)

      CALL UBLOW (8H29/09/79,M,8)
      CALL UBUNCH           (M,ND,8)
      CALL UBLOW (8H12.00.00,M,8)
      CALL UBUNCH           (M,NT,8)
      RETURN
      END
#endif
