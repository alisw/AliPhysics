#include "isajet/pilot.h"
      SUBROUTINE FORTOP
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-     add to force list forced decays for all heavy q particles
C-     if there was a request to force a heavy q decay
C-     Zero IFORCE after use
C-
C-   Created  15-DEC-1989   Serban D. Protopopescu
C-
C    Ver 7.30: Decay top quark rather than hadron, so no longer needed.
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/force.inc"
C----------------------------------------------------------------------
      RETURN
      END
