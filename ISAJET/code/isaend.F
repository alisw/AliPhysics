#include "isajet/pilot.h"
      SUBROUTINE ISAEND
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-       Terminate an ISAJET run
C-
C-   Created   4-FEB-1988   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
      CALL TIMER(2)
      CALL GETTOT(.TRUE.)
  999 RETURN
      END
