#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGXW(IERR)
C
C       SET AND CHECK X LIMITS FOR W(Z0)
C
#include "isajet/itapes.inc"
#include "isajet/jetlim.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/const.inc"
#include "isajet/dylim.inc"
#include "isajet/keys.inc"
#include "isajet/q1q2.inc"
      DATA UNDEF/-.9E9/
C
      LOGXW=.TRUE.
      FIXXW=.FALSE.
C
      IF(XWMIN.LT.UNDEF.AND.XWMAX.LT.UNDEF) THEN
        XWMIN=-1.0
        XWMAX=1.0
      ELSEIF(XWMAX.GT.UNDEF) THEN
        FIXXW=.TRUE.
        XW=XWMIN
        XWMAX=XW
C            IF XW=0 THEN YW=0
        IF(XW.NE.0) THEN
          FIXYW=.TRUE.
          YW=0
          YWMIN=0
          YWMAX=0
        ENDIF
      ENDIF
C
C            IF YW=0 THAN XW=0
      IF(YW.EQ.0) THEN
        FIXXW=.TRUE.
        XW=0
        XWMAX=0
      ENDIF
C
      RETURN
      END
