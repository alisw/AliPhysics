#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGQT(IERR)
C
C      SET AND CHECK W(Z0) PT RANGE
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
      LOGQT=.TRUE.
      FIXQT=.FALSE.
      IF(QTMIN.LT.UNDEF.AND.QTMAX.LT.UNDEF) THEN
        QTMAX=PTMAX(3)
        QTMIN=PTMIN(3)
C          SET DEFAULT QTW LIMITS IF NONE WERE SET
        IF(QTMAX.GT.0.99*HALFE) THEN
          NJET=2
          QTMIN=0.
          QTMAX=0.
          QTW=0.
          STDDY=.TRUE.
          FIXQT=.TRUE.
          PTMIN(3)=0.
          PTMAX(3)=0.
          FIXPT(3)=.TRUE.
          CALL LOGERR(0,1,IERR)
        ENDIF
      ELSEIF(FIXPT(3)) THEN
        QTW=PT(3)
        QTMIN=PTMIN(3)
        QTMAX=QTMIN
        FIXQT=.TRUE.
      ELSEIF(QTMAX.LT.UNDEF) THEN
        FIXQT=.TRUE.
        QTW=QTMIN
        QTMAX=QTMIN
        FIXPT(3)=.TRUE.
        PT(3)=QTW
        PTMIN(3)=QTMIN
        PTMAX(3)=QTMAX
      ELSE
        IF(QTMAX.LT.PTMAX(3)) PTMAX(3)=QTMAX
        IF(QTMIN.GT.PTMIN(3)) PTMIN(3)=QTMIN
        IF(QTMAX.GT.PTMAX(3)) QTMAX=PTMAX(3)
        IF(QTMIN.LT.PTMIN(3)) QTMIN=PTMIN(3)
      ENDIF
C
      RETURN
      END
