#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGTHW(IERR)
C
C        SET AND CHECK THETA LIMITS FOR W(Z0)
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
C            INVERSE HYPERBOLIC COSINE FUNCTION
      ACOSH(X)=ALOG(X+SQRT(X**2-1.0))
C            INVERSE HYPERBOLIC SINE FUNCTION
      ASINH(X)=ALOG(X+SQRT(X**2+1.0))
C
      HALFPI=PI/2.
      LOGTHW=.TRUE.
C
      IF(THWMIN.LT.UNDEF.AND.THWMAX.LT.UNDEF) THEN
          THWMIN=2.*ATAN(EXP(-YWMAX))
          THWMAX=2.*ATAN(EXP(-YWMIN))
      ELSEIF(THWMIN.GT.UNDEF) THEN
        IF(THWMAX.GT.UNDEF) THEN
          LOGTHW=.FALSE.
          CALL LOGERR(113,1,IERR)
        ELSE
          TAMIN=TAN(THWMIN)
          TAMAX=TAN(THWMAX)
          IF(THWMIN.LT.HALFPI)
     1    YWMX=ASINH(QTMAX/SQRT(QTMAX**2+QMIN**2)/TAMIN)
          IF(THWMIN.GE.HALFPI)
     1    YWMX=ASINH(QTMIN/SQRT(QTMIN**2+QMAX**2)/TAMIN)
          IF(THWMAX.GT.HALFPI)
     1    YWMN=ASINH(QTMAX/SQRT(QTMAX**2+QMIN**2)/TAMAX)
          IF(THWMAX.LT.HALFPI)
     1    YWMN=ASINH(QTMIN/SQRT(QTMIN**2+QMAX**2)/TAMAX)
          IF(YWMIN.LT.YWMN) YWMIN=YWMN
          IF(YWMAX.GT.YWMX) YWMAX=YWMX
          IF(FIXYW.AND.(YW.LT.YWMIN.OR.YW.GT.YWMAX)) THEN
            CALL LOGERR(102,1,IERR)
            LOGTHW=.FALSE.
          ENDIF
        ENDIF
      ENDIF
C
      IF(YWMIN.LT.0) THWMAX=ATAN2(QTMIN,SQRT(QTMIN**2+QMAX**2)*
     1              SINH(YWMIN))
      IF(YWMIN.GE.0) THWMAX=ATAN2(QTMAX,SQRT(QTMAX**2+QMIN**2)*
     1               SINH(YWMIN))
      IF(YWMAX.GE.0) THWMIN=ATAN2(QTMIN,SQRT(QTMIN**2+QMAX**2)*
     U               SINH(YWMAX))
      IF(YWMAX.LT.0) THWMIN=ATAN2(QTMAX,SQRT(QTMAX**2+QMIN**2)*
     1               SINH(YWMAX))
C
      RETURN
      END
