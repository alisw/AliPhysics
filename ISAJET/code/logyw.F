#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGYW(IERR)
C
C       SET AND CHECK Y LIMITS FOR W(Z0)
C
#include "isajet/itapes.inc"
#include "isajet/jetlim.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/const.inc"
#include "isajet/dylim.inc"
#include "isajet/keys.inc"
#include "isajet/q1q2.inc"
      LOGICAL COMB(2)
      DATA UNDEF/-.9E9/
C
C            INVERSE HYPERBOLIC COSINE FUNCTION
      ACOSH(X)=ALOG(X+SQRT(X**2-1.0))
C            INVERSE HYPERBOLIC SINE FUNCTION
      ASINH(X)=ALOG(X+SQRT(X**2+1.0))
      YW=1.0
      LOGYW=.TRUE.
      FIXYW=.FALSE.
C
      IF(YWMIN.LT.UNDEF.AND.YWMAX.LT.UNDEF) THEN
        YWMAX=ACOSH(HALFE/SQRT(QTMIN**2+QMIN**2))
        YWMIN=-YWMAX
      ENDIF
C
      IF(YWMAX.LT.UNDEF) THEN
        FIXYW=.TRUE.
        YW=YWMIN
        YWMAX=YWMIN
      ENDIF
C
      YWMX=ACOSH(HALFE/SQRT(QTMIN**2+QMIN**2))
      YWMN=-YWMX
      COMB(1)=YWMX.LT.YWMAX
      COMB(2)=YWMN.GT.YWMIN
C
      IF(FIXYW.AND.(COMB(1).OR.COMB(2))) THEN
        LOGYW=.FALSE.
        CALL LOGERR(102,1,IERR)
      ENDIF
C
      IF(.NOT.FIXYW) THEN
        IF(COMB(1)) YWMAX=YWMX
        IF(COMB(2)) YWMIN=YWMN
      ENDIF
C
      RETURN
      END
