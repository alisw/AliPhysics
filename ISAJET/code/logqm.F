#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGQM(IERR)
C
C     Set and check limits for gamma*/W/Z0/Higgs mass range
C     Ver 7.14: Use HMASS+-5*HGAM for MSSM default range
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/jetlim.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/const.inc"
#include "isajet/dylim.inc"
#include "isajet/keys.inc"
#include "isajet/q1q2.inc"
#include "isajet/hcon.inc"
#include "isajet/xmssm.inc"
C
      REAL UNDEF
      INTEGER IERR
      DATA UNDEF/-.9E9/
C
      LOGQM=.TRUE.
      FIXQM=.FALSE.
      IF(QMIN.LT.UNDEF.AND.QMAX.LT.UNDEF) THEN
        IF(KEYS(7).AND.GOMSSM) THEN
C          For MSSM Higgs, set default limits around Higgs
          QMAX=HMASS+5*HGAM
          QMIN=HMASS-5*HGAM
        ELSE
C          Set default QMW limits if none were set.
          QMAX=0.2*ECM
          QMIN=0.05*ECM
        ENDIF
        CALL LOGERR(0,1,IERR)
      ENDIF
      IF(QMAX.LT.UNDEF) FIXQM=.TRUE.
      IF(FIXQM) THEN
        QMW=QMIN
        QMAX=QMIN
      ENDIF
C
      RETURN
      END
