#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGPHI(IERR,DELPH)
C
C        SET AND CHECK LIMITS FOR JET PHI
C
#include "isajet/itapes.inc"
#include "isajet/jetlim.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/const.inc"
#include "isajet/dylim.inc"
#include "isajet/keys.inc"
#include "isajet/q1q2.inc"
      DIMENSION DELPH(3)
      DATA UNDEF/-.9E9/
C
      LOGPHI=.TRUE.
C
C
      DO 50 I=1,NJET
      FIXPHI(I)=.FALSE.
C
      IF(PHIMAX(I).LT.UNDEF.AND.PHIMIN(I).LT.UNDEF) THEN
        PHIMIN(I)=0.
        PHIMAX(I)=2.*PI
        DELPH(I)=PHIMAX(I)
      ELSE
        IF(PHIMAX(I).LT.UNDEF) FIXPHI(I)=.TRUE.
C
        IF(FIXPHI(I)) THEN
          PHI(I)=PHIMIN(I)
          PHIMAX(I)=PHIMIN(I)
          IF(KEYS(3).AND.I.LT.3) THEN
            LOGPHI=.FALSE.
            CALL LOGERR(105,I,IERR)
          ENDIF
          IF(I.EQ.2) THEN
            FIXPHI(1)=.TRUE.
            PHIMIN(1)=PHIMIN(2)
            PHIMAX(1)=PHIMIN(1)
          ENDIF
        ENDIF
C
        DELPH(I)=PHIMAX(I)-PHIMIN(I)
C
        IF(DELPH(I).GT.2.*PI.OR.DELPH(I).LT.0) THEN
          LOGPHI=.FALSE.
          CALL LOGERR(8,I,IERR)
        ENDIF
C
      ENDIF
C
   50 CONTINUE
C
C
      IF(KEYS(1).AND.DELPH(1).GT.DELPH(2)) THEN
        PHIMIN(1)=PHIMIN(2)+PI
        PHIMAX(1)=PHIMIN(1)+DELPH(2)
      ENDIF
C
      RETURN
      END
