#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGPHW(IERR,DELPH)
C
C       SET AND CHECK LIMITS FOR W(Z0) PHI
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
      LOGPHW=.TRUE.
      FIXPHW=.FALSE.
C
      IF(FIXPHI(3)) THEN
          FIXPHW=.TRUE.
          PHIW=AMOD(PHI(3)+PI,2.*PI)
      ELSEIF(PHWMIN.LT.UNDEF.AND.PHWMAX.LT.UNDEF) THEN
          PHWMIN=0.
          PHWMAX=2.*PI
      ELSEIF(PHWMAX.LT.UNDEF) THEN
          FIXPHW=.TRUE.
          PHW=PHWMIN
          FIXPHI(3)=.TRUE.
          PHWMAX=PHWMIN
          PHI(3)=PHIW+PI
          PHIMIN(3)=PHIW
          PHIMAX(3)=PHIW
      ENDIF
C
      DELPHW=PHWMAX-PHWMIN
C
      IF(DELPHW.LT.0.OR.DELPHW.GT.2.*PI) THEN
        CALL LOGERR(110,1,IERR)
        LOGPHW=.FALSE.
      ENDIF
C
      IF(DELPHW.LE.DELPH(3)) THEN
        PHIMIN(3)=PHWMIN+PI
        PHIMAX(3)=PHIMIN(3)+DELPHW
      ELSE
        PHWMIN=PHIMIN(3)+PI
        PHWMAX=PHWMIN+DELPH(3)
      ENDIF
C
      RETURN
      END
