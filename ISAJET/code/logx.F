#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGX(IERR)
C
C         SET AND CHECK LIMITS FOR JET FEYNMAN X
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
      HALFPI=PI/2.
      LOGX=.TRUE.
C
      DO 40 I=1,NJET
      FIXXJ(I)=.FALSE.
      IF(FIXYJ(I).AND.(FIXP(I).OR.FIXPT(I)))FIXXJ(I)=.TRUE.
      IF(FIXXJ(I)) GOTO 40
C
      IF(XJMIN(I).LT.UNDEF.AND.XJMAX(I).LT.UNDEF) THEN
        XJMAX(I)=1.0
        XJMIN(I)=-1.0
      ENDIF
C
      IF(XJMAX(I).LT.UNDEF) FIXXJ(I)=.TRUE.
      IF(FIXXJ(I)) XJMAX(I)=XJMIN(I)
C
      IF(.NOT.FIXXJ(I)) THEN
        IF(THMIN(I).LT.HALFPI) X1=PMAX(I)*COS(THMIN(I))/HALFE
        IF(THMIN(I).GE.HALFPI) X1=PMIN(I)*COS(THMIN(I))/HALFE
        IF(THMAX(I).GT.HALFPI) X2=PMAX(I)*COS(THMAX(I))/HALFE
        IF(THMAX(I).LT.HALFPI) X2=PMIN(I)*COS(THMAX(I))/HALFE
        IF(X1.LT.XJMAX(I)) XJMAX(I)=X1
        IF(X2.GT.XJMIN(I)) XJMIN(I)=X2
      ELSE
C
        XJ(I)=XJMIN(I)
C
        IF(FIXP(I)) THEN
          CTH(I)=XJ(I)*HALFE/P(I)
          IF(ABS(CTH(I)).LE.1.0) THEN
            STH(I)=SQRT(1.-CTH(I)**2)
            TH(I)=ATAN2(STH(I),CTH(I))
            YJ(I)=-ALOG(TAN(TH(I)/2.))
            FIXYJ(I)=.TRUE.
            PT(I)=P(I)*STH(I)
            FIXPT(I)=.TRUE.
            YJMIN(I)=YJ(I)
            YJMAX(I)=YJ(I)
            PTMIN(I)=PT(I)
            PTMAX(I)=PT(I)
          ELSE
            LOGX=.FALSE.
            CALL LOGERR(5,I,IERR)
          ENDIF
        ENDIF
C
        IF(FIXPT(I)) THEN
          TH(I)=ATAN(PT(I)/XJ(I)/HALFE)
          FIXYJ(I)=.TRUE.
          YJ(I)=-ALOG(TAN(TH(I)/2.))
          CTH(I)=COS(TH(I))
          STH(I)=SIN(TH(I))
          P(I)=PT(I)/STH(I)
          FIXP(I)=.TRUE.
          YJMIN(I)=YJ(I)
          YJMAX(I)=YJ(I)
          PMAX(I)=P(I)
          PMIN(I)=P(I)
        ENDIF
C
        IF(FIXYJ(I)) THEN
          FIXPT(I)=.TRUE.
          P(I)=XJ(I)*HALFE/CTH(I)
          PT(I)=P(I)*STH(I)
          FIXP(I)=.TRUE.
          PTMIN(I)=PT(I)
          PTMAX(I)=PT(I)
          PMAX(I)=P(I)
          PMIN(I)=P(I)
        ENDIF
C
      ENDIF
C
   40 CONTINUE
C
      RETURN
      END
