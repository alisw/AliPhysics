#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGYTH(IERR)
C
C       SET AND CHECK LIMITS FOR JET Y AND THETA
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
      LOGYTH=.TRUE.
C
      DO 30 I=1,NJET
      FIXYJ(I)=.FALSE.
C
      IF(FIXP(I).AND.FIXPT(I)) THEN
        STH(I)=PT(I)/P(I)
        CTHS(1,I)=SQRT(1.-STH(I)**2)
        CTHS(2,I)=-CTHS(1,I)
        THS(1,I)=ATAN2(STH(I),CTHS(1,I))
        THS(2,I)=ATAN2(STH(I),CTHS(2,I))
        YJS(1,I)=-ALOG(TAN(THS(1,I)/2.))
        YJS(2,I)=-ALOG(TAN(THS(2,I)/2.))
        XJS(1,I)=P(I)*CTHS(1,I)/HALFE
        XJS(2,I)=P(I)*CTHS(2,I)/HALFE
        YJMAX(I)=YJS(2,I)
        THMAX(I)=THS(1,I)
        THMIN(I)=THS(2,I)
        IF(YJMIN(I).EQ.YJMAX(I)) FIXYJ(I)=.TRUE.
      ENDIF
C
C
      IF(YJMIN(I).LT.UNDEF.AND.YJMAX(I).LT.UNDEF) THEN
C
        IF(THMIN(I).LT.UNDEF.AND.THMAX(I).LT.UNDEF) THEN
          YJMAX(I)=ACOSH(HALFE/PTMIN(I))
          YJMIN(I)=-YJMAX(I)
          THMIN(I)=2.*ATAN(EXP(-YJMAX(I)))
          THMAX(I)=2.*ATAN(EXP(-YJMIN(I)))
        ENDIF
C
        IF(THMAX(I).LT.UNDEF) FIXYJ(I)=.TRUE.
        IF(THMIN(I).LT.UNDEF) THMIN(I)=.001
        IF(FIXYJ(I)) THMAX(I)=THMIN(I)
        YJMIN(I)=-ALOG(TAN(THMAX(I)/2.))
        YJMAX(I)=-ALOG(TAN(THMIN(I)/2.))
        THMIN(I)=2.*ATAN(EXP(-YJMAX(I)))
        THMAX(I)=2.*ATAN(EXP(-YJMIN(I)))
      ENDIF
C
C
      IF(YJMAX(I).LT.UNDEF) FIXYJ(I)=.TRUE.
      IF(YJMIN(I).LT.UNDEF) YJMIN(I)=-YJMAX(I)
      IF(FIXYJ(I)) YJMAX(I)=YJMIN(I)
      THMIN(I)=2.*ATAN(EXP(-YJMAX(I)))
      THMAX(I)=2.*ATAN(EXP(-YJMIN(I)))
C
      IF(FIXYJ(I)) THEN
        YJ(I)=YJMIN(I)
        TH(I)=THMIN(I)
        STH(I)=SIN(TH(I))
        CTH(I)=COS(TH(I))
        IF(FIXPT(I)) P(I)=PT(I)/STH(I)
        IF(FIXP(I)) PT(I)=P(I)*STH(I)
C
        IF((FIXP(I).OR.FIXPT(I))) THEN
          XJ(I)=P(I)*CTH(I)/HALFE
          XJMIN(I)=XJ(I)
          XJMAX(I)=XJ(I)
        ENDIF
C
      ENDIF
C
C          CHECK PT LIMITS WITH P AND THETA LIMITS
      IF(.NOT.FIXPT(I)) THEN
        THETA1=AMIN1(THMIN(I),PI-THMAX(I))
        THETA2=HALFPI
        IF(THMAX(I).LT.HALFPI) THETA2=THMAX(I)
        IF(THMIN(I).GT.HALFPI) THETA2=THMIN(I)
        PT1=PMIN(I)*SIN(THETA1)
        PTMIN(I)=AMAX1(PTMIN(I),PT1)
        PT2=PMAX(I)*SIN(THETA2)
        PTMAX(I)=AMIN1(PTMAX(I),PT2)
      ENDIF
C
   30 CONTINUE
C
      RETURN
      END
