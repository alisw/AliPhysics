#include "isajet/pilot.h"
      LOGICAL FUNCTION LOGMIJ(IERR)
C
C       Set and check limits for dijet masses.
C
C       Note we use the convention that not setting an upper limit
C       gives a fixed variable, even though that currently is not
C       implemented in N-jet phase space.
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
#include "isajet/mglims.inc"
C
      REAL AMLOW,UNDEF
      INTEGER I,J,IERR
      DATA AMLOW/1.0/
      DATA UNDEF/-.9E9/
C
      LOGMIJ=.TRUE.
C
      DO 100 I=1,MXLIM
        DO 101 J=I+1,MXLIM
          FIXMIJ(I,J)=.FALSE.
          FIXMIJ(J,I)=.FALSE.
          IF(AMIJMN(I,J).LT.UNDEF.AND.AMIJMX(I,J).LT.UNDEF) THEN
            AMIJMX(I,J)=ECM
            AMIJMX(J,I)=ECM
          ENDIF
          IF(AMIJMX(I,J).GT.ECM) THEN
            AMIJMX(I,J)=ECM
            AMIJMX(J,I)=ECM
          ENDIF
          IF(AMIJMX(I,J).LT.UNDEF) THEN
            AMIJMX(I,J)=AMIJMN(I,J)
            FIXMIJ(I,J)=.TRUE.
            AMIJMX(J,I)=AMIJMN(I,J)
            FIXMIJ(J,I)=.TRUE.
          ENDIF
          IF(AMIJMN(I,J).LT.UNDEF) THEN
            AMIJMN(I,J)=AMLOW
            AMIJMN(J,I)=AMLOW
          ENDIF
101     CONTINUE
100   CONTINUE
C
      RETURN
      END
