#include "isajet/pilot.h"
      DOUBLE PRECISION FUNCTION DINCGM(A,X,EPS)
C***********************************************************************
C* Series expansion of incomplete gamma function, from Abramowitz and  *
C* Stegun. A and X are the two arguments, while EPS is the relative    *
C* precision. More accurately, if X > 0, EPS is the ratio of the last  *
C* term in the series and the sum; note that for X > 0, the series is  *
C* alternating. For X < 0, this ratio is required to be < EPS/100.     *
C***********************************************************************

      DOUBLE PRECISION A,X,EPS,SUM,TERM,XN
C      IF(DABS(A).LT.1.D-10) THEN
C        WRITE(*,*) ' Function diverges for A = 0!'
C        DINCGM = 1.D50
C        RETURN
C      ENDIF
      SUM = 1.D0/A
      TERM = 1.D0    !Term for n = 0
      H = 1.D0
      XN = 1.D0

    1 H = -H*X/XN
      TERM = H/(A+XN)
      XN = XN+1.D0
      SUM = SUM+TERM
      IF((DABS(TERM/SUM).GT.EPS.AND.X.GE.0.D0).OR.
     &   (DABS(TERM/SUM).GT.1.D-2*EPS.AND.X.LE.0.D0)) GOTO 1

      DINCGM = SUM*(X**A)
      RETURN
      END
