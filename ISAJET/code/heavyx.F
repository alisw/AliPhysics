#include "isajet/pilot.h"
      SUBROUTINE HEAVYX(X,EPS)
C
C          GENERATE X FOR HEAVY PARTICLE FRAGMENTATION ACCORDING TO
C          THE PETERSON FORM
C          D(X)=1/(X*(1-1/X-EPS/(1-X))**2)
C              =D0(X)*D1(X)*D2(X)
C          D0(X)=(1-X)**2/((1-X)**2+EPS)**2
C          D1(X)=X
C          D2(X)=(((1-X)**2+EPS)/((1-X)**2+EPS*X))**2
C          USING X=1-Y**POW
C
      DATA ALN4/1.3863/
C
C          CHOOSE POW FOR X=1-Y**POW.
C          GENERATE FLAT IN X IF EPS>1.
      IF(EPS.LT.1.) THEN
        POW=ALOG((3.+EPS)/EPS)/ALN4
        YMX=(EPS*(3.*POW-1.)/(POW+1.))**(.5/POW)
        ZMX=1-YMX**POW
        D0MX=(1-ZMX)**2/((1.-ZMX)**2+EPS)**2*POW*YMX**(POW-1.)
        D2MX=2./(2.-SQRT(EPS))
      ELSE
        POW=1.
        ZMX=0.
        D0MX=(1.-ZMX)**2/((1.-ZMX)**2+EPS)**2
        D2MX=1.+EPS
      ENDIF
C
C          GENERATE Z ACCORDING TO (1-Z)**2/((1-Z)**2+EPS*Z)**2
1     CONTINUE
      Y=RANF()
      Z=1.-Y**POW
C
      D0Z=(1.-Z)**2/((1.-Z)**2+EPS)**2*POW*Y**(POW-1.)
      IF(D0Z.LT.RANF()*D0MX) GO TO 1
C
C          CHECK REMAINING FACTORS
      D1=Z
      D2=(((1.-Z)**2+EPS)/((1.-Z)**2+EPS*Z))**2
      IF(D1*D2.LT.RANF()*D2MX) GO TO 1
C
C          GOOD X
      X=Z
      RETURN
      END
