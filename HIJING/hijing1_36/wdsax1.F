* $Id$
C
C
C
C
	FUNCTION WDSAX1(X)
C     			********THREE PARAMETER WOOD SAXON 
C                               FOR  PROJECTILE
C       HINT1(72)=R, HINT1(73)=D, HINT1(74)=W, HINT1(75)=FNORM
C
#include "hiparnt.inc"
      	WDSAX1=HINT1(75)*(1.+HINT1(74)*(X/HINT1(72))**2)/
     &       (1+EXP((X-HINT1(72))/HINT1(73)))
       	IF (HINT1(74).LT.0.) THEN
       		IF (X.GE.HINT1(72)/SQRT(ABS(HINT1(74)))) WDSAX1=0.
       	ENDIF
      	RETURN
      	END
