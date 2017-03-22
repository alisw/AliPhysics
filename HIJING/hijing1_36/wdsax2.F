* $Id$
C
C
	FUNCTION WDSAX2(X)
C     			********THREE PARAMETER WOOD SAXON 
C                               FOR  TARGET
C       HINT1(76)=R,HINT1(77)=D, HINT1(78)=W, HINT1(79)=FNORM
C
#include "hiparnt.inc"
      	WDSAX2=HINT1(79)*(1.+HINT1(78)*(X/HINT1(76))**2)/
     &       (1+EXP((X-HINT1(76))/HINT1(77)))
       	IF (HINT1(78).LT.0.) THEN
       		IF (X.GE.HINT1(76)/SQRT(ABS(HINT1(78)))) WDSAX2=0.
       	ENDIF
      	RETURN
      	END
