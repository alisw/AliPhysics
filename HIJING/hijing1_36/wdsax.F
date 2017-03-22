* $Id$
C
C
	FUNCTION WDSAX(X)
C     			********THREE PARAMETER WOOD SAXON
#include "wood_hijing.inc"
      	WDSAX=FNORM*(1.+W*(X/R)**2)/(1+EXP((X-R)/D))
       	IF (W.LT.0.) THEN
       		IF (X.GE.R/SQRT(ABS(W))) WDSAX=0.
       	ENDIF
      	RETURN
      	END
