* $Id$
C
C
C
	FUNCTION OMG0(X)
#include "hiparnt.inc"
#include "besel_hijing.inc"
	EXTERNAL BK
	X4=HIPR1(32)*SQRT(X)
	OMG0=HIPR1(32)**2*GAUSS2(BK,X4,X4+20.0,0.01)/96.0
	RETURN
	END
