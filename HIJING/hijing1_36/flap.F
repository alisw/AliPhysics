* $Id$
C
C
	FUNCTION FLAP(X)
#include "pact_hijing.inc"
#include "hiparnt.inc"
        EXTERNAL FGP1
        B1=X
        FLAP=GAUSS2(FGP1,0.0,2.0*HIPR1(40),0.01)
	RETURN
	END
