* $Id$
C
	FUNCTION FGP1(X)
#include "pact_hijing.inc"
#include "hiparnt.inc"
        EXTERNAL FGP2
        PHI=X
        FGP1=2.0*GAUSS3(FGP2,0.0,HIPR1(34),0.01)
	RETURN
	END
