* $Id$
C
C
        FUNCTION FLAP1(X)
#include "pact_hijing.inc"
        R=SQRT(BB**2+X**2)
        FLAP1=WDSAX1(R)
        RETURN
        END
