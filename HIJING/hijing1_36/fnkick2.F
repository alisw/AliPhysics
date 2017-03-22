* $Id$
C
C
	FUNCTION FNKICK2(X)
#include "hiparnt.inc"
	FNKICK2=X*EXP(-2.0*X/HIPR1(42))
	RETURN
	END
