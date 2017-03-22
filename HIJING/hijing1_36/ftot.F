* $Id$
C
C
C
C
	FUNCTION FTOT(X)
#include "hiparnt.inc"
	OMG=OMG0(X)*(HIPR1(30)+HINT1(11))/HIPR1(31)/2.0
	FTOT=2.0*(1.0-EXP(-OMG))
	RETURN
	END
