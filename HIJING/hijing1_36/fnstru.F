* $Id$
C
C
C
	FUNCTION FNSTRU(X)
#include "hiparnt.inc"
	FNSTRU=(1.0-X)**HIPR1(44)/
     &		(X**2+HIPR1(45)**2/HINT1(1)**2)**HIPR1(46)
	RETURN
	END
