* $Id$
C
C
C
	FUNCTION FNSTRUM(X)
#include "hiparnt.inc"
	FNSTRUM=1.0/((1.0-X)**2+HIPR1(45)**2/HINT1(1)**2)**HIPR1(46)
     &          /(X**2+HIPR1(45)**2/HINT1(1)**2)**HIPR1(46)
	RETURN
	END
