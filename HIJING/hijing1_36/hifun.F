* $Id$
C
C
C The next three subroutines are for Monte Carlo generation 
C according to a given function FHB. One calls first HIFUN 
C with assigned channel number I, low and up limits. Then to 
C generate the distribution one can call HIRND(I) which gives 
C you a random number generated according to the given function.
C 
	SUBROUTINE HIFUN(I,XMIN,XMAX,FHB)
#include "hijhb.inc"
	EXTERNAL FHB
	FNORM=GAUSS1(FHB,XMIN,XMAX,0.001)
	DO 100 J=1,201
		XX(I,J)=XMIN+(XMAX-XMIN)*(J-1)/200.0
		XDD=XX(I,J)
		RR(I,J)=GAUSS1(FHB,XMIN,XDD,0.001)/FNORM
100	CONTINUE
	RETURN
	END
