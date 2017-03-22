* $Id$
C
C
C
	FUNCTION GPHOTON(Y1,Y2,PT2)
	IMPLICIT REAL*8 (A-H,O-Z)
#define BLANKET_SAVE
#include "hiparnt.inc"
	DIMENSION F(2,7)
	SAVE
	XT=2.0*DSQRT(PT2)/HINT1(1)
	X1=0.50*XT*(DEXP(Y1)+DEXP(Y2))
	X2=0.50*XT*(DEXP(-Y1)+DEXP(-Y2))
	Z=DSQRT(1.D0-XT**2/X1/X2)
	SS=X1*X2*HINT1(1)**2
	T=-(1.0-Z)/2.0
	U=-(1.0+Z)/2.0
	AF=3.0
	DLAM=HIPR1(15)
	APH=12.0*3.1415926/(33.0-2.0*AF)/DLOG(PT2/DLAM**2)
	APHEM=1.0/137.0
C
	CALL PARTON(F,X1,X2,PT2)
C
	G11=-(U**2+1.0)/U/3.0*F(1,7)*(4.0*F(2,1)+4.0*F(2,2)
     &      +F(2,3)+F(2,4)+F(2,5)+F(2,6))/9.0
	G12=-(T**2+1.0)/T/3.0*F(2,7)*(4.0*F(1,1)+4.0*F(1,2)
     &      +F(1,3)+F(1,4)+F(1,5)+F(1,6))/9.0
	G2=8.0*(U**2+T**2)/U/T/9.0*(4.0*F(1,1)*F(2,2)
     &     +4.0*F(1,2)*F(2,1)+F(1,3)*F(2,4)+F(1,4)*F(2,3)
     &     +F(1,5)*F(2,6)+F(1,6)*F(2,5))/9.0
C
	GPHOTON=(G11+G12+G2)*HIPR1(23)*3.14159*APH*APHEM/SS**2
	RETURN
	END
