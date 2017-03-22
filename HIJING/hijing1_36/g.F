* $Id$
C
C
C
C
	FUNCTION G(Y1,Y2,PT2)
	IMPLICIT REAL*8 (A-H,O-Z)
#define BLANKET_SAVE
#include "hiparnt.inc"
	DIMENSION F(2,7)
	SAVE
	XT=2.0*DSQRT(PT2)/HINT1(1)
	X1=0.50*XT*(DEXP(Y1)+DEXP(Y2))
	X2=0.50*XT*(DEXP(-Y1)+DEXP(-Y2))
	DUMMY = 1.D0-XT**2/X1/X2
	if (DUMMY.LT.0.D0) DUMMY = 0.D0
	Z=DSQRT(DUMMY)
	SS=X1*X2*HINT1(1)**2
	T=-(1.0-Z)/2.0
	U=-(1.0+Z)/2.0
	AF=3.0
	DLAM=HIPR1(15)
	APH=12.0*3.1415926/(33.0-2.0*AF)/DLOG(PT2/DLAM**2)
C
	CALL PARTON(F,X1,X2,PT2)
C
	G11=( (F(1,1)+F(1,2))*(F(2,3)+F(2,4)+F(2,5)+F(2,6))
     &      +(F(1,3)+F(1,4))*(F(2,5)+F(2,6)) )*SUBCRS1(T,U)
C
	G12=( (F(2,1)+F(2,2))*(F(1,3)+F(1,4)+F(1,5)+F(1,6))
     &      +(F(2,3)+F(2,4))*(F(1,5)+F(1,6)) )*SUBCRS1(U,T)
C
	G13=(F(1,1)*F(2,1)+F(1,2)*F(2,2)+F(1,3)*F(2,3)+F(1,4)*F(2,4)
     &      +F(1,5)*F(2,5)+F(1,6)*F(2,6))*(SUBCRS1(U,T)
     &      +SUBCRS1(T,U)-8.D0/T/U/27.D0)
C
	G2=(AF-1)*(F(1,1)*F(2,2)+F(2,1)*F(1,2)+F(1,3)*F(2,4)
     &     +F(2,3)*F(1,4)+F(1,5)*F(2,6)+F(2,5)*F(1,6))*SUBCRS2(T,U)
C
	G31=(F(1,1)*F(2,2)+F(1,3)*F(2,4)+F(1,5)*F(2,6))*SUBCRS3(T,U)
	G32=(F(2,1)*F(1,2)+F(2,3)*F(1,4)+F(2,5)*F(1,6))*SUBCRS3(U,T)
C
	G4=(F(1,1)*F(2,2)+F(2,1)*F(1,2)+F(1,3)*F(2,4)+F(2,3)*F(1,4)+
     1	F(1,5)*F(2,6)+F(2,5)*F(1,6))*SUBCRS4(T,U)
C
	G5=AF*F(1,7)*F(2,7)*SUBCRS5(T,U)
C
	G61=F(1,7)*(F(2,1)+F(2,2)+F(2,3)+F(2,4)+F(2,5)
     &      +F(2,6))*SUBCRS6(T,U)
	G62=F(2,7)*(F(1,1)+F(1,2)+F(1,3)+F(1,4)+F(1,5)
     &      +F(1,6))*SUBCRS6(U,T)
C
	G7=F(1,7)*F(2,7)*SUBCRS7(T,U)
C
	G=(G11+G12+G13+G2+G31+G32+G4+G5+G61+G62+G7)*HIPR1(17)*
     1	3.14159D0*APH**2/SS**2
	RETURN
	END
