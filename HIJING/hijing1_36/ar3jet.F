* $Id$


	SUBROUTINE AR3JET(S,X1,X3,JL)
C     
#define BLANKET_SAVE
#include "hiparnt.inc"
#include "lujets_hijing.inc"
	SAVE
C     
	C=1./3.
      	IF(K(JL,2).NE.21 .AND. K(JL+1,2).NE.21) C=8./27.
      	EXP1=3
      	EXP3=3
      	IF(K(JL,2).NE.21) EXP1=2
      	IF(K(JL+1,2).NE.21) EXP3=2
      	A=0.24**2/S
      	YMA=ALOG(.5/SQRT(A)+SQRT(.25/A-1))
      	D=4.*C*YMA
      	SM1=P(JL,5)**2/S
      	SM3=P(JL+1,5)**2/S
      	XT2M=(1.-2.*SQRT(SM1)+SM1-SM3)*(1.-2.*SQRT(SM3)-SM1+SM3)
      	XT2M=MIN(.25,XT2M)
      	NTRY=0
1     	IF(NTRY.EQ.5000) THEN
	        X1=.5*(2.*SQRT(SM1)+1.+SM1-SM3)
		X3=.5*(2.*SQRT(SM3)+1.-SM1+SM3)
		RETURN
      	ENDIF
      	NTRY=NTRY+1
     
      	XT2=A*(XT2M/A)**(RLU_HIJING(0)**(1./D))
     
      	YMAX=ALOG(.5/SQRT(XT2)+SQRT(.25/XT2-1.))
      	Y=(2.*RLU_HIJING(0)-1.)*YMAX
      	X1=1.-SQRT(XT2)*EXP(Y)
      	X3=1.-SQRT(XT2)*EXP(-Y)
      	X2=2.-X1-X3
      	NEG=0
      	IF(K(JL,2).NE.21 .OR. K(JL+1,2).NE.21) THEN
        IF((1.-X1)*(1.-X2)*(1.-X3)-X2*SM1*(1.-X1)-X2*SM3*(1.-X3).
     &  LE.0..OR.X1.LE.2.*SQRT(SM1)-SM1+SM3.OR.X3.LE.2.*SQRT(SM3)
     &  -SM3+SM1) NEG=1
        X1=X1+SM1-SM3
        X3=X3-SM1+SM3
     	ENDIF
      	IF(NEG.EQ.1) GOTO 1
     
      	FG=2.*YMAX*C*(X1**EXP1+X3**EXP3)/D
      	XT2M=XT2
      	IF(FG.LT.RLU_HIJING(0)) GOTO 1
     
      	RETURN
      	END
