* $Id$
C*************************************************************


	SUBROUTINE ARORIE(S,X1,X3,JL)
C     
#define BLANKET_SAVE
#include "hiparnt.inc"
#include "lujets_hijing.inc"
	SAVE
C     
     	W=SQRT(S)
     	X2=2.-X1-X3
     	E1=.5*X1*W
     	E3=.5*X3*W
     	P1=SQRT(E1**2-P(JL,5)**2)
	P3=SQRT(E3**2-P(JL+1,5)**2)
	CBET=1.
	IF(P1.GT.0..AND.P3.GT.0.) CBET=(P(JL,5)**2
     &           +P(JL+1,5)**2+2.*E1*E3-S*(1.-X2))/(2.*P1*P3)
      	IF(ABS(CBET).GT.1.0) CBET=MAX(-1.,MIN(1.,CBET))
      	BET=ACOS(CBET)
     
C.....MINIMIZE PT1-SQUARED PLUS PT3-SQUARED.....
      	IF(P1.GE.P3) THEN
	   PSI=.5*ULANGL_HIJING(P1**2+P3**2*COS(2.*BET),-P3**2*SIN(2.
     1	*BET))
	   PT1=P1*SIN(PSI)
	   PZ1=P1*COS(PSI)
	   PT3=P3*SIN(PSI+BET)
	   PZ3=P3*COS(PSI+BET)
      	ELSE IF(P3.GT.P1) THEN
	   PSI=.5*ULANGL_HIJING(P3**2+P1**2*COS(2.*BET),-P1**2*SIN(2.
     1	*BET))
	   PT1=P1*SIN(BET+PSI)
	   PZ1=-P1*COS(BET+PSI)
	   PT3=P3*SIN(PSI)
	   PZ3=-P3*COS(PSI)
      	ENDIF
	
      	DEL=2.0*HIPR1(40)*RLU_HIJING(0)
      	P(JL,4)=E1
      	P(JL,1)=PT1*SIN(DEL)
      	P(JL,2)=-PT1*COS(DEL)
      	P(JL,3)=PZ1
      	P(JL+2,4)=E3
      	P(JL+2,1)=PT3*SIN(DEL)
      	P(JL+2,2)=-PT3*COS(DEL)
      	P(JL+2,3)=PZ3
      	P(JL+1,4)=W-E1-E3
      	P(JL+1,1)=-P(JL,1)-P(JL+2,1)
      	P(JL+1,2)=-P(JL,2)-P(JL+2,2)
      	P(JL+1,3)=-P(JL,3)-P(JL+2,3)
      	RETURN
      	END
