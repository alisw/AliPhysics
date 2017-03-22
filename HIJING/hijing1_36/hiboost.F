* $Id$
C
C
C
C
	SUBROUTINE HIBOOST
      	IMPLICIT DOUBLE PRECISION(D)  
#define BLANKET_SAVE
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "hiparnt.inc"
	SAVE
	DO 100 I=1,N
	   DBETA=P(I,3)/P(I,4)
	   IF(ABS(DBETA).GE.1.D0) THEN
	      DB=HINT1(2)
	      IF(DB.GT.0.99999999D0) THEN 
C		********Rescale boost vector if too close to unity. 
		 WRITE(6,*) '(HIBOOT:) boost vector too large' 
		 DB=0.99999999D0
	      ENDIF 
	      DGA=1D0/SQRT(1D0-DB**2)
	      DP3=P(I,3)
	      DP4=P(I,4)
	      P(I,3)=(DP3+DB*DP4)*DGA  
	      P(I,4)=(DP4+DB*DP3)*DGA  
	      GO TO 100
	   ENDIF
	   Y=0.5*DLOG((1.D0+DBETA)/(1.D0-DBETA))
	   AMT=SQRT(P(I,1)**2+P(I,2)**2+P(I,5)**2)
	   P(I,3)=AMT*SINH(Y+HINT1(3))
	   P(I,4)=AMT*COSH(Y+HINT1(3))
100	CONTINUE
	RETURN
	END
