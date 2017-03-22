* $Id$


C
C*******************************************************************
C	make  boost and rotation to entries from IMIN to IMAX
C*******************************************************************
	SUBROUTINE ATROBO(THE,PHI,BEX,BEY,BEZ,IMIN,IMAX,IERROR)
#define BLANKET_SAVE
#include "lujets_hijing.inc"
	DIMENSION ROT(3,3),PV(3)
	DOUBLE PRECISION DP(4),DBEX,DBEY,DBEZ,DGA,DGA2,DBEP,DGABEP
	SAVE
	IERROR=0
     
      	IF(IMIN.LE.0 .OR. IMAX.GT.N .OR. IMIN.GT.IMAX) RETURN

      	IF(THE**2+PHI**2.GT.1E-20) THEN
C...ROTATE (TYPICALLY FROM Z AXIS TO DIRECTION THETA,PHI)
	   ROT(1,1)=COS(THE)*COS(PHI)
	   ROT(1,2)=-SIN(PHI)
	   ROT(1,3)=SIN(THE)*COS(PHI)
	   ROT(2,1)=COS(THE)*SIN(PHI)
	   ROT(2,2)=COS(PHI)
	   ROT(2,3)=SIN(THE)*SIN(PHI)
	   ROT(3,1)=-SIN(THE)
	   ROT(3,2)=0.
	   ROT(3,3)=COS(THE)
	   DO 120 I=IMIN,IMAX
C**************	   IF(MOD(K(I,1)/10000,10).GE.6) GOTO 120
	      DO 100 J=1,3
 100		 PV(J)=P(I,J)
		 DO 110 J=1,3
 110		    P(I,J)=ROT(J,1)*PV(1)+ROT(J,2)*PV(2)
     &                     +ROT(J,3)*PV(3)
 120		 CONTINUE
	ENDIF
     
      	IF(BEX**2+BEY**2+BEZ**2.GT.1E-20) THEN
C...LORENTZ BOOST (TYPICALLY FROM REST TO MOMENTUM/ENERGY=BETA)
        	DBEX=BEX
        	DBEY=BEY
        	DBEZ=BEZ
		DGA2=1D0-DBEX**2-DBEY**2-DBEZ**2
		IF(DGA2.LE.0D0) THEN
			IERROR=1
			RETURN
		ENDIF
        	DGA=1D0/DSQRT(DGA2)
        	DO 140 I=IMIN,IMAX
C*************	   IF(MOD(K(I,1)/10000,10).GE.6) GOTO 140
        	   DO 130 J=1,4
130   		   DP(J)=P(I,J)
        	   DBEP=DBEX*DP(1)+DBEY*DP(2)+DBEZ*DP(3)
        	   DGABEP=DGA*(DGA*DBEP/(1D0+DGA)+DP(4))
        	   P(I,1)=DP(1)+DGABEP*DBEX
        	   P(I,2)=DP(2)+DGABEP*DBEY
        	   P(I,3)=DP(3)+DGABEP*DBEZ
        	   P(I,4)=DGA*(DP(4)+DBEP)
140   		CONTINUE
      	ENDIF
     
      	RETURN
      	END
