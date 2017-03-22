* $Id$
C
C
C*******************************************************************
CThis subroutine performs elastic scattering between two nucleons
C
C*******************************************************************
	SUBROUTINE HIJELS(PSC1,PSC2)
	IMPLICIT DOUBLE PRECISION(D)
	DIMENSION PSC1(5),PSC2(5)
#define BLANKET_SAVE
#include "hiparnt.inc"
	SAVE
C
	CC=1.0-HINT1(12)/HINT1(13)
	RR=(1.0-CC)*HINT1(13)/HINT1(12)/(1.0-HIPR1(33))-1.0
	BB=0.5*(3.0+RR+SQRT(9.0+10.0*RR+RR**2))
	EP=SQRT((PSC1(1)-PSC2(1))**2+(PSC1(2)-PSC2(2))**2
     &		+(PSC1(3)-PSC2(3))**2)
	IF(EP.LE.0.1) RETURN
	ELS0=98.0/EP+52.0*(1.0+RR)**2
	PCM1=PSC1(1)+PSC2(1)
	PCM2=PSC1(2)+PSC2(2)
	PCM3=PSC1(3)+PSC2(3)
	ECM=PSC1(4)+PSC2(4)
	AM1=PSC1(5)**2
	AM2=PSC2(5)**2
	AMM=ECM**2-PCM1**2-PCM2**2-PCM3**2
	IF(AMM.LE.PSC1(5)+PSC2(5)) RETURN
C		********elastic scattering only when approaching
C				to each other
	PMAX=(AMM**2+AM1**2+AM2**2-2.0*AMM*AM1-2.0*AMM*AM2
     &			-2.0*AM1*AM2)/4.0/AMM
	PMAX=ABS(PMAX)
20	TT=RLU_HIJING(0)*MIN(PMAX,1.5)
	ELS=98.0*EXP(-2.8*TT)/EP
     &		+52.0*EXP(-9.2*TT)*(1.0+RR*EXP(-4.6*(BB-1.0)*TT))**2
	IF(RLU_HIJING(0).GT.ELS/ELS0) GO TO 20
	PHI=2.0*HIPR1(40)*RLU_HIJING(0)
C
	DBX=PCM1/ECM
	DBY=PCM2/ECM
	DBZ=PCM3/ECM
        DB=SQRT(DBX**2+DBY**2+DBZ**2)
        IF(DB.GT.0.99999999D0) THEN 
          DBX=DBX*(0.99999999D0/DB) 
          DBY=DBY*(0.99999999D0/DB) 
          DBZ=DBZ*(0.99999999D0/DB) 
          DB=0.99999999D0   
	  WRITE(6,*) ' (HIJELS) boost vector too large' 
C		********Rescale boost vector if too close to unity. 
        ENDIF   
        DGA=1D0/SQRT(1D0-DB**2)      
C
	DP1=SQRT(TT)*SIN(PHI)
	DP2=SQRT(TT)*COS(PHI)
	DP3=SQRT(PMAX-TT)
	DP4=SQRT(PMAX+AM1)
        DBP=DBX*DP1+DBY*DP2+DBZ*DP3   
        DGABP=DGA*(DGA*DBP/(1D0+DGA)+DP4) 
        PSC1(1)=DP1+DGABP*DBX
        PSC1(2)=DP2+DGABP*DBY  
        PSC1(3)=DP3+DGABP*DBZ  
        PSC1(4)=DGA*(DP4+DBP)    
C	
	DP1=-SQRT(TT)*SIN(PHI)
	DP2=-SQRT(TT)*COS(PHI)
	DP3=-SQRT(PMAX-TT)
	DP4=SQRT(PMAX+AM2)
        DBP=DBX*DP1+DBY*DP2+DBZ*DP3   
        DGABP=DGA*(DGA*DBP/(1D0+DGA)+DP4) 
        PSC2(1)=DP1+DGABP*DBX
        PSC2(2)=DP2+DGABP*DBY  
        PSC2(3)=DP3+DGABP*DBZ  
        PSC2(4)=DGA*(DP4+DBP)
	RETURN
	END
