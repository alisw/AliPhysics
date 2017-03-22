* $Id$
    
C*********************************************************************  
    
      SUBROUTINE LUSPHE_HIJING(SPH,APL)    
    
C...Purpose: to perform sphericity tensor analysis to give sphericity,  
C...aplanarity and the related event axes.  
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
      DIMENSION SM(3,3),SV(3,3) 
    
C...Calculate matrix to be diagonalized.    
      NP=0  
      DO 100 J1=1,3 
      DO 100 J2=J1,3    
  100 SM(J1,J2)=0.  
      PS=0. 
      DO 120 I=1,N  
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 120  
      IF(MSTU(41).GE.2) THEN    
        KC=LUCOMP_HIJING(K(I,2))   
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.  
     &  KC.EQ.18) GOTO 120  
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE_HIJING(K(I,2))
     $       .EQ.0)GOTO 120    
      ENDIF 
      NP=NP+1   
      PA=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)    
      PWT=1.    
      IF(ABS(PARU(41)-2.).GT.0.001) PWT=MAX(1E-10,PA)**(PARU(41)-2.)    
      DO 110 J1=1,3 
      DO 110 J2=J1,3    
  110 SM(J1,J2)=SM(J1,J2)+PWT*P(I,J1)*P(I,J2)   
      PS=PS+PWT*PA**2   
  120 CONTINUE  
    
C...Very low multiplicities (0 or 1) not considered.    
      IF(NP.LE.1) THEN  
         CALL LUERRM_HIJING(8
     $        ,'(LUSPHE_HIJING:) too few particles for analysis')   
        SPH=-1. 
        APL=-1. 
        RETURN  
      ENDIF 
      DO 130 J1=1,3 
      DO 130 J2=J1,3    
  130 SM(J1,J2)=SM(J1,J2)/PS    
    
C...Find eigenvalues to matrix (third degree equation). 
      SQ=(SM(1,1)*SM(2,2)+SM(1,1)*SM(3,3)+SM(2,2)*SM(3,3)-SM(1,2)**2-   
     &SM(1,3)**2-SM(2,3)**2)/3.-1./9.   
      SR=-0.5*(SQ+1./9.+SM(1,1)*SM(2,3)**2+SM(2,2)*SM(1,3)**2+SM(3,3)*  
     &SM(1,2)**2-SM(1,1)*SM(2,2)*SM(3,3))+SM(1,2)*SM(1,3)*SM(2,3)+1./27.    
      SP=COS(ACOS(MAX(MIN(SR/SQRT(-SQ**3),1.),-1.))/3.) 
      P(N+1,4)=1./3.+SQRT(-SQ)*MAX(2.*SP,SQRT(3.*(1.-SP**2))-SP)    
      P(N+3,4)=1./3.+SQRT(-SQ)*MIN(2.*SP,-SQRT(3.*(1.-SP**2))-SP)   
      P(N+2,4)=1.-P(N+1,4)-P(N+3,4) 
      IF(P(N+2,4).LT.1E-5) THEN 
         CALL LUERRM_HIJING(8
     $        ,'(LUSPHE_HIJING:) all particles back-to-back')   
        SPH=-1. 
        APL=-1. 
        RETURN  
      ENDIF 
    
C...Find first and last eigenvector by solving equation system. 
      DO 170 I=1,3,2    
      DO 140 J1=1,3 
      SV(J1,J1)=SM(J1,J1)-P(N+I,4)  
      DO 140 J2=J1+1,3  
      SV(J1,J2)=SM(J1,J2)   
  140 SV(J2,J1)=SM(J1,J2)   
      SMAX=0.   
      DO 150 J1=1,3 
      DO 150 J2=1,3 
      IF(ABS(SV(J1,J2)).LE.SMAX) GOTO 150   
      JA=J1 
      JB=J2 
      SMAX=ABS(SV(J1,J2))   
  150 CONTINUE  
      SMAX=0.   
      DO 160 J3=JA+1,JA+2   
      J1=J3-3*((J3-1)/3)    
      RL=SV(J1,JB)/SV(JA,JB)    
      DO 160 J2=1,3 
      SV(J1,J2)=SV(J1,J2)-RL*SV(JA,J2)  
      IF(ABS(SV(J1,J2)).LE.SMAX) GOTO 160   
      JC=J1 
      SMAX=ABS(SV(J1,J2))   
  160 CONTINUE  
      JB1=JB+1-3*(JB/3) 
      JB2=JB+2-3*((JB+1)/3) 
      P(N+I,JB1)=-SV(JC,JB2)    
      P(N+I,JB2)=SV(JC,JB1) 
      P(N+I,JB)=-(SV(JA,JB1)*P(N+I,JB1)+SV(JA,JB2)*P(N+I,JB2))/ 
     &SV(JA,JB) 
      PA=SQRT(P(N+I,1)**2+P(N+I,2)**2+P(N+I,3)**2)  
      SGN=(-1.)**INT(RLU_HIJING(0)+0.5)    
      DO 170 J=1,3  
  170 P(N+I,J)=SGN*P(N+I,J)/PA  
    
C...Middle axis orthogonal to other two. Fill other codes.  
      SGN=(-1.)**INT(RLU_HIJING(0)+0.5)    
      P(N+2,1)=SGN*(P(N+1,2)*P(N+3,3)-P(N+1,3)*P(N+3,2))    
      P(N+2,2)=SGN*(P(N+1,3)*P(N+3,1)-P(N+1,1)*P(N+3,3))    
      P(N+2,3)=SGN*(P(N+1,1)*P(N+3,2)-P(N+1,2)*P(N+3,1))    
      DO 180 I=1,3  
      K(N+I,1)=31   
      K(N+I,2)=95   
      K(N+I,3)=I    
      K(N+I,4)=0    
      K(N+I,5)=0    
      P(N+I,5)=0.   
      DO 180 J=1,5  
  180 V(I,J)=0. 
    
C...Select storing option. Calculate sphericity and aplanarity. 
      MSTU(61)=N+1  
      MSTU(62)=NP   
      IF(MSTU(43).LE.1) MSTU(3)=3   
      IF(MSTU(43).GE.2) N=N+3   
      SPH=1.5*(P(N+2,4)+P(N+3,4))   
      APL=1.5*P(N+3,4)  
    
      RETURN    
      END   
