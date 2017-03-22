* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYDIFF_HIJING 
    
C...Handles diffractive and elastic scattering. 
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
    
C...Reset K, P and V vectors. Store incoming particles. 
      DO 100 JT=1,MSTP(126)+10  
      I=MINT(83)+JT 
      DO 100 J=1,5  
      K(I,J)=0  
      P(I,J)=0. 
  100 V(I,J)=0. 
      N=MINT(84)    
      MINT(3)=0 
      MINT(21)=0    
      MINT(22)=0    
      MINT(23)=0    
      MINT(24)=0    
      MINT(4)=4 
      DO 110 JT=1,2 
      I=MINT(83)+JT 
      K(I,1)=21 
      K(I,2)=MINT(10+JT)    
      P(I,5)=VINT(2+JT) 
      P(I,3)=VINT(5)*(-1)**(JT+1)   
  110 P(I,4)=SQRT(P(I,3)**2+P(I,5)**2)  
      MINT(6)=2 
    
C...Subprocess; kinematics. 
      ISUB=MINT(1)  
      SQLAM=(VINT(2)-VINT(63)-VINT(64))**2-4.*VINT(63)*VINT(64) 
      PZ=SQRT(SQLAM)/(2.*VINT(1))   
      DO 150 JT=1,2 
      I=MINT(83)+JT 
      PE=(VINT(2)+VINT(62+JT)-VINT(65-JT))/(2.*VINT(1)) 
    
C...Elastically scattered particle. 
      IF(MINT(16+JT).LE.0) THEN 
        N=N+1   
        K(N,1)=1    
        K(N,2)=K(I,2)   
        K(N,3)=I+2  
        P(N,3)=PZ*(-1)**(JT+1)  
        P(N,4)=PE   
        P(N,5)=P(I,5)   
    
C...Diffracted particle: valence quark kicked out.  
      ELSEIF(MSTP(101).EQ.1) THEN   
        N=N+2   
        K(N-1,1)=2  
        K(N,1)=1    
        K(N-1,3)=I+2    
        K(N,3)=I+2  
        CALL PYSPLI_HIJING(K(I,2),21,K(N,2),K(N-1,2))  
        P(N-1,5)=ULMASS_HIJING(K(N-1,2))   
        P(N,5)=ULMASS_HIJING(K(N,2))   
        SQLAM=(VINT(62+JT)-P(N-1,5)**2-P(N,5)**2)**2-   
     &  4.*P(N-1,5)**2*P(N,5)**2    
        P(N-1,3)=(PE*SQRT(SQLAM)+PZ*(VINT(62+JT)+P(N-1,5)**2-   
     &  P(N,5)**2))/(2.*VINT(62+JT))*(-1)**(JT+1)   
        P(N-1,4)=SQRT(P(N-1,3)**2+P(N-1,5)**2)  
        P(N,3)=PZ*(-1)**(JT+1)-P(N-1,3) 
        P(N,4)=SQRT(P(N,3)**2+P(N,5)**2)    
    
C...Diffracted particle: gluon kicked out.  
      ELSE  
        N=N+3   
        K(N-2,1)=2  
        K(N-1,1)=2  
        K(N,1)=1    
        K(N-2,3)=I+2    
        K(N-1,3)=I+2    
        K(N,3)=I+2  
        CALL PYSPLI_HIJING(K(I,2),21,K(N,2),K(N-2,2))  
        K(N-1,2)=21 
        P(N-2,5)=ULMASS_HIJING(K(N-2,2))   
        P(N-1,5)=0. 
        P(N,5)=ULMASS_HIJING(K(N,2))   
C...Energy distribution for particle into two jets. 
  120   IMB=1   
        IF(MOD(K(I,2)/1000,10).NE.0) IMB=2  
        CHIK=PARP(92+2*IMB) 
        IF(MSTP(92).LE.1) THEN  
          IF(IMB.EQ.1) CHI=RLU_HIJING(0)   
          IF(IMB.EQ.2) CHI=1.-SQRT(RLU_HIJING(0))  
        ELSEIF(MSTP(92).EQ.2) THEN  
          CHI=1.-RLU_HIJING(0)**(1./(1.+CHIK)) 
        ELSEIF(MSTP(92).EQ.3) THEN  
          CUT=2.*0.3/VINT(1)    
  130     CHI=RLU_HIJING(0)**2 
          IF((CHI**2/(CHI**2+CUT**2))**0.25*(1.-CHI)**CHIK.LT.  
     &    RLU_HIJING(0)) GOTO 130  
        ELSE    
          CUT=2.*0.3/VINT(1)    
          CUTR=(1.+SQRT(1.+CUT**2))/CUT 
  140     CHIR=CUT*CUTR**RLU_HIJING(0) 
          CHI=(CHIR**2-CUT**2)/(2.*CHIR)    
          IF((1.-CHI)**CHIK.LT.RLU_HIJING(0)) GOTO 140 
        ENDIF   
        IF(CHI.LT.P(N,5)**2/VINT(62+JT).OR.CHI.GT.1.-P(N-2,5)**2/   
     &  VINT(62+JT)) GOTO 120   
        SQM=P(N-2,5)**2/(1.-CHI)+P(N,5)**2/CHI  
        IF((SQRT(SQM)+PARJ(32))**2.GE.VINT(62+JT)) GOTO 120 
        PZI=(PE*(VINT(62+JT)-SQM)+PZ*(VINT(62+JT)+SQM))/    
     &  (2.*VINT(62+JT))    
        PEI=SQRT(PZI**2+SQM)    
        PQQP=(1.-CHI)*(PEI+PZI) 
        P(N-2,3)=0.5*(PQQP-P(N-2,5)**2/PQQP)*(-1)**(JT+1)   
        P(N-2,4)=SQRT(P(N-2,3)**2+P(N-2,5)**2)  
        P(N-1,3)=(PZ-PZI)*(-1)**(JT+1)  
        P(N-1,4)=ABS(P(N-1,3))  
        P(N,3)=PZI*(-1)**(JT+1)-P(N-2,3)    
        P(N,4)=SQRT(P(N,3)**2+P(N,5)**2)    
      ENDIF 
    
C...Documentation lines.    
      K(I+2,1)=21   
      IF(MINT(16+JT).EQ.0) K(I+2,2)=MINT(10+JT) 
      IF(MINT(16+JT).NE.0) K(I+2,2)=10*(MINT(10+JT)/10) 
      K(I+2,3)=I    
      P(I+2,3)=PZ*(-1)**(JT+1)  
      P(I+2,4)=PE   
      P(I+2,5)=SQRT(VINT(62+JT))    
  150 CONTINUE  
    
C...Rotate outgoing partons/particles using cos(theta). 
      CALL LUDBRB_HIJING(MINT(83)+3,N,ACOS(VINT(23)),VINT(24),0D0,0D0
     $     ,0D0) 
    
      RETURN    
      END   
