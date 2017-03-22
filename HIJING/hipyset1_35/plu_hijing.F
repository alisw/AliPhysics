* $Id$
    
C*********************************************************************  
    
      FUNCTION PLU_HIJING(I,J) 
    
C...Purpose: to provide various real-valued event related data. 
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
      DIMENSION PSUM(4) 
    
C...Set default value. For I = 0 sum of momenta or charges, 
C...or invariant mass of system.    
      PLU_HIJING=0.    
      IF(I.LT.0.OR.I.GT.MSTU(4).OR.J.LE.0) THEN 
      ELSEIF(I.EQ.0.AND.J.LE.4) THEN    
        DO 100 I1=1,N   
  100   IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PLU_HIJING=PLU_HIJING+P(I1,J)  
      ELSEIF(I.EQ.0.AND.J.EQ.5) THEN    
        DO 110 J1=1,4   
        PSUM(J1)=0. 
        DO 110 I1=1,N   
  110   IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PSUM(J1)=PSUM(J1)+P(I1,J1)   
        PLU_HIJING=SQRT(MAX(0.,PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)
     $       **2))   
      ELSEIF(I.EQ.0.AND.J.EQ.6) THEN    
        DO 120 I1=1,N   
 120       IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PLU_HIJING=PLU_HIJING
     $          +LUCHGE_HIJING(K(I1,2))/3.   
      ELSEIF(I.EQ.0) THEN   
    
C...Direct readout of P matrix. 
      ELSEIF(J.LE.5) THEN   
        PLU_HIJING=P(I,J)  
    
C...Charge, total momentum, transverse momentum, transverse mass.   
      ELSEIF(J.LE.12) THEN  
        IF(J.EQ.6) PLU_HIJING=LUCHGE_HIJING(K(I,2))/3.    
        IF(J.EQ.7.OR.J.EQ.8) PLU_HIJING=P(I,1)**2+P(I,2)**2+P(I,3)**2  
        IF(J.EQ.9.OR.J.EQ.10) PLU_HIJING=P(I,1)**2+P(I,2)**2   
        IF(J.EQ.11.OR.J.EQ.12) PLU_HIJING=P(I,5)**2+P(I,1)**2+P(I,2)**2    
        IF(J.EQ.8.OR.J.EQ.10.OR.J.EQ.12) PLU_HIJING=SQRT(PLU_HIJING)  
    
C...Theta and phi angle in radians or degrees.  
      ELSEIF(J.LE.16) THEN  
         IF(J.LE.14) PLU_HIJING=ULANGL_HIJING(P(I,3),SQRT(P(I,1)**2+P(I
     $        ,2)**2))    
        IF(J.GE.15) PLU_HIJING=ULANGL_HIJING(P(I,1),P(I,2))   
        IF(J.EQ.14.OR.J.EQ.16) PLU_HIJING=PLU_HIJING*180./PARU(1) 
    
C...True rapidity, rapidity with pion mass, pseudorapidity. 
      ELSEIF(J.LE.19) THEN  
        PMR=0.  
        IF(J.EQ.17) PMR=P(I,5)  
        IF(J.EQ.18) PMR=ULMASS_HIJING(211) 
        PR=MAX(1E-20,PMR**2+P(I,1)**2+P(I,2)**2)    
        PLU_HIJING=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/SQRT(PR
     $       ),1E20)),P(I,3))  
    
C...Energy and momentum fractions (only to be used in CM frame).    
      ELSEIF(J.LE.25) THEN  
         IF(J.EQ.20) PLU_HIJING=2.*SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
     $        /PARU(21) 
        IF(J.EQ.21) PLU_HIJING=2.*P(I,3)/PARU(21)  
        IF(J.EQ.22) PLU_HIJING=2.*SQRT(P(I,1)**2+P(I,2)**2)/PARU(21)   
        IF(J.EQ.23) PLU_HIJING=2.*P(I,4)/PARU(21)  
        IF(J.EQ.24) PLU_HIJING=(P(I,4)+P(I,3))/PARU(21)    
        IF(J.EQ.25) PLU_HIJING=(P(I,4)-P(I,3))/PARU(21)    
      ENDIF 
    
      RETURN    
      END   
