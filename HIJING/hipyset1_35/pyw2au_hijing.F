* $Id$
    
C***********************************************************************    
    
      FUNCTION PYW2AU_HIJING(EPS,IREIM)    
    
C...Calculates real and imaginary parts of the auxiliary function W2;   
C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,   
C...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987    
#include "ludat1_hijing.inc"
    
      ASINH(X)=LOG(X+SQRT(X**2+1.)) 
      ACOSH(X)=LOG(X+SQRT(X**2-1.)) 
    
      IF(EPS.LT.0.) THEN    
        W2RE=4.*(ASINH(SQRT(-1./EPS)))**2   
        W2IM=0. 
      ELSEIF(EPS.LT.1.) THEN    
        W2RE=4.*(ACOSH(SQRT(1./EPS)))**2-PARU(1)**2 
        W2IM=-4.*PARU(1)*ACOSH(SQRT(1./EPS))    
      ELSE  
        W2RE=-4.*(ASIN(SQRT(1./EPS)))**2    
        W2IM=0. 
      ENDIF 
    
      IF(IREIM.EQ.1) PYW2AU_HIJING=W2RE    
      IF(IREIM.EQ.2) PYW2AU_HIJING=W2IM    
    
      RETURN    
      END   
