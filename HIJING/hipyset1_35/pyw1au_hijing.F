* $Id$
    
C***********************************************************************    
    
      FUNCTION PYW1AU_HIJING(EPS,IREIM)    
    
C...Calculates real and imaginary parts of the auxiliary function W1;   
C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,   
C...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987    
#include "ludat1_hijing.inc"
    
      ASINH(X)=LOG(X+SQRT(X**2+1.)) 
      ACOSH(X)=LOG(X+SQRT(X**2-1.)) 
    
      IF(EPS.LT.0.) THEN    
        W1RE=2.*SQRT(1.-EPS)*ASINH(SQRT(-1./EPS))   
        W1IM=0. 
      ELSEIF(EPS.LT.1.) THEN    
        W1RE=2.*SQRT(1.-EPS)*ACOSH(SQRT(1./EPS))    
        W1IM=-PARU(1)*SQRT(1.-EPS)  
      ELSE  
        W1RE=2.*SQRT(EPS-1.)*ASIN(SQRT(1./EPS)) 
        W1IM=0. 
      ENDIF 
    
      IF(IREIM.EQ.1) PYW1AU_HIJING=W1RE    
      IF(IREIM.EQ.2) PYW1AU_HIJING=W1IM    
    
      RETURN    
      END   
