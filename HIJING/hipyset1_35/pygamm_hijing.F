* $Id$
    
C*********************************************************************  
    
      FUNCTION PYGAMM_HIJING(X)    
    
C...Gives ordinary Gamma function Gamma(x) for positive, real arguments;    
C...see M. Abramowitz, I. A. Stegun: Handbook of Mathematical Functions 
C...(Dover, 1965) 6.1.36.   
      DIMENSION B(8)    
      DATA B/-0.577191652,0.988205891,-0.897056937,0.918206857, 
     &-0.756704078,0.482199394,-0.193527818,0.035868343/    
    
      NX=INT(X) 
      DX=X-NX   
    
      PYGAMM_HIJING=1. 
      DO 100 I=1,8  
  100 PYGAMM_HIJING=PYGAMM_HIJING+B(I)*DX**I  
      IF(X.LT.1.) THEN  
        PYGAMM_HIJING=PYGAMM_HIJING/X 
      ELSE  
        DO 110 IX=1,NX-1    
  110   PYGAMM_HIJING=(X-IX)*PYGAMM_HIJING    
      ENDIF 
    
      RETURN    
      END   
