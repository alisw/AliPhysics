* $Id$
    
C*********************************************************************  
    
      FUNCTION ULANGL_HIJING(X,Y)  
    
C...Purpose: to reconstruct an angle from given x and y coordinates.    
#include "ludat1_hijing.inc"
    
      ULANGL_HIJING=0. 
      R=SQRT(X**2+Y**2) 
      IF(R.LT.1E-20) RETURN 
      IF(ABS(X)/R.LT.0.8) THEN  
        ULANGL_HIJING=SIGN(ACOS(X/R),Y)    
      ELSE  
        ULANGL_HIJING=ASIN(Y/R)    
        IF(X.LT.0..AND.ULANGL_HIJING.GE.0.) THEN   
          ULANGL_HIJING=PARU(1)-ULANGL_HIJING 
        ELSEIF(X.LT.0.) THEN    
          ULANGL_HIJING=-PARU(1)-ULANGL_HIJING    
        ENDIF   
      ENDIF 
    
      RETURN    
      END   
