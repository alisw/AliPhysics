* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYFRAM_HIJING(IFRAME) 
    
C...Performs transformations between different coordinate frames.   
#include "ludat1_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
    
      IF(IFRAME.LT.1.OR.IFRAME.GT.2) THEN   
        WRITE(MSTU(11),1000) IFRAME,MINT(6) 
        RETURN  
      ENDIF 
      IF(IFRAME.EQ.MINT(6)) RETURN  
    
      IF(MINT(6).EQ.1) THEN 
C...Transform from fixed target or user specified frame to  
C...CM-frame of incoming particles. 
        CALL LUROBO_HIJING(0.,0.,-VINT(8),-VINT(9),-VINT(10))  
        CALL LUROBO_HIJING(0.,-VINT(7),0.,0.,0.)   
        CALL LUROBO_HIJING(-VINT(6),0.,0.,0.,0.)   
        MINT(6)=2   
    
      ELSE  
C...Transform from particle CM-frame to fixed target or user specified  
C...frame.  
        CALL LUROBO_HIJING(VINT(6),VINT(7),VINT(8),VINT(9),VINT(10))   
        MINT(6)=1   
      ENDIF 
      MSTI(6)=MINT(6)   
    
 1000 FORMAT(1X,'Error: illegal values in subroutine PYFRAM_HIJING.',1X,   
     &'No transformation performed.'/1X,'IFRAME =',1X,I5,'; MINT(6) =', 
     &1X,I5)    
    
      RETURN    
      END   
