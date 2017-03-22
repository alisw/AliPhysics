* $Id$
    
C*********************************************************************  
    
      FUNCTION LUCHGE_HIJING(KF)   
    
C...Purpose: to give three times the charge for a particle/parton.  
#include "ludat2_hijing.inc"
    
C...Initial values. Simple case of direct readout.  
      LUCHGE_HIJING=0  
      KFA=IABS(KF)  
      KC=LUCOMP_HIJING(KFA)    
      IF(KC.EQ.0) THEN  
      ELSEIF(KFA.LE.100.OR.KC.LE.80.OR.KC.GT.100) THEN  
        LUCHGE_HIJING=KCHG(KC,1)   
    
C...Construction from quark content for heavy meson, diquark, baryon.   
      ELSEIF(MOD(KFA/1000,10).EQ.0) THEN    
        LUCHGE_HIJING=(KCHG(MOD(KFA/100,10),1)-KCHG(MOD(KFA/10,10),1))*    
     &  (-1)**MOD(KFA/100,10)   
      ELSEIF(MOD(KFA/10,10).EQ.0) THEN  
        LUCHGE_HIJING=KCHG(MOD(KFA/1000,10),1)+KCHG(MOD(KFA/100,10),1) 
      ELSE  
        LUCHGE_HIJING=KCHG(MOD(KFA/1000,10),1)+KCHG(MOD(KFA/100,10),1)+    
     &  KCHG(MOD(KFA/10,10),1)  
      ENDIF 
    
C...Add on correct sign.    
      LUCHGE_HIJING=LUCHGE_HIJING*ISIGN(1,KF) 
    
      RETURN    
      END   
