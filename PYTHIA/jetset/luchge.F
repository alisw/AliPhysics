 
C********************************************************************* 
 
      FUNCTION LUCHGE(KF) 
 
C...Purpose: to give three times the charge for a particle/parton. 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT2/ 
 
C...Initial values. Simple case of direct readout. 
      LUCHGE=0 
      KFA=IABS(KF) 
      KC=LUCOMP(KFA) 
      IF(KC.EQ.0) THEN 
      ELSEIF(KFA.LE.100.OR.KC.LE.80.OR.KC.GT.100) THEN 
        LUCHGE=KCHG(KC,1) 
 
C...Construction from quark content for heavy meson, diquark, baryon. 
      ELSEIF(MOD(KFA/1000,10).EQ.0) THEN 
        LUCHGE=(KCHG(MOD(KFA/100,10),1)-KCHG(MOD(KFA/10,10),1))* 
     &  (-1)**MOD(KFA/100,10) 
      ELSEIF(MOD(KFA/10,10).EQ.0) THEN 
        LUCHGE=KCHG(MOD(KFA/1000,10),1)+KCHG(MOD(KFA/100,10),1) 
      ELSE 
        LUCHGE=KCHG(MOD(KFA/1000,10),1)+KCHG(MOD(KFA/100,10),1)+ 
     &  KCHG(MOD(KFA/10,10),1) 
      ENDIF 
 
C...Add on correct sign. 
      LUCHGE=LUCHGE*ISIGN(1,KF) 
 
      RETURN 
      END 
