* $Id$
    
C*********************************************************************  
    
      FUNCTION KLU_HIJING(I,J) 
    
C...Purpose: to provide various integer-valued event related data.  
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
    
C...Default value. For I=0 number of entries, number of stable entries  
C...or 3 times total charge.    
      KLU_HIJING=0 
      IF(I.LT.0.OR.I.GT.MSTU(4).OR.J.LE.0) THEN 
      ELSEIF(I.EQ.0.AND.J.EQ.1) THEN    
        KLU_HIJING=N   
      ELSEIF(I.EQ.0.AND.(J.EQ.2.OR.J.EQ.6)) THEN    
        DO 100 I1=1,N   
           IF(J.EQ.2.AND.K(I1,1).GE.1.AND.K(I1,1).LE.10) KLU_HIJING
     $          =KLU_HIJING+1 
           IF(J.EQ.6.AND.K(I1,1).GE.1.AND.K(I1,1).LE.10) KLU_HIJING
     $          =KLU_HIJING+LUCHGE_HIJING(K(I1,2)) 
  100   CONTINUE    
      ELSEIF(I.EQ.0) THEN   
    
C...For I > 0 direct readout of K matrix or charge. 
      ELSEIF(J.LE.5) THEN   
        KLU_HIJING=K(I,J)  
      ELSEIF(J.EQ.6) THEN   
        KLU_HIJING=LUCHGE_HIJING(K(I,2))  
    
C...Status (existing/fragmented/decayed), parton/hadron separation. 
      ELSEIF(J.LE.8) THEN   
        IF(K(I,1).GE.1.AND.K(I,1).LE.10) KLU_HIJING=1  
        IF(J.EQ.8) KLU_HIJING=KLU_HIJING*K(I,2)   
      ELSEIF(J.LE.12) THEN  
        KFA=IABS(K(I,2))    
        KC=LUCOMP_HIJING(KFA)  
        KQ=0    
        IF(KC.NE.0) KQ=KCHG(KC,2)   
        IF(J.EQ.9.AND.KC.NE.0.AND.KQ.NE.0) KLU_HIJING=K(I,2)   
        IF(J.EQ.10.AND.KC.NE.0.AND.KQ.EQ.0) KLU_HIJING=K(I,2)  
        IF(J.EQ.11) KLU_HIJING=KC  
        IF(J.EQ.12) KLU_HIJING=KQ*ISIGN(1,K(I,2))  
    
C...Heaviest flavour in hadron/diquark. 
      ELSEIF(J.EQ.13) THEN  
        KFA=IABS(K(I,2))    
        KLU_HIJING=MOD(KFA/100,10)*(-1)**MOD(KFA/100,10)   
        IF(KFA.LT.10) KLU_HIJING=KFA   
        IF(MOD(KFA/1000,10).NE.0) KLU_HIJING=MOD(KFA/1000,10)  
        KLU_HIJING=KLU_HIJING*ISIGN(1,K(I,2)) 
    
C...Particle history: generation, ancestor, rank.   
      ELSEIF(J.LE.16) THEN  
        I2=I    
        I1=I    
  110   KLU_HIJING=KLU_HIJING+1   
        I3=I2   
        I2=I1   
        I1=K(I1,3)  
        IF(I1.GT.0.AND.K(I1,1).GT.0.AND.K(I1,1).LE.20) GOTO 110 
        IF(J.EQ.15) KLU_HIJING=I2  
        IF(J.EQ.16) THEN    
          KLU_HIJING=0 
          DO 120 I1=I2+1,I3 
 120         IF(K(I1,3).EQ.I2.AND.K(I1,1).GT.0.AND.K(I1,1).LE.20)
     $            KLU_HIJING=KLU_HIJING+1    
        ENDIF   
    
C...Particle coming from collapsing jet system or not.  
      ELSEIF(J.EQ.17) THEN  
        I1=I    
  130   KLU_HIJING=KLU_HIJING+1   
        I3=I1   
        I1=K(I1,3)  
        I0=MAX(1,I1)    
        KC=LUCOMP_HIJING(K(I0,2))  
        IF(I1.EQ.0.OR.K(I0,1).LE.0.OR.K(I0,1).GT.20.OR.KC.EQ.0) THEN    
          IF(KLU_HIJING.EQ.1) KLU_HIJING=-1   
          IF(KLU_HIJING.GT.1) KLU_HIJING=0    
          RETURN    
        ENDIF   
        IF(KCHG(KC,2).EQ.0) GOTO 130    
        IF(K(I1,1).NE.12) KLU_HIJING=0 
        IF(K(I1,1).NE.12) RETURN    
        I2=I1   
  140   I2=I2+1 
        IF(I2.LT.N.AND.K(I2,1).NE.11) GOTO 140  
        K3M=K(I3-1,3)   
        IF(K3M.GE.I1.AND.K3M.LE.I2) KLU_HIJING=0   
        K3P=K(I3+1,3)   
        IF(I3.LT.N.AND.K3P.GE.I1.AND.K3P.LE.I2) KLU_HIJING=0   
    
C...Number of decay products. Colour flow.  
      ELSEIF(J.EQ.18) THEN  
         IF(K(I,1).EQ.11.OR.K(I,1).EQ.12) KLU_HIJING=MAX(0,K(I,5)-K(I,4)
     $        +1) 
        IF(K(I,4).EQ.0.OR.K(I,5).EQ.0) KLU_HIJING=0    
      ELSEIF(J.LE.22) THEN  
        IF(K(I,1).NE.3.AND.K(I,1).NE.13.AND.K(I,1).NE.14) RETURN    
        IF(J.EQ.19) KLU_HIJING=MOD(K(I,4)/MSTU(5),MSTU(5)) 
        IF(J.EQ.20) KLU_HIJING=MOD(K(I,5)/MSTU(5),MSTU(5)) 
        IF(J.EQ.21) KLU_HIJING=MOD(K(I,4),MSTU(5)) 
        IF(J.EQ.22) KLU_HIJING=MOD(K(I,5),MSTU(5)) 
      ELSE  
      ENDIF 
    
      RETURN    
      END   
