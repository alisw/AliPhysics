 
C********************************************************************* 
 
      SUBROUTINE LUNAME(KF,CHAU) 
 
C...Purpose: to give the particle/parton name as a character string. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT4/CHAF(500) 
      CHARACTER CHAF*8 
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT4/ 
      CHARACTER CHAU*16 
 
C...Initial values. Charge. Subdivide code. 
      CHAU=' ' 
      KFA=IABS(KF) 
      KC=LUCOMP(KF) 
      IF(KC.EQ.0) RETURN 
      KQ=LUCHGE(KF) 
      KFLA=MOD(KFA/1000,10) 
      KFLB=MOD(KFA/100,10) 
      KFLC=MOD(KFA/10,10) 
      KFLS=MOD(KFA,10) 
      KFLR=MOD(KFA/10000,10) 
 
C...Read out root name and spin for simple particle. 
      IF(KFA.LE.100.OR.(KFA.GT.100.AND.KC.GT.100)) THEN 
        CHAU=CHAF(KC) 
        LEN=0 
        DO 100 LEM=1,8 
        IF(CHAU(LEM:LEM).NE.' ') LEN=LEM 
  100   CONTINUE 
 
C...Construct root name for diquark. Add on spin. 
      ELSEIF(KFLC.EQ.0) THEN 
        CHAU(1:2)=CHAF(KFLA)(1:1)//CHAF(KFLB)(1:1) 
        IF(KFLS.EQ.1) CHAU(3:4)='_0' 
        IF(KFLS.EQ.3) CHAU(3:4)='_1' 
        LEN=4 
 
C...Construct root name for heavy meson. Add on spin and heavy flavour. 
      ELSEIF(KFLA.EQ.0) THEN 
        IF(KFLB.EQ.5) CHAU(1:1)='B' 
        IF(KFLB.EQ.6) CHAU(1:1)='T' 
        IF(KFLB.EQ.7) CHAU(1:1)='L' 
        IF(KFLB.EQ.8) CHAU(1:1)='H' 
        LEN=1 
        IF(KFLR.EQ.0.AND.KFLS.EQ.1) THEN 
        ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.3) THEN 
          CHAU(2:2)='*' 
          LEN=2 
        ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.3) THEN 
          CHAU(2:3)='_1' 
          LEN=3 
        ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.1) THEN 
          CHAU(2:4)='*_0' 
          LEN=4 
        ELSEIF(KFLR.EQ.2) THEN 
          CHAU(2:4)='*_1' 
          LEN=4 
        ELSEIF(KFLS.EQ.5) THEN 
          CHAU(2:4)='*_2' 
          LEN=4 
        ENDIF 
        IF(KFLC.GE.3.AND.KFLR.EQ.0.AND.KFLS.LE.3) THEN 
          CHAU(LEN+1:LEN+2)='_'//CHAF(KFLC)(1:1) 
          LEN=LEN+2 
        ELSEIF(KFLC.GE.3) THEN 
          CHAU(LEN+1:LEN+1)=CHAF(KFLC)(1:1) 
          LEN=LEN+1 
        ENDIF 
 
C...Construct root name and spin for heavy baryon. 
      ELSE 
        IF(KFLB.LE.2.AND.KFLC.LE.2) THEN 
          CHAU='Sigma ' 
          IF(KFLC.GT.KFLB) CHAU='Lambda' 
          IF(KFLS.EQ.4) CHAU='Sigma*' 
          LEN=5 
          IF(CHAU(6:6).NE.' ') LEN=6 
        ELSEIF(KFLB.LE.2.OR.KFLC.LE.2) THEN 
          CHAU='Xi ' 
          IF(KFLA.GT.KFLB.AND.KFLB.GT.KFLC) CHAU='Xi''' 
          IF(KFLS.EQ.4) CHAU='Xi*' 
          LEN=2 
          IF(CHAU(3:3).NE.' ') LEN=3 
        ELSE 
          CHAU='Omega ' 
          IF(KFLA.GT.KFLB.AND.KFLB.GT.KFLC) CHAU='Omega''' 
          IF(KFLS.EQ.4) CHAU='Omega*' 
          LEN=5 
          IF(CHAU(6:6).NE.' ') LEN=6 
        ENDIF 
 
C...Add on heavy flavour content for heavy baryon. 
        CHAU(LEN+1:LEN+2)='_'//CHAF(KFLA)(1:1) 
        LEN=LEN+2 
        IF(KFLB.GE.KFLC.AND.KFLC.GE.4) THEN 
          CHAU(LEN+1:LEN+2)=CHAF(KFLB)(1:1)//CHAF(KFLC)(1:1) 
          LEN=LEN+2 
        ELSEIF(KFLB.GE.KFLC.AND.KFLB.GE.4) THEN 
          CHAU(LEN+1:LEN+1)=CHAF(KFLB)(1:1) 
          LEN=LEN+1 
        ELSEIF(KFLC.GT.KFLB.AND.KFLB.GE.4) THEN 
          CHAU(LEN+1:LEN+2)=CHAF(KFLC)(1:1)//CHAF(KFLB)(1:1) 
          LEN=LEN+2 
        ELSEIF(KFLC.GT.KFLB.AND.KFLC.GE.4) THEN 
          CHAU(LEN+1:LEN+1)=CHAF(KFLC)(1:1) 
          LEN=LEN+1 
        ENDIF 
      ENDIF 
 
C...Add on bar sign for antiparticle (where necessary). 
      IF(KF.GT.0.OR.LEN.EQ.0) THEN 
      ELSEIF(KFA.GT.10.AND.KFA.LE.40.AND.KQ.NE.0.AND.MOD(KQ,3).EQ.0) 
     &THEN 
      ELSEIF(KFA.EQ.89.OR.(KFA.GE.91.AND.KFA.LE.99)) THEN 
      ELSEIF(KFA.GT.100.AND.KFLA.EQ.0.AND.KQ.NE.0) THEN 
      ELSEIF(MSTU(15).LE.1) THEN 
        CHAU(LEN+1:LEN+1)='~' 
        LEN=LEN+1 
      ELSE 
        CHAU(LEN+1:LEN+3)='bar' 
        LEN=LEN+3 
      ENDIF 
 
C...Add on charge where applicable (conventional cases skipped). 
      IF(KQ.EQ.6) CHAU(LEN+1:LEN+2)='++' 
      IF(KQ.EQ.-6) CHAU(LEN+1:LEN+2)='--' 
      IF(KQ.EQ.3) CHAU(LEN+1:LEN+1)='+' 
      IF(KQ.EQ.-3) CHAU(LEN+1:LEN+1)='-' 
      IF(KQ.EQ.0.AND.(KFA.LE.22.OR.LEN.EQ.0)) THEN 
      ELSEIF(KQ.EQ.0.AND.(KFA.GE.81.AND.KFA.LE.100)) THEN 
      ELSEIF(KFA.EQ.28.OR.KFA.EQ.29) THEN 
      ELSEIF(KFA.GT.100.AND.KFLA.EQ.0.AND.KFLB.EQ.KFLC.AND. 
     &KFLB.NE.1) THEN 
      ELSEIF(KQ.EQ.0) THEN 
        CHAU(LEN+1:LEN+1)='0' 
      ENDIF 
 
      RETURN 
      END 
