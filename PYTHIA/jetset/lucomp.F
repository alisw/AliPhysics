 
C********************************************************************* 
 
      FUNCTION LUCOMP(KF) 
 
C...Purpose: to compress the standard KF codes for use in mass and decay 
C...arrays; also to check whether a given code actually is defined. 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT2/ 
      DIMENSION KFTAB(25),KCTAB(25) 
      DATA KFTAB/211,111,221,311,321,130,310,213,113,223, 
     &313,323,2112,2212,210,2110,2210,110,220,330,440,30443,30553,0,0/ 
      DATA KCTAB/101,111,112,102,103,221,222,121,131,132, 
     &122,123,332,333,281,282,283,284,285,286,287,231,235,0,0/ 
 
C...Starting values. 
      LUCOMP=0 
      KFA=IABS(KF) 
 
C...Simple cases: direct translation or table. 
      IF(KFA.EQ.0.OR.KFA.GE.100000) THEN 
        RETURN 
      ELSEIF(KFA.LE.100) THEN 
        LUCOMP=KFA 
        IF(KF.LT.0.AND.KCHG(KFA,3).EQ.0) LUCOMP=0 
        RETURN 
      ELSE 
        DO 100 IKF=1,23 
        IF(KFA.EQ.KFTAB(IKF)) THEN 
          LUCOMP=KCTAB(IKF) 
          IF(KF.LT.0.AND.KCHG(LUCOMP,3).EQ.0) LUCOMP=0 
          RETURN 
        ENDIF 
  100   CONTINUE 
      ENDIF 
 
C...Subdivide KF code into constituent pieces. 
      KFLA=MOD(KFA/1000,10) 
      KFLB=MOD(KFA/100,10) 
      KFLC=MOD(KFA/10,10) 
      KFLS=MOD(KFA,10) 
      KFLR=MOD(KFA/10000,10) 
 
C...Mesons. 
      IF(KFA-10000*KFLR.LT.1000) THEN 
        IF(KFLB.EQ.0.OR.KFLB.EQ.9.OR.KFLC.EQ.0.OR.KFLC.EQ.9) THEN 
        ELSEIF(KFLB.LT.KFLC) THEN 
        ELSEIF(KF.LT.0.AND.KFLB.EQ.KFLC) THEN 
        ELSEIF(KFLB.EQ.KFLC) THEN 
          IF(KFLR.EQ.0.AND.KFLS.EQ.1) THEN 
            LUCOMP=110+KFLB 
          ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.3) THEN 
            LUCOMP=130+KFLB 
          ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.3) THEN 
            LUCOMP=150+KFLB 
          ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.1) THEN 
            LUCOMP=170+KFLB 
          ELSEIF(KFLR.EQ.2.AND.KFLS.EQ.3) THEN 
            LUCOMP=190+KFLB 
          ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.5) THEN 
            LUCOMP=210+KFLB 
          ENDIF 
        ELSEIF(KFLB.LE.5) THEN 
          IF(KFLR.EQ.0.AND.KFLS.EQ.1) THEN 
            LUCOMP=100+((KFLB-1)*(KFLB-2))/2+KFLC 
          ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.3) THEN 
            LUCOMP=120+((KFLB-1)*(KFLB-2))/2+KFLC 
          ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.3) THEN 
            LUCOMP=140+((KFLB-1)*(KFLB-2))/2+KFLC 
          ELSEIF(KFLR.EQ.1.AND.KFLS.EQ.1) THEN 
            LUCOMP=160+((KFLB-1)*(KFLB-2))/2+KFLC 
          ELSEIF(KFLR.EQ.2.AND.KFLS.EQ.3) THEN 
            LUCOMP=180+((KFLB-1)*(KFLB-2))/2+KFLC 
          ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.5) THEN 
            LUCOMP=200+((KFLB-1)*(KFLB-2))/2+KFLC 
          ENDIF 
        ELSEIF((KFLS.EQ.1.AND.KFLR.LE.1).OR.(KFLS.EQ.3.AND.KFLR.LE.2) 
     &  .OR.(KFLS.EQ.5.AND.KFLR.EQ.0)) THEN 
          LUCOMP=80+KFLB 
        ENDIF 
 
C...Diquarks. 
      ELSEIF((KFLR.EQ.0.OR.KFLR.EQ.1).AND.KFLC.EQ.0) THEN 
        IF(KFLS.NE.1.AND.KFLS.NE.3) THEN 
        ELSEIF(KFLA.EQ.9.OR.KFLB.EQ.0.OR.KFLB.EQ.9) THEN 
        ELSEIF(KFLA.LT.KFLB) THEN 
        ELSEIF(KFLS.EQ.1.AND.KFLA.EQ.KFLB) THEN 
        ELSE 
          LUCOMP=90 
        ENDIF 
 
C...Spin 1/2 baryons. 
      ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.2) THEN 
        IF(KFLA.EQ.9.OR.KFLB.EQ.0.OR.KFLB.EQ.9.OR.KFLC.EQ.9) THEN 
        ELSEIF(KFLA.LE.KFLC.OR.KFLA.LT.KFLB) THEN 
        ELSEIF(KFLA.GE.6.OR.KFLB.GE.4.OR.KFLC.GE.4) THEN 
          LUCOMP=80+KFLA 
        ELSEIF(KFLB.LT.KFLC) THEN 
          LUCOMP=300+((KFLA+1)*KFLA*(KFLA-1))/6+(KFLC*(KFLC-1))/2+KFLB 
        ELSE 
          LUCOMP=330+((KFLA+1)*KFLA*(KFLA-1))/6+(KFLB*(KFLB-1))/2+KFLC 
        ENDIF 
 
C...Spin 3/2 baryons. 
      ELSEIF(KFLR.EQ.0.AND.KFLS.EQ.4) THEN 
        IF(KFLA.EQ.9.OR.KFLB.EQ.0.OR.KFLB.EQ.9.OR.KFLC.EQ.9) THEN 
        ELSEIF(KFLA.LT.KFLB.OR.KFLB.LT.KFLC) THEN 
        ELSEIF(KFLA.GE.6.OR.KFLB.GE.4) THEN 
          LUCOMP=80+KFLA 
        ELSE 
          LUCOMP=360+((KFLA+1)*KFLA*(KFLA-1))/6+(KFLB*(KFLB-1))/2+KFLC 
        ENDIF 
      ENDIF 
 
      RETURN 
      END 
