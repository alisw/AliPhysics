 
C********************************************************************* 
 
      SUBROUTINE LUZDIS(KFL1,KFL2,PR,Z) 
 
C...Purpose: to generate the longitudinal splitting variable z. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT1/,/LUDAT2/ 
 
C...Check if heavy flavour fragmentation. 
      KFLA=IABS(KFL1) 
      KFLB=IABS(KFL2) 
      KFLH=KFLA 
      IF(KFLA.GE.10) KFLH=MOD(KFLA/1000,10) 
 
C...Lund symmetric scaling function: determine parameters of shape. 
      IF(MSTJ(11).EQ.1.OR.(MSTJ(11).EQ.3.AND.KFLH.LE.3).OR. 
     &MSTJ(11).GE.4) THEN 
        FA=PARJ(41) 
        IF(MSTJ(91).EQ.1) FA=PARJ(43) 
        IF(KFLB.GE.10) FA=FA+PARJ(45) 
        FBB=PARJ(42) 
        IF(MSTJ(91).EQ.1) FBB=PARJ(44) 
        FB=FBB*PR 
        FC=1. 
        IF(KFLA.GE.10) FC=FC-PARJ(45) 
        IF(KFLB.GE.10) FC=FC+PARJ(45) 
        IF(MSTJ(11).GE.4.AND.KFLH.GE.4.AND.KFLH.LE.5) THEN 
          FRED=PARJ(46) 
          IF(MSTJ(11).EQ.5.AND.KFLH.EQ.5) FRED=PARJ(47) 
          FC=FC+FRED*FBB*PARF(100+KFLH)**2 
        ELSEIF(MSTJ(11).GE.4.AND.KFLH.GE.6.AND.KFLH.LE.8) THEN 
          FRED=PARJ(46) 
          IF(MSTJ(11).EQ.5) FRED=PARJ(48) 
          FC=FC+FRED*FBB*PMAS(KFLH,1)**2 
        ENDIF 
        MC=1 
        IF(ABS(FC-1.).GT.0.01) MC=2 
 
C...Determine position of maximum. Special cases for a = 0 or a = c. 
        IF(FA.LT.0.02) THEN 
          MA=1 
          ZMAX=1. 
          IF(FC.GT.FB) ZMAX=FB/FC 
        ELSEIF(ABS(FC-FA).LT.0.01) THEN 
          MA=2 
          ZMAX=FB/(FB+FC) 
        ELSE 
          MA=3 
          ZMAX=0.5*(FB+FC-SQRT((FB-FC)**2+4.*FA*FB))/(FC-FA) 
          IF(ZMAX.GT.0.9999.AND.FB.GT.100.) ZMAX=MIN(ZMAX,1.-FA/FB) 
        ENDIF 
 
C...Subdivide z range if distribution very peaked near endpoint. 
        MMAX=2 
        IF(ZMAX.LT.0.1) THEN 
          MMAX=1 
          ZDIV=2.75*ZMAX 
          IF(MC.EQ.1) THEN 
            FINT=1.-LOG(ZDIV) 
          ELSE 
            ZDIVC=ZDIV**(1.-FC) 
            FINT=1.+(1.-1./ZDIVC)/(FC-1.) 
          ENDIF 
        ELSEIF(ZMAX.GT.0.85.AND.FB.GT.1.) THEN 
          MMAX=3 
          FSCB=SQRT(4.+(FC/FB)**2) 
          ZDIV=FSCB-1./ZMAX-(FC/FB)*LOG(ZMAX*0.5*(FSCB+FC/FB)) 
          IF(MA.GE.2) ZDIV=ZDIV+(FA/FB)*LOG(1.-ZMAX) 
          ZDIV=MIN(ZMAX,MAX(0.,ZDIV)) 
          FINT=1.+FB*(1.-ZDIV) 
        ENDIF 
 
C...Choice of z, preweighted for peaks at low or high z. 
  100   Z=RLU(0) 
        FPRE=1. 
        IF(MMAX.EQ.1) THEN 
          IF(FINT*RLU(0).LE.1.) THEN 
            Z=ZDIV*Z 
          ELSEIF(MC.EQ.1) THEN 
            Z=ZDIV**Z 
            FPRE=ZDIV/Z 
          ELSE 
            Z=1./(ZDIVC+Z*(1.-ZDIVC))**(1./(1.-FC)) 
            FPRE=(ZDIV/Z)**FC 
          ENDIF 
        ELSEIF(MMAX.EQ.3) THEN 
          IF(FINT*RLU(0).LE.1.) THEN 
            Z=ZDIV+LOG(Z)/FB 
            FPRE=EXP(FB*(Z-ZDIV)) 
          ELSE 
            Z=ZDIV+Z*(1.-ZDIV) 
          ENDIF 
        ENDIF 
 
C...Weighting according to correct formula. 
        IF(Z.LE.0..OR.Z.GE.1.) GOTO 100 
        FEXP=FC*LOG(ZMAX/Z)+FB*(1./ZMAX-1./Z) 
        IF(MA.GE.2) FEXP=FEXP+FA*LOG((1.-Z)/(1.-ZMAX)) 
        FVAL=EXP(MAX(-50.,MIN(50.,FEXP))) 
        IF(FVAL.LT.RLU(0)*FPRE) GOTO 100 
 
C...Generate z according to Field-Feynman, SLAC, (1-z)**c OR z**c. 
      ELSE 
        FC=PARJ(50+MAX(1,KFLH)) 
        IF(MSTJ(91).EQ.1) FC=PARJ(59) 
  110   Z=RLU(0) 
        IF(FC.GE.0..AND.FC.LE.1.) THEN 
          IF(FC.GT.RLU(0)) Z=1.-Z**(1./3.) 
        ELSEIF(FC.GT.-1.AND.FC.LT.0.) THEN 
          IF(-4.*FC*Z*(1.-Z)**2.LT.RLU(0)*((1.-Z)**2-FC*Z)**2) GOTO 110 
        ELSE 
          IF(FC.GT.0.) Z=1.-Z**(1./FC) 
          IF(FC.LT.0.) Z=Z**(-1./FC) 
        ENDIF 
      ENDIF 
 
      RETURN 
      END 
