 
C********************************************************************* 
 
      SUBROUTINE LUXJET(ECM,NJET,CUT) 
 
C...Purpose: to select number of jets in matrix element approach. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      DIMENSION ZHUT(5) 
 
C...Relative three-jet rate in Zhu second order parametrization. 
      DATA ZHUT/3.0922, 6.2291, 7.4782, 7.8440, 8.2560/ 
 
C...Trivial result for two-jets only, including parton shower. 
      IF(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.5) THEN 
        CUT=0. 
 
C...QCD and Abelian vector gluon theory: Q^2 for jet rate and R. 
      ELSEIF(MSTJ(109).EQ.0.OR.MSTJ(109).EQ.2) THEN 
        CF=4./3. 
        IF(MSTJ(109).EQ.2) CF=1. 
        IF(MSTJ(111).EQ.0) THEN 
          Q2=ECM**2 
          Q2R=ECM**2 
        ELSEIF(MSTU(111).EQ.0) THEN 
          PARJ(169)=MIN(1.,PARJ(129)) 
          Q2=PARJ(169)*ECM**2 
          PARJ(168)=MIN(1.,MAX(PARJ(128),EXP(-12.*PARU(1)/ 
     &    ((33.-2.*MSTU(112))*PARU(111))))) 
          Q2R=PARJ(168)*ECM**2 
        ELSE 
          PARJ(169)=MIN(1.,MAX(PARJ(129),(2.*PARU(112)/ECM)**2)) 
          Q2=PARJ(169)*ECM**2 
          PARJ(168)=MIN(1.,MAX(PARJ(128),PARU(112)/ECM, 
     &    (2.*PARU(112)/ECM)**2)) 
          Q2R=PARJ(168)*ECM**2 
        ENDIF 
 
C...alpha_strong for R and R itself. 
        ALSPI=(3./4.)*CF*ULALPS(Q2R)/PARU(1) 
        IF(IABS(MSTJ(101)).EQ.1) THEN 
          RQCD=1.+ALSPI 
        ELSEIF(MSTJ(109).EQ.0) THEN 
          RQCD=1.+ALSPI+(1.986-0.115*MSTU(118))*ALSPI**2 
          IF(MSTJ(111).EQ.1) RQCD=MAX(1.,RQCD+(33.-2.*MSTU(112))/12.* 
     &    LOG(PARJ(168))*ALSPI**2) 
        ELSE 
          RQCD=1.+ALSPI-(3./32.+0.519*MSTU(118))*(4.*ALSPI/3.)**2 
        ENDIF 
 
C...alpha_strong for jet rate. Initial value for y cut. 
        ALSPI=(3./4.)*CF*ULALPS(Q2)/PARU(1) 
        CUT=MAX(0.001,PARJ(125),(PARJ(126)/ECM)**2) 
        IF(IABS(MSTJ(101)).LE.1.OR.(MSTJ(109).EQ.0.AND.MSTJ(111).EQ.0)) 
     &  CUT=MAX(CUT,EXP(-SQRT(0.75/ALSPI))/2.) 
        IF(MSTJ(110).EQ.2) CUT=MAX(0.01,MIN(0.05,CUT)) 
 
C...Parametrization of first order three-jet cross-section. 
  100   IF(MSTJ(101).EQ.0.OR.CUT.GE.0.25) THEN 
          PARJ(152)=0. 
        ELSE 
          PARJ(152)=(2.*ALSPI/3.)*((3.-6.*CUT+2.*LOG(CUT))* 
     &    LOG(CUT/(1.-2.*CUT))+(2.5+1.5*CUT-6.571)*(1.-3.*CUT)+ 
     &    5.833*(1.-3.*CUT)**2-3.894*(1.-3.*CUT)**3+ 
     &    1.342*(1.-3.*CUT)**4)/RQCD 
          IF(MSTJ(109).EQ.2.AND.(MSTJ(101).EQ.2.OR.MSTJ(101).LE.-2)) 
     &    PARJ(152)=0. 
        ENDIF 
 
C...Parametrization of second order three-jet cross-section. 
        IF(IABS(MSTJ(101)).LE.1.OR.MSTJ(101).EQ.3.OR.MSTJ(109).EQ.2.OR. 
     &  CUT.GE.0.25) THEN 
          PARJ(153)=0. 
        ELSEIF(MSTJ(110).LE.1) THEN 
          CT=LOG(1./CUT-2.) 
          PARJ(153)=ALSPI**2*CT**2*(2.419+0.5989*CT+0.6782*CT**2- 
     &    0.2661*CT**3+0.01159*CT**4)/RQCD 
 
C...Interpolation in second/first order ratio for Zhu parametrization. 
        ELSEIF(MSTJ(110).EQ.2) THEN 
          IZA=0 
          DO 110 IY=1,5 
          IF(ABS(CUT-0.01*IY).LT.0.0001) IZA=IY 
  110     CONTINUE 
          IF(IZA.NE.0) THEN 
            ZHURAT=ZHUT(IZA) 
          ELSE 
            IZ=100.*CUT 
            ZHURAT=ZHUT(IZ)+(100.*CUT-IZ)*(ZHUT(IZ+1)-ZHUT(IZ)) 
          ENDIF 
          PARJ(153)=ALSPI*PARJ(152)*ZHURAT 
        ENDIF 
 
C...Shift in second order three-jet cross-section with optimized Q^2. 
        IF(MSTJ(111).EQ.1.AND.IABS(MSTJ(101)).GE.2.AND.MSTJ(101).NE.3. 
     &  AND.CUT.LT.0.25) PARJ(153)=PARJ(153)+(33.-2.*MSTU(112))/12.* 
     &  LOG(PARJ(169))*ALSPI*PARJ(152) 
 
C...Parametrization of second order four-jet cross-section. 
        IF(IABS(MSTJ(101)).LE.1.OR.CUT.GE.0.125) THEN 
          PARJ(154)=0. 
        ELSE 
          CT=LOG(1./CUT-5.) 
          IF(CUT.LE.0.018) THEN 
            XQQGG=6.349-4.330*CT+0.8304*CT**2 
            IF(MSTJ(109).EQ.2) XQQGG=(4./3.)**2*(3.035-2.091*CT+ 
     &      0.4059*CT**2) 
            XQQQQ=1.25*(-0.1080+0.01486*CT+0.009364*CT**2) 
            IF(MSTJ(109).EQ.2) XQQQQ=8.*XQQQQ 
          ELSE 
            XQQGG=-0.09773+0.2959*CT-0.2764*CT**2+0.08832*CT**3 
            IF(MSTJ(109).EQ.2) XQQGG=(4./3.)**2*(-0.04079+0.1340*CT- 
     &      0.1326*CT**2+0.04365*CT**3) 
            XQQQQ=1.25*(0.003661-0.004888*CT-0.001081*CT**2+0.002093* 
     &      CT**3) 
            IF(MSTJ(109).EQ.2) XQQQQ=8.*XQQQQ 
          ENDIF 
          PARJ(154)=ALSPI**2*CT**2*(XQQGG+XQQQQ)/RQCD 
          PARJ(155)=XQQQQ/(XQQGG+XQQQQ) 
        ENDIF 
 
C...If negative three-jet rate, change y' optimization parameter. 
        IF(MSTJ(111).EQ.1.AND.PARJ(152)+PARJ(153).LT.0..AND. 
     &  PARJ(169).LT.0.99) THEN 
          PARJ(169)=MIN(1.,1.2*PARJ(169)) 
          Q2=PARJ(169)*ECM**2 
          ALSPI=(3./4.)*CF*ULALPS(Q2)/PARU(1) 
          GOTO 100 
        ENDIF 
 
C...If too high cross-section, use harder cuts, or fail. 
        IF(PARJ(152)+PARJ(153)+PARJ(154).GE.1) THEN 
          IF(MSTJ(110).EQ.2.AND.CUT.GT.0.0499.AND.MSTJ(111).EQ.1.AND. 
     &    PARJ(169).LT.0.99) THEN 
            PARJ(169)=MIN(1.,1.2*PARJ(169)) 
            Q2=PARJ(169)*ECM**2 
            ALSPI=(3./4.)*CF*ULALPS(Q2)/PARU(1) 
            GOTO 100 
          ELSEIF(MSTJ(110).EQ.2.AND.CUT.GT.0.0499) THEN 
            CALL LUERRM(26, 
     &      '(LUXJET:) no allowed y cut value for Zhu parametrization') 
          ENDIF 
          CUT=0.26*(4.*CUT)**(PARJ(152)+PARJ(153)+PARJ(154))**(-1./3.) 
          IF(MSTJ(110).EQ.2) CUT=MAX(0.01,MIN(0.05,CUT)) 
          GOTO 100 
        ENDIF 
 
C...Scalar gluon (first order only). 
      ELSE 
        ALSPI=ULALPS(ECM**2)/PARU(1) 
        CUT=MAX(0.001,PARJ(125),(PARJ(126)/ECM)**2,EXP(-3./ALSPI)) 
        PARJ(152)=0. 
        IF(CUT.LT.0.25) PARJ(152)=(ALSPI/3.)*((1.-2.*CUT)* 
     &  LOG((1.-2.*CUT)/CUT)+0.5*(9.*CUT**2-1.)) 
        PARJ(153)=0. 
        PARJ(154)=0. 
      ENDIF 
 
C...Select number of jets. 
      PARJ(150)=CUT 
      IF(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.5) THEN 
        NJET=2 
      ELSEIF(MSTJ(101).LE.0) THEN 
        NJET=MIN(4,2-MSTJ(101)) 
      ELSE 
        RNJ=RLU(0) 
        NJET=2 
        IF(PARJ(152)+PARJ(153)+PARJ(154).GT.RNJ) NJET=3 
        IF(PARJ(154).GT.RNJ) NJET=4 
      ENDIF 
 
      RETURN 
      END 
