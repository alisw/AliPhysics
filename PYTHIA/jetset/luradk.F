 
C********************************************************************* 
 
      SUBROUTINE LURADK(ECM,MK,PAK,THEK,PHIK,ALPK) 
 
C...Purpose: to generate initial state photon radiation. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
 
C...Function: cumulative hard photon spectrum in QFD case. 
      FXK(XX)=2.*LOG(XX)+PARJ(161)*LOG(1.-XX)+PARJ(162)*XX+ 
     &PARJ(163)*LOG((XX-SZM)**2+SZW**2)+PARJ(164)*ATAN((XX-SZM)/SZW) 
 
C...Determine whether radiative photon or not. 
      MK=0 
      PAK=0. 
      IF(PARJ(160).LT.RLU(0)) RETURN 
      MK=1 
 
C...Photon energy range. Find photon momentum in QED case. 
      XKL=PARJ(135) 
      XKU=MIN(PARJ(136),1.-(2.*PARJ(127)/ECM)**2) 
      IF(MSTJ(102).LE.1) THEN 
  100   XK=1./(1.+(1./XKL-1.)*((1./XKU-1.)/(1./XKL-1.))**RLU(0)) 
        IF(1.+(1.-XK)**2.LT.2.*RLU(0)) GOTO 100 
 
C...Ditto in QFD case, by numerical inversion of integrated spectrum. 
      ELSE 
        SZM=1.-(PARJ(123)/ECM)**2 
        SZW=PARJ(123)*PARJ(124)/ECM**2 
        FXKL=FXK(XKL) 
        FXKU=FXK(XKU) 
        FXKD=1E-4*(FXKU-FXKL) 
        FXKR=FXKL+RLU(0)*(FXKU-FXKL) 
        NXK=0 
  110   NXK=NXK+1 
        XK=0.5*(XKL+XKU) 
        FXKV=FXK(XK) 
        IF(FXKV.GT.FXKR) THEN 
          XKU=XK 
          FXKU=FXKV 
        ELSE 
          XKL=XK 
          FXKL=FXKV 
        ENDIF 
        IF(NXK.LT.15.AND.FXKU-FXKL.GT.FXKD) GOTO 110 
        XK=XKL+(XKU-XKL)*(FXKR-FXKL)/(FXKU-FXKL) 
      ENDIF 
      PAK=0.5*ECM*XK 
 
C...Photon polar and azimuthal angle. 
      PME=2.*(ULMASS(11)/ECM)**2 
  120 CTHM=PME*(2./PME)**RLU(0) 
      IF(1.-(XK**2*CTHM*(1.-0.5*CTHM)+2.*(1.-XK)*PME/MAX(PME, 
     &CTHM*(1.-0.5*CTHM)))/(1.+(1.-XK)**2).LT.RLU(0)) GOTO 120 
      CTHE=1.-CTHM 
      IF(RLU(0).GT.0.5) CTHE=-CTHE 
      STHE=SQRT(MAX(0.,(CTHM-PME)*(2.-CTHM))) 
      THEK=ULANGL(CTHE,STHE) 
      PHIK=PARU(2)*RLU(0) 
 
C...Rotation angle for hadronic system. 
      SGN=1. 
      IF(0.5*(2.-XK*(1.-CTHE))**2/((2.-XK)**2+(XK*CTHE)**2).GT. 
     &RLU(0)) SGN=-1. 
      ALPK=ASIN(SGN*STHE*(XK-SGN*(2.*SQRT(1.-XK)-2.+XK)*CTHE)/ 
     &(2.-XK*(1.-SGN*CTHE))) 
 
      RETURN 
      END 
