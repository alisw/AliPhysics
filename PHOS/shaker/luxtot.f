*CMZ :          17/07/98  15.44.35  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUXTOT(KFL,ECM,XTOT)

C...Purpose: to calculate total cross-section, including initial
C...state radiation effects.
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEND.

C...Status, (optimized) Q^2 scale, alpha_strong.
      PARJ(151)=ECM
      MSTJ(119)=10*MSTJ(102)+KFL
      IF(MSTJ(111).EQ.0) THEN
        Q2R=ECM**2
      ELSEIF(MSTU(111).EQ.0) THEN
        PARJ(168)=MIN(1.,MAX(PARJ(128),EXP(-12.*PARU(1)/
     &  ((33.-2.*MSTU(112))*PARU(111)))))
        Q2R=PARJ(168)*ECM**2
      ELSE
        PARJ(168)=MIN(1.,MAX(PARJ(128),PARU(112)/ECM,
     &  (2.*PARU(112)/ECM)**2))
        Q2R=PARJ(168)*ECM**2
      ENDIF
      ALSPI=ULALPS(Q2R)/PARU(1)

C...QCD corrections factor in R.
      IF(MSTJ(101).EQ.0.OR.MSTJ(109).EQ.1) THEN
        RQCD=1.
      ELSEIF(IABS(MSTJ(101)).EQ.1.AND.MSTJ(109).EQ.0) THEN
        RQCD=1.+ALSPI
      ELSEIF(MSTJ(109).EQ.0) THEN
        RQCD=1.+ALSPI+(1.986-0.115*MSTU(118))*ALSPI**2
        IF(MSTJ(111).EQ.1) RQCD=MAX(1.,RQCD+(33.-2.*MSTU(112))/12.*
     &  LOG(PARJ(168))*ALSPI**2)
      ELSEIF(IABS(MSTJ(101)).EQ.1) THEN
        RQCD=1.+(3./4.)*ALSPI
      ELSE
        RQCD=1.+(3./4.)*ALSPI-(3./32.+0.519*MSTU(118))*ALSPI**2
      ENDIF

C...Calculate Z0 width if default value not acceptable.
      IF(MSTJ(102).GE.3) THEN
        RVA=3.*(3.+(4.*PARU(102)-1.)**2)+6.*RQCD*(2.+(1.-8.*PARU(102)/
     &  3.)**2+(4.*PARU(102)/3.-1.)**2)
        DO 100 KFLC=5,6
        VQ=1.
        IF(MOD(MSTJ(103),2).EQ.1) VQ=SQRT(MAX(0.,1.-(2.*ULMASS(KFLC)/
     &  ECM)**2))
        IF(KFLC.EQ.5) VF=4.*PARU(102)/3.-1.
        IF(KFLC.EQ.6) VF=1.-8.*PARU(102)/3.
  100   RVA=RVA+3.*RQCD*(0.5*VQ*(3.-VQ**2)*VF**2+VQ**3)
        PARJ(124)=PARU(101)*PARJ(123)*RVA/(48.*PARU(102)*(1.-PARU(102)))
      ENDIF

C...Calculate propagator and related constants for QFD case.
      POLL=1.-PARJ(131)*PARJ(132)
      IF(MSTJ(102).GE.2) THEN
        SFF=1./(16.*PARU(102)*(1.-PARU(102)))
        SFW=ECM**4/((ECM**2-PARJ(123)**2)**2+(PARJ(123)*PARJ(124))**2)
        SFI=SFW*(1.-(PARJ(123)/ECM)**2)
        VE=4.*PARU(102)-1.
        SF1I=SFF*(VE*POLL+PARJ(132)-PARJ(131))
        SF1W=SFF**2*((VE**2+1.)*POLL+2.*VE*(PARJ(132)-PARJ(131)))
        HF1I=SFI*SF1I
        HF1W=SFW*SF1W
      ENDIF

C...Loop over different flavours: charge, velocity.
      RTOT=0.
      RQQ=0.
      RQV=0.
      RVA=0.
      DO 110 KFLC=1,MAX(MSTJ(104),KFL)
      IF(KFL.GT.0.AND.KFLC.NE.KFL) GOTO 110
      MSTJ(93)=1
      PMQ=ULMASS(KFLC)
      IF(ECM.LT.2.*PMQ+PARJ(127)) GOTO 110
      QF=KCHG(KFLC,1)/3.
      VQ=1.
      IF(MOD(MSTJ(103),2).EQ.1) VQ=SQRT(1.-(2.*PMQ/ECM)**2)

C...Calculate R and sum of charges for QED or QFD case.
      RQQ=RQQ+3.*QF**2*POLL
      IF(MSTJ(102).LE.1) THEN
        RTOT=RTOT+3.*0.5*VQ*(3.-VQ**2)*QF**2*POLL
      ELSE
        VF=SIGN(1.,QF)-4.*QF*PARU(102)
        RQV=RQV-6.*QF*VF*SF1I
        RVA=RVA+3.*(VF**2+1.)*SF1W
        RTOT=RTOT+3.*(0.5*VQ*(3.-VQ**2)*(QF**2*POLL-2.*QF*VF*HF1I+
     &  VF**2*HF1W)+VQ**3*HF1W)
      ENDIF
  110 CONTINUE
      RSUM=RQQ
      IF(MSTJ(102).GE.2) RSUM=RQQ+SFI*RQV+SFW*RVA

C...Calculate cross-section, including QCD corrections.
      PARJ(141)=RQQ
      PARJ(142)=RTOT
      PARJ(143)=RTOT*RQCD
      PARJ(144)=PARJ(143)
      PARJ(145)=PARJ(141)*86.8/ECM**2
      PARJ(146)=PARJ(142)*86.8/ECM**2
      PARJ(147)=PARJ(143)*86.8/ECM**2
      PARJ(148)=PARJ(147)
      PARJ(157)=RSUM*RQCD
      PARJ(158)=0.
      PARJ(159)=0.
      XTOT=PARJ(147)
      IF(MSTJ(107).LE.0) RETURN

C...Virtual cross-section.
      XKL=PARJ(135)
      XKU=MIN(PARJ(136),1.-(2.*PARJ(127)/ECM)**2)
      ALE=2.*LOG(ECM/ULMASS(11))-1.
      SIGV=ALE/3.+2.*LOG(ECM**2/(ULMASS(13)*ULMASS(15)))/3.-4./3.+
     &1.526*LOG(ECM**2/0.932)

C...Soft and hard radiative cross-section in QED case.
      IF(MSTJ(102).LE.1) THEN
        SIGV=1.5*ALE-0.5+PARU(1)**2/3.+2.*SIGV
        SIGS=ALE*(2.*LOG(XKL)-LOG(1.-XKL)-XKL)
        SIGH=ALE*(2.*LOG(XKU/XKL)-LOG((1.-XKU)/(1.-XKL))-(XKU-XKL))

C...Soft and hard radiative cross-section in QFD case.
      ELSE
        SZM=1.-(PARJ(123)/ECM)**2
        SZW=PARJ(123)*PARJ(124)/ECM**2
        PARJ(161)=-RQQ/RSUM
        PARJ(162)=-(RQQ+RQV+RVA)/RSUM
        PARJ(163)=(RQV*(1.-0.5*SZM-SFI)+RVA*(1.5-SZM-SFW))/RSUM
        PARJ(164)=(RQV*SZW**2*(1.-2.*SFW)+RVA*(2.*SFI+SZW**2-4.+3.*SZM-
     &  SZM**2))/(SZW*RSUM)
        SIGV=1.5*ALE-0.5+PARU(1)**2/3.+((2.*RQQ+SFI*RQV)/RSUM)*SIGV+
     &  (SZW*SFW*RQV/RSUM)*PARU(1)*20./9.
        SIGS=ALE*(2.*LOG(XKL)+PARJ(161)*LOG(1.-XKL)+PARJ(162)*XKL+
     &  PARJ(163)*LOG(((XKL-SZM)**2+SZW**2)/(SZM**2+SZW**2))+
     &  PARJ(164)*(ATAN((XKL-SZM)/SZW)-ATAN(-SZM/SZW)))
        SIGH=ALE*(2.*LOG(XKU/XKL)+PARJ(161)*LOG((1.-XKU)/(1.-XKL))+
     &  PARJ(162)*(XKU-XKL)+PARJ(163)*LOG(((XKU-SZM)**2+SZW**2)/
     &  ((XKL-SZM)**2+SZW**2))+PARJ(164)*(ATAN((XKU-SZM)/SZW)-
     &  ATAN((XKL-SZM)/SZW)))
      ENDIF

C...Total cross-section and fraction of hard photon events.
      PARJ(160)=SIGH/(PARU(1)/PARU(101)+SIGV+SIGS+SIGH)
      PARJ(157)=RSUM*(1.+(PARU(101)/PARU(1))*(SIGV+SIGS+SIGH))*RQCD
      PARJ(144)=PARJ(157)
      PARJ(148)=PARJ(144)*86.8/ECM**2
      XTOT=PARJ(148)

      RETURN
      END
