*CMZ :          17/07/98  15.44.35  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUXKFL(KFL,ECM,ECMC,KFLC)

C...Purpose: to select flavour for produced qqbar pair.
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEND.

C...Calculate maximum weight in QED or QFD case.
      IF(MSTJ(102).LE.1) THEN
        RFMAX=4./9.
      ELSE
        POLL=1.-PARJ(131)*PARJ(132)
        SFF=1./(16.*PARU(102)*(1.-PARU(102)))
        SFW=ECMC**4/((ECMC**2-PARJ(123)**2)**2+(PARJ(123)*PARJ(124))**2)
        SFI=SFW*(1.-(PARJ(123)/ECMC)**2)
        VE=4.*PARU(102)-1.
        HF1I=SFI*SFF*(VE*POLL+PARJ(132)-PARJ(131))
        HF1W=SFW*SFF**2*((VE**2+1.)*POLL+2.*VE*(PARJ(132)-PARJ(131)))
        RFMAX=MAX(4./9.*POLL-4./3.*(1.-8.*PARU(102)/3.)*HF1I+
     &  ((1.-8.*PARU(102)/3.)**2+1.)*HF1W,1./9.*POLL+2./3.*
     &  (-1.+4.*PARU(102)/3.)*HF1I+((-1.+4.*PARU(102)/3.)**2+1.)*HF1W)
      ENDIF

C...Choose flavour. Gives charge and velocity.
      NTRY=0
  100 NTRY=NTRY+1
      IF(NTRY.GT.100) THEN
        CALL LUERRM(14,'(LUXKFL:) caught in an infinite loop')
        KFLC=0
        RETURN
      ENDIF
      KFLC=KFL
      IF(KFL.LE.0) KFLC=1+INT(MSTJ(104)*RLU(0))
      MSTJ(93)=1
      PMQ=ULMASS(KFLC)
      IF(ECM.LT.2.*PMQ+PARJ(127)) GOTO 100
      QF=KCHG(KFLC,1)/3.
      VQ=1.
      IF(MOD(MSTJ(103),2).EQ.1) VQ=SQRT(MAX(0.,1.-(2.*PMQ/ECMC)**2))

C...Calculate weight in QED or QFD case.
      IF(MSTJ(102).LE.1) THEN
        RF=QF**2
        RFV=0.5*VQ*(3.-VQ**2)*QF**2
      ELSE
        VF=SIGN(1.,QF)-4.*QF*PARU(102)
        RF=QF**2*POLL-2.*QF*VF*HF1I+(VF**2+1.)*HF1W
        RFV=0.5*VQ*(3.-VQ**2)*(QF**2*POLL-2.*QF*VF*HF1I+VF**2*HF1W)+
     &  VQ**3*HF1W
      ENDIF

C...Weighting or new event (radiative photon). Cross-section update.
      IF(KFL.LE.0.AND.RF.LT.RLU(0)*RFMAX) GOTO 100
      PARJ(158)=PARJ(158)+1.
      IF(ECMC.LT.2.*PMQ+PARJ(127).OR.RFV.LT.RLU(0)*RF) KFLC=0
      IF(MSTJ(107).LE.0.AND.KFLC.EQ.0) GOTO 100
      IF(KFLC.NE.0) PARJ(159)=PARJ(159)+1.
      PARJ(144)=PARJ(157)*PARJ(159)/PARJ(158)
      PARJ(148)=PARJ(144)*86.8/ECM**2

      RETURN
      END
