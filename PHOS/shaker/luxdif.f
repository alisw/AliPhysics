*CMZ :          17/07/98  15.44.35  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUXDIF(NC,NJET,KFL,ECM,CHI,THE,PHI)

C...Purpose: to give the angular orientation of events.
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEND.

C...Charge. Factors depending on polarization for QED case.
      QF=KCHG(KFL,1)/3.
      POLL=1.-PARJ(131)*PARJ(132)
      POLD=PARJ(132)-PARJ(131)
      IF(MSTJ(102).LE.1.OR.MSTJ(109).EQ.1) THEN
        HF1=POLL
        HF2=0.
        HF3=PARJ(133)**2
        HF4=0.

C...Factors depending on flavour, energy and polarization for QFD case.
      ELSE
        SFF=1./(16.*PARU(102)*(1.-PARU(102)))
        SFW=ECM**4/((ECM**2-PARJ(123)**2)**2+(PARJ(123)*PARJ(124))**2)
        SFI=SFW*(1.-(PARJ(123)/ECM)**2)
        AE=-1.
        VE=4.*PARU(102)-1.
        AF=SIGN(1.,QF)
        VF=AF-4.*QF*PARU(102)
        HF1=QF**2*POLL-2.*QF*VF*SFI*SFF*(VE*POLL-AE*POLD)+
     &  (VF**2+AF**2)*SFW*SFF**2*((VE**2+AE**2)*POLL-2.*VE*AE*POLD)
        HF2=-2.*QF*AF*SFI*SFF*(AE*POLL-VE*POLD)+2.*VF*AF*SFW*SFF**2*
     &  (2.*VE*AE*POLL-(VE**2+AE**2)*POLD)
        HF3=PARJ(133)**2*(QF**2-2.*QF*VF*SFI*SFF*VE+(VF**2+AF**2)*
     &  SFW*SFF**2*(VE**2-AE**2))
        HF4=-PARJ(133)**2*2.*QF*VF*SFW*(PARJ(123)*PARJ(124)/ECM**2)*
     &  SFF*AE
      ENDIF

C...Mass factor. Differential cross-sections for two-jet events.
      SQ2=SQRT(2.)
      QME=0.
      IF(MSTJ(103).GE.4.AND.IABS(MSTJ(101)).LE.1.AND.MSTJ(102).LE.1.AND.
     &MSTJ(109).NE.1) QME=(2.*ULMASS(KFL)/ECM)**2
      IF(NJET.EQ.2) THEN
        SIGU=4.*SQRT(1.-QME)
        SIGL=2.*QME*SQRT(1.-QME)
        SIGT=0.
        SIGI=0.
        SIGA=0.
        SIGP=4.

C...Kinematical variables. Reduce four-jet event to three-jet one.
      ELSE
        IF(NJET.EQ.3) THEN
          X1=2.*P(NC+1,4)/ECM
          X2=2.*P(NC+3,4)/ECM
        ELSE
          ECMR=P(NC+1,4)+P(NC+4,4)+SQRT((P(NC+2,1)+P(NC+3,1))**2+
     &    (P(NC+2,2)+P(NC+3,2))**2+(P(NC+2,3)+P(NC+3,3))**2)
          X1=2.*P(NC+1,4)/ECMR
          X2=2.*P(NC+4,4)/ECMR
        ENDIF

C...Differential cross-sections for three-jet (or reduced four-jet).
        XQ=(1.-X1)/(1.-X2)
        CT12=(X1*X2-2.*X1-2.*X2+2.+QME)/SQRT((X1**2-QME)*(X2**2-QME))
        ST12=SQRT(1.-CT12**2)
        IF(MSTJ(109).NE.1) THEN
          SIGU=2.*X1**2+X2**2*(1.+CT12**2)-QME*(3.+CT12**2-X1-X2)-
     &    QME*X1/XQ+0.5*QME*((X2**2-QME)*ST12**2-2.*X2)*XQ
          SIGL=(X2*ST12)**2-QME*(3.-CT12**2-2.5*(X1+X2)+X1*X2+QME)+
     &    0.5*QME*(X1**2-X1-QME)/XQ+0.5*QME*((X2**2-QME)*CT12**2-X2)*XQ
          SIGT=0.5*(X2**2-QME-0.5*QME*(X2**2-QME)/XQ)*ST12**2
          SIGI=((1.-0.5*QME*XQ)*(X2**2-QME)*ST12*CT12+QME*(1.-X1-X2+
     &    0.5*X1*X2+0.5*QME)*ST12/CT12)/SQ2
          SIGA=X2**2*ST12/SQ2
          SIGP=2.*(X1**2-X2**2*CT12)

C...Differential cross-sect for scalar gluons (no mass or QFD effects).
        ELSE
          SIGU=2.*(2.-X1-X2)**2-(X2*ST12)**2
          SIGL=(X2*ST12)**2
          SIGT=0.5*SIGL
          SIGI=-(2.-X1-X2)*X2*ST12/SQ2
          SIGA=0.
          SIGP=0.
        ENDIF
      ENDIF

C...Upper bounds for differential cross-section.
      HF1A=ABS(HF1)
      HF2A=ABS(HF2)
      HF3A=ABS(HF3)
      HF4A=ABS(HF4)
      SIGMAX=(2.*HF1A+HF3A+HF4A)*ABS(SIGU)+2.*(HF1A+HF3A+HF4A)*
     &ABS(SIGL)+2.*(HF1A+2.*HF3A+2.*HF4A)*ABS(SIGT)+2.*SQ2*
     &(HF1A+2.*HF3A+2.*HF4A)*ABS(SIGI)+4.*SQ2*HF2A*ABS(SIGA)+
     &2.*HF2A*ABS(SIGP)

C...Generate angular orientation according to differential cross-sect.
  100 CHI=PARU(2)*RLU(0)
      CTHE=2.*RLU(0)-1.
      PHI=PARU(2)*RLU(0)
      CCHI=COS(CHI)
      SCHI=SIN(CHI)
      C2CHI=COS(2.*CHI)
      S2CHI=SIN(2.*CHI)
      THE=ACOS(CTHE)
      STHE=SIN(THE)
      C2PHI=COS(2.*(PHI-PARJ(134)))
      S2PHI=SIN(2.*(PHI-PARJ(134)))
      SIG=((1.+CTHE**2)*HF1+STHE**2*(C2PHI*HF3-S2PHI*HF4))*SIGU+
     &2.*(STHE**2*HF1-STHE**2*(C2PHI*HF3-S2PHI*HF4))*SIGL+
     &2.*(STHE**2*C2CHI*HF1+((1.+CTHE**2)*C2CHI*C2PHI-2.*CTHE*S2CHI*
     &S2PHI)*HF3-((1.+CTHE**2)*C2CHI*S2PHI+2.*CTHE*S2CHI*C2PHI)*HF4)*
     &SIGT-2.*SQ2*(2.*STHE*CTHE*CCHI*HF1-2.*STHE*(CTHE*CCHI*C2PHI-
     &SCHI*S2PHI)*HF3+2.*STHE*(CTHE*CCHI*S2PHI+SCHI*C2PHI)*HF4)*SIGI+
     &4.*SQ2*STHE*CCHI*HF2*SIGA+2.*CTHE*HF2*SIGP
      IF(SIG.LT.SIGMAX*RLU(0)) GOTO 100

      RETURN
      END
