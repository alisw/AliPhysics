*CMZ :          17/07/98  15.44.33  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUSHOW(IP1,IP2,QMAX)

C...Purpose: to generate timelike parton showers from given partons.
      IMPLICIT DOUBLE PRECISION(D)
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
      DIMENSION PMTH(5,40),PS(5),PMA(4),PMSD(4),IEP(4),IPA(4),
     &KFLA(4),KFLD(4),KFL(4),ITRY(4),ISI(4),ISL(4),DP(4),DPT(5,4)

C...Initialization of cutoff masses etc.
      IF(MSTJ(41).LE.0.OR.(MSTJ(41).EQ.1.AND.QMAX.LE.PARJ(82)).OR.
     &QMAX.LE.MIN(PARJ(82),PARJ(83)).OR.MSTJ(41).GE.3) RETURN
      PMTH(1,21)=ULMASS(21)
      PMTH(2,21)=SQRT(PMTH(1,21)**2+0.25*PARJ(82)**2)
      PMTH(3,21)=2.*PMTH(2,21)
      PMTH(4,21)=PMTH(3,21)
      PMTH(5,21)=PMTH(3,21)
      PMTH(1,22)=ULMASS(22)
      PMTH(2,22)=SQRT(PMTH(1,22)**2+0.25*PARJ(83)**2)
      PMTH(3,22)=2.*PMTH(2,22)
      PMTH(4,22)=PMTH(3,22)
      PMTH(5,22)=PMTH(3,22)
      PMQTH1=PARJ(82)
      IF(MSTJ(41).EQ.2) PMQTH1=MIN(PARJ(82),PARJ(83))
      PMQTH2=PMTH(2,21)
      IF(MSTJ(41).EQ.2) PMQTH2=MIN(PMTH(2,21),PMTH(2,22))
      DO 100 IF=1,8
      PMTH(1,IF)=ULMASS(IF)
      PMTH(2,IF)=SQRT(PMTH(1,IF)**2+0.25*PMQTH1**2)
      PMTH(3,IF)=PMTH(2,IF)+PMQTH2
      PMTH(4,IF)=SQRT(PMTH(1,IF)**2+0.25*PARJ(82)**2)+PMTH(2,21)
  100 PMTH(5,IF)=SQRT(PMTH(1,IF)**2+0.25*PARJ(83)**2)+PMTH(2,22)
      PT2MIN=MAX(0.5*PARJ(82),1.1*PARJ(81))**2
      ALAMS=PARJ(81)**2
      ALFM=LOG(PT2MIN/ALAMS)

C...Store positions of shower initiating partons.
      M3JC=0
      IF(IP1.GT.0.AND.IP1.LE.MIN(N,MSTU(4)-MSTU(32)).AND.IP2.EQ.0) THEN
        NPA=1
        IPA(1)=IP1
      ELSEIF(MIN(IP1,IP2).GT.0.AND.MAX(IP1,IP2).LE.MIN(N,MSTU(4)-
     &MSTU(32))) THEN
        NPA=2
        IPA(1)=IP1
        IPA(2)=IP2
      ELSEIF(IP1.GT.0.AND.IP1.LE.MIN(N,MSTU(4)-MSTU(32)).AND.IP2.LT.0.
     &AND.IP2.GE.-3) THEN
        NPA=IABS(IP2)
        DO 110 I=1,NPA
  110   IPA(I)=IP1+I-1
      ELSE
        CALL LUERRM(12,
     &  '(LUSHOW:) failed to reconstruct showering system')
        IF(MSTU(21).GE.1) RETURN
      ENDIF

C...Check on phase space available for emission.
      IREJ=0
      DO 120 J=1,5
  120 PS(J)=0.
      PM=0.
      DO 130 I=1,NPA
      KFLA(I)=IABS(K(IPA(I),2))
      PMA(I)=P(IPA(I),5)
      IF(KFLA(I).NE.0.AND.(KFLA(I).LE.8.OR.KFLA(I).EQ.21))
     &PMA(I)=PMTH(3,KFLA(I))
      PM=PM+PMA(I)
      IF(KFLA(I).EQ.0.OR.(KFLA(I).GT.8.AND.KFLA(I).NE.21).OR.
     &PMA(I).GT.QMAX) IREJ=IREJ+1
      DO 130 J=1,4
  130 PS(J)=PS(J)+P(IPA(I),J)
      IF(IREJ.EQ.NPA) RETURN
      PS(5)=SQRT(MAX(0.,PS(4)**2-PS(1)**2-PS(2)**2-PS(3)**2))
      IF(NPA.EQ.1) PS(5)=PS(4)
      IF(PS(5).LE.PM+PMQTH1) RETURN
      IF(NPA.EQ.2.AND.MSTJ(47).GE.1) THEN
        IF(KFLA(1).GE.1.AND.KFLA(1).LE.8.AND.KFLA(2).GE.1.AND.
     &  KFLA(2).LE.8) M3JC=1
        IF(MSTJ(47).GE.2) M3JC=1
      ENDIF

C...Define imagined single initiator of shower for parton system.
      NS=N
      IF(N.GT.MSTU(4)-MSTU(32)-5) THEN
        CALL LUERRM(11,'(LUSHOW:) no more memory left in LUJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      IF(NPA.GE.2) THEN
        K(N+1,1)=11
        K(N+1,2)=21
        K(N+1,3)=0
        K(N+1,4)=0
        K(N+1,5)=0
        P(N+1,1)=0.
        P(N+1,2)=0.
        P(N+1,3)=0.
        P(N+1,4)=PS(5)
        P(N+1,5)=PS(5)
        V(N+1,5)=PS(5)**2
        N=N+1
      ENDIF

C...Loop over partons that may branch.
      NEP=NPA
      IM=NS
      IF(NPA.EQ.1) IM=NS-1
  140 IM=IM+1
      IF(N.GT.NS) THEN
        IF(IM.GT.N) GOTO 380
        KFLM=IABS(K(IM,2))
        IF(KFLM.EQ.0.OR.(KFLM.GT.8.AND.KFLM.NE.21)) GOTO 140
        IF(P(IM,5).LT.PMTH(2,KFLM)) GOTO 140
        IGM=K(IM,3)
      ELSE
        IGM=-1
      ENDIF
      IF(N+NEP.GT.MSTU(4)-MSTU(32)-5) THEN
        CALL LUERRM(11,'(LUSHOW:) no more memory left in LUJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF

C...Position of aunt (sister to branching parton).
C...Origin and flavour of daughters.
      IAU=0
      IF(IGM.GT.0) THEN
        IF(K(IM-1,3).EQ.IGM) IAU=IM-1
        IF(N.GE.IM+1.AND.K(IM+1,3).EQ.IGM) IAU=IM+1
      ENDIF
      IF(IGM.GE.0) THEN
        K(IM,4)=N+1
        DO 150 I=1,NEP
  150   K(N+I,3)=IM
      ELSE
        K(N+1,3)=IPA(1)
      ENDIF
      IF(IGM.LE.0) THEN
        DO 160 I=1,NEP
  160   K(N+I,2)=K(IPA(I),2)
      ELSEIF(KFLM.NE.21) THEN
        K(N+1,2)=K(IM,2)
        K(N+2,2)=K(IM,5)
      ELSEIF(K(IM,5).EQ.21) THEN
        K(N+1,2)=21
        K(N+2,2)=21
      ELSE
        K(N+1,2)=K(IM,5)
        K(N+2,2)=-K(IM,5)
      ENDIF

C...Reset flags on daughers and tries made.
      DO 170 IP=1,NEP
      K(N+IP,1)=3
      K(N+IP,4)=0
      K(N+IP,5)=0
      KFLD(IP)=IABS(K(N+IP,2))
      ITRY(IP)=0
      ISL(IP)=0
      ISI(IP)=0
  170 IF(KFLD(IP).GT.0.AND.(KFLD(IP).LE.8.OR.KFLD(IP).EQ.21)) ISI(IP)=1
      ISLM=0

C...Maximum virtuality of daughters.
      IF(IGM.LE.0) THEN
        DO 180 I=1,NPA
        IF(NPA.GE.3) P(N+I,4)=(PS(4)*P(IPA(I),4)-PS(1)*P(IPA(I),1)-
     &  PS(2)*P(IPA(I),2)-PS(3)*P(IPA(I),3))/PS(5)
        P(N+I,5)=MIN(QMAX,PS(5))
        IF(NPA.GE.3) P(N+I,5)=MIN(P(N+I,5),P(N+I,4))
  180   IF(ISI(I).EQ.0) P(N+I,5)=P(IPA(I),5)
      ELSE
        IF(MSTJ(43).LE.2) PEM=V(IM,2)
        IF(MSTJ(43).GE.3) PEM=P(IM,4)
        P(N+1,5)=MIN(P(IM,5),V(IM,1)*PEM)
        P(N+2,5)=MIN(P(IM,5),(1.-V(IM,1))*PEM)
        IF(K(N+2,2).EQ.22) P(N+2,5)=PMTH(1,22)
      ENDIF
      DO 190 I=1,NEP
      PMSD(I)=P(N+I,5)
      IF(ISI(I).EQ.1) THEN
        IF(P(N+I,5).LE.PMTH(3,KFLD(I))) P(N+I,5)=PMTH(1,KFLD(I))
      ENDIF
  190 V(N+I,5)=P(N+I,5)**2

C...Choose one of the daughters for evolution.
  200 INUM=0
      IF(NEP.EQ.1) INUM=1
      DO 210 I=1,NEP
  210 IF(INUM.EQ.0.AND.ISL(I).EQ.1) INUM=I
      DO 220 I=1,NEP
      IF(INUM.EQ.0.AND.ITRY(I).EQ.0.AND.ISI(I).EQ.1) THEN
        IF(P(N+I,5).GE.PMTH(2,KFLD(I))) INUM=I
      ENDIF
  220 CONTINUE
      IF(INUM.EQ.0) THEN
        RMAX=0.
        DO 230 I=1,NEP
        IF(ISI(I).EQ.1.AND.PMSD(I).GE.PMQTH2) THEN
          RPM=P(N+I,5)/PMSD(I)
          IF(RPM.GT.RMAX.AND.P(N+I,5).GE.PMTH(2,KFLD(I))) THEN
            RMAX=RPM
            INUM=I
          ENDIF
        ENDIF
  230   CONTINUE
      ENDIF

C...Store information on choice of evolving daughter.
      INUM=MAX(1,INUM)
      IEP(1)=N+INUM
      DO 240 I=2,NEP
      IEP(I)=IEP(I-1)+1
  240 IF(IEP(I).GT.N+NEP) IEP(I)=N+1
      DO 250 I=1,NEP
  250 KFL(I)=IABS(K(IEP(I),2))
      ITRY(INUM)=ITRY(INUM)+1
      IF(ITRY(INUM).GT.200) THEN
        CALL LUERRM(14,'(LUSHOW:) caught in infinite loop')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      Z=0.5
      IF(KFL(1).EQ.0.OR.(KFL(1).GT.8.AND.KFL(1).NE.21)) GOTO 300
      IF(P(IEP(1),5).LT.PMTH(2,KFL(1))) GOTO 300

C...Calculate allowed z range.
      IF(NEP.EQ.1) THEN
        PMED=PS(4)
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN
        PMED=P(IM,5)
      ELSE
        IF(INUM.EQ.1) PMED=V(IM,1)*PEM
        IF(INUM.EQ.2) PMED=(1.-V(IM,1))*PEM
      ENDIF
      IF(MOD(MSTJ(43),2).EQ.1) THEN
        ZC=PMTH(2,21)/PMED
        ZCE=PMTH(2,22)/PMED
      ELSE
        ZC=0.5*(1.-SQRT(MAX(0.,1.-(2.*PMTH(2,21)/PMED)**2)))
        IF(ZC.LT.1E-4) ZC=(PMTH(2,21)/PMED)**2
        ZCE=0.5*(1.-SQRT(MAX(0.,1.-(2.*PMTH(2,22)/PMED)**2)))
        IF(ZCE.LT.1E-4) ZCE=(PMTH(2,22)/PMED)**2
      ENDIF
      ZC=MIN(ZC,0.491)
      ZCE=MIN(ZCE,0.491)
      IF((MSTJ(41).EQ.1.AND.ZC.GT.0.49).OR.(MSTJ(41).EQ.2.AND.
     &MIN(ZC,ZCE).GT.0.49)) THEN
        P(IEP(1),5)=PMTH(1,KFL(1))
        V(IEP(1),5)=P(IEP(1),5)**2
        GOTO 300
      ENDIF

C...Integral of Altarelli-Parisi z kernel for QCD.
      IF(MSTJ(49).EQ.0.AND.KFL(1).EQ.21) THEN
        FBR=6.*LOG((1.-ZC)/ZC)+MSTJ(45)*(0.5-ZC)
      ELSEIF(MSTJ(49).EQ.0) THEN
        FBR=(8./3.)*LOG((1.-ZC)/ZC)

C...Integral of Altarelli-Parisi z kernel for scalar gluon.
      ELSEIF(MSTJ(49).EQ.1.AND.KFL(1).EQ.21) THEN
        FBR=(PARJ(87)+MSTJ(45)*PARJ(88))*(1.-2.*ZC)
      ELSEIF(MSTJ(49).EQ.1) THEN
        FBR=(1.-2.*ZC)/3.
        IF(IGM.EQ.0.AND.M3JC.EQ.1) FBR=4.*FBR

C...Integral of Altarelli-Parisi z kernel for Abelian vector gluon.
      ELSEIF(KFL(1).EQ.21) THEN
        FBR=6.*MSTJ(45)*(0.5-ZC)
      ELSE
        FBR=2.*LOG((1.-ZC)/ZC)
      ENDIF

C...Integral of Altarelli-Parisi kernel for photon emission.
      IF(MSTJ(41).EQ.2.AND.KFL(1).GE.1.AND.KFL(1).LE.8)
     &FBRE=(KCHG(KFL(1),1)/3.)**2*2.*LOG((1.-ZCE)/ZCE)

C...Inner veto algorithm starts. Find maximum mass for evolution.
  260 PMS=V(IEP(1),5)
      IF(IGM.GE.0) THEN
        PM2=0.
        DO 270 I=2,NEP
        PM=P(IEP(I),5)
        IF(KFL(I).GT.0.AND.(KFL(I).LE.8.OR.KFL(I).EQ.21)) PM=
     &  PMTH(2,KFL(I))
  270   PM2=PM2+PM
        PMS=MIN(PMS,(P(IM,5)-PM2)**2)
      ENDIF

C...Select mass for daughter in QCD evolution.
      B0=27./6.
      DO 280 IF=4,MSTJ(45)
  280 IF(PMS.GT.4.*PMTH(2,IF)**2) B0=(33.-2.*IF)/6.
      IF(MSTJ(44).LE.0) THEN
        PMSQCD=PMS*EXP(MAX(-80.,LOG(RLU(0))*PARU(2)/(PARU(111)*FBR)))
      ELSEIF(MSTJ(44).EQ.1) THEN
        PMSQCD=4.*ALAMS*(0.25*PMS/ALAMS)**(RLU(0)**(B0/FBR))
      ELSE
        PMSQCD=PMS*RLU(0)**(ALFM*B0/FBR)
      ENDIF
      IF(ZC.GT.0.49.OR.PMSQCD.LE.PMTH(4,KFL(1))**2) PMSQCD=
     &PMTH(2,KFL(1))**2
      V(IEP(1),5)=PMSQCD
      MCE=1

C...Select mass for daughter in QED evolution.
      IF(MSTJ(41).EQ.2.AND.KFL(1).GE.1.AND.KFL(1).LE.8) THEN
        PMSQED=PMS*EXP(MAX(-80.,LOG(RLU(0))*PARU(2)/(PARU(101)*FBRE)))
        IF(ZCE.GT.0.49.OR.PMSQED.LE.PMTH(5,KFL(1))**2) PMSQED=
     &  PMTH(2,KFL(1))**2
        IF(PMSQED.GT.PMSQCD) THEN
          V(IEP(1),5)=PMSQED
          MCE=2
        ENDIF
      ENDIF

C...Check whether daughter mass below cutoff.
      P(IEP(1),5)=SQRT(V(IEP(1),5))
      IF(P(IEP(1),5).LE.PMTH(3,KFL(1))) THEN
        P(IEP(1),5)=PMTH(1,KFL(1))
        V(IEP(1),5)=P(IEP(1),5)**2
        GOTO 300
      ENDIF

C...Select z value of branching: q -> qgamma.
      IF(MCE.EQ.2) THEN
        Z=1.-(1.-ZCE)*(ZCE/(1.-ZCE))**RLU(0)
        IF(1.+Z**2.LT.2.*RLU(0)) GOTO 260
        K(IEP(1),5)=22

C...Select z value of branching: q -> qg, g -> gg, g -> qqbar.
      ELSEIF(MSTJ(49).NE.1.AND.KFL(1).NE.21) THEN
        Z=1.-(1.-ZC)*(ZC/(1.-ZC))**RLU(0)
        IF(1.+Z**2.LT.2.*RLU(0)) GOTO 260
        K(IEP(1),5)=21
      ELSEIF(MSTJ(49).EQ.0.AND.MSTJ(45)*(0.5-ZC).LT.RLU(0)*FBR) THEN
        Z=(1.-ZC)*(ZC/(1.-ZC))**RLU(0)
        IF(RLU(0).GT.0.5) Z=1.-Z
        IF((1.-Z*(1.-Z))**2.LT.RLU(0)) GOTO 260
        K(IEP(1),5)=21
      ELSEIF(MSTJ(49).NE.1) THEN
        Z=ZC+(1.-2.*ZC)*RLU(0)
        IF(Z**2+(1.-Z)**2.LT.RLU(0)) GOTO 260
        KFLB=1+INT(MSTJ(45)*RLU(0))
        PMQ=4.*PMTH(2,KFLB)**2/V(IEP(1),5)
        IF(PMQ.GE.1.) GOTO 260
        PMQ0=4.*PMTH(2,21)**2/V(IEP(1),5)
        IF(MOD(MSTJ(43),2).EQ.0.AND.(1.+0.5*PMQ)*SQRT(1.-PMQ).LT.
     &  RLU(0)*(1.+0.5*PMQ0)*SQRT(1.-PMQ0)) GOTO 260
        K(IEP(1),5)=KFLB

C...Ditto for scalar gluon model.
      ELSEIF(KFL(1).NE.21) THEN
        Z=1.-SQRT(ZC**2+RLU(0)*(1.-2.*ZC))
        K(IEP(1),5)=21
      ELSEIF(RLU(0)*(PARJ(87)+MSTJ(45)*PARJ(88)).LE.PARJ(87)) THEN
        Z=ZC+(1.-2.*ZC)*RLU(0)
        K(IEP(1),5)=21
      ELSE
        Z=ZC+(1.-2.*ZC)*RLU(0)
        KFLB=1+INT(MSTJ(45)*RLU(0))
        PMQ=4.*PMTH(2,KFLB)**2/V(IEP(1),5)
        IF(PMQ.GE.1.) GOTO 260
        K(IEP(1),5)=KFLB
      ENDIF
      IF(MCE.EQ.1.AND.MSTJ(44).GE.2) THEN
        IF(Z*(1.-Z)*V(IEP(1),5).LT.PT2MIN) GOTO 260
        IF(ALFM/LOG(V(IEP(1),5)*Z*(1.-Z)/ALAMS).LT.RLU(0)) GOTO 260
      ENDIF

C...Check if z consistent with chosen m.
      IF(KFL(1).EQ.21) THEN
        KFLGD1=IABS(K(IEP(1),5))
        KFLGD2=KFLGD1
      ELSE
        KFLGD1=KFL(1)
        KFLGD2=IABS(K(IEP(1),5))
      ENDIF
      IF(NEP.EQ.1) THEN
        PED=PS(4)
      ELSEIF(NEP.GE.3) THEN
        PED=P(IEP(1),4)
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN
        PED=0.5*(V(IM,5)+V(IEP(1),5)-PM2**2)/P(IM,5)
      ELSE
        IF(IEP(1).EQ.N+1) PED=V(IM,1)*PEM
        IF(IEP(1).EQ.N+2) PED=(1.-V(IM,1))*PEM
      ENDIF
      IF(MOD(MSTJ(43),2).EQ.1) THEN
        PMQTH3=0.5*PARJ(82)
        IF(KFLGD2.EQ.22) PMQTH3=0.5*PARJ(83)
        PMQ1=(PMTH(1,KFLGD1)**2+PMQTH3**2)/V(IEP(1),5)
        PMQ2=(PMTH(1,KFLGD2)**2+PMQTH3**2)/V(IEP(1),5)
        ZD=SQRT(MAX(0.,(1.-V(IEP(1),5)/PED**2)*((1.-PMQ1-PMQ2)**2-
     &  4.*PMQ1*PMQ2)))
        ZH=1.+PMQ1-PMQ2
      ELSE
        ZD=SQRT(MAX(0.,1.-V(IEP(1),5)/PED**2))
        ZH=1.
      ENDIF
      ZL=0.5*(ZH-ZD)
      ZU=0.5*(ZH+ZD)
      IF(Z.LT.ZL.OR.Z.GT.ZU) GOTO 260
      IF(KFL(1).EQ.21) V(IEP(1),3)=LOG(ZU*(1.-ZL)/MAX(1E-20,ZL*
     &(1.-ZU)))
      IF(KFL(1).NE.21) V(IEP(1),3)=LOG((1.-ZL)/MAX(1E-10,1.-ZU))

C...Three-jet matrix element correction.
      IF(IGM.EQ.0.AND.M3JC.EQ.1) THEN
        X1=Z*(1.+V(IEP(1),5)/V(NS+1,5))
        X2=1.-V(IEP(1),5)/V(NS+1,5)
        X3=(1.-X1)+(1.-X2)
        IF(MCE.EQ.2) THEN
          KI1=K(IPA(INUM),2)
          KI2=K(IPA(3-INUM),2)
          QF1=KCHG(IABS(KI1),1)*ISIGN(1,KI1)/3.
          QF2=KCHG(IABS(KI2),1)*ISIGN(1,KI2)/3.
          WSHOW=QF1**2*(1.-X1)/X3*(1.+(X1/(2.-X2))**2)+
     &    QF2**2*(1.-X2)/X3*(1.+(X2/(2.-X1))**2)
          WME=(QF1*(1.-X1)/X3-QF2*(1.-X2)/X3)**2*(X1**2+X2**2)
        ELSEIF(MSTJ(49).NE.1) THEN
          WSHOW=1.+(1.-X1)/X3*(X1/(2.-X2))**2+
     &    (1.-X2)/X3*(X2/(2.-X1))**2
          WME=X1**2+X2**2
        ELSE
          WSHOW=4.*X3*((1.-X1)/(2.-X2)**2+(1.-X2)/(2.-X1)**2)
          WME=X3**2
        ENDIF
        IF(WME.LT.RLU(0)*WSHOW) GOTO 260

C...Impose angular ordering by rejection of nonordered emission.
      ELSEIF(MCE.EQ.1.AND.IGM.GT.0.AND.MSTJ(42).GE.2) THEN
        MAOM=1
        ZM=V(IM,1)
        IF(IEP(1).EQ.N+2) ZM=1.-V(IM,1)
        THE2ID=Z*(1.-Z)*(ZM*P(IM,4))**2/V(IEP(1),5)
        IAOM=IM
  290   IF(K(IAOM,5).EQ.22) THEN
          IAOM=K(IAOM,3)
          IF(K(IAOM,3).LE.NS) MAOM=0
          IF(MAOM.EQ.1) GOTO 290
        ENDIF
        IF(MAOM.EQ.1) THEN
          THE2IM=V(IAOM,1)*(1.-V(IAOM,1))*P(IAOM,4)**2/V(IAOM,5)
          IF(THE2ID.LT.THE2IM) GOTO 260
        ENDIF
      ENDIF

C...Impose user-defined maximum angle at first branching.
      IF(MSTJ(48).EQ.1) THEN
        IF(NEP.EQ.1.AND.IM.EQ.NS) THEN
          THE2ID=Z*(1.-Z)*PS(4)**2/V(IEP(1),5)
          IF(THE2ID.LT.1./PARJ(85)**2) GOTO 260
        ELSEIF(NEP.EQ.2.AND.IEP(1).EQ.NS+2) THEN
          THE2ID=Z*(1.-Z)*(0.5*P(IM,4))**2/V(IEP(1),5)
          IF(THE2ID.LT.1./PARJ(85)**2) GOTO 260
        ELSEIF(NEP.EQ.2.AND.IEP(1).EQ.NS+3) THEN
          THE2ID=Z*(1.-Z)*(0.5*P(IM,4))**2/V(IEP(1),5)
          IF(THE2ID.LT.1./PARJ(86)**2) GOTO 260
        ENDIF
      ENDIF

C...End of inner veto algorithm. Check if only one leg evolved so far.
  300 V(IEP(1),1)=Z
      ISL(1)=0
      ISL(2)=0
      IF(NEP.EQ.1) GOTO 330
      IF(NEP.EQ.2.AND.P(IEP(1),5)+P(IEP(2),5).GE.P(IM,5)) GOTO 200
      DO 310 I=1,NEP
      IF(ITRY(I).EQ.0.AND.KFLD(I).GT.0.AND.(KFLD(I).LE.8.OR.KFLD(I).EQ.
     &21)) THEN
        IF(P(N+I,5).GE.PMTH(2,KFLD(I))) GOTO 200
      ENDIF
  310 CONTINUE

C...Check if chosen multiplet m1,m2,z1,z2 is physical.
      IF(NEP.EQ.3) THEN
        PA1S=(P(N+1,4)+P(N+1,5))*(P(N+1,4)-P(N+1,5))
        PA2S=(P(N+2,4)+P(N+2,5))*(P(N+2,4)-P(N+2,5))
        PA3S=(P(N+3,4)+P(N+3,5))*(P(N+3,4)-P(N+3,5))
        PTS=0.25*(2.*PA1S*PA2S+2.*PA1S*PA3S+2.*PA2S*PA3S-
     &  PA1S**2-PA2S**2-PA3S**2)/PA1S
        IF(PTS.LE.0.) GOTO 200
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2.OR.MOD(MSTJ(43),2).EQ.0) THEN
        DO 320 I1=N+1,N+2
        KFLDA=IABS(K(I1,2))
        IF(KFLDA.EQ.0.OR.(KFLDA.GT.8.AND.KFLDA.NE.21)) GOTO 320
        IF(P(I1,5).LT.PMTH(2,KFLDA)) GOTO 320
        IF(KFLDA.EQ.21) THEN
          KFLGD1=IABS(K(I1,5))
          KFLGD2=KFLGD1
        ELSE
          KFLGD1=KFLDA
          KFLGD2=IABS(K(I1,5))
        ENDIF
        I2=2*N+3-I1
        IF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN
          PED=0.5*(V(IM,5)+V(I1,5)-V(I2,5))/P(IM,5)
        ELSE
          IF(I1.EQ.N+1) ZM=V(IM,1)
          IF(I1.EQ.N+2) ZM=1.-V(IM,1)
          PML=SQRT((V(IM,5)-V(N+1,5)-V(N+2,5))**2-
     &    4.*V(N+1,5)*V(N+2,5))
          PED=PEM*(0.5*(V(IM,5)-PML+V(I1,5)-V(I2,5))+PML*ZM)/V(IM,5)
        ENDIF
        IF(MOD(MSTJ(43),2).EQ.1) THEN
          PMQTH3=0.5*PARJ(82)
          IF(KFLGD2.EQ.22) PMQTH3=0.5*PARJ(83)
          PMQ1=(PMTH(1,KFLGD1)**2+PMQTH3**2)/V(I1,5)
          PMQ2=(PMTH(1,KFLGD2)**2+PMQTH3**2)/V(I1,5)
          ZD=SQRT(MAX(0.,(1.-V(I1,5)/PED**2)*((1.-PMQ1-PMQ2)**2-
     &    4.*PMQ1*PMQ2)))
          ZH=1.+PMQ1-PMQ2
        ELSE
          ZD=SQRT(MAX(0.,1.-V(I1,5)/PED**2))
          ZH=1.
        ENDIF
        ZL=0.5*(ZH-ZD)
        ZU=0.5*(ZH+ZD)
        IF(I1.EQ.N+1.AND.(V(I1,1).LT.ZL.OR.V(I1,1).GT.ZU)) ISL(1)=1
        IF(I1.EQ.N+2.AND.(V(I1,1).LT.ZL.OR.V(I1,1).GT.ZU)) ISL(2)=1
        IF(KFLDA.EQ.21) V(I1,4)=LOG(ZU*(1.-ZL)/MAX(1E-20,ZL*(1.-ZU)))
        IF(KFLDA.NE.21) V(I1,4)=LOG((1.-ZL)/MAX(1E-10,1.-ZU))
  320   CONTINUE
        IF(ISL(1).EQ.1.AND.ISL(2).EQ.1.AND.ISLM.NE.0) THEN
          ISL(3-ISLM)=0
          ISLM=3-ISLM
        ELSEIF(ISL(1).EQ.1.AND.ISL(2).EQ.1) THEN
          ZDR1=MAX(0.,V(N+1,3)/V(N+1,4)-1.)
          ZDR2=MAX(0.,V(N+2,3)/V(N+2,4)-1.)
          IF(ZDR2.GT.RLU(0)*(ZDR1+ZDR2)) ISL(1)=0
          IF(ISL(1).EQ.1) ISL(2)=0
          IF(ISL(1).EQ.0) ISLM=1
          IF(ISL(2).EQ.0) ISLM=2
        ENDIF
        IF(ISL(1).EQ.1.OR.ISL(2).EQ.1) GOTO 200
      ENDIF
      IF(IGM.GT.0.AND.MOD(MSTJ(43),2).EQ.1.AND.(P(N+1,5).GE.
     &PMTH(2,KFLD(1)).OR.P(N+2,5).GE.PMTH(2,KFLD(2)))) THEN
        PMQ1=V(N+1,5)/V(IM,5)
        PMQ2=V(N+2,5)/V(IM,5)
        ZD=SQRT(MAX(0.,(1.-V(IM,5)/PEM**2)*((1.-PMQ1-PMQ2)**2-
     &  4.*PMQ1*PMQ2)))
        ZH=1.+PMQ1-PMQ2
        ZL=0.5*(ZH-ZD)
        ZU=0.5*(ZH+ZD)
        IF(V(IM,1).LT.ZL.OR.V(IM,1).GT.ZU) GOTO 200
      ENDIF

C...Accepted branch. Construct four-momentum for initial partons.
  330 MAZIP=0
      MAZIC=0
      IF(NEP.EQ.1) THEN
        P(N+1,1)=0.
        P(N+1,2)=0.
        P(N+1,3)=SQRT(MAX(0.,(P(IPA(1),4)+P(N+1,5))*(P(IPA(1),4)-
     &  P(N+1,5))))
        P(N+1,4)=P(IPA(1),4)
        V(N+1,2)=P(N+1,4)
      ELSEIF(IGM.EQ.0.AND.NEP.EQ.2) THEN
        PED1=0.5*(V(IM,5)+V(N+1,5)-V(N+2,5))/P(IM,5)
        P(N+1,1)=0.
        P(N+1,2)=0.
        P(N+1,3)=SQRT(MAX(0.,(PED1+P(N+1,5))*(PED1-P(N+1,5))))
        P(N+1,4)=PED1
        P(N+2,1)=0.
        P(N+2,2)=0.
        P(N+2,3)=-P(N+1,3)
        P(N+2,4)=P(IM,5)-PED1
        V(N+1,2)=P(N+1,4)
        V(N+2,2)=P(N+2,4)
      ELSEIF(NEP.EQ.3) THEN
        P(N+1,1)=0.
        P(N+1,2)=0.
        P(N+1,3)=SQRT(MAX(0.,PA1S))
        P(N+2,1)=SQRT(PTS)
        P(N+2,2)=0.
        P(N+2,3)=0.5*(PA3S-PA2S-PA1S)/P(N+1,3)
        P(N+3,1)=-P(N+2,1)
        P(N+3,2)=0.
        P(N+3,3)=-(P(N+1,3)+P(N+2,3))
        V(N+1,2)=P(N+1,4)
        V(N+2,2)=P(N+2,4)
        V(N+3,2)=P(N+3,4)

C...Construct transverse momentum for ordinary branching in shower.
      ELSE
        ZM=V(IM,1)
        PZM=SQRT(MAX(0.,(PEM+P(IM,5))*(PEM-P(IM,5))))
        PMLS=(V(IM,5)-V(N+1,5)-V(N+2,5))**2-4.*V(N+1,5)*V(N+2,5)
        IF(PZM.LE.0.) THEN
          PTS=0.
        ELSEIF(MOD(MSTJ(43),2).EQ.1) THEN
          PTS=(PEM**2*(ZM*(1.-ZM)*V(IM,5)-(1.-ZM)*V(N+1,5)-
     &    ZM*V(N+2,5))-0.25*PMLS)/PZM**2
        ELSE
          PTS=PMLS*(ZM*(1.-ZM)*PEM**2/V(IM,5)-0.25)/PZM**2
        ENDIF
        PT=SQRT(MAX(0.,PTS))

C...Find coefficient of azimuthal asymmetry due to gluon polarization.
        HAZIP=0.
        IF(MSTJ(49).NE.1.AND.MOD(MSTJ(46),2).EQ.1.AND.K(IM,2).EQ.21.
     &  AND.IAU.NE.0) THEN
          IF(K(IGM,3).NE.0) MAZIP=1
          ZAU=V(IGM,1)
          IF(IAU.EQ.IM+1) ZAU=1.-V(IGM,1)
          IF(MAZIP.EQ.0) ZAU=0.
          IF(K(IGM,2).NE.21) THEN
            HAZIP=2.*ZAU/(1.+ZAU**2)
          ELSE
            HAZIP=(ZAU/(1.-ZAU*(1.-ZAU)))**2
          ENDIF
          IF(K(N+1,2).NE.21) THEN
            HAZIP=HAZIP*(-2.*ZM*(1.-ZM))/(1.-2.*ZM*(1.-ZM))
          ELSE
            HAZIP=HAZIP*(ZM*(1.-ZM)/(1.-ZM*(1.-ZM)))**2
          ENDIF
        ENDIF

C...Find coefficient of azimuthal asymmetry due to soft gluon
C...interference.
        HAZIC=0.
        IF(MSTJ(49).NE.2.AND.MSTJ(46).GE.2.AND.(K(N+1,2).EQ.21.OR.
     &  K(N+2,2).EQ.21).AND.IAU.NE.0) THEN
          IF(K(IGM,3).NE.0) MAZIC=N+1
          IF(K(IGM,3).NE.0.AND.K(N+1,2).NE.21) MAZIC=N+2
          IF(K(IGM,3).NE.0.AND.K(N+1,2).EQ.21.AND.K(N+2,2).EQ.21.AND.
     &    ZM.GT.0.5) MAZIC=N+2
          IF(K(IAU,2).EQ.22) MAZIC=0
          ZS=ZM
          IF(MAZIC.EQ.N+2) ZS=1.-ZM
          ZGM=V(IGM,1)
          IF(IAU.EQ.IM-1) ZGM=1.-V(IGM,1)
          IF(MAZIC.EQ.0) ZGM=1.
          HAZIC=(P(IM,5)/P(IGM,5))*SQRT((1.-ZS)*(1.-ZGM)/(ZS*ZGM))
          HAZIC=MIN(0.95,HAZIC)
        ENDIF
      ENDIF

C...Construct kinematics for ordinary branching in shower.
  340 IF(NEP.EQ.2.AND.IGM.GT.0) THEN
        IF(MOD(MSTJ(43),2).EQ.1) THEN
          P(N+1,4)=PEM*V(IM,1)
        ELSE
          P(N+1,4)=PEM*(0.5*(V(IM,5)-SQRT(PMLS)+V(N+1,5)-V(N+2,5))+
     &    SQRT(PMLS)*ZM)/V(IM,5)
        ENDIF
        PHI=PARU(2)*RLU(0)
        P(N+1,1)=PT*COS(PHI)
        P(N+1,2)=PT*SIN(PHI)
        IF(PZM.GT.0.) THEN
          P(N+1,3)=0.5*(V(N+2,5)-V(N+1,5)-V(IM,5)+2.*PEM*P(N+1,4))/PZM
        ELSE
          P(N+1,3)=0.
        ENDIF
        P(N+2,1)=-P(N+1,1)
        P(N+2,2)=-P(N+1,2)
        P(N+2,3)=PZM-P(N+1,3)
        P(N+2,4)=PEM-P(N+1,4)
        IF(MSTJ(43).LE.2) THEN
          V(N+1,2)=(PEM*P(N+1,4)-PZM*P(N+1,3))/P(IM,5)
          V(N+2,2)=(PEM*P(N+2,4)-PZM*P(N+2,3))/P(IM,5)
        ENDIF
      ENDIF

C...Rotate and boost daughters.
      IF(IGM.GT.0) THEN
        IF(MSTJ(43).LE.2) THEN
          BEX=P(IGM,1)/P(IGM,4)
          BEY=P(IGM,2)/P(IGM,4)
          BEZ=P(IGM,3)/P(IGM,4)
          GA=P(IGM,4)/P(IGM,5)
          GABEP=GA*(GA*(BEX*P(IM,1)+BEY*P(IM,2)+BEZ*P(IM,3))/(1.+GA)-
     &    P(IM,4))
        ELSE
          BEX=0.
          BEY=0.
          BEZ=0.
          GA=1.
          GABEP=0.
        ENDIF
        THE=ULANGL(P(IM,3)+GABEP*BEZ,SQRT((P(IM,1)+GABEP*BEX)**2+
     &  (P(IM,2)+GABEP*BEY)**2))
        PHI=ULANGL(P(IM,1)+GABEP*BEX,P(IM,2)+GABEP*BEY)
        DO 350 I=N+1,N+2
        DP(1)=COS(THE)*COS(PHI)*P(I,1)-SIN(PHI)*P(I,2)+
     &  SIN(THE)*COS(PHI)*P(I,3)
        DP(2)=COS(THE)*SIN(PHI)*P(I,1)+COS(PHI)*P(I,2)+
     &  SIN(THE)*SIN(PHI)*P(I,3)
        DP(3)=-SIN(THE)*P(I,1)+COS(THE)*P(I,3)
        DP(4)=P(I,4)
        DBP=BEX*DP(1)+BEY*DP(2)+BEZ*DP(3)
        DGABP=GA*(GA*DBP/(1D0+GA)+DP(4))
        P(I,1)=DP(1)+DGABP*BEX
        P(I,2)=DP(2)+DGABP*BEY
        P(I,3)=DP(3)+DGABP*BEZ
  350   P(I,4)=GA*(DP(4)+DBP)
      ENDIF

C...Weight with azimuthal distribution, if required.
      IF(MAZIP.NE.0.OR.MAZIC.NE.0) THEN
        DO 360 J=1,3
        DPT(1,J)=P(IM,J)
        DPT(2,J)=P(IAU,J)
  360   DPT(3,J)=P(N+1,J)
        DPMA=DPT(1,1)*DPT(2,1)+DPT(1,2)*DPT(2,2)+DPT(1,3)*DPT(2,3)
        DPMD=DPT(1,1)*DPT(3,1)+DPT(1,2)*DPT(3,2)+DPT(1,3)*DPT(3,3)
        DPMM=DPT(1,1)**2+DPT(1,2)**2+DPT(1,3)**2
        DO 370 J=1,3
        DPT(4,J)=DPT(2,J)-DPMA*DPT(1,J)/DPMM
  370   DPT(5,J)=DPT(3,J)-DPMD*DPT(1,J)/DPMM
        DPT(4,4)=SQRT(DPT(4,1)**2+DPT(4,2)**2+DPT(4,3)**2)
        DPT(5,4)=SQRT(DPT(5,1)**2+DPT(5,2)**2+DPT(5,3)**2)
        IF(MIN(DPT(4,4),DPT(5,4)).GT.0.1*PARJ(82)) THEN
          CAD=(DPT(4,1)*DPT(5,1)+DPT(4,2)*DPT(5,2)+
     &    DPT(4,3)*DPT(5,3))/(DPT(4,4)*DPT(5,4))
          IF(MAZIP.NE.0) THEN
            IF(1.+HAZIP*(2.*CAD**2-1.).LT.RLU(0)*(1.+ABS(HAZIP)))
     &      GOTO 340
          ENDIF
          IF(MAZIC.NE.0) THEN
            IF(MAZIC.EQ.N+2) CAD=-CAD
            IF((1.-HAZIC)*(1.-HAZIC*CAD)/(1.+HAZIC**2-2.*HAZIC*CAD).
     &      LT.RLU(0)) GOTO 340
          ENDIF
        ENDIF
      ENDIF

C...Continue loop over partons that may branch, until none left.
      IF(IGM.GE.0) K(IM,1)=14
      N=N+NEP
      NEP=2
      IF(N.GT.MSTU(4)-MSTU(32)-5) THEN
        CALL LUERRM(11,'(LUSHOW:) no more memory left in LUJETS')
        IF(MSTU(21).GE.1) N=NS
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      GOTO 140

C...Set information on imagined shower initiator.
  380 IF(NPA.GE.2) THEN
        K(NS+1,1)=11
        K(NS+1,2)=94
        K(NS+1,3)=IP1
        IF(IP2.GT.0.AND.IP2.LT.IP1) K(NS+1,3)=IP2
        K(NS+1,4)=NS+2
        K(NS+1,5)=NS+1+NPA
        IIM=1
      ELSE
        IIM=0
      ENDIF

C...Reconstruct string drawing information.
      DO 390 I=NS+1+IIM,N
      IF(K(I,1).LE.10.AND.K(I,2).EQ.22) THEN
        K(I,1)=1
      ELSEIF(K(I,1).LE.10) THEN
        K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))
        K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))
      ELSEIF(K(MOD(K(I,4),MSTU(5))+1,2).NE.22) THEN
        ID1=MOD(K(I,4),MSTU(5))
        IF(K(I,2).GE.1.AND.K(I,2).LE.8) ID1=MOD(K(I,4),MSTU(5))+1
        ID2=2*MOD(K(I,4),MSTU(5))+1-ID1
        K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))+ID1
        K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))+ID2
        K(ID1,4)=K(ID1,4)+MSTU(5)*I
        K(ID1,5)=K(ID1,5)+MSTU(5)*ID2
        K(ID2,4)=K(ID2,4)+MSTU(5)*ID1
        K(ID2,5)=K(ID2,5)+MSTU(5)*I
      ELSE
        ID1=MOD(K(I,4),MSTU(5))
        ID2=ID1+1
        K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))+ID1
        K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))+ID1
        K(ID1,4)=K(ID1,4)+MSTU(5)*I
        K(ID1,5)=K(ID1,5)+MSTU(5)*I
        K(ID2,4)=0
        K(ID2,5)=0
      ENDIF
  390 CONTINUE

C...Transformation from CM frame.
      IF(NPA.GE.2) THEN
        BEX=PS(1)/PS(4)
        BEY=PS(2)/PS(4)
        BEZ=PS(3)/PS(4)
        GA=PS(4)/PS(5)
        GABEP=GA*(GA*(BEX*P(IPA(1),1)+BEY*P(IPA(1),2)+BEZ*P(IPA(1),3))
     &  /(1.+GA)-P(IPA(1),4))
      ELSE
        BEX=0.
        BEY=0.
        BEZ=0.
        GABEP=0.
      ENDIF
      THE=ULANGL(P(IPA(1),3)+GABEP*BEZ,SQRT((P(IPA(1),1)
     &+GABEP*BEX)**2+(P(IPA(1),2)+GABEP*BEY)**2))
      PHI=ULANGL(P(IPA(1),1)+GABEP*BEX,P(IPA(1),2)+GABEP*BEY)
      IF(NPA.EQ.3) THEN
        CHI=ULANGL(COS(THE)*COS(PHI)*(P(IPA(2),1)+GABEP*BEX)+COS(THE)*
     &  SIN(PHI)*(P(IPA(2),2)+GABEP*BEY)-SIN(THE)*(P(IPA(2),3)+GABEP*
     &  BEZ),-SIN(PHI)*(P(IPA(2),1)+GABEP*BEX)+COS(PHI)*(P(IPA(2),2)+
     &  GABEP*BEY))
        MSTU(33)=1
        CALL LUDBRB(NS+1,N,0.,CHI,0D0,0D0,0D0)
      ENDIF
      DBEX=DBLE(BEX)
      DBEY=DBLE(BEY)
      DBEZ=DBLE(BEZ)
      MSTU(33)=1
      CALL LUDBRB(NS+1,N,THE,PHI,DBEX,DBEY,DBEZ)

C...Decay vertex of shower.
      DO 400 I=NS+1,N
      DO 400 J=1,5
  400 V(I,J)=V(IP1,J)

C...Delete trivial shower, else connect initiators.
      IF(N.EQ.NS+NPA+IIM) THEN
        N=NS
      ELSE
        DO 410 IP=1,NPA
        K(IPA(IP),1)=14
        K(IPA(IP),4)=K(IPA(IP),4)+NS+IIM+IP
        K(IPA(IP),5)=K(IPA(IP),5)+NS+IIM+IP
        K(NS+IIM+IP,3)=IPA(IP)
        IF(IIM.EQ.1.AND.MSTU(16).NE.2) K(NS+IIM+IP,3)=NS+1
        K(NS+IIM+IP,4)=MSTU(5)*IPA(IP)+K(NS+IIM+IP,4)
  410   K(NS+IIM+IP,5)=MSTU(5)*IPA(IP)+K(NS+IIM+IP,5)
      ENDIF

      RETURN
      END
