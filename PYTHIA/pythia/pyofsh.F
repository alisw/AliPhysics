 
C***********************************************************************
 
      SUBROUTINE PYOFSH(MOFSH,KFMO,KFD1,KFD2,PMMO,RET1,RET2)
 
C...Calculates partial width and differential cross-section maxima
C...of channels/processes not allowed on mass-shell, and selects
C...masses in such channels/processes.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3)
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT5/
      DIMENSION KFD(2),MBW(2),PMD(2),PGD(2),PMG(2),PML(2),PMU(2),
     &PMH(2),ATL(2),ATU(2),ATH(2),RMG(2),INX1(100),XPT1(100),
     &FPT1(100),INX2(100),XPT2(100),FPT2(100),WDTP(0:40),
     &WDTE(0:40,0:5)
 
C...Find if particles equal, maximum mass, matrix elements, etc.
      MINT(51)=0
      ISUB=MINT(1)
      KFD(1)=IABS(KFD1)
      KFD(2)=IABS(KFD2)
      MEQL=0
      IF(KFD(1).EQ.KFD(2)) MEQL=1
      MLM=0
      IF(MOFSH.GE.2.AND.MEQL.EQ.1) MLM=INT(1.5+RLU(0))
      IF(MOFSH.LE.2.OR.MOFSH.EQ.7) THEN
        NOFF=44
        PMMX=PMMO
      ELSE
        NOFF=40
        PMMX=VINT(1)
        IF(CKIN(2).GT.CKIN(1)) PMMX=MIN(CKIN(2),VINT(1))
      ENDIF
      MMED=0
      IF((KFMO.EQ.25.OR.KFMO.EQ.35.OR.KFMO.EQ.36).AND.MEQL.EQ.1.AND.
     &(KFD(1).EQ.23.OR.KFD(1).EQ.24)) MMED=1
      IF((KFMO.EQ.32.OR.IABS(KFMO).EQ.34).AND.(KFD(1).EQ.23.OR.
     &KFD(1).EQ.24).AND.(KFD(2).EQ.23.OR.KFD(2).EQ.24)) MMED=2
      IF((KFMO.EQ.32.OR.IABS(KFMO).EQ.34).AND.(KFD(2).EQ.25.OR.
     &KFD(2).EQ.35.OR.KFD(2).EQ.36)) MMED=3
      LOOP=1
 
C...Find where Breit-Wigners are required, else select discrete masses.
  100 DO 110 I=1,2
      KFCA=KFD(I)
      IF(KFCA.GT.100) KFCA=LUCOMP(KFCA)
      IF(KFCA.GT.0) THEN
        PMD(I)=PMAS(KFCA,1)
        PGD(I)=PMAS(KFCA,2)
      ELSE
        PMD(I)=0.
        PGD(I)=0.
      ENDIF
      IF(MSTP(42).LE.0.OR.PGD(I).LT.PARP(41)) THEN
        MBW(I)=0
        PMG(I)=PMD(I)
        RMG(I)=(PMG(I)/PMMX)**2
      ELSE
        MBW(I)=1
      ENDIF
  110 CONTINUE
 
C...Find allowed mass range and Breit-Wigner parameters.
      DO 120 I=1,2
      IF(MOFSH.EQ.1.AND.LOOP.EQ.1.AND.MBW(I).EQ.1) THEN
        PML(I)=PARP(42)
        PMU(I)=PMMX-PARP(42)
        IF(MBW(3-I).EQ.0) PMU(I)=MIN(PMU(I),PMMX-PMD(3-I))
        IF(PMU(I).LT.PML(I)+PARJ(64)) MBW(I)=-1
      ELSEIF((MBW(I).EQ.1.OR.MOFSH.GE.5).AND.MOFSH.NE.7) THEN
        ILM=I
        IF(MLM.EQ.2) ILM=3-I
        PML(I)=MAX(CKIN(NOFF+2*ILM-1),PARP(42))
        IF(MOFSH.GE.5.AND.I.EQ.2) PML(I)=MAX(PML(I),2.*PMAS(KFD2,1))
        PMU(I)=PMMX-MAX(CKIN(NOFF+5-2*ILM),PARP(42))
        IF(MOFSH.GE.5.AND.I.EQ.1) PMU(I)=MIN(PMU(I),PMMX-2.*
     &  PMAS(KFD2,1))
        IF(CKIN(NOFF+2*ILM).GT.CKIN(NOFF+2*ILM-1)) PMU(I)=MIN(PMU(I),
     &  CKIN(NOFF+2*ILM))
        IF(MBW(3-I).EQ.0) PMU(I)=MIN(PMU(I),PMMX-PMD(3-I))
        IF(I.EQ.MLM) PMU(I)=MIN(PMU(I),0.5*PMMX)
        IF(MEQL.EQ.0) PMH(I)=MIN(PMU(I),0.5*PMMX)
        IF(PMU(I).LT.PML(I)+PARJ(64)) MBW(I)=-1
        IF(MBW(I).EQ.1) THEN
          ATL(I)=ATAN((PML(I)**2-PMD(I)**2)/(PMD(I)*PGD(I)))
          ATU(I)=ATAN((PMU(I)**2-PMD(I)**2)/(PMD(I)*PGD(I)))
          IF(MEQL.EQ.0) ATH(I)=ATAN((PMH(I)**2-PMD(I)**2)/(PMD(I)*
     &    PGD(I)))
        ENDIF
      ELSEIF(MBW(I).EQ.1.AND.MOFSH.EQ.7) THEN
        ILM=I
        IF(MLM.EQ.2) ILM=3-I
        PML(I)=PARP(42)
        PMU(I)=PMMX-PARP(42)
        IF(MBW(3-I).EQ.0) PMU(I)=MIN(PMU(I),PMMX-PMD(3-I))
        IF(I.EQ.MLM) PMU(I)=MIN(PMU(I),0.5*PMMX)
        IF(MEQL.EQ.0) PMH(I)=MIN(PMU(I),0.5*PMMX)
        IF(PMU(I).LT.PML(I)+PARJ(64)) MBW(I)=-1
        IF(MBW(I).EQ.1) THEN
          ATL(I)=ATAN((PML(I)**2-PMD(I)**2)/(PMD(I)*PGD(I)))
          ATU(I)=ATAN((PMU(I)**2-PMD(I)**2)/(PMD(I)*PGD(I)))
          IF(MEQL.EQ.0) ATH(I)=ATAN((PMH(I)**2-PMD(I)**2)/(PMD(I)*
     &    PGD(I)))
        ENDIF
      ENDIF
  120 CONTINUE
      IF(MBW(1).LT.0.OR.MBW(2).LT.0.OR.(MBW(1).EQ.0.AND.MBW(2).EQ.0))
     &THEN
        CALL LUERRM(13,'(PYOFSH:) no allowed decay product masses')
        MINT(51)=1
        RETURN
      ENDIF
 
C...Calculation of partial width of resonance.
      IF(MOFSH.EQ.1) THEN
 
C..If only one integration, pick that to be the inner.
        IF(MBW(1).EQ.0) THEN
          PM2=PMD(1)
          PMD(1)=PMD(2)
          PGD(1)=PGD(2)
          PML(1)=PML(2)
          PMU(1)=PMU(2)
        ELSEIF(MBW(2).EQ.0) THEN
          PM2=PMD(2)
        ENDIF
 
C...Start outer loop of integration.
        IF(MBW(1).EQ.1.AND.MBW(2).EQ.1) THEN
          ATL2=ATAN((PML(2)**2-PMD(2)**2)/(PMD(2)*PGD(2)))
          ATU2=ATAN((PMU(2)**2-PMD(2)**2)/(PMD(2)*PGD(2)))
          NPT2=1
          XPT2(1)=1.
          INX2(1)=0
          FMAX2=0.
        ENDIF
  130   IF(MBW(1).EQ.1.AND.MBW(2).EQ.1) THEN
          PM2S=PMD(2)**2+PMD(2)*PGD(2)*TAN(ATL2+XPT2(NPT2)*(ATU2-ATL2))
          PM2=MIN(PMU(2),MAX(PML(2),SQRT(MAX(0.,PM2S))))
        ENDIF
        RM2=(PM2/PMMX)**2
 
C...Start inner loop of integration.
        PML1=PML(1)
        PMU1=MIN(PMU(1),PMMX-PM2)
        IF(MEQL.EQ.1) PMU1=MIN(PMU1,PM2)
        ATL1=ATAN((PML1**2-PMD(1)**2)/(PMD(1)*PGD(1)))
        ATU1=ATAN((PMU1**2-PMD(1)**2)/(PMD(1)*PGD(1)))
        IF(PML1+PARJ(64).GE.PMU1.OR.ATL1+1E-7.GE.ATU1) THEN
          FUNC2=0.
          GOTO 180
        ENDIF
        NPT1=1
        XPT1(1)=1.
        INX1(1)=0
        FMAX1=0.
  140   PM1S=PMD(1)**2+PMD(1)*PGD(1)*TAN(ATL1+XPT1(NPT1)*(ATU1-ATL1))
        PM1=MIN(PMU1,MAX(PML1,SQRT(MAX(0.,PM1S))))
        RM1=(PM1/PMMX)**2
 
C...Evaluate function value - inner loop.
        FUNC1=SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))
        IF(MMED.EQ.1) FUNC1=FUNC1*((1.-RM1-RM2)**2+8.*RM1*RM2)
        IF(MMED.EQ.2) FUNC1=FUNC1**3*(1.+10.*RM1+10.*RM2+RM1**2+
     &  RM2**2+10.*RM1*RM2)
        IF(FUNC1.GT.FMAX1) FMAX1=FUNC1
        FPT1(NPT1)=FUNC1
 
C...Go to next position in inner loop.
        IF(NPT1.EQ.1) THEN
          NPT1=NPT1+1
          XPT1(NPT1)=0.
          INX1(NPT1)=1
          GOTO 140
        ELSEIF(NPT1.LE.8) THEN
          NPT1=NPT1+1
          IF(NPT1.LE.4.OR.NPT1.EQ.6) ISH1=1
          ISH1=ISH1+1
          XPT1(NPT1)=0.5*(XPT1(ISH1)+XPT1(INX1(ISH1)))
          INX1(NPT1)=INX1(ISH1)
          INX1(ISH1)=NPT1
          GOTO 140
        ELSEIF(NPT1.LT.100) THEN
          ISN1=ISH1
  150     ISH1=ISH1+1
          IF(ISH1.GT.NPT1) ISH1=2
          IF(ISH1.EQ.ISN1) GOTO 160
          DFPT1=ABS(FPT1(ISH1)-FPT1(INX1(ISH1)))
          IF(DFPT1.LT.PARP(43)*FMAX1) GOTO 150
          NPT1=NPT1+1
          XPT1(NPT1)=0.5*(XPT1(ISH1)+XPT1(INX1(ISH1)))
          INX1(NPT1)=INX1(ISH1)
          INX1(ISH1)=NPT1
          GOTO 140
        ENDIF
 
C...Calculate integral over inner loop.
  160   FSUM1=0.
        DO 170 IPT1=2,NPT1
        FSUM1=FSUM1+0.5*(FPT1(IPT1)+FPT1(INX1(IPT1)))*
     &  (XPT1(INX1(IPT1))-XPT1(IPT1))
  170   CONTINUE
        FUNC2=FSUM1*(ATU1-ATL1)/PARU(1)
  180   IF(MBW(1).EQ.1.AND.MBW(2).EQ.1) THEN
          IF(FUNC2.GT.FMAX2) FMAX2=FUNC2
          FPT2(NPT2)=FUNC2
 
C...Go to next position in outer loop.
          IF(NPT2.EQ.1) THEN
            NPT2=NPT2+1
            XPT2(NPT2)=0.
            INX2(NPT2)=1
            GOTO 130
          ELSEIF(NPT2.LE.8) THEN
            NPT2=NPT2+1
            IF(NPT2.LE.4.OR.NPT2.EQ.6) ISH2=1
            ISH2=ISH2+1
            XPT2(NPT2)=0.5*(XPT2(ISH2)+XPT2(INX2(ISH2)))
            INX2(NPT2)=INX2(ISH2)
            INX2(ISH2)=NPT2
            GOTO 130
          ELSEIF(NPT2.LT.100) THEN
            ISN2=ISH2
  190       ISH2=ISH2+1
            IF(ISH2.GT.NPT2) ISH2=2
            IF(ISH2.EQ.ISN2) GOTO 200
            DFPT2=ABS(FPT2(ISH2)-FPT2(INX2(ISH2)))
            IF(DFPT2.LT.PARP(43)*FMAX2) GOTO 190
            NPT2=NPT2+1
            XPT2(NPT2)=0.5*(XPT2(ISH2)+XPT2(INX2(ISH2)))
            INX2(NPT2)=INX2(ISH2)
            INX2(ISH2)=NPT2
            GOTO 130
          ENDIF
 
C...Calculate integral over outer loop.
  200     FSUM2=0.
          DO 210 IPT2=2,NPT2
          FSUM2=FSUM2+0.5*(FPT2(IPT2)+FPT2(INX2(IPT2)))*
     &    (XPT2(INX2(IPT2))-XPT2(IPT2))
  210     CONTINUE
          FSUM2=FSUM2*(ATU2-ATL2)/PARU(1)
          IF(MEQL.EQ.1) FSUM2=2.*FSUM2
        ELSE
          FSUM2=FUNC2
        ENDIF
 
C...Save result; second integration for user-selected mass range.
        IF(LOOP.EQ.1) WIDW=FSUM2
        WID2=FSUM2
        IF(LOOP.EQ.1.AND.(CKIN(46).GE.CKIN(45).OR.CKIN(48).GE.CKIN(47)
     &  .OR.MAX(CKIN(45),CKIN(47)).GE.1.01*PARP(42))) THEN
          LOOP=2
          GOTO 100
        ENDIF
        RET1=WIDW
        RET2=WID2/WIDW
 
C...Select two decay product masses of a resonance.
      ELSEIF(MOFSH.EQ.2.OR.MOFSH.EQ.7) THEN
  220   DO 230 I=1,2
        IF(MBW(I).EQ.0) GOTO 230
        PMBW=PMD(I)**2+PMD(I)*PGD(I)*TAN(ATL(I)+RLU(0)*(ATU(I)-ATL(I)))
        PMG(I)=MIN(PMU(I),MAX(PML(I),SQRT(MAX(0.,PMBW))))
        RMG(I)=(PMG(I)/PMMX)**2
  230   CONTINUE
        IF((MEQL.EQ.1.AND.PMG(MAX(1,MLM)).GT.PMG(MIN(2,3-MLM))).OR.
     &  PMG(1)+PMG(2)+PARJ(64).GT.PMMX) GOTO 220
 
C...Weight with matrix element (if none known, use beta factor).
        FLAM=SQRT(MAX(0.,(1.-RMG(1)-RMG(2))**2-4.*RMG(1)*RMG(2)))
        IF(MMED.EQ.1) THEN
          WTBE=FLAM*((1.-RMG(1)-RMG(2))**2+8.*RMG(1)*RMG(2))
        ELSEIF(MMED.EQ.2) THEN
          WTBE=FLAM**3*(1.+10.*RMG(1)+10.*RMG(2)+RMG(1)**2+
     &    RMG(2)**2+10.*RMG(1)*RMG(2))
        ELSEIF(MMED.EQ.3) THEN
          WTBE=FLAM*(RMG(1)+FLAM**2/12.)
        ELSE
          WTBE=FLAM
        ENDIF
        IF(WTBE.LT.RLU(0)) GOTO 220
        RET1=PMG(1)
        RET2=PMG(2)
 
C...Find suitable set of masses for initialization of 2 -> 2 processes.
      ELSEIF(MOFSH.EQ.3) THEN
        IF(MBW(1).NE.0.AND.MBW(2).EQ.0) THEN
          PMG(1)=MIN(PMD(1),0.5*(PML(1)+PMU(1)))
          PMG(2)=PMD(2)
        ELSEIF(MBW(2).NE.0.AND.MBW(1).EQ.0) THEN
          PMG(1)=PMD(1)
          PMG(2)=MIN(PMD(2),0.5*(PML(2)+PMU(2)))
        ELSE
          IDIV=-1
  240     IDIV=IDIV+1
          PMG(1)=MIN(PMD(1),0.1*(IDIV*PML(1)+(10-IDIV)*PMU(1)))
          PMG(2)=MIN(PMD(2),0.1*(IDIV*PML(2)+(10-IDIV)*PMU(2)))
          IF(IDIV.LE.9.AND.PMG(1)+PMG(2).GT.0.9*PMMX) GOTO 240
        ENDIF
        RET1=PMG(1)
        RET2=PMG(2)
 
C...Evaluate importance of excluded tails of Breit-Wigners.
        IF(MEQL.EQ.0.AND.MBW(1).EQ.1.AND.MBW(2).EQ.1.AND.PMD(1)+PMD(2).
     &  GT.PMMX.AND.PMH(1).GT.PML(1).AND.PMH(2).GT.PML(2)) MEQL=2
        IF(MEQL.LE.1) THEN
          VINT(80)=1.
          DO 250 I=1,2
          IF(MBW(I).NE.0) VINT(80)=VINT(80)*1.25*(ATU(I)-ATL(I))/PARU(1)
  250     CONTINUE
        ELSE
          VINT(80)=(1.25/PARU(1))**2*MAX((ATU(1)-ATL(1))*
     &    (ATH(2)-ATL(2)),(ATH(1)-ATL(1))*(ATU(2)-ATL(2)))
        ENDIF
        IF((ISUB.EQ.15.OR.ISUB.EQ.19.OR.ISUB.EQ.30.OR.ISUB.EQ.35).AND.
     &  MSTP(43).NE.2) VINT(80)=2.*VINT(80)
        IF(ISUB.EQ.22.AND.MSTP(43).NE.2) VINT(80)=4.*VINT(80)
        IF(MEQL.GE.1) VINT(80)=2.*VINT(80)
 
C...Pick one particle to be the lighter (if improves efficiency).
      ELSEIF(MOFSH.EQ.4) THEN
        IF(MEQL.EQ.0.AND.MBW(1).EQ.1.AND.MBW(2).EQ.1.AND.PMD(1)+PMD(2).
     &  GT.PMMX.AND.PMH(1).GT.PML(1).AND.PMH(2).GT.PML(2)) MEQL=2
  260   IF(MEQL.EQ.2) MLM=INT(1.5+RLU(0))
 
C...Select two masses according to Breit-Wigner + flat in s + 1/s.
        DO 270 I=1,2
        IF(MBW(I).EQ.0) GOTO 270
        PMV=PMU(I)
        IF(MEQL.EQ.2.AND.I.EQ.MLM) PMV=PMH(I)
        ATV=ATU(I)
        IF(MEQL.EQ.2.AND.I.EQ.MLM) ATV=ATH(I)
        RBR=RLU(0)
        IF((ISUB.EQ.15.OR.ISUB.EQ.19.OR.ISUB.EQ.22.OR.ISUB.EQ.30.OR.
     &  ISUB.EQ.35).AND.MSTP(43).NE.2) RBR=2.*RBR
        IF(RBR.LT.0.8) THEN
          PMSR=PMD(I)**2+PMD(I)*PGD(I)*TAN(ATL(I)+RLU(0)*(ATV-ATL(I)))
          PMG(I)=MIN(PMV,MAX(PML(I),SQRT(MAX(0.,PMSR))))
        ELSEIF(RBR.LT.0.9) THEN
          PMG(I)=SQRT(MAX(0.,PML(I)**2+RLU(0)*(PMV**2-PML(I)**2)))
        ELSEIF(RBR.LT.1.5) THEN
          PMG(I)=PML(I)*(PMV/PML(I))**RLU(0)
        ELSE
          PMG(I)=SQRT(MAX(0.,PML(I)**2*PMV**2/(PML(I)**2+RLU(0)*
     &    (PMV**2-PML(I)**2))))
        ENDIF
  270   CONTINUE
        IF((MEQL.GE.1.AND.PMG(MAX(1,MLM)).GT.PMG(MIN(2,3-MLM))).OR.
     &  PMG(1)+PMG(2)+PARJ(64).GT.PMMX) THEN
          IF(MINT(48).EQ.1) THEN
            NGEN(0,1)=NGEN(0,1)+1
            NGEN(MINT(1),1)=NGEN(MINT(1),1)+1
            GOTO 260
          ELSE
            MINT(51)=1
            RETURN
          ENDIF
        ENDIF
        RET1=PMG(1)
        RET2=PMG(2)
 
C...Give weight for selected mass distribution.
        VINT(80)=1.
        DO 280 I=1,2
        IF(MBW(I).EQ.0) GOTO 280
        PMV=PMU(I)
        IF(MEQL.EQ.2.AND.I.EQ.MLM) PMV=PMH(I)
        ATV=ATU(I)
        IF(MEQL.EQ.2.AND.I.EQ.MLM) ATV=ATH(I)
        F0=PMD(I)*PGD(I)/((PMG(I)**2-PMD(I)**2)**2+
     &  (PMD(I)*PGD(I))**2)/PARU(1)
        F1=1.
        F2=1./PMG(I)**2
        F3=1./PMG(I)**4
        FI0=(ATV-ATL(I))/PARU(1)
        FI1=PMV**2-PML(I)**2
        FI2=2.*LOG(PMV/PML(I))
        FI3=1./PML(I)**2-1./PMV**2
        IF((ISUB.EQ.15.OR.ISUB.EQ.19.OR.ISUB.EQ.22.OR.ISUB.EQ.30.OR.
     &  ISUB.EQ.35).AND.MSTP(43).NE.2) THEN
          VINT(80)=VINT(80)*20./(8.+(FI0/F0)*(F1/FI1+6.*F2/FI2+
     &    5.*F3/FI3))
        ELSE
          VINT(80)=VINT(80)*10./(8.+(FI0/F0)*(F1/FI1+F2/FI2))
        ENDIF
        VINT(80)=VINT(80)*FI0
  280   CONTINUE
        IF(MEQL.GE.1) VINT(80)=2.*VINT(80)
 
      ELSEIF(MOFSH.EQ.5) THEN
C...Find suitable set of masses for initialization of 2 -> 3 process.
        IDIV=6
  290   IDIV=IDIV-1
        IF(MBW(1).EQ.0) THEN
          PMG(1)=PMD(1)
        ELSE
          PMSR=PMD(1)**2+PMD(1)*PGD(1)*TAN(ATL(1)+0.1*IDIV*(ATU(1)-
     &    ATL(1)))
          PMG(1)=MIN(PMU(1),MAX(PML(1),SQRT(MAX(0.,PMSR))))
        ENDIF
        PMG(2)=PML(2)*(PMU(2)/PML(2))**(0.1*IDIV)
        IF(IDIV.GE.1.AND.PMG(1)+PMG(2).GT.0.9*PMMX) GOTO 290
        RET1=PMG(1)
        RET2=PMG(2)
 
C...Evaluate size of selected phase space volume.
        VINT(80)=2.*LOG(PMU(2)/PML(2))
        IF(MBW(1).NE.0) VINT(80)=VINT(80)*1.25*(ATU(1)-ATL(1))/PARU(1)
 
C...Pick decay angles.
        VINT(81)=0.
        VINT(82)=0.5*PARU(1)
        VINT(83)=1.
        VINT(84)=0.
 
C...Select flavour of resonance decays.
        KFA=KFPR(ISUB,1)
        CALL PYWIDT(KFA,PMG(1)**2,WDTP,WDTE)
        IF(KCHG(KFA,3).EQ.0) THEN
          IPM=2
        ELSE
          IPM=(5-ISIGN(1,KFA))/2
        ENDIF
        WDTE0S=WDTE(0,1)+WDTE(0,IPM)+WDTE(0,4)
        IF(WDTE0S.LE.0.) THEN
          CALL LUERRM(12,'(PYOFSH:) no allowed resonace decay channel')
          MINT(51)=1
          RETURN
        ENDIF
        WDTEC=0.
        DO 300 IDL=1,MDCY(KFA,3)
        WDTEK=WDTE(IDL,1)+WDTE(IDL,IPM)+WDTE(IDL,4)
        IF(WDTEK.GT.WDTEC) THEN
          IDC=IDL+MDCY(KFA,2)-1
          WDTEC=WDTEK
        ENDIF
  300   CONTINUE
        MINT(35)=IDC
 
C...Compensating factor for all flavours.
        KFL=IABS(KFDP(IDC,1))
        QFL=KCHG(KFL,1)/3.
        AFL=SIGN(1.,QFL+0.1)
        VFL=AFL-4.*PARU(102)*QFL
        WDTEK=VFL**2+AFL**2
        VINT(80)=VINT(80)*WDTE0S/WDTEK
 
      ELSEIF(MOFSH.EQ.6) THEN
C...Select two masses, one basically Breit-Wigner, other dm^2/m^2.
        IF(MBW(1).NE.0) THEN
          RBR=RLU(0)
          IF(RBR.LT.0.8) THEN
            PMSR=PMD(1)**2+PMD(1)*PGD(1)*TAN(ATL(1)+RLU(0)*
     &      (ATU(1)-ATL(1)))
            PMG(1)=MIN(PMU(1),MAX(PML(1),SQRT(MAX(0.,PMSR))))
          ELSEIF(RBR.LT.0.9) THEN
            PMG(1)=SQRT(MAX(0.,PML(1)**2+RLU(0)*(PMU(1)**2-PML(1)**2)))
          ELSE
            PMG(1)=PML(1)*(PMU(1)/PML(1))**RLU(0)
          ENDIF
        ENDIF
        PMG(2)=PML(2)*(PMU(2)/PML(2))**RLU(0)
        IF(SQRT(MAX(0.,1.-(PML(2)/PMG(2))**2)).LT.RLU(0).OR.
     &  PMG(1)+PMG(2)+PARJ(64).GT.PMMX) THEN
          MINT(51)=1
          RETURN
        ENDIF
        RET1=PMG(1)
        RET2=PMG(2)
 
C...Give weight for selected mass distribution.
        VINT(80)=2.*LOG(PMU(2)/PML(2))
        IF(MBW(1).NE.0) THEN
          F0=PMD(1)*PGD(1)/((PMG(1)**2-PMD(1)**2)**2+
     &    (PMD(1)*PGD(1))**2)/PARU(1)
          F1=1.
          F2=1./PMG(1)**2
          FI0=(ATU(1)-ATL(1))/PARU(1)
          FI1=PMU(1)**2-PML(1)**2
          FI2=2.*LOG(PMU(1)/PML(1))
          VINT(80)=VINT(80)*10.*FI0/(8.+(FI0/F0)*(F1/FI1+F2/FI2))
        ENDIF
 
C...Select decay angles.
        VINT(81)=2.*RLU(0)-1.
        VINT(82)=PARU(2)*RLU(0)
        VINT(83)=2.*RLU(0)-1.
        VINT(84)=PARU(2)*RLU(0)
 
C...Select flavour of resonance decays.
        KFA=KFPR(ISUB,1)
        CALL PYWIDT(KFA,PMG(1)**2,WDTP,WDTE)
        IF(KCHG(KFA,3).EQ.0) THEN
          IPM=2
        ELSE
          IPM=(5-ISIGN(1,KFA))/2
        ENDIF
        WDTE0S=WDTE(0,1)+WDTE(0,IPM)+WDTE(0,4)
        IF(WDTE0S.LE.0.) THEN
          CALL LUERRM(12,'(PYOFSH:) no allowed resonace decay channel')
          MINT(51)=1
          RETURN
        ENDIF
        RKFL=WDTE0S*RLU(0)
        IDL=0
  310   IDL=IDL+1
        IDC=IDL+MDCY(KFA,2)-1
        RKFL=RKFL-(WDTE(IDL,1)+WDTE(IDL,IPM)+WDTE(IDL,4))
        IF(IDL.LT.MDCY(KFA,3).AND.RKFL.GT.0.) GOTO 310
        MINT(35)=IDC
 
C...Compensating factor for all flavours.
        KFL=IABS(KFDP(IDC,1))
        QFL=KCHG(KFL,1)/3.
        AFL=SIGN(1.,QFL+0.1)
        VFL=AFL-4.*PARU(102)*QFL
        WDTEK=VFL**2+AFL**2
        VINT(80)=VINT(80)*WDTE0S/WDTEK
      ENDIF
 
      RETURN
      END
