*CMZ :          17/07/98  15.44.32  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUDECY(IP)

C...Purpose: to handle the decay of unstable particles.
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEEP,LUDAT3.
      COMMON /LUDAT3/ MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      SAVE /LUDAT3/
*KEND.
      DIMENSION VDCY(4),KFLO(4),KFL1(4),PV(10,5),RORD(10),UE(3),BE(3),
     &WTCOR(10)
      DATA WTCOR/2.,5.,15.,60.,250.,1500.,1.2E4,1.2E5,150.,16./

C...Functions: momentum in two-particle decays, four-product and
C...matrix element times phase space in weak decays.
      PAWT(A,B,C)=SQRT((A**2-(B+C)**2)*(A**2-(B-C)**2))/(2.*A)
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)
      HMEPS(HA)=((1.-HRQ-HA)**2+3.*HA*(1.+HRQ-HA))*
     &SQRT((1.-HRQ-HA)**2-4.*HRQ*HA)

C...Initial values.
      NTRY=0
      NSAV=N
      KFA=IABS(K(IP,2))
      KFS=ISIGN(1,K(IP,2))
      KC=LUCOMP(KFA)
      MSTJ(92)=0

C...Choose lifetime and determine decay vertex.
      IF(K(IP,1).EQ.5) THEN
        V(IP,5)=0.
      ELSEIF(K(IP,1).NE.4) THEN
        V(IP,5)=-PMAS(KC,4)*LOG(RLU(0))
      ENDIF
      DO 100 J=1,4
  100 VDCY(J)=V(IP,J)+V(IP,5)*P(IP,J)/P(IP,5)

C...Determine whether decay allowed or not.
      MOUT=0
      IF(MSTJ(22).EQ.2) THEN
        IF(PMAS(KC,4).GT.PARJ(71)) MOUT=1
      ELSEIF(MSTJ(22).EQ.3) THEN
        IF(VDCY(1)**2+VDCY(2)**2+VDCY(3)**2.GT.PARJ(72)**2) MOUT=1
      ELSEIF(MSTJ(22).EQ.4) THEN
        IF(VDCY(1)**2+VDCY(2)**2.GT.PARJ(73)**2) MOUT=1
        IF(ABS(VDCY(3)).GT.PARJ(74)) MOUT=1
      ENDIF
      IF(MOUT.EQ.1.AND.K(IP,1).NE.5) THEN
        K(IP,1)=4
        RETURN
      ENDIF

C...Check existence of decay channels. Particle/antiparticle rules.
      KCA=KC
      IF(MDCY(KC,2).GT.0) THEN
        MDMDCY=MDME(MDCY(KC,2),2)
        IF(MDMDCY.GT.80.AND.MDMDCY.LE.90) KCA=MDMDCY
      ENDIF
      IF(MDCY(KCA,2).LE.0.OR.MDCY(KCA,3).LE.0) THEN
        CALL LUERRM(9,'(LUDECY:) no decay channel defined')
        RETURN
      ENDIF
      IF(MOD(KFA/1000,10).EQ.0.AND.(KCA.EQ.85.OR.KCA.EQ.87)) KFS=-KFS
      IF(KCHG(KC,3).EQ.0) THEN
        KFSP=1
        KFSN=0
        IF(RLU(0).GT.0.5) KFS=-KFS
      ELSEIF(KFS.GT.0) THEN
        KFSP=1
        KFSN=0
      ELSE
        KFSP=0
        KFSN=1
      ENDIF

C...Sum branching ratios of allowed decay channels.
  110 NOPE=0
      BRSU=0.
      DO 120 IDL=MDCY(KCA,2),MDCY(KCA,2)+MDCY(KCA,3)-1
      IF(MDME(IDL,1).NE.1.AND.KFSP*MDME(IDL,1).NE.2.AND.
     &KFSN*MDME(IDL,1).NE.3) GOTO 120
      IF(MDME(IDL,2).GT.100) GOTO 120
      NOPE=NOPE+1
      BRSU=BRSU+BRAT(IDL)
  120 CONTINUE
      IF(NOPE.EQ.0) THEN
        CALL LUERRM(2,'(LUDECY:) all decay channels closed by user')
        RETURN
      ENDIF

C...Select decay channel among allowed ones.
  130 RBR=BRSU*RLU(0)
      IDL=MDCY(KCA,2)-1
  140 IDL=IDL+1
      IF(MDME(IDL,1).NE.1.AND.KFSP*MDME(IDL,1).NE.2.AND.
     &KFSN*MDME(IDL,1).NE.3) THEN
        IF(IDL.LT.MDCY(KCA,2)+MDCY(KCA,3)-1) GOTO 140
      ELSEIF(MDME(IDL,2).GT.100) THEN
        IF(IDL.LT.MDCY(KCA,2)+MDCY(KCA,3)-1) GOTO 140
      ELSE
        IDC=IDL
        RBR=RBR-BRAT(IDL)
        IF(IDL.LT.MDCY(KCA,2)+MDCY(KCA,3)-1.AND.RBR.GT.0.) GOTO 140
      ENDIF

C...Start readout of decay channel: matrix element, reset counters.
      MMAT=MDME(IDC,2)
  150 NTRY=NTRY+1
      IF(NTRY.GT.1000) THEN
        CALL LUERRM(14,'(LUDECY:) caught in infinite loop')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      I=N
      NP=0
      NQ=0
      MBST=0
      IF(MMAT.GE.11.AND.MMAT.NE.46.AND.P(IP,4).GT.20.*P(IP,5)) MBST=1
      DO 160 J=1,4
      PV(1,J)=0.
  160 IF(MBST.EQ.0) PV(1,J)=P(IP,J)
      IF(MBST.EQ.1) PV(1,4)=P(IP,5)
      PV(1,5)=P(IP,5)
      PS=0.
      PSQ=0.
      MREM=0

C...Read out decay products. Convert to standard flavour code.
      JTMAX=5
      IF(MDME(IDC+1,2).EQ.101) JTMAX=10
      DO 170 JT=1,JTMAX
      IF(JT.LE.5) KP=KFDP(IDC,JT)
      IF(JT.GE.6) KP=KFDP(IDC+1,JT-5)
      IF(KP.EQ.0) GOTO 170
      KPA=IABS(KP)
      KCP=LUCOMP(KPA)
      IF(KCHG(KCP,3).EQ.0.AND.KPA.NE.81.AND.KPA.NE.82) THEN
        KFP=KP
      ELSEIF(KPA.NE.81.AND.KPA.NE.82) THEN
        KFP=KFS*KP
      ELSEIF(KPA.EQ.81.AND.MOD(KFA/1000,10).EQ.0) THEN
        KFP=-KFS*MOD(KFA/10,10)
      ELSEIF(KPA.EQ.81.AND.MOD(KFA/100,10).GE.MOD(KFA/10,10)) THEN
        KFP=KFS*(100*MOD(KFA/10,100)+3)
      ELSEIF(KPA.EQ.81) THEN
        KFP=KFS*(1000*MOD(KFA/10,10)+100*MOD(KFA/100,10)+1)
      ELSEIF(KP.EQ.82) THEN
        CALL LUKFDI(-KFS*INT(1.+(2.+PARJ(2))*RLU(0)),0,KFP,KDUMP)
        IF(KFP.EQ.0) GOTO 150
        MSTJ(93)=1
        IF(PV(1,5).LT.PARJ(32)+2.*ULMASS(KFP)) GOTO 150
      ELSEIF(KP.EQ.-82) THEN
        KFP=-KFP
        IF(IABS(KFP).GT.10) KFP=KFP+ISIGN(10000,KFP)
      ENDIF
      IF(KPA.EQ.81.OR.KPA.EQ.82) KCP=LUCOMP(KFP)

C...Add decay product to event record or to quark flavour list.
      KFPA=IABS(KFP)
      KQP=KCHG(KCP,2)
      IF(MMAT.GE.11.AND.MMAT.LE.30.AND.KQP.NE.0) THEN
        NQ=NQ+1
        KFLO(NQ)=KFP
        MSTJ(93)=2
        PSQ=PSQ+ULMASS(KFLO(NQ))
      ELSEIF(MMAT.GE.42.AND.MMAT.LE.43.AND.NP.EQ.3.AND.MOD(NQ,2).EQ.1)
     &THEN
        NQ=NQ-1
        PS=PS-P(I,5)
        K(I,1)=1
        KFI=K(I,2)
        CALL LUKFDI(KFP,KFI,KFLDMP,K(I,2))
        IF(K(I,2).EQ.0) GOTO 150
        MSTJ(93)=1
        P(I,5)=ULMASS(K(I,2))
        PS=PS+P(I,5)
      ELSE
        I=I+1
        NP=NP+1
        IF(MMAT.NE.33.AND.KQP.NE.0) NQ=NQ+1
        IF(MMAT.EQ.33.AND.KQP.NE.0.AND.KQP.NE.2) NQ=NQ+1
        K(I,1)=1+MOD(NQ,2)
        IF(MMAT.EQ.4.AND.JT.LE.2.AND.KFP.EQ.21) K(I,1)=2
        IF(MMAT.EQ.4.AND.JT.EQ.3) K(I,1)=1
        K(I,2)=KFP
        K(I,3)=IP
        K(I,4)=0
        K(I,5)=0
        P(I,5)=ULMASS(KFP)
        IF(MMAT.EQ.45.AND.KFPA.EQ.89) P(I,5)=PARJ(32)
        PS=PS+P(I,5)
      ENDIF
  170 CONTINUE

C...Choose decay multiplicity in phase space model.
  180 IF(MMAT.GE.11.AND.MMAT.LE.30) THEN
        PSP=PS
        CNDE=PARJ(61)*LOG(MAX((PV(1,5)-PS-PSQ)/PARJ(62),1.1))
        IF(MMAT.EQ.12) CNDE=CNDE+PARJ(63)
  190   NTRY=NTRY+1
        IF(NTRY.GT.1000) THEN
          CALL LUERRM(14,'(LUDECY:) caught in infinite loop')
          IF(MSTU(21).GE.1) RETURN
        ENDIF
        IF(MMAT.LE.20) THEN
          GAUSS=SQRT(-2.*CNDE*LOG(MAX(1E-10,RLU(0))))*
     &    SIN(PARU(2)*RLU(0))
          ND=0.5+0.5*NP+0.25*NQ+CNDE+GAUSS
          IF(ND.LT.NP+NQ/2.OR.ND.LT.2.OR.ND.GT.10) GOTO 190
          IF(MMAT.EQ.13.AND.ND.EQ.2) GOTO 190
          IF(MMAT.EQ.14.AND.ND.LE.3) GOTO 190
          IF(MMAT.EQ.15.AND.ND.LE.4) GOTO 190
        ELSE
          ND=MMAT-20
        ENDIF

C...Form hadrons from flavour content.
        DO 200 JT=1,4
  200   KFL1(JT)=KFLO(JT)
        IF(ND.EQ.NP+NQ/2) GOTO 220
        DO 210 I=N+NP+1,N+ND-NQ/2
        JT=1+INT((NQ-1)*RLU(0))
        CALL LUKFDI(KFL1(JT),0,KFL2,K(I,2))
        IF(K(I,2).EQ.0) GOTO 190
  210   KFL1(JT)=-KFL2
  220   JT=2
        JT2=3
        JT3=4
        IF(NQ.EQ.4.AND.RLU(0).LT.PARJ(66)) JT=4
        IF(JT.EQ.4.AND.ISIGN(1,KFL1(1)*(10-IABS(KFL1(1))))*
     &  ISIGN(1,KFL1(JT)*(10-IABS(KFL1(JT)))).GT.0) JT=3
        IF(JT.EQ.3) JT2=2
        IF(JT.EQ.4) JT3=2
        CALL LUKFDI(KFL1(1),KFL1(JT),KFLDMP,K(N+ND-NQ/2+1,2))
        IF(K(N+ND-NQ/2+1,2).EQ.0) GOTO 190
        IF(NQ.EQ.4) CALL LUKFDI(KFL1(JT2),KFL1(JT3),KFLDMP,K(N+ND,2))
        IF(NQ.EQ.4.AND.K(N+ND,2).EQ.0) GOTO 190

C...Check that sum of decay product masses not too large.
        PS=PSP
        DO 230 I=N+NP+1,N+ND
        K(I,1)=1
        K(I,3)=IP
        K(I,4)=0
        K(I,5)=0
        P(I,5)=ULMASS(K(I,2))
  230   PS=PS+P(I,5)
        IF(PS+PARJ(64).GT.PV(1,5)) GOTO 190

C...Rescale energy to subtract off spectator quark mass.
      ELSEIF((MMAT.EQ.31.OR.MMAT.EQ.33.OR.MMAT.EQ.44.OR.MMAT.EQ.45).
     &AND.NP.GE.3) THEN
        PS=PS-P(N+NP,5)
        PQT=(P(N+NP,5)+PARJ(65))/PV(1,5)
        DO 240 J=1,5
        P(N+NP,J)=PQT*PV(1,J)
  240   PV(1,J)=(1.-PQT)*PV(1,J)
        IF(PS+PARJ(64).GT.PV(1,5)) GOTO 150
        ND=NP-1
        MREM=1

C...Phase space factors imposed in W decay.
      ELSEIF(MMAT.EQ.46) THEN
        MSTJ(93)=1
        PSMC=ULMASS(K(N+1,2))
        MSTJ(93)=1
        PSMC=PSMC+ULMASS(K(N+2,2))
        IF(MAX(PS,PSMC)+PARJ(32).GT.PV(1,5)) GOTO 130
        HR1=(P(N+1,5)/PV(1,5))**2
        HR2=(P(N+2,5)/PV(1,5))**2
        IF((1.-HR1-HR2)*(2.+HR1+HR2)*SQRT((1.-HR1-HR2)**2-4.*HR1*HR2).
     &  LT.2.*RLU(0)) GOTO 130
        ND=NP

C...Fully specified final state: check mass broadening effects.
      ELSE
        IF(NP.GE.2.AND.PS+PARJ(64).GT.PV(1,5)) GOTO 150
        ND=NP
      ENDIF

C...Select W mass in decay Q -> W + q, without W propagator.
      IF(MMAT.EQ.45.AND.MSTJ(25).LE.0) THEN
        HLQ=(PARJ(32)/PV(1,5))**2
        HUQ=(1.-(P(N+2,5)+PARJ(64))/PV(1,5))**2
        HRQ=(P(N+2,5)/PV(1,5))**2
  250   HW=HLQ+RLU(0)*(HUQ-HLQ)
        IF(HMEPS(HW).LT.RLU(0)) GOTO 250
        P(N+1,5)=PV(1,5)*SQRT(HW)

C...Ditto, including W propagator. Divide mass range into three regions.
      ELSEIF(MMAT.EQ.45) THEN
        HQW=(PV(1,5)/PMAS(24,1))**2
        HLW=(PARJ(32)/PMAS(24,1))**2
        HUW=((PV(1,5)-P(N+2,5)-PARJ(64))/PMAS(24,1))**2
        HRQ=(P(N+2,5)/PV(1,5))**2
        HG=PMAS(24,2)/PMAS(24,1)
        HATL=ATAN((HLW-1.)/HG)
        HM=MIN(1.,HUW-0.001)
        HMV1=HMEPS(HM/HQW)/((HM-1.)**2+HG**2)
  260   HM=HM-HG
        HMV2=HMEPS(HM/HQW)/((HM-1.)**2+HG**2)
        HSAV1=HMEPS(HM/HQW)
        HSAV2=1./((HM-1.)**2+HG**2)
        IF(HMV2.GT.HMV1.AND.HM-HG.GT.HLW) THEN
          HMV1=HMV2
          GOTO 260
        ENDIF
        HMV=MIN(2.*HMV1,HMEPS(HM/HQW)/HG**2)
        HM1=1.-SQRT(1./HMV-HG**2)
        IF(HM1.GT.HLW.AND.HM1.LT.HM) THEN
          HM=HM1
        ELSEIF(HMV2.LE.HMV1) THEN
          HM=MAX(HLW,HM-MIN(0.1,1.-HM))
        ENDIF
        HATM=ATAN((HM-1.)/HG)
        HWT1=(HATM-HATL)/HG
        HWT2=HMV*(MIN(1.,HUW)-HM)
        HWT3=0.
        IF(HUW.GT.1.) THEN
          HATU=ATAN((HUW-1.)/HG)
          HMP1=HMEPS(1./HQW)
          HWT3=HMP1*HATU/HG
        ENDIF

C...Select mass region and W mass there. Accept according to weight.
  270   HREG=RLU(0)*(HWT1+HWT2+HWT3)
        IF(HREG.LE.HWT1) THEN
          HW=1.+HG*TAN(HATL+RLU(0)*(HATM-HATL))
          HACC=HMEPS(HW/HQW)
        ELSEIF(HREG.LE.HWT1+HWT2) THEN
          HW=HM+RLU(0)*(MIN(1.,HUW)-HM)
          HACC=HMEPS(HW/HQW)/((HW-1.)**2+HG**2)/HMV
        ELSE
          HW=1.+HG*TAN(RLU(0)*HATU)
          HACC=HMEPS(HW/HQW)/HMP1
        ENDIF
        IF(HACC.LT.RLU(0)) GOTO 270
        P(N+1,5)=PMAS(24,1)*SQRT(HW)
      ENDIF

C...Determine position of grandmother, number of sisters, Q -> W sign.
      NM=0
      KFAS=0
      MSGN=0
      IF(MMAT.EQ.3.OR.MMAT.EQ.46) THEN
        IM=K(IP,3)
        IF(IM.LT.0.OR.IM.GE.IP) IM=0
        IF(IM.NE.0) KFAM=IABS(K(IM,2))
        IF(IM.NE.0.AND.MMAT.EQ.3) THEN
          DO 280 IL=MAX(IP-2,IM+1),MIN(IP+2,N)
          IF(K(IL,3).EQ.IM) NM=NM+1
  280     IF(K(IL,3).EQ.IM.AND.IL.NE.IP) ISIS=IL
          IF(NM.NE.2.OR.KFAM.LE.100.OR.MOD(KFAM,10).NE.1.OR.
     &    MOD(KFAM/1000,10).NE.0) NM=0
          IF(NM.EQ.2) THEN
            KFAS=IABS(K(ISIS,2))
            IF((KFAS.LE.100.OR.MOD(KFAS,10).NE.1.OR.
     &      MOD(KFAS/1000,10).NE.0).AND.KFAS.NE.22) NM=0
          ENDIF
        ELSEIF(IM.NE.0.AND.MMAT.EQ.46) THEN
          MSGN=ISIGN(1,K(IM,2)*K(IP,2))
          IF(KFAM.GT.100.AND.MOD(KFAM/1000,10).EQ.0) MSGN=
     &    MSGN*(-1)**MOD(KFAM/100,10)
        ENDIF
      ENDIF

C...Kinematics of one-particle decays.
      IF(ND.EQ.1) THEN
        DO 290 J=1,4
  290   P(N+1,J)=P(IP,J)
        GOTO 510
      ENDIF

C...Calculate maximum weight ND-particle decay.
      PV(ND,5)=P(N+ND,5)
      IF(ND.GE.3) THEN
        WTMAX=1./WTCOR(ND-2)
        PMAX=PV(1,5)-PS+P(N+ND,5)
        PMIN=0.
        DO 300 IL=ND-1,1,-1
        PMAX=PMAX+P(N+IL,5)
        PMIN=PMIN+P(N+IL+1,5)
  300   WTMAX=WTMAX*PAWT(PMAX,PMIN,P(N+IL,5))
      ENDIF

C...Find virtual gamma mass in Dalitz decay.
  310 IF(ND.EQ.2) THEN
      ELSEIF(MMAT.EQ.2) THEN
        PMES=4.*PMAS(11,1)**2
        PMRHO2=PMAS(131,1)**2
        PGRHO2=PMAS(131,2)**2
  320   PMST=PMES*(P(IP,5)**2/PMES)**RLU(0)
        WT=(1+0.5*PMES/PMST)*SQRT(MAX(0.,1.-PMES/PMST))*
     &  (1.-PMST/P(IP,5)**2)**3*(1.+PGRHO2/PMRHO2)/
     &  ((1.-PMST/PMRHO2)**2+PGRHO2/PMRHO2)
        IF(WT.LT.RLU(0)) GOTO 320
        PV(2,5)=MAX(2.00001*PMAS(11,1),SQRT(PMST))

C...M-generator gives weight. If rejected, try again.
      ELSE
  330   RORD(1)=1.
        DO 350 IL1=2,ND-1
        RSAV=RLU(0)
        DO 340 IL2=IL1-1,1,-1
        IF(RSAV.LE.RORD(IL2)) GOTO 350
  340   RORD(IL2+1)=RORD(IL2)
  350   RORD(IL2+1)=RSAV
        RORD(ND)=0.
        WT=1.
        DO 360 IL=ND-1,1,-1
        PV(IL,5)=PV(IL+1,5)+P(N+IL,5)+(RORD(IL)-RORD(IL+1))*(PV(1,5)-PS)
  360   WT=WT*PAWT(PV(IL,5),PV(IL+1,5),P(N+IL,5))
        IF(WT.LT.RLU(0)*WTMAX) GOTO 330
      ENDIF

C...Perform two-particle decays in respective CM frame.
  370 DO 390 IL=1,ND-1
      PA=PAWT(PV(IL,5),PV(IL+1,5),P(N+IL,5))
      UE(3)=2.*RLU(0)-1.
      PHI=PARU(2)*RLU(0)
      UE(1)=SQRT(1.-UE(3)**2)*COS(PHI)
      UE(2)=SQRT(1.-UE(3)**2)*SIN(PHI)
      DO 380 J=1,3
      P(N+IL,J)=PA*UE(J)
  380 PV(IL+1,J)=-PA*UE(J)
      P(N+IL,4)=SQRT(PA**2+P(N+IL,5)**2)
  390 PV(IL+1,4)=SQRT(PA**2+PV(IL+1,5)**2)

C...Lorentz transform decay products to lab frame.
      DO 400 J=1,4
  400 P(N+ND,J)=PV(ND,J)
      DO 430 IL=ND-1,1,-1
      DO 410 J=1,3
  410 BE(J)=PV(IL,J)/PV(IL,4)
      GA=PV(IL,4)/PV(IL,5)
      DO 430 I=N+IL,N+ND
      BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
      DO 420 J=1,3
  420 P(I,J)=P(I,J)+GA*(GA*BEP/(1.+GA)+P(I,4))*BE(J)
  430 P(I,4)=GA*(P(I,4)+BEP)

C...Matrix elements for omega and phi decays.
      IF(MMAT.EQ.1) THEN
        WT=(P(N+1,5)*P(N+2,5)*P(N+3,5))**2-(P(N+1,5)*FOUR(N+2,N+3))**2
     &  -(P(N+2,5)*FOUR(N+1,N+3))**2-(P(N+3,5)*FOUR(N+1,N+2))**2
     &  +2.*FOUR(N+1,N+2)*FOUR(N+1,N+3)*FOUR(N+2,N+3)
        IF(MAX(WT*WTCOR(9)/P(IP,5)**6,0.001).LT.RLU(0)) GOTO 310

C...Matrix elements for pi0 or eta Dalitz decay to gamma e+ e-.
      ELSEIF(MMAT.EQ.2) THEN
        FOUR12=FOUR(N+1,N+2)
        FOUR13=FOUR(N+1,N+3)
        FOUR23=0.5*PMST-0.25*PMES
        WT=(PMST-0.5*PMES)*(FOUR12**2+FOUR13**2)+
     &  PMES*(FOUR12*FOUR13+FOUR12**2+FOUR13**2)
        IF(WT.LT.RLU(0)*0.25*PMST*(P(IP,5)**2-PMST)**2) GOTO 370

C...Matrix element for S0 -> S1 + V1 -> S1 + S2 + S3 (S scalar,
C...V vector), of form cos**2(theta02) in V1 rest frame, and for
C...S0 -> gamma + V1 -> gamma + S2 + S3, of form sin**2(theta02).
      ELSEIF(MMAT.EQ.3.AND.NM.EQ.2) THEN
        FOUR10=FOUR(IP,IM)
        FOUR12=FOUR(IP,N+1)
        FOUR02=FOUR(IM,N+1)
        PMS1=P(IP,5)**2
        PMS0=P(IM,5)**2
        PMS2=P(N+1,5)**2
        IF(KFAS.NE.22) HNUM=(FOUR10*FOUR12-PMS1*FOUR02)**2
        IF(KFAS.EQ.22) HNUM=PMS1*(2.*FOUR10*FOUR12*FOUR02-
     &  PMS1*FOUR02**2-PMS0*FOUR12**2-PMS2*FOUR10**2+PMS1*PMS0*PMS2)
        HNUM=MAX(1E-6*PMS1**2*PMS0*PMS2,HNUM)
        HDEN=(FOUR10**2-PMS1*PMS0)*(FOUR12**2-PMS1*PMS2)
        IF(HNUM.LT.RLU(0)*HDEN) GOTO 370

C...Matrix element for "onium" -> g + g + g or gamma + g + g.
      ELSEIF(MMAT.EQ.4) THEN
        HX1=2.*FOUR(IP,N+1)/P(IP,5)**2
        HX2=2.*FOUR(IP,N+2)/P(IP,5)**2
        HX3=2.*FOUR(IP,N+3)/P(IP,5)**2
        WT=((1.-HX1)/(HX2*HX3))**2+((1.-HX2)/(HX1*HX3))**2+
     &  ((1.-HX3)/(HX1*HX2))**2
        IF(WT.LT.2.*RLU(0)) GOTO 310
        IF(K(IP+1,2).EQ.22.AND.(1.-HX1)*P(IP,5)**2.LT.4.*PARJ(32)**2)
     &  GOTO 310

C...Effective matrix element for nu spectrum in tau -> nu + hadrons.
      ELSEIF(MMAT.EQ.41) THEN
        HX1=2.*FOUR(IP,N+1)/P(IP,5)**2
        IF(8.*HX1*(3.-2.*HX1)/9..LT.RLU(0)) GOTO 310

C...Matrix elements for weak decays (only semileptonic for c and b)
      ELSEIF(MMAT.GE.42.AND.MMAT.LE.44.AND.ND.EQ.3) THEN
        IF(MBST.EQ.0) WT=FOUR(IP,N+1)*FOUR(N+2,N+3)
        IF(MBST.EQ.1) WT=P(IP,5)*P(N+1,4)*FOUR(N+2,N+3)
        IF(WT.LT.RLU(0)*P(IP,5)*PV(1,5)**3/WTCOR(10)) GOTO 310
      ELSEIF(MMAT.GE.42.AND.MMAT.LE.44) THEN
        DO 440 J=1,4
        P(N+NP+1,J)=0.
        DO 440 IS=N+3,N+NP
  440   P(N+NP+1,J)=P(N+NP+1,J)+P(IS,J)
        IF(MBST.EQ.0) WT=FOUR(IP,N+1)*FOUR(N+2,N+NP+1)
        IF(MBST.EQ.1) WT=P(IP,5)*P(N+1,4)*FOUR(N+2,N+NP+1)
        IF(WT.LT.RLU(0)*P(IP,5)*PV(1,5)**3/WTCOR(10)) GOTO 310

C...Angular distribution in W decay.
      ELSEIF(MMAT.EQ.46.AND.MSGN.NE.0) THEN
        IF(MSGN.GT.0) WT=FOUR(IM,N+1)*FOUR(N+2,IP+1)
        IF(MSGN.LT.0) WT=FOUR(IM,N+2)*FOUR(N+1,IP+1)
        IF(WT.LT.RLU(0)*P(IM,5)**4/WTCOR(10)) GOTO 370
      ENDIF

C...Scale back energy and reattach spectator.
      IF(MREM.EQ.1) THEN
        DO 450 J=1,5
  450   PV(1,J)=PV(1,J)/(1.-PQT)
        ND=ND+1
        MREM=0
      ENDIF

C...Low invariant mass for system with spectator quark gives particle,
C...not two jets. Readjust momenta accordingly.
      IF((MMAT.EQ.31.OR.MMAT.EQ.45).AND.ND.EQ.3) THEN
        MSTJ(93)=1
        PM2=ULMASS(K(N+2,2))
        MSTJ(93)=1
        PM3=ULMASS(K(N+3,2))
        IF(P(N+2,5)**2+P(N+3,5)**2+2.*FOUR(N+2,N+3).GE.
     &  (PARJ(32)+PM2+PM3)**2) GOTO 510
        K(N+2,1)=1
        KFTEMP=K(N+2,2)
        CALL LUKFDI(KFTEMP,K(N+3,2),KFLDMP,K(N+2,2))
        IF(K(N+2,2).EQ.0) GOTO 150
        P(N+2,5)=ULMASS(K(N+2,2))
        PS=P(N+1,5)+P(N+2,5)
        PV(2,5)=P(N+2,5)
        MMAT=0
        ND=2
        GOTO 370
      ELSEIF(MMAT.EQ.44) THEN
        MSTJ(93)=1
        PM3=ULMASS(K(N+3,2))
        MSTJ(93)=1
        PM4=ULMASS(K(N+4,2))
        IF(P(N+3,5)**2+P(N+4,5)**2+2.*FOUR(N+3,N+4).GE.
     &  (PARJ(32)+PM3+PM4)**2) GOTO 480
        K(N+3,1)=1
        KFTEMP=K(N+3,2)
        CALL LUKFDI(KFTEMP,K(N+4,2),KFLDMP,K(N+3,2))
        IF(K(N+3,2).EQ.0) GOTO 150
        P(N+3,5)=ULMASS(K(N+3,2))
        DO 460 J=1,3
  460   P(N+3,J)=P(N+3,J)+P(N+4,J)
        P(N+3,4)=SQRT(P(N+3,1)**2+P(N+3,2)**2+P(N+3,3)**2+P(N+3,5)**2)
        HA=P(N+1,4)**2-P(N+2,4)**2
        HB=HA-(P(N+1,5)**2-P(N+2,5)**2)
        HC=(P(N+1,1)-P(N+2,1))**2+(P(N+1,2)-P(N+2,2))**2+
     &  (P(N+1,3)-P(N+2,3))**2
        HD=(PV(1,4)-P(N+3,4))**2
        HE=HA**2-2.*HD*(P(N+1,4)**2+P(N+2,4)**2)+HD**2
        HF=HD*HC-HB**2
        HG=HD*HC-HA*HB
        HH=(SQRT(HG**2+HE*HF)-HG)/(2.*HF)
        DO 470 J=1,3
        PCOR=HH*(P(N+1,J)-P(N+2,J))
        P(N+1,J)=P(N+1,J)+PCOR
  470   P(N+2,J)=P(N+2,J)-PCOR
        P(N+1,4)=SQRT(P(N+1,1)**2+P(N+1,2)**2+P(N+1,3)**2+P(N+1,5)**2)
        P(N+2,4)=SQRT(P(N+2,1)**2+P(N+2,2)**2+P(N+2,3)**2+P(N+2,5)**2)
        ND=ND-1
      ENDIF

C...Check invariant mass of W jets. May give one particle or start over.
  480 IF(MMAT.GE.42.AND.MMAT.LE.44.AND.IABS(K(N+1,2)).LT.10) THEN
        PMR=SQRT(MAX(0.,P(N+1,5)**2+P(N+2,5)**2+2.*FOUR(N+1,N+2)))
        MSTJ(93)=1
        PM1=ULMASS(K(N+1,2))
        MSTJ(93)=1
        PM2=ULMASS(K(N+2,2))
        IF(PMR.GT.PARJ(32)+PM1+PM2) GOTO 490
        KFLDUM=INT(1.5+RLU(0))
        CALL LUKFDI(K(N+1,2),-ISIGN(KFLDUM,K(N+1,2)),KFLDMP,KF1)
        CALL LUKFDI(K(N+2,2),-ISIGN(KFLDUM,K(N+2,2)),KFLDMP,KF2)
        IF(KF1.EQ.0.OR.KF2.EQ.0) GOTO 150
        PSM=ULMASS(KF1)+ULMASS(KF2)
        IF(MMAT.EQ.42.AND.PMR.GT.PARJ(64)+PSM) GOTO 490
        IF(MMAT.GE.43.AND.PMR.GT.0.2*PARJ(32)+PSM) GOTO 490
        IF(ND.EQ.4.OR.KFA.EQ.15) GOTO 150
        K(N+1,1)=1
        KFTEMP=K(N+1,2)
        CALL LUKFDI(KFTEMP,K(N+2,2),KFLDMP,K(N+1,2))
        IF(K(N+1,2).EQ.0) GOTO 150
        P(N+1,5)=ULMASS(K(N+1,2))
        K(N+2,2)=K(N+3,2)
        P(N+2,5)=P(N+3,5)
        PS=P(N+1,5)+P(N+2,5)
        PV(2,5)=P(N+3,5)
        MMAT=0
        ND=2
        GOTO 370
      ENDIF

C...Phase space decay of partons from W decay.
  490 IF(MMAT.EQ.42.AND.IABS(K(N+1,2)).LT.10) THEN
        KFLO(1)=K(N+1,2)
        KFLO(2)=K(N+2,2)
        K(N+1,1)=K(N+3,1)
        K(N+1,2)=K(N+3,2)
        DO 500 J=1,5
        PV(1,J)=P(N+1,J)+P(N+2,J)
  500   P(N+1,J)=P(N+3,J)
        PV(1,5)=PMR
        N=N+1
        NP=0
        NQ=2
        PS=0.
        MSTJ(93)=2
        PSQ=ULMASS(KFLO(1))
        MSTJ(93)=2
        PSQ=PSQ+ULMASS(KFLO(2))
        MMAT=11
        GOTO 180
      ENDIF

C...Boost back for rapidly moving particle.
  510 N=N+ND
      IF(MBST.EQ.1) THEN
        DO 520 J=1,3
  520   BE(J)=P(IP,J)/P(IP,4)
        GA=P(IP,4)/P(IP,5)
        DO 540 I=NSAV+1,N
        BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
        DO 530 J=1,3
  530   P(I,J)=P(I,J)+GA*(GA*BEP/(1.+GA)+P(I,4))*BE(J)
  540   P(I,4)=GA*(P(I,4)+BEP)
      ENDIF

C...Fill in position of decay vertex.
      DO 560 I=NSAV+1,N
      DO 550 J=1,4
  550 V(I,J)=VDCY(J)
  560 V(I,5)=0.

C...Set up for parton shower evolution from jets.
      IF(MSTJ(23).GE.1.AND.MMAT.EQ.4.AND.K(NSAV+1,2).EQ.21) THEN
        K(NSAV+1,1)=3
        K(NSAV+2,1)=3
        K(NSAV+3,1)=3
        K(NSAV+1,4)=MSTU(5)*(NSAV+2)
        K(NSAV+1,5)=MSTU(5)*(NSAV+3)
        K(NSAV+2,4)=MSTU(5)*(NSAV+3)
        K(NSAV+2,5)=MSTU(5)*(NSAV+1)
        K(NSAV+3,4)=MSTU(5)*(NSAV+1)
        K(NSAV+3,5)=MSTU(5)*(NSAV+2)
        MSTJ(92)=-(NSAV+1)
      ELSEIF(MSTJ(23).GE.1.AND.MMAT.EQ.4) THEN
        K(NSAV+2,1)=3
        K(NSAV+3,1)=3
        K(NSAV+2,4)=MSTU(5)*(NSAV+3)
        K(NSAV+2,5)=MSTU(5)*(NSAV+3)
        K(NSAV+3,4)=MSTU(5)*(NSAV+2)
        K(NSAV+3,5)=MSTU(5)*(NSAV+2)
        MSTJ(92)=NSAV+2
      ELSEIF(MSTJ(23).GE.1.AND.(MMAT.EQ.32.OR.MMAT.EQ.44.OR.MMAT.EQ.46).
     &AND.IABS(K(NSAV+1,2)).LE.10.AND.IABS(K(NSAV+2,2)).LE.10) THEN
        K(NSAV+1,1)=3
        K(NSAV+2,1)=3
        K(NSAV+1,4)=MSTU(5)*(NSAV+2)
        K(NSAV+1,5)=MSTU(5)*(NSAV+2)
        K(NSAV+2,4)=MSTU(5)*(NSAV+1)
        K(NSAV+2,5)=MSTU(5)*(NSAV+1)
        MSTJ(92)=NSAV+1
      ELSEIF(MSTJ(23).GE.1.AND.MMAT.EQ.33.AND.IABS(K(NSAV+2,2)).EQ.21)
     &THEN
        K(NSAV+1,1)=3
        K(NSAV+2,1)=3
        K(NSAV+3,1)=3
        KCP=LUCOMP(K(NSAV+1,2))
        KQP=KCHG(KCP,2)*ISIGN(1,K(NSAV+1,2))
        JCON=4
        IF(KQP.LT.0) JCON=5
        K(NSAV+1,JCON)=MSTU(5)*(NSAV+2)
        K(NSAV+2,9-JCON)=MSTU(5)*(NSAV+1)
        K(NSAV+2,JCON)=MSTU(5)*(NSAV+3)
        K(NSAV+3,9-JCON)=MSTU(5)*(NSAV+2)
        MSTJ(92)=NSAV+1
      ELSEIF(MSTJ(23).GE.1.AND.MMAT.EQ.33) THEN
        K(NSAV+1,1)=3
        K(NSAV+3,1)=3
        K(NSAV+1,4)=MSTU(5)*(NSAV+3)
        K(NSAV+1,5)=MSTU(5)*(NSAV+3)
        K(NSAV+3,4)=MSTU(5)*(NSAV+1)
        K(NSAV+3,5)=MSTU(5)*(NSAV+1)
        MSTJ(92)=NSAV+1
      ENDIF

C...Mark decayed particle.
      IF(K(IP,1).EQ.5) K(IP,1)=15
      IF(K(IP,1).LE.10) K(IP,1)=11
      K(IP,4)=NSAV+1
      K(IP,5)=N

      RETURN
      END
