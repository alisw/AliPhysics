 
C*********************************************************************
 
      SUBROUTINE PYRESD
 
C...Allows resonances to decay (including parton showers for hadronic
C...channels).
      IMPLICIT DOUBLE PRECISION(D)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT4/
      DIMENSION IREF(20,8),KDCY(3),KFL1(3),KFL2(3),KEQL(3),NSD(3),
     &ILIN(6),HGZ(3,3),COUP(6,4),CORL(2,2,2),PK(6,4),PKK(6,6),
     &CTHE(3),PHI(3),WDTP(0:40),WDTE(0:40,0:5),DBEZQQ(3),DPMO(5)
      COMPLEX FGK,HA(6,6),HC(6,6)
 
C...The F, Xi and Xj functions of Gunion and Kunszt
C...(Phys. Rev. D33, 665, plus errata from the authors).
      FGK(I1,I2,I3,I4,I5,I6)=4.*HA(I1,I3)*HC(I2,I6)*(HA(I1,I5)*
     &HC(I1,I4)+HA(I3,I5)*HC(I3,I4))
      DIGK(DT,DU)=-4.*D34*D56+DT*(3.*DT+4.*DU)+DT**2*(DT*DU/(D34*D56)-
     &2.*(1./D34+1./D56)*(DT+DU)+2.*(D34/D56+D56/D34))
      DJGK(DT,DU)=8.*(D34+D56)**2-8.*(D34+D56)*(DT+DU)-6.*DT*DU-
     &2.*DT*DU*(DT*DU/(D34*D56)-2.*(1./D34+1./D56)*(DT+DU)+
     &2.*(D34/D56+D56/D34))
 
C...Some general constants.
      XW=PARU(102)
      XWV=XW
      IF(MSTP(8).GE.2) XW=1.-(PMAS(24,1)/PMAS(23,1))**2
      XW1=1.-XW
      SQMZ=PMAS(23,1)**2
      SQMW=PMAS(24,1)**2
      SH=VINT(44)
 
C...Define initial one, two or three objects.
      ISUB=MINT(1)
      DO 100 JT=1,8
      IREF(1,JT)=0
  100 CONTINUE
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN
        IREF(1,1)=MINT(84)+2+ISET(ISUB)
        IREF(1,4)=MINT(83)+6+ISET(ISUB)
      ELSEIF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN
        IREF(1,1)=MINT(84)+1+ISET(ISUB)
        IREF(1,2)=MINT(84)+2+ISET(ISUB)
        IREF(1,4)=MINT(83)+5+ISET(ISUB)
        IREF(1,5)=MINT(83)+6+ISET(ISUB)
      ELSEIF(ISET(ISUB).EQ.5) THEN
        IREF(1,1)=MINT(84)+3
        IREF(1,2)=MINT(84)+4
        IREF(1,3)=MINT(84)+5
        IREF(1,4)=MINT(83)+7
        IREF(1,5)=MINT(83)+8
        IREF(1,6)=MINT(83)+9
      ELSEIF(ISET(ISUB).EQ.6) THEN
        IREF(1,1)=MINT(84)+4
        IREF(1,2)=MINT(84)+5
        IREF(1,3)=MINT(84)+3
        IREF(1,4)=MINT(83)+8
        IREF(1,5)=MINT(83)+9
        IREF(1,6)=MINT(83)+7
      ENDIF
 
C...Check if initial resonance has been moved (in resonance + jet).
      DO 120 JT=1,3
      IF(IREF(1,JT).GT.0) THEN
        IF(K(IREF(1,JT),1).GT.10) THEN
          KFA=IABS(K(IREF(1,JT),2))
          IF(KFA.EQ.6.OR.KFA.EQ.7.OR.KFA.EQ.8.OR.KFA.EQ.39) THEN
            DO 110 I=IREF(1,JT)+1,N
            IF(K(I,1).LE.10.AND.K(I,2).EQ.K(IREF(1,JT),2)) IREF(1,JT)=I
  110       CONTINUE
          ELSE
            KDA=MOD(K(IREF(1,JT),4),MSTU(4))
            IF(KFA.GE.23.AND.KFA.LE.40.AND.KDA.GT.1) IREF(1,JT)=KDA
          ENDIF
        ENDIF
      ENDIF
  120 CONTINUE
 
C...Loop over decay history.
      NP=1
      IP=0
  130 IP=IP+1
      NINH=0
      JTMAX=2
      IF(IP.EQ.1.AND.IREF(1,2).EQ.0) JTMAX=1
      IF(IP.EQ.1.AND.IREF(1,3).NE.0) JTMAX=3
      ITLH=0
      NSAV=N
 
C...Start treatment of one or two resonances in parallel.
  140 N=NSAV
      DO 170 JT=1,JTMAX
      ID=IREF(IP,JT)
      KDCY(JT)=0
      KFL1(JT)=0
      KFL2(JT)=0
      KEQL(JT)=0
      NSD(JT)=ID
      IF(ID.EQ.0) GOTO 160
      KFA=IABS(K(ID,2))
      IF((KFA.LT.23.OR.KFA.GT.40).AND.KFA.NE.6.AND.KFA.NE.7.AND.
     &KFA.NE.8.AND.KFA.NE.17.AND.KFA.NE.18) GOTO 160
      IF(MSTP(48).LE.0.AND.KFA.EQ.6) GOTO 160
      IF(MSTP(6).NE.1.AND.MSTP(49).LE.0.AND.(KFA.EQ.7.OR.KFA.EQ.8.OR.
     &KFA.EQ.17.OR.KFA.EQ.18)) GOTO 160
      IF(K(ID,1).GT.10.OR.MDCY(KFA,1).EQ.0) GOTO 160
      IF(KFA.EQ.6.OR.(MSTP(6).NE.1.AND.(KFA.EQ.7.OR.KFA.EQ.8.OR.
     &KFA.EQ.17.OR.KFA.EQ.18))) ITLH=ITLH+1
      K(ID,4)=MSTU(5)*(K(ID,4)/MSTU(5))
      K(ID,5)=MSTU(5)*(K(ID,5)/MSTU(5))
 
C...Select decay channel.
      KFB=0
      IF(ISET(ISUB).NE.6.OR.JT.NE.3) THEN
        IF(ISUB.EQ.1.OR.ISUB.EQ.15.OR.ISUB.EQ.19.OR.ISUB.EQ.22.OR.
     &  ISUB.EQ.30.OR.ISUB.EQ.35.OR.ISUB.EQ.141) MINT(61)=1
        CALL PYWIDT(KFA,P(ID,5)**2,WDTP,WDTE)
        IF(KCHG(KFA,3).EQ.0) THEN
          IPM=2
        ELSE
          IPM=(5-ISIGN(1,K(ID,2)))/2
        ENDIF
        IF(JTMAX.GE.2.AND.JT.LE.2) KFB=IABS(K(IREF(IP,3-JT),2))
        WDTE0S=WDTE(0,1)+WDTE(0,IPM)+WDTE(0,4)
        IF(KFB.EQ.KFA) WDTE0S=WDTE0S+WDTE(0,5)
        IF(WDTE0S.LE.0.) THEN
          IF(KFA.EQ.6.OR.KFA.EQ.7.OR.KFA.EQ.8.OR.KFA.EQ.17.OR.
     &    KFA.EQ.18) THEN
            MINT(51)=1
            RETURN
          ELSE
            GOTO 160
          ENDIF
        ENDIF
        RKFL=WDTE0S*RLU(0)
        IDL=0
  150   IDL=IDL+1
        IDC=IDL+MDCY(KFA,2)-1
        RKFL=RKFL-(WDTE(IDL,1)+WDTE(IDL,IPM)+WDTE(IDL,4))
        IF(KFB.EQ.KFA) RKFL=RKFL-WDTE(IDL,5)
        IF(IDL.LT.MDCY(KFA,3).AND.RKFL.GT.0.) GOTO 150
      ELSE
        IDC=MINT(35)
      ENDIF
 
C...Read out and classify decay channel chosen.
      KFL1(JT)=KFDP(IDC,1)*ISIGN(1,K(ID,2))
      KFC1A=IABS(KFL1(JT))
      IF(KFC1A.GT.100) KFC1A=LUCOMP(KFC1A)
      IF(KCHG(KFC1A,3).EQ.0) KFL1(JT)=IABS(KFL1(JT))
      KFL2(JT)=KFDP(IDC,2)*ISIGN(1,K(ID,2))
      KFC2A=IABS(KFL2(JT))
      IF(KFC2A.GT.100) KFC2A=LUCOMP(KFC2A)
      IF(KCHG(KFC2A,3).EQ.0) KFL2(JT)=IABS(KFL2(JT))
      KDCY(JT)=2
      IF(IABS(KFL1(JT)).LE.10.OR.KFL1(JT).EQ.21) KDCY(JT)=1
      IF(IABS(KFL1(JT)).GE.23.AND.IABS(KFL1(JT)).LE.40) KDCY(JT)=3
      IF(KFB.EQ.KFA) KEQL(JT)=MDME(IDC,1)
      NSD(JT)=N
      HGZ(JT,1)=VINT(111)
      HGZ(JT,2)=VINT(112)
      HGZ(JT,3)=VINT(114)
 
C...Select masses and check that mass sum not too large.
      IF(MSTP(42).LE.0.OR.(PMAS(KFC1A,2).LT.PARP(41).AND.
     &PMAS(KFC2A,2).LT.PARP(41))) THEN
        P(N+1,5)=PMAS(KFC1A,1)
        P(N+2,5)=PMAS(KFC2A,1)
        IF(P(N+1,5)+P(N+2,5)+PARJ(64).GT.P(ID,5)) THEN
          CALL LUERRM(13,'(PYRESD:) daughter masses too large')
          MINT(51)=1
          RETURN
        ENDIF
      ELSEIF(IP.EQ.1) THEN
        CALL PYOFSH(2,KFA,KFL1(JT),KFL2(JT),P(ID,5),P(N+1,5),P(N+2,5))
        IF(MINT(51).EQ.1) RETURN
      ELSE
        CALL PYOFSH(7,KFA,KFL1(JT),KFL2(JT),P(ID,5),P(N+1,5),P(N+2,5))
        IF(MINT(51).EQ.1) RETURN
      ENDIF
 
C...Fill decay products, prepared for parton showers for quarks.
C...Special cases, done by hand, for techni-eta, t, leptoquark and q*.
      MSTU(10)=1
      IF(KFA.EQ.38.OR.KFA.EQ.39.OR.((MSTP(6).EQ.1.OR.MSTP(49).GE.1).AND.
     &(KFA.EQ.7.OR.KFA.EQ.8)).OR.KFA.EQ.6) THEN
        MSTU(19)=1
        CALL LU2ENT(N+1,KFL1(JT),KFL2(JT),P(ID,5))
        ISID=4
        IF(K(ID,2).LT.0) ISID=5
        IF(KFA.EQ.38) THEN
          IF(KFC1A.EQ.21.AND.RLU(0).GT.0.5) ISID=9-ISID
          K(N-1,1)=3
          K(N,1)=3
          K(ID,ISID)=K(ID,ISID)+(N-1)
          K(ID,9-ISID)=K(ID,9-ISID)+N
          K(N-1,ISID)=MSTU(5)*ID
          K(N-1,9-ISID)=MSTU(5)*N
          K(N,ISID)=MSTU(5)*(N-1)
          K(N,9-ISID)=MSTU(5)*ID
        ELSEIF(KFA.EQ.6.OR.(MSTP(6).NE.1.AND.(KFA.EQ.7.OR.KFA.EQ.8)))
     &  THEN
          K(N-1,1)=1
          K(N,1)=3
          K(ID,ISID)=K(ID,ISID)+N
          K(N,ISID)=MSTU(5)*ID
        ELSEIF(KFA.EQ.39) THEN
          K(N-1,1)=3
          K(N,1)=1
          K(ID,ISID)=K(ID,ISID)+(N-1)
          K(N-1,ISID)=MSTU(5)*ID
        ELSEIF(KFL1(JT).NE.21) THEN
          K(N-1,1)=1
          K(N,1)=3
          K(ID,ISID)=K(ID,ISID)+N
          K(N,ISID)=MSTU(5)*ID
        ELSE
          K(N-1,1)=3
          K(N,1)=3
          K(ID,ISID)=K(ID,ISID)+(N-1)
          K(N-1,ISID)=MSTU(5)*ID
          K(N-1,9-ISID)=MSTU(5)*N
          K(N,ISID)=MSTU(5)*(N-1)
        ENDIF
      ELSEIF(KDCY(JT).EQ.1) THEN
        CALL LU2ENT(-(N+1),KFL1(JT),KFL2(JT),P(ID,5))
      ELSE
        CALL LU2ENT(N+1,KFL1(JT),KFL2(JT),P(ID,5))
      ENDIF
      MSTU(10)=2
  160 IF(KFA.GE.23.AND.KFA.LE.40.AND.KFL1(JT).EQ.0) NINH=NINH+1
  170 CONTINUE
 
C...Check for allowed combinations. Skip if no decays.
      IF(JTMAX.GE.2) THEN
        IF(KEQL(1).EQ.4.AND.KEQL(2).EQ.4) GOTO 140
        IF(KEQL(1).EQ.5.AND.KEQL(2).EQ.5) GOTO 140
      ENDIF
      IF(JTMAX.EQ.1.AND.KDCY(1).EQ.0) GOTO 480
      IF(JTMAX.EQ.2.AND.KDCY(1).EQ.0.AND.KDCY(2).EQ.0) GOTO 480
      IF(JTMAX.EQ.3.AND.KDCY(1).EQ.0.AND.KDCY(2).EQ.0.AND.
     &KDCY(3).EQ.0) GOTO 480
 
C...Order incoming partons and outgoing resonances.
      IF(JTMAX.EQ.2.AND.MSTP(47).GE.1.AND.NINH.EQ.0) THEN
        ILIN(1)=MINT(84)+1
        IF(K(MINT(84)+1,2).GT.0) ILIN(1)=MINT(84)+2
        IF(K(ILIN(1),2).EQ.21) ILIN(1)=2*MINT(84)+3-ILIN(1)
        ILIN(2)=2*MINT(84)+3-ILIN(1)
        IMIN=1
        IF(IREF(IP,7).EQ.25.OR.IREF(IP,7).EQ.35.OR.IREF(IP,7)
     &  .EQ.36) IMIN=3
        IMAX=2
        IORD=1
        IF(K(IREF(IP,1),2).EQ.23) IORD=2
        IF(K(IREF(IP,1),2).EQ.24.AND.K(IREF(IP,2),2).EQ.-24) IORD=2
        IAKIPD=IABS(K(IREF(IP,IORD),2))
        IF(IAKIPD.EQ.25.OR.IAKIPD.EQ.35.OR.IAKIPD.EQ.36) IORD=3-IORD
        IF(KDCY(IORD).EQ.0) IORD=3-IORD
 
C...Order decay products of resonances.
        DO 180 JT=IORD,3-IORD,3-2*IORD
        IF(KDCY(JT).EQ.0) THEN
          ILIN(IMAX+1)=NSD(JT)
          IMAX=IMAX+1
        ELSEIF(K(NSD(JT)+1,2).GT.0) THEN
          ILIN(IMAX+1)=N+2*JT-1
          ILIN(IMAX+2)=N+2*JT
          IMAX=IMAX+2
          K(N+2*JT-1,2)=K(NSD(JT)+1,2)
          K(N+2*JT,2)=K(NSD(JT)+2,2)
        ELSE
          ILIN(IMAX+1)=N+2*JT
          ILIN(IMAX+2)=N+2*JT-1
          IMAX=IMAX+2
          K(N+2*JT-1,2)=K(NSD(JT)+1,2)
          K(N+2*JT,2)=K(NSD(JT)+2,2)
        ENDIF
  180   CONTINUE
 
C...Find charge, isospin, left- and righthanded couplings.
        DO 200 I=IMIN,IMAX
        DO 190 J=1,4
        COUP(I,J)=0.
  190   CONTINUE
        KFA=IABS(K(ILIN(I),2))
        IF(KFA.EQ.0.OR.KFA.GT.20) GOTO 200
        COUP(I,1)=KCHG(KFA,1)/3.
        COUP(I,2)=(-1)**MOD(KFA,2)
        COUP(I,4)=-2.*COUP(I,1)*XWV
        COUP(I,3)=COUP(I,2)+COUP(I,4)
  200   CONTINUE
 
C...Full propagator dependence and flavour correlations for 2 gamma*/Z.
        IF(ISUB.EQ.22) THEN
          DO 230 I=3,5,2
          I1=IORD
          IF(I.EQ.5) I1=3-IORD
          DO 220 J1=1,2
          DO 210 J2=1,2
          CORL(I/2,J1,J2)=COUP(1,1)**2*HGZ(I1,1)*COUP(I,1)**2/16.+
     &    COUP(1,1)*COUP(1,J1+2)*HGZ(I1,2)*COUP(I,1)*COUP(I,J2+2)/4.+
     &    COUP(1,J1+2)**2*HGZ(I1,3)*COUP(I,J2+2)**2
  210     CONTINUE
  220     CONTINUE
  230     CONTINUE
          COWT12=(CORL(1,1,1)+CORL(1,1,2))*(CORL(2,1,1)+CORL(2,1,2))+
     &    (CORL(1,2,1)+CORL(1,2,2))*(CORL(2,2,1)+CORL(2,2,2))
          COMX12=(CORL(1,1,1)+CORL(1,1,2)+CORL(1,2,1)+CORL(1,2,2))*
     &    (CORL(2,1,1)+CORL(2,1,2)+CORL(2,2,1)+CORL(2,2,2))
          IF(COWT12.LT.RLU(0)*COMX12) GOTO 140
        ENDIF
      ENDIF
 
C...Select angular orientation type - Z'/W' only.
      MZPWP=0
      IF(ISUB.EQ.141) THEN
        IF(RLU(0).LT.PARU(130)) MZPWP=1
        IF(IP.EQ.2) THEN
          IF(IABS(K(IREF(2,1),2)).EQ.37) MZPWP=2
          IAKIR=IABS(K(IREF(2,2),2))
          IF(IAKIR.EQ.25.OR.IAKIR.EQ.35.OR.IAKIR.EQ.36) MZPWP=2
        ENDIF
        IF(IP.GE.3) MZPWP=2
      ELSEIF(ISUB.EQ.142) THEN
        IF(RLU(0).LT.PARU(136)) MZPWP=1
        IF(IP.EQ.2) THEN
          IAKIR=IABS(K(IREF(2,2),2))
          IF(IAKIR.EQ.25.OR.IAKIR.EQ.35.OR.IAKIR.EQ.36) MZPWP=2
        ENDIF
        IF(IP.GE.3) MZPWP=2
      ENDIF
 
C...Select random angles (begin of weighting procedure).
  240 DO 250 JT=1,JTMAX
      IF(KDCY(JT).EQ.0) GOTO 250
      IF(ISET(ISUB).EQ.6.AND.JT.EQ.3) GOTO 250
      IF(JTMAX.EQ.1) THEN
        CTHE(JT)=VINT(13)+(VINT(33)-VINT(13)+VINT(34)-VINT(14))*RLU(0)
        IF(CTHE(JT).GT.VINT(33)) CTHE(JT)=CTHE(JT)+VINT(14)-VINT(33)
        PHI(JT)=VINT(24)
      ELSE
        CTHE(JT)=2.*RLU(0)-1.
        PHI(JT)=PARU(2)*RLU(0)
      ENDIF
  250 CONTINUE
 
      IF(JTMAX.EQ.2.AND.MSTP(47).GE.1.AND.NINH.EQ.0) THEN
C...Construct massless four-vectors.
        DO 270 I=N+1,N+4
        K(I,1)=1
        DO 260 J=1,5
        P(I,J)=0.
        V(I,J)=0.
  260   CONTINUE
  270   CONTINUE
        DO 280 JT=1,JTMAX
        IF(KDCY(JT).EQ.0) GOTO 280
        ID=IREF(IP,JT)
        P(N+2*JT-1,3)=0.5*P(ID,5)
        P(N+2*JT-1,4)=0.5*P(ID,5)
        P(N+2*JT,3)=-0.5*P(ID,5)
        P(N+2*JT,4)=0.5*P(ID,5)
        CALL LUDBRB(N+2*JT-1,N+2*JT,ACOS(CTHE(JT)),PHI(JT),DBLE(P(ID,1)/
     &  P(ID,4)),DBLE(P(ID,2)/P(ID,4)),DBLE(P(ID,3)/P(ID,4)))
  280   CONTINUE
 
C...Store incoming and outgoing momenta, with random rotation to
C...avoid accidental zeroes in HA expressions.
        DO 300 I=1,IMAX
        K(N+4+I,1)=1
        P(N+4+I,4)=SQRT(P(ILIN(I),1)**2+P(ILIN(I),2)**2+P(ILIN(I),3)**2+
     &  P(ILIN(I),5)**2)
        P(N+4+I,5)=P(ILIN(I),5)
        DO 290 J=1,3
        P(N+4+I,J)=P(ILIN(I),J)
  290   CONTINUE
  300   CONTINUE
  310   THERR=ACOS(2.*RLU(0)-1.)
        PHIRR=PARU(2)*RLU(0)
        CALL LUDBRB(N+5,N+4+IMAX,THERR,PHIRR,0D0,0D0,0D0)
        DO 330 I=1,IMAX
        IF(P(N+4+I,1)**2+P(N+4+I,2)**2.LT.1E-4*P(N+4+I,4)**2) GOTO 310
        DO 320 J=1,4
        PK(I,J)=P(N+4+I,J)
  320   CONTINUE
  330   CONTINUE
 
C...Calculate internal products.
        IF(ISUB.EQ.22.OR.ISUB.EQ.23.OR.ISUB.EQ.25.OR.ISUB.EQ.141.OR.
     &  ISUB.EQ.142) THEN
          DO 350 I1=IMIN,IMAX-1
          DO 340 I2=I1+1,IMAX
          HA(I1,I2)=SQRT((PK(I1,4)-PK(I1,3))*(PK(I2,4)+PK(I2,3))/
     &    (1E-20+PK(I1,1)**2+PK(I1,2)**2))*CMPLX(PK(I1,1),PK(I1,2))-
     &    SQRT((PK(I1,4)+PK(I1,3))*(PK(I2,4)-PK(I2,3))/
     &    (1E-20+PK(I2,1)**2+PK(I2,2)**2))*CMPLX(PK(I2,1),PK(I2,2))
          HC(I1,I2)=CONJG(HA(I1,I2))
          IF(I1.LE.2) HA(I1,I2)=CMPLX(0.,1.)*HA(I1,I2)
          IF(I1.LE.2) HC(I1,I2)=CMPLX(0.,1.)*HC(I1,I2)
          HA(I2,I1)=-HA(I1,I2)
          HC(I2,I1)=-HC(I1,I2)
  340     CONTINUE
  350     CONTINUE
        ENDIF
        DO 370 I=1,2
        DO 360 J=1,4
        PK(I,J)=-PK(I,J)
  360   CONTINUE
  370   CONTINUE
        DO 390 I1=IMIN,IMAX-1
        DO 380 I2=I1+1,IMAX
        PKK(I1,I2)=2.*(PK(I1,4)*PK(I2,4)-PK(I1,1)*PK(I2,1)-
     &  PK(I1,2)*PK(I2,2)-PK(I1,3)*PK(I2,3))
        PKK(I2,I1)=PKK(I1,I2)
  380   CONTINUE
  390   CONTINUE
      ENDIF
 
      KFAGM=IABS(IREF(IP,7))
      IF(MSTP(47).LE.0.OR.NINH.NE.0) THEN
C...Isotropic decay selected by user.
        WT=1.
        WTMAX=1.
 
      ELSEIF(ITLH.GE.1) THEN
C... Isotropic decay t -> b + W etc for 4th generation q and l.
        WT=1.
        WTMAX=1.
 
      ELSEIF(IREF(IP,7).EQ.25.OR.IREF(IP,7).EQ.35.OR.
     &IREF(IP,7).EQ.36) THEN
C...Angular weight for H0 -> Z0 + Z0 or W+ + W- -> 4 quarks/leptons.
        WT=16.*PKK(3,5)*PKK(4,6)
        IF(IP.EQ.1) WTMAX=SH**2
        IF(IP.GE.2) WTMAX=P(IREF(IP,8),5)**4
        KFA=IABS(K(IREF(IP,1),2))
        IF(KFA.NE.23.AND.KFA.NE.24) WT=WTMAX
 
      ELSEIF((KFAGM.EQ.6.OR.(MSTP(6).NE.1.AND.(KFAGM.EQ.7.OR.
     &KFAGM.EQ.8.OR.KFAGM.EQ.17.OR.KFAGM.EQ.18))).AND.
     &IABS(K(IREF(IP,1),2)).EQ.24) THEN
C...Angular correlation in f -> f' + W -> f' + 2 quarks/leptons.
        I1=IREF(IP,8)
        IF(MOD(KFAGM,2).EQ.0) THEN
          I2=N+1
          I3=N+2
        ELSE
          I2=N+2
          I3=N+1
        ENDIF
        I4=IREF(IP,2)
        WT=(P(I1,4)*P(I2,4)-P(I1,1)*P(I2,1)-P(I1,2)*P(I2,2)-
     &  P(I1,3)*P(I2,3))*(P(I3,4)*P(I4,4)-P(I3,1)*P(I4,1)-
     &  P(I3,2)*P(I4,2)-P(I3,3)*P(I4,3))
        WTMAX=(P(I1,5)**4-P(IREF(IP,1),5)**4)/8.
        IF(KFAGM.EQ.6.AND.MSTP(48).LE.1) WT=WTMAX
        IF(KFAGM.NE.6.AND.MSTP(49).LE.1) WT=WTMAX
 
      ELSEIF(ISUB.EQ.1) THEN
C...Angular weight for gamma*/Z0 -> 2 quarks/leptons.
        EI=KCHG(IABS(MINT(15)),1)/3.
        AI=SIGN(1.,EI+0.1)
        VI=AI-4.*EI*XWV
        EF=KCHG(IABS(KFL1(1)),1)/3.
        AF=SIGN(1.,EF+0.1)
        VF=AF-4.*EF*XWV
        ASYM=2.*(EI*AI*VINT(112)*EF*AF+4.*VI*AI*VINT(114)*VF*AF)/
     &  (EI**2*VINT(111)*EF**2+EI*VI*VINT(112)*EF*VF+
     &  (VI**2+AI**2)*VINT(114)*(VF**2+AF**2))
        WT=1.+ASYM*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1))+CTHE(1)**2
        WTMAX=2.+ABS(ASYM)
 
      ELSEIF(ISUB.EQ.2) THEN
C...Angular weight for W+/- -> 2 quarks/leptons.
        WT=(1.+CTHE(1)*ISIGN(1,MINT(15)*KFL1(1)))**2
        WTMAX=4.
 
      ELSEIF(ISUB.EQ.15.OR.ISUB.EQ.19) THEN
C...Angular weight for f + f~ -> gluon/gamma + (gamma*/Z0) ->
C...-> gluon/gamma + 2 quarks/leptons.
        CLILF=COUP(1,1)**2*HGZ(2,1)*COUP(3,1)**2/16.+
     &  COUP(1,1)*COUP(1,3)*HGZ(2,2)*COUP(3,1)*COUP(3,3)/4.+
     &  COUP(1,3)**2*HGZ(2,3)*COUP(3,3)**2
        CLIRF=COUP(1,1)**2*HGZ(2,1)*COUP(3,1)**2/16.+
     &  COUP(1,1)*COUP(1,3)*HGZ(2,2)*COUP(3,1)*COUP(3,4)/4.+
     &  COUP(1,3)**2*HGZ(2,3)*COUP(3,4)**2
        CRILF=COUP(1,1)**2*HGZ(2,1)*COUP(3,1)**2/16.+
     &  COUP(1,1)*COUP(1,4)*HGZ(2,2)*COUP(3,1)*COUP(3,3)/4.+
     &  COUP(1,4)**2*HGZ(2,3)*COUP(3,3)**2
        CRIRF=COUP(1,1)**2*HGZ(2,1)*COUP(3,1)**2/16.+
     &  COUP(1,1)*COUP(1,4)*HGZ(2,2)*COUP(3,1)*COUP(3,4)/4.+
     &  COUP(1,4)**2*HGZ(2,3)*COUP(3,4)**2
        WT=(CLILF+CRIRF)*(PKK(1,3)**2+PKK(2,4)**2)+
     &  (CLIRF+CRILF)*(PKK(1,4)**2+PKK(2,3)**2)
        WTMAX=(CLILF+CLIRF+CRILF+CRIRF)*
     &  ((PKK(1,3)+PKK(1,4))**2+(PKK(2,3)+PKK(2,4))**2)
 
      ELSEIF(ISUB.EQ.16.OR.ISUB.EQ.20) THEN
C...Angular weight for f + f~' -> gluon/gamma + W+/- ->
C...-> gluon/gamma + 2 quarks/leptons.
        WT=PKK(1,3)**2+PKK(2,4)**2
        WTMAX=(PKK(1,3)+PKK(1,4))**2+(PKK(2,3)+PKK(2,4))**2
 
      ELSEIF(ISUB.EQ.22) THEN
C...Angular weight for f + f~ -> Z0 + Z0 -> 4 quarks/leptons.
        S34=P(IREF(IP,IORD),5)**2
        S56=P(IREF(IP,3-IORD),5)**2
        TI=PKK(1,3)+PKK(1,4)+S34
        UI=PKK(1,5)+PKK(1,6)+S56
        FGK135=ABS(FGK(1,2,3,4,5,6)/TI+FGK(1,2,5,6,3,4)/UI)**2
        FGK145=ABS(FGK(1,2,4,3,5,6)/TI+FGK(1,2,5,6,4,3)/UI)**2
        FGK136=ABS(FGK(1,2,3,4,6,5)/TI+FGK(1,2,6,5,3,4)/UI)**2
        FGK146=ABS(FGK(1,2,4,3,6,5)/TI+FGK(1,2,6,5,4,3)/UI)**2
        FGK253=ABS(FGK(2,1,5,6,3,4)/TI+FGK(2,1,3,4,5,6)/UI)**2
        FGK263=ABS(FGK(2,1,6,5,3,4)/TI+FGK(2,1,3,4,6,5)/UI)**2
        FGK254=ABS(FGK(2,1,5,6,4,3)/TI+FGK(2,1,4,3,5,6)/UI)**2
        FGK264=ABS(FGK(2,1,6,5,4,3)/TI+FGK(2,1,4,3,6,5)/UI)**2
        WT=
     &  CORL(1,1,1)*CORL(2,1,1)*FGK135+CORL(1,1,2)*CORL(2,1,1)*FGK145+
     &  CORL(1,1,1)*CORL(2,1,2)*FGK136+CORL(1,1,2)*CORL(2,1,2)*FGK146+
     &  CORL(1,2,1)*CORL(2,2,1)*FGK253+CORL(1,2,2)*CORL(2,2,1)*FGK263+
     &  CORL(1,2,1)*CORL(2,2,2)*FGK254+CORL(1,2,2)*CORL(2,2,2)*FGK264
        WTMAX=16.*((CORL(1,1,1)+CORL(1,1,2))*(CORL(2,1,1)+CORL(2,1,2))+
     &  (CORL(1,2,1)+CORL(1,2,2))*(CORL(2,2,1)+CORL(2,2,2)))*S34*S56*
     &  ((TI**2+UI**2+2.*SH*(S34+S56))/(TI*UI)-S34*S56*(1./TI**2+
     &  1./UI**2))
 
      ELSEIF(ISUB.EQ.23) THEN
C...Angular weight for f + f~' -> Z0 + W+/- -> 4 quarks/leptons.
        D34=P(IREF(IP,IORD),5)**2
        D56=P(IREF(IP,3-IORD),5)**2
        DT=PKK(1,3)+PKK(1,4)+D34
        DU=PKK(1,5)+PKK(1,6)+D56
        FACBW=1./((SH-SQMW)**2+SQMW*PMAS(24,2)**2)
        CAWZ=COUP(2,3)/SNGL(DT)-2.*XW1*COUP(1,2)*(SH-SQMW)*FACBW
        CBWZ=COUP(1,3)/SNGL(DU)+2.*XW1*COUP(1,2)*(SH-SQMW)*FACBW
        FGK135=ABS(CAWZ*FGK(1,2,3,4,5,6)+CBWZ*FGK(1,2,5,6,3,4))
        FGK136=ABS(CAWZ*FGK(1,2,3,4,6,5)+CBWZ*FGK(1,2,6,5,3,4))
        WT=(COUP(5,3)*FGK135)**2+(COUP(5,4)*FGK136)**2
        WTMAX=4.*D34*D56*(COUP(5,3)**2+COUP(5,4)**2)*(CAWZ**2*
     &  DIGK(DT,DU)+CBWZ**2*DIGK(DU,DT)+CAWZ*CBWZ*DJGK(DT,DU))
 
      ELSEIF(ISUB.EQ.24.OR.ISUB.EQ.171.OR.ISUB.EQ.176) THEN
C...Angular weight for f + f~ -> Z0 + H0 -> 2 quarks/leptons + H0
C...(or H'0, or A0).
        WT=((COUP(1,3)*COUP(3,3))**2+(COUP(1,4)*COUP(3,4))**2)*
     &  PKK(1,3)*PKK(2,4)+((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)*
     &  COUP(3,3))**2)*PKK(1,4)*PKK(2,3)
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*
     &  (PKK(1,3)+PKK(1,4))*(PKK(2,3)+PKK(2,4))
 
      ELSEIF(ISUB.EQ.25) THEN
C...Angular weight for f + f~ -> W+ + W- -> 4 quarks/leptons.
        D34=P(IREF(IP,IORD),5)**2
        D56=P(IREF(IP,3-IORD),5)**2
        DT=PKK(1,3)+PKK(1,4)+D34
        DU=PKK(1,5)+PKK(1,6)+D56
        FACBW=1./((SH-SQMZ)**2+SQMZ*PMAS(23,2)**2)
        CDWW=(COUP(1,3)*SQMZ*(SH-SQMZ)*FACBW+COUP(1,2))/SH
        CAWW=CDWW+0.5*(COUP(1,2)+1.)/SNGL(DT)
        CBWW=CDWW+0.5*(COUP(1,2)-1.)/SNGL(DU)
        CCWW=COUP(1,4)*SQMZ*(SH-SQMZ)*FACBW/SH
        FGK135=ABS(CAWW*FGK(1,2,3,4,5,6)-CBWW*FGK(1,2,5,6,3,4))
        FGK253=ABS(FGK(2,1,5,6,3,4)-FGK(2,1,3,4,5,6))
        WT=FGK135**2+(CCWW*FGK253)**2
        WTMAX=4.*D34*D56*(CAWW**2*DIGK(DT,DU)+CBWW**2*DIGK(DU,DT)-CAWW*
     &  CBWW*DJGK(DT,DU)+CCWW**2*(DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU)))
 
      ELSEIF(ISUB.EQ.26.OR.ISUB.EQ.172.OR.ISUB.EQ.177) THEN
C...Angular weight for f + f~' -> W+/- + H0 -> 2 quarks/leptons + H0
C...(or H'0, or A0).
        WT=PKK(1,3)*PKK(2,4)
        WTMAX=(PKK(1,3)+PKK(1,4))*(PKK(2,3)+PKK(2,4))
 
      ELSEIF(ISUB.EQ.30.OR.ISUB.EQ.35) THEN
C...Angular weight for f + g/gamma -> f + (gamma*/Z0)
C...-> f + 2 quarks/leptons.
        CLILF=COUP(1,1)**2*HGZ(2,1)*COUP(3,1)**2/16.+
     &  COUP(1,1)*COUP(1,3)*HGZ(2,2)*COUP(3,1)*COUP(3,3)/4.+
     &  COUP(1,3)**2*HGZ(2,3)*COUP(3,3)**2
        CLIRF=COUP(1,1)**2*HGZ(2,1)*COUP(3,1)**2/16.+
     &  COUP(1,1)*COUP(1,3)*HGZ(2,2)*COUP(3,1)*COUP(3,4)/4.+
     &  COUP(1,3)**2*HGZ(2,3)*COUP(3,4)**2
        CRILF=COUP(1,1)**2*HGZ(2,1)*COUP(3,1)**2/16.+
     &  COUP(1,1)*COUP(1,4)*HGZ(2,2)*COUP(3,1)*COUP(3,3)/4.+
     &  COUP(1,4)**2*HGZ(2,3)*COUP(3,3)**2
        CRIRF=COUP(1,1)**2*HGZ(2,1)*COUP(3,1)**2/16.+
     &  COUP(1,1)*COUP(1,4)*HGZ(2,2)*COUP(3,1)*COUP(3,4)/4.+
     &  COUP(1,4)**2*HGZ(2,3)*COUP(3,4)**2
        IF(K(ILIN(1),2).GT.0) WT=(CLILF+CRIRF)*(PKK(1,4)**2+
     &  PKK(3,5)**2)+(CLIRF+CRILF)*(PKK(1,3)**2+PKK(4,5)**2)
        IF(K(ILIN(1),2).LT.0) WT=(CLILF+CRIRF)*(PKK(1,3)**2+
     &  PKK(4,5)**2)+(CLIRF+CRILF)*(PKK(1,4)**2+PKK(3,5)**2)
        WTMAX=(CLILF+CLIRF+CRILF+CRIRF)*
     &  ((PKK(1,3)+PKK(1,4))**2+(PKK(3,5)+PKK(4,5))**2)
 
      ELSEIF(ISUB.EQ.31) THEN
C...Angular weight for f + g -> f' + W+/- -> f' + 2 quarks/leptons.
        IF(K(ILIN(1),2).GT.0) WT=PKK(1,4)**2+PKK(3,5)**2
        IF(K(ILIN(1),2).LT.0) WT=PKK(1,3)**2+PKK(4,5)**2
        WTMAX=(PKK(1,3)+PKK(1,4))**2+(PKK(3,5)+PKK(4,5))**2
 
      ELSEIF(ISUB.EQ.71.OR.ISUB.EQ.72.OR.ISUB.EQ.73.OR.ISUB.EQ.76.OR.
     &ISUB.EQ.77) THEN
C...Angular weight for V_L1 + V_L2 -> V_L3 + V_L4 (V = Z/W).
        WT=16.*PKK(3,5)*PKK(4,6)
        WTMAX=SH**2
 
      ELSEIF(ISUB.EQ.110) THEN
C...Angular weight for f + f~ -> gamma + H0 -> gamma + X is isotropic.
        WT=1.
        WTMAX=1.
 
      ELSEIF(ISUB.EQ.141) THEN
        IF(IP.EQ.1.AND.(KDCY(1).EQ.1.OR.KDCY(1).EQ.2)) THEN
C...Angular weight for f + f~ -> gamma*/Z0/Z'0 -> 2 quarks/leptons.
C...Couplings of incoming flavour.
          KFAI=IABS(MINT(15))
          EI=KCHG(KFAI,1)/3.
          AI=SIGN(1.,EI+0.1)
          VI=AI-4.*EI*XWV
          KFAIC=1
          IF(KFAI.LE.10.AND.MOD(KFAI,2).EQ.0) KFAIC=2
          IF(KFAI.GT.10.AND.MOD(KFAI,2).NE.0) KFAIC=3
          IF(KFAI.GT.10.AND.MOD(KFAI,2).EQ.0) KFAIC=4
          VPI=PARU(119+2*KFAIC)
          API=PARU(120+2*KFAIC)
C...Couplings of final flavour.
          KFAF=IABS(KFL1(1))
          EF=KCHG(KFAF,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
          KFAFC=1
          IF(KFAF.LE.10.AND.MOD(KFAF,2).EQ.0) KFAFC=2
          IF(KFAF.GT.10.AND.MOD(KFAF,2).NE.0) KFAFC=3
          IF(KFAF.GT.10.AND.MOD(KFAF,2).EQ.0) KFAFC=4
          VPF=PARU(119+2*KFAFC)
          APF=PARU(120+2*KFAFC)
C...Asymmetry and weight.
          ASYM=2.*(EI*AI*VINT(112)*EF*AF+EI*API*VINT(113)*EF*APF+
     &    4.*VI*AI*VINT(114)*VF*AF+(VI*API+VPI*AI)*VINT(115)*
     &    (VF*APF+VPF*AF)+4.*VPI*API*VINT(116)*VPF*APF)/
     &    (EI**2*VINT(111)*EF**2+EI*VI*VINT(112)*EF*VF+
     &    EI*VPI*VINT(113)*EF*VPF+(VI**2+AI**2)*VINT(114)*
     &    (VF**2+AF**2)+(VI*VPI+AI*API)*VINT(115)*(VF*VPF+AF*APF)+
     &    (VPI**2+API**2)*VINT(116)*(VPF**2+APF**2))
          WT=1.+ASYM*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1))+CTHE(1)**2
          WTMAX=2.+ABS(ASYM)
        ELSEIF(IP.EQ.1.AND.IABS(KFL1(1)).EQ.24) THEN
C...Angular weight for f + f~ -> Z' -> W+ + W-.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          CCOS2=-(1./16.)*((1.-RM1-RM2)**2-4.*RM1*RM2)*
     &    (1.-2.*RM1-2.*RM2+RM1**2+RM2**2+10.*RM1*RM2)
          CFLAT=-CCOS2+0.5*(RM1+RM2)*(1.-2.*RM1-2.*RM2+(RM2-RM1)**2)
          WT=CFLAT+CCOS2*CTHE(1)**2
          WTMAX=CFLAT+MAX(0.,CCOS2)
        ELSEIF(IP.EQ.1.AND.(KFL1(1).EQ.25.OR.KFL1(1).EQ.35.OR.
     &  IABS(KFL1(1)).EQ.37)) THEN
C...Angular weight for f + f~ -> Z' -> H0 + A0, H'0 + A0, H+ + H-.
          WT=1.-CTHE(1)**2
          WTMAX=1.
        ELSEIF(IP.EQ.1.AND.KDCY(1).EQ.3) THEN
C...Angular weight for f + f~ -> Z' -> Z0 + H0.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          FLAM2=MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2)
          WT=1.+FLAM2*(1.-CTHE(1)**2)/(8.*RM1)
          WTMAX=1.+FLAM2/(8.*RM1)
        ELSEIF(MZPWP.EQ.0) THEN
C...Angular weight for f + f~ -> Z' -> W+ + W- -> 4 quarks/leptons
C...(W:s like if intermediate Z).
          D34=P(IREF(IP,IORD),5)**2
          D56=P(IREF(IP,3-IORD),5)**2
          DT=PKK(1,3)+PKK(1,4)+D34
          DU=PKK(1,5)+PKK(1,6)+D56
          FGK135=ABS(FGK(1,2,3,4,5,6)-FGK(1,2,5,6,3,4))
          FGK253=ABS(FGK(2,1,5,6,3,4)-FGK(2,1,3,4,5,6))
          WT=(COUP(1,3)*FGK135)**2+(COUP(1,4)*FGK253)**2
          WTMAX=4.*D34*D56*(COUP(1,3)**2+COUP(1,4)**2)*
     &    (DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU))
        ELSEIF(MZPWP.EQ.1) THEN
C...Angular weight for f + f~ -> Z' -> W+ + W- -> 4 quarks/leptons
C...(W:s approximately longitudinal, like if intermediate H).
          WT=16.*PKK(3,5)*PKK(4,6)
          WTMAX=SH**2
        ELSE
C...Angular weight for f + f~ -> Z' -> H+ + H-, Z0 + H0, H0 + A0,
C...H'0 + A0 -> 4 quarks/leptons.
          WT=1.
          WTMAX=1.
        ENDIF
 
      ELSEIF(ISUB.EQ.142) THEN
        IF(IP.EQ.1.AND.(KDCY(1).EQ.1.OR.KDCY(1).EQ.2)) THEN
C...Angular weight for f + f~' -> W'+/- -> 2 quarks/leptons.
          KFAI=IABS(MINT(15))
          KFAIC=1
          IF(KFAI.GT.10) KFAIC=2
          VI=PARU(129+2*KFAIC)
          AI=PARU(130+2*KFAIC)
          KFAF=IABS(KFL1(1))
          KFAFC=1
          IF(KFAF.GT.10) KFAFC=2
          VF=PARU(129+2*KFAFC)
          AF=PARU(130+2*KFAFC)
          ASYM=8.*VI*AI*VF*AF/((VI**2+AI**2)*(VF**2+AF**2))
          WT=1.+ASYM*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1))+CTHE(1)**2
          WTMAX=2.+ABS(ASYM)
        ELSEIF(IP.EQ.1.AND.IABS(KFL2(1)).EQ.23) THEN
C...Angular weight for f + f~' -> W'+/- -> W+/- + Z0.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          CCOS2=-(1./16.)*((1.-RM1-RM2)**2-4.*RM1*RM2)*
     &    (1.-2.*RM1-2.*RM2+RM1**2+RM2**2+10.*RM1*RM2)
          CFLAT=-CCOS2+0.5*(RM1+RM2)*(1.-2.*RM1-2.*RM2+(RM2-RM1)**2)
          WT=CFLAT+CCOS2*CTHE(1)**2
          WTMAX=CFLAT+MAX(0.,CCOS2)
        ELSEIF(IP.EQ.1.AND.KDCY(1).EQ.3) THEN
C...Angular weight for f + f~ -> W'+/- -> W+/- + H0.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          FLAM2=MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2)
          WT=1.+FLAM2*(1.-CTHE(1)**2)/(8.*RM1)
          WTMAX=1.+FLAM2/(8.*RM1)
        ELSEIF(MZPWP.EQ.0) THEN
C...Angular weight for f + f~' -> W' -> W + Z0 -> 4 quarks/leptons
C...(W/Z like if intermediate W).
          D34=P(IREF(IP,IORD),5)**2
          D56=P(IREF(IP,3-IORD),5)**2
          DT=PKK(1,3)+PKK(1,4)+D34
          DU=PKK(1,5)+PKK(1,6)+D56
          FGK135=ABS(FGK(1,2,3,4,5,6)-FGK(1,2,5,6,3,4))
          FGK136=ABS(FGK(1,2,3,4,6,5)-FGK(1,2,6,5,3,4))
          WT=(COUP(5,3)*FGK135)**2+(COUP(5,4)*FGK136)**2
          WTMAX=4.*D34*D56*(COUP(5,3)**2+COUP(5,4)**2)*
     &    (DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU))
        ELSEIF(MZPWP.EQ.1) THEN
C...Angular weight for f + f~' -> W' -> W + Z0 -> 4 quarks/leptons
C...(W/Z approximately longitudinal, like if intermediate H).
          WT=16.*PKK(3,5)*PKK(4,6)
          WTMAX=SH**2
        ELSE
C...Angular weight for f + f~ -> W' -> W + H0 -> whatever.
          WT=1.
          WTMAX=1.
        ENDIF
 
      ELSEIF(ISUB.EQ.145.OR.ISUB.EQ.162.OR.ISUB.EQ.163.OR.ISUB.EQ.164)
     &THEN
C...Isotropic decay of leptoquarks (assumed spin 0).
        WT=1.
        WTMAX=1.
 
      ELSEIF(ISUB.EQ.147.OR.ISUB.EQ.148) THEN
C...Decays of (spin 1/2) q* -> q + (g,gamma) or (Z0,W+-).
        SIDE=1.
        IF(MINT(16).EQ.21) SIDE=-1.
        IF(IP.EQ.1.AND.(KFL1(1).EQ.21.OR.KFL1(1).EQ.22)) THEN
          WT=1.+SIDE*CTHE(1)
          WTMAX=2.
        ELSEIF(IP.EQ.1) THEN
          RM1=P(NSD(1)+1,5)**2/SH
          WT=1.+SIDE*CTHE(1)*(1.-0.5*RM1)/(1.+0.5*RM1)
          WTMAX=1.+(1.-0.5*RM1)/(1.+0.5*RM1)
        ELSE
C...W/Z decay assumed isotropic, since not known.
          WT=1.
          WTMAX=1.
        ENDIF
 
      ELSEIF(ISUB.EQ.149) THEN
C...Isotropic decay of techni-eta.
        WT=1.
        WTMAX=1.
 
C...Obtain correct angular distribution by rejection techniques.
      ELSE
        WT=1.
        WTMAX=1.
      ENDIF
      IF(WT.LT.RLU(0)*WTMAX) GOTO 240
 
C...Construct massive four-vectors using angles chosen. Mark decayed
C...resonances, add documentation lines. Shower evolution.
  400 DO 470 JT=1,JTMAX
      IF(KDCY(JT).EQ.0) GOTO 470
      ID=IREF(IP,JT)
      IF(ISET(ISUB).NE.6.OR.JT.NE.3) THEN
        DO 410 J=1,5
        DPMO(J)=P(ID,J)
  410   CONTINUE
        DPMO(4)=SQRT(DPMO(1)**2+DPMO(2)**2+DPMO(3)**2+DPMO(5)**2)
        CALL LUDBRB(NSD(JT)+1,NSD(JT)+2,ACOS(CTHE(JT)),PHI(JT),
     &  DPMO(1)/DPMO(4),DPMO(2)/DPMO(4),DPMO(3)/DPMO(4))
      ELSE
C...Z + q + q~ : angles already known, in rest frame of system.
        DO 420 J=1,3
        DBEZQQ(J)=(P(ID,J)+P(ID+1,J)+P(ID+2,J))/(P(ID,4)+P(ID+1,4)+
     &  P(ID+2,4))
  420   CONTINUE
        K(N+1,1)=1
        DO 430 J=1,5
        P(N+1,J)=P(ID,J)
  430   CONTINUE
        CALL LUDBRB(N+1,N+1,0.,0.,-DBEZQQ(1),-DBEZQQ(2),-DBEZQQ(3))
        PHIZQQ=ULANGL(P(N+1,1),P(N+1,2))
        THEZQQ=ULANGL(P(N+1,3),SQRT(P(N+1,1)**2+P(N+1,2)**2))
        CALL LUDBRB(NSD(JT)+1,NSD(JT)+2,ACOS(VINT(81)),VINT(82),
     &  0D0,0D0,DBLE(SQRT(P(N+1,1)**2+P(N+1,2)**2+P(N+1,3)**2)/
     &  P(N+1,4)))
        CALL LUDBRB(NSD(JT)+1,NSD(JT)+2,THEZQQ,PHIZQQ,DBEZQQ(1),
     &  DBEZQQ(2),DBEZQQ(3))
      ENDIF
      K(ID,1)=K(ID,1)+10
      KFA=IABS(K(ID,2))
      IF(KFA.EQ.38.OR.KFA.EQ.39.OR.((MSTP(6).EQ.1.OR.MSTP(49).GE.1).AND.
     &(KFA.EQ.7.OR.KFA.EQ.8)).OR.(MSTP(48).GE.1.AND.KFA.EQ.6)) THEN
C...Do not kill colour flow through techni-eta, t, leptoquark or q*!
      ELSE
        K(ID,4)=NSD(JT)+1
        K(ID,5)=NSD(JT)+2
      ENDIF
      IDOC=MINT(83)+MINT(4)
      DO 450 I=NSD(JT)+1,NSD(JT)+2
      I1=MINT(83)+MINT(4)+1
      K(I,3)=I1
      IF(MSTP(128).GE.1) K(I,3)=ID
      IF(MSTP(128).LE.1.AND.MINT(4).LT.MSTP(126)) THEN
        MINT(4)=MINT(4)+1
        K(I1,1)=21
        K(I1,2)=K(I,2)
        K(I1,3)=IREF(IP,JT+3)
        DO 440 J=1,5
        P(I1,J)=P(I,J)
  440   CONTINUE
      ENDIF
  450 CONTINUE
C...Shower - top currently special case.
      NSHBEF=N
      IF(MSTP(71).GE.1.AND.(KDCY(JT).LE.2.OR.KFA.EQ.6.OR.KFA.EQ.7.OR.
     &KFA.EQ.8)) CALL LUSHOW(NSD(JT)+1,NSD(JT)+2,P(ID,5))
      NSHAFT=N
 
C...Check if new resonances were produced.
      KNSDA1=IABS(K(NSD(JT)+1,2))
      KNSDA2=IABS(K(NSD(JT)+2,2))
      IF(KNSDA1.EQ.6.OR.KNSDA2.EQ.6.OR.KNSDA1.EQ.7.OR.KNSDA2.EQ.7.OR.
     &KNSDA1.EQ.8.OR.KNSDA2.EQ.8) THEN
        NSD1=0
        NSD2=0
        DO 460 I=NSD(JT)+1,N
        IF(K(I,1).LT.10.AND.K(I,2).EQ.K(NSD(JT)+1,2)) NSD1=I
        IF(K(I,1).LT.10.AND.K(I,2).EQ.K(NSD(JT)+2,2)) NSD2=I
  460   CONTINUE
        IF(NSD1.NE.0.AND.NSD2.NE.0) THEN
          NP=NP+1
          IREF(NP,1)=NSD1
          IREF(NP,2)=NSD2
          IREF(NP,3)=0
          IREF(NP,4)=IDOC+1
          IREF(NP,5)=IDOC+2
          IREF(NP,6)=0
          IREF(NP,7)=K(IREF(IP,JT),2)
          IREF(NP,8)=IREF(IP,JT)
        ENDIF
      ELSEIF(KDCY(JT).EQ.3.OR.KFA.EQ.6.OR.KFA.EQ.7.OR.KFA.EQ.8) THEN
        NP=NP+1
        IREF(NP,1)=NSD(JT)+1
        IREF(NP,2)=NSD(JT)+2
        IF(NSHAFT-NSHBEF.GT.0) THEN
          IREF(NP,1)=NSHBEF+2
          IREF(NP,2)=NSHBEF+3
        ENDIF
        IREF(NP,3)=0
        IREF(NP,4)=IDOC+1
        IREF(NP,5)=IDOC+2
        IREF(NP,6)=0
        IREF(NP,7)=K(IREF(IP,JT),2)
        IREF(NP,8)=IREF(IP,JT)
      ENDIF
  470 CONTINUE
 
C...Fill information for 2 -> 1 -> 2. Loop back if needed.
      IF(JTMAX.EQ.1.AND.KDCY(1).NE.0) THEN
        MINT(7)=MINT(83)+6+2*ISET(ISUB)
        MINT(8)=MINT(83)+7+2*ISET(ISUB)
        MINT(25)=KFL1(1)
        MINT(26)=KFL2(1)
        VINT(23)=CTHE(1)
        RM3=P(N-1,5)**2/SH
        RM4=P(N,5)**2/SH
        BE34=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4))
        VINT(45)=-0.5*SH*(1.-RM3-RM4-BE34*CTHE(1))
        VINT(46)=-0.5*SH*(1.-RM3-RM4+BE34*CTHE(1))
        VINT(48)=0.25*SH*BE34**2*MAX(0.,1.-CTHE(1)**2)
        VINT(47)=SQRT(VINT(48))
      ENDIF
  480 IF(IP.LT.NP) GOTO 130
 
      RETURN
      END
