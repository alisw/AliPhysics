 
C*********************************************************************
 
      SUBROUTINE PYSCAT
 
C...Finds outgoing flavours and event type; sets up the kinematics
C...and colour flow of the hard scattering.
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3)
      COMMON/PYUPPR/NUP,KUP(20,7),PUP(20,5),NFUP,IFUP(10,2),Q2UP(0:10)
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT3/,/PYINT4/,
     &/PYINT5/,/PYUPPR/
      DIMENSION WDTP(0:40),WDTE(0:40,0:5),PMQ(2),Z(2),CTHE(2),PHI(2),
     &KUPPO(20),VINTSV(41:66)
      SAVE VINTSV
 
C...Read out process.
      ISUB=MINT(1)
      ISUBSV=ISUB
 
C...Restore information for low-pT processes.
      IF(ISUB.EQ.95.AND.MINT(57).GE.1) THEN
        DO 100 J=41,66
  100   VINT(J)=VINTSV(J)
      ENDIF
 
C...Convert H' or A process into equivalent H one.
      IHIGG=1
      KFHIGG=25
      IF((ISUB.GE.151.AND.ISUB.LE.160).OR.(ISUB.GE.171.AND.
     &ISUB.LE.190)) THEN
        IHIGG=2
        IF(MOD(ISUB-1,10).GE.5) IHIGG=3
        KFHIGG=33+IHIGG
        IF(ISUB.EQ.151.OR.ISUB.EQ.156) ISUB=3
        IF(ISUB.EQ.152.OR.ISUB.EQ.157) ISUB=102
        IF(ISUB.EQ.153.OR.ISUB.EQ.158) ISUB=103
        IF(ISUB.EQ.171.OR.ISUB.EQ.176) ISUB=24
        IF(ISUB.EQ.172.OR.ISUB.EQ.177) ISUB=26
        IF(ISUB.EQ.173.OR.ISUB.EQ.178) ISUB=123
        IF(ISUB.EQ.174.OR.ISUB.EQ.179) ISUB=124
        IF(ISUB.EQ.181.OR.ISUB.EQ.186) ISUB=121
        IF(ISUB.EQ.182.OR.ISUB.EQ.187) ISUB=122
      ENDIF
 
C...Choice of subprocess, number of documentation lines.
      IDOC=6+ISET(ISUB)
      IF(ISUB.EQ.95) IDOC=8
      IF(ISET(ISUB).EQ.5.OR.ISET(ISUB).EQ.6) IDOC=9
      IF(ISET(ISUB).EQ.11) IDOC=4+NUP
      MINT(3)=IDOC-6
      IF(IDOC.GE.9.AND.ISET(ISUB).LE.4) IDOC=IDOC+2
      MINT(4)=IDOC
      IPU1=MINT(84)+1
      IPU2=MINT(84)+2
      IPU3=MINT(84)+3
      IPU4=MINT(84)+4
      IPU5=MINT(84)+5
      IPU6=MINT(84)+6
 
C...Reset K, P and V vectors. Store incoming particles.
      DO 120 JT=1,MSTP(126)+20
      I=MINT(83)+JT
      DO 110 J=1,5
      K(I,J)=0
      P(I,J)=0.
      V(I,J)=0.
  110 CONTINUE
  120 CONTINUE
      DO 140 JT=1,2
      I=MINT(83)+JT
      K(I,1)=21
      K(I,2)=MINT(10+JT)
      DO 130 J=1,5
      P(I,J)=VINT(285+5*JT+J)
  130 CONTINUE
  140 CONTINUE
      MINT(6)=2
      KFRES=0
 
C...Store incoming partons in their CM-frame.
      SH=VINT(44)
      SHR=SQRT(SH)
      SHP=VINT(26)*VINT(2)
      SHPR=SQRT(SHP)
      SHUSER=SHR
      IF(ISET(ISUB).GE.3.AND.ISET(ISUB).LE.5) SHUSER=SHPR
      DO 150 JT=1,2
      I=MINT(84)+JT
      K(I,1)=14
      K(I,2)=MINT(14+JT)
      K(I,3)=MINT(83)+2+JT
      P(I,3)=0.5*SHUSER*(-1.)**(JT-1)
      P(I,4)=0.5*SHUSER
  150 CONTINUE
 
C...Copy incoming partons to documentation lines.
      DO 170 JT=1,2
      I1=MINT(83)+4+JT
      I2=MINT(84)+JT
      K(I1,1)=21
      K(I1,2)=K(I2,2)
      K(I1,3)=I1-2
      DO 160 J=1,5
      P(I1,J)=P(I2,J)
  160 CONTINUE
  170 CONTINUE
 
C...Choose new quark/lepton flavour for relevant annihilation graphs.
      IF(ISUB.EQ.12.OR.ISUB.EQ.53.OR.ISUB.EQ.54.OR.ISUB.EQ.58) THEN
        IGLGA=21
        IF(ISUB.EQ.58) IGLGA=22
        CALL PYWIDT(IGLGA,SH,WDTP,WDTE)
  180   RKFL=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))*RLU(0)
        DO 190 I=1,MDCY(IGLGA,3)
        KFLF=KFDP(I+MDCY(IGLGA,2)-1,1)
        RKFL=RKFL-(WDTE(I,1)+WDTE(I,2)+WDTE(I,4))
        IF(RKFL.LE.0.) GOTO 200
  190   CONTINUE
  200   CONTINUE
        IF(ISUB.EQ.12.AND.MSTP(5).EQ.1.AND.IABS(MINT(15)).LE.2.AND.
     &  IABS(KFLF).GE.3) THEN
          FACQQB=VINT(58)**2*4./9.*(VINT(45)**2+VINT(46)**2)/
     &    VINT(44)**2
          FACCIB=VINT(46)**2/PARU(155)**4
          IF(FACQQB/(FACQQB+FACCIB).LT.RLU(0)) GOTO 180
        ELSEIF(ISUB.EQ.54) THEN
          IF((KCHG(IABS(KFLF),1)/2.)**2.LT.RLU(0)) GOTO 180
        ELSEIF(ISUB.EQ.58) THEN
          IF((KCHG(IABS(KFLF),1)/3.)**2.LT.RLU(0)) GOTO 180
        ENDIF
      ENDIF
 
C...Final state flavours and colour flow: default values.
      JS=1
      MINT(21)=MINT(15)
      MINT(22)=MINT(16)
      MINT(23)=0
      MINT(24)=0
      KCC=20
      KCS=ISIGN(1,MINT(15))
 
      IF(ISET(ISUB).EQ.11) THEN
C...User-defined processes: find products.
        IRUP=0
        DO 210 IUP=3,NUP
        IF(KUP(IUP,1).NE.1) THEN
        ELSEIF(IRUP.LE.5) THEN
          IRUP=IRUP+1
          MINT(20+IRUP)=KUP(IUP,2)
        ENDIF
  210   CONTINUE
 
      ELSEIF(ISUB.LE.10) THEN
      IF(ISUB.EQ.1) THEN
C...f + f~ -> gamma*/Z0.
        KFRES=23
 
      ELSEIF(ISUB.EQ.2) THEN
C...f + f~' -> W+/- .
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        KFRES=ISIGN(24,KCH1+KCH2)
 
      ELSEIF(ISUB.EQ.3) THEN
C...f + f~ -> H0 (or H'0, or A0).
        KFRES=KFHIGG
 
      ELSEIF(ISUB.EQ.4) THEN
C...gamma + W+/- -> W+/-.
 
      ELSEIF(ISUB.EQ.5) THEN
C...Z0 + Z0 -> H0.
        XH=SH/SHP
        MINT(21)=MINT(15)
        MINT(22)=MINT(16)
        PMQ(1)=ULMASS(MINT(21))
        PMQ(2)=ULMASS(MINT(22))
  220   JT=INT(1.5+RLU(0))
        ZMIN=2.*PMQ(JT)/SHPR
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        Z(JT)=ZMIN+(ZMAX-ZMIN)*RLU(0)
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.
     &  (1.-XH)**2/(4.*XH)*RLU(0)) GOTO 220
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 220
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(JT)**2)/(Z(JT)*SHP)
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))
        Z(3-JT)=1.-XH/(1.-Z(JT))
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 220
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP)
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))
        PHIR=PARU(2)*RLU(0)
        CPHI=COS(PHIR)
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI
        Z1=2.-Z(JT)
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP)
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*
     &  PMQ(3-JT)**2/SHP))
        ZMIN=2.*PMQ(3-JT)/SHPR
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 220
        KCC=22
        KFRES=25
 
      ELSEIF(ISUB.EQ.6) THEN
C...Z0 + W+/- -> W+/-.
 
      ELSEIF(ISUB.EQ.7) THEN
C...W+ + W- -> Z0.
 
      ELSEIF(ISUB.EQ.8) THEN
C...W+ + W- -> H0.
        XH=SH/SHP
  230   DO 260 JT=1,2
        I=MINT(14+JT)
        IA=IABS(I)
        IF(IA.LE.10) THEN
          RVCKM=VINT(180+I)*RLU(0)
          DO 240 J=1,MSTP(1)
          IB=2*J-1+MOD(IA,2)
          IPM=(5-ISIGN(1,I))/2
          IDC=J+MDCY(IA,2)+2
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 240
          MINT(20+JT)=ISIGN(IB,I)
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
          IF(RVCKM.LE.0.) GOTO 250
  240     CONTINUE
        ELSE
          IB=2*((IA+1)/2)-1+MOD(IA,2)
          MINT(20+JT)=ISIGN(IB,I)
        ENDIF
  250   PMQ(JT)=ULMASS(MINT(20+JT))
  260   CONTINUE
        JT=INT(1.5+RLU(0))
        ZMIN=2.*PMQ(JT)/SHPR
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(ZMIN.GE.ZMAX) GOTO 230
        Z(JT)=ZMIN+(ZMAX-ZMIN)*RLU(0)
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.
     &  (1.-XH)**2/(4.*XH)*RLU(0)) GOTO 230
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 230
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(24,1)**2-PMQ(JT)**2)/(Z(JT)*SHP)
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))
        Z(3-JT)=1.-XH/(1.-Z(JT))
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 230
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(24,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP)
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))
        PHIR=PARU(2)*RLU(0)
        CPHI=COS(PHIR)
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI
        Z1=2.-Z(JT)
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP)
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*
     &  PMQ(3-JT)**2/SHP))
        ZMIN=2.*PMQ(3-JT)/SHPR
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 230
        KCC=22
        KFRES=25
 
      ELSEIF(ISUB.EQ.10) THEN
C...f + f' -> f + f' (gamma/Z/W exchange); th = (p(f)-p(f))**2.
        IF(MINT(2).EQ.1) THEN
          KCC=22
        ELSE
C...W exchange: need to mix flavours according to CKM matrix.
          DO 280 JT=1,2
          I=MINT(14+JT)
          IA=IABS(I)
          IF(IA.LE.10) THEN
            RVCKM=VINT(180+I)*RLU(0)
            DO 270 J=1,MSTP(1)
            IB=2*J-1+MOD(IA,2)
            IPM=(5-ISIGN(1,I))/2
            IDC=J+MDCY(IA,2)+2
            IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 270
            MINT(20+JT)=ISIGN(IB,I)
            RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
            IF(RVCKM.LE.0.) GOTO 280
  270       CONTINUE
          ELSE
            IB=2*((IA+1)/2)-1+MOD(IA,2)
            MINT(20+JT)=ISIGN(IB,I)
          ENDIF
  280     CONTINUE
          KCC=22
        ENDIF
      ENDIF
 
      ELSEIF(ISUB.LE.20) THEN
      IF(ISUB.EQ.11) THEN
C...f + f' -> f + f' (g exchange); th = (p(f)-p(f))**2.
        KCC=MINT(2)
        IF(MINT(15)*MINT(16).LT.0) KCC=KCC+2
 
      ELSEIF(ISUB.EQ.12) THEN
C...f + f~ -> f' + f~'; th = (p(f)-p(f'))**2.
        MINT(21)=ISIGN(KFLF,MINT(15))
        MINT(22)=-MINT(21)
        KCC=4
 
      ELSEIF(ISUB.EQ.13) THEN
C...f + f~ -> g + g; th arbitrary.
        MINT(21)=21
        MINT(22)=21
        KCC=MINT(2)+4
 
      ELSEIF(ISUB.EQ.14) THEN
C...f + f~ -> g + gamma; th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(20+JS)=21
        MINT(23-JS)=22
        KCC=17+JS
 
      ELSEIF(ISUB.EQ.15) THEN
C...f + f~ -> g + Z0; th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(20+JS)=21
        MINT(23-JS)=23
        KCC=17+JS
 
      ELSEIF(ISUB.EQ.16) THEN
C...f + f~' -> g + W+/-; th = (p(f)-p(W-))**2 or (p(f~')-p(W+))**2.
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        IF(MINT(15)*(KCH1+KCH2).LT.0) JS=2
        MINT(20+JS)=21
        MINT(23-JS)=ISIGN(24,KCH1+KCH2)
        KCC=17+JS
 
      ELSEIF(ISUB.EQ.17) THEN
C...f + f~ -> g + H0; th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(20+JS)=21
        MINT(23-JS)=25
        KCC=17+JS
 
      ELSEIF(ISUB.EQ.18) THEN
C...f + f~ -> gamma + gamma; th arbitrary.
        MINT(21)=22
        MINT(22)=22
 
      ELSEIF(ISUB.EQ.19) THEN
C...f + f~ -> gamma + Z0; th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(20+JS)=22
        MINT(23-JS)=23
 
      ELSEIF(ISUB.EQ.20) THEN
C...f + f~' -> gamma + W+/-; th = (p(f)-p(W-))**2 or (p(f~')-p(W+))**2.
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        IF(MINT(15)*(KCH1+KCH2).LT.0) JS=2
        MINT(20+JS)=22
        MINT(23-JS)=ISIGN(24,KCH1+KCH2)
      ENDIF
 
      ELSEIF(ISUB.LE.30) THEN
      IF(ISUB.EQ.21) THEN
C...f + f~ -> gamma + H0; th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(20+JS)=22
        MINT(23-JS)=25
 
      ELSEIF(ISUB.EQ.22) THEN
C...f + f~ -> Z0 + Z0; th arbitrary.
        MINT(21)=23
        MINT(22)=23
 
      ELSEIF(ISUB.EQ.23) THEN
C...f + f~' -> Z0 + W+/-; th = (p(f)-p(W-))**2 or (p(f~')-p(W+))**2.
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        IF(MINT(15)*(KCH1+KCH2).LT.0) JS=2
        MINT(20+JS)=23
        MINT(23-JS)=ISIGN(24,KCH1+KCH2)
 
      ELSEIF(ISUB.EQ.24) THEN
C...f + f~ -> Z0 + H0 (or H'0, or A0); th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(20+JS)=23
        MINT(23-JS)=KFHIGG
 
      ELSEIF(ISUB.EQ.25) THEN
C...f + f~ -> W+ + W-; th = (p(f)-p(W-))**2.
        MINT(21)=-ISIGN(24,MINT(15))
        MINT(22)=-MINT(21)
 
      ELSEIF(ISUB.EQ.26) THEN
C...f + f~' -> W+/- + H0 (or H'0, or A0);
C...th = (p(f)-p(W-))**2 or (p(f~')-p(W+))**2.
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        IF(MINT(15)*(KCH1+KCH2).GT.0) JS=2
        MINT(20+JS)=ISIGN(24,KCH1+KCH2)
        MINT(23-JS)=KFHIGG
 
      ELSEIF(ISUB.EQ.27) THEN
C...f + f~ -> H0 + H0.
 
      ELSEIF(ISUB.EQ.28) THEN
C...f + g -> f + g; th = (p(f)-p(f))**2.
        KCC=MINT(2)+6
        IF(MINT(15).EQ.21) KCC=KCC+2
        IF(MINT(15).NE.21) KCS=ISIGN(1,MINT(15))
        IF(MINT(16).NE.21) KCS=ISIGN(1,MINT(16))
 
      ELSEIF(ISUB.EQ.29) THEN
C...f + g -> f + gamma; th = (p(f)-p(f))**2.
        IF(MINT(15).EQ.21) JS=2
        MINT(23-JS)=22
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))
 
      ELSEIF(ISUB.EQ.30) THEN
C...f + g -> f + Z0; th = (p(f)-p(f))**2.
        IF(MINT(15).EQ.21) JS=2
        MINT(23-JS)=23
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))
      ENDIF
 
      ELSEIF(ISUB.LE.40) THEN
      IF(ISUB.EQ.31) THEN
C...f + g -> f' + W+/-; th = (p(f)-p(f'))**2; choose flavour f'.
        IF(MINT(15).EQ.21) JS=2
        I=MINT(14+JS)
        IA=IABS(I)
        MINT(23-JS)=ISIGN(24,KCHG(IA,1)*I)
        RVCKM=VINT(180+I)*RLU(0)
        DO 290 J=1,MSTP(1)
        IB=2*J-1+MOD(IA,2)
        IPM=(5-ISIGN(1,I))/2
        IDC=J+MDCY(IA,2)+2
        IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 290
        MINT(20+JS)=ISIGN(IB,I)
        RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
        IF(RVCKM.LE.0.) GOTO 300
  290   CONTINUE
  300   KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))
 
      ELSEIF(ISUB.EQ.32) THEN
C...f + g -> f + H0; th = (p(f)-p(f))**2.
        IF(MINT(15).EQ.21) JS=2
        MINT(23-JS)=25
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))
 
      ELSEIF(ISUB.EQ.33) THEN
C...f + gamma -> f + g; th=(p(f)-p(f))**2.
        IF(MINT(15).EQ.22) JS=2
        MINT(23-JS)=21
        KCC=24+JS
        KCS=ISIGN(1,MINT(14+JS))
 
      ELSEIF(ISUB.EQ.34) THEN
C...f + gamma -> f + gamma; th=(p(f)-p(f))**2.
        IF(MINT(15).EQ.22) JS=2
        KCC=22
        KCS=ISIGN(1,MINT(14+JS))
 
      ELSEIF(ISUB.EQ.35) THEN
C...f + gamma -> f + Z0; th=(p(f)-p(f))**2.
        IF(MINT(15).EQ.22) JS=2
        MINT(23-JS)=23
        KCC=22
 
      ELSEIF(ISUB.EQ.36) THEN
C...f + gamma -> f' + W+/-; th=(p(f)-p(f'))**2.
        IF(MINT(15).EQ.22) JS=2
        I=MINT(14+JS)
        IA=IABS(I)
        MINT(23-JS)=ISIGN(24,KCHG(IA,1)*I)
        IF(IA.LE.10) THEN
          RVCKM=VINT(180+I)*RLU(0)
          DO 310 J=1,MSTP(1)
          IB=2*J-1+MOD(IA,2)
          IPM=(5-ISIGN(1,I))/2
          IDC=J+MDCY(IA,2)+2
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 310
          MINT(20+JS)=ISIGN(IB,I)
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
          IF(RVCKM.LE.0.) GOTO 320
  310     CONTINUE
        ELSE
          IB=2*((IA+1)/2)-1+MOD(IA,2)
          MINT(20+JS)=ISIGN(IB,I)
        ENDIF
  320   KCC=22
 
      ELSEIF(ISUB.EQ.37) THEN
C...f + gamma -> f + H0.
 
      ELSEIF(ISUB.EQ.38) THEN
C...f + Z0 -> f + g.
 
      ELSEIF(ISUB.EQ.39) THEN
C...f + Z0 -> f + gamma.
 
      ELSEIF(ISUB.EQ.40) THEN
C...f + Z0 -> f + Z0.
      ENDIF
 
      ELSEIF(ISUB.LE.50) THEN
      IF(ISUB.EQ.41) THEN
C...f + Z0 -> f' + W+/-.
 
      ELSEIF(ISUB.EQ.42) THEN
C...f + Z0 -> f + H0.
 
      ELSEIF(ISUB.EQ.43) THEN
C...f + W+/- -> f' + g.
 
      ELSEIF(ISUB.EQ.44) THEN
C...f + W+/- -> f' + gamma.
 
      ELSEIF(ISUB.EQ.45) THEN
C...f + W+/- -> f' + Z0.
 
      ELSEIF(ISUB.EQ.46) THEN
C...f + W+/- -> f' + W+/-.
 
      ELSEIF(ISUB.EQ.47) THEN
C...f + W+/- -> f' + H0.
 
      ELSEIF(ISUB.EQ.48) THEN
C...f + H0 -> f + g.
 
      ELSEIF(ISUB.EQ.49) THEN
C...f + H0 -> f + gamma.
 
      ELSEIF(ISUB.EQ.50) THEN
C...f + H0 -> f + Z0.
      ENDIF
 
      ELSEIF(ISUB.LE.60) THEN
      IF(ISUB.EQ.51) THEN
C...f + H0 -> f' + W+/-.
 
      ELSEIF(ISUB.EQ.52) THEN
C...f + H0 -> f + H0.
 
      ELSEIF(ISUB.EQ.53) THEN
C...g + g -> f + f~; th arbitrary.
        KCS=(-1)**INT(1.5+RLU(0))
        MINT(21)=ISIGN(KFLF,KCS)
        MINT(22)=-MINT(21)
        KCC=MINT(2)+10
 
      ELSEIF(ISUB.EQ.54) THEN
C...g + gamma -> f + f~; th arbitrary.
        KCS=(-1)**INT(1.5+RLU(0))
        MINT(21)=ISIGN(KFLF,KCS)
        MINT(22)=-MINT(21)
        KCC=27
        IF(MINT(16).EQ.21) KCC=28
 
      ELSEIF(ISUB.EQ.55) THEN
C...g + Z0 -> f + f~.
 
      ELSEIF(ISUB.EQ.56) THEN
C...g + W+/- -> f + f~'.
 
      ELSEIF(ISUB.EQ.57) THEN
C...g + H0 -> f + f~.
 
      ELSEIF(ISUB.EQ.58) THEN
C...gamma + gamma -> f + f~; th arbitrary.
        KCS=(-1)**INT(1.5+RLU(0))
        MINT(21)=ISIGN(KFLF,KCS)
        MINT(22)=-MINT(21)
        KCC=21
 
      ELSEIF(ISUB.EQ.59) THEN
C...gamma + Z0 -> f + f~.
 
      ELSEIF(ISUB.EQ.60) THEN
C...gamma + W+/- -> f + f~'.
      ENDIF
 
      ELSEIF(ISUB.LE.70) THEN
      IF(ISUB.EQ.61) THEN
C...gamma + H0 -> f + f~.
 
      ELSEIF(ISUB.EQ.62) THEN
C...Z0 + Z0 -> f + f~.
 
      ELSEIF(ISUB.EQ.63) THEN
C...Z0 + W+/- -> f + f~'.
 
      ELSEIF(ISUB.EQ.64) THEN
C...Z0 + H0 -> f + f~.
 
      ELSEIF(ISUB.EQ.65) THEN
C...W+ + W- -> f + f~.
 
      ELSEIF(ISUB.EQ.66) THEN
C...W+/- + H0 -> f + f~'.
 
      ELSEIF(ISUB.EQ.67) THEN
C...H0 + H0 -> f + f~.
 
      ELSEIF(ISUB.EQ.68) THEN
C...g + g -> g + g; th arbitrary.
        KCC=MINT(2)+12
        KCS=(-1)**INT(1.5+RLU(0))
 
      ELSEIF(ISUB.EQ.69) THEN
C...gamma + gamma -> W+ + W-; th arbitrary.
        MINT(21)=24
        MINT(22)=-24
        KCC=21
 
      ELSEIF(ISUB.EQ.70) THEN
C...gamma + W+/- -> Z0 + W+/-; th=(p(W)-p(W))**2.
        IF(MINT(15).EQ.22) MINT(21)=23
        IF(MINT(16).EQ.22) MINT(22)=23
        KCC=21
      ENDIF
 
      ELSEIF(ISUB.LE.80) THEN
      IF(ISUB.EQ.71.OR.ISUB.EQ.72) THEN
C...Z0 + Z0 -> Z0 + Z0; Z0 + Z0 -> W+ + W-.
        XH=SH/SHP
        MINT(21)=MINT(15)
        MINT(22)=MINT(16)
        PMQ(1)=ULMASS(MINT(21))
        PMQ(2)=ULMASS(MINT(22))
  330   JT=INT(1.5+RLU(0))
        ZMIN=2.*PMQ(JT)/SHPR
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        Z(JT)=ZMIN+(ZMAX-ZMIN)*RLU(0)
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.
     &  (1.-XH)**2/(4.*XH)*RLU(0)) GOTO 330
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 330
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(JT)**2)/(Z(JT)*SHP)
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))
        Z(3-JT)=1.-XH/(1.-Z(JT))
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 330
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP)
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))
        PHIR=PARU(2)*RLU(0)
        CPHI=COS(PHIR)
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI
        Z1=2.-Z(JT)
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP)
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*
     &  PMQ(3-JT)**2/SHP))
        ZMIN=2.*PMQ(3-JT)/SHPR
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 330
        KCC=22
 
      ELSEIF(ISUB.EQ.73) THEN
C...Z0 + W+/- -> Z0 + W+/-.
        JS=MINT(2)
        XH=SH/SHP
  340   JT=3-MINT(2)
        I=MINT(14+JT)
        IA=IABS(I)
        IF(IA.LE.10) THEN
          RVCKM=VINT(180+I)*RLU(0)
          DO 350 J=1,MSTP(1)
          IB=2*J-1+MOD(IA,2)
          IPM=(5-ISIGN(1,I))/2
          IDC=J+MDCY(IA,2)+2
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 350
          MINT(20+JT)=ISIGN(IB,I)
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
          IF(RVCKM.LE.0.) GOTO 360
  350     CONTINUE
        ELSE
          IB=2*((IA+1)/2)-1+MOD(IA,2)
          MINT(20+JT)=ISIGN(IB,I)
        ENDIF
  360   PMQ(JT)=ULMASS(MINT(20+JT))
        MINT(23-JT)=MINT(17-JT)
        PMQ(3-JT)=ULMASS(MINT(23-JT))
        JT=INT(1.5+RLU(0))
        ZMIN=2.*PMQ(JT)/SHPR
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(ZMIN.GE.ZMAX) GOTO 340
        Z(JT)=ZMIN+(ZMAX-ZMIN)*RLU(0)
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.
     &  (1.-XH)**2/(4.*XH)*RLU(0)) GOTO 340
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 340
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(JT)**2)/(Z(JT)*SHP)
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))
        Z(3-JT)=1.-XH/(1.-Z(JT))
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 340
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(23,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP)
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))
        PHIR=PARU(2)*RLU(0)
        CPHI=COS(PHIR)
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI
        Z1=2.-Z(JT)
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP)
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*
     &  PMQ(3-JT)**2/SHP))
        ZMIN=2.*PMQ(3-JT)/SHPR
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 340
        KCC=22
 
      ELSEIF(ISUB.EQ.74) THEN
C...Z0 + H0 -> Z0 + H0.
 
      ELSEIF(ISUB.EQ.75) THEN
C...W+ + W- -> gamma + gamma.
 
      ELSEIF(ISUB.EQ.76.OR.ISUB.EQ.77) THEN
C...W+ + W- -> Z0 + Z0; W+ + W- -> W+ + W-.
        XH=SH/SHP
  370   DO 400 JT=1,2
        I=MINT(14+JT)
        IA=IABS(I)
        IF(IA.LE.10) THEN
          RVCKM=VINT(180+I)*RLU(0)
          DO 380 J=1,MSTP(1)
          IB=2*J-1+MOD(IA,2)
          IPM=(5-ISIGN(1,I))/2
          IDC=J+MDCY(IA,2)+2
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 380
          MINT(20+JT)=ISIGN(IB,I)
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
          IF(RVCKM.LE.0.) GOTO 390
  380     CONTINUE
        ELSE
          IB=2*((IA+1)/2)-1+MOD(IA,2)
          MINT(20+JT)=ISIGN(IB,I)
        ENDIF
  390   PMQ(JT)=ULMASS(MINT(20+JT))
  400   CONTINUE
        JT=INT(1.5+RLU(0))
        ZMIN=2.*PMQ(JT)/SHPR
        ZMAX=1.-PMQ(3-JT)/SHPR-(SH-PMQ(JT)**2)/(SHPR*(SHPR-PMQ(3-JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(ZMIN.GE.ZMAX) GOTO 370
        Z(JT)=ZMIN+(ZMAX-ZMIN)*RLU(0)
        IF(-1.+(1.+XH)/(1.-Z(JT))-XH/(1.-Z(JT))**2.LT.
     &  (1.-XH)**2/(4.*XH)*RLU(0)) GOTO 370
        SQC1=1.-4.*PMQ(JT)**2/(Z(JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 370
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(24,1)**2-PMQ(JT)**2)/(Z(JT)*SHP)
        CTHE(JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1
        CTHE(JT)=MIN(1.,MAX(-1.,CTHE(JT)))
        Z(3-JT)=1.-XH/(1.-Z(JT))
        SQC1=1.-4.*PMQ(3-JT)**2/(Z(3-JT)**2*SHP)
        IF(SQC1.LT.1.E-8) GOTO 370
        C1=SQRT(SQC1)
        C2=1.+2.*(PMAS(24,1)**2-PMQ(3-JT)**2)/(Z(3-JT)*SHP)
        CTHE(3-JT)=(C2-(C2**2-C1**2)/(C2+(2.*RLU(0)-1.)*C1))/C1
        CTHE(3-JT)=MIN(1.,MAX(-1.,CTHE(3-JT)))
        PHIR=PARU(2)*RLU(0)
        CPHI=COS(PHIR)
        ANG=CTHE(1)*CTHE(2)-SQRT(1.-CTHE(1)**2)*SQRT(1.-CTHE(2)**2)*CPHI
        Z1=2.-Z(JT)
        Z2=ANG*SQRT(Z(JT)**2-4.*PMQ(JT)**2/SHP)
        Z3=1.-Z(JT)-XH+(PMQ(1)**2+PMQ(2)**2)/SHP
        Z(3-JT)=2./(Z1**2-Z2**2)*(Z1*Z3+Z2*SQRT(Z3**2-(Z1**2-Z2**2)*
     &  PMQ(3-JT)**2/SHP))
        ZMIN=2.*PMQ(3-JT)/SHPR
        ZMAX=1.-PMQ(JT)/SHPR-(SH-PMQ(3-JT)**2)/(SHPR*(SHPR-PMQ(JT)))
        ZMAX=MIN(1.-XH,ZMAX)
        IF(Z(3-JT).LT.ZMIN.OR.Z(3-JT).GT.ZMAX) GOTO 370
        KCC=22
 
      ELSEIF(ISUB.EQ.78) THEN
C...W+/- + H0 -> W+/- + H0.
 
      ELSEIF(ISUB.EQ.79) THEN
C...H0 + H0 -> H0 + H0.
 
      ELSEIF(ISUB.EQ.80) THEN
C...q + gamma -> q' + pi+/-; th=(p(q)-p(q'))**2.
        IF(MINT(15).EQ.22) JS=2
        I=MINT(14+JS)
        IA=IABS(I)
        MINT(23-JS)=ISIGN(211,KCHG(IA,1)*I)
        IB=3-IA
        MINT(20+JS)=ISIGN(IB,I)
        KCC=22
      ENDIF
 
      ELSEIF(ISUB.LE.90) THEN
      IF(ISUB.EQ.81) THEN
C...q + q~ -> Q + Q~; th = (p(q)-p(Q))**2.
        MINT(21)=ISIGN(MINT(55),MINT(15))
        MINT(22)=-MINT(21)
        KCC=4
 
      ELSEIF(ISUB.EQ.82) THEN
C...g + g -> Q + Q~; th arbitrary.
        KCS=(-1)**INT(1.5+RLU(0))
        MINT(21)=ISIGN(MINT(55),KCS)
        MINT(22)=-MINT(21)
        KCC=MINT(2)+10
 
      ELSEIF(ISUB.EQ.83) THEN
C...f + q -> f' + Q; th = (p(f) - p(f'))**2.
        KFOLD=MINT(16)
        IF(MINT(2).EQ.2) KFOLD=MINT(15)
        KFAOLD=IABS(KFOLD)
        IF(KFAOLD.GT.10) THEN
          KFANEW=KFAOLD+2*MOD(KFAOLD,2)-1
        ELSE
          RCKM=VINT(180+KFOLD)*RLU(0)
          IPM=(5-ISIGN(1,KFOLD))/2
          KFANEW=-MOD(KFAOLD+1,2)
  410     KFANEW=KFANEW+2
          IDC=MDCY(KFAOLD,2)+(KFANEW+1)/2+2
          IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.IPM) THEN
            IF(MOD(KFAOLD,2).EQ.0) RCKM=RCKM-VCKM(KFAOLD/2,(KFANEW+1)/2)
            IF(MOD(KFAOLD,2).EQ.1) RCKM=RCKM-VCKM(KFANEW/2,(KFAOLD+1)/2)
          ENDIF
          IF(KFANEW.LE.6.AND.RCKM.GT.0.) GOTO 410
        ENDIF
        IF(MINT(2).EQ.1) THEN
          MINT(21)=ISIGN(MINT(55),MINT(15))
          MINT(22)=ISIGN(KFANEW,MINT(16))
        ELSE
          MINT(21)=ISIGN(KFANEW,MINT(15))
          MINT(22)=ISIGN(MINT(55),MINT(16))
          JS=2
        ENDIF
        KCC=22
 
      ELSEIF(ISUB.EQ.84) THEN
C...g + gamma -> Q + Q~; th arbitary.
        KCS=(-1)**INT(1.5+RLU(0))
        MINT(21)=ISIGN(MINT(55),KCS)
        MINT(22)=-MINT(21)
        KCC=27
        IF(MINT(16).EQ.21) KCC=28
 
      ELSEIF(ISUB.EQ.85) THEN
C...gamma + gamma -> F + F~; th arbitary.
        KCS=(-1)**INT(1.5+RLU(0))
        MINT(21)=ISIGN(MINT(56),KCS)
        MINT(22)=-MINT(21)
        KCC=21
 
      ELSEIF(ISUB.GE.86.AND.ISUB.LE.89) THEN
C...g + g -> (J/Psi, chi_0c, chi_1c or chi_2c) + g
        MINT(21)=KFPR(ISUB,1)
        MINT(22)=KFPR(ISUB,2)
        KCC=24
        KCS=(-1)**INT(1.5+RLU(0))
      ENDIF
 
      ELSEIF(ISUB.LE.100) THEN
      IF(ISUB.EQ.95) THEN
C...Low-pT ( = energyless g + g -> g + g).
        KCC=MINT(2)+12
        KCS=(-1)**INT(1.5+RLU(0))
 
      ELSEIF(ISUB.EQ.96) THEN
C...Multiple interactions (should be reassigned to QCD process).
      ENDIF
 
      ELSEIF(ISUB.LE.110) THEN
      IF(ISUB.EQ.101) THEN
C...g + g -> gamma*/Z0.
        KCC=21
        KFRES=22
 
      ELSEIF(ISUB.EQ.102) THEN
C...g + g -> H0 (or H'0, or A0).
        KCC=21
        KFRES=KFHIGG
 
      ELSEIF(ISUB.EQ.103) THEN
C...gamma + gamma -> H0 (or H'0, or A0).
        KCC=21
        KFRES=KFHIGG
 
      ELSEIF(ISUB.EQ.110) THEN
C...f + f~ -> gamma + H0; th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(20+JS)=22
        MINT(23-JS)=KFHIGG
      ENDIF
 
      ELSEIF(ISUB.LE.120) THEN
      IF(ISUB.EQ.111) THEN
C...f + f~ -> g + H0; th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(20+JS)=21
        MINT(23-JS)=25
        KCC=17+JS
 
      ELSEIF(ISUB.EQ.112) THEN
C...f + g -> f + H0; th = (p(f) - p(f))**2.
        IF(MINT(15).EQ.21) JS=2
        MINT(23-JS)=25
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))
 
      ELSEIF(ISUB.EQ.113) THEN
C...g + g -> g + H0; th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(23-JS)=25
        KCC=22+JS
        KCS=(-1)**INT(1.5+RLU(0))
 
      ELSEIF(ISUB.EQ.114) THEN
C...g + g -> gamma + gamma; th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(21)=22
        MINT(22)=22
        KCC=21
 
      ELSEIF(ISUB.EQ.115) THEN
C...g + g -> g + gamma; th arbitrary.
        IF(RLU(0).GT.0.5) JS=2
        MINT(23-JS)=22
        KCC=22+JS
        KCS=(-1)**INT(1.5+RLU(0))
 
      ELSEIF(ISUB.EQ.116) THEN
C...g + g -> gamma + Z0.
 
      ELSEIF(ISUB.EQ.117) THEN
C...g + g -> Z0 + Z0.
 
      ELSEIF(ISUB.EQ.118) THEN
C...g + g -> W+ + W-.
      ENDIF
 
      ELSEIF(ISUB.LE.140) THEN
      IF(ISUB.EQ.121) THEN
C...g + g -> Q + Q~ + H0.
        KCS=(-1)**INT(1.5+RLU(0))
        MINT(21)=ISIGN(KFPR(ISUBSV,2),KCS)
        MINT(22)=-MINT(21)
        KCC=11+INT(0.5+RLU(0))
        KFRES=KFHIGG
 
      ELSEIF(ISUB.EQ.122) THEN
C...q + q~ -> Q + Q~ + H0.
        MINT(21)=ISIGN(KFPR(ISUBSV,2),MINT(15))
        MINT(22)=-MINT(21)
        KCC=4
        KFRES=KFHIGG
 
      ELSEIF(ISUB.EQ.123) THEN
C...f + f' -> f + f' + H0 (or H'0, or A0) (Z0 + Z0 -> H0 as
C...inner process).
        KCC=22
        KFRES=KFHIGG
 
      ELSEIF(ISUB.EQ.124) THEN
C...f + f' -> f" + f"' + H0 (or H'0, or A) (W+ + W- -> H0 as
C...inner process).
        DO 430 JT=1,2
        I=MINT(14+JT)
        IA=IABS(I)
        IF(IA.LE.10) THEN
          RVCKM=VINT(180+I)*RLU(0)
          DO 420 J=1,MSTP(1)
          IB=2*J-1+MOD(IA,2)
          IPM=(5-ISIGN(1,I))/2
          IDC=J+MDCY(IA,2)+2
          IF(MDME(IDC,1).NE.1.AND.MDME(IDC,1).NE.IPM) GOTO 420
          MINT(20+JT)=ISIGN(IB,I)
          RVCKM=RVCKM-VCKM((IA+1)/2,(IB+1)/2)
          IF(RVCKM.LE.0.) GOTO 430
  420     CONTINUE
        ELSE
          IB=2*((IA+1)/2)-1+MOD(IA,2)
          MINT(20+JT)=ISIGN(IB,I)
        ENDIF
  430   CONTINUE
        KCC=22
        KFRES=KFHIGG
 
      ELSEIF(ISUB.EQ.131) THEN
C...g + g -> Z0 + q + q~.
        MINT(21)=KFPR(131,1)
        MINT(22)=KFPR(131,2)
        MINT(23)=-MINT(22)
        KCC=MINT(2)+10
        KCS=1
      ENDIF
 
      ELSEIF(ISUB.LE.160) THEN
      IF(ISUB.EQ.141) THEN
C...f + f~ -> gamma*/Z0/Z'0.
        KFRES=32
 
      ELSEIF(ISUB.EQ.142) THEN
C...f + f~' -> W'+/- .
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        KFRES=ISIGN(34,KCH1+KCH2)
 
      ELSEIF(ISUB.EQ.143) THEN
C...f + f~' -> H+/-.
        KCH1=KCHG(IABS(MINT(15)),1)*ISIGN(1,MINT(15))
        KCH2=KCHG(IABS(MINT(16)),1)*ISIGN(1,MINT(16))
        KFRES=ISIGN(37,KCH1+KCH2)
 
      ELSEIF(ISUB.EQ.144) THEN
C...f + f~' -> R.
        KFRES=ISIGN(40,MINT(15)+MINT(16))
 
      ELSEIF(ISUB.EQ.145) THEN
C...q + l -> LQ (leptoquark).
        IF(IABS(MINT(16)).LE.8) JS=2
        KFRES=ISIGN(39,MINT(14+JS))
        KCC=28+JS
        KCS=ISIGN(1,MINT(14+JS))
 
      ELSEIF(ISUB.EQ.147.OR.ISUB.EQ.148) THEN
C...q + g -> q* (excited quark).
        IF(MINT(15).EQ.21) JS=2
        KFRES=MINT(14+JS)+ISIGN(6,MINT(14+JS))
        KCC=30+JS
        KCS=ISIGN(1,MINT(14+JS))
 
      ELSEIF(ISUB.EQ.149) THEN
C...g + g -> eta_techni.
        KFRES=38
        KCC=23
        KCS=(-1)**INT(1.5+RLU(0))
      ENDIF
 
      ELSE
      IF(ISUB.EQ.161) THEN
C...f + g -> f' + H+/-; th = (p(f)-p(f'))**2.
        IF(MINT(15).EQ.21) JS=2
        I=MINT(14+JS)
        IA=IABS(I)
        MINT(23-JS)=ISIGN(37,KCHG(IA,1)*I)
        IB=IA+MOD(IA,2)-MOD(IA+1,2)
        MINT(20+JS)=ISIGN(IB,I)
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))
 
      ELSEIF(ISUB.EQ.162) THEN
C...q + g -> LQ + l~; LQ=leptoquark; th=(p(q)-p(LQ))^2.
        IF(MINT(15).EQ.21) JS=2
        MINT(20+JS)=ISIGN(39,MINT(14+JS))
        KFLQL=KFDP(MDCY(39,2),2)
        MINT(23-JS)=-ISIGN(KFLQL,MINT(14+JS))
        KCC=15+JS
        KCS=ISIGN(1,MINT(14+JS))
 
      ELSEIF(ISUB.EQ.163) THEN
C...g + g -> LQ + LQ~; LQ=leptoquark; th arbitrary.
        KCS=(-1)**INT(1.5+RLU(0))
        MINT(21)=ISIGN(39,KCS)
        MINT(22)=-MINT(21)
        KCC=MINT(2)+10
 
      ELSEIF(ISUB.EQ.164) THEN
C...q + q~ -> LQ + LQ~; LQ=leptoquark; th=(p(q)-p(LQ))**2.
        MINT(21)=ISIGN(39,MINT(15))
        MINT(22)=-MINT(21)
        KCC=4
 
      ELSEIF(ISUB.EQ.165) THEN
C...q + q~ -> l- + l+; th=(p(q)-p(l-))**2.
        MINT(21)=ISIGN(KFPR(ISUB,1),MINT(15))
        MINT(22)=-MINT(21)
 
      ELSEIF(ISUB.EQ.166) THEN
C...q + q~' -> l + nu; th=(p(u)-p(nu))**2 or (p(u~)-p(nu~))**2.
        IF(MOD(MINT(15),2).EQ.0) THEN
          MINT(21)=ISIGN(KFPR(ISUB,1)+1,MINT(15))
          MINT(22)=ISIGN(KFPR(ISUB,1),MINT(16))
        ELSE
          MINT(21)=ISIGN(KFPR(ISUB,1),MINT(15))
          MINT(22)=ISIGN(KFPR(ISUB,1)+1,MINT(16))
        ENDIF
 
      ELSEIF(ISUB.EQ.167.OR.ISUB.EQ.168) THEN
C...q + q' -> q" + q* (excited quark).
        KFQEXC=ISUB-166
        KFQSTR=ISUB-160
        JS=MINT(2)
        MINT(20+JS)=ISIGN(KFQSTR,MINT(14+JS))
        IF(IABS(MINT(15)).NE.KFQEXC.AND.IABS(MINT(16)).NE.KFQEXC)
     &  MINT(23-JS)=ISIGN(KFQEXC,MINT(17-JS))
        KCC=22
      ENDIF
      ENDIF
 
      IF(ISET(ISUB).EQ.11) THEN
C...Store documentation for user-defined processes.
        BEZUP=(PUP(1,4)-PUP(2,4))/(PUP(1,4)+PUP(2,4))
        KUPPO(1)=MINT(83)+5
        KUPPO(2)=MINT(83)+6
        I=MINT(83)+6
        DO 450 IUP=3,NUP
        KUPPO(IUP)=0
        IF(MSTP(128).GE.2.AND.KUP(IUP,3).NE.0) THEN
          IDOC=IDOC-1
          MINT(4)=MINT(4)-1
          GOTO 450
        ENDIF
        I=I+1
        KUPPO(IUP)=I
        K(I,1)=21
        K(I,2)=KUP(IUP,2)
        K(I,3)=0
        IF(KUP(IUP,3).NE.0) K(I,3)=KUPPO(KUP(IUP,3))
        K(I,4)=0
        K(I,5)=0
        DO 440 J=1,5
        P(I,J)=PUP(IUP,J)
  440   CONTINUE
  450   CONTINUE
        CALL LUDBRB(MINT(83)+7,MINT(83)+4+NUP,0.,VINT(24),0D0,0D0,
     &  -DBLE(BEZUP))
 
C...Store final state partons for user-defined processes.
        N=IPU2
        DO 470 IUP=3,NUP
        N=N+1
        K(N,1)=1
        IF(KUP(IUP,1).NE.1) K(N,1)=11
        K(N,2)=KUP(IUP,2)
        IF(MSTP(128).LE.0.OR.KUP(IUP,3).EQ.0) THEN
          K(N,3)=KUPPO(IUP)
        ELSE
          K(N,3)=MINT(84)+KUP(IUP,3)
        ENDIF
        K(N,4)=0
        K(N,5)=0
        DO 460 J=1,5
        P(N,J)=PUP(IUP,J)
  460   CONTINUE
  470   CONTINUE
        CALL LUDBRB(IPU3,N,0.,VINT(24),0D0,0D0,-DBLE(BEZUP))
 
C...Arrange colour flow for user-defined processes.
        N=MINT(84)
        DO 480 IUP=1,NUP
        N=N+1
        IF(KCHG(LUCOMP(K(N,2)),2).EQ.0) GOTO 480
        IF(K(N,1).EQ.1) K(N,1)=3
        IF(K(N,1).EQ.11) K(N,1)=14
        IF(KUP(IUP,4).NE.0) K(N,4)=K(N,4)+MSTU(5)*(KUP(IUP,4)+MINT(84))
        IF(KUP(IUP,5).NE.0) K(N,5)=K(N,5)+MSTU(5)*(KUP(IUP,5)+MINT(84))
        IF(KUP(IUP,6).NE.0) K(N,4)=K(N,4)+KUP(IUP,6)+MINT(84)
        IF(KUP(IUP,7).NE.0) K(N,5)=K(N,5)+KUP(IUP,7)+MINT(84)
  480   CONTINUE
 
      ELSEIF(IDOC.EQ.7) THEN
C...Resonance not decaying; store kinematics.
        I=MINT(83)+7
        K(IPU3,1)=1
        K(IPU3,2)=KFRES
        K(IPU3,3)=I
        P(IPU3,4)=SHUSER
        P(IPU3,5)=SHUSER
        K(I,1)=21
        K(I,2)=KFRES
        P(I,4)=SHUSER
        P(I,5)=SHUSER
        N=IPU3
        MINT(21)=KFRES
        MINT(22)=0
 
C...Special cases: colour flow in g + g -> eta_techni, q + l -> LQ
C...and q + g -> q*.
        IF(KFRES.EQ.38.OR.IABS(KFRES).EQ.39.OR.(MSTP(6).EQ.1.AND.
     &  (IABS(KFRES).EQ.7.OR.IABS(KFRES).EQ.8))) THEN
          K(IPU3,1)=3
          DO 490 J=1,2
          JC=J
          IF(KCS.EQ.-1) JC=3-J
          IF(ICOL(KCC,1,JC).NE.0.AND.K(IPU1,1).EQ.14) K(IPU1,J+3)=
     &    MINT(84)+ICOL(KCC,1,JC)
          IF(ICOL(KCC,2,JC).NE.0.AND.K(IPU2,1).EQ.14) K(IPU2,J+3)=
     &    MINT(84)+ICOL(KCC,2,JC)
          IF(ICOL(KCC,3,JC).NE.0.AND.K(IPU3,1).EQ.3) K(IPU3,J+3)=
     &    MSTU(5)*(MINT(84)+ICOL(KCC,3,JC))
  490     CONTINUE
        ELSE
          K(IPU1,4)=IPU2
          K(IPU1,5)=IPU2
          K(IPU2,4)=IPU1
          K(IPU2,5)=IPU1
        ENDIF
 
      ELSEIF(IDOC.EQ.8) THEN
C...2 -> 2 processes: store outgoing partons in their CM-frame.
        DO 500 JT=1,2
        I=MINT(84)+2+JT
        K(I,1)=1
        IF(IABS(MINT(20+JT)).LE.100) THEN
          IF(KCHG(IABS(MINT(20+JT)),2).NE.0) K(I,1)=3
        ENDIF
        K(I,2)=MINT(20+JT)
        K(I,3)=MINT(83)+IDOC+JT-2
        KFAA=IABS(K(I,2))
        IF(KFAA.GE.23.OR.(KFAA.EQ.6.AND.KFPR(ISUBSV,1).NE.0.AND.
     &  MSTP(48).GE.1).OR.((KFAA.EQ.7.OR.KFAA.EQ.8.OR.KFAA.EQ.17.OR.
     &  KFAA.EQ.18).AND.KFPR(ISUBSV,1).NE.0.AND.MSTP(49).GE.1)) THEN
          P(I,5)=SQRT(VINT(63+MOD(JS+JT,2)))
        ELSEIF((KFAA.EQ.7.OR.KFAA.EQ.8).AND.MSTP(6).EQ.1.AND.
     &  KFPR(ISUBSV,2).NE.0) THEN
          P(I,5)=SQRT(VINT(64))
        ELSE
          P(I,5)=ULMASS(K(I,2))
        ENDIF
  500   CONTINUE
        IF(P(IPU3,5)+P(IPU4,5).GE.SHR) THEN
          KFA1=IABS(MINT(21))
          KFA2=IABS(MINT(22))
          IF((KFA1.GT.3.AND.KFA1.NE.21).OR.(KFA2.GT.3.AND.KFA2.NE.21))
     &    THEN
            MINT(51)=1
            RETURN
          ENDIF
          P(IPU3,5)=0.
          P(IPU4,5)=0.
        ENDIF
        P(IPU3,4)=0.5*(SHR+(P(IPU3,5)**2-P(IPU4,5)**2)/SHR)
        P(IPU3,3)=SQRT(MAX(0.,P(IPU3,4)**2-P(IPU3,5)**2))
        P(IPU4,4)=SHR-P(IPU3,4)
        P(IPU4,3)=-P(IPU3,3)
        N=IPU4
        MINT(7)=MINT(83)+7
        MINT(8)=MINT(83)+8
 
C...Rotate outgoing partons using cos(theta)=(th-uh)/lam(sh,sqm3,sqm4).
        CALL LUDBRB(IPU3,IPU4,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0)
 
      ELSEIF(IDOC.EQ.9.AND.ISET(ISUB).EQ.5) THEN
C...2 -> 3 processes (alt 1): store outgoing partons in their CM frame.
        DO 510 JT=1,2
        I=MINT(84)+2+JT
        K(I,1)=1
        IF(IABS(MINT(20+JT)).LE.100) THEN
          IF(KCHG(IABS(MINT(20+JT)),2).NE.0) K(I,1)=3
        ENDIF
        K(I,2)=MINT(20+JT)
        K(I,3)=MINT(83)+IDOC+JT-3
        IF(IABS(K(I,2)).LE.22) THEN
          P(I,5)=ULMASS(K(I,2))
        ELSE
          P(I,5)=SQRT(VINT(63+MOD(JS+JT,2)))
        ENDIF
        PT=SQRT(MAX(0.,VINT(197+5*JT)-P(I,5)**2+VINT(196+5*JT)**2))
        P(I,1)=PT*COS(VINT(198+5*JT))
        P(I,2)=PT*SIN(VINT(198+5*JT))
  510   CONTINUE
        K(IPU5,1)=1
        K(IPU5,2)=KFRES
        K(IPU5,3)=MINT(83)+IDOC
        P(IPU5,5)=SHR
        P(IPU5,1)=-P(IPU3,1)-P(IPU4,1)
        P(IPU5,2)=-P(IPU3,2)-P(IPU4,2)
        PMS1=P(IPU3,5)**2+P(IPU3,1)**2+P(IPU3,2)**2
        PMS2=P(IPU4,5)**2+P(IPU4,1)**2+P(IPU4,2)**2
        PMS3=P(IPU5,5)**2+P(IPU5,1)**2+P(IPU5,2)**2
        PMT3=SQRT(PMS3)
        P(IPU5,3)=PMT3*SINH(VINT(211))
        P(IPU5,4)=PMT3*COSH(VINT(211))
        PMS12=(SHPR-P(IPU5,4))**2-P(IPU5,3)**2
        SQL12=(PMS12-PMS1-PMS2)**2-4.*PMS1*PMS2
        IF(SQL12.LE.0.) THEN
          MINT(51)=1
          RETURN
        ENDIF
        P(IPU3,3)=(-P(IPU5,3)*(PMS12+PMS1-PMS2)+
     &  VINT(213)*(SHPR-P(IPU5,4))*SQRT(SQL12))/(2.*PMS12)
        P(IPU4,3)=-P(IPU3,3)-P(IPU5,3)
        P(IPU3,4)=SQRT(PMS1+P(IPU3,3)**2)
        P(IPU4,4)=SQRT(PMS2+P(IPU4,3)**2)
        MINT(23)=KFRES
        N=IPU5
        MINT(7)=MINT(83)+7
        MINT(8)=MINT(83)+8
 
      ELSEIF(IDOC.EQ.9) THEN
C...2 -> 3 processes: store outgoing partons in their CM frame.
        DO 520 JT=1,3
        I=MINT(84)+2+JT
        K(I,1)=1
        IF(IABS(MINT(20+JT)).LE.10.OR.MINT(20+JT).EQ.21) K(I,1)=3
        K(I,2)=MINT(20+JT)
        K(I,3)=MINT(83)+IDOC+JT-3
        IF(JT.EQ.1) THEN
          P(I,5)=SQRT(VINT(63))
        ELSE
          P(I,5)=PMAS(KFPR(ISUB,2),1)
        ENDIF
  520   CONTINUE
        P(IPU3,4)=0.5*(SHR+(VINT(63)-VINT(64))/SHR)
        P(IPU3,3)=SQRT(MAX(0.,P(IPU3,4)**2-P(IPU3,5)**2))
        P(IPU4,4)=0.5*SQRT(VINT(64))
        P(IPU4,3)=SQRT(MAX(0.,P(IPU4,4)**2-P(IPU4,5)**2))
        P(IPU5,4)=P(IPU4,4)
        P(IPU5,3)=-P(IPU4,3)
        N=IPU5
        MINT(7)=MINT(83)+7
        MINT(8)=MINT(83)+9
 
C...Rotate and boost outgoing partons.
        CALL LUDBRB(IPU4,IPU5,ACOS(VINT(83)),VINT(84),0D0,0D0,0D0)
        CALL LUDBRB(IPU4,IPU5,0.,0.,0D0,0D0,
     &  -DBLE(P(IPU3,3)/(SHR-P(IPU3,4))))
        CALL LUDBRB(IPU3,IPU5,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0)
 
      ELSEIF(IDOC.EQ.11) THEN
C...Z0 + Z0 -> H0, W+ + W- -> H0: store Higgs and outgoing partons.
        PHI(1)=PARU(2)*RLU(0)
        PHI(2)=PHI(1)-PHIR
        DO 530 JT=1,2
        I=MINT(84)+2+JT
        K(I,1)=1
        IF(IABS(MINT(20+JT)).LE.10.OR.MINT(20+JT).EQ.21) K(I,1)=3
        K(I,2)=MINT(20+JT)
        K(I,3)=MINT(83)+IDOC+JT-2
        P(I,5)=ULMASS(K(I,2))
        IF(0.5*SHPR*Z(JT).LE.P(I,5)) P(I,5)=0.
        PABS=SQRT(MAX(0.,(0.5*SHPR*Z(JT))**2-P(I,5)**2))
        PTABS=PABS*SQRT(MAX(0.,1.-CTHE(JT)**2))
        P(I,1)=PTABS*COS(PHI(JT))
        P(I,2)=PTABS*SIN(PHI(JT))
        P(I,3)=PABS*CTHE(JT)*(-1)**(JT+1)
        P(I,4)=0.5*SHPR*Z(JT)
        IZW=MINT(83)+6+JT
        K(IZW,1)=21
        K(IZW,2)=23
        IF(ISUB.EQ.8) K(IZW,2)=ISIGN(24,LUCHGE(MINT(14+JT)))
        K(IZW,3)=IZW-2
        P(IZW,1)=-P(I,1)
        P(IZW,2)=-P(I,2)
        P(IZW,3)=(0.5*SHPR-PABS*CTHE(JT))*(-1)**(JT+1)
        P(IZW,4)=0.5*SHPR*(1.-Z(JT))
        P(IZW,5)=-SQRT(MAX(0.,P(IZW,3)**2+PTABS**2-P(IZW,4)**2))
  530   CONTINUE
        I=MINT(83)+9
        K(IPU5,1)=1
        K(IPU5,2)=KFRES
        K(IPU5,3)=I
        P(IPU5,5)=SHR
        P(IPU5,1)=-P(IPU3,1)-P(IPU4,1)
        P(IPU5,2)=-P(IPU3,2)-P(IPU4,2)
        P(IPU5,3)=-P(IPU3,3)-P(IPU4,3)
        P(IPU5,4)=SHPR-P(IPU3,4)-P(IPU4,4)
        K(I,1)=21
        K(I,2)=KFRES
        DO 540 J=1,5
        P(I,J)=P(IPU5,J)
  540   CONTINUE
        N=IPU5
        MINT(23)=KFRES
 
      ELSEIF(IDOC.EQ.12) THEN
C...Z0 and W+/- scattering: store bosons and outgoing partons.
        PHI(1)=PARU(2)*RLU(0)
        PHI(2)=PHI(1)-PHIR
        JTRAN=INT(1.5+RLU(0))
        DO 550 JT=1,2
        I=MINT(84)+2+JT
        K(I,1)=1
        IF(IABS(MINT(20+JT)).LE.10.OR.MINT(20+JT).EQ.21) K(I,1)=3
        K(I,2)=MINT(20+JT)
        K(I,3)=MINT(83)+IDOC+JT-2
        P(I,5)=ULMASS(K(I,2))
        IF(0.5*SHPR*Z(JT).LE.P(I,5)) P(I,5)=0.
        PABS=SQRT(MAX(0.,(0.5*SHPR*Z(JT))**2-P(I,5)**2))
        PTABS=PABS*SQRT(MAX(0.,1.-CTHE(JT)**2))
        P(I,1)=PTABS*COS(PHI(JT))
        P(I,2)=PTABS*SIN(PHI(JT))
        P(I,3)=PABS*CTHE(JT)*(-1)**(JT+1)
        P(I,4)=0.5*SHPR*Z(JT)
        IZW=MINT(83)+6+JT
        K(IZW,1)=21
        IF(MINT(14+JT).EQ.MINT(20+JT)) THEN
          K(IZW,2)=23
        ELSE
          K(IZW,2)=ISIGN(24,LUCHGE(MINT(14+JT))-LUCHGE(MINT(20+JT)))
        ENDIF
        K(IZW,3)=IZW-2
        P(IZW,1)=-P(I,1)
        P(IZW,2)=-P(I,2)
        P(IZW,3)=(0.5*SHPR-PABS*CTHE(JT))*(-1)**(JT+1)
        P(IZW,4)=0.5*SHPR*(1.-Z(JT))
        P(IZW,5)=-SQRT(MAX(0.,P(IZW,3)**2+PTABS**2-P(IZW,4)**2))
        IPU=MINT(84)+4+JT
        K(IPU,1)=3
        K(IPU,2)=KFPR(ISUB,JT)
        IF(ISUB.EQ.72.AND.JT.EQ.JTRAN) K(IPU,2)=-K(IPU,2)
        IF(ISUB.EQ.73.OR.ISUB.EQ.77) K(IPU,2)=K(IZW,2)
        K(IPU,3)=MINT(83)+8+JT
        IF(IABS(K(IPU,2)).LE.10.OR.K(IPU,2).EQ.21) THEN
          P(IPU,5)=ULMASS(K(IPU,2))
        ELSE
          P(IPU,5)=SQRT(VINT(63+MOD(JS+JT,2)))
        ENDIF
        MINT(22+JT)=K(IPU,2)
  550   CONTINUE
C...Find rotation and boost for hard scattering subsystem.
        I1=MINT(83)+7
        I2=MINT(83)+8
        BEXCM=(P(I1,1)+P(I2,1))/(P(I1,4)+P(I2,4))
        BEYCM=(P(I1,2)+P(I2,2))/(P(I1,4)+P(I2,4))
        BEZCM=(P(I1,3)+P(I2,3))/(P(I1,4)+P(I2,4))
        GAMCM=(P(I1,4)+P(I2,4))/SHR
        BEPCM=BEXCM*P(I1,1)+BEYCM*P(I1,2)+BEZCM*P(I1,3)
        PX=P(I1,1)+GAMCM*(GAMCM/(1.+GAMCM)*BEPCM-P(I1,4))*BEXCM
        PY=P(I1,2)+GAMCM*(GAMCM/(1.+GAMCM)*BEPCM-P(I1,4))*BEYCM
        PZ=P(I1,3)+GAMCM*(GAMCM/(1.+GAMCM)*BEPCM-P(I1,4))*BEZCM
        THECM=ULANGL(PZ,SQRT(PX**2+PY**2))
        PHICM=ULANGL(PX,PY)
C...Store hard scattering subsystem. Rotate and boost it.
        SQLAM=(SH-P(IPU5,5)**2-P(IPU6,5)**2)**2-4.*P(IPU5,5)**2*
     &  P(IPU6,5)**2
        PABS=SQRT(MAX(0.,SQLAM/(4.*SH)))
        CTHWZ=VINT(23)
        STHWZ=SQRT(MAX(0.,1.-CTHWZ**2))
        PHIWZ=VINT(24)-PHICM
        P(IPU5,1)=PABS*STHWZ*COS(PHIWZ)
        P(IPU5,2)=PABS*STHWZ*SIN(PHIWZ)
        P(IPU5,3)=PABS*CTHWZ
        P(IPU5,4)=SQRT(PABS**2+P(IPU5,5)**2)
        P(IPU6,1)=-P(IPU5,1)
        P(IPU6,2)=-P(IPU5,2)
        P(IPU6,3)=-P(IPU5,3)
        P(IPU6,4)=SQRT(PABS**2+P(IPU6,5)**2)
        CALL LUDBRB(IPU5,IPU6,THECM,PHICM,DBLE(BEXCM),DBLE(BEYCM),
     &  DBLE(BEZCM))
        DO 570 JT=1,2
        I1=MINT(83)+8+JT
        I2=MINT(84)+4+JT
        K(I1,1)=21
        K(I1,2)=K(I2,2)
        DO 560 J=1,5
        P(I1,J)=P(I2,J)
  560   CONTINUE
  570   CONTINUE
        N=IPU6
        MINT(7)=MINT(83)+9
        MINT(8)=MINT(83)+10
      ENDIF
 
      IF(ISET(ISUB).EQ.11) THEN
      ELSEIF(IDOC.GE.8.AND.ISET(ISUB).NE.6) THEN
C...Store colour connection indices.
        DO 580 J=1,2
        JC=J
        IF(KCS.EQ.-1) JC=3-J
        IF(ICOL(KCC,1,JC).NE.0.AND.K(IPU1,1).EQ.14) K(IPU1,J+3)=
     &  K(IPU1,J+3)+MINT(84)+ICOL(KCC,1,JC)
        IF(ICOL(KCC,2,JC).NE.0.AND.K(IPU2,1).EQ.14) K(IPU2,J+3)=
     &  K(IPU2,J+3)+MINT(84)+ICOL(KCC,2,JC)
        IF(ICOL(KCC,3,JC).NE.0.AND.K(IPU3,1).EQ.3) K(IPU3,J+3)=
     &  MSTU(5)*(MINT(84)+ICOL(KCC,3,JC))
        IF(ICOL(KCC,4,JC).NE.0.AND.K(IPU4,1).EQ.3) K(IPU4,J+3)=
     &  MSTU(5)*(MINT(84)+ICOL(KCC,4,JC))
  580   CONTINUE
 
C...Copy outgoing partons to documentation lines.
        IMAX=2
        IF(IDOC.EQ.9) IMAX=3
        DO 600 I=1,IMAX
        I1=MINT(83)+IDOC-IMAX+I
        I2=MINT(84)+2+I
        K(I1,1)=21
        K(I1,2)=K(I2,2)
        IF(IDOC.LE.9) K(I1,3)=0
        IF(IDOC.GE.11) K(I1,3)=MINT(83)+2+I
        DO 590 J=1,5
        P(I1,J)=P(I2,J)
  590   CONTINUE
  600   CONTINUE
 
      ELSEIF(IDOC.EQ.9) THEN
C...Store colour connection indices.
        DO 610 J=1,2
        JC=J
        IF(KCS.EQ.-1) JC=3-J
        IF(ICOL(KCC,1,JC).NE.0.AND.K(IPU1,1).EQ.14) K(IPU1,J+3)=
     &  K(IPU1,J+3)+MINT(84)+ICOL(KCC,1,JC)+
     &  MAX(0,MIN(1,ICOL(KCC,1,JC)-2))
        IF(ICOL(KCC,2,JC).NE.0.AND.K(IPU2,1).EQ.14) K(IPU2,J+3)=
     &  K(IPU2,J+3)+MINT(84)+ICOL(KCC,2,JC)+
     &  MAX(0,MIN(1,ICOL(KCC,2,JC)-2))
        IF(ICOL(KCC,3,JC).NE.0.AND.K(IPU4,1).EQ.3) K(IPU4,J+3)=
     &  MSTU(5)*(MINT(84)+ICOL(KCC,3,JC))
        IF(ICOL(KCC,4,JC).NE.0.AND.K(IPU5,1).EQ.3) K(IPU5,J+3)=
     &  MSTU(5)*(MINT(84)+ICOL(KCC,4,JC))
  610   CONTINUE
 
C...Copy outgoing partons to documentation lines.
        DO 630 I=1,3
        I1=MINT(83)+IDOC-3+I
        I2=MINT(84)+2+I
        K(I1,1)=21
        K(I1,2)=K(I2,2)
        K(I1,3)=0
        DO 620 J=1,5
        P(I1,J)=P(I2,J)
  620   CONTINUE
  630   CONTINUE
      ENDIF
 
C...Low-pT events: remove gluons used for string drawing purposes.
      IF(ISUB.EQ.95) THEN
        K(IPU3,1)=K(IPU3,1)+10
        K(IPU4,1)=K(IPU4,1)+10
        DO 640 J=41,66
        VINTSV(J)=VINT(J)
        VINT(J)=0.
  640   CONTINUE
        DO 660 I=MINT(83)+5,MINT(83)+8
        DO 650 J=1,5
        P(I,J)=0.
  650   CONTINUE
  660   CONTINUE
      ENDIF
 
      RETURN
      END
