 
C***********************************************************************
 
      SUBROUTINE PYSIGH(NCHN,SIGS)
 
C...Differential matrix elements for all included subprocesses.
C...Note that what is coded is (disregarding the COMFAC factor)
C...1) for 2 -> 1 processes: s-hat/pi*d(sigma-hat), where,
C...when d(sigma-hat) is given in the zero-width limit, the delta
C...function in tau is replaced by a (modified) Breit-Wigner:
C...1/pi*s*H_res/((s*tau-m_res^2)^2+H_res^2),
C...where H_res = s-hat/m_res*Gamma_res(s-hat);
C...2) for 2 -> 2 processes: (s-hat)**2/pi*d(sigma-hat)/d(t-hat);
C...i.e., dimensionless quantities.
C...3) for 2 -> 3 processes: abs(M)^2, where the total cross-section is
C...Integral abs(M)^2/(2shat') * (prod_(i=1)^3 d^3p_i/((2pi)^3*2E_i)) *
C...(2pi)^4 delta^4(P - sum p_i).
C...COMFAC contains the factor pi/s (or equivalent) and
C...the conversion factor from GeV^-2 to mb.
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
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT3/,/PYINT4/,
     &/PYINT5/,/PYINT7/
      DIMENSION X(2),XPQ(-25:25),KFAC(2,-40:40),WDTP(0:40),
     &WDTE(0:40,0:5),HGZ(6,3),HL3(3),HR3(3),HL4(3),HR4(3)
      COMPLEX A004,A204,A114,A00U,A20U,A11U
      COMPLEX CIGTOT,CIZTOT,F0ALP,F1ALP,F2ALP,F0BET,F1BET,F2BET,FIF,
     &COULCK,COULCP,COULCD,COULCR,COULCS
 
C...The following gives an interface for process 131, gg -> Zqq,
C...to the matrix element package of Ronald Kleiss.
      COMMON/RKBBVC/RKMQ,RKMZ,RKGZ,RKVQ,RKAQ,RKVL,RKAL
      SAVE /RKBBVC/
      DIMENSION RKG1(0:3),RKG2(0:3),RKQ1(0:3),RKQ2(0:3),RKL1(0:3),
     &RKL2(0:3)
 
C...Reset number of channels and cross-section.
      NCHN=0
      SIGS=0.
 
C...Convert H' or A process into equivalent H one.
      ISUB=MINT(1)
      ISUBSV=ISUB
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
 
C...Read kinematical variables and limits.
      ISTSB=ISET(ISUBSV)
      TAUMIN=VINT(11)
      YSTMIN=VINT(12)
      CTNMIN=VINT(13)
      CTPMIN=VINT(14)
      TAUPMN=VINT(16)
      TAU=VINT(21)
      YST=VINT(22)
      CTH=VINT(23)
      XT2=VINT(25)
      TAUP=VINT(26)
      TAUMAX=VINT(31)
      YSTMAX=VINT(32)
      CTNMAX=VINT(33)
      CTPMAX=VINT(34)
      TAUPMX=VINT(36)
 
C...Derive kinematical quantities.
      TAUE=TAU
      IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=TAUP
      X(1)=SQRT(TAUE)*EXP(YST)
      X(2)=SQRT(TAUE)*EXP(-YST)
      IF(MINT(45).EQ.2.AND.ISTSB.GE.1) THEN
        IF(X(1).GT.0.9999) RETURN
      ELSEIF(MINT(45).EQ.3) THEN
        X(1)=MIN(0.9999989,X(1))
      ENDIF
      IF(MINT(46).EQ.2.AND.ISTSB.GE.1) THEN
        IF(X(2).GT.0.9999) RETURN
      ELSEIF(MINT(46).EQ.3) THEN
        X(2)=MIN(0.9999989,X(2))
      ENDIF
      SH=TAU*VINT(2)
      SQM3=VINT(63)
      SQM4=VINT(64)
      RM3=SQM3/SH
      RM4=SQM4/SH
      BE34=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4))
      RPTS=4.*VINT(71)**2/SH
      BE34L=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4-RPTS))
      RM34=MAX(1E-20,2.*RM3*RM4)
      RSQM=1.+RM34
      IF(2.*VINT(71)**2/(VINT(21)*VINT(2)).LT.0.0001) RM34=MAX(RM34,
     &2.*VINT(71)**2/(VINT(21)*VINT(2)))
      RTHM=(4.*RM3*RM4+RPTS)/(1.-RM3-RM4+BE34L)
      IF(ISTSB.EQ.0) THEN
        TH=VINT(45)
        UH=-0.5*SH*MAX(RTHM,1.-RM3-RM4+BE34*CTH)
        SQPTH=MAX(VINT(71)**2,0.25*SH*BE34**2*VINT(59)**2)
      ELSE
        TH=-0.5*SH*MAX(RTHM,1.-RM3-RM4-BE34*CTH)
        UH=-0.5*SH*MAX(RTHM,1.-RM3-RM4+BE34*CTH)
        SQPTH=MAX(VINT(71)**2,0.25*SH*BE34**2*(1.-CTH**2))
      ENDIF
      SH2=SH**2
      TH2=TH**2
      UH2=UH**2
 
C...Choice of Q2 scale: hard, structure functions, parton showers.
      IF(ISTSB.EQ.1.OR.ISTSB.EQ.3.OR.ISTSB.EQ.5) THEN
        Q2=SH
      ELSEIF(MOD(ISTSB,2).EQ.0.OR.ISTSB.EQ.9) THEN
        IF(MSTP(32).EQ.1) THEN
          Q2=2.*SH*TH*UH/(SH**2+TH**2+UH**2)
        ELSEIF(MSTP(32).EQ.2) THEN
          Q2=SQPTH+0.5*(SQM3+SQM4)
        ELSEIF(MSTP(32).EQ.3) THEN
          Q2=MIN(-TH,-UH)
        ELSEIF(MSTP(32).EQ.4) THEN
          Q2=SH
        ELSEIF(MSTP(32).EQ.5) THEN
          Q2=-TH
        ENDIF
        IF(ISTSB.EQ.9) Q2=SQPTH
        IF((ISTSB.EQ.9.AND.MSTP(82).GE.2).OR.(ISTSB.NE.9.AND.
     &  MSTP(85).EQ.1)) Q2=Q2+PARP(82)**2
      ENDIF
      Q2SF=Q2
      IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
        Q2SF=PMAS(23,1)**2
        IF(ISUB.EQ.8.OR.ISUB.EQ.76.OR.ISUB.EQ.77.OR.ISUB.EQ.124)
     &  Q2SF=PMAS(24,1)**2
        IF(ISUB.EQ.121.OR.ISUB.EQ.122) THEN
          Q2SF=PMAS(KFPR(ISUBSV,2),1)**2
          IF(MSTP(39).EQ.2) Q2SF=Q2SF+MAX(VINT(202),VINT(207))
          IF(MSTP(39).EQ.3) Q2SF=SH
          IF(MSTP(39).EQ.4) Q2SF=VINT(26)*VINT(2)
        ENDIF
      ENDIF
      Q2PS=Q2SF
      Q2SF=Q2SF*PARP(34)
      IF(MSTP(22).GE.1.AND.(ISUB.EQ.10.OR.ISUB.EQ.83).AND.
     &(MINT(43).EQ.2.OR.MINT(43).EQ.3)) THEN
        XBJ=X(2)
        IF(MINT(43).EQ.3) XBJ=X(1)
        IF(MSTP(22).EQ.1) THEN
          Q2PS=-TH
        ELSEIF(MSTP(22).EQ.2) THEN
          Q2PS=((1.-XBJ)/XBJ)*(-TH)
        ELSEIF(MSTP(22).EQ.3) THEN
          Q2PS=SQRT((1.-XBJ)/XBJ)*(-TH)
        ELSE
          Q2PS=(1.-XBJ)*MAX(1.,-LOG(XBJ))*(-TH)
        ENDIF
      ENDIF
 
C...Store derived kinematical quantities.
      VINT(41)=X(1)
      VINT(42)=X(2)
      VINT(44)=SH
      VINT(43)=SQRT(SH)
      VINT(45)=TH
      VINT(46)=UH
      VINT(48)=SQPTH
      VINT(47)=SQRT(SQPTH)
      VINT(50)=TAUP*VINT(2)
      VINT(49)=SQRT(MAX(0.,VINT(50)))
      VINT(52)=Q2
      VINT(51)=SQRT(Q2)
      VINT(54)=Q2SF
      VINT(53)=SQRT(Q2SF)
      VINT(56)=Q2PS
      VINT(55)=SQRT(Q2PS)
 
C...Calculate parton structure functions.
      IF(ISTSB.LE.0) GOTO 160
      IF(MINT(47).GE.2) THEN
        DO 110 I=3-MIN(2,MINT(45)),MIN(2,MINT(46))
        XSF=X(I)
        IF(ISTSB.EQ.9) XSF=X(I)/VINT(142+I)
        MINT(105)=MINT(102+I)
        MINT(109)=MINT(106+I)
        IF(MSTP(57).LE.1) THEN
          CALL PYSTFU(MINT(10+I),XSF,Q2SF,XPQ)
        ELSE
          CALL PYSTFL(MINT(10+I),XSF,Q2SF,XPQ)
        ENDIF
        DO 100 KFL=-25,25
        XSFX(I,KFL)=XPQ(KFL)
  100   CONTINUE
  110   CONTINUE
      ENDIF
 
C...Calculate alpha_em, alpha_strong and K-factor.
      XW=PARU(102)
      XWV=XW
      IF(MSTP(8).GE.2.OR.(ISUB.GE.71.AND.ISUB.LE.77)) XW=
     &1.-(PMAS(24,1)/PMAS(23,1))**2
      XW1=1.-XW
      XWC=1./(16.*XW*XW1)
      AEM=ULALEM(Q2)
      IF(MSTP(8).GE.1) AEM=SQRT(2.)*PARU(105)*PMAS(24,1)**2*XW/PARU(1)
      IF(MSTP(33).NE.3) AS=ULALPS(PARP(34)*Q2)
      FACK=1.
      FACA=1.
      IF(MSTP(33).EQ.1) THEN
        FACK=PARP(31)
      ELSEIF(MSTP(33).EQ.2) THEN
        FACK=PARP(31)
        FACA=PARP(32)/PARP(31)
      ELSEIF(MSTP(33).EQ.3) THEN
        Q2AS=PARP(33)*Q2
        IF(ISTSB.EQ.9.AND.MSTP(82).GE.2) Q2AS=Q2AS+
     &  PARU(112)*PARP(82)
        AS=ULALPS(Q2AS)
      ENDIF
      VINT(138)=1.
      VINT(57)=AEM
      VINT(58)=AS
 
C...Set flags for allowed reacting partons/leptons.
      DO 140 I=1,2
      DO 120 J=-25,25
      KFAC(I,J)=0
  120 CONTINUE
      IF(MINT(44+I).EQ.1) THEN
        KFAC(I,MINT(10+I))=1
      ELSEIF(MINT(40+I).EQ.1.AND.MSTP(12).EQ.0) THEN
        KFAC(I,MINT(10+I))=1
        KFAC(I,22)=1
        KFAC(I,24)=1
        KFAC(I,-24)=1
      ELSE
        DO 130 J=-25,25
        KFAC(I,J)=KFIN(I,J)
        IF(IABS(J).GT.MSTP(58).AND.IABS(J).LE.10) KFAC(I,J)=0
        IF(XSFX(I,J).LT.1E-10) KFAC(I,J)=0
  130   CONTINUE
      ENDIF
  140 CONTINUE
 
C...Lower and upper limit for fermion flavour loops.
      MMIN1=0
      MMAX1=0
      MMIN2=0
      MMAX2=0
      DO 150 J=-20,20
      IF(KFAC(1,-J).EQ.1) MMIN1=-J
      IF(KFAC(1,J).EQ.1) MMAX1=J
      IF(KFAC(2,-J).EQ.1) MMIN2=-J
      IF(KFAC(2,J).EQ.1) MMAX2=J
  150 CONTINUE
      MMINA=MIN(MMIN1,MMIN2)
      MMAXA=MAX(MMAX1,MMAX2)
 
C...Common conversion factors (including Jacobian) for subprocesses.
      SQMZ=PMAS(23,1)**2
      SQMW=PMAS(24,1)**2
      SQMH=PMAS(KFHIGG,1)**2
      GMMH=PMAS(KFHIGG,1)*PMAS(KFHIGG,2)
      SQMZP=PMAS(32,1)**2
      SQMWP=PMAS(34,1)**2
      SQMHC=PMAS(37,1)**2
      SQMLQ=PMAS(39,1)**2
      SQMR=PMAS(40,1)**2
 
C...Phase space integral in tau.
      COMFAC=PARU(1)*PARU(5)/VINT(2)
      IF(MINT(41).EQ.2.AND.MINT(42).EQ.2) COMFAC=COMFAC*FACK
      IF((MINT(47).GE.2.OR.(ISTSB.GE.3.AND.ISTSB.LE.5)).AND.
     &ISTSB.NE.9) THEN
        ATAU1=LOG(TAUMAX/TAUMIN)
        ATAU2=(TAUMAX-TAUMIN)/(TAUMAX*TAUMIN)
        H1=COEF(ISUBSV,1)+(ATAU1/ATAU2)*COEF(ISUBSV,2)/TAU
        IF(MINT(72).GE.1) THEN
          TAUR1=VINT(73)
          GAMR1=VINT(74)
          ATAUD=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR1)/(TAUMAX+TAUR1))
          ATAU3=ATAUD/TAUR1
          IF(ATAUD.GT.1E-6) H1=H1+
     &    (ATAU1/ATAU3)*COEF(ISUBSV,3)/(TAU+TAUR1)
          ATAUD=ATAN((TAUMAX-TAUR1)/GAMR1)-ATAN((TAUMIN-TAUR1)/GAMR1)
          ATAU4=ATAUD/GAMR1
          IF(ATAUD.GT.1E-6) H1=H1+
     &    (ATAU1/ATAU4)*COEF(ISUBSV,4)*TAU/((TAU-TAUR1)**2+GAMR1**2)
        ENDIF
        IF(MINT(72).EQ.2) THEN
          TAUR2=VINT(75)
          GAMR2=VINT(76)
          ATAUD=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR2)/(TAUMAX+TAUR2))
          ATAU5=ATAUD/TAUR2
          IF(ATAUD.GT.1E-6) H1=H1+
     &    (ATAU1/ATAU5)*COEF(ISUBSV,5)/(TAU+TAUR2)
          ATAUD=ATAN((TAUMAX-TAUR2)/GAMR2)-ATAN((TAUMIN-TAUR2)/GAMR2)
          ATAU6=ATAUD/GAMR2
          IF(ATAUD.GT.1E-6) H1=H1+
     &    (ATAU1/ATAU6)*COEF(ISUBSV,6)*TAU/((TAU-TAUR2)**2+GAMR2**2)
        ENDIF
        IF(MINT(47).EQ.5.AND.(ISTSB.LE.2.OR.ISTSB.GE.6)) THEN
          ATAU7=LOG(MAX(2E-6,1.-TAUMIN)/MAX(2E-6,1.-TAUMAX))
          IF(ATAU7.GT.1E-6) H1=H1+(ATAU1/ATAU7)*COEF(ISUBSV,7)*TAU/
     &    MAX(2E-6,1.-TAU)
        ENDIF
        COMFAC=COMFAC*ATAU1/(TAU*H1)
      ENDIF
 
C...Phase space integral in y*.
      IF(MINT(47).GE.4.AND.ISTSB.NE.9) THEN
        AYST0=YSTMAX-YSTMIN
        IF(AYST0.LT.1E-6) THEN
          COMFAC=0.
        ELSE
          AYST1=0.5*(YSTMAX-YSTMIN)**2
          AYST2=AYST1
          AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))
          H2=(AYST0/AYST1)*COEF(ISUBSV,8)*(YST-YSTMIN)+
     &    (AYST0/AYST2)*COEF(ISUBSV,9)*(YSTMAX-YST)+
     &    (AYST0/AYST3)*COEF(ISUBSV,10)/COSH(YST)
          IF(MINT(45).EQ.3) THEN
            YST0=-0.5*LOG(TAUE)
            AYST4=LOG(MAX(1E-6,EXP(YST0-YSTMIN)-1.)/
     &      MAX(1E-6,EXP(YST0-YSTMAX)-1.))
            IF(AYST4.GT.1E-6) H2=H2+(AYST0/AYST4)*COEF(ISUBSV,11)/
     &      MAX(1E-6,1.-EXP(YST-YST0))
          ENDIF
          IF(MINT(46).EQ.3) THEN
            YST0=-0.5*LOG(TAUE)
            AYST5=LOG(MAX(1E-6,EXP(YST0+YSTMAX)-1.)/
     &      MAX(1E-6,EXP(YST0+YSTMIN)-1.))
            IF(AYST5.GT.1E-6) H2=H2+(AYST0/AYST5)*COEF(ISUBSV,12)/
     &      MAX(1E-6,1.-EXP(-YST-YST0))
          ENDIF
          COMFAC=COMFAC*AYST0/H2
        ENDIF
      ENDIF
 
C...2 -> 1 processes: reduction in angular part of phase space integral
C...for case of decaying resonance.
      ACTH0=CTNMAX-CTNMIN+CTPMAX-CTPMIN
      IF((ISTSB.EQ.1.OR.ISTSB.EQ.3.OR.ISTSB.EQ.5)) THEN
        IF(MDCY(KFPR(ISUBSV,1),1).EQ.1) THEN
          IF(KFPR(ISUB,1).EQ.25.OR.KFPR(ISUB,1).EQ.37.OR.
     &    KFPR(ISUB,1).EQ.39) THEN
            COMFAC=COMFAC*0.5*ACTH0
          ELSE
            COMFAC=COMFAC*0.125*(3.*ACTH0+CTNMAX**3-CTNMIN**3+
     &      CTPMAX**3-CTPMIN**3)
          ENDIF
        ENDIF
 
C...2 -> 2 processes: angular part of phase space integral.
      ELSEIF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) THEN
        ACTH1=LOG((MAX(RM34,RSQM-CTNMIN)*MAX(RM34,RSQM-CTPMIN))/
     &  (MAX(RM34,RSQM-CTNMAX)*MAX(RM34,RSQM-CTPMAX)))
        ACTH2=LOG((MAX(RM34,RSQM+CTNMAX)*MAX(RM34,RSQM+CTPMAX))/
     &  (MAX(RM34,RSQM+CTNMIN)*MAX(RM34,RSQM+CTPMIN)))
        ACTH3=1./MAX(RM34,RSQM-CTNMAX)-1./MAX(RM34,RSQM-CTNMIN)+
     &  1./MAX(RM34,RSQM-CTPMAX)-1./MAX(RM34,RSQM-CTPMIN)
        ACTH4=1./MAX(RM34,RSQM+CTNMIN)-1./MAX(RM34,RSQM+CTNMAX)+
     &  1./MAX(RM34,RSQM+CTPMIN)-1./MAX(RM34,RSQM+CTPMAX)
        H3=COEF(ISUBSV,13)+
     &  (ACTH0/ACTH1)*COEF(ISUBSV,14)/MAX(RM34,RSQM-CTH)+
     &  (ACTH0/ACTH2)*COEF(ISUBSV,15)/MAX(RM34,RSQM+CTH)+
     &  (ACTH0/ACTH3)*COEF(ISUBSV,16)/MAX(RM34,RSQM-CTH)**2+
     &  (ACTH0/ACTH4)*COEF(ISUBSV,17)/MAX(RM34,RSQM+CTH)**2
        COMFAC=COMFAC*ACTH0*0.5*BE34/H3
 
C...2 -> 2 processes: take into account final state Breit-Wigners.
        COMFAC=COMFAC*VINT(80)
      ENDIF
 
C...2 -> 3, 4 processes: phace space integral in tau'.
      IF(MINT(47).GE.2.AND.ISTSB.GE.3.AND.ISTSB.LE.5) THEN
        ATAUP1=LOG(TAUPMX/TAUPMN)
        ATAUP2=((1.-TAU/TAUPMX)**4-(1.-TAU/TAUPMN)**4)/(4.*TAU)
        H4=COEF(ISUBSV,18)+
     &  (ATAUP1/ATAUP2)*COEF(ISUBSV,19)*(1.-TAU/TAUP)**3/TAUP
        IF(MINT(47).EQ.5) THEN
          ATAUP3=LOG(MAX(2E-6,1.-TAUPMN)/MAX(2E-6,1.-TAUPMX))
          H4=H4+(ATAUP1/ATAUP3)*COEF(ISUBSV,20)*TAUP/MAX(2E-6,1.-TAUP)
        ENDIF
        COMFAC=COMFAC*ATAUP1/H4
      ENDIF
 
C...2 -> 3, 4 processes: effective W/Z structure functions.
      IF(ISTSB.EQ.3.OR.ISTSB.EQ.4) THEN
        IF(1.-TAU/TAUP.GT.1.E-4) THEN
          FZW=(1.+TAU/TAUP)*LOG(TAUP/TAU)-2.*(1.-TAU/TAUP)
        ELSE
          FZW=1./6.*(1.-TAU/TAUP)**3*TAU/TAUP
        ENDIF
        COMFAC=COMFAC*FZW
      ENDIF
 
C...2 -> 3 processes: phase space integrals for pT1, pT2, y3, mirror.
      IF(ISTSB.EQ.5) THEN
        COMFAC=COMFAC*VINT(205)*VINT(210)*VINT(212)*VINT(214)/
     &  (128.*PARU(1)**4*VINT(220))*(TAU**2/TAUP)
      ENDIF
 
C...2 -> 2 processes: optional dampening by pT^4/(pT0^2+pT^2)^2.
      IF(MSTP(85).EQ.1.AND.MOD(ISTSB,2).EQ.0) COMFAC=COMFAC*
     &SQPTH**2/(PARP(82)**2+SQPTH)**2
 
C...gamma + gamma: include factor 2 when different nature.
      IF(MINT(11).EQ.22.AND.MINT(12).EQ.22.AND.MINT(123).GE.4)
     &COMFAC=2.*COMFAC
 
C...Phase space integral for low-pT and multiple interactions.
      IF(ISTSB.EQ.9) THEN
        COMFAC=PARU(1)*PARU(5)*FACK*0.5*VINT(2)/SH2
        ATAU1=LOG(2.*(1.+SQRT(1.-XT2))/XT2-1.)
        ATAU2=2.*ATAN(1./XT2-1.)/SQRT(XT2)
        H1=COEF(ISUBSV,1)+(ATAU1/ATAU2)*COEF(ISUBSV,2)/SQRT(TAU)
        COMFAC=COMFAC*ATAU1/H1
        AYST0=YSTMAX-YSTMIN
        AYST1=0.5*(YSTMAX-YSTMIN)**2
        AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))
        H2=(AYST0/AYST1)*COEF(ISUBSV,8)*(YST-YSTMIN)+
     &  (AYST0/AYST1)*COEF(ISUBSV,9)*(YSTMAX-YST)+
     &  (AYST0/AYST3)*COEF(ISUBSV,10)/COSH(YST)
        COMFAC=COMFAC*AYST0/H2
        IF(MSTP(82).LE.1) COMFAC=COMFAC*XT2**2*(1./VINT(149)-1.)
C...For MSTP(82)>=2 an additional factor (xT2/(xT2+VINT(149))**2 is
C...introduced to make cross-section finite for xT2 -> 0.
        IF(MSTP(82).GE.2) COMFAC=COMFAC*XT2**2/(VINT(149)*
     &  (1.+VINT(149)))
      ENDIF
 
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
      IF((MSTP(46).GE.3.AND.MSTP(46).LE.6).AND.(ISUB.EQ.71.OR.ISUB.EQ.
     &72.OR.ISUB.EQ.73.OR.ISUB.EQ.76.OR.ISUB.EQ.77)) THEN
C...Calculate M_R and N_R functions for Higgs-like and QCD-like models.
        IF(MSTP(46).LE.4) THEN
          HDTLH=LOG(PMAS(25,1)/PARP(44))
          HDTMR=(4.5*PARU(1)/SQRT(3.)-74./9.)/8.+HDTLH/12.
          HDTNR=-1./18.+HDTLH/6.
        ELSE
          HDTNM=0.125*(1./(288.*PARU(1)**2)+(PARP(47)/PARP(45))**2)
          HDTLQ=LOG(PARP(45)/PARP(44))
          HDTMR=-(4.*PARU(1))**2*0.5*HDTNM+HDTLQ/12.
          HDTNR=(4.*PARU(1))**2*HDTNM+HDTLQ/6.
        ENDIF
 
C...Calculate lowest and next-to-lowest order partial wave amplitudes.
        HDTV=1./(16.*PARU(1)*PARP(47)**2)
        A00L=HDTV*SH
        A20L=-0.5*A00L
        A11L=A00L/6.
        HDTLS=LOG(SH/PARP(44)**2)
        A004=(HDTV*SH)**2/(4.*PARU(1))*CMPLX((176.*HDTMR+112.*HDTNR)/3.+
     &  11./27.-(50./9.)*HDTLS,4.*PARU(1))
        A204=(HDTV*SH)**2/(4.*PARU(1))*CMPLX(32.*(HDTMR+2.*HDTNR)/3.+
     &  25./54.-(20./9.)*HDTLS,PARU(1))
        A114=(HDTV*SH)**2/(6.*PARU(1))*CMPLX(4.*(-2.*HDTMR+HDTNR)-
     &  1./18.,PARU(1)/6.)
 
C...Unitarize partial wave amplitudes with Pade or K-matrix method.
        IF(MSTP(46).EQ.3.OR.MSTP(46).EQ.5) THEN
          A00U=A00L/(1.-A004/A00L)
          A20U=A20L/(1.-A204/A20L)
          A11U=A11L/(1.-A114/A11L)
        ELSE
          A00U=(A00L+REAL(A004))/(1.-CMPLX(0.,A00L+REAL(A004)))
          A20U=(A20L+REAL(A204))/(1.-CMPLX(0.,A20L+REAL(A204)))
          A11U=(A11L+REAL(A114))/(1.-CMPLX(0.,A11L+REAL(A114)))
        ENDIF
      ENDIF
 
C...A: 2 -> 1, tree diagrams.
 
  160 IF(ISUB.LE.10) THEN
      IF(ISUB.EQ.1) THEN
C...f + f~ -> gamma*/Z0.
        MINT(61)=2
        CALL PYWIDT(23,SH,WDTP,WDTE)
        HP0=AEM/3.*SH
        HP1=AEM/3.*XWC*SH
        HS=HP1*WDTP(0)
        FACZ=4.*COMFAC*3.
        DO 170 I=MMINA,MMAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 170
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        HI0=HP0
        IF(IABS(I).LE.10) HI0=HI0*FACA/3.
        HI1=HP1
        IF(IABS(I).LE.10) HI1=HI1*FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZ*(EI**2/SH2*HI0*HP0*VINT(111)+EI*VI*(1.-SQMZ/SH)/
     &  ((SH-SQMZ)**2+HS**2)*(HI0*HP1+HI1*HP0)*VINT(112)+
     &  (VI**2+AI**2)/((SH-SQMZ)**2+HS**2)*HI1*HP1*VINT(114))
  170   CONTINUE
 
      ELSEIF(ISUB.EQ.2) THEN
C...f + f~' -> W+/-.
        CALL PYWIDT(24,SH,WDTP,WDTE)
        HP=AEM/(24.*XW)*SH
        HS=HP*WDTP(0)
        FACBW=4.*COMFAC/((SH-SQMW)**2+HS**2)*3.
        DO 190 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 190
        IA=IABS(I)
        DO 180 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 180
        JA=IABS(J)
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 180
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 180
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        HI=HP*2.
        IF(IA.LE.10) HI=HI*VCKM((IA+1)/2,(JA+1)/2)*FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        HF=HP*(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))
        SIGH(NCHN)=HI*FACBW*HF
  180   CONTINUE
  190   CONTINUE
 
      ELSEIF(ISUB.EQ.3) THEN
C...f + f~ -> H0 (or H'0, or A0).
        CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=4.*COMFAC/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        DO 200 I=MMINA,MMAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 200
        IA=IABS(I)
        RMQ=PMAS(IA,1)**2/SH
        HI=HP*RMQ
        IF(IA.LE.10) HI=HP*RMQ*FACA/3.
        IF(IA.LE.10.AND.MSTP(37).EQ.1.AND.MSTP(2).GE.1) HI=HI*
     &  (LOG(MAX(4.,PARP(37)**2*RMQ*SH/PARU(117)**2))/
     &  LOG(MAX(4.,SH/PARU(117)**2)))**(24./(33.-2.*MSTU(118)))
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
          IKFI=1
          IF(IA.LE.10.AND.MOD(IA,2).EQ.0) IKFI=2
          IF(IA.GT.10) IKFI=3
          HI=HI*PARU(150+10*IHIGG+IKFI)**2
        ENDIF
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
  200   CONTINUE
 
      ELSEIF(ISUB.EQ.4) THEN
C...gamma + W+/- -> W+/-.
 
      ELSEIF(ISUB.EQ.5) THEN
C...Z0 + Z0 -> H0.
        CALL PYWIDT(25,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=4.*COMFAC/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        HI=HP/4.
        FACI=8./(PARU(1)**2*XW1)*(AEM*XWC)**2
        DO 220 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 220
        DO 210 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 210
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        EJ=KCHG(IABS(J),1)/3.
        AJ=SIGN(1.,EJ)
        VJ=AJ-4.*EJ*XWV
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACI*(VI**2+AI**2)*(VJ**2+AJ**2)*HI*FACBW*HF
  210   CONTINUE
  220   CONTINUE
 
      ELSEIF(ISUB.EQ.6) THEN
C...Z0 + W+/- -> W+/-.
 
      ELSEIF(ISUB.EQ.7) THEN
C...W+ + W- -> Z0.
 
      ELSEIF(ISUB.EQ.8) THEN
C...W+ + W- -> H0.
        CALL PYWIDT(25,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=4.*COMFAC/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        HI=HP/2.
        FACI=1./(4.*PARU(1)**2)*(AEM/XW)**2
        DO 240 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 240
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)
        DO 230 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 230
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)
        IF(EI*EJ.GT.0.) GOTO 230
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACI*VINT(180+I)*VINT(180+J)*HI*FACBW*HF
  230   CONTINUE
  240   CONTINUE
 
C...B: 2 -> 2, tree diagrams.
 
      ELSEIF(ISUB.EQ.10) THEN
C...f + f' -> f + f' (gamma/Z/W exchange).
        FACGGF=COMFAC*AEM**2*2.*(SH2+UH2)/TH2
        FACGZF=COMFAC*AEM**2*XWC*4.*SH2/(TH*(TH-SQMZ))
        FACZZF=COMFAC*(AEM*XWC)**2*2.*SH2/(TH-SQMZ)**2
        FACWWF=COMFAC*(0.5*AEM/XW)**2*SH2/(TH-SQMW)**2
        DO 260 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 260
        IA=IABS(I)
        DO 250 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 250
        JA=IABS(J)
C...Electroweak couplings.
        EI=KCHG(IA,1)*ISIGN(1,I)/3.
        AI=SIGN(1.,KCHG(IA,1)+0.5)*ISIGN(1,I)
        VI=AI-4.*EI*XWV
        EJ=KCHG(JA,1)*ISIGN(1,J)/3.
        AJ=SIGN(1.,KCHG(JA,1)+0.5)*ISIGN(1,J)
        VJ=AJ-4.*EJ*XWV
        EPSIJ=ISIGN(1,I*J)
C...gamma/Z exchange, only gamma exchange, or only Z exchange.
        IF(MSTP(21).GE.1.AND.MSTP(21).LE.4) THEN
          IF(MSTP(21).EQ.1.OR.MSTP(21).EQ.4) THEN
            FACNCF=FACGGF*EI**2*EJ**2+FACGZF*EI*EJ*
     &      (VI*VJ*(1.+UH2/SH2)+AI*AJ*EPSIJ*(1.-UH2/SH2))+
     &      FACZZF*((VI**2+AI**2)*(VJ**2+AJ**2)*(1.+UH2/SH2)+
     &      4.*VI*VJ*AI*AJ*EPSIJ*(1.-UH2/SH2))
          ELSEIF(MSTP(21).EQ.2) THEN
            FACNCF=FACGGF*EI**2*EJ**2
          ELSE
            FACNCF=FACZZF*((VI**2+AI**2)*(VJ**2+AJ**2)*(1.+UH2/SH2)+
     &      4.*VI*VJ*AI*AJ*EPSIJ*(1.-UH2/SH2))
          ENDIF
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACNCF
        ENDIF
C...W exchange.
        IF((MSTP(21).EQ.1.OR.MSTP(21).EQ.5).AND.AI*AJ.LT.0.) THEN
          FACCCF=FACWWF*VINT(180+I)*VINT(180+J)
          IF(EPSIJ.LT.0.) FACCCF=FACCCF*UH2/SH2
          IF(IA.GT.10.AND.MOD(IA,2).EQ.0) FACCCF=2.*FACCCF
          IF(JA.GT.10.AND.MOD(JA,2).EQ.0) FACCCF=2.*FACCCF
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=2
          SIGH(NCHN)=FACCCF
        ENDIF
  250   CONTINUE
  260   CONTINUE
      ENDIF
 
      ELSEIF(ISUB.LE.20) THEN
      IF(ISUB.EQ.11) THEN
C...f + f' -> f + f' (g exchange).
        FACQQ1=COMFAC*AS**2*4./9.*(SH2+UH2)/TH2
        FACQQB=COMFAC*AS**2*4./9.*((SH2+UH2)/TH2*FACA-
     &  MSTP(34)*2./3.*UH2/(SH*TH))
        FACQQ2=COMFAC*AS**2*4./9.*((SH2+TH2)/UH2-
     &  MSTP(34)*2./3.*SH2/(TH*UH))
        IF(MSTP(5).GE.1) THEN
C...Modifications from contact interactions (compositeness).
          FACCI1=FACQQ1+COMFAC*(SH2/PARU(155)**4)
          FACCIB=FACQQB+COMFAC*(8./9.)*(AS*PARU(156)/PARU(155)**2)*
     &    (UH2/TH+UH2/SH)+COMFAC*(5./3.)*(UH2/PARU(155)**4)
          FACCI2=FACQQ2+COMFAC*(8./9.)*(AS*PARU(156)/PARU(155)**2)*
     &    (SH2/TH+SH2/UH)+COMFAC*(5./3.)*(SH2/PARU(155)**4)
          FACCI3=FACQQ1+COMFAC*(UH2/PARU(155)**4)
        ENDIF
        DO 280 I=MMIN1,MMAX1
        IA=IABS(I)
        IF(I.EQ.0.OR.IA.GT.MSTP(58).OR.KFAC(1,I).EQ.0) GOTO 280
        DO 270 J=MMIN2,MMAX2
        JA=IABS(J)
        IF(J.EQ.0.OR.JA.GT.MSTP(58).OR.KFAC(2,J).EQ.0) GOTO 270
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        IF(MSTP(5).LE.0.OR.(MSTP(5).EQ.1.AND.(IA.GE.3.OR.JA.GE.3)))
     &  THEN
          SIGH(NCHN)=FACQQ1
          IF(I.EQ.-J) SIGH(NCHN)=FACQQB
        ELSE
          SIGH(NCHN)=FACCI1
          IF(I*J.LT.0) SIGH(NCHN)=FACCI3
          IF(I.EQ.-J) SIGH(NCHN)=FACCIB
        ENDIF
        IF(I.EQ.J) THEN
          SIGH(NCHN)=0.5*SIGH(NCHN)
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=2
          IF(MSTP(5).LE.0.OR.(MSTP(5).EQ.1.AND.IA.GE.3)) THEN
            SIGH(NCHN)=0.5*FACQQ2
          ELSE
            SIGH(NCHN)=0.5*FACCI2
          ENDIF
        ENDIF
  270   CONTINUE
  280   CONTINUE
 
      ELSEIF(ISUB.EQ.12) THEN
C...f + f~ -> f' + f~' (q + q~ -> q' + q~' only).
        CALL PYWIDT(21,SH,WDTP,WDTE)
        FACQQB=COMFAC*AS**2*4./9.*(TH2+UH2)/SH2*(WDTE(0,1)+WDTE(0,2)+
     &  WDTE(0,4))
        IF(MSTP(5).EQ.1) THEN
C...Modifications from contact interactions (compositeness).
          FACCIB=FACQQB
          DO 290 I=1,2
          FACCIB=FACCIB+COMFAC*(UH2/PARU(155)**4)*(WDTE(I,1)+WDTE(I,2)+
     &    WDTE(I,4))
  290     CONTINUE
        ELSEIF(MSTP(5).GE.2) THEN
          FACCIB=FACQQB+COMFAC*(UH2/PARU(155)**4)*(WDTE(0,1)+WDTE(0,2)+
     &    WDTE(0,4))
        ENDIF
        DO 300 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 300
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        IF(MSTP(5).LE.0.OR.(MSTP(5).EQ.1.AND.IABS(I).GE.3)) THEN
          SIGH(NCHN)=FACQQB
        ELSE
          SIGH(NCHN)=FACCIB
        ENDIF
  300   CONTINUE
 
      ELSEIF(ISUB.EQ.13) THEN
C...f + f~ -> g + g (q + q~ -> g + g only).
        FACGG1=COMFAC*AS**2*32./27.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)
        FACGG2=COMFAC*AS**2*32./27.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)
        DO 310 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 310
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACGG1
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=2
        SIGH(NCHN)=0.5*FACGG2
  310   CONTINUE
 
      ELSEIF(ISUB.EQ.14) THEN
C...f + f~ -> g + gamma (q + q~ -> g + gamma only).
        FACGG=COMFAC*AS*AEM*8./9.*(TH2+UH2)/(TH*UH)
        DO 320 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 320
        EI=KCHG(IABS(I),1)/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGG*EI**2
  320   CONTINUE
 
      ELSEIF(ISUB.EQ.15) THEN
C...f + f~ -> g + (gamma*/Z0) (q + q~ -> g + (gamma*/Z0) only).
        FACZG=COMFAC*AS*AEM*(8./9.)*(TH2+UH2+2.*SQM4*SH)/(TH*UH)
C...gamma, gamma/Z interference and Z couplings to final fermion pairs.
        HFGG=0.
        HFGZ=0.
        HFZZ=0.
        HBW4=0.
        RADC4=1.+ULALPS(SQM4)/PARU(1)
        DO 330 I=1,MIN(16,MDCY(23,3))
        IDC=I+MDCY(23,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 330
        IMDM=0
        IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2.OR.MDME(IDC,1).EQ.4)
     &  IMDM=1
        IF(I.LE.8) THEN
          EF=KCHG(I,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
        ELSEIF(I.LE.16) THEN
          EF=KCHG(I+2,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
        ENDIF
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
        IF(4.*RM1.LT.1.) THEN
          FCOF=1.
          IF(I.LE.8) FCOF=3.*RADC4
          BE34=SQRT(MAX(0.,1.-4.*RM1))
          IF(IMDM.EQ.1) THEN
            HFGG=HFGG+FCOF*EF**2*(1.+2.*RM1)*BE34
            HFGZ=HFGZ+FCOF*EF*VF*(1.+2.*RM1)*BE34
            HFZZ=HFZZ+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
          ENDIF
          HBW4=HBW4+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
        ENDIF
  330   CONTINUE
C...Propagators: as simulated in PYOFSH and as desired.
        GMMZ=PMAS(23,1)*PMAS(23,2)
        HBW4=HBW4*XWC*SQMZ/((SQM4-SQMZ)**2+GMMZ**2)
        MINT(15)=1
        MINT(61)=1
        CALL PYWIDT(23,SQM4,WDTP,WDTE)
        HFGG=HFGG*VINT(111)/SQM4
        HFGZ=HFGZ*VINT(112)/SQM4
        HFZZ=HFZZ*VINT(114)/SQM4
C...Loop over flavours; consider full gamma/Z structure.
        DO 340 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 340
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZG*(EI**2*HFGG+EI*VI*HFGZ+
     &  (VI**2+AI**2)*HFZZ)/HBW4
  340   CONTINUE
 
      ELSEIF(ISUB.EQ.16) THEN
C...f + f~' -> g + W+/- (q + q~' -> g + W+/- only).
        FACWG=COMFAC*AS*AEM/XW*2./9.*(TH2+UH2+2.*SQM4*SH)/(TH*UH)
C...Propagators: as simulated in PYOFSH and as desired.
        GMMW=PMAS(24,1)*PMAS(24,2)
        HBW4=GMMW/((SQM4-SQMW)**2+GMMW**2)
        CALL PYWIDT(24,SQM4,WDTP,WDTE)
        AEMC=ULALEM(SQM4)
        IF(MSTP(8).GE.1) AEMC=AEM
        GMMWC=SQM4*WDTP(0)*AEMC/(24.*XW)
        HBW4C=GMMWC/((SQM4-SQMW)**2+GMMWC**2)
        FACWG=FACWG*HBW4C/HBW4
        DO 360 I=MMIN1,MMAX1
        IA=IABS(I)
        IF(I.EQ.0.OR.IA.GT.10.OR.KFAC(1,I).EQ.0) GOTO 360
        DO 350 J=MMIN2,MMAX2
        JA=IABS(J)
        IF(J.EQ.0.OR.JA.GT.10.OR.KFAC(2,J).EQ.0) GOTO 350
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 350
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        WIDSC=(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))/WDTP(0)
        FCKM=VCKM((IA+1)/2,(JA+1)/2)
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWG*FCKM*WIDSC
  350   CONTINUE
  360   CONTINUE
 
      ELSEIF(ISUB.EQ.17) THEN
C...f + f~ -> g + H0 (q + q~ -> g + H0 only).
 
      ELSEIF(ISUB.EQ.18) THEN
C...f + f~ -> gamma + gamma.
        FACGG=COMFAC*AEM**2*2.*(TH2+UH2)/(TH*UH)
        DO 370 I=MMINA,MMAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 370
        EI=KCHG(IABS(I),1)/3.
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACGG*FCOI*EI**4
  370   CONTINUE
 
      ELSEIF(ISUB.EQ.19) THEN
C...f + f~ -> gamma + (gamma*/Z0).
        FACGZ=COMFAC*2.*AEM**2*(TH2+UH2+2.*SQM4*SH)/(TH*UH)
C...gamma, gamma/Z interference and Z couplings to final fermion pairs.
        HFGG=0.
        HFGZ=0.
        HFZZ=0.
        HBW4=0.
        RADC4=1.+ULALPS(SQM4)/PARU(1)
        DO 380 I=1,MIN(16,MDCY(23,3))
        IDC=I+MDCY(23,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 380
        IMDM=0
        IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2.OR.MDME(IDC,1).EQ.4)
     &  IMDM=1
        IF(I.LE.8) THEN
          EF=KCHG(I,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
        ELSEIF(I.LE.16) THEN
          EF=KCHG(I+2,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
        ENDIF
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
        IF(4.*RM1.LT.1.) THEN
          FCOF=1.
          IF(I.LE.8) FCOF=3.*RADC4
          BE34=SQRT(MAX(0.,1.-4.*RM1))
          IF(IMDM.EQ.1) THEN
            HFGG=HFGG+FCOF*EF**2*(1.+2.*RM1)*BE34
            HFGZ=HFGZ+FCOF*EF*VF*(1.+2.*RM1)*BE34
            HFZZ=HFZZ+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
          ENDIF
          HBW4=HBW4+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
        ENDIF
  380   CONTINUE
C...Propagators: as simulated in PYOFSH and as desired.
        GMMZ=PMAS(23,1)*PMAS(23,2)
        HBW4=HBW4*XWC*SQMZ/((SQM4-SQMZ)**2+GMMZ**2)
        MINT(15)=1
        MINT(61)=1
        CALL PYWIDT(23,SQM4,WDTP,WDTE)
        HFGG=HFGG*VINT(111)/SQM4
        HFGZ=HFGZ*VINT(112)/SQM4
        HFZZ=HFZZ*VINT(114)/SQM4
C...Loop over flavours; consider full gamma/Z structure.
        DO 390 I=MMINA,MMAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 390
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGZ*FCOI*EI**2*(EI**2*HFGG+EI*VI*HFGZ+
     &  (VI**2+AI**2)*HFZZ)/HBW4
  390   CONTINUE
 
      ELSEIF(ISUB.EQ.20) THEN
C...f + f~' -> gamma + W+/-.
        FACGW=COMFAC*0.5*AEM**2/XW
C...Propagators: as simulated in PYOFSH and as desired.
        GMMW=PMAS(24,1)*PMAS(24,2)
        HBW4=GMMW/((SQM4-SQMW)**2+GMMW**2)
        CALL PYWIDT(24,SQM4,WDTP,WDTE)
        AEMC=ULALEM(SQM4)
        IF(MSTP(8).GE.1) AEMC=AEM
        GMMWC=SQM4*WDTP(0)*AEMC/(24.*XW)
        HBW4C=GMMWC/((SQM4-SQMW)**2+GMMWC**2)
        FACGW=FACGW*HBW4C/HBW4
C...Anomalous couplings.
        TERM1=(TH2+UH2+2.*SQM4*SH)/(TH*UH)
        TERM2=0.
        TERM3=0.
        IF(MSTP(5).GE.1) THEN
          TERM2=PARU(153)*(TH-UH)/(TH+UH)
          TERM3=0.5*PARU(153)**2*(TH*UH+(TH2+UH2)*SH/
     &    (4.*PMAS(24,1)**2))/(TH+UH)**2
        ENDIF
        DO 410 I=MMIN1,MMAX1
        IA=IABS(I)
        IF(I.EQ.0.OR.IA.GT.20.OR.KFAC(1,I).EQ.0) GOTO 410
        DO 400 J=MMIN2,MMAX2
        JA=IABS(J)
        IF(J.EQ.0.OR.JA.GT.20.OR.KFAC(2,J).EQ.0) GOTO 400
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 400
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 400
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        WIDSC=(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))/WDTP(0)
        IF(IA.LE.10) THEN
          FACWR=UH/(TH+UH)-1./3.
          FCKM=VCKM((IA+1)/2,(JA+1)/2)
          FCOI=FACA/3.
        ELSE
          FACWR=-TH/(TH+UH)
          FCKM=1.
          FCOI=1.
        ENDIF
        FACWK=TERM1*FACWR**2+TERM2*FACWR+TERM3
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGW*FACWK*FCOI*FCKM*WIDSC
  400   CONTINUE
  410   CONTINUE
      ENDIF
 
      ELSEIF(ISUB.LE.30) THEN
      IF(ISUB.EQ.21) THEN
C...f + f~ -> gamma + H0.
 
      ELSEIF(ISUB.EQ.22) THEN
C...f + f~ -> (gamma*/Z0) + (gamma*/Z0).
C...Kinematics dependence.
        FACZZ=COMFAC*AEM**2*((TH2+UH2+2.*(SQM3+SQM4)*SH)/(TH*UH)-
     &  SQM3*SQM4*(1./TH2+1./UH2))
C...gamma, gamma/Z interference and Z couplings to final fermion pairs.
        DO 430 I=1,6
        DO 420 J=1,3
        HGZ(I,J)=0.
  420   CONTINUE
  430   CONTINUE
        HBW3=0.
        HBW4=0.
        RADC3=1.+ULALPS(SQM3)/PARU(1)
        RADC4=1.+ULALPS(SQM4)/PARU(1)
        DO 440 I=1,MIN(16,MDCY(23,3))
        IDC=I+MDCY(23,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 440
        IMDM=0
        IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2) IMDM=1
        IF(MDME(IDC,1).EQ.4.OR.MDME(IDC,1).EQ.5) IMDM=MDME(IDC,1)-2
        IF(I.LE.8) THEN
          EF=KCHG(I,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
        ELSEIF(I.LE.16) THEN
          EF=KCHG(I+2,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
        ENDIF
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM3
        IF(4.*RM1.LT.1.) THEN
          FCOF=1.
          IF(I.LE.8) FCOF=3.*RADC3
          BE34=SQRT(MAX(0.,1.-4.*RM1))
          IF(IMDM.GE.1) THEN
            HGZ(1,IMDM)=HGZ(1,IMDM)+FCOF*EF**2*(1.+2.*RM1)*BE34
            HGZ(2,IMDM)=HGZ(2,IMDM)+FCOF*EF*VF*(1.+2.*RM1)*BE34
            HGZ(3,IMDM)=HGZ(3,IMDM)+FCOF*(VF**2*(1.+2.*RM1)+
     &      AF**2*(1.-4.*RM1))*BE34
          ENDIF
          HBW3=HBW3+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
        ENDIF
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
        IF(4.*RM1.LT.1.) THEN
          FCOF=1.
          IF(I.LE.8) FCOF=3.*RADC4
          BE34=SQRT(MAX(0.,1.-4.*RM1))
          IF(IMDM.GE.1) THEN
            HGZ(4,IMDM)=HGZ(4,IMDM)+FCOF*EF**2*(1.+2.*RM1)*BE34
            HGZ(5,IMDM)=HGZ(5,IMDM)+FCOF*EF*VF*(1.+2.*RM1)*BE34
            HGZ(6,IMDM)=HGZ(6,IMDM)+FCOF*(VF**2*(1.+2.*RM1)+
     &      AF**2*(1.-4.*RM1))*BE34
          ENDIF
          HBW4=HBW4+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
        ENDIF
  440   CONTINUE
C...Propagators: as simulated in PYOFSH and as desired.
        GMMZ=PMAS(23,1)*PMAS(23,2)
        HBW3=HBW3*XWC*SQMZ/((SQM3-SQMZ)**2+GMMZ**2)
        HBW4=HBW4*XWC*SQMZ/((SQM4-SQMZ)**2+GMMZ**2)
        MINT(15)=1
        MINT(61)=1
        CALL PYWIDT(23,SQM3,WDTP,WDTE)
        DO 450 J=1,3
        HGZ(1,J)=HGZ(1,J)*VINT(111)/SQM3
        HGZ(2,J)=HGZ(2,J)*VINT(112)/SQM3
        HGZ(3,J)=HGZ(3,J)*VINT(114)/SQM3
  450   CONTINUE
        MINT(61)=1
        CALL PYWIDT(23,SQM4,WDTP,WDTE)
        DO 460 J=1,3
        HGZ(4,J)=HGZ(4,J)*VINT(111)/SQM4
        HGZ(5,J)=HGZ(5,J)*VINT(112)/SQM4
        HGZ(6,J)=HGZ(6,J)*VINT(114)/SQM4
  460   CONTINUE
C...Loop over flavours; separate left- and right-handed couplings.
        DO 480 I=MMINA,MMAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 480
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        VALI=VI-AI
        VARI=VI+AI
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        DO 470 J=1,3
        HL3(J)=EI**2*HGZ(1,J)+EI*VALI*HGZ(2,J)+VALI**2*HGZ(3,J)
        HR3(J)=EI**2*HGZ(1,J)+EI*VARI*HGZ(2,J)+VARI**2*HGZ(3,J)
        HL4(J)=EI**2*HGZ(4,J)+EI*VALI*HGZ(5,J)+VALI**2*HGZ(6,J)
        HR4(J)=EI**2*HGZ(4,J)+EI*VARI*HGZ(5,J)+VARI**2*HGZ(6,J)
  470   CONTINUE
        FACLR=HL3(1)*HL4(1)+HL3(1)*(HL4(2)+HL4(3))+
     &  HL4(1)*(HL3(2)+HL3(3))+HL3(2)*HL4(3)+HL4(2)*HL3(3)+
     &  HR3(1)*HR4(1)+HR3(1)*(HR4(2)+HR4(3))+
     &  HR4(1)*(HR3(2)+HR3(3))+HR3(2)*HR4(3)+HR4(2)*HR3(3)
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACZZ*FCOI*FACLR/(HBW3*HBW4)
  480   CONTINUE
 
      ELSEIF(ISUB.EQ.23) THEN
C...f + f~' -> Z0 + W+/-.
        FACZW=COMFAC*0.5*(AEM/XW)**2
        FACZW=FACZW*WIDS(23,2)
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
        FACBW=1./((SH-SQMW)**2+SQMW*PMAS(24,2)**2)
        DO 500 I=MMIN1,MMAX1
        IA=IABS(I)
        IF(I.EQ.0.OR.IA.GT.20.OR.KFAC(1,I).EQ.0) GOTO 500
        DO 490 J=MMIN2,MMAX2
        JA=IABS(J)
        IF(J.EQ.0.OR.JA.GT.20.OR.KFAC(2,J).EQ.0) GOTO 490
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 490
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 490
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        EI=KCHG(IA,1)/3.
        AI=SIGN(1.,EI+0.1)
        VI=AI-4.*EI*XWV
        EJ=KCHG(JA,1)/3.
        AJ=SIGN(1.,EJ+0.1)
        VJ=AJ-4.*EJ*XWV
        IF(VI+AI.GT.0) THEN
          VISAV=VI
          AISAV=AI
          VI=VJ
          AI=AJ
          VJ=VISAV
          AJ=AISAV
        ENDIF
        FCKM=1.
        IF(IA.LE.10) FCKM=VCKM((IA+1)/2,(JA+1)/2)
        FCOI=1.
        IF(IA.LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZW*FCOI*FCKM*(FACBW*((9.-8.*XW)/4.*THUH+
     &  (8.*XW-6.)/4.*SH*(SQM3+SQM4))+(THUH-SH*(SQM3+SQM4))*
     &  (SH-SQMW)*FACBW*0.5*((VJ+AJ)/TH-(VI+AI)/UH)+
     &  THUH/(16.*XW1)*((VJ+AJ)**2/TH2+(VI+AI)**2/UH2)+
     &  SH*(SQM3+SQM4)/(8.*XW1)*(VI+AI)*(VJ+AJ)/(TH*UH))*
     &  WIDS(24,(5-KCHW)/2)
  490   CONTINUE
  500   CONTINUE
 
      ELSEIF(ISUB.EQ.24) THEN
C...f + f~ -> Z0 + H0 (or H'0, or A0).
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
        FACHZ=COMFAC*8.*(AEM*XWC)**2*
     &  (THUH+2.*SH*SQM3)/((SH-SQMZ)**2+SQMZ*PMAS(23,2)**2)
        FACHZ=FACHZ*WIDS(23,2)*WIDS(KFHIGG,2)
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACHZ=FACHZ*
     &  PARU(154+10*IHIGG)**2
        DO 510 I=MMINA,MMAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 510
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACHZ*FCOI*(VI**2+AI**2)
  510   CONTINUE
 
      ELSEIF(ISUB.EQ.25) THEN
C...f + f~ -> W+ + W-.
C...Propagators: Z0, W+- as simulated in PYOFSH and as desired.
        CALL PYWIDT(23,SH,WDTP,WDTE)
        GMMZC=AEM/(48.*XW*XW1)*SH*WDTP(0)
        HBWZC=SH**2/((SH-SQMZ)**2+GMMZC**2)
        GMMW=PMAS(24,1)*PMAS(24,2)
        HBW3=GMMW/((SQM3-SQMW)**2+GMMW**2)
        AEM3=ULALEM(SQM3)
        IF(MSTP(8).GE.1) AEM3=AEM
        CALL PYWIDT(24,SQM3,WDTP,WDTE)
        GMMW3=AEM3/(24.*XW)*SQM3*WDTP(0)
        HBW3C=GMMW3/((SQM3-SQMW)**2+GMMW3**2)
        HBW4=GMMW/((SQM4-SQMW)**2+GMMW**2)
        AEM4=ULALEM(SQM4)
        IF(MSTP(8).GE.1) AEM4=AEM
        CALL PYWIDT(24,SQM4,WDTP,WDTE)
        GMMW4=AEM4/(24.*XW)*SQM4*WDTP(0)
        HBW4C=GMMW4/((SQM4-SQMW)**2+GMMW4**2)
C...Kinematical functions.
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
        THUH34=(2.*SH*(SQM3+SQM4)+THUH)/(SQM3*SQM4)
        GS=(((SH-SQM3-SQM4)**2-4.*SQM3*SQM4)*THUH34+12.*THUH)/SH2
        GT=THUH34+4.*THUH/TH2
        GST=((SH-SQM3-SQM4)*THUH34+4.*(SH*(SQM3+SQM4)-THUH)/TH)/SH
        GU=THUH34+4.*THUH/UH2
        GSU=((SH-SQM3-SQM4)*THUH34+4.*(SH*(SQM3+SQM4)-THUH)/UH)/SH
C...Common factors and couplings.
        FACWW=COMFAC*(HBW3C/HBW3)*(HBW4C/HBW4)
        FACWW=FACWW*WIDS(24,1)
        CGG=AEM**2/2.
        CGZ=AEM**2/(4.*XW)*HBWZC*(1.-SQMZ/SH)
        CZZ=AEM**2/(32.*XW**2)*HBWZC
        CNG=AEM**2/(4.*XW)
        CNZ=AEM**2/(16.*XW**2)*HBWZC*(1.-SQMZ/SH)
        CNN=AEM**2/(16.*XW**2)
C...Coulomb factor for W+W- pair.
        IF(MSTP(40).GE.1.AND.MSTP(40).LE.3) THEN
          COULE=(SH-4.*SQMW)/(4.*PMAS(24,1))
          COULP=MAX(1E-10,0.5*BE34*SQRT(SH))
          IF(COULE.LT.100.*PMAS(24,2)) THEN
            COULP1=SQRT(0.5*PMAS(24,1)*(SQRT(COULE**2+PMAS(24,2)**2)-
     &      COULE))
          ELSE
            COULP1=SQRT(0.5*PMAS(24,1)*(0.5*PMAS(24,2)**2/COULE))
          ENDIF
          IF(COULE.GT.-100.*PMAS(24,2)) THEN
            COULP2=SQRT(0.5*PMAS(24,1)*(SQRT(COULE**2+PMAS(24,2)**2)+
     &      COULE))
          ELSE
            COULP2=SQRT(0.5*PMAS(24,1)*(0.5*PMAS(24,2)**2/ABS(COULE)))
          ENDIF
          IF(MSTP(40).EQ.1) THEN
            COULDC=PARU(1)-2.*ATAN((COULP1**2+COULP2**2-COULP**2)/
     &      MAX(1E-10,2.*COULP*COULP1)) 
            FACCOU=1.+0.5*PARU(101)*COULDC/MAX(1E-5,BE34)
          ELSEIF(MSTP(40).EQ.2) THEN
            COULCK=CMPLX(COULP1,COULP2)
            COULCP=CMPLX(0.,COULP)
            COULCD=(COULCK+COULCP)/(COULCK-COULCP)
            COULCR=1.+(PARU(101)*SQRT(SH))/(4.*COULCP)*LOG(COULCD)
            COULCS=CMPLX(0.,0.)
            NSTP=100
            DO 515 ISTP=1,NSTP
            COULXX=(ISTP-0.5)/NSTP
            COULCS=COULCS+(1./COULXX)*LOG((1.+COULXX*COULCD)/
     &      (1.+COULXX/COULCD))
  515       CONTINUE
            COULCR=COULCR+(PARU(101)**2*SH)/(16.*COULCP*COULCK)*
     &      (COULCS/NSTP)
            FACCOU=ABS(COULCR)**2
          ELSEIF(MSTP(40).EQ.3) THEN
            COULDC=PARU(1)-2.*(1.-BE34)**2*ATAN((COULP1**2+COULP2**2-
     &      COULP**2)/MAX(1E-10,2.*COULP*COULP1)) 
            FACCOU=1.+0.5*PARU(101)*COULDC/MAX(1E-5,BE34)
          ENDIF
        ELSEIF(MSTP(40).EQ.4) THEN
          FACCOU=1.+0.5*PARU(101)*PARU(1)/MAX(1E-5,BE34)
        ELSE
          FACCOU=1.
        ENDIF
        VINT(95)=FACCOU
        FACWW=FACWW*FACCOU
C...Loop over allowed flavours.
        DO 520 I=MMINA,MMAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 520
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI+0.1)
        VI=AI-4.*EI*XWV
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        IF(AI.LT.0.) THEN
          DSIGWW=(CGG*EI**2+CGZ*VI*EI+CZZ*(VI**2+AI**2))*GS+
     &    (CNG*EI+CNZ*(VI+AI))*GST+CNN*GT
        ELSE
          DSIGWW=(CGG*EI**2+CGZ*VI*EI+CZZ*(VI**2+AI**2))*GS-
     &    (CNG*EI+CNZ*(VI+AI))*GSU+CNN*GU
        ENDIF
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWW*FCOI*DSIGWW
  520   CONTINUE
 
      ELSEIF(ISUB.EQ.26) THEN
C...f + f~' -> W+/- + H0 (or H'0, or A0).
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
        FACHW=COMFAC*0.125*(AEM/XW)**2*(THUH+2.*SH*SQM3)/
     &  ((SH-SQMW)**2+SQMW*PMAS(24,2)**2)
        FACHW=FACHW*WIDS(KFHIGG,2)
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACHW=FACHW*
     &  PARU(155+10*IHIGG)**2
        DO 540 I=MMIN1,MMAX1
        IA=IABS(I)
        IF(I.EQ.0.OR.IA.GT.20.OR.KFAC(1,I).EQ.0) GOTO 540
        DO 530 J=MMIN2,MMAX2
        JA=IABS(J)
        IF(J.EQ.0.OR.JA.GT.20.OR.KFAC(1,J).EQ.0) GOTO 530
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 530
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 530
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        FCKM=1.
        IF(IA.LE.10) FCKM=VCKM((IA+1)/2,(JA+1)/2)
        FCOI=1.
        IF(IA.LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACHW*FCOI*FCKM*WIDS(24,(5-KCHW)/2)
  530   CONTINUE
  540   CONTINUE
 
      ELSEIF(ISUB.EQ.27) THEN
C...f + f~ -> H0 + H0.
 
      ELSEIF(ISUB.EQ.28) THEN
C...f + g -> f + g (q + g -> q + g only).
        FACQG1=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*UH2/TH2-UH/SH)*
     &  FACA
        FACQG2=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*SH2/TH2-SH/UH)
        DO 560 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.10) GOTO 560
        DO 550 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 550
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 550
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQG1
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=2
        SIGH(NCHN)=FACQG2
  550   CONTINUE
  560   CONTINUE
 
      ELSEIF(ISUB.EQ.29) THEN
C...f + g -> f + gamma (q + g -> q + gamma only).
        FGQ=COMFAC*FACA*AS*AEM*1./3.*(SH2+UH2)/(-SH*UH)
        DO 580 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 580
        EI=KCHG(IABS(I),1)/3.
        FACGQ=FGQ*EI**2
        DO 570 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 570
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 570
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGQ
  570   CONTINUE
  580   CONTINUE
 
      ELSEIF(ISUB.EQ.30) THEN
C...f + g -> f + (gamma*/Z0) (q + g -> q + (gamma*/Z0) only).
        FZQ=COMFAC*FACA*AS*AEM*(1./3.)*(SH2+UH2+2.*SQM4*TH)/(-SH*UH)
C...gamma, gamma/Z interference and Z couplings to final fermion pairs.
        HFGG=0.
        HFGZ=0.
        HFZZ=0.
        HBW4=0.
        RADC4=1.+ULALPS(SQM4)/PARU(1)
        DO 590 I=1,MIN(16,MDCY(23,3))
        IDC=I+MDCY(23,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 590
        IMDM=0
        IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2.OR.MDME(IDC,1).EQ.4)
     &  IMDM=1
        IF(I.LE.8) THEN
          EF=KCHG(I,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
        ELSEIF(I.LE.16) THEN
          EF=KCHG(I+2,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
        ENDIF
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
        IF(4.*RM1.LT.1.) THEN
          FCOF=1.
          IF(I.LE.8) FCOF=3.*RADC4
          BE34=SQRT(MAX(0.,1.-4.*RM1))
          IF(IMDM.EQ.1) THEN
            HFGG=HFGG+FCOF*EF**2*(1.+2.*RM1)*BE34
            HFGZ=HFGZ+FCOF*EF*VF*(1.+2.*RM1)*BE34
            HFZZ=HFZZ+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
          ENDIF
          HBW4=HBW4+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
        ENDIF
  590   CONTINUE
C...Propagators: as simulated in PYOFSH and as desired.
        GMMZ=PMAS(23,1)*PMAS(23,2)
        HBW4=HBW4*XWC*SQMZ/((SQM4-SQMZ)**2+GMMZ**2)
        MINT(15)=1
        MINT(61)=1
        CALL PYWIDT(23,SQM4,WDTP,WDTE)
        HFGG=HFGG*VINT(111)/SQM4
        HFGZ=HFGZ*VINT(112)/SQM4
        HFZZ=HFZZ*VINT(114)/SQM4
C...Loop over flavours; consider full gamma/Z structure.
        DO 610 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 610
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        FACZQ=FZQ*(EI**2*HFGG+EI*VI*HFGZ+
     &  (VI**2+AI**2)*HFZZ)/HBW4
        DO 600 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 600
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 600
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZQ
  600   CONTINUE
  610   CONTINUE
      ENDIF
 
      ELSEIF(ISUB.LE.40) THEN
      IF(ISUB.EQ.31) THEN
C...f + g -> f' + W+/- (q + g -> q' + W+/- only).
        FACWQ=COMFAC*FACA*AS*AEM/XW*1./12.*
     &  (SH2+UH2+2.*SQM4*TH)/(-SH*UH)
C...Propagators: as simulated in PYOFSH and as desired.
        GMMW=PMAS(24,1)*PMAS(24,2)
        HBW4=GMMW/((SQM4-SQMW)**2+GMMW**2)
        CALL PYWIDT(24,SQM4,WDTP,WDTE)
        AEMC=ULALEM(SQM4)
        IF(MSTP(8).GE.1) AEMC=AEM
        GMMWC=SQM4*WDTP(0)*AEMC/(24.*XW)
        HBW4C=GMMWC/((SQM4-SQMW)**2+GMMWC**2)
        FACWQ=FACWQ*HBW4C/HBW4
        DO 630 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 630
        IA=IABS(I)
        KCHW=ISIGN(1,KCHG(IA,1)*ISIGN(1,I))
        WIDSC=(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))/WDTP(0)
        DO 620 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 620
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 620
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWQ*VINT(180+I)*WIDSC
  620   CONTINUE
  630   CONTINUE
 
      ELSEIF(ISUB.EQ.32) THEN
C...f + g -> f + H0 (q + g -> q + H0 only).
 
      ELSEIF(ISUB.EQ.33) THEN
C...f + gamma -> f + g (q + gamma -> q + g only).
        FGQ=COMFAC*AS*AEM*8./3.*(SH2+UH2)/(-SH*UH)
        DO 650 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 650
        EI=KCHG(IABS(I),1)/3.
        FACGQ=FGQ*EI**2
        DO 640 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 640
        IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 640
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGQ
  640   CONTINUE
  650   CONTINUE
 
      ELSEIF(ISUB.EQ.34) THEN
C...f + gamma -> f + gamma.
        FGQ=COMFAC*AEM**2*2.*(SH2+UH2)/(-SH*UH)
        DO 670 I=MMINA,MMAXA
        IF(I.EQ.0) GOTO 670
        EI=KCHG(IABS(I),1)/3.
        FACGQ=FGQ*EI**4
        DO 660 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 660
        IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 660
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGQ
  660   CONTINUE
  670   CONTINUE
 
      ELSEIF(ISUB.EQ.35) THEN
C...f + gamma -> f + (gamma*/Z0).
        FZQN=COMFAC*2.*AEM**2*(SH2+UH2+2.*SQM4*TH)
        FZQD=SQPTH*SQM4-SH*UH
C...gamma, gamma/Z interference and Z couplings to final fermion pairs.
        HFGG=0.
        HFGZ=0.
        HFZZ=0.
        HBW4=0.
        RADC4=1.+ULALPS(SQM4)/PARU(1)
        DO 680 I=1,MIN(16,MDCY(23,3))
        IDC=I+MDCY(23,2)-1
        IF(MDME(IDC,1).LT.0) GOTO 680
        IMDM=0
        IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2.OR.MDME(IDC,1).EQ.4)
     &  IMDM=1
        IF(I.LE.8) THEN
          EF=KCHG(I,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
        ELSEIF(I.LE.16) THEN
          EF=KCHG(I+2,1)/3.
          AF=SIGN(1.,EF+0.1)
          VF=AF-4.*EF*XWV
        ENDIF
        RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
        IF(4.*RM1.LT.1.) THEN
          FCOF=1.
          IF(I.LE.8) FCOF=3.*RADC4
          BE34=SQRT(MAX(0.,1.-4.*RM1))
          IF(IMDM.EQ.1) THEN
            HFGG=HFGG+FCOF*EF**2*(1.+2.*RM1)*BE34
            HFGZ=HFGZ+FCOF*EF*VF*(1.+2.*RM1)*BE34
            HFZZ=HFZZ+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
          ENDIF
          HBW4=HBW4+FCOF*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*BE34
        ENDIF
  680   CONTINUE
C...Propagators: as simulated in PYOFSH and as desired.
        GMMZ=PMAS(23,1)*PMAS(23,2)
        HBW4=HBW4*XWC*SQMZ/((SQM4-SQMZ)**2+GMMZ**2)
        MINT(15)=1
        MINT(61)=1
        CALL PYWIDT(23,SQM4,WDTP,WDTE)
        HFGG=HFGG*VINT(111)/SQM4
        HFGZ=HFGZ*VINT(112)/SQM4
        HFZZ=HFZZ*VINT(114)/SQM4
C...Loop over flavours; consider full gamma/Z structure.
        DO 700 I=MMINA,MMAXA
        IF(I.EQ.0) GOTO 700
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        FACZQ=EI**2*(EI**2*HFGG+EI*VI*HFGZ+
     &  (VI**2+AI**2)*HFZZ)/HBW4
        DO 690 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 690
        IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 690
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZQ*FZQN/MAX(PMAS(IABS(I),1)**2*SQM4,FZQD)
  690   CONTINUE
  700   CONTINUE
 
      ELSEIF(ISUB.EQ.36) THEN
C...f + gamma -> f' + W+/-.
        FWQ=COMFAC*AEM**2/(2.*XW)*
     &  (SH2+UH2+2.*SQM4*TH)/(SQPTH*SQM4-SH*UH)
C...Propagators: as simulated in PYOFSH and as desired.
        GMMW=PMAS(24,1)*PMAS(24,2)
        HBW4=GMMW/((SQM4-SQMW)**2+GMMW**2)
        CALL PYWIDT(24,SQM4,WDTP,WDTE)
        AEMC=ULALEM(SQM4)
        IF(MSTP(8).GE.1) AEMC=AEM
        GMMWC=SQM4*WDTP(0)*AEMC/(24.*XW)
        HBW4C=GMMWC/((SQM4-SQMW)**2+GMMWC**2)
        FWQ=FWQ*HBW4C/HBW4
        DO 720 I=MMINA,MMAXA
        IF(I.EQ.0) GOTO 720
        IA=IABS(I)
        EIA=ABS(KCHG(IABS(I),1)/3.)
        FACWQ=FWQ*(EIA-SH/(SH+UH))**2
        KCHW=ISIGN(1,KCHG(IA,1)*ISIGN(1,I))
        WIDSC=(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))/WDTP(0)
        DO 710 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 710
        IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 710
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWQ*VINT(180+I)*WIDSC
  710   CONTINUE
  720   CONTINUE
 
      ELSEIF(ISUB.EQ.37) THEN
C...f + gamma -> f + H0.
 
      ELSEIF(ISUB.EQ.38) THEN
C...f + Z0 -> f + g (q + Z0 -> q + g only).
 
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
C...f + W+/- -> f' + g (q + W+/- -> q' + g only).
 
      ELSEIF(ISUB.EQ.44) THEN
C...f + W+/- -> f' + gamma.
 
      ELSEIF(ISUB.EQ.45) THEN
C...f + W+/- -> f' + Z0.
 
      ELSEIF(ISUB.EQ.46) THEN
C...f + W+/- -> f' + W+/-.
 
      ELSEIF(ISUB.EQ.47) THEN
C...f + W+/- -> f' + H0.
 
      ELSEIF(ISUB.EQ.48) THEN
C...f + H0 -> f + g (q + H0 -> q + g only).
 
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
C...g + g -> f + f~ (g + g -> q + q~ only).
        CALL PYWIDT(21,SH,WDTP,WDTE)
        FACQQ1=COMFAC*AS**2*1./6.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)*
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA
        FACQQ2=COMFAC*AS**2*1./6.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)*
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 730
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQQ1
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=2
        SIGH(NCHN)=FACQQ2
  730   CONTINUE
 
      ELSEIF(ISUB.EQ.54) THEN
C...g + gamma -> f + f~ (g + gamma -> q + q~ only).
        CALL PYWIDT(21,SH,WDTP,WDTE)
        WDTESU=0.
        DO 740 I=1,MIN(8,MDCY(21,3))
        EF=KCHG(I,1)/3.
        WDTESU=WDTESU+EF**2*(WDTE(I,1)+WDTE(I,2)+WDTE(I,3)+WDTE(I,4))
  740   CONTINUE
        FACQQ=COMFAC*AEM*AS*WDTESU*(TH2+UH2)/(TH*UH)
        IF(KFAC(1,21)*KFAC(2,22).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=22
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ
        ENDIF
        IF(KFAC(1,22)*KFAC(2,21).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=22
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ
        ENDIF
 
      ELSEIF(ISUB.EQ.55) THEN
C...g + Z -> f + f~ (g + Z -> q + q~ only).
 
      ELSEIF(ISUB.EQ.56) THEN
C...g + W -> f + f'~ (g + W -> q + q'~ only).
 
      ELSEIF(ISUB.EQ.57) THEN
C...g + H0 -> f + f~ (g + H0 -> q + q~ only).
 
      ELSEIF(ISUB.EQ.58) THEN
C...gamma + gamma -> f + f~.
        CALL PYWIDT(22,SH,WDTP,WDTE)
        WDTESU=0.
        DO 750 I=1,MIN(12,MDCY(22,3))
        IF(I.LE.8) EF= KCHG(I,1)/3.
        IF(I.GE.9) EF= KCHG(9+2*(I-8),1)/3.
        WDTESU=WDTESU+EF**2*(WDTE(I,1)+WDTE(I,2)+WDTE(I,3)+WDTE(I,4))
  750   CONTINUE
        FACFF=COMFAC*AEM**2*WDTESU*2.*(TH2+UH2)/(TH*UH)
        IF(KFAC(1,22)*KFAC(2,22).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=22
          ISIG(NCHN,2)=22
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACFF
        ENDIF
 
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
C...g + g -> g + g.
        FACGG1=COMFAC*AS**2*9./4.*(SH2/TH2+2.*SH/TH+3.+2.*TH/SH+
     &  TH2/SH2)*FACA
        FACGG2=COMFAC*AS**2*9./4.*(UH2/SH2+2.*UH/SH+3.+2.*SH/UH+
     &  SH2/UH2)*FACA
        FACGG3=COMFAC*AS**2*9./4.*(TH2/UH2+2.*TH/UH+3.+2.*UH/TH+
     &  UH2/TH2)
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 760
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACGG1
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=2
        SIGH(NCHN)=0.5*FACGG2
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=3
        SIGH(NCHN)=0.5*FACGG3
  760   CONTINUE
 
      ELSEIF(ISUB.EQ.69) THEN
C...gamma + gamma -> W+ + W-.
        SQMWE=MAX(0.5*SQMW,SQRT(SQM3*SQM4))
        FPROP=SH2/((SQMWE-TH)*(SQMWE-UH))
        FACWW=COMFAC*6.*AEM**2*(1.-FPROP*(4./3.+2.*SQMWE/SH)+
     &  FPROP**2*(2./3.+2.*(SQMWE/SH)**2))*WIDS(24,1)
        IF(KFAC(1,22)*KFAC(2,22).EQ.0) GOTO 770
        NCHN=NCHN+1
        ISIG(NCHN,1)=22
        ISIG(NCHN,2)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWW
  770   CONTINUE
 
      ELSEIF(ISUB.EQ.70) THEN
C...gamma + W+/- -> Z0 + W+/-.
        SQMWE=MAX(0.5*SQMW,SQRT(SQM3*SQM4))
        FPROP=(TH-SQMWE)**2/(-SH*(SQMWE-UH))
        FACZW=COMFAC*6.*AEM**2*(XW1/XW)*
     &  (1.-FPROP*(4./3.+2.*SQMWE/(TH-SQMWE))+
     &  FPROP**2*(2./3.+2.*(SQMWE/(TH-SQMWE))**2))*WIDS(23,2)
        DO 790 KCHW=1,-1,-2
        DO 780 ISDE=1,2
        IF(KFAC(ISDE,22)*KFAC(3-ISDE,24*KCHW).EQ.0) GOTO 780
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=22
        ISIG(NCHN,3-ISDE)=24*KCHW
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZW*WIDS(24,(5-KCHW)/2)
  780   CONTINUE
  790   CONTINUE
      ENDIF
 
      ELSEIF(ISUB.LE.80) THEN
      IF(ISUB.EQ.71) THEN
C...Z0 + Z0 -> Z0 + Z0.
        IF(SH.LE.4.01*SQMZ) GOTO 820
 
        IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons.
          BE2=1.-4.*SQMZ/SH
          TH=-0.5*SH*BE2*(1.-CTH)
          UH=-0.5*SH*BE2*(1.+CTH)
          IF(MAX(TH,UH).GT.-1.) GOTO 820
          SHANG=1./XW1*SQMW/SQMZ*(1.+BE2)**2
          ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
          ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
          THANG=1./XW1*SQMW/SQMZ*(BE2-CTH)**2
          ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG
          ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG
          UHANG=1./XW1*SQMW/SQMZ*(BE2+CTH)**2
          AUHRE=(UH-SQMH)/((UH-SQMH)**2+GMMH**2)*UHANG
          AUHIM=-GMMH/((UH-SQMH)**2+GMMH**2)*UHANG
          FACZZ=COMFAC*1./(4096.*PARU(1)**2*16.*XW1**2)*
     &    (AEM/XW)**4*(SH/SQMW)**2*(SQMZ/SQMW)*SH2
          IF(MSTP(46).LE.0) FACZZ=FACZZ*(ASHRE**2+ASHIM**2)
          IF(MSTP(46).EQ.1) FACZZ=FACZZ*((ASHRE+ATHRE+AUHRE)**2+
     &    (ASHIM+ATHIM+AUHIM)**2)
          IF(MSTP(46).EQ.2) FACZZ=0.
 
        ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
          FACZZ=COMFAC*(AEM/(16.*PARU(1)*XW*XW1))**2*(64./9.)*
     &    ABS(A00U+2.*A20U)**2
        ENDIF
        FACZZ=FACZZ*WIDS(23,1)
 
        DO 810 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 810
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        AVI=AI**2+VI**2
        DO 800 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 800
        EJ=KCHG(IABS(J),1)/3.
        AJ=SIGN(1.,EJ)
        VJ=AJ-4.*EJ*XWV
        AVJ=AJ**2+VJ**2
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACZZ*AVI*AVJ
  800   CONTINUE
  810   CONTINUE
  820   CONTINUE
 
      ELSEIF(ISUB.EQ.72) THEN
C...Z0 + Z0 -> W+ + W-.
        IF(SH.LE.4.01*SQMZ) GOTO 850
 
        IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons.
          BE2=SQRT((1.-4.*SQMW/SH)*(1.-4.*SQMZ/SH))
          CTH2=CTH**2
          TH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH-BE2*CTH)
          UH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH+BE2*CTH)
          IF(MAX(TH,UH).GT.-1.) GOTO 850
          SHANG=4.*SQRT(SQMW/(SQMZ*XW1))*(1.-2.*SQMW/SH)*
     &    (1.-2.*SQMZ/SH)
          ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
          ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
          ATWRE=XW1/SQMZ*SH/(TH-SQMW)*((CTH-BE2)**2*(3./2.+BE2/2.*
     &    CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/
     &    SH*(1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+
     &    4.*(SQMW**2+SQMZ**2)/SH2*CTH2+2.*(SQMW+SQMZ)/SH*BE2*CTH))
          ATWIM=0.
          AUWRE=XW1/SQMZ*SH/(UH-SQMW)*((CTH+BE2)**2*(3./2.-BE2/2.*
     &    CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/
     &    SH*(1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+
     &    4.*(SQMW**2+SQMZ**2)/SH2*CTH2-2.*(SQMW+SQMZ)/SH*BE2*CTH))
          AUWIM=0.
          A4RE=2.*XW1/SQMZ*(3.-CTH2-4.*(SQMW+SQMZ)/SH)
          A4IM=0.
          FACWW=COMFAC*1./(4096.*PARU(1)**2*16.*XW1**2)*
     &    (AEM/XW)**4*(SH/SQMW)**2*(SQMZ/SQMW)*SH2
          IF(MSTP(46).LE.0) FACWW=FACWW*(ASHRE**2+ASHIM**2)
          IF(MSTP(46).EQ.1) FACWW=FACWW*((ASHRE+ATWRE+AUWRE+A4RE)**2+
     &    (ASHIM+ATWIM+AUWIM+A4IM)**2)
          IF(MSTP(46).EQ.2) FACWW=FACWW*((ATWRE+AUWRE+A4RE)**2+
     &    (ATWIM+AUWIM+A4IM)**2)
 
        ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
          FACWW=COMFAC*(AEM/(16.*PARU(1)*XW*XW1))**2*(64./9.)*
     &    ABS(A00U-A20U)**2
        ENDIF
        FACWW=FACWW*WIDS(24,1)
 
        DO 840 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 840
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        AVI=AI**2+VI**2
        DO 830 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 830
        EJ=KCHG(IABS(J),1)/3.
        AJ=SIGN(1.,EJ)
        VJ=AJ-4.*EJ*XWV
        AVJ=AJ**2+VJ**2
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWW*AVI*AVJ
  830   CONTINUE
  840   CONTINUE
  850   CONTINUE
 
      ELSEIF(ISUB.EQ.73) THEN
C...Z0 + W+/- -> Z0 + W+/-.
        IF(SH.LE.2.*SQMZ+2.*SQMW) GOTO 880
 
        IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons.
          BE2=1.-2.*(SQMZ+SQMW)/SH+((SQMZ-SQMW)/SH)**2
          EP1=1.-(SQMZ-SQMW)/SH
          EP2=1.+(SQMZ-SQMW)/SH
          TH=-0.5*SH*BE2*(1.-CTH)
          UH=(SQMZ-SQMW)**2/SH-0.5*SH*BE2*(1.+CTH)
          IF(MAX(TH,UH).GT.-1.) GOTO 880
          THANG=(BE2-EP1*CTH)*(BE2-EP2*CTH)
          ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG
          ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG
          ASWRE=-XW1/SQMZ*SH/(SH-SQMW)*(-BE2*(EP1+EP2)**4*CTH+
     &    1./4.*(BE2+EP1*EP2)**2*((EP1-EP2)**2-4.*BE2*CTH)+
     &    2.*BE2*(BE2+EP1*EP2)*(EP1+EP2)**2*CTH-
     &    1./16.*SH/SQMW*(EP1**2-EP2**2)**2*(BE2+EP1*EP2)**2)
          ASWIM=0.
          AUWRE=XW1/SQMZ*SH/(UH-SQMW)*(-BE2*(EP2+EP1*CTH)*
     &    (EP1+EP2*CTH)*(BE2+EP1*EP2)+BE2*(EP2+EP1*CTH)*
     &    (BE2+EP1*EP2*CTH)*(2.*EP2-EP2*CTH+EP1)-BE2*(EP2+EP1*CTH)**2*
     &    (BE2-EP2**2*CTH)-1./8.*(BE2+EP1*EP2*CTH)**2*((EP1+EP2)**2+
     &    2.*BE2*(1.-CTH))+1./32.*SH/SQMW*(BE2+EP1*EP2*CTH)**2*
     &    (EP1**2-EP2**2)**2-BE2*(EP1+EP2*CTH)*(EP2+EP1*CTH)*
     &    (BE2+EP1*EP2)+BE2*(EP1+EP2*CTH)*(BE2+EP1*EP2*CTH)*
     &    (2.*EP1-EP1*CTH+EP2)-BE2*(EP1+EP2*CTH)**2*(BE2-EP1**2*CTH)-
     &    1./8.*(BE2+EP1*EP2*CTH)**2*((EP1+EP2)**2+2.*BE2*(1.-CTH))+
     &    1./32.*SH/SQMW*(BE2+EP1*EP2*CTH)**2*(EP1**2-EP2**2)**2)
          AUWIM=0.
          A4RE=XW1/SQMZ*(EP1**2*EP2**2*(CTH**2-1.)-
     &    2.*BE2*(EP1**2+EP2**2+EP1*EP2)*CTH-2.*BE2*EP1*EP2)
          A4IM=0.
          FACZW=COMFAC*1./(4096.*PARU(1)**2*4.*XW1)*(AEM/XW)**4*
     &    (SH/SQMW)**2*SQRT(SQMZ/SQMW)*SH2
          IF(MSTP(46).LE.0) FACZW=0.
          IF(MSTP(46).EQ.1) FACZW=FACZW*((ATHRE+ASWRE+AUWRE+A4RE)**2+
     &    (ATHIM+ASWIM+AUWIM+A4IM)**2)
          IF(MSTP(46).EQ.2) FACZW=FACZW*((ASWRE+AUWRE+A4RE)**2+
     &    (ASWIM+AUWIM+A4IM)**2)
 
        ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
          FACZW=COMFAC*AEM**2/(64.*PARU(1)**2*XW**2*XW1)*16.*
     &    ABS(A20U+3.*A11U*CTH)**2
        ENDIF
        FACZW=FACZW*WIDS(23,2)
 
        DO 870 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 870
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        AVI=AI**2+VI**2
        KCHWI=ISIGN(1,KCHG(IABS(I),1)*ISIGN(1,I))
        DO 860 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 860
        EJ=KCHG(IABS(J),1)/3.
        AJ=SIGN(1.,EJ)
        VJ=AI-4.*EJ*XWV
        AVJ=AJ**2+VJ**2
        KCHWJ=ISIGN(1,KCHG(IABS(J),1)*ISIGN(1,J))
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZW*AVI*VINT(180+J)*WIDS(24,(5-KCHWJ)/2)
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=2
        SIGH(NCHN)=FACZW*VINT(180+I)*WIDS(24,(5-KCHWI)/2)*AVJ
  860   CONTINUE
  870   CONTINUE
  880   CONTINUE
 
      ELSEIF(ISUB.EQ.75) THEN
C...W+ + W- -> gamma + gamma.
 
      ELSEIF(ISUB.EQ.76) THEN
C...W+ + W- -> Z0 + Z0.
        IF(SH.LE.4.01*SQMZ) GOTO 910
 
        IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons.
          BE2=SQRT((1.-4.*SQMW/SH)*(1.-4.*SQMZ/SH))
          CTH2=CTH**2
          TH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH-BE2*CTH)
          UH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH+BE2*CTH)
          IF(MAX(TH,UH).GT.-1.) GOTO 910
          SHANG=4.*SQRT(SQMW/(SQMZ*XW1))*(1.-2.*SQMW/SH)*
     &    (1.-2.*SQMZ/SH)
          ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
          ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
          ATWRE=XW1/SQMZ*SH/(TH-SQMW)*((CTH-BE2)**2*(3./2.+BE2/2.*
     &    CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/
     &    SH*(1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+
     &    4.*(SQMW**2+SQMZ**2)/SH2*CTH2+2.*(SQMW+SQMZ)/SH*BE2*CTH))
          ATWIM=0.
          AUWRE=XW1/SQMZ*SH/(UH-SQMW)*((CTH+BE2)**2*(3./2.-BE2/2.*
     &    CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/
     &    SH*(1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+
     &    4.*(SQMW**2+SQMZ**2)/SH2*CTH2-2.*(SQMW+SQMZ)/SH*BE2*CTH))
          AUWIM=0.
          A4RE=2.*XW1/SQMZ*(3.-CTH2-4.*(SQMW+SQMZ)/SH)
          A4IM=0.
          FACZZ=COMFAC*1./(4096.*PARU(1)**2)*(AEM/XW)**4*
     &    (SH/SQMW)**2*SH2
          IF(MSTP(46).LE.0) FACZZ=FACZZ*(ASHRE**2+ASHIM**2)
          IF(MSTP(46).EQ.1) FACZZ=FACZZ*((ASHRE+ATWRE+AUWRE+A4RE)**2+
     &    (ASHIM+ATWIM+AUWIM+A4IM)**2)
          IF(MSTP(46).EQ.2) FACZZ=FACZZ*((ATWRE+AUWRE+A4RE)**2+
     &    (ATWIM+AUWIM+A4IM)**2)
 
        ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
          FACZZ=COMFAC*(AEM/(4.*PARU(1)*XW))**2*(64./9.)*
     &    ABS(A00U-A20U)**2
        ENDIF
        FACZZ=FACZZ*WIDS(23,1)
 
        DO 900 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 900
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)
        DO 890 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 890
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)
        IF(EI*EJ.GT.0.) GOTO 890
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=0.5*FACZZ*VINT(180+I)*VINT(180+J)
  890   CONTINUE
  900   CONTINUE
  910   CONTINUE
 
      ELSEIF(ISUB.EQ.77) THEN
C...W+/- + W+/- -> W+/- + W+/-.
        IF(SH.LE.4.01*SQMW) GOTO 940
 
        IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons.
          BE2=1.-4.*SQMW/SH
          BE4=BE2**2
          CTH2=CTH**2
          CTH3=CTH**3
          TH=-0.5*SH*BE2*(1.-CTH)
          UH=-0.5*SH*BE2*(1.+CTH)
          IF(MAX(TH,UH).GT.-1.) GOTO 940
          SHANG=(1.+BE2)**2
          ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
          ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
          THANG=(BE2-CTH)**2
          ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG
          ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG
          UHANG=(BE2+CTH)**2
          AUHRE=(UH-SQMH)/((UH-SQMH)**2+GMMH**2)*UHANG
          AUHIM=-GMMH/((UH-SQMH)**2+GMMH**2)*UHANG
          SGZANG=1./SQMW*BE2*(3.-BE2)**2*CTH
          ASGRE=XW*SGZANG
          ASGIM=0.
          ASZRE=XW1*SH/(SH-SQMZ)*SGZANG
          ASZIM=0.
          TGZANG=1./SQMW*(BE2*(4.-2.*BE2+BE4)+BE2*(4.-10.*BE2+BE4)*CTH+
     &    (2.-11.*BE2+10.*BE4)*CTH2+BE2*CTH3)
          ATGRE=0.5*XW*SH/TH*TGZANG
          ATGIM=0.
          ATZRE=0.5*XW1*SH/(TH-SQMZ)*TGZANG
          ATZIM=0.
          UGZANG=1./SQMW*(BE2*(4.-2.*BE2+BE4)-BE2*(4.-10.*BE2+BE4)*CTH+
     &    (2.-11.*BE2+10.*BE4)*CTH2-BE2*CTH3)
          AUGRE=0.5*XW*SH/UH*UGZANG
          AUGIM=0.
          AUZRE=0.5*XW1*SH/(UH-SQMZ)*UGZANG
          AUZIM=0.
          A4ARE=1./SQMW*(1.+2.*BE2-6.*BE2*CTH-CTH2)
          A4AIM=0.
          A4SRE=2./SQMW*(1.+2.*BE2-CTH2)
          A4SIM=0.
          FWW=COMFAC*1./(4096.*PARU(1)**2)*(AEM/XW)**4*(SH/SQMW)**2*SH2
          IF(MSTP(46).LE.0) THEN
            AWWARE=ASHRE
            AWWAIM=ASHIM
            AWWSRE=0.
            AWWSIM=0.
          ELSEIF(MSTP(46).EQ.1) THEN
            AWWARE=ASHRE+ATHRE+ASGRE+ASZRE+ATGRE+ATZRE+A4ARE
            AWWAIM=ASHIM+ATHIM+ASGIM+ASZIM+ATGIM+ATZIM+A4AIM
            AWWSRE=-ATHRE-AUHRE+ATGRE+ATZRE+AUGRE+AUZRE+A4SRE
            AWWSIM=-ATHIM-AUHIM+ATGIM+ATZIM+AUGIM+AUZIM+A4SIM
          ELSE
            AWWARE=ASGRE+ASZRE+ATGRE+ATZRE+A4ARE
            AWWAIM=ASGIM+ASZIM+ATGIM+ATZIM+A4AIM
            AWWSRE=ATGRE+ATZRE+AUGRE+AUZRE+A4SRE
            AWWSIM=ATGIM+ATZIM+AUGIM+AUZIM+A4SIM
          ENDIF
          AWWA2=AWWARE**2+AWWAIM**2
          AWWS2=AWWSRE**2+AWWSIM**2
 
        ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron.
          FWWA=COMFAC*(AEM/(4.*PARU(1)*XW))**2*(64./9.)*
     &    ABS(A00U+0.5*A20U+4.5*A11U*CTH)**2
          FWWS=COMFAC*(AEM/(4.*PARU(1)*XW))**2*64.*ABS(A20U)**2
        ENDIF
 
        DO 930 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 930
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)
        DO 920 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 920
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)
        IF(EI*EJ.LT.0.) THEN
C...W+W-
          IF(MSTP(45).EQ.1) GOTO 920
          IF(MSTP(46).LE.2) FACWW=FWW*AWWA2*WIDS(24,1)
          IF(MSTP(46).GE.3) FACWW=FWWA*WIDS(24,1)
        ELSE
C...W+W+/W-W-
          IF(MSTP(45).EQ.2) GOTO 920
          IF(MSTP(46).LE.2) FACWW=FWW*AWWS2
          IF(MSTP(46).GE.3) FACWW=FWWS
          IF(EI.GT.0.) FACWW=FACWW*VINT(91)
          IF(EI.LT.0.) FACWW=FACWW*VINT(92)
        ENDIF
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACWW*VINT(180+I)*VINT(180+J)
        IF(EI*EJ.GT.0.) SIGH(NCHN)=0.5*SIGH(NCHN)
  920   CONTINUE
  930   CONTINUE
  940   CONTINUE
 
      ELSEIF(ISUB.EQ.78) THEN
C...W+/- + H0 -> W+/- + H0.
 
      ELSEIF(ISUB.EQ.79) THEN
C...H0 + H0 -> H0 + H0.
 
      ELSEIF(ISUB.EQ.80) THEN
C...q + gamma -> q' + pi+/-.
        FQPI=COMFAC*(2.*AEM/9.)*(-SH/TH)*(1./SH2+1./TH2)
        ASSH=ULALPS(MAX(0.5,0.5*SH))
        Q2FPSH=0.55/LOG(MAX(2.,2.*SH))
        DELSH=UH*SQRT(ASSH*Q2FPSH)
        ASUH=ULALPS(MAX(0.5,-0.5*UH))
        Q2FPUH=0.55/LOG(MAX(2.,-2.*UH))
        DELUH=SH*SQRT(ASUH*Q2FPUH)
        DO 960 I=MAX(-2,MMINA),MIN(2,MMAXA)
        IF(I.EQ.0) GOTO 960
        EI=KCHG(IABS(I),1)/3.
        EJ=SIGN(1.-ABS(EI),EI)
        DO 950 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 950
        IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 950
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FQPI*(EI*DELSH+EJ*DELUH)**2
  950   CONTINUE
  960   CONTINUE
 
      ENDIF
 
C...C: 2 -> 2, tree diagrams with masses.
 
      ELSEIF(ISUB.LE.90) THEN
      IF(ISUB.EQ.81) THEN
C...q + q~ -> Q + Q~.
        FACQQB=COMFAC*AS**2*4./9.*(((TH-SQM3)**2+
     &  (UH-SQM3)**2)/SH2+2.*SQM3/SH)
        IF(MSTP(35).GE.1) FACQQB=FACQQB*PYHFTH(SH,SQM3,0.)
        WID2=1.
        IF(MINT(55).EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,1)
        IF((MINT(55).EQ.7.OR.MINT(55).EQ.8).AND.MSTP(49).GE.1)
     &  WID2=WIDS(MINT(55)+20,1)
        FACQQB=FACQQB*WID2
        DO 970 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 970
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQQB
  970   CONTINUE
 
      ELSEIF(ISUB.EQ.82) THEN
C...g + g -> Q + Q~.
        IF(MSTP(34).EQ.0) THEN
          FACQQ1=COMFAC*FACA*AS**2*(1./6.)*((UH-SQM3)/(TH-SQM3)-
     &    2.*(UH-SQM3)**2/SH2+4.*(SQM3/SH)*(TH*UH-SQM3**2)/
     &    (TH-SQM3)**2)
          FACQQ2=COMFAC*FACA*AS**2*(1./6.)*((TH-SQM3)/(UH-SQM3)-
     &    2.*(TH-SQM3)**2/SH2+4.*(SQM3/SH)*(TH*UH-SQM3**2)/
     &    (UH-SQM3)**2)
        ELSE
          FACQQ1=COMFAC*FACA*AS**2*(1./6.)*((UH-SQM3)/(TH-SQM3)-
     &    2.25*(UH-SQM3)**2/SH2+4.5*(SQM3/SH)*(TH*UH-SQM3**2)/
     &    (TH-SQM3)**2+0.5*SQM3*TH/(TH-SQM3)**2-SQM3**2/(SH*(TH-SQM3)))
          FACQQ2=COMFAC*FACA*AS**2*(1./6.)*((TH-SQM3)/(UH-SQM3)-
     &    2.25*(TH-SQM3)**2/SH2+4.5*(SQM3/SH)*(TH*UH-SQM3**2)/
     &    (UH-SQM3)**2+0.5*SQM3*UH/(UH-SQM3)**2-SQM3**2/(SH*(UH-SQM3)))
        ENDIF
        IF(MSTP(35).GE.1) THEN
          FATRE=PYHFTH(SH,SQM3,2./7.)
          FACQQ1=FACQQ1*FATRE
          FACQQ2=FACQQ2*FATRE
        ENDIF
        WID2=1.
        IF(MINT(55).EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,1)
        IF((MINT(55).EQ.7.OR.MINT(55).EQ.8).AND.MSTP(49).GE.1)
     &  WID2=WIDS(MINT(55)+20,1)
        FACQQ1=FACQQ1*WID2
        FACQQ2=FACQQ2*WID2
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 980
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQQ1
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=2
        SIGH(NCHN)=FACQQ2
  980   CONTINUE
 
      ELSEIF(ISUB.EQ.83) THEN
C...f + q -> f' + Q.
        FACQQS=COMFAC*(0.5*AEM/XW)**2*SH*(SH-SQM3)/(SQMW-TH)**2
        FACQQU=COMFAC*(0.5*AEM/XW)**2*UH*(UH-SQM3)/(SQMW-TH)**2
        DO 1000 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1000
        DO 990 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 990
        IF(I*J.GT.0.AND.MOD(IABS(I+J),2).EQ.0) GOTO 990
        IF(I*J.LT.0.AND.MOD(IABS(I+J),2).EQ.1) GOTO 990
        IF(IABS(I).LT.MINT(55).AND.MOD(IABS(I+MINT(55)),2).EQ.1) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=1
          IF(MOD(MINT(55),2).EQ.0) FACCKM=VCKM(MINT(55)/2,
     &    (IABS(I)+1)/2)*VINT(180+J)
          IF(MOD(MINT(55),2).EQ.1) FACCKM=VCKM(IABS(I)/2,
     &    (MINT(55)+1)/2)*VINT(180+J)
          WID2=1.
          IF(I.GT.0) THEN
            IF(MINT(55).EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,2)
            IF((MINT(55).EQ.7.OR.MINT(55).EQ.8).AND.MSTP(49).GE.1)
     &      WID2=WIDS(MINT(55)+20,2)
          ELSE
            IF(MINT(55).EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,3)
            IF((MINT(55).EQ.7.OR.MINT(55).EQ.8).AND.MSTP(49).GE.1)
     &      WID2=WIDS(MINT(55)+20,3)
          ENDIF
          IF(I*J.GT.0) SIGH(NCHN)=FACQQS*FACCKM*WID2
          IF(I*J.LT.0) SIGH(NCHN)=FACQQU*FACCKM*WID2
        ENDIF
        IF(IABS(J).LT.MINT(55).AND.MOD(IABS(J+MINT(55)),2).EQ.1) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=2
          IF(MOD(MINT(55),2).EQ.0) FACCKM=VCKM(MINT(55)/2,
     &    (IABS(J)+1)/2)*VINT(180+I)
          IF(MOD(MINT(55),2).EQ.1) FACCKM=VCKM(IABS(J)/2,
     &    (MINT(55)+1)/2)*VINT(180+I)
          IF(J.GT.0) THEN
            IF(MINT(55).EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,2)
            IF((MINT(55).EQ.7.OR.MINT(55).EQ.8).AND.MSTP(49).GE.1)
     &      WID2=WIDS(MINT(55)+20,2)
          ELSE
            IF(MINT(55).EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,3)
            IF((MINT(55).EQ.7.OR.MINT(55).EQ.8).AND.MSTP(49).GE.1)
     &      WID2=WIDS(MINT(55)+20,3)
          ENDIF
          IF(I*J.GT.0) SIGH(NCHN)=FACQQS*FACCKM*WID2
          IF(I*J.LT.0) SIGH(NCHN)=FACQQU*FACCKM*WID2
        ENDIF
  990   CONTINUE
 1000   CONTINUE
 
      ELSEIF(ISUB.EQ.84) THEN
C...g + gamma -> Q + Q~.
        FMTU=SQM3/(SQM3-TH)+SQM3/(SQM3-UH)
        FACQQ=COMFAC*AS*AEM*(KCHG(IABS(MINT(55)),1)/3.)**2*
     &  ((SQM3-TH)/(SQM3-UH)+(SQM3-UH)/(SQM3-TH)+4.*FMTU*(1.-FMTU))
        IF(MSTP(35).GE.1) FACQQ=FACQQ*PYHFTH(SH,SQM3,0.)
        WID2=1.
        IF(MINT(55).EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,1)
        IF((MINT(55).EQ.7.OR.MINT(55).EQ.8).AND.MSTP(49).GE.1)
     &  WID2=WIDS(MINT(55)+20,1)
        FACQQ=FACQQ*WID2
        IF(KFAC(1,21)*KFAC(2,22).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=22
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ
        ENDIF
        IF(KFAC(1,22)*KFAC(2,21).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=22
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ
        ENDIF
 
      ELSEIF(ISUB.EQ.85) THEN
C...gamma + gamma -> F + F~ (heavy fermion, quark or lepton).
        FMTU=SQM3/(SQM3-TH)+SQM3/(SQM3-UH)
        FACFF=COMFAC*AEM**2*(KCHG(IABS(MINT(56)),1)/3.)**4*2.*
     &  ((SQM3-TH)/(SQM3-UH)+(SQM3-UH)/(SQM3-TH)+4.*FMTU*(1.-FMTU))
        IF(IABS(MINT(56)).LT.10) FACFF=3.*FACFF
        IF(IABS(MINT(56)).LT.10.AND.MSTP(35).GE.1)
     &  FACFF=FACFF*PYHFTH(SH,SQM3,1.)
        WID2=1.
        IF(MINT(56).EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,1)
        IF((MINT(56).EQ.7.OR.MINT(56).EQ.8).AND.MSTP(49).GE.1)
     &  WID2=WIDS(MINT(56)+20,1)
        IF(MINT(56).EQ.17.AND.MSTP(49).GE.1) WID2=WIDS(29,1)
        FACFF=FACFF*WID2
        IF(KFAC(1,22)*KFAC(2,22).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=22
          ISIG(NCHN,2)=22
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACFF
        ENDIF
 
      ELSEIF(ISUB.EQ.86) THEN
C...g + g -> J/Psi + g.
        FACQQG=COMFAC*AS**3*(5./9.)*PARP(38)*SQRT(SQM3)*
     &  (((SH*(SH-SQM3))**2+(TH*(TH-SQM3))**2+(UH*(UH-SQM3))**2)/
     &  ((TH-SQM3)*(UH-SQM3))**2)/(SH-SQM3)**2
        IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQG
        ENDIF
 
      ELSEIF(ISUB.EQ.87) THEN
C...g + g -> chi_0c + g.
        PGTW=(SH*TH+TH*UH+UH*SH)/SH2
        QGTW=(SH*TH*UH)/SH**3
        RGTW=SQM3/SH
        FACQQG=COMFAC*AS**3*4.*(PARP(39)/SQRT(SQM3))*(1./SH)*
     &  (9.*RGTW**2*PGTW**4*(RGTW**4-2.*RGTW**2*PGTW+PGTW**2)-
     &  6.*RGTW*PGTW**3*QGTW*(2.*RGTW**4-5.*RGTW**2*PGTW+PGTW**2)-
     &  PGTW**2*QGTW**2*(RGTW**4+2.*RGTW**2*PGTW-PGTW**2)+
     &  2.*RGTW*PGTW*QGTW**3*(RGTW**2-PGTW)+6.*RGTW**2*QGTW**4)/
     &  (QGTW*(QGTW-RGTW*PGTW)**4)
        IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQG
        ENDIF
 
      ELSEIF(ISUB.EQ.88) THEN
C...g + g -> chi_1c + g.
        PGTW=(SH*TH+TH*UH+UH*SH)/SH2
        QGTW=(SH*TH*UH)/SH**3
        RGTW=SQM3/SH
        FACQQG=COMFAC*AS**3*12.*(PARP(39)/SQRT(SQM3))*(1./SH)*
     &  PGTW**2*(RGTW*PGTW**2*(RGTW**2-4.*PGTW)+2.*QGTW*(-RGTW**4+
     &  5.*RGTW**2*PGTW+PGTW**2)-15.*RGTW*QGTW**2)/
     &  (QGTW-RGTW*PGTW)**4
        IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQG
        ENDIF
 
      ELSEIF(ISUB.EQ.89) THEN
C...g + g -> chi_2c + g.
        PGTW=(SH*TH+TH*UH+UH*SH)/SH2
        QGTW=(SH*TH*UH)/SH**3
        RGTW=SQM3/SH
        FACQQG=COMFAC*AS**3*4.*(PARP(39)/SQRT(SQM3))*(1./SH)*
     &  (12.*RGTW**2*PGTW**4*(RGTW**4-2.*RGTW**2*PGTW+PGTW**2)-
     &  3.*RGTW*PGTW**3*QGTW*(8.*RGTW**4-RGTW**2*PGTW+4.*PGTW**2)+
     &  2.*PGTW**2*QGTW**2*(-7.*RGTW**4+43.*RGTW**2*PGTW+PGTW**2)+
     &  RGTW*PGTW*QGTW**3*(16.*RGTW**2-61.*PGTW)+12.*RGTW**2*QGTW**4)/
     &  (QGTW*(QGTW-RGTW*PGTW)**4)
        IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQG
        ENDIF
      ENDIF
 
C...D: Mimimum bias processes.
 
      ELSEIF(ISUB.LE.100) THEN
      IF(ISUB.EQ.91) THEN
C...Elastic scattering.
        SIGS=SIGT(0,0,1)
 
      ELSEIF(ISUB.EQ.92) THEN
C...Single diffractive scattering (first side, i.e. XB).
        SIGS=SIGT(0,0,2)
 
      ELSEIF(ISUB.EQ.93) THEN
C...Single diffractive scattering (second side, i.e. AX).
        SIGS=SIGT(0,0,3)
 
      ELSEIF(ISUB.EQ.94) THEN
C...Double diffractive scattering.
        SIGS=SIGT(0,0,4)
 
      ELSEIF(ISUB.EQ.95) THEN
C...Low-pT scattering.
        SIGS=SIGT(0,0,5)
 
      ELSEIF(ISUB.EQ.96) THEN
C...Multiple interactions: sum of QCD processes.
        CALL PYWIDT(21,SH,WDTP,WDTE)
 
C...q + q' -> q + q'.
        FACQQ1=COMFAC*AS**2*4./9.*(SH2+UH2)/TH2
        FACQQB=COMFAC*AS**2*4./9.*((SH2+UH2)/TH2*FACA-
     &  MSTP(34)*2./3.*UH2/(SH*TH))
        FACQQ2=COMFAC*AS**2*4./9.*((SH2+TH2)/UH2-
     &  MSTP(34)*2./3.*SH2/(TH*UH))
        DO 1020 I=-3,3
        IF(I.EQ.0) GOTO 1020
        DO 1010 J=-3,3
        IF(J.EQ.0) GOTO 1010
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=111
        SIGH(NCHN)=FACQQ1
        IF(I.EQ.-J) SIGH(NCHN)=FACQQB
        IF(I.EQ.J) THEN
          SIGH(NCHN)=0.5*SIGH(NCHN)
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=112
          SIGH(NCHN)=0.5*FACQQ2
        ENDIF
 1010   CONTINUE
 1020   CONTINUE
 
C...q + q~ -> q' + q~' or g + g.
        FACQQB=COMFAC*AS**2*4./9.*(TH2+UH2)/SH2*(WDTE(0,1)+WDTE(0,2)+
     &  WDTE(0,3)+WDTE(0,4))
        FACGG1=COMFAC*AS**2*32./27.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)
        FACGG2=COMFAC*AS**2*32./27.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)
        DO 1030 I=-3,3
        IF(I.EQ.0) GOTO 1030
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=121
        SIGH(NCHN)=FACQQB
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=131
        SIGH(NCHN)=0.5*FACGG1
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=132
        SIGH(NCHN)=0.5*FACGG2
 1030   CONTINUE
 
C...q + g -> q + g.
        FACQG1=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*UH2/TH2-UH/SH)*
     &  FACA
        FACQG2=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*SH2/TH2-SH/UH)
        DO 1050 I=-3,3
        IF(I.EQ.0) GOTO 1050
        DO 1040 ISDE=1,2
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=281
        SIGH(NCHN)=FACQG1
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=282
        SIGH(NCHN)=FACQG2
 1040   CONTINUE
 1050   CONTINUE
 
C...g + g -> q + q~ or g + g.
        FACQQ1=COMFAC*AS**2*1./6.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)*
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA
        FACQQ2=COMFAC*AS**2*1./6.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)*
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA
        FACGG1=COMFAC*AS**2*9./4.*(SH2/TH2+2.*SH/TH+3.+2.*TH/SH+
     &  TH2/SH2)*FACA
        FACGG2=COMFAC*AS**2*9./4.*(UH2/SH2+2.*UH/SH+3.+2.*SH/UH+
     &  SH2/UH2)*FACA
        FACGG3=COMFAC*AS**2*9./4.*(TH2/UH2+2.*TH/UH+3+2.*UH/TH+UH2/TH2)
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=531
        SIGH(NCHN)=FACQQ1
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=532
        SIGH(NCHN)=FACQQ2
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=681
        SIGH(NCHN)=0.5*FACGG1
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=682
        SIGH(NCHN)=0.5*FACGG2
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=683
        SIGH(NCHN)=0.5*FACGG3
      ENDIF
 
C...E: 2 -> 1, loop diagrams.
 
      ELSEIF(ISUB.LE.110) THEN
      IF(ISUB.EQ.101) THEN
C...g + g -> gamma*/Z0.
 
      ELSEIF(ISUB.EQ.102) THEN
C...g + g -> H0 (or H'0, or A0).
        CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=4.*COMFAC/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        HI=HP*WDTP(13)/32.
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1060
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
 1060   CONTINUE
 
      ELSEIF(ISUB.EQ.103) THEN
C...gamma + gamma -> H0 (or H'0, or A0).
        CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=4.*COMFAC/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        HI=HP*WDTP(14)*2.
        IF(KFAC(1,22)*KFAC(2,22).EQ.0) GOTO 1070
        NCHN=NCHN+1
        ISIG(NCHN,1)=22
        ISIG(NCHN,2)=22
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
 1070   CONTINUE
 
C...F: 2 -> 2, box diagrams.
 
      ELSEIF(ISUB.EQ.110) THEN
C...f + f~ -> gamma + H0.
        THUH=MAX(TH*UH,SH*CKIN(3)**2)
        FACHG=COMFAC*(3.*AEM**4)/(2.*PARU(1)**2*XW*SQMW)*SH*THUH
        FACHG=FACHG*WIDS(KFHIGG,2)
C...Calculate loop contributions for intermediate gamma* and Z0.
        CIGTOT=CMPLX(0.,0.)
        CIZTOT=CMPLX(0.,0.)
        JMAX=3*MSTP(1)+1
        DO 1080 J=1,JMAX
        IF(J.LE.2*MSTP(1)) THEN
          FNC=1.
          EJ=KCHG(J,1)/3.
          AJ=SIGN(1.,EJ+0.1)
          VJ=AJ-4.*EJ*XWV
          BALP=SQM4/(2.*PMAS(J,1))**2
          BBET=SH/(2.*PMAS(J,1))**2
        ELSEIF(J.LE.3*MSTP(1)) THEN
          FNC=3.
          JL=2*(J-2*MSTP(1))-1
          EJ=KCHG(10+JL,1)/3.
          AJ=SIGN(1.,EJ+0.1)
          VJ=AJ-4.*EJ*XWV
          BALP=SQM4/(2.*PMAS(10+JL,1))**2
          BBET=SH/(2.*PMAS(10+JL,1))**2
        ELSE
          BALP=SQM4/(2.*PMAS(24,1))**2
          BBET=SH/(2.*PMAS(24,1))**2
        ENDIF
        BABI=1./(BALP-BBET)
        IF(BALP.LT.1.) THEN
          F0ALP=CMPLX(ASIN(SQRT(BALP)),0.)
          F1ALP=F0ALP**2
        ELSE
          F0ALP=CMPLX(LOG(SQRT(BALP)+SQRT(BALP-1.)),-0.5*PARU(1))
          F1ALP=-F0ALP**2
        ENDIF
        F2ALP=SQRT(ABS(BALP-1.)/BALP)*F0ALP
        IF(BBET.LT.1.) THEN
          F0BET=CMPLX(ASIN(SQRT(BBET)),0.)
          F1BET=F0BET**2
        ELSE
          F0BET=CMPLX(LOG(SQRT(BBET)+SQRT(BBET-1.)),-0.5*PARU(1))
          F1BET=-F0BET**2
        ENDIF
        F2BET=SQRT(ABS(BBET-1.)/BBET)*F0BET
        IF(J.LE.3*MSTP(1)) THEN
          FIF=0.5*BABI+BABI**2*(0.5*(1.-BALP+BBET)*(F1BET-F1ALP)+
     &    BBET*(F2BET-F2ALP))
          CIGTOT=CIGTOT+FNC*EJ**2*FIF
          CIZTOT=CIZTOT+FNC*EJ*VJ*FIF
        ELSE
          TXW=XW/XW1
          CIGTOT=CIGTOT-0.5*(BABI*(1.5+BALP)+BABI**2*((1.5-3.*BALP+
     &    4.*BBET)*(F1BET-F1ALP)+BBET*(2.*BALP+3.)*(F2BET-F2ALP)))
          CIZTOT=CIZTOT-0.5*BABI*XW1*((5.-TXW+2.*BALP*(1.-TXW))*
     &    (1.+2.*BABI*BBET*(F2BET-F2ALP))+BABI*(4.*BBET*(3.-TXW)-
     &    (2.*BALP-1.)*(5.-TXW))*(F1BET-F1ALP))
        ENDIF
 1080   CONTINUE
        GMMZ=PMAS(23,1)*PMAS(23,2)
        CIGTOT=CIGTOT/SH
        CIZTOT=CIZTOT*XWC/CMPLX(SH-SQMZ,GMMZ)
C...Loop over initial flavours.
        DO 1090 I=MMINA,MMAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 1090
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACHG*FCOI*(ABS(EI*CIGTOT+VI*CIZTOT)**2+
     &  ABS(AI*CIZTOT)**2)
 1090   CONTINUE
 
      ENDIF
 
      ELSEIF(ISUB.LE.120) THEN
      IF(ISUB.EQ.111) THEN
C...f + f~ -> g + H0 (q + q~ -> g + H0 only).
        A5STUR=0.
        A5STUI=0.
        DO 1100 I=1,2*MSTP(1)
        SQMQ=PMAS(I,1)**2
        EPSS=4.*SQMQ/SH
        EPSH=4.*SQMQ/SQMH
        CALL PYWAUX(1,EPSS,W1SR,W1SI)
        CALL PYWAUX(1,EPSH,W1HR,W1HI)
        CALL PYWAUX(2,EPSS,W2SR,W2SI)
        CALL PYWAUX(2,EPSH,W2HR,W2HI)
        A5STUR=A5STUR+EPSH*(1.+SH/(TH+UH)*(W1SR-W1HR)+
     &  (0.25-SQMQ/(TH+UH))*(W2SR-W2HR))
        A5STUI=A5STUI+EPSH*(SH/(TH+UH)*(W1SI-W1HI)+
     &  (0.25-SQMQ/(TH+UH))*(W2SI-W2HI))
 1100   CONTINUE
        FACGH=COMFAC*FACA/(144.*PARU(1)**2)*AEM/XW*AS**3*SQMH/SQMW*
     &  SQMH/SH*(UH**2+TH**2)/(UH+TH)**2*(A5STUR**2+A5STUI**2)
        FACGH=FACGH*WIDS(25,2)
        DO 1110 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 1110
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGH
 1110   CONTINUE
 
      ELSEIF(ISUB.EQ.112) THEN
C...f + g -> f + H0 (q + g -> q + H0 only).
        A5TSUR=0.
        A5TSUI=0.
        DO 1120 I=1,2*MSTP(1)
        SQMQ=PMAS(I,1)**2
        EPST=4.*SQMQ/TH
        EPSH=4.*SQMQ/SQMH
        CALL PYWAUX(1,EPST,W1TR,W1TI)
        CALL PYWAUX(1,EPSH,W1HR,W1HI)
        CALL PYWAUX(2,EPST,W2TR,W2TI)
        CALL PYWAUX(2,EPSH,W2HR,W2HI)
        A5TSUR=A5TSUR+EPSH*(1.+TH/(SH+UH)*(W1TR-W1HR)+
     &  (0.25-SQMQ/(SH+UH))*(W2TR-W2HR))
        A5TSUI=A5TSUI+EPSH*(TH/(SH+UH)*(W1TI-W1HI)+
     &  (0.25-SQMQ/(SH+UH))*(W2TI-W2HI))
 1120   CONTINUE
        FACQH=COMFAC*FACA/(384.*PARU(1)**2)*AEM/XW*AS**3*SQMH/SQMW*
     &  SQMH/(-TH)*(UH**2+SH**2)/(UH+SH)**2*(A5TSUR**2+A5TSUI**2)
        FACQH=FACQH*WIDS(25,2)
        DO 1140 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 1140
        DO 1130 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 1130
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 1130
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQH
 1130   CONTINUE
 1140   CONTINUE
 
      ELSEIF(ISUB.EQ.113) THEN
C...g + g -> g + H0.
        A2STUR=0.
        A2STUI=0.
        A2USTR=0.
        A2USTI=0.
        A2TUSR=0.
        A2TUSI=0.
        A4STUR=0.
        A4STUI=0.
        DO 1150 I=1,2*MSTP(1)
        SQMQ=PMAS(I,1)**2
        EPSS=4.*SQMQ/SH
        EPST=4.*SQMQ/TH
        EPSU=4.*SQMQ/UH
        EPSH=4.*SQMQ/SQMH
        IF(EPSH.LT.1.E-6) GOTO 1150
        CALL PYWAUX(1,EPSS,W1SR,W1SI)
        CALL PYWAUX(1,EPST,W1TR,W1TI)
        CALL PYWAUX(1,EPSU,W1UR,W1UI)
        CALL PYWAUX(1,EPSH,W1HR,W1HI)
        CALL PYWAUX(2,EPSS,W2SR,W2SI)
        CALL PYWAUX(2,EPST,W2TR,W2TI)
        CALL PYWAUX(2,EPSU,W2UR,W2UI)
        CALL PYWAUX(2,EPSH,W2HR,W2HI)
        CALL PYI3AU(EPSS,TH/UH,Y3STUR,Y3STUI)
        CALL PYI3AU(EPSS,UH/TH,Y3SUTR,Y3SUTI)
        CALL PYI3AU(EPST,SH/UH,Y3TSUR,Y3TSUI)
        CALL PYI3AU(EPST,UH/SH,Y3TUSR,Y3TUSI)
        CALL PYI3AU(EPSU,SH/TH,Y3USTR,Y3USTI)
        CALL PYI3AU(EPSU,TH/SH,Y3UTSR,Y3UTSI)
        CALL PYI3AU(EPSH,SQMH/SH*TH/UH,YHSTUR,YHSTUI)
        CALL PYI3AU(EPSH,SQMH/SH*UH/TH,YHSUTR,YHSUTI)
        CALL PYI3AU(EPSH,SQMH/TH*SH/UH,YHTSUR,YHTSUI)
        CALL PYI3AU(EPSH,SQMH/TH*UH/SH,YHTUSR,YHTUSI)
        CALL PYI3AU(EPSH,SQMH/UH*SH/TH,YHUSTR,YHUSTI)
        CALL PYI3AU(EPSH,SQMH/UH*TH/SH,YHUTSR,YHUTSI)
        W3STUR=YHSTUR-Y3STUR-Y3UTSR
        W3STUI=YHSTUI-Y3STUI-Y3UTSI
        W3SUTR=YHSUTR-Y3SUTR-Y3TUSR
        W3SUTI=YHSUTI-Y3SUTI-Y3TUSI
        W3TSUR=YHTSUR-Y3TSUR-Y3USTR
        W3TSUI=YHTSUI-Y3TSUI-Y3USTI
        W3TUSR=YHTUSR-Y3TUSR-Y3SUTR
        W3TUSI=YHTUSI-Y3TUSI-Y3SUTI
        W3USTR=YHUSTR-Y3USTR-Y3TSUR
        W3USTI=YHUSTI-Y3USTI-Y3TSUI
        W3UTSR=YHUTSR-Y3UTSR-Y3STUR
        W3UTSI=YHUTSI-Y3UTSI-Y3STUI
        B2STUR=SQMQ/SQMH**2*(SH*(UH-SH)/(SH+UH)+2.*TH*UH*(UH+2.*SH)/
     &  (SH+UH)**2*(W1TR-W1HR)+(SQMQ-SH/4.)*(0.5*W2SR+0.5*W2HR-W2TR+
     &  W3STUR)+SH2*(2.*SQMQ/(SH+UH)**2-0.5/(SH+UH))*(W2TR-W2HR)+
     &  0.5*TH*UH/SH*(W2HR-2.*W2TR)+0.125*(SH-12.*SQMQ-4.*TH*UH/SH)*
     &  W3TSUR)
        B2STUI=SQMQ/SQMH**2*(2.*TH*UH*(UH+2.*SH)/(SH+UH)**2*
     &  (W1TI-W1HI)+(SQMQ-SH/4.)*(0.5*W2SI+0.5*W2HI-W2TI+W3STUI)+
     &  SH2*(2.*SQMQ/(SH+UH)**2-0.5/(SH+UH))*(W2TI-W2HI)+0.5*TH*UH/SH*
     &  (W2HI-2.*W2TI)+0.125*(SH-12.*SQMQ-4.*TH*UH/SH)*W3TSUI)
        B2SUTR=SQMQ/SQMH**2*(SH*(TH-SH)/(SH+TH)+2.*UH*TH*(TH+2.*SH)/
     &  (SH+TH)**2*(W1UR-W1HR)+(SQMQ-SH/4.)*(0.5*W2SR+0.5*W2HR-W2UR+
     &  W3SUTR)+SH2*(2.*SQMQ/(SH+TH)**2-0.5/(SH+TH))*(W2UR-W2HR)+
     &  0.5*UH*TH/SH*(W2HR-2.*W2UR)+0.125*(SH-12.*SQMQ-4.*UH*TH/SH)*
     &  W3USTR)
        B2SUTI=SQMQ/SQMH**2*(2.*UH*TH*(TH+2.*SH)/(SH+TH)**2*
     &  (W1UI-W1HI)+(SQMQ-SH/4.)*(0.5*W2SI+0.5*W2HI-W2UI+W3SUTI)+
     &  SH2*(2.*SQMQ/(SH+TH)**2-0.5/(SH+TH))*(W2UI-W2HI)+0.5*UH*TH/SH*
     &  (W2HI-2.*W2UI)+0.125*(SH-12.*SQMQ-4.*UH*TH/SH)*W3USTI)
        B2TSUR=SQMQ/SQMH**2*(TH*(UH-TH)/(TH+UH)+2.*SH*UH*(UH+2.*TH)/
     &  (TH+UH)**2*(W1SR-W1HR)+(SQMQ-TH/4.)*(0.5*W2TR+0.5*W2HR-W2SR+
     &  W3TSUR)+TH2*(2.*SQMQ/(TH+UH)**2-0.5/(TH+UH))*(W2SR-W2HR)+
     &  0.5*SH*UH/TH*(W2HR-2.*W2SR)+0.125*(TH-12.*SQMQ-4.*SH*UH/TH)*
     &  W3STUR)
        B2TSUI=SQMQ/SQMH**2*(2.*SH*UH*(UH+2.*TH)/(TH+UH)**2*
     &  (W1SI-W1HI)+(SQMQ-TH/4.)*(0.5*W2TI+0.5*W2HI-W2SI+W3TSUI)+
     &  TH2*(2.*SQMQ/(TH+UH)**2-0.5/(TH+UH))*(W2SI-W2HI)+0.5*SH*UH/TH*
     &  (W2HI-2.*W2SI)+0.125*(TH-12.*SQMQ-4.*SH*UH/TH)*W3STUI)
        B2TUSR=SQMQ/SQMH**2*(TH*(SH-TH)/(TH+SH)+2.*UH*SH*(SH+2.*TH)/
     &  (TH+SH)**2*(W1UR-W1HR)+(SQMQ-TH/4.)*(0.5*W2TR+0.5*W2HR-W2UR+
     &  W3TUSR)+TH2*(2.*SQMQ/(TH+SH)**2-0.5/(TH+SH))*(W2UR-W2HR)+
     &  0.5*UH*SH/TH*(W2HR-2.*W2UR)+0.125*(TH-12.*SQMQ-4.*UH*SH/TH)*
     &  W3UTSR)
        B2TUSI=SQMQ/SQMH**2*(2.*UH*SH*(SH+2.*TH)/(TH+SH)**2*
     &  (W1UI-W1HI)+(SQMQ-TH/4.)*(0.5*W2TI+0.5*W2HI-W2UI+W3TUSI)+
     &  TH2*(2.*SQMQ/(TH+SH)**2-0.5/(TH+SH))*(W2UI-W2HI)+0.5*UH*SH/TH*
     &  (W2HI-2.*W2UI)+0.125*(TH-12.*SQMQ-4.*UH*SH/TH)*W3UTSI)
        B2USTR=SQMQ/SQMH**2*(UH*(TH-UH)/(UH+TH)+2.*SH*TH*(TH+2.*UH)/
     &  (UH+TH)**2*(W1SR-W1HR)+(SQMQ-UH/4.)*(0.5*W2UR+0.5*W2HR-W2SR+
     &  W3USTR)+UH2*(2.*SQMQ/(UH+TH)**2-0.5/(UH+TH))*(W2SR-W2HR)+
     &  0.5*SH*TH/UH*(W2HR-2.*W2SR)+0.125*(UH-12.*SQMQ-4.*SH*TH/UH)*
     &  W3SUTR)
        B2USTI=SQMQ/SQMH**2*(2.*SH*TH*(TH+2.*UH)/(UH+TH)**2*
     &  (W1SI-W1HI)+(SQMQ-UH/4.)*(0.5*W2UI+0.5*W2HI-W2SI+W3USTI)+
     &  UH2*(2.*SQMQ/(UH+TH)**2-0.5/(UH+TH))*(W2SI-W2HI)+0.5*SH*TH/UH*
     &  (W2HI-2.*W2SI)+0.125*(UH-12.*SQMQ-4.*SH*TH/UH)*W3SUTI)
        B2UTSR=SQMQ/SQMH**2*(UH*(SH-UH)/(UH+SH)+2.*TH*SH*(SH+2.*UH)/
     &  (UH+SH)**2*(W1TR-W1HR)+(SQMQ-UH/4.)*(0.5*W2UR+0.5*W2HR-W2TR+
     &  W3UTSR)+UH2*(2.*SQMQ/(UH+SH)**2-0.5/(UH+SH))*(W2TR-W2HR)+
     &  0.5*TH*SH/UH*(W2HR-2.*W2TR)+0.125*(UH-12.*SQMQ-4.*TH*SH/UH)*
     &  W3TUSR)
        B2UTSI=SQMQ/SQMH**2*(2.*TH*SH*(SH+2.*UH)/(UH+SH)**2*
     &  (W1TI-W1HI)+(SQMQ-UH/4.)*(0.5*W2UI+0.5*W2HI-W2TI+W3UTSI)+
     &  UH2*(2.*SQMQ/(UH+SH)**2-0.5/(UH+SH))*(W2TI-W2HI)+0.5*TH*SH/UH*
     &  (W2HI-2.*W2TI)+0.125*(UH-12.*SQMQ-4.*TH*SH/UH)*W3TUSI)
        B4STUR=0.25*EPSH*(-2./3.+0.25*(EPSH-1.)*(W2SR-W2HR+W3STUR))
        B4STUI=0.25*EPSH*0.25*(EPSH-1.)*(W2SI-W2HI+W3STUI)
        B4TUSR=0.25*EPSH*(-2./3.+0.25*(EPSH-1.)*(W2TR-W2HR+W3TUSR))
        B4TUSI=0.25*EPSH*0.25*(EPSH-1.)*(W2TI-W2HI+W3TUSI)
        B4USTR=0.25*EPSH*(-2./3.+0.25*(EPSH-1.)*(W2UR-W2HR+W3USTR))
        B4USTI=0.25*EPSH*0.25*(EPSH-1.)*(W2UI-W2HI+W3USTI)
        A2STUR=A2STUR+B2STUR+B2SUTR
        A2STUI=A2STUI+B2STUI+B2SUTI
        A2USTR=A2USTR+B2USTR+B2UTSR
        A2USTI=A2USTI+B2USTI+B2UTSI
        A2TUSR=A2TUSR+B2TUSR+B2TSUR
        A2TUSI=A2TUSI+B2TUSI+B2TSUI
        A4STUR=A4STUR+B4STUR+B4USTR+B4TUSR
        A4STUI=A4STUI+B4STUI+B4USTI+B4TUSI
 1150   CONTINUE
        FACGH=COMFAC*FACA*3./(128.*PARU(1)**2)*AEM/XW*AS**3*
     &  SQMH/SQMW*SQMH**3/(SH*TH*UH)*(A2STUR**2+A2STUI**2+A2USTR**2+
     &  A2USTI**2+A2TUSR**2+A2TUSI**2+A4STUR**2+A4STUI**2)
        FACGH=FACGH*WIDS(25,2)
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1160
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACGH
 1160   CONTINUE
 
      ELSEIF(ISUB.EQ.114.OR.ISUB.EQ.115) THEN
C...g + g -> gamma + gamma or g + g -> g + gamma.
        A0STUR=0.
        A0STUI=0.
        A0TSUR=0.
        A0TSUI=0.
        A0UTSR=0.
        A0UTSI=0.
        A1STUR=0.
        A1STUI=0.
        A2STUR=0.
        A2STUI=0.
        ALST=LOG(-SH/TH)
        ALSU=LOG(-SH/UH)
        ALTU=LOG(TH/UH)
        IMAX=2*MSTP(1)
        IF(MSTP(38).GE.1.AND.MSTP(38).LE.8) IMAX=MSTP(38)
        DO 1170 I=1,IMAX
        EI=KCHG(IABS(I),1)/3.
        EIWT=EI**2
        IF(ISUB.EQ.115) EIWT=EI
        SQMQ=PMAS(I,1)**2
        EPSS=4.*SQMQ/SH
        EPST=4.*SQMQ/TH
        EPSU=4.*SQMQ/UH
        IF((MSTP(38).GE.1.AND.MSTP(38).LE.8).OR.EPSS.LT.1.E-4) THEN
          B0STUR=1.+(TH-UH)/SH*ALTU+0.5*(TH2+UH2)/SH2*(ALTU**2+
     &    PARU(1)**2)
          B0STUI=0.
          B0TSUR=1.+(SH-UH)/TH*ALSU+0.5*(SH2+UH2)/TH2*ALSU**2
          B0TSUI=-PARU(1)*((SH-UH)/TH+(SH2+UH2)/TH2*ALSU)
          B0UTSR=1.+(SH-TH)/UH*ALST+0.5*(SH2+TH2)/UH2*ALST**2
          B0UTSI=-PARU(1)*((SH-TH)/UH+(SH2+TH2)/UH2*ALST)
          B1STUR=-1.
          B1STUI=0.
          B2STUR=-1.
          B2STUI=0.
        ELSE
          CALL PYWAUX(1,EPSS,W1SR,W1SI)
          CALL PYWAUX(1,EPST,W1TR,W1TI)
          CALL PYWAUX(1,EPSU,W1UR,W1UI)
          CALL PYWAUX(2,EPSS,W2SR,W2SI)
          CALL PYWAUX(2,EPST,W2TR,W2TI)
          CALL PYWAUX(2,EPSU,W2UR,W2UI)
          CALL PYI3AU(EPSS,TH/UH,Y3STUR,Y3STUI)
          CALL PYI3AU(EPSS,UH/TH,Y3SUTR,Y3SUTI)
          CALL PYI3AU(EPST,SH/UH,Y3TSUR,Y3TSUI)
          CALL PYI3AU(EPST,UH/SH,Y3TUSR,Y3TUSI)
          CALL PYI3AU(EPSU,SH/TH,Y3USTR,Y3USTI)
          CALL PYI3AU(EPSU,TH/SH,Y3UTSR,Y3UTSI)
          B0STUR=1.+(1.+2.*TH/SH)*W1TR+(1.+2.*UH/SH)*W1UR+
     &    0.5*((TH2+UH2)/SH2-EPSS)*(W2TR+W2UR)-
     &    0.25*EPST*(1.-0.5*EPSS)*(Y3SUTR+Y3TUSR)-
     &    0.25*EPSU*(1.-0.5*EPSS)*(Y3STUR+Y3UTSR)+
     &    0.25*(-2.*(TH2+UH2)/SH2+4.*EPSS+EPST+EPSU+0.5*EPST*EPSU)*
     &    (Y3TSUR+Y3USTR)
          B0STUI=(1.+2.*TH/SH)*W1TI+(1.+2.*UH/SH)*W1UI+
     &    0.5*((TH2+UH2)/SH2-EPSS)*(W2TI+W2UI)-
     &    0.25*EPST*(1.-0.5*EPSS)*(Y3SUTI+Y3TUSI)-
     &    0.25*EPSU*(1.-0.5*EPSS)*(Y3STUI+Y3UTSI)+
     &    0.25*(-2.*(TH2+UH2)/SH2+4.*EPSS+EPST+EPSU+0.5*EPST*EPSU)*
     &    (Y3TSUI+Y3USTI)
          B0TSUR=1.+(1.+2.*SH/TH)*W1SR+(1.+2.*UH/TH)*W1UR+
     &    0.5*((SH2+UH2)/TH2-EPST)*(W2SR+W2UR)-
     &    0.25*EPSS*(1.-0.5*EPST)*(Y3TUSR+Y3SUTR)-
     &    0.25*EPSU*(1.-0.5*EPST)*(Y3TSUR+Y3USTR)+
     &    0.25*(-2.*(SH2+UH2)/TH2+4.*EPST+EPSS+EPSU+0.5*EPSS*EPSU)*
     &    (Y3STUR+Y3UTSR)
          B0TSUI=(1.+2.*SH/TH)*W1SI+(1.+2.*UH/TH)*W1UI+
     &    0.5*((SH2+UH2)/TH2-EPST)*(W2SI+W2UI)-
     &    0.25*EPSS*(1.-0.5*EPST)*(Y3TUSI+Y3SUTI)-
     &    0.25*EPSU*(1.-0.5*EPST)*(Y3TSUI+Y3USTI)+
     &    0.25*(-2.*(SH2+UH2)/TH2+4.*EPST+EPSS+EPSU+0.5*EPSS*EPSU)*
     &    (Y3STUI+Y3UTSI)
          B0UTSR=1.+(1.+2.*TH/UH)*W1TR+(1.+2.*SH/UH)*W1SR+
     &    0.5*((TH2+SH2)/UH2-EPSU)*(W2TR+W2SR)-
     &    0.25*EPST*(1.-0.5*EPSU)*(Y3USTR+Y3TSUR)-
     &    0.25*EPSS*(1.-0.5*EPSU)*(Y3UTSR+Y3STUR)+
     &    0.25*(-2.*(TH2+SH2)/UH2+4.*EPSU+EPST+EPSS+0.5*EPST*EPSS)*
     &    (Y3TUSR+Y3SUTR)
          B0UTSI=(1.+2.*TH/UH)*W1TI+(1.+2.*SH/UH)*W1SI+
     &    0.5*((TH2+SH2)/UH2-EPSU)*(W2TI+W2SI)-
     &    0.25*EPST*(1.-0.5*EPSU)*(Y3USTI+Y3TSUI)-
     &    0.25*EPSS*(1.-0.5*EPSU)*(Y3UTSI+Y3STUI)+
     &    0.25*(-2.*(TH2+SH2)/UH2+4.*EPSU+EPST+EPSS+0.5*EPST*EPSS)*
     &    (Y3TUSI+Y3SUTI)
          B1STUR=-1.-0.25*(EPSS+EPST+EPSU)*(W2SR+W2TR+W2UR)+
     &    0.25*(EPSU+0.5*EPSS*EPST)*(Y3SUTR+Y3TUSR)+
     &    0.25*(EPST+0.5*EPSS*EPSU)*(Y3STUR+Y3UTSR)+
     &    0.25*(EPSS+0.5*EPST*EPSU)*(Y3TSUR+Y3USTR)
          B1STUI=-0.25*(EPSS+EPST+EPSU)*(W2SI+W2TI+W2UI)+
     &    0.25*(EPSU+0.5*EPSS*EPST)*(Y3SUTI+Y3TUSI)+
     &    0.25*(EPST+0.5*EPSS*EPSU)*(Y3STUI+Y3UTSI)+
     &    0.25*(EPSS+0.5*EPST*EPSU)*(Y3TSUI+Y3USTI)
          B2STUR=-1.+0.125*EPSS*EPST*(Y3SUTR+Y3TUSR)+
     &    0.125*EPSS*EPSU*(Y3STUR+Y3UTSR)+
     &    0.125*EPST*EPSU*(Y3TSUR+Y3USTR)
          B2STUI=0.125*EPSS*EPST*(Y3SUTI+Y3TUSI)+
     &    0.125*EPSS*EPSU*(Y3STUI+Y3UTSI)+
     &    0.125*EPST*EPSU*(Y3TSUI+Y3USTI)
        ENDIF
        A0STUR=A0STUR+EIWT*B0STUR
        A0STUI=A0STUI+EIWT*B0STUI
        A0TSUR=A0TSUR+EIWT*B0TSUR
        A0TSUI=A0TSUI+EIWT*B0TSUI
        A0UTSR=A0UTSR+EIWT*B0UTSR
        A0UTSI=A0UTSI+EIWT*B0UTSI
        A1STUR=A1STUR+EIWT*B1STUR
        A1STUI=A1STUI+EIWT*B1STUI
        A2STUR=A2STUR+EIWT*B2STUR
        A2STUI=A2STUI+EIWT*B2STUI
 1170   CONTINUE
        ASQSUM=A0STUR**2+A0STUI**2+A0TSUR**2+A0TSUI**2+A0UTSR**2+
     &  A0UTSI**2+4.*A1STUR**2+4.*A1STUI**2+A2STUR**2+A2STUI**2
        FACGG=COMFAC*FACA/(16.*PARU(1)**2)*AS**2*AEM**2*ASQSUM
        FACGP=COMFAC*FACA*5./(192.*PARU(1)**2)*AS**3*AEM*ASQSUM
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1180
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        IF(ISUB.EQ.114) SIGH(NCHN)=0.5*FACGG
        IF(ISUB.EQ.115) SIGH(NCHN)=FACGP
 1180   CONTINUE
 
      ELSEIF(ISUB.EQ.116) THEN
C...g + g -> gamma + Z0.
 
      ELSEIF(ISUB.EQ.117) THEN
C...g + g -> Z0 + Z0.
 
      ELSEIF(ISUB.EQ.118) THEN
C...g + g -> W+ + W-.
 
      ENDIF
 
C...G: 2 -> 3, tree diagrams.
 
      ELSEIF(ISUB.LE.140) THEN
      IF(ISUB.EQ.121) THEN
C...g + g -> Q + Q~ + H0.
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1190
        IA=KFPR(ISUBSV,2)
        PMF=PMAS(IA,1)
        FACQQH=COMFAC*(4.*PARU(1)*AEM/XW)*(4.*PARU(1)*AS)**2*
     &  (0.5*PMF/PMAS(24,1))**2
        IF(IA.LE.10.AND.MSTP(37).EQ.1.AND.MSTP(2).GE.1) FACQQH=
     &  FACQQH*(LOG(MAX(4.,PARP(37)**2*PMF**2/PARU(117)**2))/
     &  LOG(MAX(4.,SH/PARU(117)**2)))**(24./(33.-2.*MSTU(118)))
        WID2=1.
        IF(IA.EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,1)
        IF((IA.EQ.7.OR.IA.EQ.8).AND.MSTP(49).GE.1) WID2=WIDS(IA+20,1)
        FACQQH=FACQQH*WID2
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
          IKFI=1
          IF(IA.LE.10.AND.MOD(IA,2).EQ.0) IKFI=2
          IF(IA.GT.10) IKFI=3
          FACQQH=FACQQH*PARU(150+10*IHIGG+IKFI)**2
        ENDIF
        CALL PYQQBH(WTQQBH)
        CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=(1./PARU(1))*VINT(2)*HF/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQQH*WTQQBH*FACBW
 1190   CONTINUE
 
      ELSEIF(ISUB.EQ.122) THEN
C...q + q~ -> Q + Q~ + H0.
        IA=KFPR(ISUBSV,2)
        PMF=PMAS(IA,1)
        FACQQH=COMFAC*(4.*PARU(1)*AEM/XW)*(4.*PARU(1)*AS)**2*
     &  (0.5*PMF/PMAS(24,1))**2
        IF(IA.LE.10.AND.MSTP(37).EQ.1.AND.MSTP(2).GE.1) FACQQH=
     &  FACQQH*(LOG(MAX(4.,PARP(37)**2*PMF**2/PARU(117)**2))/
     &  LOG(MAX(4.,SH/PARU(117)**2)))**(24./(33.-2.*MSTU(118)))
        WID2=1.
        IF(IA.EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,1)
        IF((IA.EQ.7.OR.IA.EQ.8).AND.MSTP(49).GE.1) WID2=WIDS(IA+20,1)
        FACQQH=FACQQH*WID2
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
          IKFI=1
          IF(IA.LE.10.AND.MOD(IA,2).EQ.0) IKFI=2
          IF(IA.GT.10) IKFI=3
          FACQQH=FACQQH*PARU(150+10*IHIGG+IKFI)**2
        ENDIF
        CALL PYQQBH(WTQQBH)
        CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=(1./PARU(1))*VINT(2)*HF/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        DO 1200 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 1200
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACQQH*WTQQBH*FACBW
 1200   CONTINUE
 
      ELSEIF(ISUB.EQ.123) THEN
C...f + f' -> f + f' + H0 (or H'0, or A0) (Z0 + Z0 -> H0 as
C...inner process).
        FACNOR=COMFAC*(4.*PARU(1)*AEM/(XW*XW1))**3*SQMZ/32.
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACNOR=FACNOR*
     &  PARU(154+10*IHIGG)**2
        FACPRP=1./((VINT(215)-VINT(204)**2)*(VINT(216)-VINT(209)**2))**2
        FACZZ1=FACNOR*FACPRP*(0.5*TAUP*VINT(2))*VINT(219)
        FACZZ2=FACNOR*FACPRP*VINT(217)*VINT(218)
        CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=(1./PARU(1))*VINT(2)*HF/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        DO 1220 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1220
        IA=IABS(I)
        DO 1210 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1210
        JA=IABS(J)
        EI=KCHG(IA,1)*ISIGN(1,I)/3.
        AI=SIGN(1.,KCHG(IA,1)+0.5)*ISIGN(1,I)
        VI=AI-4.*EI*XWV
        EJ=KCHG(JA,1)*ISIGN(1,J)/3.
        AJ=SIGN(1.,KCHG(JA,1)+0.5)*ISIGN(1,J)
        VJ=AJ-4.*EJ*XWV
        FACLR1=(VI**2+AI**2)*(VJ**2+AJ**2)+4.*VI*AI*VJ*AJ
        FACLR2=(VI**2+AI**2)*(VJ**2+AJ**2)-4.*VI*AI*VJ*AJ
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=(FACLR1*FACZZ1+FACLR2*FACZZ2)*FACBW
 1210   CONTINUE
 1220   CONTINUE
 
      ELSEIF(ISUB.EQ.124) THEN
C...f + f' -> f" + f"' + H0 (or H'0, or A0) (W+ + W- -> H0 as
C...inner process).
        FACNOR=COMFAC*(4.*PARU(1)*AEM/XW)**3*SQMW
        IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACNOR=FACNOR*
     &  PARU(155+10*IHIGG)**2
        FACPRP=1./((VINT(215)-VINT(204)**2)*(VINT(216)-VINT(209)**2))**2
        FACWW=FACNOR*FACPRP*(0.5*TAUP*VINT(2))*VINT(219)
        CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        FACBW=(1./PARU(1))*VINT(2)*HF/((SH-SQMH)**2+HS**2)
        IF(ABS(SH-SQMH).GT.100.*HS) FACBW=0.
        DO 1240 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1240
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)
        DO 1230 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1230
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)
        IF(EI*EJ.GT.0.) GOTO 1230
        FACLR=VINT(180+I)*VINT(180+J)
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACLR*FACWW*FACBW
 1230   CONTINUE
 1240   CONTINUE
 
      ELSEIF(ISUB.EQ.131) THEN
C...g + g -> Z0 + q + qbar.
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1280
 
C...Read out information on flavours, masses, couplings.
        KFQ=KFPR(131,2)
        KFL=IABS(KFDP(MINT(35),1))
        PMH=SQRT(SH)
        PMQQ=SQRT(VINT(64))
        PMLL=SQRT(VINT(63))
        PMQ=PMAS(KFQ,1)
        QFQ=KCHG(KFQ,1)/3.
        AFQ=SIGN(1.,QFQ+0.1)
        VFQ=AFQ-4.*XWV*QFQ
        QFL=KCHG(KFL,1)/3.
        AFL=SIGN(1.,QFL+0.1)
        VFL=AFL-4.*XWV*QFL
        WID2=1.
        IF(KFQ.EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,1)
        IF((KFQ.EQ.7.OR.KFQ.EQ.8).AND.MSTP(49).GE.1) WID2=WIDS(KFQ+20,1)
 
C...Set line numbers for particles.
        IG1=MINT(84)+1
        IG2=MINT(84)+2
        IQ1=MINT(84)+3
        IQ2=MINT(84)+4
        IL1=MINT(84)+5
        IL2=MINT(84)+6
        IZ=MINT(84)+7
 
C...Reconstruct decay kinematics.
        DO 1260 I=MINT(84)+1,MINT(84)+7
        K(I,1)=1
        DO 1250 J=1,5
        P(I,J)=0.
 1250   CONTINUE
 1260   CONTINUE
        P(IG1,4)=0.5*PMH
        P(IG1,3)=P(IG1,4)
        P(IG2,4)=P(IG1,4)
        P(IG2,3)=-P(IG1,3)
        P(IQ1,5)=PMQ
        P(IQ1,4)=0.5*PMQQ
        P(IQ1,3)=SQRT(MAX(0.,P(IQ1,4)**2-PMQ**2))
        P(IQ2,5)=PMQ
        P(IQ2,4)=P(IQ1,4)
        P(IQ2,3)=-P(IQ1,3)
        P(IL1,4)=0.5*PMLL
        P(IL1,3)=P(IL1,4)
        P(IL2,4)=P(IL1,4)
        P(IL2,3)=-P(IL1,3)
        P(IZ,5)=PMLL
        P(IZ,4)=0.5*(PMH+(PMLL**2-PMQQ**2)/PMH)
        P(IZ,3)=SQRT(MAX(0.,P(IZ,4)**2-PMLL**2))
        CALL LUDBRB(IQ1,IQ2,ACOS(VINT(83)),VINT(84),0D0,0D0,
     &  -DBLE(P(IZ,3)/(PMH-P(IZ,4))))
        CALL LUDBRB(IL1,IL2,ACOS(VINT(81)),VINT(82),0D0,0D0,
     &  DBLE(P(IZ,3)/P(IZ,4)))
        CALL LUDBRB(IQ1,IZ,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0)
 
C...Interface information to program of Ronald Kleiss.
        RKMQ=PMQ
        RKMZ=PMAS(23,1)
        RKGZ=PMAS(23,2)
        RKVQ=VFQ
        RKAQ=AFQ
        RKVL=VFL
        RKAL=AFL
        RKG1(0)=P(IG1,4)
        RKG2(0)=P(IG2,4)
        RKQ1(0)=P(IQ1,4)
        RKQ2(0)=P(IQ2,4)
        RKL1(0)=P(IL1,4)
        RKL2(0)=P(IL2,4)
        DO 1270 J=1,3
        RKG1(J)=P(IG1,J)
        RKG2(J)=P(IG2,J)
        RKQ1(J)=P(IQ1,J)
        RKQ2(J)=P(IQ2,J)
        RKL1(J)=P(IL1,J)
        RKL2(J)=P(IL2,J)
 1270   CONTINUE
        CALL RKBBV(RKG1,RKG2,RKQ1,RKQ2,RKL1,RKL2,1,RKRES)
 
C...Multiply with normalization factors.
        WTMEP=1./(2.*SH*PARU(2)**8)
        WTCOU=AS**2*(4.*PARU(1)*AEM*XWC)**2
        WTZQQ=WTMEP*WTCOU*RKRES
        WTPHS=(PARU(1)/2.)**2*PMQQ**2*
     &  (PARU(1)*((PMLL**2-PMAS(23,1)**2)**2+(PMAS(23,1)*
     &  PMAS(23,2))**2)/(PMAS(23,1)*PMAS(23,2)))*0.5*SH
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=INT(1.5+RLU(0))
        SIGH(NCHN)=COMFAC*WTPHS*WTZQQ*WID2
 1280   CONTINUE
      ENDIF
 
C...H: 2 -> 1, tree diagrams, non-standard model processes.
 
      ELSEIF(ISUB.LE.160) THEN
      IF(ISUB.EQ.141) THEN
C...f + f~ -> gamma*/Z0/Z'0.
        MINT(61)=2
        CALL PYWIDT(32,SH,WDTP,WDTE)
        HP0=AEM/3.*SH
        HP1=AEM/3.*XWC*SH
        HP2=HP1
        HS=HP1*VINT(117)
        HSP=HP2*WDTP(0)
        FACZP=4.*COMFAC*3.
        DO 1290 I=MMINA,MMAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 1290
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI)
        VI=AI-4.*EI*XWV
        IF(IABS(I).LT.10) THEN
          VPI=PARU(123-2*MOD(IABS(I),2))
          API=PARU(124-2*MOD(IABS(I),2))
        ELSE
          VPI=PARU(127-2*MOD(IABS(I),2))
          API=PARU(128-2*MOD(IABS(I),2))
        ENDIF
        HI0=HP0
        IF(IABS(I).LE.10) HI0=HI0*FACA/3.
        HI1=HP1
        IF(IABS(I).LE.10) HI1=HI1*FACA/3.
        HI2=HP2
        IF(IABS(I).LE.10) HI2=HI2*FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACZP*(EI**2/SH2*HI0*HP0*VINT(111)+EI*VI*
     &  (1.-SQMZ/SH)/((SH-SQMZ)**2+HS**2)*(HI0*HP1+HI1*HP0)*VINT(112)+
     &  EI*VPI*(1.-SQMZP/SH)/((SH-SQMZP)**2+HSP**2)*(HI0*HP2+HI2*HP0)*
     &  VINT(113)+(VI**2+AI**2)/((SH-SQMZ)**2+HS**2)*HI1*HP1*VINT(114)+
     &  (VI*VPI+AI*API)*((SH-SQMZ)*(SH-SQMZP)+HS*HSP)/(((SH-SQMZ)**2+
     &  HS**2)*((SH-SQMZP)**2+HSP**2))*(HI1*HP2+HI2*HP1)*VINT(115)+
     &  (VPI**2+API**2)/((SH-SQMZP)**2+HSP**2)*HI2*HP2*VINT(116))
 1290   CONTINUE
 
      ELSEIF(ISUB.EQ.142) THEN
C...f + f~' -> W'+/-.
        CALL PYWIDT(34,SH,WDTP,WDTE)
        HP=AEM/(24.*XW)*SH
        HS=HP*WDTP(0)
        FACBW=4.*COMFAC/((SH-SQMWP)**2+HS**2)*3.
        DO 1310 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1310
        IA=IABS(I)
        DO 1300 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1300
        JA=IABS(J)
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 1300
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 1300
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        HI=HP*(PARU(133)**2+PARU(134)**2)
        IF(IA.LE.10) HI=HP*(PARU(131)**2+PARU(132)**2)*
     &  VCKM((IA+1)/2,(JA+1)/2)*FACA/3.
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        HF=HP*(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))
        SIGH(NCHN)=HI*FACBW*HF
 1300   CONTINUE
 1310   CONTINUE
 
      ELSEIF(ISUB.EQ.143) THEN
C...f + f~' -> H+/-.
        CALL PYWIDT(37,SH,WDTP,WDTE)
        HP=AEM/(8.*XW)*SH/SQMW*SH
        HS=HP*WDTP(0)
        FACBW=4.*COMFAC/((SH-SQMHC)**2+HS**2)
        DO 1330 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1330
        IA=IABS(I)
        IM=(MOD(IA,10)+1)/2
        DO 1320 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1320
        JA=IABS(J)
        JM=(MOD(JA,10)+1)/2
        IF(I*J.GT.0.OR.IA.EQ.JA.OR.IM.NE.JM) GOTO 1320
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 1320
        IF(MOD(IA,2).EQ.0) THEN
          IU=IA
          IL=JA
        ELSE
          IU=JA
          IL=IA
        ENDIF
        RML=PMAS(IL,1)**2/SH
        RMU=PMAS(IU,1)**2/SH
        IF(IL.LE.10.AND.MSTP(37).EQ.1.AND.MSTP(2).GE.1) RML=RML*
     &  (LOG(MAX(4.,PARP(37)**2*RML*SH/PARU(117)**2))/
     &  LOG(MAX(4.,SH/PARU(117)**2)))**(24./(33.-2.*MSTU(118)))
        HI=HP*(RML*PARU(141)**2+RMU/PARU(141)**2)
        IF(IA.LE.10) HI=HI*FACA/3.
        KCHHC=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
        HF=HP*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
 1320   CONTINUE
 1330   CONTINUE
 
      ELSEIF(ISUB.EQ.144) THEN
C...f + f~' -> R.
        CALL PYWIDT(40,SH,WDTP,WDTE)
        HP=AEM/(12.*XW)*SH
        HS=HP*WDTP(0)
        FACBW=4.*COMFAC/((SH-SQMR)**2+HS**2)*3.
        DO 1350 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1350
        IA=IABS(I)
        DO 1340 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1340
        JA=IABS(J)
        IF(I*J.GT.0.OR.IABS(IA-JA).NE.2) GOTO 1340
        HI=HP
        IF(IA.LE.10) HI=HI*FACA/3.
        HF=HP*(WDTE(0,1)+WDTE(0,(10-(I+J))/4)+WDTE(0,4))
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
 1340   CONTINUE
 1350   CONTINUE
 
      ELSEIF(ISUB.EQ.145) THEN
C...q + l -> LQ (leptoquark).
        CALL PYWIDT(39,SH,WDTP,WDTE)
        HP=AEM/4.*SH
        HS=HP*WDTP(0)
        FACBW=4.*COMFAC/((SH-SQMLQ)**2+HS**2)
        IF(ABS(SH-SQMLQ).GT.100.*HS) FACBW=0.
        KFLQQ=KFDP(MDCY(39,2),1)
        KFLQL=KFDP(MDCY(39,2),2)
        DO 1370 I=MMIN1,MMAX1
        IF(KFAC(1,I).EQ.0) GOTO 1370
        IA=IABS(I)
        IF(IA.NE.KFLQQ.AND.IA.NE.KFLQL) GOTO 1370
        DO 1360 J=MMIN2,MMAX2
        IF(KFAC(2,J).EQ.0) GOTO 1360
        JA=IABS(J)
        IF(JA.NE.KFLQQ.AND.JA.NE.KFLQL) GOTO 1360
        IF(I*J.NE.KFLQQ*KFLQL) GOTO 1360
        IF(IA.EQ.KFLQQ) KCHLQ=ISIGN(1,I)
        IF(JA.EQ.KFLQQ) KCHLQ=ISIGN(1,J)
        HI=HP*PARU(151)
        HF=HP*(WDTE(0,1)+WDTE(0,(5-KCHLQ)/2)+WDTE(0,4))
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
 1360   CONTINUE
 1370   CONTINUE
 
      ELSEIF(ISUB.EQ.147.OR.ISUB.EQ.148) THEN
C...d + g -> d* and u + g -> u* (excited quarks).
        KFQEXC=ISUB-146
        KFQSTR=ISUB-140
        CALL PYWIDT(KFQSTR,SH,WDTP,WDTE)
        HP=SH
        HS=HP*WDTP(0)
        FACBW=COMFAC/((SH-PMAS(KFQSTR,1)**2)**2+HS**2)
        FACBW=FACBW*AS*PARU(159)**2*SH/(3.*PARU(155)**2)
        IF(ABS(SH-PMAS(KFQSTR,1)**2).GT.100.*HS) FACBW=0.
        DO 1390 I=-KFQEXC,KFQEXC,2*KFQEXC
        DO 1380 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 1380
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 1380
        HI=HP
        IF(I.GT.0) HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        IF(I.LT.0) HF=HP*(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
 1380   CONTINUE
 1390   CONTINUE
 
      ELSEIF(ISUB.EQ.149) THEN
C...g + g -> eta_techni.
        CALL PYWIDT(38,SH,WDTP,WDTE)
        HP=SH
        HS=HP*WDTP(0)
        FACBW=COMFAC*0.5/((SH-PMAS(38,1)**2)**2+HS**2)
        IF(ABS(SH-PMAS(38,1)**2).GT.100.*HS) FACBW=0.
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1400
        HI=HP*WDTP(3)
        HF=HP*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=HI*FACBW*HF
 1400   CONTINUE
 
      ENDIF
 
C...I: 2 -> 2, tree diagrams, non-standard model processes.
 
      ELSE
      IF(ISUB.EQ.161) THEN
C...f + g -> f' + H+/- (b + g -> t + H+/- only)
C...(choice of only b and t to avoid kinematics problems).
        FHCQ=COMFAC*FACA*AS*AEM/XW*1./24
        DO 1420 I=MMINA,MMAXA
        IA=IABS(I)
        IF(IA.NE.5) GOTO 1420
        SQML=PMAS(IA,1)**2
        IF(IA.LE.10.AND.MSTP(37).EQ.1.AND.MSTP(2).GE.1) SQML=SQML*
     &  (LOG(MAX(4.,PARP(37)**2*SQML/PARU(117)**2))/
     &  LOG(MAX(4.,SH/PARU(117)**2)))**(24./(33.-2.*MSTU(118)))
        IUA=IA+MOD(IA,2)
        SQMQ=PMAS(IUA,1)**2
        FACHCQ=FHCQ*(SQML*PARU(141)**2+SQMQ/PARU(141)**2)/SQMW*
     &  (SH/(SQMQ-UH)+2.*SQMQ*(SQMHC-UH)/(SQMQ-UH)**2+(SQMQ-UH)/SH+
     &  2.*SQMQ/(SQMQ-UH)+2.*(SQMHC-UH)/(SQMQ-UH)*(SQMHC-SQMQ-SH)/SH)
        KCHHC=ISIGN(1,KCHG(IA,1)*ISIGN(1,I))
        DO 1410 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 1410
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,1).EQ.0) GOTO 1410
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACHCQ*WIDS(37,(5-KCHHC)/2)
 1410   CONTINUE
 1420   CONTINUE
 
      ELSEIF(ISUB.EQ.162) THEN
C...q + g -> LQ + l~; LQ=leptoquark.
        FACLQ=COMFAC*FACA*PARU(151)*(AS*AEM/6.)*(-TH/SH)*
     &  (UH2+SQMLQ**2)/(UH-SQMLQ)**2
        KFLQQ=KFDP(MDCY(39,2),1)
        DO 1440 I=MMINA,MMAXA
        IF(IABS(I).NE.KFLQQ) GOTO 1440
        KCHLQ=ISIGN(1,I)
        DO 1430 ISDE=1,2
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 1430
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 1430
        NCHN=NCHN+1
        ISIG(NCHN,ISDE)=I
        ISIG(NCHN,3-ISDE)=21
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACLQ*WIDS(39,(5-KCHLQ)/2)
 1430   CONTINUE
 1440   CONTINUE
 
      ELSEIF(ISUB.EQ.163) THEN
C...g + g -> LQ + LQ~; LQ=leptoquark.
        FACLQ=COMFAC*FACA*WIDS(39,1)*(AS**2/2.)*
     &  (7./48.+3.*(UH-TH)**2/(16.*SH2))*(1.+2.*SQMLQ*TH/(TH-SQMLQ)**2+
     &  2.*SQMLQ*UH/(UH-SQMLQ)**2+4.*SQMLQ**2/((TH-SQMLQ)*(UH-SQMLQ)))
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 1450
        NCHN=NCHN+1
        ISIG(NCHN,1)=21
        ISIG(NCHN,2)=21
C...Since don't know proper colour flow, randomize between alternatives.
        ISIG(NCHN,3)=INT(1.5+RLU(0))
        SIGH(NCHN)=FACLQ
 1450   CONTINUE
 
      ELSEIF(ISUB.EQ.164) THEN
C...q + q~ -> LQ + LQ~; LQ=leptoquark.
        FACLQA=COMFAC*WIDS(39,1)*(AS**2/9.)*
     &  (SH*(SH-4.*SQMLQ)-(UH-TH)**2)/SH2
        FACLQS=COMFAC*WIDS(39,1)*((PARU(151)**2*AEM**2/8.)*
     &  (-SH*TH-(SQMLQ-TH)**2)/TH2+(PARU(151)*AEM*AS/18.)*
     &  ((SQMLQ-TH)*(UH-TH)+SH*(SQMLQ+TH))/(SH*TH))
        KFLQQ=KFDP(MDCY(39,2),1)
        DO 1460 I=MMINA,MMAXA
        IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &  KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 1460
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=FACLQA
        IF(IABS(I).EQ.KFLQQ) SIGH(NCHN)=FACLQA+FACLQS
 1460   CONTINUE
 
      ELSEIF(ISUB.EQ.165) THEN
C...q + q~ -> l+ + l- (including contact term for compositeness).
        ZRATR=XWC*SH*(SH-SQMZ)/((SH-SQMZ)**2+SQMZ*PMAS(23,2)**2)
        ZRATI=XWC*SH*PMAS(23,1)*PMAS(23,2)/
     &  ((SH-SQMZ)**2+SQMZ*PMAS(23,2)**2)
        KFF=IABS(KFPR(ISUB,1))
        EF=KCHG(KFF,1)/3.
        AF=SIGN(1.,EF+0.1)
        VF=AF-4.*EF*XWV
        VALF=VF+AF
        VARF=VF-AF
        FCOF=1.
        IF(KFF.LE.10) FCOF=3.
        WID2=1.
        IF(KFF.EQ.6.AND.MSTP(48).GE.1) WID2=WIDS(26,1)
        IF((KFF.EQ.7.OR.KFF.EQ.8).AND.MSTP(49).GE.1) WID2=WIDS(KFF+20,1)
        IF((KFF.EQ.17.OR.KFF.EQ.18).AND.MSTP(49).GE.1) WID2=
     &  WIDS(KFF+12,1)
        DO 1470 I=MMINA,MMAXA
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 1470
        EI=KCHG(IABS(I),1)/3.
        AI=SIGN(1.,EI+0.1)
        VI=AI-4.*EI*XWV
        VALI=VI+AI
        VARI=VI-AI
        FCOI=1.
        IF(IABS(I).LE.10) FCOI=FACA/3.
        IF((MSTP(5).EQ.1.AND.IABS(I).LE.2).OR.MSTP(5).EQ.2) THEN
          FGZA=(EI*EF+VALI*VALF*ZRATR+PARU(156)*SH/
     &    (AEM*PARU(155)**2))**2+(VALI*VALF*ZRATI)**2+
     &    (EI*EF+VARI*VARF*ZRATR)**2+(VARI*VARF*ZRATI)**2
        ELSE
          FGZA=(EI*EF+VALI*VALF*ZRATR)**2+(VALI*VALF*ZRATI)**2+
     &    (EI*EF+VARI*VARF*ZRATR)**2+(VARI*VARF*ZRATI)**2
        ENDIF
        FGZB=(EI*EF+VALI*VARF*ZRATR)**2+(VALI*VARF*ZRATI)**2+
     &  (EI*EF+VARI*VALF*ZRATR)**2+(VARI*VALF*ZRATI)**2
        FGZAB=AEM**2*(FGZA*UH2/SH2+FGZB*TH2/SH2)
        IF((MSTP(5).EQ.3.AND.IABS(I).EQ.2).OR.(MSTP(5).EQ.4.AND.
     &  MOD(IABS(I),2).EQ.0)) FGZAB=FGZAB+SH2/(2.*PARU(155)**4)
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=-I
        ISIG(NCHN,3)=1
        SIGH(NCHN)=COMFAC*FCOI*FCOF*FGZAB*WID2
 1470   CONTINUE
 
      ELSEIF(ISUB.EQ.166) THEN
C...q + q'~ -> l + nu_l (including contact term for compositeness).
        WFAC=(1./4.)*(AEM/XW)**2*UH2/((SH-SQMW)**2+SQMW*PMAS(24,2)**2)
        WCIFAC=WFAC+SH2/(4.*PARU(155)**4)
        KFF=IABS(KFPR(ISUB,1))
        FCOF=1.
        IF(KFF.LE.10) FCOF=3.
        DO 1490 I=MMIN1,MMAX1
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 1490
        IA=IABS(I)
        DO 1480 J=MMIN2,MMAX2
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 1480
        JA=IABS(J)
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 1480
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 1480
        FCOI=1.
        IF(IA.LE.10) FCOI=VCKM((IA+1)/2,(JA+1)/2)*FACA/3.
        WID2=1.
        IF((I.GT.0.AND.MOD(I,2).EQ.0).OR.(J.GT.0.AND.MOD(J,2).EQ.0))
     &  THEN
          IF(KFF.EQ.5.AND.MSTP(48).GE.1) WID2=WIDS(26,2)
          IF(KFF.EQ.7.AND.MSTP(49).GE.1) WID2=WIDS(28,2)*WIDS(27,3)
          IF(KFF.EQ.17.AND.MSTP(49).GE.1) WID2=WIDS(30,2)*WIDS(29,3)
        ELSE
          IF(KFF.EQ.5.AND.MSTP(48).GE.1) WID2=WIDS(26,3)
          IF(KFF.EQ.7.AND.MSTP(49).GE.1) WID2=WIDS(28,3)*WIDS(27,2)
          IF(KFF.EQ.17.AND.MSTP(49).GE.1) WID2=WIDS(30,3)*WIDS(29,2)
        ENDIF
        NCHN=NCHN+1
        ISIG(NCHN,1)=I
        ISIG(NCHN,2)=J
        ISIG(NCHN,3)=1
        SIGH(NCHN)=COMFAC*FCOI*FCOF*WFAC*WID2
        IF((MSTP(5).EQ.3.AND.IA.LE.2.AND.JA.LE.2).OR.MSTP(5).EQ.4)
     &  SIGH(NCHN)=COMFAC*FCOI*FCOF*WCIFAC*WID2
 1480   CONTINUE
 1490   CONTINUE
 
      ELSEIF(ISUB.EQ.167.OR.ISUB.EQ.168) THEN
C...d + g -> d* and u + g -> u* (excited quarks).
        KFQEXC=ISUB-166
        KFQSTR=ISUB-160
        FACQSA=COMFAC*(SH/PARU(155)**2)**2*(1.-SQM4/SH)
        FACQSB=COMFAC*0.25*(SH/PARU(155)**2)**2*(1.-SQM4/SH)*
     &  (1.+SQM4/SH)*(1.+CTH)*(1.+((SH-SQM4)/(SH+SQM4))*CTH)
C...Propagators: as simulated in PYOFSH and as desired.
        GMMQ=PMAS(KFQSTR,1)*PMAS(KFQSTR,2)
        HBW4=GMMQ/((SQM4-PMAS(KFQSTR,1)**2)**2+GMMQ**2)
        CALL PYWIDT(KFQSTR,SQM4,WDTP,WDTE)
        GMMQC=SQM4*WDTP(0)
        HBW4C=GMMQC/((SQM4-PMAS(KFQSTR,1)**2)**2+GMMQC**2)
        FACQSA=FACQSA*HBW4C/HBW4
        FACQSB=FACQSB*HBW4C/HBW4
        DO 1510 I=MMIN1,MMAX1
        IA=IABS(I)
        IF(I.EQ.0.OR.IA.GT.6.OR.KFAC(1,I).EQ.0) GOTO 1510
        DO 1500 J=MMIN2,MMAX2
        JA=IABS(J)
        IF(J.EQ.0.OR.JA.GT.6.OR.KFAC(2,J).EQ.0) GOTO 1500
        IF(IA.EQ.KFQEXC.AND.I.EQ.J) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=1
          SIGH(NCHN)=(4./3.)*FACQSA
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=2
          SIGH(NCHN)=(4./3.)*FACQSA
        ELSEIF((IA.EQ.KFQEXC.OR.JA.EQ.KFQEXC).AND.I*J.GT.0) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=1
          IF(JA.EQ.KFQEXC) ISIG(NCHN,3)=2
          SIGH(NCHN)=FACQSA
        ELSEIF(IA.EQ.KFQEXC.AND.I.EQ.-J) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=1
          SIGH(NCHN)=(8./3.)*FACQSB
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=2
          SIGH(NCHN)=(8./3.)*FACQSB
        ELSEIF(I.EQ.-J) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQSB
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=2
          SIGH(NCHN)=FACQSB
        ELSEIF(IA.EQ.KFQEXC.OR.JA.EQ.KFQEXC) THEN
          NCHN=NCHN+1
          ISIG(NCHN,1)=I
          ISIG(NCHN,2)=J
          ISIG(NCHN,3)=1
          IF(JA.EQ.KFQEXC) ISIG(NCHN,3)=2
          SIGH(NCHN)=FACQSB
        ENDIF
 1500   CONTINUE
 1510   CONTINUE
 
      ENDIF
      ENDIF
 
C...Multiply with structure functions.
      IF(ISUB.LE.90.OR.ISUB.GE.96) THEN
        DO 1520 ICHN=1,NCHN
        IF(MINT(45).GE.2) THEN
          KFL1=ISIG(ICHN,1)
          SIGH(ICHN)=SIGH(ICHN)*XSFX(1,KFL1)
        ENDIF
        IF(MINT(46).GE.2) THEN
          KFL2=ISIG(ICHN,2)
          SIGH(ICHN)=SIGH(ICHN)*XSFX(2,KFL2)
        ENDIF
        SIGS=SIGS+SIGH(ICHN)
 1520   CONTINUE
      ENDIF
 
      RETURN
      END
