 
C*********************************************************************
 
      SUBROUTINE PYMAXI
 
C...Finds optimal set of coefficients for kinematical variable selection
C...and the maximum of the part of the differential cross-section used
C...in the event weighting.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3)
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3)
      COMMON/PYINT6/PROC(0:200)
      CHARACTER PROC*28
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT3/,/PYINT4/,
     &/PYINT5/,/PYINT6/,/PYINT7/
      CHARACTER CVAR(4)*4
      DIMENSION NPTS(4),MVARPT(500,4),VINTPT(500,30),SIGSPT(500),
     &NAREL(7),WTREL(7),WTMAT(7,7),WTRELN(7),COEFU(7),COEFO(7),
     &IACCMX(4),SIGSMX(4),SIGSSM(3)
      DATA CVAR/'tau ','tau''','y*  ','cth '/
      DATA SIGSSM/3*0./
 
C...Select subprocess to study: skip cases not applicable.
      NPOSI=0
      VINT(143)=1.
      VINT(144)=1.
      XSEC(0,1)=0.
      DO 440 ISUB=1,200
      MINT(51)=0
      IF(ISET(ISUB).EQ.11) THEN
        XSEC(ISUB,1)=1.00001*COEF(ISUB,1)
        NPOSI=NPOSI+1
        GOTO 430
      ELSEIF(ISUB.GE.91.AND.ISUB.LE.95) THEN
        XSEC(ISUB,1)=SIGT(0,0,ISUB-90)
        IF(MSUB(ISUB).NE.1) GOTO 440
        NPOSI=NPOSI+1
        GOTO 430
      ELSEIF(ISUB.EQ.96) THEN
        IF(MINT(50).EQ.0) GOTO 440
        IF(MSUB(95).NE.1.AND.MSTP(81).LE.0.AND.MSTP(131).LE.0) GOTO 440
        IF(MINT(49).EQ.0.AND.MSTP(131).EQ.0) GOTO 440
      ELSEIF(ISUB.EQ.11.OR.ISUB.EQ.12.OR.ISUB.EQ.13.OR.ISUB.EQ.28.OR.
     &ISUB.EQ.53.OR.ISUB.EQ.68) THEN
        IF(MSUB(ISUB).NE.1.OR.MSUB(95).EQ.1) GOTO 440
      ELSE
        IF(MSUB(ISUB).NE.1) GOTO 440
      ENDIF
      MINT(1)=ISUB
      ISTSB=ISET(ISUB)
      IF(ISUB.EQ.96) ISTSB=2
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5000) ISUB
      MWTXS=0
      IF(MSTP(142).GE.1.AND.ISUB.NE.96.AND.MSUB(91)+MSUB(92)+MSUB(93)+
     &MSUB(94)+MSUB(95).EQ.0) MWTXS=1
 
C...Find resonances (explicit or implicit in cross-section).
      MINT(72)=0
      KFR1=0
      IF(ISTSB.EQ.1.OR.ISTSB.EQ.3.OR.ISTSB.EQ.5) THEN
        KFR1=KFPR(ISUB,1)
      ELSEIF(ISUB.EQ.24.OR.ISUB.EQ.25.OR.ISUB.EQ.110.OR.ISUB.EQ.165.OR.
     &ISUB.EQ.171.OR.ISUB.EQ.176) THEN
        KFR1=23
      ELSEIF(ISUB.EQ.23.OR.ISUB.EQ.26.OR.ISUB.EQ.166.OR.ISUB.EQ.172.OR.
     &ISUB.EQ.177) THEN
        KFR1=24
      ELSEIF(ISUB.GE.71.AND.ISUB.LE.77) THEN
        KFR1=25
        IF(MSTP(46).EQ.5) THEN
          KFR1=30
          PMAS(30,1)=PARP(45)
          PMAS(30,2)=PARP(45)**3/(96.*PARU(1)*PARP(47)**2)
        ENDIF
      ENDIF
      CKMX=CKIN(2)
      IF(CKMX.LE.0.) CKMX=VINT(1)
      IF(KFR1.NE.0) THEN
        IF(CKIN(1).GT.PMAS(KFR1,1)+20.*PMAS(KFR1,2).OR.
     &  CKMX.LT.PMAS(KFR1,1)-20.*PMAS(KFR1,2)) KFR1=0
      ENDIF
      IF(KFR1.NE.0) THEN
        TAUR1=PMAS(KFR1,1)**2/VINT(2)
        GAMR1=PMAS(KFR1,1)*PMAS(KFR1,2)/VINT(2)
        MINT(72)=1
        MINT(73)=KFR1
        VINT(73)=TAUR1
        VINT(74)=GAMR1
      ENDIF
      IF(ISUB.EQ.141) THEN
        KFR2=23
        TAUR2=PMAS(KFR2,1)**2/VINT(2)
        GAMR2=PMAS(KFR2,1)*PMAS(KFR2,2)/VINT(2)
        IF(CKIN(1).GT.PMAS(KFR2,1)+20.*PMAS(KFR2,2).OR.
     &  CKMX.LT.PMAS(KFR2,1)-20.*PMAS(KFR2,2)) KFR2=0
        IF(KFR2.NE.0.AND.KFR1.NE.0) THEN
          MINT(72)=2
          MINT(74)=KFR2
          VINT(75)=TAUR2
          VINT(76)=GAMR2
        ELSEIF(KFR2.NE.0) THEN
          KFR1=KFR2
          TAUR1=TAUR2
          GAMR1=GAMR2
          MINT(72)=1
          MINT(73)=KFR1
          VINT(73)=TAUR1
          VINT(74)=GAMR1
        ENDIF
      ENDIF
 
C...Find product masses and minimum pT of process.
      SQM3=0.
      SQM4=0.
      MINT(71)=0
      VINT(71)=CKIN(3)
      VINT(80)=1.
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN
        NBW=0
        DO 100 I=1,2
        IF(KFPR(ISUB,I).EQ.0) THEN
        ELSEIF(MSTP(42).LE.0.OR.PMAS(LUCOMP(KFPR(ISUB,I)),2).LT.
     &  PARP(41)) THEN
          IF(I.EQ.1) SQM3=PMAS(LUCOMP(KFPR(ISUB,I)),1)**2
          IF(I.EQ.2) SQM4=PMAS(LUCOMP(KFPR(ISUB,I)),1)**2
        ELSE
          NBW=NBW+1
        ENDIF
  100   CONTINUE
        IF(NBW.GE.1) THEN
          CALL PYOFSH(3,0,KFPR(ISUB,1),KFPR(ISUB,2),0.,PQM3,PQM4)
          IF(MINT(51).EQ.1) THEN
            WRITE(MSTU(11),5100) ISUB
            MSUB(ISUB)=0
            GOTO 440
          ENDIF
          SQM3=PQM3**2
          SQM4=PQM4**2
        ENDIF
        IF(MIN(SQM3,SQM4).LT.CKIN(6)**2) MINT(71)=1
        IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5))
        IF(ISUB.EQ.96.AND.MSTP(82).LE.1) VINT(71)=PARP(81)
        IF(ISUB.EQ.96.AND.MSTP(82).GE.2) VINT(71)=0.08*PARP(82)
      ELSEIF(ISTSB.EQ.6) THEN
        CALL PYOFSH(5,0,KFPR(ISUB,1),KFPR(ISUB,2),0.,PQM3,PQM4)
        IF(MINT(51).EQ.1) THEN
          WRITE(MSTU(11),5100) ISUB
          MSUB(ISUB)=0
          GOTO 440
        ENDIF
        SQM3=PQM3**2
        SQM4=PQM4**2
      ENDIF
      VINT(63)=SQM3
      VINT(64)=SQM4
 
C...Prepare for additional variable choices in 2 -> 3.
      IF(ISTSB.EQ.5) THEN
        VINT(201)=0.
        IF(KFPR(ISUB,2).GT.0) VINT(201)=PMAS(KFPR(ISUB,2),1)
        VINT(206)=VINT(201)
        VINT(204)=PMAS(23,1)
        IF(ISUB.EQ.124) VINT(204)=PMAS(24,1)
        IF(ISUB.EQ.121.OR.ISUB.EQ.122.OR.ISUB.EQ.181.OR.ISUB.EQ.182.OR.
     &  ISUB.EQ.186.OR.ISUB.EQ.187) VINT(204)=VINT(201)
        VINT(209)=VINT(204)
      ENDIF
 
C...Number of points for each variable: tau, tau', y*, cos(theta-hat).
      NPTS(1)=2+2*MINT(72)
      IF(MINT(47).EQ.1) THEN
        IF(ISTSB.EQ.1.OR.ISTSB.EQ.2.OR.ISTSB.EQ.6) NPTS(1)=1
      ELSEIF(MINT(47).EQ.5) THEN
        IF(ISTSB.LE.2.OR.ISTSB.GE.6) NPTS(1)=NPTS(1)+1
      ENDIF
      NPTS(2)=1
      IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
        IF(MINT(47).GE.2) NPTS(2)=2
        IF(MINT(47).EQ.5) NPTS(2)=3
      ENDIF
      NPTS(3)=1
      IF(MINT(47).GE.4) NPTS(3)=3
      IF(MINT(45).EQ.3) NPTS(3)=NPTS(3)+1
      IF(MINT(46).EQ.3) NPTS(3)=NPTS(3)+1
      NPTS(4)=1
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) NPTS(4)=5
      NTRY=NPTS(1)*NPTS(2)*NPTS(3)*NPTS(4)
 
C...Reset coefficients of cross-section weighting.
      DO 110 J=1,20
      COEF(ISUB,J)=0.
  110 CONTINUE
      COEF(ISUB,1)=1.
      COEF(ISUB,8)=0.5
      COEF(ISUB,9)=0.5
      COEF(ISUB,13)=1.
      COEF(ISUB,18)=1.
      MCTH=0
      MTAUP=0
      METAUP=0
      VINT(23)=0.
      VINT(26)=0.
      SIGSAM=0.
 
C...Find limits and select tau, y*, cos(theta-hat) and tau' values,
C...in grid of phase space points.
      CALL PYKLIM(1)
      METAU=MINT(51)
      NACC=0
      DO 140 ITRY=1,NTRY
      MINT(51)=0
      IF(METAU.EQ.1) GOTO 140
      IF(MOD(ITRY-1,NPTS(2)*NPTS(3)*NPTS(4)).EQ.0) THEN
        MTAU=1+(ITRY-1)/(NPTS(2)*NPTS(3)*NPTS(4))
        IF(MTAU.GT.2+2*MINT(72)) MTAU=7
        CALL PYKMAP(1,MTAU,0.5)
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) CALL PYKLIM(4)
        METAUP=MINT(51)
      ENDIF
      IF(METAUP.EQ.1) GOTO 140
      IF(ISTSB.GE.3.AND.ISTSB.LE.5.AND.MOD(ITRY-1,NPTS(3)*NPTS(4))
     &.EQ.0) THEN
        MTAUP=1+MOD((ITRY-1)/(NPTS(3)*NPTS(4)),NPTS(2))
        CALL PYKMAP(4,MTAUP,0.5)
      ENDIF
      IF(MOD(ITRY-1,NPTS(3)*NPTS(4)).EQ.0) THEN
        CALL PYKLIM(2)
        MEYST=MINT(51)
      ENDIF
      IF(MEYST.EQ.1) GOTO 140
      IF(MOD(ITRY-1,NPTS(4)).EQ.0) THEN
        MYST=1+MOD((ITRY-1)/NPTS(4),NPTS(3))
        IF(MYST.EQ.4.AND.MINT(45).NE.3) MYST=5
        CALL PYKMAP(2,MYST,0.5)
        CALL PYKLIM(3)
        MECTH=MINT(51)
      ENDIF
      IF(MECTH.EQ.1) GOTO 140
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) THEN
        MCTH=1+MOD(ITRY-1,NPTS(4))
        CALL PYKMAP(3,MCTH,0.5)
      ENDIF
      IF(ISUB.EQ.96) VINT(25)=VINT(21)*(1.-VINT(23)**2)
 
C...Store position and limits.
      MINT(51)=0
      CALL PYKLIM(0)
      IF(MINT(51).EQ.1) GOTO 140
      NACC=NACC+1
      MVARPT(NACC,1)=MTAU
      MVARPT(NACC,2)=MTAUP
      MVARPT(NACC,3)=MYST
      MVARPT(NACC,4)=MCTH
      DO 120 J=1,30
      VINTPT(NACC,J)=VINT(10+J)
  120 CONTINUE
 
C...Normal case: calculate cross-section.
      IF(ISTSB.NE.5) THEN
        CALL PYSIGH(NCHN,SIGS)
        IF(MWTXS.EQ.1) THEN
          CALL PYEVWT(WTXS)
          SIGS=WTXS*SIGS
        ENDIF
 
C..2 -> 3: find highest value out of a number of tries.
      ELSE
        SIGS=0.
        DO 130 IKIN3=1,MSTP(129)
        CALL PYKMAP(5,0,0.)
        IF(MINT(51).EQ.1) GOTO 130
        CALL PYSIGH(NCHN,SIGTMP)
        IF(MWTXS.EQ.1) THEN
          CALL PYEVWT(WTXS)
          SIGTMP=WTXS*SIGTMP
        ENDIF
        IF(SIGTMP.GT.SIGS) SIGS=SIGTMP
  130   CONTINUE
      ENDIF
 
C...Store cross-section.
      SIGSPT(NACC)=SIGS
      IF(SIGS.GT.SIGSAM) SIGSAM=SIGS
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5200) MTAU,MYST,MCTH,MTAUP,
     &VINT(21),VINT(22),VINT(23),VINT(26),SIGS
  140 CONTINUE
      IF(NACC.EQ.0) THEN
        WRITE(MSTU(11),5100) ISUB
        MSUB(ISUB)=0
        GOTO 440
      ELSEIF(SIGSAM.EQ.0.) THEN
        WRITE(MSTU(11),5300) ISUB
        MSUB(ISUB)=0
        GOTO 440
      ENDIF
      IF(ISUB.NE.96) NPOSI=NPOSI+1
 
C...Calculate integrals in tau over maximal phase space limits.
      TAUMIN=VINT(11)
      TAUMAX=VINT(31)
      ATAU1=LOG(TAUMAX/TAUMIN)
      IF(NPTS(1).GE.2) THEN
        ATAU2=(TAUMAX-TAUMIN)/(TAUMAX*TAUMIN)
      ENDIF
      IF(NPTS(1).GE.4) THEN
        ATAU3=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR1)/(TAUMAX+TAUR1))/TAUR1
        ATAU4=(ATAN((TAUMAX-TAUR1)/GAMR1)-ATAN((TAUMIN-TAUR1)/GAMR1))/
     &  GAMR1
      ENDIF
      IF(NPTS(1).GE.6) THEN
        ATAU5=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR2)/(TAUMAX+TAUR2))/TAUR2
        ATAU6=(ATAN((TAUMAX-TAUR2)/GAMR2)-ATAN((TAUMIN-TAUR2)/GAMR2))/
     &  GAMR2
      ENDIF
      IF(NPTS(1).GT.2+2*MINT(72)) THEN
        ATAU7=LOG(MAX(2E-6,1.-TAUMIN)/MAX(2E-6,1.-TAUMAX))
      ENDIF
 
C...Reset. Sum up cross-sections in points calculated.
      DO 300 IVAR=1,4
      IF(NPTS(IVAR).EQ.1) GOTO 300
      IF(ISUB.EQ.96.AND.IVAR.EQ.4) GOTO 300
      NBIN=NPTS(IVAR)
      DO 160 J1=1,NBIN
      NAREL(J1)=0
      WTREL(J1)=0.
      COEFU(J1)=0.
      DO 150 J2=1,NBIN
      WTMAT(J1,J2)=0.
  150 CONTINUE
  160 CONTINUE
      DO 170 IACC=1,NACC
      IBIN=MVARPT(IACC,IVAR)
      IF(IVAR.EQ.1.AND.IBIN.EQ.7) IBIN=3+2*MINT(72)
      IF(IVAR.EQ.3.AND.IBIN.EQ.5.AND.MINT(45).NE.3) IBIN=4
      NAREL(IBIN)=NAREL(IBIN)+1
      WTREL(IBIN)=WTREL(IBIN)+SIGSPT(IACC)
 
C...Sum up tau cross-section pieces in points used.
      IF(IVAR.EQ.1) THEN
        TAU=VINTPT(IACC,11)
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ATAU1/ATAU2)/TAU
        IF(NBIN.GE.4) THEN
          WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ATAU1/ATAU3)/(TAU+TAUR1)
          WTMAT(IBIN,4)=WTMAT(IBIN,4)+(ATAU1/ATAU4)*TAU/
     &    ((TAU-TAUR1)**2+GAMR1**2)
        ENDIF
        IF(NBIN.GE.6) THEN
          WTMAT(IBIN,5)=WTMAT(IBIN,5)+(ATAU1/ATAU5)/(TAU+TAUR2)
          WTMAT(IBIN,6)=WTMAT(IBIN,6)+(ATAU1/ATAU6)*TAU/
     &    ((TAU-TAUR2)**2+GAMR2**2)
        ENDIF
        IF(NBIN.GT.2+2*MINT(72)) THEN
          WTMAT(IBIN,NBIN)=WTMAT(IBIN,NBIN)+(ATAU1/ATAU7)*
     &    TAU/MAX(2E-6,1.-TAU)
        ENDIF
 
C...Sum up tau' cross-section pieces in points used.
      ELSEIF(IVAR.EQ.2) THEN
        TAU=VINTPT(IACC,11)
        TAUP=VINTPT(IACC,16)
        TAUPMN=VINTPT(IACC,6)
        TAUPMX=VINTPT(IACC,26)
        ATAUP1=LOG(TAUPMX/TAUPMN)
        ATAUP2=((1.-TAU/TAUPMX)**4-(1.-TAU/TAUPMN)**4)/(4.*TAU)
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ATAUP1/ATAUP2)*(1.-TAU/TAUP)**3/
     &  TAUP
        IF(NBIN.GE.3) THEN
          ATAUP3=LOG(MAX(2E-6,1.-TAUPMN)/MAX(2E-6,1.-TAUPMX))
          WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ATAUP1/ATAUP3)*
     &    TAUP/MAX(2E-6,1.-TAUP)
        ENDIF
 
C...Sum up y* cross-section pieces in points used.
      ELSEIF(IVAR.EQ.3) THEN
        YST=VINTPT(IACC,12)
        YSTMIN=VINTPT(IACC,2)
        YSTMAX=VINTPT(IACC,22)
        AYST0=YSTMAX-YSTMIN
        AYST1=0.5*(YSTMAX-YSTMIN)**2
        AYST2=AYST1
        AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+(AYST0/AYST1)*(YST-YSTMIN)
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(AYST0/AYST2)*(YSTMAX-YST)
        WTMAT(IBIN,3)=WTMAT(IBIN,3)+(AYST0/AYST3)/COSH(YST)
        IF(MINT(45).EQ.3) THEN
          TAUE=VINTPT(IACC,11)
          IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=VINTPT(IACC,16)
          YST0=-0.5*LOG(TAUE)
          AYST4=LOG(MAX(1E-6,EXP(YST0-YSTMIN)-1.)/
     &    MAX(1E-6,EXP(YST0-YSTMAX)-1.))
          WTMAT(IBIN,4)=WTMAT(IBIN,4)+(AYST0/AYST4)/
     &    MAX(1E-6,1.-EXP(YST-YST0))
        ENDIF
        IF(MINT(46).EQ.3) THEN
          TAUE=VINTPT(IACC,11)
          IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=VINTPT(IACC,16)
          YST0=-0.5*LOG(TAUE)
          AYST5=LOG(MAX(1E-6,EXP(YST0+YSTMAX)-1.)/
     &    MAX(1E-6,EXP(YST0+YSTMIN)-1.))
          WTMAT(IBIN,NBIN)=WTMAT(IBIN,NBIN)+(AYST0/AYST5)/
     &    MAX(1E-6,1.-EXP(-YST-YST0))
        ENDIF
 
C...Sum up cos(theta-hat) cross-section pieces in points used.
      ELSE
        RM34=MAX(1E-20,2.*SQM3*SQM4/(VINTPT(IACC,11)*VINT(2))**2)
        RSQM=1.+RM34
        CTHMAX=SQRT(1.-4.*VINT(71)**2/(TAUMAX*VINT(2)))
        CTHMIN=-CTHMAX
        IF(CTHMAX.GT.0.9999) RM34=MAX(RM34,2.*VINT(71)**2/
     &  (TAUMAX*VINT(2)))
        ACTH1=CTHMAX-CTHMIN
        ACTH2=LOG(MAX(RM34,RSQM-CTHMIN)/MAX(RM34,RSQM-CTHMAX))
        ACTH3=LOG(MAX(RM34,RSQM+CTHMAX)/MAX(RM34,RSQM+CTHMIN))
        ACTH4=1./MAX(RM34,RSQM-CTHMAX)-1./MAX(RM34,RSQM-CTHMIN)
        ACTH5=1./MAX(RM34,RSQM+CTHMIN)-1./MAX(RM34,RSQM+CTHMAX)
        CTH=VINTPT(IACC,13)
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ACTH1/ACTH2)/MAX(RM34,RSQM-CTH)
        WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ACTH1/ACTH3)/MAX(RM34,RSQM+CTH)
        WTMAT(IBIN,4)=WTMAT(IBIN,4)+(ACTH1/ACTH4)/MAX(RM34,RSQM-CTH)**2
        WTMAT(IBIN,5)=WTMAT(IBIN,5)+(ACTH1/ACTH5)/MAX(RM34,RSQM+CTH)**2
      ENDIF
  170 CONTINUE
 
C...Check that equation system solvable; else trivial way out.
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5400) CVAR(IVAR)
      MSOLV=1
      WTRELS=0.
      DO 180 IBIN=1,NBIN
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5500) (WTMAT(IBIN,IRED),
     &IRED=1,NBIN),WTREL(IBIN)
      IF(NAREL(IBIN).EQ.0) MSOLV=0
      WTRELS=WTRELS+WTREL(IBIN)
  180 CONTINUE
      IF(MSOLV.EQ.0) THEN
        DO 190 IBIN=1,NBIN
        COEFU(IBIN)=1.
        WTRELN(IBIN)=0.1
        IF(WTRELS.GT.0.) WTRELN(IBIN)=MAX(0.1,WTREL(IBIN)/WTRELS)
  190   CONTINUE
 
C...Solve to find relative importance of cross-section pieces.
      ELSE
        DO 200 IBIN=1,NBIN
        WTRELN(IBIN)=MAX(0.1,WTREL(IBIN)/WTRELS)
  200   CONTINUE
        DO 230 IRED=1,NBIN-1
        DO 220 IBIN=IRED+1,NBIN
        RQT=WTMAT(IBIN,IRED)/WTMAT(IRED,IRED)
        WTREL(IBIN)=WTREL(IBIN)-RQT*WTREL(IRED)
        DO 210 ICOE=IRED,NBIN
        WTMAT(IBIN,ICOE)=WTMAT(IBIN,ICOE)-RQT*WTMAT(IRED,ICOE)
  210   CONTINUE
  220   CONTINUE
  230   CONTINUE
        DO 250 IRED=NBIN,1,-1
        DO 240 ICOE=IRED+1,NBIN
        WTREL(IRED)=WTREL(IRED)-WTMAT(IRED,ICOE)*COEFU(ICOE)
  240   CONTINUE
        COEFU(IRED)=WTREL(IRED)/WTMAT(IRED,IRED)
  250   CONTINUE
      ENDIF
 
C...Normalize coefficients, with piece shared democratically.
      COEFSU=0.
      WTRELS=0.
      DO 260 IBIN=1,NBIN
      COEFU(IBIN)=MAX(0.,COEFU(IBIN))
      COEFSU=COEFSU+COEFU(IBIN)
      WTRELS=WTRELS+WTRELN(IBIN)
  260 CONTINUE
      IF(COEFSU.GT.0.) THEN
        DO 270 IBIN=1,NBIN
        COEFO(IBIN)=PARP(122)/NBIN+(1.-PARP(122))*0.5*
     &  (COEFU(IBIN)/COEFSU+WTRELN(IBIN)/WTRELS)
  270   CONTINUE
      ELSE
        DO 280 IBIN=1,NBIN
        COEFO(IBIN)=1./NBIN
  280   CONTINUE
      ENDIF
      IF(IVAR.EQ.1) IOFF=0
      IF(IVAR.EQ.2) IOFF=17
      IF(IVAR.EQ.3) IOFF=7
      IF(IVAR.EQ.4) IOFF=12
      DO 290 IBIN=1,NBIN
      ICOF=IOFF+IBIN
      IF(IVAR.EQ.1.AND.IBIN.GT.2+2*MINT(72)) ICOF=7
      IF(IVAR.EQ.3.AND.IBIN.EQ.4.AND.MINT(45).NE.3) ICOF=ICOF+1
      COEF(ISUB,ICOF)=COEFO(IBIN)
  290 CONTINUE
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5600) CVAR(IVAR),
     &(COEFO(IBIN),IBIN=1,NBIN)
  300 CONTINUE
 
C...Find two most promising maxima among points previously determined.
      DO 310 J=1,4
      IACCMX(J)=0
      SIGSMX(J)=0.
  310 CONTINUE
      NMAX=0
      DO 370 IACC=1,NACC
      DO 320 J=1,30
      VINT(10+J)=VINTPT(IACC,J)
  320 CONTINUE
      IF(ISTSB.NE.5) THEN
        CALL PYSIGH(NCHN,SIGS)
        IF(MWTXS.EQ.1) THEN
          CALL PYEVWT(WTXS)
          SIGS=WTXS*SIGS
        ENDIF
      ELSE
        SIGS=0.
        DO 330 IKIN3=1,MSTP(129)
        CALL PYKMAP(5,0,0.)
        IF(MINT(51).EQ.1) GOTO 330
        CALL PYSIGH(NCHN,SIGTMP)
        IF(MWTXS.EQ.1) THEN
          CALL PYEVWT(WTXS)
          SIGTMP=WTXS*SIGTMP
        ENDIF
        IF(SIGTMP.GT.SIGS) SIGS=SIGTMP
  330   CONTINUE
      ENDIF
      IEQ=0
      DO 340 IMV=1,NMAX
      IF(ABS(SIGS-SIGSMX(IMV)).LT.1E-4*(SIGS+SIGSMX(IMV))) IEQ=IMV
  340 CONTINUE
      IF(IEQ.EQ.0) THEN
        DO 350 IMV=NMAX,1,-1
        IIN=IMV+1
        IF(SIGS.LE.SIGSMX(IMV)) GOTO 360
        IACCMX(IMV+1)=IACCMX(IMV)
        SIGSMX(IMV+1)=SIGSMX(IMV)
  350   CONTINUE
        IIN=1
  360   IACCMX(IIN)=IACC
        SIGSMX(IIN)=SIGS
        IF(NMAX.LE.1) NMAX=NMAX+1
      ENDIF
  370 CONTINUE
 
C...Read out starting position for search.
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5700)
      SIGSAM=SIGSMX(1)
      DO 420 IMAX=1,NMAX
      IACC=IACCMX(IMAX)
      MTAU=MVARPT(IACC,1)
      MTAUP=MVARPT(IACC,2)
      MYST=MVARPT(IACC,3)
      MCTH=MVARPT(IACC,4)
      VTAU=0.5
      VYST=0.5
      VCTH=0.5
      VTAUP=0.5
 
C...Starting point and step size in parameter space.
      DO 410 IRPT=1,2
      DO 400 IVAR=1,4
      IF(NPTS(IVAR).EQ.1) GOTO 400
      IF(IVAR.EQ.1) VVAR=VTAU
      IF(IVAR.EQ.2) VVAR=VTAUP
      IF(IVAR.EQ.3) VVAR=VYST
      IF(IVAR.EQ.4) VVAR=VCTH
      IF(IVAR.EQ.1) MVAR=MTAU
      IF(IVAR.EQ.2) MVAR=MTAUP
      IF(IVAR.EQ.3) MVAR=MYST
      IF(IVAR.EQ.4) MVAR=MCTH
      IF(IRPT.EQ.1) VDEL=0.1
      IF(IRPT.EQ.2) VDEL=MAX(0.01,MIN(0.05,VVAR-0.02,0.98-VVAR))
      IF(IRPT.EQ.1) VMAR=0.02
      IF(IRPT.EQ.2) VMAR=0.002
      IMOV0=1
      IF(IRPT.EQ.1.AND.IVAR.EQ.1) IMOV0=0
      DO 390 IMOV=IMOV0,8
 
C...Define new point in parameter space.
      IF(IMOV.EQ.0) THEN
        INEW=2
        VNEW=VVAR
      ELSEIF(IMOV.EQ.1) THEN
        INEW=3
        VNEW=VVAR+VDEL
      ELSEIF(IMOV.EQ.2) THEN
        INEW=1
        VNEW=VVAR-VDEL
      ELSEIF(SIGSSM(3).GE.MAX(SIGSSM(1),SIGSSM(2)).AND.
     &VVAR+2.*VDEL.LT.1.-VMAR) THEN
        VVAR=VVAR+VDEL
        SIGSSM(1)=SIGSSM(2)
        SIGSSM(2)=SIGSSM(3)
        INEW=3
        VNEW=VVAR+VDEL
      ELSEIF(SIGSSM(1).GE.MAX(SIGSSM(2),SIGSSM(3)).AND.
     &VVAR-2.*VDEL.GT.VMAR) THEN
        VVAR=VVAR-VDEL
        SIGSSM(3)=SIGSSM(2)
        SIGSSM(2)=SIGSSM(1)
        INEW=1
        VNEW=VVAR-VDEL
      ELSEIF(SIGSSM(3).GE.SIGSSM(1)) THEN
        VDEL=0.5*VDEL
        VVAR=VVAR+VDEL
        SIGSSM(1)=SIGSSM(2)
        INEW=2
        VNEW=VVAR
      ELSE
        VDEL=0.5*VDEL
        VVAR=VVAR-VDEL
        SIGSSM(3)=SIGSSM(2)
        INEW=2
        VNEW=VVAR
      ENDIF
 
C...Convert to relevant variables and find derived new limits.
      IF(IVAR.EQ.1) THEN
        VTAU=VNEW
        CALL PYKMAP(1,MTAU,VTAU)
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) CALL PYKLIM(4)
      ENDIF
      IF(IVAR.LE.2.AND.ISTSB.GE.3.AND.ISTSB.LE.5) THEN
        IF(IVAR.EQ.2) VTAUP=VNEW
        CALL PYKMAP(4,MTAUP,VTAUP)
      ENDIF
      IF(IVAR.LE.2) CALL PYKLIM(2)
      IF(IVAR.LE.3) THEN
        IF(IVAR.EQ.3) VYST=VNEW
        CALL PYKMAP(2,MYST,VYST)
        CALL PYKLIM(3)
      ENDIF
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) THEN
        IF(IVAR.EQ.4) VCTH=VNEW
        CALL PYKMAP(3,MCTH,VCTH)
      ENDIF
      IF(ISUB.EQ.96) VINT(25)=VINT(21)*(1.-VINT(23)**2)
 
C...Evaluate cross-section. Save new maximum. Final maximum.
      IF(ISTSB.NE.5) THEN
        CALL PYSIGH(NCHN,SIGS)
        IF(MWTXS.EQ.1) THEN
          CALL PYEVWT(WTXS)
          SIGS=WTXS*SIGS
        ENDIF
      ELSE
        SIGS=0.
        DO 380 IKIN3=1,MSTP(129)
        CALL PYKMAP(5,0,0.)
        IF(MINT(51).EQ.1) GOTO 380
        CALL PYSIGH(NCHN,SIGTMP)
        IF(MWTXS.EQ.1) THEN
          CALL PYEVWT(WTXS)
          SIGTMP=WTXS*SIGTMP
        ENDIF
        IF(SIGTMP.GT.SIGS) SIGS=SIGTMP
  380   CONTINUE
      ENDIF
      SIGSSM(INEW)=SIGS
      IF(SIGS.GT.SIGSAM) SIGSAM=SIGS
      IF(MSTP(122).GE.2) WRITE(MSTU(11),5800) IMAX,IVAR,MVAR,IMOV,
     &VNEW,VINT(21),VINT(22),VINT(23),VINT(26),SIGS
  390 CONTINUE
  400 CONTINUE
  410 CONTINUE
  420 CONTINUE
      IF(MSTP(121).EQ.1) SIGSAM=PARP(121)*SIGSAM
      XSEC(ISUB,1)=1.05*SIGSAM
  430 CONTINUE
      IF(MSTP(173).EQ.1.AND.ISUB.NE.96) XSEC(ISUB,1)=
     &PARP(174)*XSEC(ISUB,1)
      IF(ISUB.NE.96) XSEC(0,1)=XSEC(0,1)+XSEC(ISUB,1)
  440 CONTINUE
      MINT(51)=0
 
C...Print summary table.
      IF(NPOSI.EQ.0) THEN
        WRITE(MSTU(11),5900)
        STOP
      ENDIF
      IF(MSTP(122).GE.1) THEN
        WRITE(MSTU(11),6000)
        WRITE(MSTU(11),6100)
        DO 450 ISUB=1,200
        IF(MSUB(ISUB).NE.1.AND.ISUB.NE.96) GOTO 450
        IF(ISUB.EQ.96.AND.MINT(50).EQ.0) GOTO 450
        IF(ISUB.EQ.96.AND.MSUB(95).NE.1.AND.MSTP(81).LE.0) GOTO 450
        IF(ISUB.EQ.96.AND.MINT(49).EQ.0.AND.MSTP(131).EQ.0) GOTO 450
        IF(MSUB(95).EQ.1.AND.(ISUB.EQ.11.OR.ISUB.EQ.12.OR.ISUB.EQ.13.OR.
     &  ISUB.EQ.28.OR.ISUB.EQ.53.OR.ISUB.EQ.68)) GOTO 450
        WRITE(MSTU(11),6200) ISUB,PROC(ISUB),XSEC(ISUB,1)
  450   CONTINUE
        WRITE(MSTU(11),6300)
      ENDIF
 
C...Format statements for maximization results.
 5000 FORMAT(/1X,'Coefficient optimization and maximum search for ',
     &'subprocess no',I4/1X,'Coefficient modes     tau',10X,'y*',9X,
     &'cth',9X,'tau''',7X,'sigma')
 5100 FORMAT(1X,'Warning: requested subprocess ',I3,' has no allowed ',
     &'phase space.'/1X,'Process switched off!')
 5200 FORMAT(1X,4I4,F12.8,F12.6,F12.7,F12.8,1P,E12.4)
 5300 FORMAT(1X,'Warning: requested subprocess ',I3,' has vanishing ',
     &'cross-section.'/1X,'Process switched off!')
 5400 FORMAT(1X,'Coefficients of equation system to be solved for ',A4)
 5500 FORMAT(1X,1P,8E11.3)
 5600 FORMAT(1X,'Result for ',A4,':',7F9.4)
 5700 FORMAT(1X,'Maximum search for given coefficients'/2X,'MAX VAR ',
     &'MOD MOV   VNEW',7X,'tau',7X,'y*',8X,'cth',7X,'tau''',7X,'sigma')
 5800 FORMAT(1X,4I4,F8.4,F11.7,F9.3,F11.6,F11.7,1P,E12.4)
 5900 FORMAT(1X,'Error: no requested process has non-vanishing ',
     &'cross-section.'/1X,'Execution stopped!')
 6000 FORMAT(/1X,8('*'),1X,'PYMAXI: summary of differential ',
     &'cross-section maximum search',1X,8('*'))
 6100 FORMAT(/11X,58('=')/11X,'I',38X,'I',17X,'I'/11X,'I  ISUB  ',
     &'Subprocess name',15X,'I  Maximum value  I'/11X,'I',38X,'I',
     &17X,'I'/11X,58('=')/11X,'I',38X,'I',17X,'I')
 6200 FORMAT(11X,'I',2X,I3,3X,A28,2X,'I',2X,1P,E12.4,3X,'I')
 6300 FORMAT(11X,'I',38X,'I',17X,'I'/11X,58('='))
 
      RETURN
      END
