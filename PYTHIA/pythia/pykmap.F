 
C*********************************************************************
 
      SUBROUTINE PYKMAP(IVAR,MVAR,VVAR)
 
C...Maps a uniform distribution into a distribution of a kinematical
C...variable according to one of the possibilities allowed. It is
C...assumed that kinematical limits have been set by a PYKLIM call.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/
 
C...Convert VVAR to tau variable.
      ISUB=MINT(1)
      ISTSB=ISET(ISUB)
      IF(IVAR.EQ.1) THEN
        TAUMIN=VINT(11)
        TAUMAX=VINT(31)
        IF(MVAR.EQ.3.OR.MVAR.EQ.4) THEN
          TAURE=VINT(73)
          GAMRE=VINT(74)
        ELSEIF(MVAR.EQ.5.OR.MVAR.EQ.6) THEN
          TAURE=VINT(75)
          GAMRE=VINT(76)
        ENDIF
        IF(MINT(47).EQ.1.AND.(ISTSB.EQ.1.OR.ISTSB.EQ.2.OR.ISTSB.EQ.6))
     &  THEN
          TAU=1.
        ELSEIF(MVAR.EQ.1) THEN
          TAU=TAUMIN*(TAUMAX/TAUMIN)**VVAR
        ELSEIF(MVAR.EQ.2) THEN
          TAU=TAUMAX*TAUMIN/(TAUMIN+(TAUMAX-TAUMIN)*VVAR)
        ELSEIF(MVAR.EQ.3.OR.MVAR.EQ.5) THEN
          RATGEN=(TAURE+TAUMAX)/(TAURE+TAUMIN)*TAUMIN/TAUMAX
          TAU=TAURE*TAUMIN/((TAURE+TAUMIN)*RATGEN**VVAR-TAUMIN)
        ELSEIF(MVAR.EQ.4.OR.MVAR.EQ.6) THEN
          AUPP=ATAN((TAUMAX-TAURE)/GAMRE)
          ALOW=ATAN((TAUMIN-TAURE)/GAMRE)
          TAU=TAURE+GAMRE*TAN(ALOW+(AUPP-ALOW)*VVAR)
        ELSE
          AUPP=LOG(MAX(2E-6,1.-TAUMAX))
          ALOW=LOG(MAX(2E-6,1.-TAUMIN))
          TAU=1.-EXP(AUPP+VVAR*(ALOW-AUPP))
        ENDIF
        VINT(21)=MIN(TAUMAX,MAX(TAUMIN,TAU))
 
C...Convert VVAR to y* variable.
      ELSEIF(IVAR.EQ.2) THEN
        YSTMIN=VINT(12)
        YSTMAX=VINT(32)
        TAUE=VINT(21)
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=VINT(26)
        IF(MINT(47).EQ.1) THEN
          YST=0.
        ELSEIF(MINT(47).EQ.2) THEN
          YST=-0.5*LOG(TAUE)
        ELSEIF(MINT(47).EQ.3) THEN
          YST=0.5*LOG(TAUE)
        ELSEIF(MVAR.EQ.1) THEN
          YST=YSTMIN+(YSTMAX-YSTMIN)*SQRT(VVAR)
        ELSEIF(MVAR.EQ.2) THEN
          YST=YSTMAX-(YSTMAX-YSTMIN)*SQRT(1.-VVAR)
        ELSEIF(MVAR.EQ.3) THEN
          AUPP=ATAN(EXP(YSTMAX))
          ALOW=ATAN(EXP(YSTMIN))
          YST=LOG(TAN(ALOW+(AUPP-ALOW)*VVAR))
        ELSEIF(MVAR.EQ.4) THEN
          YST0=-0.5*LOG(TAUE)
          AUPP=LOG(MAX(1E-6,EXP(YST0-YSTMIN)-1.))
          ALOW=LOG(MAX(1E-6,EXP(YST0-YSTMAX)-1.))
          YST=YST0-LOG(1.+EXP(ALOW+VVAR*(AUPP-ALOW)))
        ELSE
          YST0=-0.5*LOG(TAUE)
          AUPP=LOG(MAX(1E-6,EXP(YST0+YSTMIN)-1.))
          ALOW=LOG(MAX(1E-6,EXP(YST0+YSTMAX)-1.))
          YST=LOG(1.+EXP(AUPP+VVAR*(ALOW-AUPP)))-YST0
        ENDIF
        VINT(22)=MIN(YSTMAX,MAX(YSTMIN,YST))
 
C...Convert VVAR to cos(theta-hat) variable.
      ELSEIF(IVAR.EQ.3) THEN
        RM34=MAX(1E-20,2.*VINT(63)*VINT(64)/(VINT(21)*VINT(2))**2)
        RSQM=1.+RM34
        IF(2.*VINT(71)**2/(VINT(21)*VINT(2)).LT.0.0001) RM34=MAX(RM34,
     &  2.*VINT(71)**2/(VINT(21)*VINT(2)))
        CTNMIN=VINT(13)
        CTNMAX=VINT(33)
        CTPMIN=VINT(14)
        CTPMAX=VINT(34)
        IF(MVAR.EQ.1) THEN
          ANEG=CTNMAX-CTNMIN
          APOS=CTPMAX-CTPMIN
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN
            VCTN=VVAR*(ANEG+APOS)/ANEG
            CTH=CTNMIN+(CTNMAX-CTNMIN)*VCTN
          ELSE
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS
            CTH=CTPMIN+(CTPMAX-CTPMIN)*VCTP
          ENDIF
        ELSEIF(MVAR.EQ.2) THEN
          RMNMIN=MAX(RM34,RSQM-CTNMIN)
          RMNMAX=MAX(RM34,RSQM-CTNMAX)
          RMPMIN=MAX(RM34,RSQM-CTPMIN)
          RMPMAX=MAX(RM34,RSQM-CTPMAX)
          ANEG=LOG(RMNMIN/RMNMAX)
          APOS=LOG(RMPMIN/RMPMAX)
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN
            VCTN=VVAR*(ANEG+APOS)/ANEG
            CTH=RSQM-RMNMIN*(RMNMAX/RMNMIN)**VCTN
          ELSE
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS
            CTH=RSQM-RMPMIN*(RMPMAX/RMPMIN)**VCTP
          ENDIF
        ELSEIF(MVAR.EQ.3) THEN
          RMNMIN=MAX(RM34,RSQM+CTNMIN)
          RMNMAX=MAX(RM34,RSQM+CTNMAX)
          RMPMIN=MAX(RM34,RSQM+CTPMIN)
          RMPMAX=MAX(RM34,RSQM+CTPMAX)
          ANEG=LOG(RMNMAX/RMNMIN)
          APOS=LOG(RMPMAX/RMPMIN)
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN
            VCTN=VVAR*(ANEG+APOS)/ANEG
            CTH=RMNMIN*(RMNMAX/RMNMIN)**VCTN-RSQM
          ELSE
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS
            CTH=RMPMIN*(RMPMAX/RMPMIN)**VCTP-RSQM
          ENDIF
        ELSEIF(MVAR.EQ.4) THEN
          RMNMIN=MAX(RM34,RSQM-CTNMIN)
          RMNMAX=MAX(RM34,RSQM-CTNMAX)
          RMPMIN=MAX(RM34,RSQM-CTPMIN)
          RMPMAX=MAX(RM34,RSQM-CTPMAX)
          ANEG=1./RMNMAX-1./RMNMIN
          APOS=1./RMPMAX-1./RMPMIN
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN
            VCTN=VVAR*(ANEG+APOS)/ANEG
            CTH=RSQM-1./(1./RMNMIN+ANEG*VCTN)
          ELSE
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS
            CTH=RSQM-1./(1./RMPMIN+APOS*VCTP)
          ENDIF
        ELSEIF(MVAR.EQ.5) THEN
          RMNMIN=MAX(RM34,RSQM+CTNMIN)
          RMNMAX=MAX(RM34,RSQM+CTNMAX)
          RMPMIN=MAX(RM34,RSQM+CTPMIN)
          RMPMAX=MAX(RM34,RSQM+CTPMAX)
          ANEG=1./RMNMIN-1./RMNMAX
          APOS=1./RMPMIN-1./RMPMAX
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN
            VCTN=VVAR*(ANEG+APOS)/ANEG
            CTH=1./(1./RMNMIN-ANEG*VCTN)-RSQM
          ELSE
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS
            CTH=1./(1./RMPMIN-APOS*VCTP)-RSQM
          ENDIF
        ENDIF
        IF(CTH.LT.0.) CTH=MIN(CTNMAX,MAX(CTNMIN,CTH))
        IF(CTH.GT.0.) CTH=MIN(CTPMAX,MAX(CTPMIN,CTH))
        VINT(23)=CTH
 
C...Convert VVAR to tau' variable.
      ELSEIF(IVAR.EQ.4) THEN
        TAU=VINT(21)
        TAUPMN=VINT(16)
        TAUPMX=VINT(36)
        IF(MINT(47).EQ.1) THEN
          TAUP=1.
        ELSEIF(MVAR.EQ.1) THEN
          TAUP=TAUPMN*(TAUPMX/TAUPMN)**VVAR
        ELSEIF(MVAR.EQ.2) THEN
          AUPP=(1.-TAU/TAUPMX)**4
          ALOW=(1.-TAU/TAUPMN)**4
          TAUP=TAU/MAX(1E-7,1.-(ALOW+(AUPP-ALOW)*VVAR)**0.25)
        ELSE
          AUPP=LOG(MAX(2E-6,1.-TAUPMX))
          ALOW=LOG(MAX(2E-6,1.-TAUPMN))
          TAUP=1.-EXP(AUPP+VVAR*(ALOW-AUPP))
        ENDIF
        VINT(26)=MIN(TAUPMX,MAX(TAUPMN,TAUP))
 
C...Selection of extra variables needed in 2 -> 3 process:
C...pT1, pT2, phi1, phi2, y3 for three outgoing particles.
C...Since no options are available, the functions of PYKLIM
C...and PYKMAP are joint for these choices.
      ELSEIF(IVAR.EQ.5) THEN
 
C...Read out total energy and particle masses.
        MINT(51)=0
        MPTPK=1
        IF(ISUB.EQ.123.OR.ISUB.EQ.124.OR.ISUB.EQ.173.OR.ISUB.EQ.174
     &  .OR.ISUB.EQ.178.OR.ISUB.EQ.179) MPTPK=2
        SHP=VINT(26)*VINT(2)
        SHPR=SQRT(SHP)
        PM1=VINT(201)
        PM2=VINT(206)
        PM3=SQRT(VINT(21))*VINT(1)
        IF(PM1+PM2+PM3.GT.0.9999*SHPR) THEN
          MINT(51)=1
          RETURN
        ENDIF
        PMRS1=VINT(204)**2
        PMRS2=VINT(209)**2
 
C...Specify coefficients of pT choice; upper and lower limits.
        IF(MPTPK.EQ.1) THEN
          HWT1=0.4
          HWT2=0.4
        ELSE
          HWT1=0.05
          HWT2=0.05
        ENDIF
        HWT3=1.-HWT1-HWT2
        PTSMX1=((SHP-PM1**2-(PM2+PM3)**2)**2-(2.*PM1*(PM2+PM3))**2)/
     &  (4.*SHP)
        IF(CKIN(52).GT.0.) PTSMX1=MIN(PTSMX1,CKIN(52)**2)
        PTSMN1=CKIN(51)**2
        PTSMX2=((SHP-PM2**2-(PM1+PM3)**2)**2-(2.*PM2*(PM1+PM3))**2)/
     &  (4.*SHP)
        IF(CKIN(54).GT.0.) PTSMX2=MIN(PTSMX2,CKIN(54)**2)
        PTSMN2=CKIN(53)**2
 
C...Select transverse momenta according to
C...dp_T^2 * (a + b/(M^2 + p_T^2) + c/(M^2 + p_T^2)^2).
        HMX=PMRS1+PTSMX1
        HMN=PMRS1+PTSMN1
        IF(HMX.LT.1.0001*HMN) THEN
          MINT(51)=1
          RETURN
        ENDIF
        HDE=PTSMX1-PTSMN1
        RPT=RLU(0)
        IF(RPT.LT.HWT1) THEN
          PTS1=PTSMN1+RLU(0)*HDE
        ELSEIF(RPT.LT.HWT1+HWT2) THEN
          PTS1=MAX(PTSMN1,HMN*(HMX/HMN)**RLU(0)-PMRS1)
        ELSE
          PTS1=MAX(PTSMN1,HMN*HMX/(HMN+RLU(0)*HDE)-PMRS1)
        ENDIF
        WTPTS1=HDE/(HWT1+HWT2*HDE/(LOG(HMX/HMN)*(PMRS1+PTS1))+
     &  HWT3*HMN*HMX/(PMRS1+PTS1)**2)
        HMX=PMRS2+PTSMX2
        HMN=PMRS2+PTSMN2
        IF(HMX.LT.1.0001*HMN) THEN
          MINT(51)=1
          RETURN
        ENDIF
        HDE=PTSMX2-PTSMN2
        RPT=RLU(0)
        IF(RPT.LT.HWT1) THEN
          PTS2=PTSMN2+RLU(0)*HDE
        ELSEIF(RPT.LT.HWT1+HWT2) THEN
          PTS2=MAX(PTSMN2,HMN*(HMX/HMN)**RLU(0)-PMRS2)
        ELSE
          PTS2=MAX(PTSMN2,HMN*HMX/(HMN+RLU(0)*HDE)-PMRS2)
        ENDIF
        WTPTS2=HDE/(HWT1+HWT2*HDE/(LOG(HMX/HMN)*(PMRS2+PTS2))+
     &  HWT3*HMN*HMX/(PMRS2+PTS2)**2)
 
C...Select azimuthal angles and check pT choice.
        PHI1=PARU(2)*RLU(0)
        PHI2=PARU(2)*RLU(0)
        PHIR=PHI2-PHI1
        PTS3=MAX(0.,PTS1+PTS2+2.*SQRT(PTS1*PTS2)*COS(PHIR))
        IF(PTS3.LT.CKIN(55)**2.OR.(CKIN(56).GT.0..AND.PTS3.GT.
     &  CKIN(56)**2)) THEN
          MINT(51)=1
          RETURN
        ENDIF
 
C...Calculate transverse masses and check phase space not closed.
        PMS1=PM1**2+PTS1
        PMS2=PM2**2+PTS2
        PMS3=PM3**2+PTS3
        PMT1=SQRT(PMS1)
        PMT2=SQRT(PMS2)
        PMT3=SQRT(PMS3)
        PM12=(PMT1+PMT2)**2
        IF(PMT1+PMT2+PMT3.GT.0.9999*SHPR) THEN
          MINT(51)=1
          RETURN
        ENDIF
 
C...Select rapidity for particle 3 and check phase space not closed.
        Y3MAX=LOG((SHP+PMS3-PM12+SQRT(MAX(0.,(SHP-PMS3-PM12)**2-
     &  4.*PMS3*PM12)))/(2.*SHPR*PMT3))
        IF(Y3MAX.LT.1E-6) THEN
          MINT(51)=1
          RETURN
        ENDIF
        Y3=(2.*RLU(0)-1.)*0.999999*Y3MAX
        PZ3=PMT3*SINH(Y3)
        PE3=PMT3*COSH(Y3)
 
C...Find momentum transfers in two mirror solutions (in 1-2 frame).
        PZ12=-PZ3
        PE12=SHPR-PE3
        PMS12=PE12**2-PZ12**2
        SQL12=SQRT(MAX(0.,(PMS12-PMS1-PMS2)**2-4.*PMS1*PMS2))
        IF(SQL12.LT.1E-6*SHP) THEN
          MINT(51)=1
          RETURN
        ENDIF
        PMM1=PMS12+PMS1-PMS2
        PMM2=PMS12+PMS2-PMS1
        TFAC=-SHPR/(2.*PMS12)
        T1P=TFAC*(PE12-PZ12)*(PMM1-SQL12)
        T1N=TFAC*(PE12-PZ12)*(PMM1+SQL12)
        T2P=TFAC*(PE12+PZ12)*(PMM2-SQL12)
        T2N=TFAC*(PE12+PZ12)*(PMM2+SQL12)
 
C...Construct relative mirror weights and make choice.
        IF(MPTPK.EQ.1) THEN
          WTPU=1.
          WTNU=1.
        ELSE
          WTPU=1./((T1P-PMRS1)*(T2P-PMRS2))**2
          WTNU=1./((T1N-PMRS1)*(T2N-PMRS2))**2
        ENDIF
        WTP=WTPU/(WTPU+WTNU)
        WTN=WTNU/(WTPU+WTNU)
        EPS=1.
        IF(WTN.GT.RLU(0)) EPS=-1.
 
C...Store result of variable choice and associated weights.
        VINT(202)=PTS1
        VINT(207)=PTS2
        VINT(203)=PHI1
        VINT(208)=PHI2
        VINT(205)=WTPTS1
        VINT(210)=WTPTS2
        VINT(211)=Y3
        VINT(212)=Y3MAX
        VINT(213)=EPS
        IF(EPS.GT.0.) THEN
          VINT(214)=1./WTP
          VINT(215)=T1P
          VINT(216)=T2P
        ELSE
          VINT(214)=1./WTN
          VINT(215)=T1N
          VINT(216)=T2N
        ENDIF
        VINT(217)=-0.5*TFAC*(PE12-PZ12)*(PMM2+EPS*SQL12)
        VINT(218)=-0.5*TFAC*(PE12+PZ12)*(PMM1+EPS*SQL12)
        VINT(219)=0.5*(PMS12-PTS3)
        VINT(220)=SQL12
      ENDIF
 
      RETURN
      END
