 
C*********************************************************************
 
      SUBROUTINE PYSTEL(X,Q2,XPEL)
 
C...Gives electron structure function.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/,/LUDAT2/
      SAVE /PYPARS/,/PYINT1/
      DIMENSION XPEL(-25:25),XPGA(-6:6),SXP(0:6)
 
C...Interface to PDFLIB.
      COMMON/W50513/XMIN,XMAX,Q2MIN,Q2MAX
      SAVE /W50513/
      DOUBLE PRECISION XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     &VALUE(20),XMIN,XMAX,Q2MIN,Q2MAX
      CHARACTER*20 PARM(20)
      DATA VALUE/20*0D0/,PARM/20*' '/
 
C...Some common constants.
      DO 100 KFL=-25,25
      XPEL(KFL)=0.
  100 CONTINUE
      AEM=PARU(101)
      PME=PMAS(11,1)
      XL=LOG(MAX(1E-10,X))
      X1L=LOG(MAX(1E-10,1.-X))
      HLE=LOG(MAX(3.,Q2/PME**2))
      HBE2=(AEM/PARU(1))*(HLE-1.)
 
C...Electron inside electron, see R. Kleiss et al., in Z physics at
C...LEP 1, CERN 89-08, p. 34
      IF(MSTP(59).LE.1) THEN
        HDE=1.+(AEM/PARU(1))*(1.5*HLE+1.289868)+(AEM/PARU(1))**2*
     &  (-2.164868*HLE**2+9.840808*HLE-10.130464)
        HEE=HBE2*(1.-X)**(HBE2-1.)*SQRT(MAX(0.,HDE))-
     &  0.5*HBE2*(1.+X)+HBE2**2/8.*((1.+X)*(-4.*X1L+3.*XL)-
     &  4.*XL/(1.-X)-5.-X)
      ELSE
        HEE=HBE2*(1.-X)**(HBE2-1.)*EXP(0.172784*HBE2)/PYGAMM(1.+HBE2)-
     &  0.5*HBE2*(1.+X)+HBE2**2/8.*((1.+X)*(-4.*X1L+3.*XL)-
     &  4.*XL/(1.-X)-5.-X)        
      ENDIF
      IF(X.GT.0.9999.AND.X.LE.0.999999) THEN
        HEE=HEE*100.**HBE2/(100.**HBE2-1.)
      ELSEIF(X.GT.0.999999) THEN
        HEE=0.
      ENDIF
      XPEL(11)=X*HEE
 
C...Photon and (transverse) W- inside electron.
      AEMP=ULALEM(PME*SQRT(MAX(0.,Q2)))/PARU(2)
      IF(MSTP(13).LE.1) THEN
        HLG=HLE
      ELSE
        HLG=LOG(MAX(1.,(PARP(13)/PME**2)*(1.-X)/X**2))
      ENDIF
      XPEL(22)=AEMP*HLG*(1.+(1.-X)**2)
      HLW=LOG(1.+Q2/PMAS(24,1)**2)/(4.*PARU(102))
      XPEL(-24)=AEMP*HLW*(1.+(1.-X)**2)
 
C...Electron or positron inside photon inside electron.
      IF(MSTP(12).EQ.1) THEN
        XFSEA=0.5*(AEMP*(HLE-1.))**2*(4./3.+X-X**2-4.*X**3/3.+
     &  2.*X*(1.+X)*XL)
        XPEL(11)=XPEL(11)+XFSEA
        XPEL(-11)=XFSEA
 
C...Initialize PDFLIB photon structure functions.
        IF(MSTP(56).EQ.2) THEN
          PARM(1)='NPTYPE'
          VALUE(1)=3
          PARM(2)='NGROUP'
          VALUE(2)=MSTP(55)/1000
          PARM(3)='NSET'
          VALUE(3)=MOD(MSTP(55),1000)
          IF(MINT(93).NE.3000000+MSTP(55)) THEN
            CALL PDFSET(PARM,VALUE)
            MINT(93)=3000000+MSTP(55)
          ENDIF
        ENDIF
 
C...Quarks and gluons inside photon inside electron:
C...numerical convolution required.
        DO 110 KFL=0,6
        SXP(KFL)=0.
  110   CONTINUE
        SUMXPP=0.
        ITER=-1
  120   ITER=ITER+1
        SUMXP=SUMXPP
        NSTP=2**(ITER-1)
        IF(ITER.EQ.0) NSTP=2
        DO 130 KFL=0,6
        SXP(KFL)=0.5*SXP(KFL)
  130   CONTINUE
        WTSTP=0.5/NSTP
        IF(ITER.EQ.0) WTSTP=0.5
C...Pick grid of x_{gamma} values logarithmically even.
        DO 150 ISTP=1,NSTP
        IF(ITER.EQ.0) THEN
          XLE=XL*(ISTP-1)
        ELSE
          XLE=XL*(ISTP-0.5)/NSTP
        ENDIF
        XE=MIN(0.999999,EXP(XLE))
        XG=MIN(0.999999,X/XE)
C...Evaluate photon inside electron structure function for convolution.
        XPGP=1.+(1.-XE)**2
        IF(MSTP(13).LE.1) THEN
          XPGP=XPGP*HLE
        ELSE
          XPGP=XPGP*LOG(MAX(1.,(PARP(13)/PME**2)*(1.-XE)/XE**2))
        ENDIF
C...Evaluate photon structure functions for convolution.
        IF(MSTP(56).EQ.1) THEN
          CALL PYSTGA(XG,Q2,XPGA)
          DO 140 KFL=0,5
          SXP(KFL)=SXP(KFL)+WTSTP*XPGP*XPGA(KFL)
  140     CONTINUE
        ELSEIF(MSTP(56).EQ.2) THEN
C...Call PDFLIB structure functions.
          XX=XG
          QQ=SQRT(MAX(0.,SNGL(Q2MIN),Q2))
          IF(MSTP(57).EQ.0) QQ=SQRT(Q2MIN)
          CALL STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
          SXP(0)=SXP(0)+WTSTP*XPGP*GLU
          SXP(1)=SXP(1)+WTSTP*XPGP*DNV
          SXP(2)=SXP(2)+WTSTP*XPGP*UPV
          SXP(3)=SXP(3)+WTSTP*XPGP*STR
          SXP(4)=SXP(4)+WTSTP*XPGP*CHM
          SXP(5)=SXP(5)+WTSTP*XPGP*BOT
          SXP(6)=SXP(6)+WTSTP*XPGP*TOP
        ENDIF
  150   CONTINUE
        SUMXPP=SXP(0)+2.*SXP(1)+2.*SXP(2)
        IF(ITER.LE.2.OR.(ITER.LE.7.AND.ABS(SUMXPP-SUMXP).GT.
     &  PARP(14)*(SUMXPP+SUMXP))) GOTO 120
 
C...Put convolution into output arrays.
        FCONV=AEMP*(-XL)
        XPEL(0)=FCONV*SXP(0)
        DO 160 KFL=1,6
        XPEL(KFL)=FCONV*SXP(KFL)
        XPEL(-KFL)=XPEL(KFL)
  160   CONTINUE
      ENDIF
 
      RETURN
      END
