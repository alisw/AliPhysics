 
C*********************************************************************
 
      SUBROUTINE PYXTOT
 
C...Parametrizes total, elastic and diffractive cross-sections
C...for different energies and beams. Donnachie-Landshoff for
C...total and Schuler-Sjostrand for elastic and diffractive.
C...Process code IPROC:
C...=  1 : p + p;
C...=  2 : pbar + p;
C...=  3 : pi+ + p;
C...=  4 : pi- + p;
C...=  5 : pi0 + p;
C...=  6 : phi + p;
C...=  7 : J/psi + p;
C...= 11 : rho + rho;
C...= 12 : rho + phi;
C...= 13 : rho + J/psi;
C...= 14 : phi + phi;
C...= 15 : phi + J/psi;
C...= 16 : J/psi + J/psi;
C...= 21 : gamma + p (DL);
C...= 22 : gamma + p (VDM).
C...= 23 : gamma + pi (DL);
C...= 24 : gamma + pi (VDM);
C...= 25 : gamma + gamma (DL);
C...= 26 : gamma + gamma (VDM).
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3)
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      SAVE /LUDAT1/
      SAVE /PYPARS/,/PYINT1/,/PYINT5/,/PYINT7/
      DIMENSION NPROC(30),XPAR(30),YPAR(30),IHADA(20),IHADB(20),
     &PMHAD(4),BHAD(4),BETP(4),IFITSD(20),IFITDD(20),CEFFS(10,10),
     &CEFFD(10,10),SIGTMP(6,0:5)
 
C...Common constants.
      DATA EPS/0.0808/, ETA/-0.4525/, ALP/0.25/, CRES/2./, PMRC/1.062/,
     &SMP/0.880/, FACEL/0.0511/, FACSD/0.0336/, FACDD/0.0084/
 
C...Number of multiple processes to be evaluated (= 0 : undefined).
      DATA NPROC/7*1,3*0,6*1,4*0,4*3,2*6,4*0/
C...X and Y parameters of sigmatot = X * s**epsilon + Y * s**(-eta).
      DATA XPAR/2*21.70,3*13.63,10.01,0.970,3*0.,
     &8.56,6.29,0.609,4.62,0.447,0.0434,4*0.,
     &0.0677,0.0534,0.0425,0.0335,2.11E-4,1.31E-4,4*0./
      DATA YPAR/56.08,98.39,27.56,36.02,31.79,-1.51,-0.146,3*0.,
     &18.02,-0.86,-0.083,0.041,-0.0039,0.00038,4*0.,
     &0.129,0.115,0.081,0.072,2.97E-4,2.36E-4,4*0./
 
C...Beam and target hadron class:
C...= 1 : p/n ; = 2 : pi/rho/omega; = 3 : phi; = 4 : J/psi.
      DATA IHADA/2*1,3*2,3,4,3*0,3*2,2*3,4,4*0/
      DATA IHADB/7*1,3*0,2,3,4,3,2*4,4*0/
C...Characteristic class masses, slope parameters, beta = sqrt(X).
      DATA PMHAD/0.938,0.770,1.020,3.097/
      DATA BHAD/2.3,1.4,1.4,0.23/
      DATA BETP/4.658,2.926,2.149,0.208/
 
C...Fitting constants used in parametrizations of diffractive results.
      DATA IFITSD/2*1,3*2,3,4,3*0,5,6,7,8,9,10,4*0/
      DATA IFITDD/2*1,3*2,3,4,3*0,5,6,7,8,9,10,4*0/
      DATA ((CEFFS(J1,J2),J2=1,10),J1=1,10)/
     & 0.213, 0.0, -0.47, 150., 0.213, 0.0, -0.47, 150., 0., 0.,
     & 0.213, 0.0, -0.47, 150., 0.267, 0.0, -0.47, 100., 0., 0.,
     & 0.213, 0.0, -0.47, 150., 0.232, 0.0, -0.47, 110., 0., 0.,
     & 0.213, 7.0, -0.55, 800., 0.115, 0.0, -0.47, 110., 0., 0.,
     & 0.267, 0.0, -0.46,  75., 0.267, 0.0, -0.46,  75., 0., 0.,
     & 0.232, 0.0, -0.46,  85., 0.267, 0.0, -0.48, 100., 0., 0.,
     & 0.115, 0.0, -0.50,  90., 0.267, 6.0, -0.56, 420., 0., 0.,
     & 0.232, 0.0, -0.48, 110., 0.232, 0.0, -0.48, 110., 0., 0.,
     & 0.115, 0.0, -0.52, 120., 0.232, 6.0, -0.56, 470., 0., 0.,
     & 0.115, 5.5, -0.58, 570., 0.115, 5.5, -0.58, 570., 0., 0./
      DATA ((CEFFD(J1,J2),J2=1,10),J1=1,10)/
     & 3.11, -7.34,  9.71, 0.068, -0.42, 1.31, -1.37,  35.0,  118., 0.,
     & 3.11, -7.10,  10.6, 0.073, -0.41, 1.17, -1.41,  31.6,   95., 0.,
     & 3.12, -7.43,  9.21, 0.067, -0.44, 1.41, -1.35,  36.5,  132., 0.,
     & 3.13, -8.18, -4.20, 0.056, -0.71, 3.12, -1.12,  55.2, 1298., 0.,
     & 3.11, -6.90,  11.4, 0.078, -0.40, 1.05, -1.40,  28.4,   78., 0.,
     & 3.11, -7.13,  10.0, 0.071, -0.41, 1.23, -1.34,  33.1,  105., 0.,
     & 3.12, -7.90, -1.49, 0.054, -0.64, 2.72, -1.13,  53.1,  995., 0.,
     & 3.11, -7.39,  8.22, 0.065, -0.44, 1.45, -1.36,  38.1,  148., 0.,
     & 3.18, -8.95, -3.37, 0.057, -0.76, 3.32, -1.12,  55.6, 1472., 0.,
     & 4.18, -29.2,  56.2, 0.074, -1.36, 6.67, -1.14, 116.2, 6532., 0./
 
C...Parameters. Combinations of the energy.
      AEM=PARU(101)
      PMTH=PARP(102)
      S=VINT(2)
      SRT=VINT(1)
      SEPS=S**EPS
      SETA=S**ETA
      SLOG=LOG(S)
 
C...Ratio of gamma/pi (for rescaling in structure functions).
      VINT(281)=(XPAR(22)*SEPS+YPAR(22)*SETA)/
     &(XPAR(5)*SEPS+YPAR(5)*SETA)
      IF(MINT(50).NE.1) RETURN
 
C...Order flavours of incoming particles: KF1 < KF2.
      IF(IABS(MINT(11)).LE.IABS(MINT(12))) THEN
        KF1=IABS(MINT(11))
        KF2=IABS(MINT(12))
        IORD=1
      ELSE
        KF1=IABS(MINT(12))
        KF2=IABS(MINT(11))
        IORD=2
      ENDIF
      ISGN12=ISIGN(1,MINT(11)*MINT(12))
 
C...Find process number (for lookup tables).
      IF(KF1.GT.1000) THEN
        IPROC=1
        IF(ISGN12.LT.0) IPROC=2
      ELSEIF(KF1.GT.100.AND.KF2.GT.1000) THEN
        IPROC=3
        IF(ISGN12.LT.0) IPROC=4
        IF(KF1.EQ.111) IPROC=5
      ELSEIF(KF1.GT.100) THEN
        IPROC=11
      ELSEIF(KF2.GT.1000) THEN
        IPROC=21
        IF(MINT(123).EQ.2) IPROC=22
      ELSEIF(KF2.GT.100) THEN
        IPROC=23
        IF(MINT(123).EQ.2) IPROC=24
      ELSE
        IPROC=25
        IF(MINT(123).EQ.2) IPROC=26
      ENDIF
 
C... Number of multiple processes to be stored; beam/target side.
      NPR=NPROC(IPROC)
      MINT(101)=1
      MINT(102)=1
      IF(NPR.EQ.3) THEN
        MINT(100+IORD)=4
      ELSEIF(NPR.EQ.6) THEN
        MINT(101)=4
        MINT(102)=4
      ENDIF
      N1=0
      IF(MINT(101).EQ.4) N1=4
      N2=0
      IF(MINT(102).EQ.4) N2=4
 
C...Do not do any more for user-set or undefined cross-sections.
      IF(MSTP(31).LE.0) RETURN
      IF(NPR.EQ.0) CALL LUERRM(26,
     &'(PYXTOT:) cross section for this process not yet implemented')
 
C...Parameters. Combinations of the energy.
      AEM=PARU(101)
      PMTH=PARP(102)
      S=VINT(2)
      SRT=VINT(1)
      SEPS=S**EPS
      SETA=S**ETA
      SLOG=LOG(S)
 
C...Loop over multiple processes (for VDM).
      DO 110 I=1,NPR
      IF(NPR.EQ.1) THEN
        IPR=IPROC
      ELSEIF(NPR.EQ.3) THEN
        IPR=I+4
        IF(KF2.LT.1000) IPR=I+10
      ELSEIF(NPR.EQ.6) THEN
        IPR=I+10
      ENDIF
 
C...Evaluate hadron species, mass, slope contribution and fit number.
      IHA=IHADA(IPR)
      IHB=IHADB(IPR)
      PMA=PMHAD(IHA)
      PMB=PMHAD(IHB)
      BHA=BHAD(IHA)
      BHB=BHAD(IHB)
      ISD=IFITSD(IPR)
      IDD=IFITDD(IPR)
 
C...Skip if energy too low relative to masses.
      DO 100 J=0,5
      SIGTMP(I,J)=0.
  100 CONTINUE
      IF(SRT.LT.1.5*(PMA+PMB)) GOTO 110
 
C...Total cross-section. Elastic slope parameter and cross-section.
      SIGTMP(I,0)=XPAR(IPR)*SEPS+YPAR(IPR)*SETA
      BEL=2.*BHA+2.*BHB+4.*SEPS-4.2
      SIGTMP(I,1)=FACEL*SIGTMP(I,0)**2/BEL
 
C...Diffractive scattering A + B -> X + B.
      BSD=2.*BHB
      SQML=(PMA+PMTH)**2
      SQMU=S*CEFFS(ISD,1)+CEFFS(ISD,2)
      SUM1=LOG((BSD+2.*ALP*LOG(S/SQML))/
     &(BSD+2.*ALP*LOG(S/SQMU)))/(2.*ALP)
      BXB=CEFFS(ISD,3)+CEFFS(ISD,4)/S
      SUM2=CRES*LOG(1.+((PMA+PMRC)/(PMA+PMTH))**2)/
     &(BSD+2.*ALP*LOG(S/((PMA+PMTH)*(PMA+PMRC)))+BXB)
      SIGTMP(I,2)=FACSD*XPAR(IPR)*BETP(IHB)*MAX(0.,SUM1+SUM2)
 
C...Diffractive scattering A + B -> A + X.
      BSD=2.*BHA
      SQML=(PMB+PMTH)**2
      SQMU=S*CEFFS(ISD,5)+CEFFS(ISD,6)
      SUM1=LOG((BSD+2.*ALP*LOG(S/SQML))/
     &(BSD+2.*ALP*LOG(S/SQMU)))/(2.*ALP)
      BAX=CEFFS(ISD,7)+CEFFS(ISD,8)/S
      SUM2=CRES*LOG(1.+((PMB+PMRC)/(PMB+PMTH))**2)/
     &(BSD+2.*ALP*LOG(S/((PMB+PMTH)*(PMB+PMRC)))+BAX)
      SIGTMP(I,3)=FACSD*XPAR(IPR)*BETP(IHA)*MAX(0.,SUM1+SUM2)
 
C...Order single diffractive correctly.
      IF(IORD.EQ.2) THEN
        SIGSAV=SIGTMP(I,2)
        SIGTMP(I,2)=SIGTMP(I,3)
        SIGTMP(I,3)=SIGSAV
      ENDIF
 
C...Double diffractive scattering A + B -> X1 + X2.
      YEFF=LOG(S*SMP/((PMA+PMTH)*(PMB+PMTH))**2)
      DEFF=CEFFD(IDD,1)+CEFFD(IDD,2)/SLOG+CEFFD(IDD,3)/SLOG**2
      SUM1=DEFF+YEFF*(LOG(MAX(1E-10,YEFF/DEFF))-1.)/(2.*ALP)
      IF(YEFF.LE.0) SUM1=0.
      SQMU=S*(CEFFD(IDD,4)+CEFFD(IDD,5)/SLOG+CEFFD(IDD,6)/SLOG**2)
      SLUP=LOG(MAX(1.1,S/(ALP*(PMA+PMTH)**2*(PMB+PMTH)*(PMB+PMRC))))
      SLDN=LOG(MAX(1.1,S/(ALP*SQMU*(PMB+PMTH)*(PMB+PMRC))))
      SUM2=CRES*LOG(1.+((PMB+PMRC)/(PMB+PMTH))**2)*LOG(SLUP/SLDN)/
     &(2.*ALP)
      SLUP=LOG(MAX(1.1,S/(ALP*(PMB+PMTH)**2*(PMA+PMTH)*(PMA+PMRC))))
      SLDN=LOG(MAX(1.1,S/(ALP*SQMU*(PMA+PMTH)*(PMA+PMRC))))
      SUM3=CRES*LOG(1.+((PMA+PMRC)/(PMA+PMTH))**2)*LOG(SLUP/SLDN)/
     &(2.*ALP)
      BXX=CEFFD(IDD,7)+CEFFD(IDD,8)/SRT+CEFFD(IDD,9)/S
      SLRR=LOG(S/(ALP*(PMA+PMTH)*(PMA+PMRC)*(PMB+PMTH)*(PMB*PMRC)))
      SUM4=CRES**2*LOG(1.+((PMA+PMRC)/(PMA+PMTH))**2)*
     &LOG(1.+((PMB+PMRC)/(PMB+PMTH))**2)/MAX(0.1,2.*ALP*SLRR+BXX)
      SIGTMP(I,4)=FACDD*XPAR(IPR)*MAX(0.,SUM1+SUM2+SUM3+SUM4)
 
C...Non-diffractive by unitarity.
      SIGTMP(I,5)=SIGTMP(I,0)-SIGTMP(I,1)-SIGTMP(I,2)-SIGTMP(I,3)-
     &SIGTMP(I,4)
  110 CONTINUE
 
C...Put temporary results in output array: only one process.
      IF(MINT(101).EQ.1.AND.MINT(102).EQ.1) THEN
        DO 120 J=0,5
        SIGT(0,0,J)=SIGTMP(1,J)
  120   CONTINUE
 
C...Beam multiple processes.
      ELSEIF(MINT(101).EQ.4.AND.MINT(102).EQ.1) THEN
        DO 140 I=1,4
        CONV=AEM/PARP(160+I)
        I1=MAX(1,I-1)
        DO 130 J=0,5
        SIGT(I,0,J)=CONV*SIGTMP(I1,J)
  130   CONTINUE
  140   CONTINUE
        DO 150 J=0,5
        SIGT(0,0,J)=SIGT(1,0,J)+SIGT(2,0,J)+SIGT(3,0,J)+SIGT(4,0,J)
  150   CONTINUE
 
C...Target multiple processes.
      ELSEIF(MINT(101).EQ.1.AND.MINT(102).EQ.4) THEN
        DO 170 I=1,4
        CONV=AEM/PARP(160+I)
        IV=MAX(1,I-1)
        DO 160 J=0,5
        SIGT(0,I,J)=CONV*SIGTMP(IV,J)
  160   CONTINUE
  170   CONTINUE
        DO 180 J=0,5
        SIGT(0,0,J)=SIGT(0,1,J)+SIGT(0,2,J)+SIGT(0,3,J)+SIGT(0,4,J)
  180   CONTINUE
 
C...Both beam and target multiple processes.
      ELSE
        DO 210 I1=1,4
        DO 200 I2=1,4
        CONV=AEM**2/(PARP(160+I1)*PARP(160+I2))
        IF(I1.LE.2) THEN
          IV=MAX(1,I2-1)
        ELSEIF(I2.LE.2) THEN
          IV=MAX(1,I1-1)
        ELSEIF(I1.EQ.I2) THEN
          IV=2*I1-2
        ELSE
          IV=5
        ENDIF
        DO 190 J=0,5
        JV=J
        IF(I2.GT.I1.AND.(J.EQ.2.OR.J.EQ.3)) JV=5-J
        SIGT(I1,I2,J)=CONV*SIGTMP(IV,JV)
  190   CONTINUE
  200   CONTINUE
  210   CONTINUE
        DO 230 J=0,5
        DO 220 I=1,4
        SIGT(I,0,J)=SIGT(I,1,J)+SIGT(I,2,J)+SIGT(I,3,J)+SIGT(I,4,J)
        SIGT(0,I,J)=SIGT(1,I,J)+SIGT(2,I,J)+SIGT(3,I,J)+SIGT(4,I,J)
  220   CONTINUE
        SIGT(0,0,J)=SIGT(1,0,J)+SIGT(2,0,J)+SIGT(3,0,J)+SIGT(4,0,J)
  230   CONTINUE
      ENDIF
 
C...Scale up uniformly for Donnachie-Landshoff parametrization.
      IF(IPROC.EQ.21.OR.IPROC.EQ.23.OR.IPROC.EQ.25) THEN
        RFAC=(XPAR(IPROC)*SEPS+YPAR(IPROC)*SETA)/SIGT(0,0,0)
        DO 260 I1=0,N1
        DO 250 I2=0,N2
        DO 240 J=0,5
        SIGT(I1,I2,J)=RFAC*SIGT(I1,I2,J)
  240   CONTINUE
  250   CONTINUE
  260   CONTINUE
      ENDIF
 
      RETURN
      END
