C*********************************************************************
C*********************************************************************
C*                                                                  **
C*                                                 December 1993    **
C*                                                                  **
C*           The Lund Monte Carlo for Hadronic Processes            **
C*                                                                  **
C*                        PYTHIA version 5.7                        **
C*                                                                  **
C*                        Torbjorn Sjostrand                        **
C*                Department of theoretical physics 2               **
C*                        University of Lund                        **
C*               Solvegatan 14A, S-223 62 Lund, Sweden              **
C*                    E-mail torbjorn@thep.lu.se                    **
C*                    phone +46 - 46 - 222 48 16                    ** 
C*                                                                  **
C*         Several parts are written by Hans-Uno Bengtsson          **
C*     CTEQ 2 parton distributions are by the CTEQ collaboration    **
C*   SaS photon parton distributions together with Gerhard Schuler  **
C*    g + g -> Z + b + bbar matrix element code by Ronald Kleiss    **
C*     g + g and q + qbar -> t + tbar + H code by Zoltan Kunszt     **
C*                                                                  **
C*   The latest program version and documentation is found on WWW   **
C*         http://thep.lu.se/tf2/staff/torbjorn/Welcome.html        **
C*                                                                  **
C*        Copyright Torbjorn Sjostrand and CERN, Geneva 1993        **
C*                                                                  **
C*********************************************************************
C*********************************************************************
C                                                                    *
C  List of subprograms in order of appearance, with main purpose     *
C  (S = subroutine, F = function, B = block data)                    *
C                                                                    *
C  S   PYINIT   to administer the initialization procedure           *
C  S   PYEVNT   to administer the generation of an event             *
C  S   PYSTAT   to print cross-section and other information         *
C  S   PYINRE   to initialize treatment of resonances                *
C  S   PYINBM   to read in beam, target and frame choices            *
C  S   PYINKI   to initialize kinematics of incoming particles       *
C  S   PYINPR   to set up the selection of included processes        *
C  S   PYXTOT   to give total, elastic and diffractive cross-sect.   *
C  S   PYMAXI   to find differential cross-section maxima            *
C  S   PYPILE   to select multiplicity of pileup events              *
C  S   PYSAVE   to save alternatives for gamma-p and gamma-gamma     *
C  S   PYRAND   to select subprocess and kinematics for event        *
C  S   PYSCAT   to set up kinematics and colour flow of event        *
C  S   PYSSPA   to simulate initial state spacelike showers          *
C  S   PYRESD   to perform resonance decays                          *
C  S   PYMULT   to generate multiple interactions                    *
C  S   PYREMN   to add on target remnants                            *
C  S   PYDIFF   to set up kinematics for diffractive events          *
C  S   PYDOCU   to compute cross-sections and handle documentation   *
C  S   PYFRAM   to perform boosts between different frames           *
C  S   PYWIDT   to calculate full and partial widths of resonances   *
C  S   PYOFSH   to calculate partial width into off-shell channels   *
C  S   PYKLIM   to calculate borders of allowed kinematical region   *
C  S   PYKMAP   to construct value of kinematical variable           *
C  S   PYSIGH   to calculate differential cross-sections             *
C  S   PYSTFU   to evaluate structure functions                      *
C  S   PYSTFL   to evaluate structure functions at low x and Q^2     *
C  S   PYSTEL   to evaluate electron structure function              *
C  S   PYSTGA   to evaluate photon structure function (generic)      *
C  S   PYGGAM   to evaluate photon structure function (SaS sets)     *
C  S   PYGVMD   to evaluate VMD part of photon structure functions   *
C  S   PYGANO   to evaluate anomalous part of photon str. func.      *
C  S   PYGBEH   to evaluate Bethe-Heitler part of photon str. func.  *
C  S   PYGDIR   to evaluate direct contribution to photon str. func. *
C  S   PYSTPI   to evaluate pion structure function                  *
C  S   PYSTPR   to evaluate proton structure function                *
C  F   PYCTQ2   to evaluate the CTEQ 2 proton structure function     *
C  F   PYHFTH   to evaluate threshold factor for heavy flavour       *
C  S   PYSPLI   to find flavours left in hadron when one removed     *
C  F   PYGAMM   to evaluate ordinary Gamma function Gamma(x)         *
C  S   PYWAUX   to evaluate auxiliary functions W1(s) and W2(s)      *
C  S   PYI3AU   to evaluate auxiliary function I3(s,t,u,v)           *
C  F   PYSPEN   to evaluate Spence (dilogarithm) function Sp(x)      *
C  S   PYQQBH   to evaluate matrix element for g + g -> Q + Q~ + H   *
C  S   PYTEST   to test the proper functioning of the package        *
C  B   PYDATA   to contain all default values                        *
C  S   PYKCUT   to provide dummy routine for user kinematical cuts   *
C  S   PYEVWT   to provide dummy routine for weighting events        *
C  S   PYUPIN   to initialize a user process                         *
C  S   PYUPEV   to generate a user process event (dummy routine)     *
C  S   PDFSET   dummy routine to be removed when using PDFLIB        *
C  S   STRUCTM  dummy routine to be removed when using PDFLIB        *
C                                                                    *
C*********************************************************************
 
      SUBROUTINE PYINIT(FRAME,BEAM,TARGET,WIN)
 
C...Initializes the generation procedure; finds maxima of the
C...differential cross-sections to be used for weighting.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/LUDAT4/CHAF(500)
      CHARACTER CHAF*8
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3)
      COMMON/PYINT9/DXSEC(0:200)
      DOUBLE PRECISION DXSEC
      SAVE /LUDAT1/,/LUDAT2/,/LUDAT3/,/LUDAT4/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT5/,/PYINT9/
      DIMENSION ALAMIN(20),NFIN(20)
      CHARACTER*(*) FRAME,BEAM,TARGET
      CHARACTER CHFRAM*8,CHBEAM*8,CHTARG*8,CHLH(2)*6
 
C...Interface to PDFLIB.
      COMMON/W50512/QCDL4,QCDL5
      SAVE /W50512/
      DOUBLE PRECISION VALUE(20),QCDL4,QCDL5
      CHARACTER*20 PARM(20)
      DATA VALUE/20*0D0/,PARM/20*' '/
 
C...Data:Lambda and n_f values for structure functions; months.
      DATA ALAMIN/0.20,0.29,0.20,0.40,0.213,0.208,0.208,0.322,
     &0.190,0.235,10*0.2/,NFIN/20*4/
      DATA CHLH/'lepton','hadron'/
 
C...Reset MINT and VINT arrays. Write headers.
      DO 100 J=1,400
      MINT(J)=0
      VINT(J)=0.
  100 CONTINUE
      IF(MSTU(12).GE.1) CALL LULIST(0)
      IF(MSTP(122).GE.1) WRITE(MSTU(11),5100)
 
C...Maximum 4 generations; set maximum number of allowed flavours.
      MSTP(1)=MIN(4,MSTP(1))
      MSTU(114)=MIN(MSTU(114),2*MSTP(1))
      MSTP(58)=MIN(MSTP(58),2*MSTP(1))
 
C...Sum up Cabibbo-Kobayashi-Maskawa factors for each quark/lepton.
      DO 120 I=-20,20
      VINT(180+I)=0.
      IA=IABS(I)
      IF(IA.GE.1.AND.IA.LE.2*MSTP(1)) THEN
        DO 110 J=1,MSTP(1)
        IB=2*J-1+MOD(IA,2)
        IPM=(5-ISIGN(1,I))/2
        IDC=J+MDCY(IA,2)+2
        IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.IPM) VINT(180+I)=
     &  VINT(180+I)+VCKM((IA+1)/2,(IB+1)/2)
  110   CONTINUE
      ELSEIF(IA.GE.11.AND.IA.LE.10+2*MSTP(1)) THEN
        VINT(180+I)=1.
      ENDIF
  120 CONTINUE
 
C...Initialize structure functions: PDFLIB.
      IF(MSTP(52).EQ.2) THEN
        PARM(1)='NPTYPE'
        VALUE(1)=1
        PARM(2)='NGROUP'
        VALUE(2)=MSTP(51)/1000
        PARM(3)='NSET'
        VALUE(3)=MOD(MSTP(51),1000)
        PARM(4)='TMAS'
        VALUE(4)=PMAS(6,1)
        CALL PDFSET(PARM,VALUE)
        MINT(93)=1000000+MSTP(51)
      ENDIF
 
C...Choose Lambda value to use in alpha-strong.
      MSTU(111)=MSTP(2)
      IF(MSTP(3).GE.2) THEN
        ALAM=0.2
        NF=4
        IF(MSTP(52).EQ.1.AND.MSTP(51).GE.1.AND.MSTP(51).LE.10) THEN
          ALAM=ALAMIN(MSTP(51))
          NF=NFIN(MSTP(51))
        ELSEIF(MSTP(52).EQ.2) THEN
          ALAM=QCDL4
          NF=4
        ENDIF
        PARP(1)=ALAM
        PARP(61)=ALAM
        PARP(72)=ALAM
        PARU(112)=ALAM
        MSTU(112)=NF
        IF(MSTP(3).EQ.3) PARJ(81)=ALAM      
      ENDIF
 
C...Initialize widths and partial widths for resonances.
      CALL PYINRE
 
C...Identify beam and target particles and frame of process.
      CHFRAM=FRAME//' '
      CHBEAM=BEAM//' '
      CHTARG=TARGET//' '
      CALL PYINBM(CHFRAM,CHBEAM,CHTARG,WIN)
      IF(MINT(65).EQ.1) GOTO 170
 
C...For gamma-p or gamma-gamma allow many (3 or 6) alternatives.
C...For e-gamma allow 2 alternatives.
      MINT(121)=1
      MINT(123)=MSTP(14)
      IF(MSTP(14).EQ.10.AND.(MSEL.EQ.1.OR.MSEL.EQ.2)) THEN
        IF((MINT(11).EQ.22.OR.MINT(12).EQ.22).AND.
     &  (IABS(MINT(11)).GE.28.OR.IABS(MINT(12)).GE.28)) MINT(121)=3
        IF(MINT(11).EQ.22.AND.MINT(12).EQ.22) MINT(121)=6
        IF((MINT(11).EQ.22.OR.MINT(12).EQ.22).AND.
     &  (IABS(MINT(11)).EQ.11.OR.IABS(MINT(12)).EQ.11)) MINT(121)=2
      ENDIF
 
C...Set up kinematics of process.
      CALL PYINKI(0)
 
C...Loop over gamma-p or gamma-gamma alternatives.
      DO 160 IGA=1,MINT(121)
      MINT(122)=IGA
 
C...Select partonic subprocesses to be included in the simulation.
      CALL PYINPR
 
C...Count number of subprocesses on.
      MINT(48)=0
      DO 130 ISUB=1,200
      IF(MINT(50).EQ.0.AND.ISUB.GE.91.AND.ISUB.LE.96.AND.
     &MSUB(ISUB).EQ.1) THEN
        WRITE(MSTU(11),5200) ISUB,CHLH(MINT(41)),CHLH(MINT(42))
        STOP
      ELSEIF(MSUB(ISUB).EQ.1.AND.ISET(ISUB).EQ.-1) THEN
        WRITE(MSTU(11),5300) ISUB
        STOP
      ELSEIF(MSUB(ISUB).EQ.1.AND.ISET(ISUB).LE.-2) THEN
        WRITE(MSTU(11),5400) ISUB
        STOP
      ELSEIF(MSUB(ISUB).EQ.1) THEN
        MINT(48)=MINT(48)+1
      ENDIF
  130 CONTINUE
      IF(MINT(48).EQ.0) THEN
        WRITE(MSTU(11),5500)
        STOP
      ENDIF
      MINT(49)=MINT(48)-MSUB(91)-MSUB(92)-MSUB(93)-MSUB(94)
 
C...Reset variables for cross-section calculation.
      DO 150 I=0,200
      DO 140 J=1,3
      NGEN(I,J)=0
      XSEC(I,J)=0.
  140 CONTINUE
      DXSEC(I)=0D0
  150 CONTINUE
 
C...Find parametrized total cross-sections.
      CALL PYXTOT
 
C...Maxima of differential cross-sections.
      IF(MSTP(121).LE.1) CALL PYMAXI
 
C...Initialize possibility of pileup events.
      IF(MINT(121).GT.1) MSTP(131)=0
      IF(MSTP(131).NE.0) CALL PYPILE(1)
 
C...Initialize multiple interactions with variable impact parameter.
      IF(MINT(50).EQ.1.AND.(MINT(49).NE.0.OR.MSTP(131).NE.0).AND.
     &MSTP(82).GE.2) CALL PYMULT(1)
 
C...Save results for gamma-p and gamma-gamma alternatives.
      IF(MINT(121).GT.1) CALL PYSAVE(1,IGA)
  160 CONTINUE
 
C...Initialization finished.
  170 IF(MSTP(122).GE.1) WRITE(MSTU(11),5600)
 
C...Formats for initialization information.
 5100 FORMAT('1',18('*'),1X,'PYINIT: initialization of PYTHIA ',
     &'routines',1X,17('*'))
 5200 FORMAT(1X,'Error: process number ',I3,' not meaningful for ',A6,
     &'-',A6,' interactions.'/1X,'Execution stopped!')
 5300 FORMAT(1X,'Error: requested subprocess',I4,' not implemented.'/
     &1X,'Execution stopped!')
 5400 FORMAT(1X,'Error: requested subprocess',I4,' not existing.'/
     &1X,'Execution stopped!')
 5500 FORMAT(1X,'Error: no subprocess switched on.'/
     &1X,'Execution stopped.')
 5600 FORMAT(/1X,22('*'),1X,'PYINIT: initialization completed',1X,
     &22('*'))
 
      RETURN
      END
