 
C*********************************************************************
 
      SUBROUTINE PYINPR
 
C...Selects partonic subprocesses to be included in the simulation.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      SAVE /LUDAT1/,/LUDAT3/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/
 
C...Reset processes to be included.
      IF(MSEL.NE.0) THEN
        DO 100 I=1,200
        MSUB(I)=0
  100   CONTINUE
      ENDIF
 
C...For e-gamma witn MSTP(14)=10 allow mixture of VMD and anomalous.
      IF(MINT(121).EQ.2) THEN
        MSUB(10)=1
        MINT(123)=MINT(122)+1
 
C...For gamma-p or gamma-gamma with MSTP(14)=10 allow mixture.
C...Here also set a few parameters otherwise normally not touched.
      ELSEIF(MINT(121).GT.1) THEN
 
C...Structure functions dampened at small Q2; go to low energies,
C...alpha_s <1; no minimum pT cut-off a priori.
        MSTP(57)=3
        MSTP(85)=0
        PARP(2)=2.
        PARU(115)=1.
        CKIN(5)=0.2
        CKIN(6)=0.2
 
C...Define pT cut-off parameters and whether run involves low-pT.
        IF(MSTP(82).LE.1) THEN
          PTMVMD=1.30+0.15*LOG(VINT(1)/200.)/LOG(900./200.)
        ELSE
          PTMVMD=1.25+0.15*LOG(VINT(1)/200.)/LOG(900./200.)
        ENDIF
        PTMDIR=PARP(15)
        PTMANO=PTMVMD
        IF(MSTP(15).EQ.5) PTMANO=0.70+0.17*LOG(1.+0.05*VINT(1))**2
        IPTL=1
        IF(VINT(285).GT.MAX(PTMVMD,PTMDIR,PTMANO)) IPTL=0
        IF(MSEL.EQ.2) IPTL=1
 
C...Set up for p/VMD * VMD.
        IF(MINT(122).EQ.1) THEN
          MINT(123)=2
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(IPTL.EQ.1) MSUB(95)=1
          IF(MSEL.EQ.2) THEN
            MSUB(91)=1
            MSUB(92)=1
            MSUB(93)=1
            MSUB(94)=1
          ENDIF
          PARP(81)=PTMVMD
          PARP(82)=PTMVMD
          IF(IPTL.EQ.1) CKIN(3)=0.
 
C...Set up for p/VMD * direct gamma.
        ELSEIF(MINT(122).EQ.2) THEN
          MINT(123)=0
          IF(MINT(121).EQ.6) MINT(123)=5
          MSUB(33)=1
          MSUB(54)=1
          IF(IPTL.EQ.1) CKIN(3)=PTMDIR
 
C...Set up for p/VMD * anomalous gamma.
        ELSEIF(MINT(122).EQ.3) THEN
          MINT(123)=3
          IF(MINT(121).EQ.6) MINT(123)=7
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(MSTP(82).GE.2) MSTP(85)=1
          IF(IPTL.EQ.1) CKIN(3)=PTMANO
 
C...Set up for direct * direct gamma (switch off leptons).
        ELSEIF(MINT(122).EQ.4) THEN
          MINT(123)=0
          MSUB(58)=1
          DO 110 II=MDCY(22,2),MDCY(22,2)+MDCY(22,3)-1
          IF(IABS(KFDP(II,1)).GE.10) MDME(II,1)=MIN(0,MDME(II,1))
  110     CONTINUE
          IF(IPTL.EQ.1) CKIN(3)=PTMDIR
 
C...Set up for direct * anomalous gamma.
        ELSEIF(MINT(122).EQ.5) THEN
          MINT(123)=6
          MSUB(33)=1
          MSUB(54)=1
          IF(IPTL.EQ.1) CKIN(3)=PTMANO
 
C...Set up for anomalous * anomalous gamma.
        ELSEIF(MINT(122).EQ.6) THEN
          MINT(123)=3
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(MSTP(82).GE.2) MSTP(85)=1
          IF(IPTL.EQ.1) CKIN(3)=PTMANO
        ENDIF
 
C...End of special set up for gamma-p and gamma-gamma.
        CKIN(1)=2.*CKIN(3)
      ENDIF
 
C...Flavour information for individual beams.
      DO 120 I=1,2
      MINT(40+I)=1
      IF(MINT(123).GE.1.AND.MINT(10+I).EQ.22) MINT(40+I)=2
      IF(IABS(MINT(10+I)).GT.100) MINT(40+I)=2
      IF(MINT(10+I).EQ.28.OR.MINT(10+I).EQ.29) MINT(40+I)=2
      MINT(44+I)=MINT(40+I)
      IF(MSTP(11).GE.1.AND.IABS(MINT(10+I)).EQ.11) MINT(44+I)=3
  120 CONTINUE
 
C...If two gammas, whereof one direct, pick the first.
      IF(MINT(11).EQ.22.AND.MINT(12).EQ.22) THEN
        IF(MINT(123).GE.4.AND.MINT(123).LE.6) THEN
          MINT(41)=1
          MINT(45)=1
        ENDIF
      ELSEIF(MINT(11).EQ.22.OR.MINT(12).EQ.22) THEN
        IF(MINT(123).GE.4) CALL LUERRM(26,
     &  '(PYINPR:) unallowed MSTP(14) code for single photon')
      ENDIF
 
C...Flavour information on combination of incoming particles.
      MINT(43)=2*MINT(41)+MINT(42)-2
      MINT(44)=MINT(43)
      IF(MINT(123).LE.0) THEN
        IF(MINT(11).EQ.22) MINT(43)=MINT(43)+2
        IF(MINT(12).EQ.22) MINT(43)=MINT(43)+1
      ELSEIF(MINT(123).LE.3) THEN
        IF(MINT(11).EQ.22) MINT(44)=MINT(44)-2
        IF(MINT(12).EQ.22) MINT(44)=MINT(44)-1
      ELSEIF(MINT(11).EQ.22.AND.MINT(12).EQ.22) THEN
        MINT(43)=4
        MINT(44)=1
      ENDIF
      MINT(47)=2*MIN(2,MINT(45))+MIN(2,MINT(46))-2
      IF(MIN(MINT(45),MINT(46)).EQ.3) MINT(47)=5
      MINT(50)=0
      IF(MINT(41).EQ.2.AND.MINT(42).EQ.2) MINT(50)=1
      IF((MINT(11).EQ.22.OR.MINT(12).EQ.22).AND.MINT(123).GE.3)
     &MINT(50)=0
      MINT(107)=0
      IF(MINT(11).EQ.22) THEN
        MINT(107)=MINT(123)
        IF(MINT(123).GE.4) MINT(107)=0
        IF(MINT(123).EQ.7) MINT(107)=2
      ENDIF
      MINT(108)=0
      IF(MINT(12).EQ.22) THEN
        MINT(108)=MINT(123)
        IF(MINT(123).GE.4) MINT(108)=MINT(123)-3
        IF(MINT(123).EQ.7) MINT(108)=3
      ENDIF
 
C...Select default processes according to incoming beams
C...(already done for gamma-p and gamma-gamma with MSTP(14)=10).
      IF(MINT(121).GT.1) THEN
      ELSEIF(MSEL.EQ.1.OR.MSEL.EQ.2) THEN
 
        IF(MINT(43).EQ.1) THEN
C...Lepton + lepton -> gamma/Z0 or W.
          IF(MINT(11)+MINT(12).EQ.0) MSUB(1)=1
          IF(MINT(11)+MINT(12).NE.0) MSUB(2)=1
 
        ELSEIF(MINT(43).LE.3.AND.MINT(123).EQ.0.AND.
     &  (MINT(11).EQ.22.OR.MINT(12).EQ.22)) THEN
C...Unresolved photon + lepton: Compton scattering.
          MSUB(34)=1
 
        ELSEIF(MINT(43).LE.3) THEN
C...Lepton + hadron: deep inelastic scattering.
          MSUB(10)=1
 
        ELSEIF(MINT(123).EQ.0.AND.MINT(11).EQ.22.AND.
     &  MINT(12).EQ.22) THEN
C...Two unresolved photons: fermion pair production.
          MSUB(58)=1
 
        ELSEIF((MINT(123).EQ.0.AND.(MINT(11).EQ.22.OR.MINT(12).EQ.22))
     &  .OR.(MINT(123).GE.4.AND.MINT(123).LE.6.AND.MINT(11).EQ.22.AND.
     &   MINT(12).EQ.22)) THEN
C...Unresolved photon + hadron: photon-parton scattering.
          MSUB(33)=1
          MSUB(34)=1
          MSUB(54)=1
 
        ELSEIF(MSEL.EQ.1) THEN
C...High-pT QCD processes:
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(MSTP(82).LE.1.AND.CKIN(3).LT.PARP(81)) MSUB(95)=1
          IF(MSTP(82).GE.2.AND.CKIN(3).LT.PARP(82)) MSUB(95)=1
          IF(MSUB(95).EQ.1.AND.MINT(50).EQ.0) MSUB(95)=0
 
        ELSE
C...All QCD processes:
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          MSUB(91)=1
          MSUB(92)=1
          MSUB(93)=1
          MSUB(94)=1
          MSUB(95)=1
        ENDIF
 
      ELSEIF(MSEL.GE.4.AND.MSEL.LE.8) THEN
C...Heavy quark production.
        MSUB(81)=1
        MSUB(82)=1
        MSUB(84)=1
        DO 130 J=1,MIN(8,MDCY(21,3))
        MDME(MDCY(21,2)+J-1,1)=0
  130   CONTINUE
        MDME(MDCY(21,2)+MSEL-1,1)=1
        MSUB(85)=1
        DO 140 J=1,MIN(12,MDCY(22,3))
        MDME(MDCY(22,2)+J-1,1)=0
  140   CONTINUE
        MDME(MDCY(22,2)+MSEL-1,1)=1
 
      ELSEIF(MSEL.EQ.10) THEN
C...Prompt photon production:
        MSUB(14)=1
        MSUB(18)=1
        MSUB(29)=1
 
      ELSEIF(MSEL.EQ.11) THEN
C...Z0/gamma* production:
        MSUB(1)=1
 
      ELSEIF(MSEL.EQ.12) THEN
C...W+/- production:
        MSUB(2)=1
 
      ELSEIF(MSEL.EQ.13) THEN
C...Z0 + jet:
        MSUB(15)=1
        MSUB(30)=1
 
      ELSEIF(MSEL.EQ.14) THEN
C...W+/- + jet:
        MSUB(16)=1
        MSUB(31)=1
 
      ELSEIF(MSEL.EQ.15) THEN
C...Z0 & W+/- pair production:
        MSUB(19)=1
        MSUB(20)=1
        MSUB(22)=1
        MSUB(23)=1
        MSUB(25)=1
 
      ELSEIF(MSEL.EQ.16) THEN
C...H0 production:
        MSUB(3)=1
        MSUB(102)=1
        MSUB(103)=1
        MSUB(123)=1
        MSUB(124)=1
 
      ELSEIF(MSEL.EQ.17) THEN
C...H0 & Z0 or W+/- pair production:
        MSUB(24)=1
        MSUB(26)=1
 
      ELSEIF(MSEL.EQ.18) THEN
C...H0 production; interesting processes in e+e-.
        MSUB(24)=1
        MSUB(103)=1
        MSUB(123)=1
        MSUB(124)=1
 
      ELSEIF(MSEL.EQ.19) THEN
C...H0, H'0 and A0 production; interesting processes in e+e-.
        MSUB(24)=1
        MSUB(103)=1
        MSUB(123)=1
        MSUB(124)=1
        MSUB(153)=1
        MSUB(171)=1
        MSUB(173)=1
        MSUB(174)=1
        MSUB(158)=1
        MSUB(176)=1
        MSUB(178)=1
        MSUB(179)=1
 
      ELSEIF(MSEL.EQ.21) THEN
C...Z'0 production:
        MSUB(141)=1
 
      ELSEIF(MSEL.EQ.22) THEN
C...W'+/- production:
        MSUB(142)=1
 
      ELSEIF(MSEL.EQ.23) THEN
C...H+/- production:
        MSUB(143)=1
 
      ELSEIF(MSEL.EQ.24) THEN
C...R production:
        MSUB(144)=1
 
      ELSEIF(MSEL.EQ.25) THEN
C...LQ (leptoquark) production.
        MSUB(145)=1
        MSUB(162)=1
        MSUB(163)=1
        MSUB(164)=1
 
      ELSEIF(MSEL.GE.35.AND.MSEL.LE.38) THEN
C...Production of one heavy quark (W exchange):
        MSUB(83)=1
        DO 150 J=1,MIN(8,MDCY(21,3))
        MDME(MDCY(21,2)+J-1,1)=0
  150   CONTINUE
        MDME(MDCY(21,2)+MSEL-31,1)=1
      ENDIF
 
C...Find heaviest new quark flavour allowed in processes 81-84.
      KFLQM=1
      DO 160 I=1,MIN(8,MDCY(21,3))
      IDC=I+MDCY(21,2)-1
      IF(MDME(IDC,1).LE.0) GOTO 160
      KFLQM=I
  160 CONTINUE
      IF(MSTP(7).GE.1.AND.MSTP(7).LE.8.AND.(MSEL.LE.3.OR.MSEL.GE.9))
     &KFLQM=MSTP(7)
      MINT(55)=KFLQM
      KFPR(81,1)=KFLQM
      KFPR(81,2)=KFLQM
      KFPR(82,1)=KFLQM
      KFPR(82,2)=KFLQM
      KFPR(83,1)=KFLQM
      KFPR(84,1)=KFLQM
      KFPR(84,2)=KFLQM
 
C...Find heaviest new fermion flavour allowed in process 85.
      KFLFM=1
      DO 170 I=1,MIN(12,MDCY(22,3))
      IDC=I+MDCY(22,2)-1
      IF(MDME(IDC,1).LE.0) GOTO 170
      KFLFM=KFDP(IDC,1)
  170 CONTINUE
      IF(((MSTP(7).GE.1.AND.MSTP(7).LE.8).OR.(MSTP(7).GE.11.AND.
     &MSTP(7).LE.18)).AND.(MSEL.LE.3.OR.MSEL.GE.9)) KFLFM=MSTP(7)
      MINT(56)=KFLFM
      KFPR(85,1)=KFLFM
      KFPR(85,2)=KFLFM
 
      RETURN
      END
