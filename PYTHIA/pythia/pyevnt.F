 
C*********************************************************************
 
      SUBROUTINE PYEVNT
 
C...Administers the generation of a high-pT event via calls to
C...a number of subroutines.
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3)
      COMMON/PYUPPR/NUP,KUP(20,7),PUP(20,5),NFUP,IFUP(10,2),Q2UP(0:10)
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/
      SAVE /PYPARS/,/PYINT1/,/PYINT2/,/PYINT5/,/PYUPPR/
      DIMENSION VTX(4)
 
C...Initial values for some counters.
      N=0
      MINT(5)=MINT(5)+1
      MINT(7)=0
      MINT(8)=0
      MINT(83)=0
      MINT(84)=MSTP(126)
      MSTU(24)=0
      MSTU70=0
      MSTJ14=MSTJ(14)
 
C...If variable energies: redo incoming kinematics and cross-section.
      MSTI(61)=0
      IF(MSTP(171).EQ.1) THEN
        CALL PYINKI(1)
        IF(MSTI(61).EQ.1) THEN
          MINT(5)=MINT(5)-1
          RETURN
        ENDIF
        IF(MINT(121).GT.1) CALL PYSAVE(3,1)
        CALL PYXTOT
      ENDIF
 
C...Loop over number of pileup events; check space left.
      IF(MSTP(131).LE.0) THEN
        NPILE=1
      ELSE
        CALL PYPILE(2)
        NPILE=MINT(81)
      ENDIF
      DO 250 IPILE=1,NPILE
      IF(MINT(84)+100.GE.MSTU(4)) THEN
        CALL LUERRM(11,
     &  '(PYEVNT:) no more space in LUJETS for pileup events')
        IF(MSTU(21).GE.1) GOTO 260
      ENDIF
      MINT(82)=IPILE
 
C...Generate variables of hard scattering.
      MINT(51)=0
      MSTI(52)=0
  100 CONTINUE
      IF(MINT(51).NE.0.OR.MSTU(24).NE.0) MSTI(52)=MSTI(52)+1
      MINT(31)=0
      MINT(51)=0
      MINT(57)=0
      CALL PYRAND
      IF(MSTI(61).EQ.1) THEN
        MINT(5)=MINT(5)-1
        RETURN
      ENDIF
      IF(MINT(51).EQ.2) RETURN
      ISUB=MINT(1)
      IF(MSTP(111).EQ.-1) GOTO 240
 
      IF(ISUB.LE.90.OR.ISUB.GE.95) THEN
C...Hard scattering (including low-pT):
C...reconstruct kinematics and colour flow of hard scattering.
  110   MINT(51)=0
        CALL PYSCAT
        IF(MINT(51).EQ.1) GOTO 100
        IPU1=MINT(84)+1
        IPU2=MINT(84)+2
        IF(ISUB.EQ.95) GOTO 130
 
C...Showering of initial state partons (optional).
        ALAMSV=PARJ(81)
        PARJ(81)=PARP(72)
        IF(MSTP(61).GE.1.AND.MINT(47).GE.2) CALL PYSSPA(IPU1,IPU2)
        PARJ(81)=ALAMSV
        IF(MINT(51).EQ.1) GOTO 100
 
C...Showering of final state partons (optional).
        ALAMSV=PARJ(81)
        PARJ(81)=PARP(72)
        IF(MSTP(71).GE.1.AND.ISET(ISUB).GE.2.AND.ISET(ISUB).LE.10) THEN
          IPU3=MINT(84)+3
          IPU4=MINT(84)+4
          IF(ISET(ISUB).EQ.5.OR.ISET(ISUB).EQ.6) IPU4=-3
          QMAX=VINT(55)
          IF(ISET(ISUB).EQ.2) QMAX=SQRT(PARP(71))*VINT(55)
          CALL LUSHOW(IPU3,IPU4,QMAX)
        ELSEIF(MSTP(71).GE.1.AND.ISET(ISUB).EQ.11.AND.NFUP.GE.1) THEN
          DO 120 IUP=1,NFUP
          IPU3=IFUP(IUP,1)+MINT(84)
          IPU4=IFUP(IUP,2)+MINT(84)
          QMAX=SQRT(MAX(0.,Q2UP(IUP)))
          CALL LUSHOW(IPU3,IPU4,QMAX)
  120     CONTINUE
        ENDIF
        PARJ(81)=ALAMSV
 
C...Decay of final state resonances.
        IF(MSTP(41).GE.1.AND.ISET(ISUB).LE.10) CALL PYRESD
        IF(MINT(51).EQ.1) GOTO 100
        MINT(52)=N
 
C...Multiple interactions.
        IF(MSTP(81).GE.1.AND.MINT(50).EQ.1) CALL PYMULT(6)
        MINT(53)=N
 
C...Hadron remnants and primordial kT.
  130   CALL PYREMN(IPU1,IPU2)
        IF(MINT(51).EQ.1.AND.MINT(57).GE.1.AND.MINT(57).LE.5) GOTO 110
        IF(MINT(51).EQ.1) GOTO 100
 
      ELSE
C...Diffractive and elastic scattering.
        CALL PYDIFF
      ENDIF
 
C...Recalculate energies from momenta and masses (if desired).
      IF(MSTP(113).GE.1) THEN
        DO 140 I=MINT(83)+1,N
        IF(K(I,1).GT.0.AND.K(I,1).LE.10) P(I,4)=SQRT(P(I,1)**2+
     &  P(I,2)**2+P(I,3)**2+P(I,5)**2)
  140   CONTINUE
        NRECAL=N
      ENDIF
 
C...Rearrange partons along strings, check invariant mass cuts.
      MSTU(28)=0
      IF(MSTP(111).LE.0) MSTJ(14)=-1
      CALL LUPREP(MINT(84)+1)
      MSTJ(14)=MSTJ14
      IF(MSTP(112).EQ.1.AND.MSTU(28).EQ.3) GOTO 100
      IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) THEN
        DO 170 I=MINT(84)+1,N
        IF(K(I,2).EQ.94) THEN
          DO 160 I1=I+1,MIN(N,I+3)
          IF(K(I1,3).EQ.I) THEN
            K(I1,3)=MOD(K(I1,4)/MSTU(5),MSTU(5))
            IF(K(I1,3).EQ.0) THEN
              DO 150 II=MINT(84)+1,I-1
              IF(K(II,2).EQ.K(I1,2)) THEN
                IF(MOD(K(II,4),MSTU(5)).EQ.I1.OR.
     &          MOD(K(II,5),MSTU(5)).EQ.I1) K(I1,3)=II
              ENDIF
  150         CONTINUE
              IF(K(I+1,3).EQ.0) K(I+1,3)=K(I,3)
            ENDIF
          ENDIF
  160     CONTINUE
        ENDIF
  170   CONTINUE
        CALL LUEDIT(12)
        CALL LUEDIT(14)
        IF(MSTP(125).EQ.0) CALL LUEDIT(15)
        IF(MSTP(125).EQ.0) MINT(4)=0
        DO 190 I=MINT(83)+1,N
        IF(K(I,1).EQ.11.AND.K(I,4).EQ.0.AND.K(I,5).EQ.0) THEN
          DO 180 I1=I+1,N
          IF(K(I1,3).EQ.I.AND.K(I,4).EQ.0) K(I,4)=I1
          IF(K(I1,3).EQ.I) K(I,5)=I1
  180     CONTINUE
        ENDIF
  190   CONTINUE
      ENDIF
 
C...Introduce separators between sections in LULIST event listing.
      IF(IPILE.EQ.1.AND.MSTP(125).LE.0) THEN
        MSTU70=1
        MSTU(71)=N
      ELSEIF(IPILE.EQ.1) THEN
        MSTU70=3
        MSTU(71)=2
        MSTU(72)=MINT(4)
        MSTU(73)=N
      ENDIF
 
C...Go back to lab frame (needed for vertices, also in fragmentation).
      CALL PYFRAM(1)
 
C...Set nonvanishing production vertex (optional).
      IF(MSTP(151).EQ.1) THEN
        DO 200 J=1,4
        VTX(J)=PARP(150+J)*SQRT(-2.*LOG(MAX(1E-10,RLU(0))))*
     &  SIN(PARU(2)*RLU(0))
  200   CONTINUE
        DO 220 I=MINT(83)+1,N
        DO 210 J=1,4
        V(I,J)=V(I,J)+VTX(J)
  210   CONTINUE
  220   CONTINUE
      ENDIF
 
C...Perform hadronization (if desired).
      IF(MSTP(111).GE.1) THEN
        CALL LUEXEC
        IF(MSTU(24).NE.0) GOTO 100
      ENDIF
      IF(MSTP(113).GE.1) THEN
        DO 230 I=NRECAL,N
        IF(P(I,5).GT.0.) P(I,4)=SQRT(P(I,1)**2+
     &  P(I,2)**2+P(I,3)**2+P(I,5)**2)
  230   CONTINUE
      ENDIF
      IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) CALL LUEDIT(14)
 
C...Store event information and calculate Monte Carlo estimates of
C...subprocess cross-sections.
  240 IF(IPILE.EQ.1) CALL PYDOCU
 
C...Set counters for current pileup event and loop to next one.
      MSTI(41)=IPILE
      IF(IPILE.GE.2.AND.IPILE.LE.10) MSTI(40+IPILE)=ISUB
      IF(MSTU70.LT.10) THEN
        MSTU70=MSTU70+1
        MSTU(70+MSTU70)=N
      ENDIF
      MINT(83)=N
      MINT(84)=N+MSTP(126)
      IF(IPILE.LT.NPILE) CALL PYFRAM(2)
  250 CONTINUE
 
C...Generic information on pileup events. Reconstruct missing history.
      IF(MSTP(131).EQ.1.AND.MSTP(133).GE.1) THEN
        PARI(91)=VINT(132)
        PARI(92)=VINT(133)
        PARI(93)=VINT(134)
        IF(MSTP(133).GE.2) PARI(93)=PARI(93)*XSEC(0,3)/VINT(131)
      ENDIF
      CALL LUEDIT(16)
 
C...Transform to the desired coordinate frame.
  260 CALL PYFRAM(MSTP(124))
      MSTU(70)=MSTU70
      PARU(21)=VINT(1)
 
      RETURN
      END
