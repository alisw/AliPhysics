 
C*********************************************************************
 
      SUBROUTINE PYINKI(MODKI)
 
C...Sets up kinematics, including rotations and boosts to/from CM frame.
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/
 
C...Set initial flavour state.
      N=2
      DO 100 I=1,2
      K(I,1)=1
      K(I,2)=MINT(10+I)
  100 CONTINUE
 
C...Reset boost. Do kinematics for various cases.
      DO 110 J=6,10
      VINT(J)=0.
  110 CONTINUE
 
C...Set up kinematics for events defined in CM frame.
      IF(MINT(111).EQ.1) THEN
        WIN=VINT(290)
        IF(MODKI.EQ.1) WIN=PARP(171)*VINT(290)
        S=WIN**2
        P(1,5)=VINT(3)
        P(2,5)=VINT(4)
        P(1,1)=0.
        P(1,2)=0.
        P(2,1)=0.
        P(2,2)=0.
        P(1,3)=SQRT(((S-P(1,5)**2-P(2,5)**2)**2-(2.*P(1,5)*P(2,5))**2)/
     &  (4.*S))
        P(2,3)=-P(1,3)
        P(1,4)=SQRT(P(1,3)**2+P(1,5)**2)
        P(2,4)=SQRT(P(2,3)**2+P(2,5)**2)
 
C...Set up kinematics for fixed target events.
      ELSEIF(MINT(111).EQ.2) THEN
        WIN=VINT(290)
        IF(MODKI.EQ.1) WIN=PARP(171)*VINT(290)
        P(1,5)=VINT(3)
        P(2,5)=VINT(4)
        P(1,1)=0.
        P(1,2)=0.
        P(2,1)=0.
        P(2,2)=0.
        P(1,3)=WIN
        P(1,4)=SQRT(P(1,3)**2+P(1,5)**2)
        P(2,3)=0.
        P(2,4)=P(2,5)
        S=P(1,5)**2+P(2,5)**2+2.*P(2,4)*P(1,4)
        VINT(10)=P(1,3)/(P(1,4)+P(2,4))
        CALL LUROBO(0.,0.,0.,0.,-VINT(10))
 
C...Set up kinematics for events in user-defined frame.
      ELSEIF(MINT(111).EQ.3) THEN
        P(1,5)=VINT(3)
        P(2,5)=VINT(4)
        P(1,4)=SQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)
        P(2,4)=SQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)
        DO 120 J=1,3
        VINT(7+J)=(DBLE(P(1,J))+DBLE(P(2,J)))/DBLE(P(1,4)+P(2,4))
  120   CONTINUE
        CALL LUROBO(0.,0.,-VINT(8),-VINT(9),-VINT(10))
        VINT(7)=ULANGL(P(1,1),P(1,2))
        CALL LUROBO(0.,-VINT(7),0.,0.,0.)
        VINT(6)=ULANGL(P(1,3),P(1,1))
        CALL LUROBO(-VINT(6),0.,0.,0.,0.)
        S=P(1,5)**2+P(2,5)**2+2.*(P(1,4)*P(2,4)-P(1,3)*P(2,3))
 
C...Set up kinematics for events with user-defined four-vectors.
      ELSEIF(MINT(111).EQ.4) THEN
        PMS1=P(1,4)**2-P(1,1)**2-P(1,2)**2-P(1,3)**2
        P(1,5)=SIGN(SQRT(ABS(PMS1)),PMS1)
        PMS2=P(2,4)**2-P(2,1)**2-P(2,2)**2-P(2,3)**2
        P(2,5)=SIGN(SQRT(ABS(PMS2)),PMS2)
        DO 130 J=1,3
        VINT(7+J)=(DBLE(P(1,J))+DBLE(P(2,J)))/DBLE(P(1,4)+P(2,4))
  130   CONTINUE
        CALL LUROBO(0.,0.,-VINT(8),-VINT(9),-VINT(10))
        VINT(7)=ULANGL(P(1,1),P(1,2))
        CALL LUROBO(0.,-VINT(7),0.,0.,0.)
        VINT(6)=ULANGL(P(1,3),P(1,1))
        CALL LUROBO(-VINT(6),0.,0.,0.,0.)
        S=(P(1,4)+P(2,4))**2
 
C...Set up kinematics for events with user-defined five-vectors.
      ELSEIF(MINT(111).EQ.5) THEN
        DO 140 J=1,3
        VINT(7+J)=(DBLE(P(1,J))+DBLE(P(2,J)))/DBLE(P(1,4)+P(2,4))
  140   CONTINUE
        CALL LUROBO(0.,0.,-VINT(8),-VINT(9),-VINT(10))
        VINT(7)=ULANGL(P(1,1),P(1,2))
        CALL LUROBO(0.,-VINT(7),0.,0.,0.)
        VINT(6)=ULANGL(P(1,3),P(1,1))
        CALL LUROBO(-VINT(6),0.,0.,0.,0.)
        S=(P(1,4)+P(2,4))**2
      ENDIF
 
C...Return or error for too low CM energy.
      IF(MODKI.EQ.1.AND.S.LT.PARP(2)**2) THEN
        IF(MSTP(172).LE.1) THEN
          CALL LUERRM(23,
     &    '(PYINKI:) too low invariant mass in this event')
        ELSE
          MSTI(61)=1
          RETURN
        ENDIF
      ENDIF
 
C...Save information on incoming particles.
      VINT(1)=SQRT(S)
      VINT(2)=S
      IF(MINT(111).GE.4) VINT(3)=P(1,5)
      IF(MINT(111).GE.4) VINT(4)=P(2,5)
      VINT(5)=P(1,3)
      IF(MODKI.EQ.0) VINT(289)=S
      DO 150 J=1,5
      V(1,J)=0.
      V(2,J)=0.
      VINT(290+J)=P(1,J)
      VINT(295+J)=P(2,J)
  150 CONTINUE
 
C...Store pT cut-off and related constants to be used in generation.
      IF(MODKI.EQ.0) VINT(285)=CKIN(3)
      IF(MSTP(82).LE.1) THEN
        IF(MINT(121).GT.1) PARP(81)=1.30+0.15*LOG(VINT(1)/200.)/
     &  LOG(900./200.)
        PTMN=PARP(81)
      ELSE
        IF(MINT(121).GT.1) PARP(82)=1.25+0.15*LOG(VINT(1)/200.)/
     &  LOG(900./200.)
        PTMN=PARP(82)
      ENDIF
      VINT(149)=4.*PTMN**2/S
 
      RETURN
      END
