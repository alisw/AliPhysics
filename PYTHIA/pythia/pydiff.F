 
C*********************************************************************
 
      SUBROUTINE PYDIFF
 
C...Handles diffractive and elastic scattering.
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /LUJETS/,/LUDAT1/
      SAVE /PYPARS/,/PYINT1/
      DOUBLE PRECISION DBETAZ
 
C...Reset K, P and V vectors. Store incoming particles.
      DO 110 JT=1,MSTP(126)+10
      I=MINT(83)+JT
      DO 100 J=1,5
      K(I,J)=0
      P(I,J)=0.
      V(I,J)=0.
  100 CONTINUE
  110 CONTINUE
      N=MINT(84)
      MINT(3)=0
      MINT(21)=0
      MINT(22)=0
      MINT(23)=0
      MINT(24)=0
      MINT(4)=4
      DO 130 JT=1,2
      I=MINT(83)+JT
      K(I,1)=21
      K(I,2)=MINT(10+JT)
      DO 120 J=1,5
      P(I,J)=VINT(285+5*JT+J)
  120 CONTINUE
  130 CONTINUE
      MINT(6)=2
 
C...Subprocess; kinematics.
      SQLAM=(VINT(2)-VINT(63)-VINT(64))**2-4.*VINT(63)*VINT(64)
      PZ=SQRT(SQLAM)/(2.*VINT(1))
      DO 200 JT=1,2
      I=MINT(83)+JT
      PE=(VINT(2)+VINT(62+JT)-VINT(65-JT))/(2.*VINT(1))
      KFH=MINT(102+JT)
 
C...Elastically scattered particle.
      IF(MINT(16+JT).LE.0) THEN
        N=N+1
        K(N,1)=1
        K(N,2)=KFH
        K(N,3)=I+2
        P(N,3)=PZ*(-1)**(JT+1)
        P(N,4)=PE
        P(N,5)=SQRT(VINT(62+JT))
 
C...Decay rho from elastic scattering of gamma with sin**2(theta)
C...distribution of decay products (in rho rest frame).
        IF(KFH.EQ.113.AND.MINT(10+JT).EQ.22.AND.MSTP(102).EQ.1) THEN
          NSAV=N
          DBETAZ=DBLE(P(N,3))/SQRT(DBLE(P(N,3))**2+DBLE(P(N,5))**2)
          P(N,3)=0.
          P(N,4)=P(N,5)
          CALL LUDECY(NSAV)
          IF(N.EQ.NSAV+2.AND.IABS(K(NSAV+1,2)).EQ.211) THEN
            PHI=ULANGL(P(NSAV+1,1),P(NSAV+1,2))
            CALL LUDBRB(NSAV+1,NSAV+2,0.,-PHI,0D0,0D0,0D0)
            THE=ULANGL(P(NSAV+1,3),P(NSAV+1,1))
            CALL LUDBRB(NSAV+1,NSAV+2,-THE,0.,0D0,0D0,0D0)
  140       CTHE=2.*RLU(0)-1.
            IF(1.-CTHE**2.LT.RLU(0)) GOTO 140
            CALL LUDBRB(NSAV+1,NSAV+2,ACOS(CTHE),PHI,0D0,0D0,0D0)
          ENDIF
          CALL LUDBRB(NSAV,NSAV+2,0.,0.,0D0,0D0,DBETAZ)
        ENDIF
 
C...Diffracted particle: low-mass system to two particles.
      ELSEIF(VINT(62+JT).LT.(VINT(66+JT)+PARP(103))**2) THEN
        N=N+2
        K(N-1,1)=1
        K(N,1)=1
        K(N-1,3)=I+2
        K(N,3)=I+2
        PMMAS=SQRT(VINT(62+JT))
        NTRY=0
  150   NTRY=NTRY+1
        IF(NTRY.LT.20) THEN
          MINT(105)=MINT(102+JT)
          MINT(109)=MINT(106+JT)
          CALL PYSPLI(KFH,21,KFL1,KFL2)
          CALL LUKFDI(KFL1,0,KFL3,KF1)
          IF(KF1.EQ.0) GOTO 150
          CALL LUKFDI(KFL2,-KFL3,KFLDUM,KF2)
          IF(KF2.EQ.0) GOTO 150
        ELSE
          KF1=KFH
          KF2=111
        ENDIF
        PM1=ULMASS(KF1)
        PM2=ULMASS(KF2)
        IF(PM1+PM2+PARJ(64).GT.PMMAS) GOTO 150
        K(N-1,2)=KF1
        K(N,2)=KF2
        P(N-1,5)=PM1
        P(N,5)=PM2
        PZP=SQRT(MAX(0.,(PMMAS**2-PM1**2-PM2**2)**2-4.*PM1**2*PM2**2))/
     &  (2.*PMMAS)
        P(N-1,3)=PZP
        P(N,3)=-PZP
        P(N-1,4)=SQRT(PM1**2+PZP**2)
        P(N,4)=SQRT(PM2**2+PZP**2)
        CALL LUDBRB(N-1,N,ACOS(2.*RLU(0)-1.),PARU(2)*RLU(0),0D0,0D0,0D0)
        DBETAZ=DBLE(PZ)*(-1)**(JT+1)/SQRT(DBLE(PZ)**2+DBLE(PMMAS)**2)
        CALL LUDBRB(N-1,N,0.,0.,0D0,0D0,DBETAZ)
 
C...Diffracted particle: valence quark kicked out.
      ELSEIF(MSTP(101).EQ.1.OR.(MSTP(101).EQ.3.AND.RLU(0).LT.PARP(101)))
     &THEN
        N=N+2
        K(N-1,1)=2
        K(N,1)=1
        K(N-1,3)=I+2
        K(N,3)=I+2
        MINT(105)=MINT(102+JT)
        MINT(109)=MINT(106+JT)
        CALL PYSPLI(KFH,21,K(N,2),K(N-1,2))
        P(N-1,5)=ULMASS(K(N-1,2))
        P(N,5)=ULMASS(K(N,2))
        SQLAM=(VINT(62+JT)-P(N-1,5)**2-P(N,5)**2)**2-
     &  4.*P(N-1,5)**2*P(N,5)**2
        P(N-1,3)=(PE*SQRT(SQLAM)+PZ*(VINT(62+JT)+P(N-1,5)**2-
     &  P(N,5)**2))/(2.*VINT(62+JT))*(-1)**(JT+1)
        P(N-1,4)=SQRT(P(N-1,3)**2+P(N-1,5)**2)
        P(N,3)=PZ*(-1)**(JT+1)-P(N-1,3)
        P(N,4)=SQRT(P(N,3)**2+P(N,5)**2)
 
C...Diffracted particle: gluon kicked out.
      ELSE
        N=N+3
        K(N-2,1)=2
        K(N-1,1)=2
        K(N,1)=1
        K(N-2,3)=I+2
        K(N-1,3)=I+2
        K(N,3)=I+2
        MINT(105)=MINT(102+JT)
        MINT(109)=MINT(106+JT)
        CALL PYSPLI(KFH,21,K(N,2),K(N-2,2))
        K(N-1,2)=21
        P(N-2,5)=ULMASS(K(N-2,2))
        P(N-1,5)=0.
        P(N,5)=ULMASS(K(N,2))
C...Energy distribution for particle into two jets.
  160   IMB=1
        IF(MOD(KFH/1000,10).NE.0) IMB=2
        CHIK=PARP(92+2*IMB)
        IF(MSTP(92).LE.1) THEN
          IF(IMB.EQ.1) CHI=RLU(0)
          IF(IMB.EQ.2) CHI=1.-SQRT(RLU(0))
        ELSEIF(MSTP(92).EQ.2) THEN
          CHI=1.-RLU(0)**(1./(1.+CHIK))
        ELSEIF(MSTP(92).EQ.3) THEN
          CUT=2.*0.3/VINT(1)
  170     CHI=RLU(0)**2
          IF((CHI**2/(CHI**2+CUT**2))**0.25*(1.-CHI)**CHIK.LT.
     &    RLU(0)) GOTO 170
        ELSEIF(MSTP(92).EQ.4) THEN
          CUT=2.*0.3/VINT(1)
          CUTR=(1.+SQRT(1.+CUT**2))/CUT
  180     CHIR=CUT*CUTR**RLU(0)
          CHI=(CHIR**2-CUT**2)/(2.*CHIR)
          IF((1.-CHI)**CHIK.LT.RLU(0)) GOTO 180
        ELSE
          CUT=2.*0.3/VINT(1)
          CUTA=CUT**(1.-PARP(98))
          CUTB=(1.+CUT)**(1.-PARP(98))
  190     CHI=(CUTA+RLU(0)*(CUTB-CUTA))**(1./(1.-PARP(98)))
          IF(((CHI+CUT)**2/(2.*(CHI**2+CUT**2)))**
     &    (0.5*PARP(98))*(1.-CHI)**CHIK.LT.RLU(0)) GOTO 190
        ENDIF
        IF(CHI.LT.P(N,5)**2/VINT(62+JT).OR.CHI.GT.1.-P(N-2,5)**2/
     &  VINT(62+JT)) GOTO 160
        SQM=P(N-2,5)**2/(1.-CHI)+P(N,5)**2/CHI
        IF((SQRT(SQM)+PARJ(32))**2.GE.VINT(62+JT)) GOTO 160
        PZI=(PE*(VINT(62+JT)-SQM)+PZ*(VINT(62+JT)+SQM))/
     &  (2.*VINT(62+JT))
        PEI=SQRT(PZI**2+SQM)
        PQQP=(1.-CHI)*(PEI+PZI)
        P(N-2,3)=0.5*(PQQP-P(N-2,5)**2/PQQP)*(-1)**(JT+1)
        P(N-2,4)=SQRT(P(N-2,3)**2+P(N-2,5)**2)
        P(N-1,4)=0.5*(VINT(62+JT)-SQM)/(PEI+PZI)
        P(N-1,3)=P(N-1,4)*(-1)**JT
        P(N,3)=PZI*(-1)**(JT+1)-P(N-2,3)
        P(N,4)=SQRT(P(N,3)**2+P(N,5)**2)
      ENDIF
 
C...Documentation lines.
      K(I+2,1)=21
      IF(MINT(16+JT).EQ.0) K(I+2,2)=KFH
      IF(MINT(16+JT).NE.0) K(I+2,2)=10*(KFH/10)
      K(I+2,3)=I
      P(I+2,3)=PZ*(-1)**(JT+1)
      P(I+2,4)=PE
      P(I+2,5)=SQRT(VINT(62+JT))
  200 CONTINUE
 
C...Rotate outgoing partons/particles using cos(theta).
      IF(VINT(23).LT.0.9) THEN
        CALL LUDBRB(MINT(83)+3,N,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0)
      ELSE
        CALL LUDBRB(MINT(83)+3,N,ASIN(VINT(59)),VINT(24),0D0,0D0,0D0)
      ENDIF
 
      RETURN
      END
