*CMZ :          17/07/98  15.44.31  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUPREP(IP)

C...Purpose: to rearrange partons along strings, to allow small systems
C...to collapse into one or two particles and to check flavours.
      IMPLICIT DOUBLE PRECISION(D)
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEEP,LUDAT3.
      COMMON /LUDAT3/ MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      SAVE /LUDAT3/
*KEND.
      DIMENSION DPS(5),DPC(5),UE(3)

C...Rearrange parton shower product listing along strings: begin loop.
      I1=N
      DO 130 MQGST=1,2
      DO 120 I=MAX(1,IP),N
      IF(K(I,1).NE.3) GOTO 120
      KC=LUCOMP(K(I,2))
      IF(KC.EQ.0) GOTO 120
      KQ=KCHG(KC,2)
      IF(KQ.EQ.0.OR.(MQGST.EQ.1.AND.KQ.EQ.2)) GOTO 120

C...Pick up loose string end.
      KCS=4
      IF(KQ*ISIGN(1,K(I,2)).LT.0) KCS=5
      IA=I
      NSTP=0
  100 NSTP=NSTP+1
      IF(NSTP.GT.4*N) THEN
        CALL LUERRM(14,'(LUPREP:) caught in infinite loop')
        RETURN
      ENDIF

C...Copy undecayed parton.
      IF(K(IA,1).EQ.3) THEN
        IF(I1.GE.MSTU(4)-MSTU(32)-5) THEN
          CALL LUERRM(11,'(LUPREP:) no more memory left in LUJETS')
          RETURN
        ENDIF
        I1=I1+1
        K(I1,1)=2
        IF(NSTP.GE.2.AND.IABS(K(IA,2)).NE.21) K(I1,1)=1
        K(I1,2)=K(IA,2)
        K(I1,3)=IA
        K(I1,4)=0
        K(I1,5)=0
        DO 110 J=1,5
        P(I1,J)=P(IA,J)
  110   V(I1,J)=V(IA,J)
        K(IA,1)=K(IA,1)+10
        IF(K(I1,1).EQ.1) GOTO 120
      ENDIF

C...Go to next parton in colour space.
      IB=IA
      IF(MOD(K(IB,KCS)/MSTU(5)**2,2).EQ.0.AND.MOD(K(IB,KCS),MSTU(5)).
     &NE.0) THEN
        IA=MOD(K(IB,KCS),MSTU(5))
        K(IB,KCS)=K(IB,KCS)+MSTU(5)**2
        MREV=0
      ELSE
        IF(K(IB,KCS).GE.2*MSTU(5)**2.OR.MOD(K(IB,KCS)/MSTU(5),MSTU(5)).
     &  EQ.0) KCS=9-KCS
        IA=MOD(K(IB,KCS)/MSTU(5),MSTU(5))
        K(IB,KCS)=K(IB,KCS)+2*MSTU(5)**2
        MREV=1
      ENDIF
      IF(IA.LE.0.OR.IA.GT.N) THEN
        CALL LUERRM(12,'(LUPREP:) colour rearrangement failed')
        RETURN
      ENDIF
      IF(MOD(K(IA,4)/MSTU(5),MSTU(5)).EQ.IB.OR.MOD(K(IA,5)/MSTU(5),
     &MSTU(5)).EQ.IB) THEN
        IF(MREV.EQ.1) KCS=9-KCS
        IF(MOD(K(IA,KCS)/MSTU(5),MSTU(5)).NE.IB) KCS=9-KCS
        K(IA,KCS)=K(IA,KCS)+2*MSTU(5)**2
      ELSE
        IF(MREV.EQ.0) KCS=9-KCS
        IF(MOD(K(IA,KCS),MSTU(5)).NE.IB) KCS=9-KCS
        K(IA,KCS)=K(IA,KCS)+MSTU(5)**2
      ENDIF
      IF(IA.NE.I) GOTO 100
      K(I1,1)=1
  120 CONTINUE
  130 CONTINUE
      N=I1

C...Find lowest-mass colour singlet jet system, OK if above threshold.
      IF(MSTJ(14).LE.0) GOTO 320
      NS=N
  140 NSIN=N-NS
      PDM=1.+PARJ(32)
      IC=0
      DO 190 I=MAX(1,IP),NS
      IF(K(I,1).NE.1.AND.K(I,1).NE.2) THEN
      ELSEIF(K(I,1).EQ.2.AND.IC.EQ.0) THEN
        NSIN=NSIN+1
        IC=I
        DO 150 J=1,4
  150   DPS(J)=P(I,J)
        MSTJ(93)=1
        DPS(5)=ULMASS(K(I,2))
      ELSEIF(K(I,1).EQ.2) THEN
        DO 160 J=1,4
  160   DPS(J)=DPS(J)+P(I,J)
      ELSEIF(IC.NE.0.AND.KCHG(LUCOMP(K(I,2)),2).NE.0) THEN
        DO 170 J=1,4
  170   DPS(J)=DPS(J)+P(I,J)
        MSTJ(93)=1
        DPS(5)=DPS(5)+ULMASS(K(I,2))
        PD=SQRT(MAX(0D0,DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2))-DPS(5)
        IF(PD.LT.PDM) THEN
          PDM=PD
          DO 180 J=1,5
  180     DPC(J)=DPS(J)
          IC1=IC
          IC2=I
        ENDIF
        IC=0
      ELSE
        NSIN=NSIN+1
      ENDIF
  190 CONTINUE
      IF(PDM.GE.PARJ(32)) GOTO 320

C...Fill small-mass system as cluster.
      NSAV=N
      PECM=SQRT(MAX(0D0,DPC(4)**2-DPC(1)**2-DPC(2)**2-DPC(3)**2))
      K(N+1,1)=11
      K(N+1,2)=91
      K(N+1,3)=IC1
      K(N+1,4)=N+2
      K(N+1,5)=N+3
      P(N+1,1)=DPC(1)
      P(N+1,2)=DPC(2)
      P(N+1,3)=DPC(3)
      P(N+1,4)=DPC(4)
      P(N+1,5)=PECM

C...Form two particles from flavours of lowest-mass system, if feasible.
      K(N+2,1)=1
      K(N+3,1)=1
      IF(MSTU(16).NE.2) THEN
        K(N+2,3)=N+1
        K(N+3,3)=N+1
      ELSE
        K(N+2,3)=IC1
        K(N+3,3)=IC2
      ENDIF
      K(N+2,4)=0
      K(N+3,4)=0
      K(N+2,5)=0
      K(N+3,5)=0
      IF(IABS(K(IC1,2)).NE.21) THEN
        KC1=LUCOMP(K(IC1,2))
        KC2=LUCOMP(K(IC2,2))
        IF(KC1.EQ.0.OR.KC2.EQ.0) GOTO 320
        KQ1=KCHG(KC1,2)*ISIGN(1,K(IC1,2))
        KQ2=KCHG(KC2,2)*ISIGN(1,K(IC2,2))
        IF(KQ1+KQ2.NE.0) GOTO 320
  200   CALL LUKFDI(K(IC1,2),0,KFLN,K(N+2,2))
        CALL LUKFDI(K(IC2,2),-KFLN,KFLDMP,K(N+3,2))
        IF(K(N+2,2).EQ.0.OR.K(N+3,2).EQ.0) GOTO 200
      ELSE
        IF(IABS(K(IC2,2)).NE.21) GOTO 320
  210   CALL LUKFDI(1+INT((2.+PARJ(2))*RLU(0)),0,KFLN,KFDMP)
        CALL LUKFDI(KFLN,0,KFLM,K(N+2,2))
        CALL LUKFDI(-KFLN,-KFLM,KFLDMP,K(N+3,2))
        IF(K(N+2,2).EQ.0.OR.K(N+3,2).EQ.0) GOTO 210
      ENDIF
      P(N+2,5)=ULMASS(K(N+2,2))
      P(N+3,5)=ULMASS(K(N+3,2))
      IF(P(N+2,5)+P(N+3,5)+PARJ(64).GE.PECM.AND.NSIN.EQ.1) GOTO 320
      IF(P(N+2,5)+P(N+3,5)+PARJ(64).GE.PECM) GOTO 260

C...Perform two-particle decay of jet system, if possible.
      IF(PECM.GE.0.02*DPC(4)) THEN
        PA=SQRT((PECM**2-(P(N+2,5)+P(N+3,5))**2)*(PECM**2-
     &  (P(N+2,5)-P(N+3,5))**2))/(2.*PECM)
        UE(3)=2.*RLU(0)-1.
        PHI=PARU(2)*RLU(0)
        UE(1)=SQRT(1.-UE(3)**2)*COS(PHI)
        UE(2)=SQRT(1.-UE(3)**2)*SIN(PHI)
        DO 220 J=1,3
        P(N+2,J)=PA*UE(J)
  220   P(N+3,J)=-PA*UE(J)
        P(N+2,4)=SQRT(PA**2+P(N+2,5)**2)
        P(N+3,4)=SQRT(PA**2+P(N+3,5)**2)
        MSTU(33)=1
        CALL LUDBRB(N+2,N+3,0.,0.,DPC(1)/DPC(4),DPC(2)/DPC(4),
     &  DPC(3)/DPC(4))
      ELSE
        NP=0
        DO 230 I=IC1,IC2
  230   IF(K(I,1).EQ.1.OR.K(I,1).EQ.2) NP=NP+1
        HA=P(IC1,4)*P(IC2,4)-P(IC1,1)*P(IC2,1)-P(IC1,2)*P(IC2,2)-
     &  P(IC1,3)*P(IC2,3)
        IF(NP.GE.3.OR.HA.LE.1.25*P(IC1,5)*P(IC2,5)) GOTO 260
        HD1=0.5*(P(N+2,5)**2-P(IC1,5)**2)
        HD2=0.5*(P(N+3,5)**2-P(IC2,5)**2)
        HR=SQRT(MAX(0.,((HA-HD1-HD2)**2-(P(N+2,5)*P(N+3,5))**2)/
     &  (HA**2-(P(IC1,5)*P(IC2,5))**2)))-1.
        HC=P(IC1,5)**2+2.*HA+P(IC2,5)**2
        HK1=((P(IC2,5)**2+HA)*HR+HD1-HD2)/HC
        HK2=((P(IC1,5)**2+HA)*HR+HD2-HD1)/HC
        DO 240 J=1,4
        P(N+2,J)=(1.+HK1)*P(IC1,J)-HK2*P(IC2,J)
  240   P(N+3,J)=(1.+HK2)*P(IC2,J)-HK1*P(IC1,J)
      ENDIF
      DO 250 J=1,4
      V(N+1,J)=V(IC1,J)
      V(N+2,J)=V(IC1,J)
  250 V(N+3,J)=V(IC2,J)
      V(N+1,5)=0.
      V(N+2,5)=0.
      V(N+3,5)=0.
      N=N+3
      GOTO 300

C...Else form one particle from the flavours available, if possible.
  260 K(N+1,5)=N+2
      IF(IABS(K(IC1,2)).GT.100.AND.IABS(K(IC2,2)).GT.100) THEN
        GOTO 320
      ELSEIF(IABS(K(IC1,2)).NE.21) THEN
        CALL LUKFDI(K(IC1,2),K(IC2,2),KFLDMP,K(N+2,2))
      ELSE
        KFLN=1+INT((2.+PARJ(2))*RLU(0))
        CALL LUKFDI(KFLN,-KFLN,KFLDMP,K(N+2,2))
      ENDIF
      IF(K(N+2,2).EQ.0) GOTO 260
      P(N+2,5)=ULMASS(K(N+2,2))

C...Find parton/particle which combines to largest extra mass.
      IR=0
      HA=0.
      HSM=0.
      DO 280 MCOMB=1,3
      IF(IR.NE.0) GOTO 280
      DO 270 I=MAX(1,IP),N
      IF(K(I,1).LE.0.OR.K(I,1).GT.10.OR.(I.GE.IC1.AND.I.LE.IC2.
     &AND.K(I,1).GE.1.AND.K(I,1).LE.2)) GOTO 270
      IF(MCOMB.EQ.1) KCI=LUCOMP(K(I,2))
      IF(MCOMB.EQ.1.AND.KCI.EQ.0) GOTO 270
      IF(MCOMB.EQ.1.AND.KCHG(KCI,2).EQ.0.AND.I.LE.NS) GOTO 270
      IF(MCOMB.EQ.2.AND.IABS(K(I,2)).GT.10.AND.IABS(K(I,2)).LE.100)
     &GOTO 270
      HCR=DPC(4)*P(I,4)-DPC(1)*P(I,1)-DPC(2)*P(I,2)-DPC(3)*P(I,3)
      HSR=2.*HCR+PECM**2-P(N+2,5)**2-2.*P(N+2,5)*P(I,5)
      IF(HSR.GT.HSM) THEN
        IR=I
        HA=HCR
        HSM=HSR
      ENDIF
  270 CONTINUE
  280 CONTINUE

C...Shuffle energy and momentum to put new particle on mass shell.
      IF(IR.NE.0) THEN
        HB=PECM**2+HA
        HC=P(N+2,5)**2+HA
        HD=P(IR,5)**2+HA
        HK2=0.5*(HB*SQRT(MAX(0.,((HB+HC)**2-4.*(HB+HD)*P(N+2,5)**2)/
     &  (HA**2-(PECM*P(IR,5))**2)))-(HB+HC))/(HB+HD)
        HK1=(0.5*(P(N+2,5)**2-PECM**2)+HD*HK2)/HB
        DO 290 J=1,4
        P(N+2,J)=(1.+HK1)*DPC(J)-HK2*P(IR,J)
        P(IR,J)=(1.+HK2)*P(IR,J)-HK1*DPC(J)
        V(N+1,J)=V(IC1,J)
  290   V(N+2,J)=V(IC1,J)
        V(N+1,5)=0.
        V(N+2,5)=0.
        N=N+2
      ELSE
        CALL LUERRM(3,'(LUPREP:) no match for collapsing cluster')
        RETURN
      ENDIF

C...Mark collapsed system and store daughter pointers. Iterate.
  300 DO 310 I=IC1,IC2
      IF((K(I,1).EQ.1.OR.K(I,1).EQ.2).AND.KCHG(LUCOMP(K(I,2)),2).NE.0)
     &THEN
        K(I,1)=K(I,1)+10
        IF(MSTU(16).NE.2) THEN
          K(I,4)=NSAV+1
          K(I,5)=NSAV+1
        ELSE
          K(I,4)=NSAV+2
          K(I,5)=N
        ENDIF
      ENDIF
  310 CONTINUE
      IF(N.LT.MSTU(4)-MSTU(32)-5) GOTO 140

C...Check flavours and invariant masses in parton systems.
  320 NP=0
      KFN=0
      KQS=0
      DO 330 J=1,5
  330 DPS(J)=0.
      DO 360 I=MAX(1,IP),N
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 360
      KC=LUCOMP(K(I,2))
      IF(KC.EQ.0) GOTO 360
      KQ=KCHG(KC,2)*ISIGN(1,K(I,2))
      IF(KQ.EQ.0) GOTO 360
      NP=NP+1
      IF(KQ.NE.2) THEN
        KFN=KFN+1
        KQS=KQS+KQ
        MSTJ(93)=1
        DPS(5)=DPS(5)+ULMASS(K(I,2))
      ENDIF
      DO 340 J=1,4
  340 DPS(J)=DPS(J)+P(I,J)
      IF(K(I,1).EQ.1) THEN
        IF(NP.NE.1.AND.(KFN.EQ.1.OR.KFN.GE.3.OR.KQS.NE.0)) CALL
     &  LUERRM(2,'(LUPREP:) unphysical flavour combination')
        IF(NP.NE.1.AND.DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2.LT.
     &  (0.9*PARJ(32)+DPS(5))**2) CALL LUERRM(3,
     &  '(LUPREP:) too small mass in jet system')
        NP=0
        KFN=0
        KQS=0
        DO 350 J=1,5
  350   DPS(J)=0.
      ENDIF
  360 CONTINUE

      RETURN
      END
