 
C********************************************************************* 
 
      SUBROUTINE LUCELL(NJET) 
 
C...Purpose: to provide a simple way of jet finding in an eta-phi-ET 
C...coordinate frame, as used for calorimeters at hadron colliders. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
 
C...Loop over all particles. Find cell that was hit by given particle. 
      PTLRAT=1./SINH(PARU(51))**2 
      NP=0 
      NC=N 
      DO 110 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 110 
      IF(P(I,1)**2+P(I,2)**2.LE.PTLRAT*P(I,3)**2) GOTO 110 
      IF(MSTU(41).GE.2) THEN 
        KC=LUCOMP(K(I,2)) 
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &  KC.EQ.18) GOTO 110 
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &  GOTO 110 
      ENDIF 
      NP=NP+1 
      PT=SQRT(P(I,1)**2+P(I,2)**2) 
      ETA=SIGN(LOG((SQRT(PT**2+P(I,3)**2)+ABS(P(I,3)))/PT),P(I,3)) 
      IETA=MAX(1,MIN(MSTU(51),1+INT(MSTU(51)*0.5*(ETA/PARU(51)+1.)))) 
      PHI=ULANGL(P(I,1),P(I,2)) 
      IPHI=MAX(1,MIN(MSTU(52),1+INT(MSTU(52)*0.5*(PHI/PARU(1)+1.)))) 
      IETPH=MSTU(52)*IETA+IPHI 
 
C...Add to cell already hit, or book new cell. 
      DO 100 IC=N+1,NC 
      IF(IETPH.EQ.K(IC,3)) THEN 
        K(IC,4)=K(IC,4)+1 
        P(IC,5)=P(IC,5)+PT 
        GOTO 110 
      ENDIF 
  100 CONTINUE 
      IF(NC.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUCELL:) no more memory left in LUJETS') 
        NJET=-2 
        RETURN 
      ENDIF 
      NC=NC+1 
      K(NC,3)=IETPH 
      K(NC,4)=1 
      K(NC,5)=2 
      P(NC,1)=(PARU(51)/MSTU(51))*(2*IETA-1-MSTU(51)) 
      P(NC,2)=(PARU(1)/MSTU(52))*(2*IPHI-1-MSTU(52)) 
      P(NC,5)=PT 
  110 CONTINUE 
 
C...Smear true bin content by calorimeter resolution. 
      IF(MSTU(53).GE.1) THEN 
        DO 130 IC=N+1,NC 
        PEI=P(IC,5) 
        IF(MSTU(53).EQ.2) PEI=P(IC,5)*COSH(P(IC,1)) 
  120   PEF=PEI+PARU(55)*SQRT(-2.*LOG(MAX(1E-10,RLU(0)))*PEI)* 
     &  COS(PARU(2)*RLU(0)) 
        IF(PEF.LT.0..OR.PEF.GT.PARU(56)*PEI) GOTO 120 
        P(IC,5)=PEF 
        IF(MSTU(53).EQ.2) P(IC,5)=PEF/COSH(P(IC,1)) 
  130   CONTINUE 
      ENDIF 
 
C...Remove cells below threshold. 
      IF(PARU(58).GT.0.) THEN 
        NCC=NC 
        NC=N 
        DO 140 IC=N+1,NCC 
        IF(P(IC,5).GT.PARU(58)) THEN 
          NC=NC+1 
          K(NC,3)=K(IC,3) 
          K(NC,4)=K(IC,4) 
          K(NC,5)=K(IC,5) 
          P(NC,1)=P(IC,1) 
          P(NC,2)=P(IC,2) 
          P(NC,5)=P(IC,5) 
        ENDIF 
  140   CONTINUE 
      ENDIF 
 
C...Find initiator cell: the one with highest pT of not yet used ones. 
      NJ=NC 
  150 ETMAX=0. 
      DO 160 IC=N+1,NC 
      IF(K(IC,5).NE.2) GOTO 160 
      IF(P(IC,5).LE.ETMAX) GOTO 160 
      ICMAX=IC 
      ETA=P(IC,1) 
      PHI=P(IC,2) 
      ETMAX=P(IC,5) 
  160 CONTINUE 
      IF(ETMAX.LT.PARU(52)) GOTO 220 
      IF(NJ.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUCELL:) no more memory left in LUJETS') 
        NJET=-2 
        RETURN 
      ENDIF 
      K(ICMAX,5)=1 
      NJ=NJ+1 
      K(NJ,4)=0 
      K(NJ,5)=1 
      P(NJ,1)=ETA 
      P(NJ,2)=PHI 
      P(NJ,3)=0. 
      P(NJ,4)=0. 
      P(NJ,5)=0. 
 
C...Sum up unused cells within required distance of initiator. 
      DO 170 IC=N+1,NC 
      IF(K(IC,5).EQ.0) GOTO 170 
      IF(ABS(P(IC,1)-ETA).GT.PARU(54)) GOTO 170 
      DPHIA=ABS(P(IC,2)-PHI) 
      IF(DPHIA.GT.PARU(54).AND.DPHIA.LT.PARU(2)-PARU(54)) GOTO 170 
      PHIC=P(IC,2) 
      IF(DPHIA.GT.PARU(1)) PHIC=PHIC+SIGN(PARU(2),PHI) 
      IF((P(IC,1)-ETA)**2+(PHIC-PHI)**2.GT.PARU(54)**2) GOTO 170 
      K(IC,5)=-K(IC,5) 
      K(NJ,4)=K(NJ,4)+K(IC,4) 
      P(NJ,3)=P(NJ,3)+P(IC,5)*P(IC,1) 
      P(NJ,4)=P(NJ,4)+P(IC,5)*PHIC 
      P(NJ,5)=P(NJ,5)+P(IC,5) 
  170 CONTINUE 
 
C...Reject cluster below minimum ET, else accept. 
      IF(P(NJ,5).LT.PARU(53)) THEN 
        NJ=NJ-1 
        DO 180 IC=N+1,NC 
        IF(K(IC,5).LT.0) K(IC,5)=-K(IC,5) 
  180   CONTINUE 
      ELSEIF(MSTU(54).LE.2) THEN 
        P(NJ,3)=P(NJ,3)/P(NJ,5) 
        P(NJ,4)=P(NJ,4)/P(NJ,5) 
        IF(ABS(P(NJ,4)).GT.PARU(1)) P(NJ,4)=P(NJ,4)-SIGN(PARU(2), 
     &  P(NJ,4)) 
        DO 190 IC=N+1,NC 
        IF(K(IC,5).LT.0) K(IC,5)=0 
  190   CONTINUE 
      ELSE 
        DO 200 J=1,4 
        P(NJ,J)=0. 
  200   CONTINUE 
        DO 210 IC=N+1,NC 
        IF(K(IC,5).GE.0) GOTO 210 
        P(NJ,1)=P(NJ,1)+P(IC,5)*COS(P(IC,2)) 
        P(NJ,2)=P(NJ,2)+P(IC,5)*SIN(P(IC,2)) 
        P(NJ,3)=P(NJ,3)+P(IC,5)*SINH(P(IC,1)) 
        P(NJ,4)=P(NJ,4)+P(IC,5)*COSH(P(IC,1)) 
        K(IC,5)=0 
  210   CONTINUE 
      ENDIF 
      GOTO 150 
 
C...Arrange clusters in falling ET sequence. 
  220 DO 250 I=1,NJ-NC 
      ETMAX=0. 
      DO 230 IJ=NC+1,NJ 
      IF(K(IJ,5).EQ.0) GOTO 230 
      IF(P(IJ,5).LT.ETMAX) GOTO 230 
      IJMAX=IJ 
      ETMAX=P(IJ,5) 
  230 CONTINUE 
      K(IJMAX,5)=0 
      K(N+I,1)=31 
      K(N+I,2)=98 
      K(N+I,3)=I 
      K(N+I,4)=K(IJMAX,4) 
      K(N+I,5)=0 
      DO 240 J=1,5 
      P(N+I,J)=P(IJMAX,J) 
      V(N+I,J)=0. 
  240 CONTINUE 
  250 CONTINUE 
      NJET=NJ-NC 
 
C...Convert to massless or massive four-vectors. 
      IF(MSTU(54).EQ.2) THEN 
        DO 260 I=N+1,N+NJET 
        ETA=P(I,3) 
        P(I,1)=P(I,5)*COS(P(I,4)) 
        P(I,2)=P(I,5)*SIN(P(I,4)) 
        P(I,3)=P(I,5)*SINH(ETA) 
        P(I,4)=P(I,5)*COSH(ETA) 
        P(I,5)=0. 
  260   CONTINUE 
      ELSEIF(MSTU(54).GE.3) THEN 
        DO 270 I=N+1,N+NJET 
        P(I,5)=SQRT(MAX(0.,P(I,4)**2-P(I,1)**2-P(I,2)**2-P(I,3)**2)) 
  270   CONTINUE 
      ENDIF 
 
C...Information about storage. 
      MSTU(61)=N+1 
      MSTU(62)=NP 
      MSTU(63)=NC-N 
      IF(MSTU(43).LE.1) MSTU(3)=NJET 
      IF(MSTU(43).GE.2) N=N+NJET 
 
      RETURN 
      END 
