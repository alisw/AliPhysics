 
C********************************************************************* 
 
      SUBROUTINE LUEXEC 
 
C...Purpose: to administrate the fragmentation and decay chain. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/ 
      DIMENSION PS(2,6) 
 
C...Initialize and reset. 
      MSTU(24)=0 
      IF(MSTU(12).GE.1) CALL LULIST(0) 
      MSTU(31)=MSTU(31)+1 
      MSTU(1)=0 
      MSTU(2)=0 
      MSTU(3)=0 
      IF(MSTU(17).LE.0) MSTU(90)=0 
      MCONS=1 
 
C...Sum up momentum, energy and charge for starting entries. 
      NSAV=N 
      DO 110 I=1,2 
      DO 100 J=1,6 
      PS(I,J)=0. 
  100 CONTINUE 
  110 CONTINUE 
      DO 130 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 130 
      DO 120 J=1,4 
      PS(1,J)=PS(1,J)+P(I,J) 
  120 CONTINUE 
      PS(1,6)=PS(1,6)+LUCHGE(K(I,2)) 
  130 CONTINUE 
      PARU(21)=PS(1,4) 
 
C...Prepare system for subsequent fragmentation/decay. 
      CALL LUPREP(0) 
 
C...Loop through jet fragmentation and particle decays. 
      MBE=0 
  140 MBE=MBE+1 
      IP=0 
  150 IP=IP+1 
      KC=0 
      IF(K(IP,1).GT.0.AND.K(IP,1).LE.10) KC=LUCOMP(K(IP,2)) 
      IF(KC.EQ.0) THEN 
 
C...Particle decay if unstable and allowed. Save long-lived particle 
C...decays until second pass after Bose-Einstein effects. 
      ELSEIF(KCHG(KC,2).EQ.0) THEN 
        IF(MSTJ(21).GE.1.AND.MDCY(KC,1).GE.1.AND.(MSTJ(51).LE.0.OR.MBE 
     &  .EQ.2.OR.PMAS(KC,2).GE.PARJ(91).OR.IABS(K(IP,2)).EQ.311)) 
     &  CALL LUDECY(IP) 
 
C...Decay products may develop a shower. 
        IF(MSTJ(92).GT.0) THEN 
          IP1=MSTJ(92) 
          QMAX=SQRT(MAX(0.,(P(IP1,4)+P(IP1+1,4))**2-(P(IP1,1)+P(IP1+1, 
     &    1))**2-(P(IP1,2)+P(IP1+1,2))**2-(P(IP1,3)+P(IP1+1,3))**2)) 
          CALL LUSHOW(IP1,IP1+1,QMAX) 
          CALL LUPREP(IP1) 
          MSTJ(92)=0 
        ELSEIF(MSTJ(92).LT.0) THEN 
          IP1=-MSTJ(92) 
          CALL LUSHOW(IP1,-3,P(IP,5)) 
          CALL LUPREP(IP1) 
          MSTJ(92)=0 
        ENDIF 
 
C...Jet fragmentation: string or independent fragmentation. 
      ELSEIF(K(IP,1).EQ.1.OR.K(IP,1).EQ.2) THEN 
        MFRAG=MSTJ(1) 
        IF(MFRAG.GE.1.AND.K(IP,1).EQ.1) MFRAG=2 
        IF(MSTJ(21).GE.2.AND.K(IP,1).EQ.2.AND.N.GT.IP) THEN 
          IF(K(IP+1,1).EQ.1.AND.K(IP+1,3).EQ.K(IP,3).AND. 
     &    K(IP,3).GT.0.AND.K(IP,3).LT.IP) THEN 
            IF(KCHG(LUCOMP(K(K(IP,3),2)),2).EQ.0) MFRAG=MIN(1,MFRAG) 
          ENDIF 
        ENDIF 
        IF(MFRAG.EQ.1) CALL LUSTRF(IP) 
        IF(MFRAG.EQ.2) CALL LUINDF(IP) 
        IF(MFRAG.EQ.2.AND.K(IP,1).EQ.1) MCONS=0 
        IF(MFRAG.EQ.2.AND.(MSTJ(3).LE.0.OR.MOD(MSTJ(3),5).EQ.0)) MCONS=0 
      ENDIF 
 
C...Loop back if enough space left in LUJETS and no error abort. 
      IF(MSTU(24).NE.0.AND.MSTU(21).GE.2) THEN 
      ELSEIF(IP.LT.N.AND.N.LT.MSTU(4)-20-MSTU(32)) THEN 
        GOTO 150 
      ELSEIF(IP.LT.N) THEN 
        CALL LUERRM(11,'(LUEXEC:) no more memory left in LUJETS') 
      ENDIF 
 
C...Include simple Bose-Einstein effect parametrization if desired. 
      IF(MBE.EQ.1.AND.MSTJ(51).GE.1) THEN 
        CALL LUBOEI(NSAV) 
        GOTO 140 
      ENDIF 
 
C...Check that momentum, energy and charge were conserved. 
      DO 170 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 170 
      DO 160 J=1,4 
      PS(2,J)=PS(2,J)+P(I,J) 
  160 CONTINUE 
      PS(2,6)=PS(2,6)+LUCHGE(K(I,2)) 
  170 CONTINUE 
      PDEV=(ABS(PS(2,1)-PS(1,1))+ABS(PS(2,2)-PS(1,2))+ABS(PS(2,3)- 
     &PS(1,3))+ABS(PS(2,4)-PS(1,4)))/(1.+ABS(PS(2,4))+ABS(PS(1,4))) 
      IF(MCONS.EQ.1.AND.PDEV.GT.PARU(11)) CALL LUERRM(15, 
     &'(LUEXEC:) four-momentum was not conserved') 
      IF(MCONS.EQ.1.AND.ABS(PS(2,6)-PS(1,6)).GT.0.1) CALL LUERRM(15, 
     &'(LUEXEC:) charge was not conserved') 
 
      RETURN 
      END 
