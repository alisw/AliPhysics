 
C********************************************************************* 
 
      SUBROUTINE LUCLUS(NJET) 
 
C...Purpose: to subdivide the particle content of an event into 
C...jets/clusters. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION PS(5) 
      SAVE NSAV,NP,PS,PSS,RINIT,NPRE,NREM 
 
C...Functions: distance measure in pT or (pseudo)mass. 
      R2T(I1,I2)=(P(I1,5)*P(I2,5)-P(I1,1)*P(I2,1)-P(I1,2)*P(I2,2)- 
     &P(I1,3)*P(I2,3))*2.*P(I1,5)*P(I2,5)/(0.0001+P(I1,5)+P(I2,5))**2 
      R2M(I1,I2)=2.*P(I1,4)*P(I2,4)*(1.-(P(I1,1)*P(I2,1)+P(I1,2)* 
     &P(I2,2)+P(I1,3)*P(I2,3))/(P(I1,5)*P(I2,5))) 
 
C...If first time, reset. If reentering, skip preliminaries. 
      IF(MSTU(48).LE.0) THEN 
        NP=0 
        DO 100 J=1,5 
        PS(J)=0. 
  100   CONTINUE 
        PSS=0. 
      ELSE 
        NJET=NSAV 
        IF(MSTU(43).GE.2) N=N-NJET 
        DO 110 I=N+1,N+NJET 
        P(I,5)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
  110   CONTINUE 
        IF(MSTU(46).LE.3) R2ACC=PARU(44)**2 
        IF(MSTU(46).GE.4) R2ACC=PARU(45)*PS(5)**2 
        NLOOP=0 
        GOTO 300 
      ENDIF 
 
C...Find which particles are to be considered in cluster search. 
      DO 140 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 140 
      IF(MSTU(41).GE.2) THEN 
        KC=LUCOMP(K(I,2)) 
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &  KC.EQ.18) GOTO 140 
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &  GOTO 140 
      ENDIF 
      IF(N+2*NP.GE.MSTU(4)-MSTU(32)-5) THEN 
        CALL LUERRM(11,'(LUCLUS:) no more memory left in LUJETS') 
        NJET=-1 
        RETURN 
      ENDIF 
 
C...Take copy of these particles, with space left for jets later on. 
      NP=NP+1 
      K(N+NP,3)=I 
      DO 120 J=1,5 
      P(N+NP,J)=P(I,J) 
  120 CONTINUE 
      IF(MSTU(42).EQ.0) P(N+NP,5)=0. 
      IF(MSTU(42).EQ.1.AND.K(I,2).NE.22) P(N+NP,5)=PMAS(101,1) 
      P(N+NP,4)=SQRT(P(N+NP,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2) 
      P(N+NP,5)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
      DO 130 J=1,4 
      PS(J)=PS(J)+P(N+NP,J) 
  130 CONTINUE 
      PSS=PSS+P(N+NP,5) 
  140 CONTINUE 
      DO 160 I=N+1,N+NP 
      K(I+NP,3)=K(I,3) 
      DO 150 J=1,5 
      P(I+NP,J)=P(I,J) 
  150 CONTINUE 
  160 CONTINUE 
      PS(5)=SQRT(MAX(0.,PS(4)**2-PS(1)**2-PS(2)**2-PS(3)**2)) 
 
C...Very low multiplicities not considered. 
      IF(NP.LT.MSTU(47)) THEN 
        CALL LUERRM(8,'(LUCLUS:) too few particles for analysis') 
        NJET=-1 
        RETURN 
      ENDIF 
 
C...Find precluster configuration. If too few jets, make harder cuts. 
      NLOOP=0 
      IF(MSTU(46).LE.3) R2ACC=PARU(44)**2 
      IF(MSTU(46).GE.4) R2ACC=PARU(45)*PS(5)**2 
      RINIT=1.25*PARU(43) 
      IF(NP.LE.MSTU(47)+2) RINIT=0. 
  170 RINIT=0.8*RINIT 
      NPRE=0 
      NREM=NP 
      DO 180 I=N+NP+1,N+2*NP 
      K(I,4)=0 
  180 CONTINUE 
 
C...Sum up small momentum region. Jet if enough absolute momentum. 
      IF(MSTU(46).LE.2) THEN 
        DO 190 J=1,4 
        P(N+1,J)=0. 
  190   CONTINUE 
        DO 210 I=N+NP+1,N+2*NP 
        IF(P(I,5).GT.2.*RINIT) GOTO 210 
        NREM=NREM-1 
        K(I,4)=1 
        DO 200 J=1,4 
        P(N+1,J)=P(N+1,J)+P(I,J) 
  200   CONTINUE 
  210   CONTINUE 
        P(N+1,5)=SQRT(P(N+1,1)**2+P(N+1,2)**2+P(N+1,3)**2) 
        IF(P(N+1,5).GT.2.*RINIT) NPRE=1 
        IF(RINIT.GE.0.2*PARU(43).AND.NPRE+NREM.LT.MSTU(47)) GOTO 170 
        IF(NREM.EQ.0) GOTO 170 
      ENDIF 
 
C...Find fastest remaining particle. 
  220 NPRE=NPRE+1 
      PMAX=0. 
      DO 230 I=N+NP+1,N+2*NP 
      IF(K(I,4).NE.0.OR.P(I,5).LE.PMAX) GOTO 230 
      IMAX=I 
      PMAX=P(I,5) 
  230 CONTINUE 
      DO 240 J=1,5 
      P(N+NPRE,J)=P(IMAX,J) 
  240 CONTINUE 
      NREM=NREM-1 
      K(IMAX,4)=NPRE 
 
C...Sum up precluster around it according to pT separation. 
      IF(MSTU(46).LE.2) THEN 
        DO 260 I=N+NP+1,N+2*NP 
        IF(K(I,4).NE.0) GOTO 260 
        R2=R2T(I,IMAX) 
        IF(R2.GT.RINIT**2) GOTO 260 
        NREM=NREM-1 
        K(I,4)=NPRE 
        DO 250 J=1,4 
        P(N+NPRE,J)=P(N+NPRE,J)+P(I,J) 
  250   CONTINUE 
  260   CONTINUE 
        P(N+NPRE,5)=SQRT(P(N+NPRE,1)**2+P(N+NPRE,2)**2+P(N+NPRE,3)**2) 
 
C...Sum up precluster around it according to mass separation. 
      ELSE 
  270   IMIN=0 
        R2MIN=RINIT**2 
        DO 280 I=N+NP+1,N+2*NP 
        IF(K(I,4).NE.0) GOTO 280 
        R2=R2M(I,N+NPRE) 
        IF(R2.GE.R2MIN) GOTO 280 
        IMIN=I 
        R2MIN=R2 
  280   CONTINUE 
        IF(IMIN.NE.0) THEN 
          DO 290 J=1,4 
          P(N+NPRE,J)=P(N+NPRE,J)+P(IMIN,J) 
  290     CONTINUE 
          P(N+NPRE,5)=SQRT(P(N+NPRE,1)**2+P(N+NPRE,2)**2+P(N+NPRE,3)**2) 
          NREM=NREM-1 
          K(IMIN,4)=NPRE 
          GOTO 270 
        ENDIF 
      ENDIF 
 
C...Check if more preclusters to be found. Start over if too few. 
      IF(RINIT.GE.0.2*PARU(43).AND.NPRE+NREM.LT.MSTU(47)) GOTO 170 
      IF(NREM.GT.0) GOTO 220 
      NJET=NPRE 
 
C...Reassign all particles to nearest jet. Sum up new jet momenta. 
  300 TSAV=0. 
      PSJT=0. 
  310 IF(MSTU(46).LE.1) THEN 
        DO 330 I=N+1,N+NJET 
        DO 320 J=1,4 
        V(I,J)=0. 
  320   CONTINUE 
  330 CONTINUE 
        DO 360 I=N+NP+1,N+2*NP 
        R2MIN=PSS**2 
        DO 340 IJET=N+1,N+NJET 
        IF(P(IJET,5).LT.RINIT) GOTO 340 
        R2=R2T(I,IJET) 
        IF(R2.GE.R2MIN) GOTO 340 
        IMIN=IJET 
        R2MIN=R2 
  340   CONTINUE 
        K(I,4)=IMIN-N 
        DO 350 J=1,4 
        V(IMIN,J)=V(IMIN,J)+P(I,J) 
  350   CONTINUE 
  360   CONTINUE 
        PSJT=0. 
        DO 380 I=N+1,N+NJET 
        DO 370 J=1,4 
        P(I,J)=V(I,J) 
  370   CONTINUE 
        P(I,5)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
        PSJT=PSJT+P(I,5) 
  380   CONTINUE 
      ENDIF 
 
C...Find two closest jets. 
      R2MIN=2.*MAX(R2ACC,PS(5)**2) 
      DO 400 ITRY1=N+1,N+NJET-1 
      DO 390 ITRY2=ITRY1+1,N+NJET 
      IF(MSTU(46).LE.2) R2=R2T(ITRY1,ITRY2) 
      IF(MSTU(46).GE.3) R2=R2M(ITRY1,ITRY2) 
      IF(R2.GE.R2MIN) GOTO 390 
      IMIN1=ITRY1 
      IMIN2=ITRY2 
      R2MIN=R2 
  390 CONTINUE 
  400 CONTINUE 
 
C...If allowed, join two closest jets and start over. 
      IF(NJET.GT.MSTU(47).AND.R2MIN.LT.R2ACC) THEN 
        IREC=MIN(IMIN1,IMIN2) 
        IDEL=MAX(IMIN1,IMIN2) 
        DO 410 J=1,4 
        P(IREC,J)=P(IMIN1,J)+P(IMIN2,J) 
  410   CONTINUE 
        P(IREC,5)=SQRT(P(IREC,1)**2+P(IREC,2)**2+P(IREC,3)**2) 
        DO 430 I=IDEL+1,N+NJET 
        DO 420 J=1,5 
        P(I-1,J)=P(I,J) 
  420   CONTINUE 
  430 CONTINUE 
        IF(MSTU(46).GE.2) THEN 
          DO 440 I=N+NP+1,N+2*NP 
          IORI=N+K(I,4) 
          IF(IORI.EQ.IDEL) K(I,4)=IREC-N 
          IF(IORI.GT.IDEL) K(I,4)=K(I,4)-1 
  440     CONTINUE 
        ENDIF 
        NJET=NJET-1 
        GOTO 300 
 
C...Divide up broad jet if empty cluster in list of final ones. 
      ELSEIF(NJET.EQ.MSTU(47).AND.MSTU(46).LE.1.AND.NLOOP.LE.2) THEN 
        DO 450 I=N+1,N+NJET 
        K(I,5)=0 
  450   CONTINUE 
        DO 460 I=N+NP+1,N+2*NP 
        K(N+K(I,4),5)=K(N+K(I,4),5)+1 
  460   CONTINUE 
        IEMP=0 
        DO 470 I=N+1,N+NJET 
        IF(K(I,5).EQ.0) IEMP=I 
  470   CONTINUE 
        IF(IEMP.NE.0) THEN 
          NLOOP=NLOOP+1 
          ISPL=0 
          R2MAX=0. 
          DO 480 I=N+NP+1,N+2*NP 
          IF(K(N+K(I,4),5).LE.1.OR.P(I,5).LT.RINIT) GOTO 480 
          IJET=N+K(I,4) 
          R2=R2T(I,IJET) 
          IF(R2.LE.R2MAX) GOTO 480 
          ISPL=I 
          R2MAX=R2 
  480     CONTINUE 
          IF(ISPL.NE.0) THEN 
            IJET=N+K(ISPL,4) 
            DO 490 J=1,4 
            P(IEMP,J)=P(ISPL,J) 
            P(IJET,J)=P(IJET,J)-P(ISPL,J) 
  490       CONTINUE 
            P(IEMP,5)=P(ISPL,5) 
            P(IJET,5)=SQRT(P(IJET,1)**2+P(IJET,2)**2+P(IJET,3)**2) 
            IF(NLOOP.LE.2) GOTO 300 
          ENDIF 
        ENDIF 
      ENDIF 
 
C...If generalized thrust has not yet converged, continue iteration. 
      IF(MSTU(46).LE.1.AND.NLOOP.LE.2.AND.PSJT/PSS.GT.TSAV+PARU(48)) 
     &THEN 
        TSAV=PSJT/PSS 
        GOTO 310 
      ENDIF 
 
C...Reorder jets according to energy. 
      DO 510 I=N+1,N+NJET 
      DO 500 J=1,5 
      V(I,J)=P(I,J) 
  500 CONTINUE 
  510 CONTINUE 
      DO 540 INEW=N+1,N+NJET 
      PEMAX=0. 
      DO 520 ITRY=N+1,N+NJET 
      IF(V(ITRY,4).LE.PEMAX) GOTO 520 
      IMAX=ITRY 
      PEMAX=V(ITRY,4) 
  520 CONTINUE 
      K(INEW,1)=31 
      K(INEW,2)=97 
      K(INEW,3)=INEW-N 
      K(INEW,4)=0 
      DO 530 J=1,5 
      P(INEW,J)=V(IMAX,J) 
  530 CONTINUE 
      V(IMAX,4)=-1. 
      K(IMAX,5)=INEW 
  540 CONTINUE 
 
C...Clean up particle-jet assignments and jet information. 
      DO 550 I=N+NP+1,N+2*NP 
      IORI=K(N+K(I,4),5) 
      K(I,4)=IORI-N 
      IF(K(K(I,3),1).NE.3) K(K(I,3),4)=IORI-N 
      K(IORI,4)=K(IORI,4)+1 
  550 CONTINUE 
      IEMP=0 
      PSJT=0. 
      DO 570 I=N+1,N+NJET 
      K(I,5)=0 
      PSJT=PSJT+P(I,5) 
      P(I,5)=SQRT(MAX(P(I,4)**2-P(I,5)**2,0.)) 
      DO 560 J=1,5 
      V(I,J)=0. 
  560 CONTINUE 
      IF(K(I,4).EQ.0) IEMP=I 
  570 CONTINUE 
 
C...Select storing option. Output variables. Check for failure. 
      MSTU(61)=N+1 
      MSTU(62)=NP 
      MSTU(63)=NPRE 
      PARU(61)=PS(5) 
      PARU(62)=PSJT/PSS 
      PARU(63)=SQRT(R2MIN) 
      IF(NJET.LE.1) PARU(63)=0. 
      IF(IEMP.NE.0) THEN 
        CALL LUERRM(8,'(LUCLUS:) failed to reconstruct as requested') 
        NJET=-1 
      ENDIF 
      IF(MSTU(43).LE.1) MSTU(3)=NJET 
      IF(MSTU(43).GE.2) N=N+NJET 
      NSAV=NJET 
 
      RETURN 
      END 
