*CMZ :          17/07/98  15.44.34  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUJMAS(PMH,PML)

C...Purpose: to determine, approximately, the two jet masses that
C...minimize the sum m_H^2 + m_L^2, a la Clavelli and Wyler.
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEND.
      DIMENSION SM(3,3),SAX(3),PS(3,5)

C...Reset.
      NP=0
      DO 110 J1=1,3
      DO 100 J2=J1,3
  100 SM(J1,J2)=0.
      DO 110 J2=1,4
  110 PS(J1,J2)=0.
      PSS=0.

C...Take copy of particles that are to be considered in mass analysis.
      DO 150 I=1,N
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 150
      IF(MSTU(41).GE.2) THEN
        KC=LUCOMP(K(I,2))
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.
     &  KC.EQ.18) GOTO 150
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0)
     &  GOTO 150
      ENDIF
      IF(N+NP+1.GE.MSTU(4)-MSTU(32)-5) THEN
        CALL LUERRM(11,'(LUJMAS:) no more memory left in LUJETS')
        PMH=-2.
        PML=-2.
        RETURN
      ENDIF
      NP=NP+1
      DO 120 J=1,5
  120 P(N+NP,J)=P(I,J)
      IF(MSTU(42).EQ.0) P(N+NP,5)=0.
      IF(MSTU(42).EQ.1.AND.K(I,2).NE.22) P(N+NP,5)=PMAS(101,1)
      P(N+NP,4)=SQRT(P(N+NP,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2)

C...Fill information in sphericity tensor and total momentum vector.
      DO 130 J1=1,3
      DO 130 J2=J1,3
  130 SM(J1,J2)=SM(J1,J2)+P(I,J1)*P(I,J2)
      PSS=PSS+(P(I,1)**2+P(I,2)**2+P(I,3)**2)
      DO 140 J=1,4
  140 PS(3,J)=PS(3,J)+P(N+NP,J)
  150 CONTINUE

C...Very low multiplicities (0 or 1) not considered.
      IF(NP.LE.1) THEN
        CALL LUERRM(8,'(LUJMAS:) too few particles for analysis')
        PMH=-1.
        PML=-1.
        RETURN
      ENDIF
      PARU(61)=SQRT(MAX(0.,PS(3,4)**2-PS(3,1)**2-PS(3,2)**2-PS(3,3)**2))

C...Find largest eigenvalue to matrix (third degree equation).
      DO 160 J1=1,3
      DO 160 J2=J1,3
  160 SM(J1,J2)=SM(J1,J2)/PSS
      SQ=(SM(1,1)*SM(2,2)+SM(1,1)*SM(3,3)+SM(2,2)*SM(3,3)-SM(1,2)**2-
     &SM(1,3)**2-SM(2,3)**2)/3.-1./9.
      SR=-0.5*(SQ+1./9.+SM(1,1)*SM(2,3)**2+SM(2,2)*SM(1,3)**2+SM(3,3)*
     &SM(1,2)**2-SM(1,1)*SM(2,2)*SM(3,3))+SM(1,2)*SM(1,3)*SM(2,3)+1./27.
      SP=COS(ACOS(MAX(MIN(SR/SQRT(-SQ**3),1.),-1.))/3.)
      SMA=1./3.+SQRT(-SQ)*MAX(2.*SP,SQRT(3.*(1.-SP**2))-SP)

C...Find largest eigenvector by solving equation system.
      DO 170 J1=1,3
      SM(J1,J1)=SM(J1,J1)-SMA
      DO 170 J2=J1+1,3
  170 SM(J2,J1)=SM(J1,J2)
      SMAX=0.
      DO 180 J1=1,3
      DO 180 J2=1,3
      IF(ABS(SM(J1,J2)).LE.SMAX) GOTO 180
      JA=J1
      JB=J2
      SMAX=ABS(SM(J1,J2))
  180 CONTINUE
      SMAX=0.
      DO 190 J3=JA+1,JA+2
      J1=J3-3*((J3-1)/3)
      RL=SM(J1,JB)/SM(JA,JB)
      DO 190 J2=1,3
      SM(J1,J2)=SM(J1,J2)-RL*SM(JA,J2)
      IF(ABS(SM(J1,J2)).LE.SMAX) GOTO 190
      JC=J1
      SMAX=ABS(SM(J1,J2))
  190 CONTINUE
      JB1=JB+1-3*(JB/3)
      JB2=JB+2-3*((JB+1)/3)
      SAX(JB1)=-SM(JC,JB2)
      SAX(JB2)=SM(JC,JB1)
      SAX(JB)=-(SM(JA,JB1)*SAX(JB1)+SM(JA,JB2)*SAX(JB2))/SM(JA,JB)

C...Divide particles into two initial clusters by hemisphere.
      DO 200 I=N+1,N+NP
      PSAX=P(I,1)*SAX(1)+P(I,2)*SAX(2)+P(I,3)*SAX(3)
      IS=1
      IF(PSAX.LT.0.) IS=2
      K(I,3)=IS
      DO 200 J=1,4
  200 PS(IS,J)=PS(IS,J)+P(I,J)
      PMS=MAX(1E-10,PS(1,4)**2-PS(1,1)**2-PS(1,2)**2-PS(1,3)**2)+
     &MAX(1E-10,PS(2,4)**2-PS(2,1)**2-PS(2,2)**2-PS(2,3)**2)

C...Reassign one particle at a time; find maximum decrease of m^2 sum.
  210 PMD=0.
      IM=0
      DO 220 J=1,4
  220 PS(3,J)=PS(1,J)-PS(2,J)
      DO 230 I=N+1,N+NP
      PPS=P(I,4)*PS(3,4)-P(I,1)*PS(3,1)-P(I,2)*PS(3,2)-P(I,3)*PS(3,3)
      IF(K(I,3).EQ.1) PMDI=2.*(P(I,5)**2-PPS)
      IF(K(I,3).EQ.2) PMDI=2.*(P(I,5)**2+PPS)
      IF(PMDI.LT.PMD) THEN
        PMD=PMDI
        IM=I
      ENDIF
  230 CONTINUE

C...Loop back if significant reduction in sum of m^2.
      IF(PMD.LT.-PARU(48)*PMS) THEN
        PMS=PMS+PMD
        IS=K(IM,3)
        DO 240 J=1,4
        PS(IS,J)=PS(IS,J)-P(IM,J)
  240   PS(3-IS,J)=PS(3-IS,J)+P(IM,J)
        K(IM,3)=3-IS
        GOTO 210
      ENDIF

C...Final masses and output.
      MSTU(61)=N+1
      MSTU(62)=NP
      PS(1,5)=SQRT(MAX(0.,PS(1,4)**2-PS(1,1)**2-PS(1,2)**2-PS(1,3)**2))
      PS(2,5)=SQRT(MAX(0.,PS(2,4)**2-PS(2,1)**2-PS(2,2)**2-PS(2,3)**2))
      PMH=MAX(PS(1,5),PS(2,5))
      PML=MIN(PS(1,5),PS(2,5))

      RETURN
      END
