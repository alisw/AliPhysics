 
C********************************************************************* 
 
      SUBROUTINE LUSPHE(SPH,APL) 
 
C...Purpose: to perform sphericity tensor analysis to give sphericity, 
C...aplanarity and the related event axes. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION SM(3,3),SV(3,3) 
 
C...Calculate matrix to be diagonalized. 
      NP=0 
      DO 110 J1=1,3 
      DO 100 J2=J1,3 
      SM(J1,J2)=0. 
  100 CONTINUE 
  110 CONTINUE 
      PS=0. 
      DO 140 I=1,N 
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 140 
      IF(MSTU(41).GE.2) THEN 
        KC=LUCOMP(K(I,2)) 
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &  KC.EQ.18) GOTO 140 
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &  GOTO 140 
      ENDIF 
      NP=NP+1 
      PA=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) 
      PWT=1. 
      IF(ABS(PARU(41)-2.).GT.0.001) PWT=MAX(1E-10,PA)**(PARU(41)-2.) 
      DO 130 J1=1,3 
      DO 120 J2=J1,3 
      SM(J1,J2)=SM(J1,J2)+PWT*P(I,J1)*P(I,J2) 
  120 CONTINUE 
  130 CONTINUE 
      PS=PS+PWT*PA**2 
  140 CONTINUE 
 
C...Very low multiplicities (0 or 1) not considered. 
      IF(NP.LE.1) THEN 
        CALL LUERRM(8,'(LUSPHE:) too few particles for analysis') 
        SPH=-1. 
        APL=-1. 
        RETURN 
      ENDIF 
      DO 160 J1=1,3 
      DO 150 J2=J1,3 
      SM(J1,J2)=SM(J1,J2)/PS 
  150 CONTINUE 
  160 CONTINUE 
 
C...Find eigenvalues to matrix (third degree equation). 
      SQ=(SM(1,1)*SM(2,2)+SM(1,1)*SM(3,3)+SM(2,2)*SM(3,3)-SM(1,2)**2- 
     &SM(1,3)**2-SM(2,3)**2)/3.-1./9. 
      SR=-0.5*(SQ+1./9.+SM(1,1)*SM(2,3)**2+SM(2,2)*SM(1,3)**2+SM(3,3)* 
     &SM(1,2)**2-SM(1,1)*SM(2,2)*SM(3,3))+SM(1,2)*SM(1,3)*SM(2,3)+1./27. 
      SP=COS(ACOS(MAX(MIN(SR/SQRT(-SQ**3),1.),-1.))/3.) 
      P(N+1,4)=1./3.+SQRT(-SQ)*MAX(2.*SP,SQRT(3.*(1.-SP**2))-SP) 
      P(N+3,4)=1./3.+SQRT(-SQ)*MIN(2.*SP,-SQRT(3.*(1.-SP**2))-SP) 
      P(N+2,4)=1.-P(N+1,4)-P(N+3,4) 
      IF(P(N+2,4).LT.1E-5) THEN 
        CALL LUERRM(8,'(LUSPHE:) all particles back-to-back') 
        SPH=-1. 
        APL=-1. 
        RETURN 
      ENDIF 
 
C...Find first and last eigenvector by solving equation system. 
      DO 240 I=1,3,2 
      DO 180 J1=1,3 
      SV(J1,J1)=SM(J1,J1)-P(N+I,4) 
      DO 170 J2=J1+1,3 
      SV(J1,J2)=SM(J1,J2) 
      SV(J2,J1)=SM(J1,J2) 
  170 CONTINUE 
  180 CONTINUE 
      SMAX=0. 
      DO 200 J1=1,3 
      DO 190 J2=1,3 
      IF(ABS(SV(J1,J2)).LE.SMAX) GOTO 190 
      JA=J1 
      JB=J2 
      SMAX=ABS(SV(J1,J2)) 
  190 CONTINUE 
  200 CONTINUE 
      SMAX=0. 
      DO 220 J3=JA+1,JA+2 
      J1=J3-3*((J3-1)/3) 
      RL=SV(J1,JB)/SV(JA,JB) 
      DO 210 J2=1,3 
      SV(J1,J2)=SV(J1,J2)-RL*SV(JA,J2) 
      IF(ABS(SV(J1,J2)).LE.SMAX) GOTO 210 
      JC=J1 
      SMAX=ABS(SV(J1,J2)) 
  210 CONTINUE 
  220 CONTINUE 
      JB1=JB+1-3*(JB/3) 
      JB2=JB+2-3*((JB+1)/3) 
      P(N+I,JB1)=-SV(JC,JB2) 
      P(N+I,JB2)=SV(JC,JB1) 
      P(N+I,JB)=-(SV(JA,JB1)*P(N+I,JB1)+SV(JA,JB2)*P(N+I,JB2))/ 
     &SV(JA,JB) 
      PA=SQRT(P(N+I,1)**2+P(N+I,2)**2+P(N+I,3)**2) 
      SGN=(-1.)**INT(RLU(0)+0.5) 
      DO 230 J=1,3 
      P(N+I,J)=SGN*P(N+I,J)/PA 
  230 CONTINUE 
  240 CONTINUE 
 
C...Middle axis orthogonal to other two. Fill other codes. 
      SGN=(-1.)**INT(RLU(0)+0.5) 
      P(N+2,1)=SGN*(P(N+1,2)*P(N+3,3)-P(N+1,3)*P(N+3,2)) 
      P(N+2,2)=SGN*(P(N+1,3)*P(N+3,1)-P(N+1,1)*P(N+3,3)) 
      P(N+2,3)=SGN*(P(N+1,1)*P(N+3,2)-P(N+1,2)*P(N+3,1)) 
      DO 260 I=1,3 
      K(N+I,1)=31 
      K(N+I,2)=95 
      K(N+I,3)=I 
      K(N+I,4)=0 
      K(N+I,5)=0 
      P(N+I,5)=0. 
      DO 250 J=1,5 
      V(I,J)=0. 
  250 CONTINUE 
  260 CONTINUE 
 
C...Calculate sphericity and aplanarity. Select storing option. 
      SPH=1.5*(P(N+2,4)+P(N+3,4)) 
      APL=1.5*P(N+3,4) 
      MSTU(61)=N+1 
      MSTU(62)=NP 
      IF(MSTU(43).LE.1) MSTU(3)=3 
      IF(MSTU(43).GE.2) N=N+3 
 
      RETURN 
      END 
