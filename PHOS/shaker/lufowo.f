*CMZ :          17/07/98  15.44.34  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUFOWO(H10,H20,H30,H40)

C...Purpose: to calculate the first few Fox-Wolfram moments.
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

C...Copy momenta for particles and calculate H0.
      NP=0
      H0=0.
      HD=0.
      DO 110 I=1,N
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 110
      IF(MSTU(41).GE.2) THEN
        KC=LUCOMP(K(I,2))
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.
     &  KC.EQ.18) GOTO 110
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0)
     &  GOTO 110
      ENDIF
      IF(N+NP.GE.MSTU(4)-MSTU(32)-5) THEN
        CALL LUERRM(11,'(LUFOWO:) no more memory left in LUJETS')
        H10=-1.
        H20=-1.
        H30=-1.
        H40=-1.
        RETURN
      ENDIF
      NP=NP+1
      DO 100 J=1,3
  100 P(N+NP,J)=P(I,J)
      P(N+NP,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
      H0=H0+P(N+NP,4)
      HD=HD+P(N+NP,4)**2
  110 CONTINUE
      H0=H0**2

C...Very low multiplicities (0 or 1) not considered.
      IF(NP.LE.1) THEN
        CALL LUERRM(8,'(LUFOWO:) too few particles for analysis')
        H10=-1.
        H20=-1.
        H30=-1.
        H40=-1.
        RETURN
      ENDIF

C...Calculate H1 - H4.
      H10=0.
      H20=0.
      H30=0.
      H40=0.
      DO 120 I1=N+1,N+NP
      DO 120 I2=I1+1,N+NP
      CTHE=(P(I1,1)*P(I2,1)+P(I1,2)*P(I2,2)+P(I1,3)*P(I2,3))/
     &(P(I1,4)*P(I2,4))
      H10=H10+P(I1,4)*P(I2,4)*CTHE
      H20=H20+P(I1,4)*P(I2,4)*(1.5*CTHE**2-0.5)
      H30=H30+P(I1,4)*P(I2,4)*(2.5*CTHE**3-1.5*CTHE)
      H40=H40+P(I1,4)*P(I2,4)*(4.375*CTHE**4-3.75*CTHE**2+0.375)
  120 CONTINUE

C...Calculate H1/H0 - H4/H0. Output.
      MSTU(61)=N+1
      MSTU(62)=NP
      H10=(HD+2.*H10)/H0
      H20=(HD+2.*H20)/H0
      H30=(HD+2.*H30)/H0
      H40=(HD+2.*H40)/H0

      RETURN
      END
