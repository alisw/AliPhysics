*CMZ :          17/07/98  15.44.31  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LU2ENT(IP,KF1,KF2,PECM)

C...Purpose: to store two partons/particles in their CM frame,
C...with the first along the +z axis.
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

C...Standard checks.
      MSTU(28)=0
      IF(MSTU(12).GE.1) CALL LULIST(0)
      IPA=MAX(1,IABS(IP))
      IF(IPA.GT.MSTU(4)-1) CALL LUERRM(21,
     &'(LU2ENT:) writing outside LUJETS memory')
      KC1=LUCOMP(KF1)
      KC2=LUCOMP(KF2)
      IF(KC1.EQ.0.OR.KC2.EQ.0) CALL LUERRM(12,
     &'(LU2ENT:) unknown flavour code')

C...Find masses. Reset K, P and V vectors.
      PM1=0.
      IF(MSTU(10).EQ.1) PM1=P(IPA,5)
      IF(MSTU(10).GE.2) PM1=ULMASS(KF1)
      PM2=0.
      IF(MSTU(10).EQ.1) PM2=P(IPA+1,5)
      IF(MSTU(10).GE.2) PM2=ULMASS(KF2)
      DO 100 I=IPA,IPA+1
      DO 100 J=1,5
      K(I,J)=0
      P(I,J)=0.
  100 V(I,J)=0.

C...Check flavours.
      KQ1=KCHG(KC1,2)*ISIGN(1,KF1)
      KQ2=KCHG(KC2,2)*ISIGN(1,KF2)
      IF(KQ1+KQ2.NE.0.AND.KQ1+KQ2.NE.4) CALL LUERRM(2,
     &'(LU2ENT:) unphysical flavour combination')
      K(IPA,2)=KF1
      K(IPA+1,2)=KF2

C...Store partons/particles in K vectors for normal case.
      IF(IP.GE.0) THEN
        K(IPA,1)=1
        IF(KQ1.NE.0.AND.KQ2.NE.0) K(IPA,1)=2
        K(IPA+1,1)=1

C...Store partons in K vectors for parton shower evolution.
      ELSE
        IF(KQ1.EQ.0.OR.KQ2.EQ.0) CALL LUERRM(2,
     &  '(LU2ENT:) requested flavours can not develop parton shower')
        K(IPA,1)=3
        K(IPA+1,1)=3
        K(IPA,4)=MSTU(5)*(IPA+1)
        K(IPA,5)=K(IPA,4)
        K(IPA+1,4)=MSTU(5)*IPA
        K(IPA+1,5)=K(IPA+1,4)
      ENDIF

C...Check kinematics and store partons/particles in P vectors.
      IF(PECM.LE.PM1+PM2) CALL LUERRM(13,
     &'(LU2ENT:) energy smaller than sum of masses')
      PA=SQRT(MAX(0.,(PECM**2-PM1**2-PM2**2)**2-(2.*PM1*PM2)**2))/
     &(2.*PECM)
      P(IPA,3)=PA
      P(IPA,4)=SQRT(PM1**2+PA**2)
      P(IPA,5)=PM1
      P(IPA+1,3)=-PA
      P(IPA+1,4)=SQRT(PM2**2+PA**2)
      P(IPA+1,5)=PM2

C...Set N. Optionally fragment/decay.
      N=IPA+1
      IF(IP.EQ.0) CALL LUEXEC

      RETURN
      END
