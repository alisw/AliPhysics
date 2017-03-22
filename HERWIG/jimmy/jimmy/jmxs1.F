      DOUBLE PRECISION FUNCTION JMXS1( X1, X2, PT2, FLAG1, FLAG2 )
c ------------------------------------------------------------------
c  Purpose: The fully differential cross section dsigma/dx1dx2dt
c           for E-P resolved scattering WITHOUT eikonal integral.
c
c   Inputs: x1, x2  -> "Bjorken" x's of photon and proton,
c                       respectively
c           pt2     ->  pt2 of the partons
c           flag1   ->  If flag1=0 then no YGAMMA included,
c           flag2   ->  If flag2=0 calculate dxsec/d(pt2)
c
c Author: JRF/JMB
c   Mods: JRF Oct '93 - Include different parton processes.
c ------------------------------------------------------------------
#include "herwig65.inc"
#include "jimmy.inc"


      INTEGER I, J, NF, FLAG1, FLAG2, PID(4)
      DOUBLE PRECISION Y, PT2, X1, X2, S, T, U, Q2, G1
      DOUBLE PRECISION U1, UU1, D1, DD1, S1, SS1, G2, U2, UU2
      DOUBLE PRECISION D2, DD2, S2, SS2, SCL2, ALPHA, F, F1, F2
      DOUBLE PRECISION X, C1, CC1, C2, CC2, JMPTOT
      DOUBLE PRECISION SMIN, SCALE
      DOUBLE PRECISION DIST(13), HWUALF
      LOGICAL VETO

      IF (JMBUG.GT.3) THEN
        WRITE(JMOUT,8900) 'JMXS1 called'
      ENDIF

c     Centre-of-mass energy for this parton-parton collision.
      S = 2.D0*(EBEAM1*EBEAM2+PBEAM1*PBEAM2)*YGAMMA*X1*X2

c     Appropriate kinematics
      Y = 4.*PT2/S

c     Define Mandelstam variables in parton subprocess
      IF (Y .GE. .99D0) THEN
        IF (FLAG2.EQ.0) THEN
c         We are calculating dxsec/d(pt2) and are at the pole.
	  IF (JMBUG.GT.3) THEN
            WRITE(JMOUT,8900) 'JMXS1 at the pole. ZERO'
          ENDIF
          JMXS1 = 0.D0
          RETURN
        ELSE
          T = -S/2.D0
        ENDIF
      ELSE
        T = -S/2.*(1.0D0-(1.0D0-Y)**.5)
      ENDIF

      U=-S-T
c     Select Scale
      EMSCA = SQRT(TWO*S*T*U/(S*S+T*T+U*U))
      SCALE = EMSCA
      Q2    = EMSCA**2

      IF ((X1 .GE. 0.99D0).OR.(X2 .GE. 0.99D0)) THEN
	IF (JMBUG.GT.3) THEN
          WRITE(JMOUT,8900) 'JMXS1 zero: x out of upper range'
	ENDIF
        JMXS1=0.D0
        RETURN
      ENDIF
      IF ((X1 .EQ. 0.0D0).OR.(X2 .EQ. 0.D0)) THEN
	IF (JMBUG.GT.3) THEN
          WRITE(JMOUT,8900) 'JMXS1 zero: x out of lower range'
	ENDIF
        JMXS1=0.D0
        RETURN
      ENDIF

c     Apply detailed phase space cutoffs.
      CALL JMKNIF(X1,X2,PT2,VETO)
      IF (VETO) THEN
	IF (JMBUG.GT.3) THEN
          WRITE(JMOUT,8900) 'JMXS1 vetoed by jmknif'
	ENDIF
        JMXS1 = 0.D0
        RETURN
      ENDIF

c     BEAM1 pdfs, with applied cutoffs on x1
      X = X1
      IF (JMBUG.GT.3) THEN
        write(JMOUT,8900) 'x1,x,scale'
        write(JMOUT,*) x1,x,scale
      ENDIF
      IF (IPART1.EQ.121.OR.IPART1.EQ.127) THEN
C       photon PDF
        CALL HWSFUN(X,SCALE,59,NSTRU,DIST,1)
      ELSE
        CALL HWSFUN(X,SCALE,IPART1,NSTRU,DIST,1)
      ENDIF
      IF (JMBUG.GT.3) THEN
        write(JMOUT,8900) 'x1,x,scale,dist'
        write(JMOUT,*) x2,x,scale,dist
      ENDIF
      G1  = DIST(13)/X1*JMVETO(1,13)
      U1  = DIST(2)/X1*JMVETO(1,2)
      UU1 = DIST(8)/X1*JMVETO(1,8)
      D1  = DIST(1)/X1*JMVETO(1,1)
      DD1 = DIST(7)/X1*JMVETO(1,7)
      S1  = DIST(3)/X1*JMVETO(1,3)
      SS1 = DIST(9)/X1*JMVETO(1,9)
      C1  = DIST(4)/X1*JMVETO(1,4)
      CC1 = DIST(10)/X1*JMVETO(1,10)

c     BEAM2 pdfs, with applied cutoffs on x2
      X = X2
*      write(*,*) 'x2,x,scale',x2,x,scale
      IF (IPART2.EQ.121.OR.IPART2.EQ.127) THEN
C       photon PDF
        CALL HWSFUN(X,SCALE,59,NSTRU,DIST,2)
      ELSE
C       proton PDF
        CALL HWSFUN(X,SCALE,IPART2,NSTRU,DIST,2)
      ENDIF
*      write(*,*) 'x2,x,scale,dist',x2,x,scale,dist
      G2  = DIST(13)/X2*JMVETO(2,13)
      U2  = DIST(2)/X2*JMVETO(2,2)
      UU2 = DIST(8)/X2*JMVETO(2,8)
      D2  = DIST(1)/X2*JMVETO(2,1)
      DD2 = DIST(7)/X2*JMVETO(2,7)
      S2  = DIST(3)/X2*JMVETO(2,3)
      SS2 = DIST(9)/X2*JMVETO(2,9)
      C2  = DIST(4)/X2*JMVETO(2,4)
      CC2 = DIST(10)/X2*JMVETO(2,10)

C     Define LAMBDA_QCD, ALPHA_S and number of flavours.
c      NF = NFLAV
c      SCL2  = QCDLAM**2
c      ALPHA = 4.D0*PIFAC/(11.D0-2.D0/3.D0*NF)/DLOG(Q2/SCL2)
c      write(*,*) 'alpha 1',alpha
      ALPHA = HWUALF(1,EMSCA)
c      write(*,*) 'alpha 2',alpha

C********************************************************************
C**** JMPROC(1:117) CONTAINS THE PROBABILITIES FOR ALL 117         **
C**** SUBPROCESSES.                                                **
C**** NOTATION TO BE USED IN REST OF THIS COMMENT IS:              **
C**** P(i) = ij-kl means that P(i) is the probability for ij=>kl   **
C**** Where i,j,k and l are parton types (u,d,s,c for quarks,      **
C**** U,D,S,C for antiquarks and g for gluons) and i refers to     **
C**** partons originating from the photon whilst j refers to those **
C**** originating from the proton. Each probability is colour and  **
C**** spin averaged and the Mandelstam variable, t, is defined by  **
C**** the momenta of i and k (or j and l). The hard scale must be  **
C**** t-u symmetric since we always use t-u symmetrised cross      **
C**** sections (nb: p_T is t-u symmetric)                          **
C********************************************************************
C
C**** JMPROC(1) = gg-gg
C
C**** JMPROC(2) to JMPROC(5)  =  gg-qQ  (all equal)
C
C**** JMPROC(6) to JMPROC(13) = qq-qq, QQ-QQ (i.e. uu-uu,dd-dd,..,UU-UU,..,CC-CC)
C
C**** JMPROC(14) to JMPROC(61) = qq'-qq', qQ'-qQ', Qq'-Qq', QQ'-QQ'
C
C**** JMPROC(62) to JMPROC(93) = qQ-(uU, dD, sS, cC), Qq-(uU, dD, sS, cC)
C
C**** JMPROC(94) to JMPROC(101)  = qQ-gg, Qq-gg
C
C**** JMPROC(102) to JMPROC(117) = gq-gq, qg-gq, gQ-gQ, Qg-gQ
C
C**** THE SUM OF THE JMPROC(i)'s RESIDES IN JMPTOT
C********************************************************************
C       gg-gg
      JMPROC(1)=9./2.*(3.0D0-t*u/s/s-s*u/t/t-s*t/u/u)*g1*g2
C       gg-qQ
      JMPROC(2)=(t/u+u/t-9./4.*(t*t+u*u)/s/s)/3.d0*g1*g2
      JMPROC(3)=JMPROC(2)
      JMPROC(4)=JMPROC(2)
      JMPROC(5)=JMPROC(2)
C       qq-qq,QQ-QQ
        f=4./9.*((s*s+u*u)/t/t+(s*s+t*t)/u/u)-8./27.*s*s/t/u
      JMPROC(6)=u1*u2*f
      JMPROC(7)=d1*d2*f
      JMPROC(8)=s1*s2*f
      JMPROC(9)=c1*c2*f
      JMPROC(10)=uu1*uu2*f
      JMPROC(11)=dd1*dd2*f
      JMPROC(12)=ss1*ss2*f
      JMPROC(13)=cc1*cc2*f
C       qq'-qq', qQ'-qQ', Qq'-Qq', QQ'-QQ'
        f=4./9.*((s*s+u*u)/t/t+(s*s+t*t)/u/u)
      JMPROC(14)=u1*d2*f
      JMPROC(15)=u1*s2*f
      JMPROC(16)=u1*c2*f
      JMPROC(17)=d1*u2*f
      JMPROC(18)=d1*s2*f
      JMPROC(19)=d1*c2*f
      JMPROC(20)=s1*u2*f
      JMPROC(21)=s1*d2*f
      JMPROC(22)=s1*c2*f
      JMPROC(23)=c1*u2*f
      JMPROC(24)=c1*d2*f
      JMPROC(25)=c1*s2*f
      JMPROC(26)=u1*dd2*f
      JMPROC(27)=u1*ss2*f
      JMPROC(28)=u1*cc2*f
      JMPROC(29)=d1*uu2*f
      JMPROC(30)=d1*ss2*f
      JMPROC(31)=d1*cc2*f
      JMPROC(32)=s1*uu2*f
      JMPROC(33)=s1*dd2*f
      JMPROC(34)=s1*cc2*f
      JMPROC(35)=c1*uu2*f
      JMPROC(36)=c1*dd2*f
      JMPROC(37)=c1*ss2*f
      JMPROC(38)=uu1*d2*f
      JMPROC(39)=uu1*s2*f
      JMPROC(40)=uu1*c2*f
      JMPROC(41)=dd1*u2*f
      JMPROC(42)=dd1*s2*f
      JMPROC(43)=dd1*c2*f
      JMPROC(44)=ss1*u2*f
      JMPROC(45)=ss1*d2*f
      JMPROC(46)=ss1*c2*f
      JMPROC(47)=cc1*u2*f
      JMPROC(48)=cc1*d2*f
      JMPROC(49)=cc1*s2*f
      JMPROC(50)=uu1*dd2*f
      JMPROC(51)=uu1*ss2*f
      JMPROC(52)=uu1*cc2*f
      JMPROC(53)=dd1*uu2*f
      JMPROC(54)=dd1*ss2*f
      JMPROC(55)=dd1*cc2*f
      JMPROC(56)=ss1*uu2*f
      JMPROC(57)=ss1*dd2*f
      JMPROC(58)=ss1*cc2*f
      JMPROC(59)=cc1*uu2*f
      JMPROC(60)=cc1*dd2*f
      JMPROC(61)=cc1*ss2*f
C       qQ-(uU, dD, sS, cC), Qq-(uU, dD, sS, cC)
        f1=4./9.*(u*u+t*t)/s/s*2.d0
        f2=f1+4./9.*((s*s+u*u)/t/t+(s*s+t*t)/u/u)-8./27.*
     *     (u*u/s/t+t*t/s/u)
      JMPROC(62)=u1*uu2*f2
      JMPROC(63)=u1*uu2*f1
      JMPROC(64)=JMPROC(63)
      JMPROC(65)=JMPROC(63)
      JMPROC(66)=d1*dd2*f1
      JMPROC(67)=d1*dd2*f2
      JMPROC(68)=JMPROC(66)
      JMPROC(69)=JMPROC(66)
      JMPROC(70)=s1*ss2*f1
      JMPROC(71)=JMPROC(70)
      JMPROC(72)=s1*ss2*f2
      JMPROC(73)=JMPROC(70)
      JMPROC(74)=c1*cc2*f1
      JMPROC(75)=JMPROC(74)
      JMPROC(76)=JMPROC(74)
      JMPROC(77)=c1*cc2*f2
      JMPROC(78)=uu1*u2*f2
      JMPROC(79)=uu1*u2*f1
      JMPROC(80)=JMPROC(79)
      JMPROC(81)=JMPROC(79)
      JMPROC(82)=dd1*d2*f1
      JMPROC(83)=dd1*d2*f2
      JMPROC(84)=JMPROC(82)
      JMPROC(85)=JMPROC(82)
      JMPROC(86)=ss1*s2*f1
      JMPROC(87)=JMPROC(86)
      JMPROC(88)=ss1*s2*f2
      JMPROC(89)=JMPROC(86)
      JMPROC(90)=cc1*c2*f1
      JMPROC(91)=JMPROC(90)
      JMPROC(92)=JMPROC(90)
      JMPROC(93)=cc1*c2*f2
C       qQ-gg, Qq-gg
        f=32./27.*(t/u+u/t)-8./3.*(t*t+u*u)/s/s
      JMPROC(94)=u1*uu2*f
      JMPROC(95)=d1*dd2*f
      JMPROC(96)=s1*ss2*f
      JMPROC(97)=c1*cc2*f
      JMPROC(98)=uu1*u2*f
      JMPROC(99)=dd1*d2*f
      JMPROC(100)=ss1*s2*f
      JMPROC(101)=cc1*c2*f
C       gq-gq, qg-gq, gQ-gQ, Qg-gQ
        f=(s*s+u*u)/t/t+(s*s+t*t)/u/u-4./9.*(s/u+u/s+s/t+t/s)
      JMPROC(102)=g1*u2*f
      JMPROC(103)=g1*d2*f
      JMPROC(104)=g1*s2*f
      JMPROC(105)=g1*c2*f
      JMPROC(106)=g2*u1*f
      JMPROC(107)=g2*d1*f
      JMPROC(108)=g2*s1*f
      JMPROC(109)=g2*c1*f
      JMPROC(110)=g1*uu2*f
      JMPROC(111)=g1*dd2*f
      JMPROC(112)=g1*ss2*f
      JMPROC(113)=g1*cc2*f
      JMPROC(114)=g2*uu1*f
      JMPROC(115)=g2*dd1*f
      JMPROC(116)=g2*ss1*f
      JMPROC(117)=g2*cc1*f

      JMPTOT=0.D0
      DO I=1,117
C       Parton type from photon.
*        PID(1) = JMPTYP(I)/1000000
C       Parton type from proton.
*        PID(2) = MOD(JMPTYP(I)/10000,100)
C       Scattered partons
        PID(3) = MOD(JMPTYP(I)/100,100)
        PID(4) = MOD(JMPTYP(I),100)
        DO J=3,4
          IF (PID(J).EQ.0) THEN
            PID(J)=13
          ELSE IF (PID(J).GT.10) THEN
            PID(J) = PID(J)-4
          ENDIF
        ENDDO
        SMIN = RMASS(PID(3))**2 + RMASS(PID(4))**2+ TWO*PTJIM**2+TWO
     $       * SQRT(PTJIM**2+RMASS(PID(3))**2)*SQRT(PTJIM**2+RMASS(PID(4
     $       ))**2)
        IF (S.LT.SMIN) THEN
          JMPROC(I)=0.0
        ELSE
          JMPTOT=JMPTOT+JMPROC(I)
        ENDIF
      ENDDO

      IF (FLAG2.EQ.0) THEN
C       Calculate dxsec/d(pt2)
        JMXS1=JMPTOT*PIFAC*ALPHA**2/(S**2 * SQRT(1.D0-Y))
      ELSE
        WRITE(JMOUT,8900) 'JMXS1 Called with illegal flag2'
        WRITE(JMOUT,*) flag2
        JMXS1=0.D0
      ENDIF

C     WEISZACKER-WILLIAMS flux factor.
C     Only Want This When NOT integrating (i.e. never again)
C     (flag=0 otherwise)
      IF (FLAG1.NE.0) THEN
        WRITE(JMOUT,8900) 'JMXS1 ERROR'
      ENDIF

 8900 FORMAT(A)

      RETURN
      END
