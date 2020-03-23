c         1         2         3         4         5         6         7         8
c---------|---------|---------|---------|---------|---------|---------|---------|
*-- Author :    R.Lednicky      20/01/95
      SUBROUTINE FSIW(J,WEIF,WEI,WEIN)

C=======================================================================
C    Calculates final state interaction (FSI) weights
C    WEIF = weight due to particle - (effective) nucleus FSI (p-N)
C    WEI  = weight due to p-p-N FSI
C    WEIN = weight due to p-p FSI; note that WEIN=WEI if I3C=0;
C                                  note that if I3C=1 the calculation of
C                                  WEIN can be skipped by putting J=0
C.......................................................................
C    Correlation Functions:
C    CF(p-p-N)   = sum(WEI)/sum(WEIF)
C    CF(p-p)     = sum(WEIN)/sum(1); here the nucleus is completely
C                                    inactive
C    CF(p-p-"N") = sum(WEIN*WEIF')/sum(WEIF'), where WEIN and WEIF'
C                  are not correlated (calculated at different emission
C                  points, e.g., for different events);
C                  thus here the nucleus affects one-particle
C                  spectra but not the correlation
C.......................................................................
C    User must supply data file <fn> on unit NUNIT (e.g. =11) specifying
C    LL   : particle pair
C    NS   : approximation used to calculate Bethe-Salpeter amplitude
C    ITEST: test switch
C           If ITEST=1 then also following parameters are required
C    ICH  : 1(0) Coulomb interaction between the two particles ON (OFF)
C    IQS  : 1(0) quantum statistics for the two particles ON (OFF)
C    ISI  : 1(0) strong interaction between the two particles ON (OFF)
C    I3C  : 1(0) Coulomb interaction with residual nucleus ON (OFF)
C    This data file can contain other information useful for the user.
C    It is read by subroutines READINT4 and READREA8(4) (or READ_FILE).
C----------------------------------------------------------------------
C-   LL       1  2  3  4  5   6   7   8  9 10  11  12  13  14 15 16 17
C-   part. 1: n  p  n  a  pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-
C-   part. 2: n  p  p  a  pi- pi0 pi+ d  d  K-  K+  p   p  K- K+ p  p
C   NS=1 y/n: +  +  +  +  +   -   -   -  -  -   -   -   -  -  -  -  -
C----------------------------------------------------------------------
C-   LL       18 19 20 21 22 23  24 25 26 27 28 29 30
C-   part. 1: d  d  t  t  K0 K0  d  p  p  p  n  lam p
C-   part. 2: d  a  t  a  K0 K0b t  t  a  lam lam lam pb
C   NS=1 y/n: -  -  -  -  -  -   -  -  -  +  +  +  -
C----------------------------------------------------------------------
C   NS=1  Square well potential,
C   NS=3  not used
C   NS=4  scattered wave approximated by the spherical wave,
C   NS=2  same as NS=4 but the approx. of equal emission times in PRF
C         not required (t=0 approx. used in all other cases).
C   Note: if NS=2,4, the B-S amplitude diverges at zero distance r* in
C         the two-particle c.m.s.; user can specify a cutoff AA in
C         SUBROUTINE FSIINI, for example:
C         IF(NS.EQ.2.OR.NS.EQ.4)AA=5.D0 !! in 1/GeV --> AA=1. fm
C---------------------------------------------------------------------
C    ITEST=1 any values of parameters ICH, IQS, ISI, I3C are allowed
C            and should be given in data file <fn>
C    ITEST=0 physical values of these parameters are put automatically
C            in FSIINI (their values are not required in data file)
C=====================================================================
C    At the beginning of calculation user should call FSIINI,
C    which reads LL, NS, ITEST (and eventually ICH, IQS, ISI, I3C)
C    and initializes various parameters.
C    In particular the constants in
C      COMMON/FSI_CONS/PI,PI2,SPI,DR,W
C    may be useful for the user:
C     W=1/.1973D0    ! from fm to 1/GeV
C     PI=4*DATAN(1.D0)
C     PI2=2*PI
C     SPI=DSQRT(PI)
C     DR=180.D0/PI   ! from radian to degree
C      _______________________________________________________
C  !! |Important note: all real quantities are assumed REAL*8 | !!
C      -------------------------------------------------------
C    For each event user should fill in the following information
C    in COMMONs (all COMMONs in FSI calculation start with FSI_):
C    ...................................................................
C     COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
C    Only
C         AMN  = mass of the effective nucleus   [GeV/c**2]
C         CN   = charge of the effective nucleus [elem. charge units]
C    are required
C    ...................................................................
C     COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1, !part. momenta in the rest frame
C    1               P2X,P2Y,P2Z,E2,P2  !of effective nucleus (NRF)
C    Only the components
C                        PiX,PiY,PiZ  [GeV/c]
C    in NRF are required.
C    To make the corresponding Lorentz transformation user can use the
C    subroutines LTRAN and LTRANB
C    ...................................................................
C    COMMON/FSI_COOR/X1,Y1,Z1,T1,R1,     ! 4-coord. of emission
C    1               X2,Y2,Z2,T2,R2      ! points in NRF
C    The componets
C                       Xi,Yi,Zi  [fm]
C    and emission times
C                          Ti   [fm/c]
C    should be given in NRF with the origin assumed at the center
C    of the effective nucleus. If the effect of residual nucleus is
C    not calculated within FSIW, the NRF can be any fixed frame.
C-----------------------------------------------------------------------
C    Before calling FSIW the user must call
C     CALL LTRAN12
C    Besides Lorentz transformation to pair rest frame:
C    (p1-p2)/2 --> k* it also transforms 4-coordinates of
C    emission points from fm to 1/GeV and calculates Ei,Pi and Ri.
C    Note that |k*|=AK in COMMON/FSI_PRF/
C-----------------------------------------------------------------------
C    After making some additional filtering using k* (say k* < k*max)
C    or direction of vector k*,
C    user can finally call FSIW to calculate the FSI weights
C    to be used to construct the correlation function
C=======================================================================

      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_JR/JRAT
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  ! particle momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS, ! k*=(p1-p2)/2 and x1-x2
     1               X,Y,Z,T,RP,RPS      ! in pair rest frame (PRF)
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, !4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_FFPN/FF12,FF21
      COMPLEX*16 FF12,FF21
C------------------------------------------------------------------
C==> AC1,2 = "relativistic" Bohr radii for particle-nucleus systems
      C1N=C1*CN
      IF(C1N.NE.0.D0)AC1=137.036D0/(C1N*E1) !m1-->E1
      C2N=C2*CN
      IF(C2N.NE.0.D0)AC2=137.036D0/(C2N*E2) !m2-->E2

C-----------------------------------------------------------
      CALL FSIPN(WEIF) !weight due to particle-nucleus FSI
      JRAT=0
      call BoostToPrf()
      CALL FSIWF(WEI)  !weight due to particle-particle-nucleus FSI
      WEIN=WEI
      IF(I3C*J.NE.0) THEN
         FF12=DCMPLX(1.D0,0.D0)
         FF21=DCMPLX(1.D0,0.D0)
         JRAT=1
         CALL VZ(WEIN)          ! weight due to particle-particle FSI
         ENDIF
      RETURN
      END
c_________________________________________________________________
      subroutine FsiWFast(WEI,aDX,aDY,aDZ,aDPx,aDPY,aDPZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS, ! k*=(p1-p2)/2 and x1-x2
     1               X,Y,Z,T,RP,RPS      ! in pair rest frame (PRF)

      X=aDx
      Y=aDy
      Z=aDz
      T=0
      RPS=X*X+Y*Y+Z*Z
      RP=DSQRT(RPS)
      PPX=aDPx
      PPY=aDPy
      PPZ=aDPz
      AKS=PPX*PPX+PPY*PPY+PPZ*PPZ
      AK=DSQRT(AKS)
      CALL FSIWF(WEI)  !weight due to particle-particle-nucleus FSI
      RETURN
      END
C=======================================================================


      SUBROUTINE LTRAN(P0,P,PS)
C==>calculating particle 4-momentum PS={PSX,PSY,PSZ,ES}
C   in rest frame of a system 0 with 4-momentum P0={P0X,P0Y,P0Z,E0}
C   from its 4-momentum P={PX,PY,PZ,E}

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P0(4),P(4),PS(4)
C-----------------------------------------------------------------------
      P0S=P0(1)**2+P0(2)**2+P0(3)**2
      AM0=DSQRT(P0(4)**2-P0S)
      EPM=P0(4)+AM0
      PP0=P(1)*P0(1)+P(2)*P0(2)+P(3)*P0(3)
      H=(PP0/EPM-P(4))/AM0
      PS(1)=P(1)+P0(1)*H
      PS(2)=P(2)+P0(2)*H
      PS(3)=P(3)+P0(3)*H
      PS(4)=(P0(4)*P(4)-PP0)/AM0
      RETURN
      END

      SUBROUTINE LTRANB(P0,PS,P)
C==>calculating particle 4-momentum P={PX,PY,PZ,E}
C   from its 4-momentum PS={PSX,PSY,PSZ,ES}
C   in rest frame of a system 0 with 4-momentum P0={P0X,P0Y,P0Z,E0}

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P0(4),P(4),PS(4)
C-----------------------------------------------------------------------
      P0S=P0(1)**2+P0(2)**2+P0(3)**2
      AM0=DSQRT(P0(4)**2-P0S)
      EPM=P0(4)+AM0
      PSP0=PS(1)*P0(1)+PS(2)*P0(2)+PS(3)*P0(3)
      HS=(PSP0/EPM+PS(4))/AM0
      P(1)=PS(1)+P0(1)*HS
      P(2)=PS(2)+P0(2)*HS
      P(3)=PS(3)+P0(3)*HS
      P(4)=(P0(4)*PS(4)+PSP0)/AM0
      RETURN
      END

      SUBROUTINE LTRAN12
C==>calculating particle momentum in PRF {EE,PPX,PPY,PPZ} from
C-  the momentum of the first particle {E1,P1X,P1Y,P1Z) in NRF
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  !part. momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_P12/P12X,P12Y,P12Z,E12,P12,AM12,EPM
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, !4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_CONS/PI,PI2,SPI,DR,W

C   fm --> 1/GeV
      X1=X1*W
      Y1=Y1*W
      Z1=Z1*W
      T1=T1*W
      X2=X2*W
      Y2=Y2*W
      Z2=Z2*W
      T2=T2*W
C   calculating Ri, Pi and Ei
      R1=DSQRT(X1*X1+Y1*Y1+Z1*Z1)
      R2=DSQRT(X2*X2+Y2*Y2+Z2*Z2)
      P1S=P1X*P1X+P1Y*P1Y+P1Z*P1Z
      P2S=P2X*P2X+P2Y*P2Y+P2Z*P2Z
      P1=DSQRT(P1S)
      P2=DSQRT(P2S)
      E1=DSQRT(AM1*AM1+P1S)
      E2=DSQRT(AM2*AM2+P2S)
C-----------------------------------------------------------------------
      E12=E1+E2
      P12X=P1X+P2X
      P12Y=P1Y+P2Y
      P12Z=P1Z+P2Z
      P12S=P12X**2+P12Y**2+P12Z**2
      AM12=DSQRT(E12**2-P12S)
      EPM=E12+AM12
      P12=DSQRT(P12S)
      P112=P1X*P12X+P1Y*P12Y+P1Z*P12Z
      H1=(P112/EPM-E1)/AM12
      PPX=P1X+P12X*H1
      PPY=P1Y+P12Y*H1
      PPZ=P1Z+P12Z*H1
      EE=(E12*E1-P112)/AM12
      AKS=EE**2-AM1**2
      AK=DSQRT(AKS)

CW      WRITE(6,38)'AK ',AK,'K ',PPX,PPY,PPZ,EE
38    FORMAT(A7,E11.4,A7,4E11.4)
      RETURN
      END

      SUBROUTINE FSIPN(WEIF)
C  calculating particle-nucleus Coulomb Wave functions FFij
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  !part. momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, ! 4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_ICH1/ICH1
      COMMON/FSI_ETA/ETA
      COMMON/FSI_FFPN/FF12,FF21
      COMPLEX*16 FF1,FF12,FF21
      FF12=DCMPLX(1.D0,0.D0)
      FF21=DCMPLX(1.D0,0.D0)
      ACH=1.D0
      WEIF=1.D0
      IF(I3C.EQ.0)RETURN
      ICH1=IDINT(C1)
      IF(ICH1.EQ.0)GOTO 11
      XH=AC1*P1
      ACH=ACP(XH)
      ACHR=DSQRT(ACH)
      ETA=0.D0
      IF(XH.NE.0.D0)ETA=1/XH
      RHOS=P1*R1
      HS=X1*P1X+Y1*P1Y+Z1*P1Z
      FF12=FF12*FF1(RHOS,HS)
      IF(IQS.EQ.0)GOTO 11
      RHOS=P1*R2
      HS=X2*P1X+Y2*P1Y+Z2*P1Z
      FF21=FF21*FF1(RHOS,HS)
 11   ICH1=IDINT(C2)
      IF(ICH1.EQ.0)GOTO 10
      XH=AC2*P2
      ACH=ACP(XH)
      ACHR=DSQRT(ACH)
      ETA=0.D0
      IF(XH.NE.0.D0)ETA=1/XH
      RHOS=P2*R2
      HS=X2*P2X+Y2*P2Y+Z2*P2Z
      FF12=FF12*FF1(RHOS,HS)
CW      WRITE(6,41)'AC2 ',AC2,'ACH ',ACH,'ETA ',ETA,'RHOS ',RHOS,'HS ',HS
41    FORMAT(5(A5,E11.4))
CW      WRITE(6,40)'FF12 ',DREAL(FF12),DIMAG(FF12)
      IF(IQS.EQ.0)GOTO 10
      RHOS=P2*R1
      HS=X1*P2X+Y1*P2Y+Z1*P2Z
      FF21=FF21*FF1(RHOS,HS)
CW      WRITE(6,41)'AC1 ',AC1,'ACH ',ACH,'ETA ',ETA,'RHOS ',RHOS,'HS ',HS
CW      WRITE(6,40)'FF21 ',DREAL(FF21),DIMAG(FF21)
40    FORMAT(A7,2E12.4)
 10   CONTINUE

C  WEIF = the weight due to the Coulomb particle-nucleus interaction
      WEIF=DREAL(FF12)**2+DIMAG(FF12)**2
      IF(IQS.EQ.1)WEIF=0.5D0*(WEIF+DREAL(FF21)**2+DIMAG(FF21)**2)
      RETURN
      END

      FUNCTION GPIPI(X,J)
C--- GPIPI = k*COTG(DELTA), X=k^2
C--  J=1(2) corresponds to isospin=0(2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AAPI/AAPI(20,2)
      COMMON/FSI_C/HELP(20),AM,AMS,DM
      OMS=X+AMS
      OM=DSQRT(OMS)
      XX=X/AMS
      XXS=XX*XX
c      GOTO 1 ! Nagels (79)
C--- Colangelo(01)
      GPIPI=OM*(4*OMS-AAPI(5,J))/(4*AMS-AAPI(5,J))
      GPIPI=GPIPI/(AAPI(1,J)+AAPI(2,J)*XX+AAPI(3,J)*XXS+
     +             AAPI(4,J)*XXS*XX)
      RETURN
C--- Nagels (79)
 1    continue
      GPIPI=OM/AAPI(1,J)
      GPIPI=GPIPI*(1+(AAPI(3,J)-AAPI(1,J)**2)*XX+AAPI(4,J)*XXS)
      GPIPI=GPIPI/(1+(AAPI(3,J)+AAPI(2,J)/AAPI(1,J))*XX)
      RETURN
      END

      FUNCTION GPIN(X,J)
C--- GPIN = k*COTG(DELTA), X=k^2
C--  J=1(2) corresponds to piN isospin=1/2(3/2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AAPIN/AAPIN(20,2)
      GPIN=1/AAPIN(1,J)+.5D0*AAPIN(2,J)*X
      RETURN
      END

      FUNCTION GND(X,J)
C--- GND = k*COTG(DELTA), X=k^2
C--- J=1(2) corresp. to nd(pd), S=1/2,
C--- J=3(4) corresp. to nd(pd), S=3/2
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AAND/AAND(20,4)
      XX=X
      GND=1/AAND(1,J)+.5D0*AAND(2,J)*X
      DO 1 I=4,4
      XX=XX*X
   1  GND=GND+AAND(I,J)*XX
      GND=GND/(1+AAND(3,J)*X)
      RETURN
      END
      FUNCTION GDD(X,J)
C--- GDD = k*COTG(DELTA), X=k^2
C--- J=1,2,3 corresp. to S=0,1,2
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AADD/AADD(20,3)
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMMON/FSI_CONS/PI,PI2,SPI,DR,W
      COMPLEX*16 C
      E=X/2/AM
      ER=DSQRT(E)
      IF(J.EQ.1)THEN
       GDD=ER*(AADD(1,1)*DEXP(-E/AADD(2,1))-AADD(3,1))
       GDD=GDD/DR   ! from degree to radian
       TAND=DTAN(GDD)
       IF(TAND.EQ.0.D0)TAND=1.D-10
       GDD=DSQRT(X)/TAND
      END IF
      IF(J.EQ.2)THEN
       GDD=1.D10
      END IF
      IF(J.EQ.3)THEN
       GDD=ER*(AADD(1,3)+AADD(2,3)*E)
       GDD=GDD/DR    ! from degree to radian
       TAND=DTAN(GDD)
       IF(TAND.EQ.0.D0)TAND=1.D-10
       GDD=DSQRT(X)/TAND
      END IF
      RETURN
      END
      BLOCK DATA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AAPI/AAPI(20,2)/FSI_AAND/AAND(20,4)
      COMMON/FSI_AADD/AADD(20,3)/FSI_AAPIN/AAPIN(20,2)
      COMMON/FSI_AAKK/AAKK(9)/FSI_AAPAP/AAPAPR(3,2),AAPAPI(3,2)
C---Parameters for GPIPI (I,J), J=1,2 -> isospin=0,2
C--- Nagels (79)
c      DATA AAPI/.2600D00, .2500D00, .3480D00,-.0637D00, 16*0.D0,
c     1         -.0280D00,-.0820D00, .2795D00,-.0086D00, 16*0.D0/
C--- Colangelo(01)
      DATA AAPI/.220D0, .268D0,-.139D-1,-.139D-2, 36.77D0, 15*0.D0,
     1         -.444D-1,-.857D-1,-.221D-2,-.129D-3,-21.62D0, 15*0.D0/
C***
C---Parameters for GPIN (I,J), J=1,2 -> piN isospin=1/2,3/2
      DATA AAPIN/ .12265D1, .1563D2,     18*0.D0,
     1            -.750D0,  .7688D2,     18*0.D0/
C---Parameters for GND (I,J), J=1(2) corresp. to nd(pd), S=1/2,
C---                          J=3(4) corresp. to nd(pd), S=3/2
      DATA AAND/-.3295D1,-.8819D3, .4537D4, .28645D5, 16*0.D0,
     1          -.13837D2, .11505D2, .0D0, .10416D2,  16*0.D0,
     2          -.32180D2, .1014D2,  .0D0, .0D0,      16*0.D0,
     3          -.60213D2, .1333D2,  .0D0,-.70309D2,  16*0.D0/
      DATA AADD/ .10617D4, .3194D-2, .56849D3, 17*0.D0,
     1           20*0.D0,
     2          -.1085D4, .1987D5, 18*0.D0/
C--- AAKK= m_K^2, m_pi^2, m_eta^2,
C---       m_delta^2, gam(S*-->KK-b),
C---       gam(S*-->pipi), gam(delta-->KK-b), gam(delta-->pi eta)
      DATA AAKK/.247677D0,.01947977D0,.2997015D0,
     1 .9565D0,.9487D0,.792D0,.199D0,.333D00,.222D0/ ! Martin'77
c     1 .9467D0,.9698D0,2.763D0,.5283D0,.4038D0,.3711D0/ ! Antonelli'02
c     1 .9920D0,.9841D0,1.305D0,.2684D0,.5555D0,.4401D0/ ! Achasov1'01,03
c     1 .9920D0,1.0060D0,1.305D0,.2684D0,.8365D0,.4580D0/! Achasov2'01,03
C---Parameters for PAP (I,J), j=1,2 -> isospin I=0,2
C---                          i=1-3 -> a_singlet, a_triplet, d [fm]
C---    Im a_IS (I=isospin, S=spin) are fixed by atomic data and
C       n-bar survival time up to one free parameter, e.g. Im a_00
C---    Batty (89), Kerbikov (93):
C--- Ima_10=1.96-Ima_00, Ima_01=0.367-Ima_00/3, Ima_11=0.453+Ima_00/3
C---       In DATA we put Ima_00=0.3.
C---    Re a_IS are fixed by atomic data up to three free parameters
C---    Batty (89):
C---         Rea_aver(pp-bar)=Re[(a_00+a_10)+3(a_01+a_11)]/8=-0.9
C---       In DATA we used Rea_IS from Paris potential Pignone (94)
C---       rescaled by 1.67 to satisfy the atomic constraint.
C---    Effective radius is is taken independent of IS from the phase
C---    shift fit by Pirner et al. (91).
      DATA AAPAPR/-0.94D0, -1.98D0,  .1D0,
     1            -1.40D0,  0.37D0,  .1D0/ ! Re
      DATA AAPAPI/ 0.3 D0,  .267D0,-.01D0,
     1             1.66D0,  .553D0,-.01D0/ ! Im
      END

      SUBROUTINE CPIPI ! calculates pi+pi- scattering amplitude
                       ! accounting for pi0pi0->pi+pi- transition
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMMON/FSI_2CHA/AK2,AK2S,AMU2_AMU1 ! k* (kappa) for 2-nd channel
      COMPLEX*16 C
      DATA AM2C/.13498D0/ ! pi0 mass
      AMU2_AMU1=AM2C/AM   ! AM2C=2*mu(pi0pi0), AM=2*mu(pi+pi-)
      AK2S0=AMU2_AMU1*AKS
      AK2S =AK2S0-2*AM2C*(AM2C-AM)
      AK2=DSQRT(AK2S)  !  k2
      GPI2=GPIPI(AKS,2)
      GPI1=GPIPI(AKS,1)
      RK11=.6667D0*GPI1+.3333D0*GPI2
      RK22=.6667D0*GPI2+.3333D0*GPI1
      RK12=-0.47140452D0*(GPI1-GPI2)
      C(3)=RK11+DCMPLX(-HCP2,-AAK)      ! (1/f)11
      C(5)=RK22
      C(5)=C(5)-DCMPLX(0.D0,AK2)        ! (1/f)22
      C(7)=RK12                         ! (1/f)12
      C(10)=C(3)*C(5)-C(7)*C(7)
      C(1)=C(5)/C(10)                   ! f11
      C(2)=-C(7)/C(10)                  ! f12
      RETURN
      END

      SUBROUTINE CKKB  ! calculates KK-b scattering amplitude,
                       ! saturated by S*(980) and delta(982) resonances
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_AAKK/AAKK(9)
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMPLEX*16 C
      S4=AKS+AAKK(1)
      S=4*S4
      AKPIPI=DSQRT(S4-AAKK(2))
      EETA2=(S+AAKK(3)-AAKK(2))**2/4/S
      AKPIETA=DSQRT(EETA2-AAKK(3))
      C(1)=AAKK(6)/2/DCMPLX(AAKK(4)-S,
     ,-AK*AAKK(6)-AKPIPI*AAKK(7))
      C(1)=C(1)+AAKK(8)/2/DCMPLX(AAKK(5)-S,
     ,-AK*AAKK(8)-AKPIETA*AAKK(9))
      RETURN
      END

      SUBROUTINE CPAP  ! calculates pp-bar scattering amplitude
                       ! accounting for nn-bar->pp-bar channel
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1     X,Y,Z,T,RP,RPS
      COMMON/FSI_AAPAP/AAPAPR(3,2),AAPAPI(3,2)
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMMON/FSI_2CHA/AK2,AK2S,AMU2_AMU1 ! k* (kappa) for 2-nd channel
      COMPLEX*16 C
        DATA AM2C/.93956563D0/ ! neutron mass
      AMU2_AMU1=AM2C/AM   ! AM2C=2*mu(nn-bar), AM=2*mu(pp-bar)
      AK2S0=AMU2_AMU1*AKS
        AK2S =AK2S0-2*AM2C*(AM2C-AM)
      IF(AK2S.GE.0.D0)THEN
           AK2=DCMPLX(DSQRT(AK2S),0.D0)  !  k2
       ELSE
           AK2=DCMPLX(0.D0,DSQRT(-AK2S)) !  kappa2
      ENDIF
      C(10)=C(6+(ISPIN-1)*2)+
     + DCMPLX(AAPAPR(3,ISPIN)*AKS/2-0.016D0-HCP2,
     ,        AAPAPI(3,ISPIN)*AKS/2-AAK)         ! (1/f)11
      C(5)=C(6+(ISPIN-1)*2)+
     + DCMPLX(AAPAPR(3,ISPIN)*AK2S0/2,
     ,        AAPAPI(3,ISPIN)*AK2S0/2)
      IF(AK2S.GE.0.D0)THEN
           C(5)=C(5)-DCMPLX(0.D0,AK2)
       ELSE
           C(5)=C(5)+DCMPLX(AK2,0.D0)              ! (1/f)22
        ENDIF
      C(10)=C(10)*C(5)-C(7+(ISPIN-1)*2)*C(7+(ISPIN-1)*2)
      C(ISPIN)=C(5)/C(10)                ! f11
      C(ISPIN+2)=-C(7+(ISPIN-1)*2)/C(10) ! f12
      RETURN
      END

c      SUBROUTINE FSIINI(I_ITEST,I_ICH,I_IQS,I_ISI,I_I3C)
      SUBROUTINE FSIINI(I_ITEST,I_LL,I_NS,I_ICH,I_IQS,I_ISI,I_I3C)
C          SUBROUTINE FSIINI
C---Note:
C-- ICH= 0 (1) if the Coulomb interaction is absent (present);
C-- ISPIN= JJ= 1,2,..,MSPIN denote increasing values of the pair
C-- total spin S.
C-- To calculate the CF of two particles (with masses m1, m2 and
C-- charges C1, C2) the following information is required:
C-- AM= twice the reduced mass= 2*m1*m2/(m1+m2) in GeV/c^2,
C-- DM= (m1-m2)/(m1+m2), required if NS=2;
C-- AC= Bohr radius= 2*137.036*0.1973/(C1*C2*AMH) in fm;
C-- AC > 1.D9 if C1*C2= 0, AC < 0 if C1*C2 < 0;
C-- MSPIN= MSPINH(LL)= number of the values of the total pair spin S;
C-- FD= FDH(LL,JJ), RD= RDH(LL,JJ)= scattering length and effective
C-- radius for each value of the total pair spin S, JJ= 1,..,MSPIN;     ;
C-- the corresponding square well parameters EB= EBH(LL,JJ), RB=
C-- RBH(LL,JJ) (required if NS=1) may be calculated by sear.f;
C-- if the effective range approximation is not valid (as is the case,
C-- e.g., for two-pion system) a code for calculation of the scattering
C-- amplitude should be supplemented;
C-- RHO= RHOH(LL,JJ), SF= SFH(LL,JJ), SE= SEH(LL) are spin factors;
C-- RHO= the probability that the spins j1 and j2 of the two particles
C-- will combine in a total spin S;
C-- RHO= (2*S+1)/[(2j1+1)*(2j2+1)] for unpolarized particles;
C-- RHO= (1-P1*P2)/4 and (3+P1*P2)/4 correspond to S=0 and 1 in the
C-- case of spin-1/2 particles with polarizations P1 and P2;
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_SPIN/RHO(10)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_FD/FD(10),RD(10)
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMMON/FSI_CONS/PI,PI2,SPI,DR,W
      COMPLEX*16 C
      COMMON/FSI_AA/AA
      COMMON/FSI_AAPI/AAPI(20,2)/FSI_AAND/AAND(20,4)
      COMMON/FSI_AAPIN/AAPIN(20,2)
      COMMON/FSI_AAPAP/AAPAPR(3,2),AAPAPI(3,2)
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
      COMMON/FSI_FDH/FDH(30,10),RDH(30,10),EBH(30,10),RBH(30,10)
      COMMON/FSI_RHOH/RHOH(30,10)
      COMMON/FSI_AMCH/AM1H(30),AM2H(30),C1H(30),C2H(30),MSPINH(30)
C============= declarations pour l'appel de READ_FILE()============
c      CHARACTER*10 KEY
c      CHARACTER*8  CH8
c      INTEGER*4    INT4
c      REAL*8       REAL8
c      INTEGER*4    IERR
C
C--- mass of the first and second particle
      DATA AM1H/.93956563D0,.93827231D0,.93956563D0,3.72737978D0,
     C          .13957D0,.13498D0,.13957D0, .93956563D0, .93827231D0,
     C          4*.13957D0,4*.493677D0,
     C          2*1.87561339D0,2*2.80892165D0,2*.497672D0,
     C          1.87561339D0,3*.93827231D0,.93956563D0,
     C          1.115684D0,.93827231D0/
      DATA AM2H/.93956563D0,.93827231D0,.93827231D0,3.72737978D0,
     C          .13957D0,.13498D0,.13957D0, 2*1.87561339D0,
     C          2*.493677D0,2*.93827231D0,
     C          2*.493677D0,2*.93827231D0,
     C          1.87561339D0,3.72737978D0,2.80892165D0,3.72737978D0,
     C          2*.497672D0,2*2.80892165D0,3.72737978D0,
     C          3*1.115684D0,.93827231D0/
c--------|---------|---------|---------|---------|---------|---------|----------
C---  charge of the first and second particle
      DATA C1H/0.D0,1.D0,0.D0,2.D0, 1.D0,0.D0,1.D0,0.D0,1.D0,
     C         3*1.D0,-1.D0,3*1.D0,-1.D0,
     C         4*1.D0,2*0.D0,4*1.D0,2*0.D0, 1.D0/
      DATA C2H/0.D0,1.D0,1.D0,2.D0,-1.D0,0.D0,3*1.D0,
     C         -1.D0,3*1.D0,-1.D0,3*1.D0,
     C         1.D0,2.D0,1.D0,2.D0,2*0.D0,2*1.D0,2.D0,3*0.D0,-1.D0/
C---MSPIN vs (LL)
      DATA MSPINH/3*2,4*1,2*2,8*1,3,1,2,1,2*1,2*2,1,3*2, 2/
C---Spin factors <RHO vs (LL,ISPIN)
      DATA RHOH/3*.25D0, 4*1.D0, 2*.3333D0, 8*1.D0,
     1          .1111D0,1.D0,.25D0,1.D0,2*1.D0,
     1          .3333D0,.25D0,1.D0,3*.25D0, .25D0,
     2          3*.75D0, 4*0.D0, 2*.6667D0, 8*0.D0,
     2          .3333D0,.0D0,.75D0,.0D0,2*0.D0,
     2          .6667D0,.75D0,0.D0,3*.75D0, .75D0,
     3          17*.0D0,.5556D0,3*0.D0, 8*0.D0,1*0.D0,210*0.D0/
C---Scattering length FD and effective radius RD in fm vs (LL,ISPIN)
      DATA FDH/17.0D0,7.77D0,23.7D0,2230.1218D0,.225D0,.081D0,-.063D0,
     1     -.65D0,-2.73D0,
     1     .137D0,-.071D0,-.148D0,.112D0,2*1.D-6,-.360D0,
     1     2*1.D-6,1.344D0,6*1.D-6,-5.628D0,2.18D0,2.40D0,
     1     2.81D0,              ! ND potential
C     1     0.50D0,              ! NSC97e potential lam-lam
     1     1*0.001D0,
C     c     2 -10.8D0,2*-5.4D0,4*0.D0,-6.35D0,-11.88D0,8*0.D0,9*0.D0,
     2     3*-5.4D0,4*0.D0,-6.35D0,-11.88D0,8*0.D0,9*0.D0,
     2     1.93D0,1.84D0,
     2     0.50D0,              ! triplet f0 lam-lam=singlet f0 ND
                                ! not contributing in s-wave FSI approx.
     2     1*0.001D0,
     3     240*0.D0/
c--------|---------|---------|---------|---------|---------|---------|----------
      DATA RDH/2.7D0,2.77D0,2.7D0,1.12139906D0,-44.36D0,64.0D0,784.9D0,
     1     477.9D0, 2.27D0, 9*0.D0,-69.973D0, 6*0.D0,3.529D0,
     1     3.19D0,3.15D0,
     1     2.95D0,              ! ND potential lam-lam
C     1  10.6D0, ! NSC97e potential lam-lam
     1     1*0.D0,
     2     3*1.7D0,4*0.D0,2.0D0,2.63D0, 17*0.D0,3.35D0,3.37D0,
     2     2.95D0,              ! triplet d0 lam-lam=singlet d0 ND
                                ! not contributing in s-wave approx.
     2     1*0.D0,
     3     240*0.D0/
C---  Corresponding square well parameters RB (width in fm) and
C--   EB =SQRT(-AM*U) (in GeV/c); U is the well height
      DATA RBH/2.545739D0,   2.779789D0, 2.585795D0, 5.023544D0,
     1     .124673D0, .3925180D0,.09D0, 2.D0, 4.058058D0, 17*0.D0,
     1     2.252623D0, 2.278575D0,
     1     2.234089D0,          ! ND potential lam-lam
C     1  3.065796D0, ! NSC97e potential lam-lam
     1     1*0.001D0,
     2     3*2.003144D0,
     2     4*0.D0, 2.D0, 4.132163D0, 17*0.D0,
     2     2.272703D0, 2.256355D0,
     2     2.234089D0,          ! triplet potential lam-lam=singlet ND
                                ! not contributing in s-wave FSI approx.
     2     1*0.001D0,
     3     240*0.D0/
      DATA EBH/.1149517D0,    .1046257D0,   .1148757D0, .1186010D0,
     1     .7947389D0,2.281208D0,8.7D0,.4D0,.1561219D0,17*0.D0,
     1     .1013293D0, .1020966D0,
     1     .1080476D0,          ! ND potential lam-lam
C     1    .04115994D0, ! NSC97e potential lam-lam
     1     1*0.001D0,
     2     3*.1847221D0,
     2     4*0.D0, .4D0, .1150687D0, 17*0.D0,
     2     .09736083D0, .09708310D0,
     2     .1080476D0,          ! triplet potential lam-lam= singlet ND
                                ! not contributing in s-wave FSI approx.
     2     1*0.001D0,
     3     240*0.D0/

C++----- add to be able to call several time-------
      integer ifirst
      data ifirst/0/
      ifirst=ifirst+1

C=======< constants >========================
      W=1/.1973D0    ! from fm to 1/GeV
      PI=4*DATAN(1.D0)
      PI2=2*PI
      SPI=DSQRT(PI)
      DR=180.D0/PI   ! from radian to degree
        AC1=1.D10
        AC2=1.D10
C=======< condition de calculs >=============
      NUNIT=11 ! for IBM in Nantes
C      NUNIT=4 ! for SUN in Prague
C++      CALL readint4(NUNIT,'ITEST     ',ITEST)
C++      CALL readint4(NUNIT,'NS        ',NS)
      ITEST=I_ITEST
C      LL=I_LL            ! New parameter init for C++ interface
C      NS=I_NS
      call LLINI(I_LL,I_NS,I_ITEST)
      IF(ITEST.EQ.1)THEN
C++      CALL readint4(NUNIT,'ICH       ',ICH)
C++      CALL readint4(NUNIT,'IQS       ',IQS)
C++      CALL readint4(NUNIT,'ISI       ',ISI)
C++      CALL readint4(NUNIT,'I3C       ',I3C)
       ICH=I_ICH
       IQS=I_IQS
       ISI=I_ISI
       I3C=I_I3C
c      print *,'FSINI parameters = LL, NS : ',LL,NS
c      print *,'                 = ICH    : ',ICH
c      print *,'                 = IQS,ISI: ',IQS,ISI
c      print *,'                 = I3C    : ',I3C
      ENDIF
C============================================
      IF (IFIRST.LE.1) THEN
         DO 3 J1=1,30
            DO 3 J2=1,10
               FDH(J1,J2)=FDH(J1,J2)*W
               RDH(J1,J2)=RDH(J1,J2)*W
 3             RBH(J1,J2)=RBH(J1,J2)*W
c      print *,"FD,RD: ",FDH(J1,J2),RDH(J1,J2)
c      print *,"EB,RB: ",EBH(J1,J2),RBH(J1,J2)

      ENDIF
C===================================
      RETURN
      END
C
      SUBROUTINE LLINI(lll,I_NS,I_ITEST)
C===> Initialisation for a given LL value.
C     ===========================================
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_SPIN/RHO(10)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_FD/FD(10),RD(10)
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMMON/FSI_CONS/PI,PI2,SPI,DR,W
      COMPLEX*16 C
      COMMON/FSI_AA/AA
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_AAPI/AAPI(20,2)/FSI_AAND/AAND(20,4)
      COMMON/FSI_AAPIN/AAPIN(20,2)
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
      COMMON/FSI_FDH/FDH(30,10),RDH(30,10),EBH(30,10),RBH(30,10)
      COMMON/FSI_RHOH/RHOH(30,10)
      COMMON/FSI_AMCH/AM1H(30),AM2H(30),C1H(30),C2H(30),MSPINH(30)
      COMMON/FSI_AAKK/AAKK(9)/FSI_AAPAP/AAPAPR(3,2),AAPAPI(3,2)

C++----- add to be able to call several time-------
      integer ifirst
      data ifirst/0/
      ifirst=ifirst+1

C===> LL - Initialisation ========================================
C---- Setting particle masses and charges
      LL=lll
      NS=I_NS
      AM1=AM1H(LL)
      AM2=AM2H(LL)
      C1=C1H(LL)
      C2=C2H(LL)
C      print *,"llini: ",LL,AM1,AM2,C1,C2
C---> Switches:
C     ISI=1(0)  the strong interaction between the two particles ON (OFF)
C     IQS=1(0)  the quantum statistics ON (OFF);
C               should be OFF for nonidentical particles
C     I3C=1(0)  the Coulomb interaction with the nucleus ON (OFF)
C     I3S=1(0)  the strong interaction with the nucleus ON (OFF)
C     ICH=1(0)  if C1*C2 is different from 0 (is equal to 0)
C-    To switch off the Coulomb force between the two particles
C     put ICH=0 and substitute the strong amplitude parameters by
C     the ones not affected by Coulomb interaction
C---- --------------------------------------------------------------------
      IF (I_ITEST.NE.1) THEN
         ICH=0
         IF(C1*C2.NE.0.D0) ICH=1
         IQS=0
         IF(C1+AM1.EQ.C2+AM2) IQS=1
         I3S=0                  ! only this option is available
         ISI=1
         I3C=1
      ENDIF
C---> Calcul. twice the reduced mass (AM), the relative
C     mass difference (DM) and the Bohr radius (AC)
      AM=2*AM1*AM2/(AM1+AM2)
      AMS=AM*AM
      DM=(AM1-AM2)/(AM1+AM2)
      AC=1.D10
      C12=C1*C2
      IF(C12.NE.0.D0)AC=2*137.036D0/(C12*AM)
C---Setting spin factors
      MSPIN=MSPINH(LL)
c      print *,"LLINI MSPIN: ",MSPIN
      MSP=MSPIN
      DO 91 ISPIN=1,10
 91      RHO(ISPIN)=RHOH(LL,ISPIN)
c         print *,"RHO: ",ISPIN,RHO(ISPIN)
C---> Integration limit AA in the spherical wave approximation
      AA=0.D0
cc      IF(NS.EQ.2.OR.NS.EQ.4)AA=.5D0 !!in 1/GeV --> 0.1 fm
      IF(NS.EQ.2.OR.NS.EQ.4)AA=6.D0 !!in 1/GeV --> 1.2 fm
C---> Setting scatt. length (FD), eff. radius (RD) and, if possible,
C--   also the corresp. square well parameters (EB, RB)
      DO 55 JJ=1,MSP
         ISPIN=JJ
         FD(JJ)=FDH(LL,JJ)
         RD(JJ)=RDH(LL,JJ)
         EB(JJ)=EBH(LL,JJ)
         RB(JJ)=RBH(LL,JJ)
c         print *,"FD,RD,EB,RB: ",FD(JJ),RD(JJ),EB(JJ),RB(JJ)
C---Resets FD and RD for a nucleon-deuteron system (LL=8,9)
      IF(LL.EQ.8.OR.LL.EQ.9)THEN
       JH=LL-7+2*JJ-2
       FD(JJ)=AAND(1,JH)
       RD(JJ)=AAND(2,JH)-2*AAND(3,JH)/AAND(1,JH)
      ENDIF
C---Resets FD and RD for a pion-pion system (LL=5,6,7)
      IF(LL.EQ.5.OR.LL.EQ.6.OR.LL.EQ.7)THEN
       IF(LL.EQ.7)FD(JJ)=AAPI(1,2)/AM
       IF(LL.EQ.5)FD(JJ)=(.6667D0*AAPI(1,1)+.3333D0*AAPI(1,2))/AM
       IF(LL.EQ.6)FD(JJ)=(.3333D0*AAPI(1,1)+.6667D0*AAPI(1,2))/AM
       AKS=0.D0
       DAKS=1.D-5
       AKSH=AKS+DAKS
       AKH=DSQRT(AKSH)
       GPI1H=GPIPI(AKSH,1)
       GPI2H=GPIPI(AKSH,2)
       H=1/FD(JJ)
       IF(LL.EQ.7)C(JJ)=1/DCMPLX(GPI2H,-AKH)
       IF(LL.EQ.5)
     + C(JJ)=.6667D0/DCMPLX(GPI1H,-AKH)+.3333D0/DCMPLX(GPI2H,-AKH)
       IF(LL.EQ.6)
     + C(JJ)=.3333D0/DCMPLX(GPI1H,-AKH)+.6667D0/DCMPLX(GPI2H,-AKH)
       HH=DREAL(1/C(JJ))
       RD(JJ)=2*(HH-H)/DAKS
      ENDIF
C---Resets FD and RD for a pion-nucleon system (LL=12,13)
C *** SWITCH OFF SPECIAL PiN ***
C      IF(LL.EQ.12.OR.LL.EQ.13)THEN
      IF(LL.EQ.32.OR.LL.EQ.33)THEN
       IF(LL.EQ.12)FD(JJ)=AAPIN(1,2)
       IF(LL.EQ.13)FD(JJ)=(.6667D0*AAPIN(1,1)+.3333D0*AAPIN(1,2))
       AKS=0.D0
       DAKS=1.D-5
       AKSH=AKS+DAKS
       AKH=DSQRT(AKSH)
       GPI1H=GPIN(AKSH,1)
       GPI2H=GPIN(AKSH,2)
       H=1/FD(JJ)
       IF(LL.EQ.12)C(JJ)=1/DCMPLX(GPI2H,-AKH)
       IF(LL.EQ.13)
     + C(JJ)=.6667D0/DCMPLX(GPI1H,-AKH)+.3333D0/DCMPLX(GPI2H,-AKH)
       HH=DREAL(1/C(JJ))
       RD(JJ)=2*(HH-H)/DAKS
c       print *,"FD,RD,EB,RB: ",FD(JJ),RD(JJ),EB(JJ),RB(JJ)
      ENDIF
C---fm to 1/GeV for pp-bar system
      IF(LL.EQ.30)THEN
      DO 4 I3=1,3
      AAPAPR(I3,JJ)=AAPAPR(I3,JJ)*W
 4    AAPAPI(I3,JJ)=AAPAPI(I3,JJ)*W
C---Calculates complex elements M11=M22=C(6), M12=M21=C(7) for I=0
C---   at k*=0                  M11=M22=C(8), M12=M21=C(9) for I=1
        C(7+(JJ-1)*2)=2*DCMPLX(AAPAPR(1,JJ),AAPAPI(1,JJ))*
     *                DCMPLX(AAPAPR(2,JJ),AAPAPI(2,JJ)) ! 2a_0Sa_1S
        C(6+(JJ-1)*2)=DCMPLX(AAPAPR(1,JJ)+AAPAPR(2,JJ),
     ,                     AAPAPI(1,JJ)+AAPAPI(2,JJ))/
     /                                  C(7+(JJ-1)*2)   ! M11=M22
      C(7+(JJ-1)*2)=-DCMPLX(AAPAPR(1,JJ)-AAPAPR(2,JJ),
     ,                      AAPAPI(1,JJ)-AAPAPI(2,JJ))/
     /                                   C(7+(JJ-1)*2)  ! M12=M21
      ENDIF
C---Calculation continues for any system (any LL)
 55   CONTINUE
      RETURN
      END
C=======================================================
C

c++  This routine is used to init mass and charge of the nucleus.

      SUBROUTINE FSINUCL(R_AMN,R_CN)

      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2

      AMN=R_AMN
      CN=R_CN

      RETURN
      END

C======================================================
C

      SUBROUTINE FSIMOMENTUM(PP1,PP2)

      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  ! particle momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2


      REAL*8 PP1(3),PP2(3)
c      Print *,"momentum",pp1,pp2
      P1X=PP1(1)
      P1Y=PP1(2)
      P1Z=PP1(3)
      P2X=PP2(1)
      P2Y=PP2(2)
      P2Z=PP2(3)
      RETURN
      END



C======================================================
C

      SUBROUTINE FSIPOSITION(XT1,XT2)

      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, !4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2

      REAL*8 XT1(4),XT2(4)
clc      print *,'fsi',xt1,xt2
      X1=XT1(1)
      Y1=XT1(2)
      Z1=XT1(3)
      T1=XT1(4)
      X2=XT2(1)
      Y2=XT2(2)
      Z2=XT2(3)
      T2=XT2(4)
      RETURN
      END


C======================================================
C======================================================
C
      subroutine BoostToPrf()
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_CVK/V,CVK
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  !part. momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, ! 4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_P12/P12X,P12Y,P12Z,E12,P12,AM12,EPM

      XS=X1-X2
      YS=Y1-Y2
      ZS=Z1-Z2
      TS=T1-T2
      RS12=XS*P12X+YS*P12Y+ZS*P12Z
      H1=(RS12/EPM-TS)/AM12
      X=XS+P12X*H1
      Y=YS+P12Y*H1
      Z=ZS+P12Z*H1
      T=(E12*TS-RS12)/AM12
      RPS=X*X+Y*Y+Z*Z
      RP=DSQRT(RPS)
CW      WRITE(6,38)'RP ',RP,'X ',X,Y,Z,T
38    FORMAT(A7,E11.4,A7,4E11.4)

      CVK=(P12X*PPX+P12Y*PPY+P12Z*PPZ)/(P12*AK)
      V=P12/E12
      return
      end

      SUBROUTINE FSIWF(WEI)
C==> Prepares necessary quantities and call VZ(WEI) to calculate
C    the weight due to FSI
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_CVK/V,CVK
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  !part. momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, ! 4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_SPIN/RHO(10)
      COMMON/FSI_BP/B,P
      COMMON/FSI_ETA/ETA
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_RR/F(10)
      COMMON/FSI_FD/FD(10),RD(10)
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMPLEX*16 C,F
      COMMON/FSI_AA/AA
      COMMON/FSI_SHH/SH,CHH
      COMMON/FSI_AAPI/AAPI(20,2)/FSI_AAND/AAND(20,4)
      COMMON/FSI_AAPIN/AAPIN(20,2)
      COMMON/FSI_P12/P12X,P12Y,P12Z,E12,P12,AM12,EPM
      COMMON/FSI_2CHA/AK2,AK2S,AMU2_AMU1 ! k* (kappa) for 2-nd channel
C==>calculating relative 4-coordinates of the particles in PRF
C-  {T,X,Y,Z} from the relative coordinates {TS,XS,YS,ZS} in NRF
      XS=X1-X2
      YS=Y1-Y2
      ZS=Z1-Z2
      TS=T1-T2
      RS12=XS*P12X+YS*P12Y+ZS*P12Z
      H1=(RS12/EPM-TS)/AM12
      X=XS+P12X*H1
      Y=YS+P12Y*H1
      Z=ZS+P12Z*H1
      T=(E12*TS-RS12)/AM12
      RPS=X*X+Y*Y+Z*Z
      RP=DSQRT(RPS)
CW      WRITE(6,38)'RP ',RP,'X ',X,Y,Z,T
38    FORMAT(A7,E11.4,A7,4E11.4)

      CVK=(P12X*PPX+P12Y*PPY+P12Z*PPZ)/(P12*AK)
      V=P12/E12

      ACH=1.D0
      IF(ICH.EQ.0)GOTO 21
      XH=AC*AK
      ACH=ACP(XH)
      ACHR=DSQRT(ACH)
      ETA=0.D0
      IF(XH.NE.0.D0)ETA=1/XH
C---HCP, HPR needed (e.g. in GST) if ICH=1
      HCP=HC(XH)
c      PRINT *,"HC ", HCP
      HPR=HCP+.1544313298D0
c      PRINT *,"HPR ", HPR
      AAK=ACH*AK    !
      HCP2=2*HCP/AC ! needed to calculate C(JJ) for charged particles
  21  CONTINUE
      MSP=MSPIN
      DO 30 JJ=1,MSP
      ISPIN=JJ
      IF(NS.NE.1)GOTO22
C---Calc. quantities for the square well potential;
C-- for LL=6-26 the square well potential is not possible or available
      IF(LL.EQ.4)GOTO 22
      BK(JJ)=DSQRT(EB(JJ)**2+AKS)
      XRA=2*RB(JJ)/AC
      HRA=BK(JJ)*RB(JJ)
      CALL SEQ(XRA,HRA)
      SBKRB(JJ)=HRA*B
      HRA=AK*RB(JJ)
      CALL GST(XRA,HRA)
      SDK(JJ)=SH
      CDK(JJ)=CHH
c      PRINT *,"G ",CHH SH
      SDKK(JJ)=RB(JJ)
      IF(AK.NE.0.D0)SDKK(JJ)=SH/AK
      IF(ICH.EQ.1)SDK(JJ)=ACH*SDK(JJ)
  22  CONTINUE
C-----------------------------------------------------------------------
C---Calc. the strong s-wave scattering amplitude = C(JJ)
C-- divided by Coulomb penetration factor squared (if ICH=1)
      IF(NS.NE.1)GOTO 230
      IF(LL.NE.4)GOTO 230 ! SW scat. amplitude used for alfa-alfa only
      GAK=G(AK)
      AKACH=AK
      IF(ICH.EQ.1)AKACH=AK*ACH
      C(JJ)=1/DCMPLX(GAK,-AKACH) ! amplitude for the SW-potential
      GOTO 30
 230  IF(LL.EQ.5.OR.LL.EQ.6.OR.LL.EQ.7)GOTO20    ! pipi
C *** SWITCH OFF SPECIAL PiN ***
C      IF(LL.EQ.12.OR.LL.EQ.13)GOTO20             ! piN
      IF(LL.EQ.32.OR.LL.EQ.33)GOTO20             ! piN
      IF(LL.EQ.8.OR.LL.EQ.9.OR.LL.EQ.18)GOTO20   ! Nd, dd
      IF(LL.EQ.14.OR.LL.EQ.17.OR.LL.EQ.23)GOTO27 ! K+K-, K-p, K0K0-b
	IF(LL.EQ.30)GOTO 28                        ! pp-bar
      A1=RD(JJ)*FD(JJ)*AKS
      A2=1+.5D0*A1
      IF(ICH.EQ.1)A2=A2-2*HCP*FD(JJ)/AC
      AKF=AK*FD(JJ)
      IF(ICH.EQ.1)AKF=AKF*ACH
      C(JJ)=FD(JJ)/DCMPLX(A2,-AKF)
      GOTO30
 20   CONTINUE
C---Calc. scatt. ampl. C(JJ) for pipi, piN and Nd, dd
      IF(LL.EQ.5)THEN
       CALL CPIPI  ! pi+pi- C(1)=f(pi+pi- -> pi+pi-)
c                           C(2)=f(pi0pi0 -> pi+pi-)
       GOTO 30
      ENDIF
      JH=LL-7+2*JJ-2
      IF(LL.EQ.8.OR.LL.EQ.9)GPI2=GND(AKS,JH)
      IF(LL.EQ.18)GPI2=GDD(AKS,JJ)
      IF(LL.EQ.5.OR.LL.EQ.6.OR.LL.EQ.7)GPI2=GPIPI(AKS,2)
C *** SWITCH OFF SPECIAL PiN ***
C      IF(LL.EQ.12.OR.LL.EQ.13)GPI2=GPIN(AKS,2)
      IF(LL.EQ.32.OR.LL.EQ.33)GPI2=GPIN(AKS,2)
      C(JJ)=1.D0/DCMPLX(GPI2,-AK) !pi+pi+, nd, pd, pi+p, dd
C *** SWITCH OFF SPECIAL PiN ***
C      IF(LL.NE.5.AND.LL.NE.6.AND.LL.NE.13)GOTO27
      IF(LL.NE.5.AND.LL.NE.6.AND.LL.NE.33)GOTO27
      IF(LL.EQ.5.OR.LL.EQ.6)GPI1=GPIPI(AKS,1)
C *** SWITCH OFF SPECIAL PiN ***
C      IF(LL.EQ.13)GPI1=GPIN(AKS,1)
C      IF(LL.EQ.5.OR.LL.EQ.13)
      IF(LL.EQ.5.OR.LL.EQ.33)
     c           C(JJ)=.6667D0/DCMPLX(GPI1,-AK)+.3333D0*C(JJ) !pi+pi-,pi-p
      IF(LL.EQ.6)C(JJ)=.3333D0/DCMPLX(GPI1,-AK)+.6667D0*C(JJ) !pi0pi0
 27   CONTINUE
C---Calc. K+K-, K0K0-b or K-p s-wave scatt. ampl.
      IF(LL.EQ.14.OR.LL.EQ.23)CALL CKKB
c      IF(LL.EQ.17)C(JJ)=DCMPLX(-3.29D0,3.55D0)     ! Martin'76 (K-p)
c      IF(LL.EQ.17)C(JJ)=DCMPLX(3.29D0,3.55D0)     ! Martin'76 (K-p) WRONG SIGN!!!
c      IF(LL.EQ.17)C(JJ)=DCMPLX(-2.585D0,4.156D0)   ! Borasoy'04 (K-p)
      IF(LL.EQ.17)C(JJ)=DCMPLX(-5.321D0,3.8D0)   ! Borasoy'06 (K-p) New, better fit
c      IF(LL.EQ.17)C(JJ)=DCMPLX(-3.371D0,3.244D0)   ! Martin'81 (K-p)
C---Calc. pi+pi-, pi+pi+, pd, pi+p, pi-p, K+K- or K-p s-wave scatt. ampl.
C-- divided by Coulomb penetration factor squared (if ICH=1)
      IF(ICH.EQ.0)GOTO 30
      AAK=ACH*AK    !
      HCP2=2*HCP/AC ! needed to calculate C(JJ) for charged particles
      C(JJ)=1/(1/C(JJ)-HCP2+DCMPLX(0.D0,AK-AAK))
 28   CONTINUE
C---Calc. pp-bar s-wave scatt. ampl.
      CALL CPAP
 30   CONTINUE
C***********************************************************************
      CALL VZ(WEI)
      RETURN
      END
      SUBROUTINE VZ(WEI)
C==> Calculates the weight WEI due to FSI
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_JR/JRAT
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_SPIN/RHO(10)
      COMMON/FSI_ETA/ETA
      COMMON/FSI_AA/AA
      COMMON/FSI_FFF/F12,F21
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_FD/FD(10),RD(10)
      COMMON/FSI_RR/F(10)
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMMON/FSI_COULPH/EIDC
      COMPLEX*16 F,C,G,PSI12,PSI21
      COMPLEX*16 F12,F21
      COMPLEX*16 EIDC
      COMPLEX*8 Z8,CGAMMA
      COMMON/FSI_FFPN/FF12,FF21
      COMPLEX*16 FF12,FF21
      COMMON/FSI_2CHA/AK2,AK2S,AMU2_AMU1 ! k* (kappa) for 2-nd channel
      WEI=0.D0
      IF(JRAT.EQ.1)GOTO 11
      RHOS=AK*RP
      HS=X*PPX+Y*PPY+Z*PPZ
      IF(RHOS.LT.15.D0.AND.RHOS+DABS(HS).LT.20.D0)GOTO 2
C---Calc. EIDC=exp(i*Coul.Ph.);
C-- used in calc. of hypergeom. f-s in SEQA, FAS at k*R > 15, 20
      Z8=CMPLX(1.,SNGL(ETA))
      Z8=CGAMMA(Z8)
      EIDC=Z8/CABS(Z8)
 2    CALL FF(RHOS,HS)
 11   MSP=MSPIN
      IF(ISI.EQ.0)GOTO 4  ! the strong interaction ON (OFF) if ISI=1(0)
      IF(RP.LT.AA)GOTO 4
      IF(JRAT.NE.1) CALL FIRT
      IF(IQS.EQ.0)GOTO 5  ! the quantum statistics ON (OFF) if IQS=1(0)
      JSIGN=-1
      DO 1 JJ=1,MSP
      JSIGN=-JSIGN
      G=F(JJ)*C(JJ)
      IF(ICH.EQ.1)G=G*ACHR
      PSI12=FF12*(F12+G)
      PSI21=FF21*(F21+G)
      G=PSI12+JSIGN*PSI21
 1    WEI=WEI+RHO(JJ)*(DREAL(G)**2+DIMAG(G)**2)
      GOTO 8
 5    DO 6 JJ=1,MSP
      G=F(JJ)*C(JJ)
      IF(ICH.EQ.1)G=G*ACHR
CW      WRITE(6,38)'JJ ',JJ,'F ',DREAL(F(JJ)),DIMAG(F(JJ))
CW      WRITE(6,38)'JJ ',JJ,'C ',DREAL(C(JJ)),DIMAG(C(JJ))
CW      WRITE(6,38)'JJ ',JJ,'G ',DREAL(G),DIMAG(G)
CW      WRITE(6,38)'JJ ',JJ,'F12+G ',DREAL(F12+G),DIMAG(F12+G)
CW      WRITE(6,38)'JJ ',JJ,'F21+G ',DREAL(F21+G),DIMAG(F21+G)
38    FORMAT(A7,I3,A7,2E11.4)
      PSI12=FF12*(F12+G)
 6    WEI=WEI+RHO(JJ)*(DREAL(PSI12)**2+DIMAG(PSI12)**2)
c--- Account for pi0pi0 -> pi+pi- or nn-bar -> pp-bar channel -----
      IF(LL.EQ.5.OR.LL.EQ.30)THEN
      DO 61 JJ=1,MSP
        HH=RHO(JJ)*(DREAL(C(JJ+MSP))**2+DIMAG(C(JJ+MSP))**2)*
     *  AMU2_AMU1*ACH/RPS
        IF(AK2S.LT.0)HH=HH*DEXP(-2*RP*AK2)
 61   WEI=WEI+HH
	ENDIF
c------------------------------------------------------------------
      RETURN
 4    PSI12=FF12*F12
C      PRINT *,"FF12 ",FF12
C      PRINT *,"PSI12 ",PSI12
      IF(IQS.EQ.0)GOTO 50 ! the quantum statistics ON (OFF) if IQS=1(0)
      PSI21=FF21*F21
C      PRINT *,"FF21 ",FF21
C      PRINT *,"PSI21 ",PSI21
      JSIGN=-1
      DO 3 JJ=1,MSP
      JSIGN=-JSIGN
      G=PSI12+JSIGN*PSI21
C      PRINT *,"G",G
 3    WEI=WEI+RHO(JJ)*(DREAL(G)**2+DIMAG(G)**2)
C      PRINT *,"WEI",WEI/2
C      PRINT *,"1/Ka",(1/(AK*AC))
      GOTO 8
 50   WEI=DREAL(PSI12)**2+DIMAG(PSI12)**2
      RETURN
 8    WEI=WEI/2
      RETURN
      END
      SUBROUTINE FIRT
C---CALC. THE F(JJ)
C-- F(JJ)*C(JJ)= DEVIATION OF THE BETHE-SALPETER AMPL. FROM PLANE WAVE
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_SHH/SH,CHH
      COMMON/FSI_BP/B,P
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_RR/F(10)
      EQUIVALENCE(RSS,RP),(TSS,T)
      COMPLEX*16 F,C,CH1
      MSP=MSPIN
      DO 10 JJ=1,MSP
      IF(JJ.GT.1)GOTO 3
      XRA=2*RSS/AC
      IF(AK.NE.0.D0)GOTO2
      SHK=1.D0
      SH=.0D0
      SHH=SH
      CHH=1/RSS
      GOTO3
  2   H=AK*RSS
      CALL GST(XRA,H)
      SH=SH/RSS
      CHH=CHH/RSS
      SHH=SH
c      PRINT *,"G ",CHH, SH
      IF(ICH.EQ.1) SHH=ACH*SH
  3   IF(NS.EQ.2)GOTO1
C---F= ASYMPTOTIC FORMULA (T= 0 APPROX.); NS= 4
  6   F(JJ)=DCMPLX(CHH,SHH)
      IF(NS.NE.1)GOTO 10
C---F INSIDE THE SQUARE-WELL (T= 0 APPROX.); NS= 1
      IF(RSS.GE.RB(JJ)) GOTO 10
      IF(AK.NE.0.D0.AND.JJ.EQ.1)SHK=SH/AK
      H=BK(JJ)*RSS
      CALL GST(XRA,H)
      SKR=B*BK(JJ)
      F(JJ)=DCMPLX(CDK(JJ),SDK(JJ))*SKR
      CH1=(SDKK(JJ)*SKR-SHK*SBKRB(JJ))/C(JJ)
      F(JJ)=(F(JJ)+CH1)/SBKRB(JJ)
      GOTO 10
  1   CONTINUE
C---F= ASYMPTOTIC FORMULA (T= 0 NOT REQUIRED); NS= 2
      IF(JJ.GT.1)GOTO 8
      IF(TSS.EQ.0.D0)GOTO6
      TSSA=DABS(TSS)
      IF(DM.NE.0.D0)GOTO 11
      H=AM*.5D0/TSSA
      IF(AK.NE.0.D0)GOTO4
      HM=H*RPS
      IF(HM.GE.3.D15)GOTO6
      FS1=DFRSIN(HM)
      FC1=DFRCOS(HM)
      FC2=FC1
      FS2=FS1
      GOTO5
  4   CONTINUE
      H1=AK*TSSA/AM
      HM=H*(RSS-H1)**2
      HP=H*(RSS+H1)**2
      IF(HP.GE.3.D15)GOTO6
      FS1=DFRSIN(HM)
      FC1=DFRCOS(HM)
      FS2=DFRSIN(HP)
      FC2=DFRCOS(HP)
      GOTO 5
  11  CONTINUE
      FS1=0.D0
      FS2=0.D0
      FC1=0.D0
      FC2=0.D0
      DO 13 I=1,2
      IF(I.EQ.1)TSSH=TSSA*(1+DM)
      IF(I.EQ.2)TSSH=TSSA*(1-DM)
      H=AM*.5D0/TSSH
      IF(AK.NE.0.D0)GOTO 12
      HM=H*RPS
      IF(HM.GE.3.D15)GOTO6
      FS1=FS1+DFRSIN(HM)/2
      FC1=FC1+DFRCOS(HM)/2
      IF(I.EQ.1)GOTO 13
      FC2=FC1
      FS2=FS1
      GOTO 13
  12  CONTINUE
      H1=AK*TSSH/AM
      HM=H*(RSS-H1)**2
      HP=H*(RSS+H1)**2
      IF(HP.GE.3.D15)GOTO6
      FS1=FS1+DFRSIN(HM)/2
      FC1=FC1+DFRCOS(HM)/2
      FS2=FS2+DFRSIN(HP)/2
      FC2=FC2+DFRCOS(HP)/2
  13  CONTINUE
  5   C12=FC1+FS2
      S12=FS1+FC2
      A12=FS1-FC2
      A21=FS2-FC1
      A2=.5D0*(CHH*(A12+A21)+SH*(A12-A21))+SHH
      A1=.5D0*(CHH*(C12+S12)+SH*(C12-S12))
      F(JJ)=.3989422D0*DCMPLX(A1,A2)
      GOTO 10
  8   F(JJ)=F(1)
 10   CONTINUE
      RETURN
      END
      FUNCTION EXF(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X.LT.-15.D0) GO TO 1
      EXF=DEXP(X)
      RETURN
  1   EXF=.0D0
      RETURN
      END
      SUBROUTINE SEQ(X,H)
C---CALC. FUNCTIONS B, P (EQS. (17) OF G-K-L-L);
C-- NEEDED TO CALC. THE CONFLUENT HYPERGEOMETRIC FUNCTION GST.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_BP/B,P
      DIMENSION BH(3),PH(3)
      DATA ERR/1.D-7/
      BH(1)=1.D0
      PH(1)=1.D0
      PH(2)=.0D0
      BH(2)=.5D0*X
      B=1+BH(2)
      P=1.D0
      HS=H*H
      J=0
  2   J=J+1
      BH(3)=(X*BH(2)-HS*BH(1))/((J+1)*(J+2))
      PH(3)=(X*PH(2)-HS*PH(1)-(2*J+1)*X*BH(2))/(J*(J+1))
      B=B+BH(3)
      P=P+PH(3)
      Z=DABS(BH(2))+DABS(BH(3))+DABS(PH(2))+DABS(PH(3))
      IF(Z.LT.ERR)RETURN
      BH(1)=BH(2)
      BH(2)=BH(3)
      PH(1)=PH(2)
      PH(2)=PH(3)
      GOTO 2
      END
      SUBROUTINE SEQA(X,H)
C---CALC. FUNCTIONS CHH=REAL(GST), SH=IMAG(GST)/ACH, B=SH/H
C-- IN THE ASYMPTOTIC REGION H=K*R >> 1.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_BP/B,P
      COMMON/FSI_SHH/SH,CHH
      COMMON/FSI_ETA/ETA
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_COULPH/EIDC
      COMPLEX*16 EIDC,GST
      ARG=H-ETA*DLOG(2*H)
      GST=DCMPLX(DCOS(ARG),DSIN(ARG))
      GST=ACHR*EIDC*GST
      CHH=DREAL(GST)
      SH=DIMAG(GST)/ACH
      B=SH/H
      RETURN
      END
      SUBROUTINE FF(RHO,H)
C---Calc. F12, F21;
C-- F12= FF0* plane wave,  FF0=F*ACHR,
C---F is the confluent hypergeometric function,
C-- ACHR=sqrt(ACH), where ACH is the Coulomb factor
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ETA/ETA
      COMMON/FSI_FFF/F12,F21
      COMPLEX*16 FF0,F12,F21
C      PRINT *,"FF for ",RHO,H
      C=DCOS(H)
      S=DSIN(H)
C      PRINT *,"Cos Sin ",C,S
      F12=DCMPLX(C,-S)
      F21=DCMPLX(C,S)
      IF(ICH.EQ.0)RETURN
      RHOP=RHO+H
      RHOM=RHO-H
      F12=FF0(RHO,H)*F12
      F21=FF0(RHO,-H)*F21
c      PRINT *,"F12 ",F12
c      PRINT *,"F21 ",F21
      RETURN
      END
      FUNCTION FAS(RKS)
C-- FAS=F*ACHR
C---F is the confluent hypergeometric function at k*r >> 1
C-- ACHR=sqrt(ACH), where ACH is the Coulomb factor
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 FAS,EIDC,ZZ1
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_ETA/ETA
      COMMON/FSI_COULPH/EIDC
      D1=DLOG(RKS)*ETA
      D2=ETA*ETA/RKS
      ZZ1=DCMPLX(DCOS(D1),DSIN(D1))/EIDC
      FAS=DCMPLX(1.D0,-D2)*ZZ1
      FAS=FAS-DCMPLX(DCOS(RKS),DSIN(RKS))*ETA/RKS/ZZ1
      RETURN
      END
      FUNCTION FF0(RHO,H)
C-- FF0=F*ACHR
C-- F is the confluent hypergeometric function
C-- (Eq. (15) of G-K-L-L), F= 1 at r* << AC
C-- ACHR=sqrt(ACH), where ACH is the Coulomb factor
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ETA/ETA
      COMPLEX*16 ALF,ALF1,Z,S,A,FF0,FAS
      DATA ERR/1.D-7/
      S=DCMPLX(1.D0,0.D0)
      FF0=S
      RHOP=RHO+H
CC      GOTO 5 ! rejects the approx. calcul. of hyperg. f-ion F
      IF(RHOP.LT.20.D0)GOTO5
      FF0=FAS(RHOP) ! approx. calc.
      RETURN
  5   ALF=DCMPLX(.0D0,-ETA)
      ALF1=ALF-1
      Z=DCMPLX(.0D0,RHOP)
      J=0
  3   J=J+1
      A=(ALF1+J)/(J*J)
      S=S*A*Z
      FF0=FF0+S
      ZR=DABS(DREAL(S))
      ZI=DABS(DIMAG(S))
      IF((ZR+ZI).GT.ERR)GOTO3
C      PRINT *,"FF0 ",FF0
C      PRINT *,"ACHR ",ACHR
      FF0=FF0*ACHR
      RETURN
      END
      FUNCTION HC(XA)
C---HC = h-function of Landau-Lifshitz: h(x)=Re[psi(1-i/x)]+ln(x)
C-- psi(x) is the digamma function (the logarithmic derivative of
C-- the gamma function)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BN(15)
      DATA BN/.8333333333D-1,.8333333333D-2,.396825396825D-2,
     1        .4166666667D-2,.7575757576D-2,.2109279609D-1,
     2        .8333333333D-1,.4432598039D0 ,.305395433D1,
     3        .2645621212D2, .2814601449D3, .3607510546D4,
     4        .5482758333D5, .9749368235D6, .200526958D8/
      X=DABS(XA)
      IF(X.LT..33D0) GOTO 1
CC      IF(X.GE.3.5D0) GO TO 2
      S=.0D0
      N=0
   3  N=N+1
      DS=1.D0/N/((N*X)**2+1)
      S=S+DS
      IF(DS.GT.0.1D-12) GOTO 3
C---Provides 7 digit accuracy
      HC=S-.5772156649D0+DLOG(X)
      RETURN
CC   2  HC=1.2D0/X**2+DLOG(X)-.5772156649 D0
CC      RETURN
   1  X2=X*X
      XP=X2
      HC=0.D0
      IMA=9
      IF(X.LT.0.1D0)IMA=3
      DO 4 I=1,IMA
      HC=HC+XP*BN(I)
   4  XP=X2*XP
      RETURN
      END
      FUNCTION ACP(X)
C--- ACP = COULOMB PENETRATION FACTOR
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X.LT.0.05D0.AND.X.GE.0.D0) GO TO 1
      Y=6.2831853D0/X
      ACP=Y/(EXF(Y)-1)
      RETURN
   1  ACP=1.D-6
      RETURN
      END
      SUBROUTINE GST(X,H)
C---CALC. THE CONFL. HYPERGEOM. F-N = CHH+i*SH
C-- AND THE COULOMB F-S B, P (CALLS SEQ OR SEQA).
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_SHH/SH,CHH
      COMMON/FSI_BP/B,P
  1   IF(ICH.EQ.1)GOTO 2
  3   SH=DSIN(H)
      CHH=DCOS(H)
      B=1.D0
      IF(H.NE.0.D0)B=SH/H
      P=CHH
      RETURN
  2   CONTINUE
      IF(H.GT.15.D0)GOTO4 ! comment out if you want to reject
                   ! the approximate calculation of hyperg. f-ion G
      CALL SEQ(X,H) ! exact calculation
      SH=H*B
      CHH=P+B*X*(DLOG(DABS(X))+HPR)
c      PRINT *,"CH SH ",CHH, SH
c      PRINT *,"X H ",X, H
c      PRINT *,"P B ",P, B
      RETURN
  4   CALL SEQA(X,H)
      RETURN
      END
      FUNCTION FF1(RHO,H)
C---FF1=FF0; used for particle-nucleus system
C-- FF0=F12*ACHR
C-- F12 is the confluent hypergeometric function
C-- (Eq. (15) of G-K-L-L), F12= 1 at r* << AC
C-- ACHR=sqrt(ACH), where ACH is the Coulomb factor
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ETA/ETA
      COMMON/FSI_COULPH/EIDC
      COMMON/FSI_ICH1/ICH1
      COMPLEX*16 FF0,FF1
      COMPLEX*16 EIDC
      COMPLEX*8 Z8,CGAMMA
      FF1=DCMPLX(1.D0,0.D0)
      IF(ICH1.EQ.0)GOTO 2
      IF(RHO.LT.15.D0.AND.RHO+H.LT.20.D0)GOTO 2
C---Calc. EIDC=exp(i*Coul.Ph.);
C-- used in calc. of hypergeom. f-s in SEQA, FAS at k*R > 15, 20
      Z8=CMPLX(1.,SNGL(ETA))
      Z8=CGAMMA(Z8)
      EIDC=Z8/CABS(Z8)
 2    FF1=FF0(RHO,H)
      RETURN
      END

      FUNCTION G(AK)
C---Used to calculate SW scattering amplitude for alpa-alpha system
C-- and for sear.f (square well potential search)
C---NOTE THAT SCATT. AMPL.= 1/CMPLX(G(AK),-AK*ACH)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,HCP2,AAK,JJ,MSPIN
      COMMON/FSI_BP/B,P/FSI_DERIV/BPR,PPR/FSI_SHH/SH,CHH
      COMMON/FSI_DAK/DAK,IFUN
      HCP2=.0D0
      ACH=1.D0
      IF(ICH.EQ.0)GOTO 1
      XH=AC*AK
      HCP=HC(XH)
      HPR=HCP+.1544313298D0
      ACH=ACP(XH)
      HCP2=2*HCP/AC
   1  AKS=AK**2
      BK(JJ)=DSQRT(AKS+EB(JJ)**2) ! kappa=kp
      X=2*RB(JJ)/AC
      H=BK(JJ)*RB(JJ)             ! kp*d
      CALL GST(X,H)
      BRHO=B                      ! B(kp,d)
      SBKRB(JJ)=SH                ! kp*d*B(kp,d)
      CALL DERIW(X,H)
      BRHOP=BPR                   ! B'(kp,d)= dB(kp,r)/dln(r) at r=d
      H=AK*RB(JJ)
      CALL GST(X,H)
      CDK(JJ)=CHH                 ! ReG(k,d)
      BRHOS=B                     !  B(k,d)
      PRHOS=P                     !  P(k,d)
      SDK(JJ)=SH
      IF(ICH.EQ.0)GOTO 2
      SDK(JJ)=ACH*SH              ! ImG(k,d)
      IF(AK.EQ.0.D0.AND.AC.LT.0.D0)SDK(JJ)=3.14159*X*B
   2  SDKK(JJ)=RB(JJ)
      IF(AK.NE.0.D0)SDKK(JJ)=SH/AK ! d*B(k,d)
      CALL DERIW(X,H)              ! PPR=P'(k,d)= dP(k,r)/dln(r) at r=d
      ZZ=PPR-PRHOS
      IF(ICH.EQ.1)ZZ=ZZ+X*(BRHOS+BPR*(DLOG(DABS(X))+HPR))
C   ZZ= P'(k,d)-P(k,d)+x*{B(k,d)+B'(k,d)*[ln!x!+2*C-1+h(k*ac)]}
      GG=(BRHOP*CDK(JJ)-BRHO*ZZ)/RB(JJ)
C   GG= [B'(kp,d)*ReG(k,d)-B(kp,d)*ZZ]/d
      G=GG/(BRHO*BPR-BRHOP*BRHOS)
C    G= GG/[B(kp,d)*B'(k,d)-B'(kp,d)*B(k,d)]
      RETURN
      END
      SUBROUTINE DERIW(X,H)
C---CALLED BY F-N G(AK)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_BP/B,P/FSI_DERIV/BPR,PPR
      HH=.1D-3
      CALL GST(X,H-HH)
      Q1=P
      B1=B
      CALL GST(X,H+HH)
      HHH=HH+HH
      BPR=H*(B-B1)/HHH
      PPR=H*(P-Q1)/HHH
      IF(ICH.EQ.0)RETURN
      CALL GST(X-HH,H)
      Q1=P
      B1=B
      CALL GST(X+HH,H)
      BPR=BPR+X*(B-B1)/HHH
      PPR=PPR+X*(P-Q1)/HHH
      RETURN
      END
C================================================================

c      SUBROUTINE READ_FILE(KEY,CH8,INT4,REAL8,IERR,NUNIT)
C     ==========
C
C     Routine to read one parameter of the program in the file
C     DATA NUNIT defined in FSI3B EXEC
C     NUNIT=11 for IBM in Nantes, 4 for SUN in Prague
C
C     INPUT  : KEY (CHARACTER*10) :
C     OUTPUT : case of KEY : CH8   : (CHARACTER*8)
C                            INT4  : (INTEGER*4)
C                            REAL8 : (REAL*8)
C                     (only one of them)
C              IERR (INTEGER) : 0 : no error
C                               1 : key not found
c
c      CHARACTER*10 KEY,TEST
c      CHARACTER*4  TYPE
c      CHARACTER*8  CH8
c      INTEGER*4    INT4
c      REAL*8       REAL8
c      INTEGER*4    IERR
c
c      IERR=0
c      REWIND(NUNIT)
c1     READ(NUNIT,FMT='(A10,2X,A4)')TEST,TYPE
c      IF (TEST.EQ.KEY) THEN
c        BACKSPACE(NUNIT)
c        IF (TYPE.EQ.'CHAR') READ(NUNIT,FMT='(18X,A8,54X)')CH8
c        IF (TYPE.EQ.'INT4') READ(NUNIT,FMT='(18X,I8,54X)')INT4
c        IF (TYPE.EQ.'REA8') READ(NUNIT,FMT='(18X,F10.5,52X)')REAL8
c      ELSE
c        IF (TEST.NE.'* E.O.F. *') THEN
c          GOTO 1
c        ELSE
c          IERR=1
c        ENDIF
c      ENDIF
cc      IF(IERR.EQ.1)STOP
c      RETURN
c      END

