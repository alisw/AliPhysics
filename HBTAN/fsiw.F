c         1         2         3         4         5         6         7         8
c---------|---------|---------|---------|---------|---------|---------|---------|
*-- Author :    R.Lednicky      20/01/95
      SUBROUTINE fsiw 
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
C-   LL       1  2  3  4   5    6   7  8  9 10  11  12  13  14 15 16 17
C-   part. 1: n  p  n alfa pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-
C-   part. 2: n  p  p alfa pi- pi0 pi+ d  d  K-  K+  p   p  K- K+ p  p
C   NS=1 y/n: +  +  +  +   +    -   -  -  -  -   -   -   -  -  -  -  -
C----------------------------------------------------------------------
C-   LL       18  19 20  21  22 23 24 25 26    27     28
C-   part. 1:  d  d   t  t   K0 K0  d p  p      p      n
C-   part. 2:  d alfa t alfa K0 K0b t t alfa lambda lambda
C   NS=1 y/n:  -  -   -  -   -  -   - -  -      +      +
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
      COMMON/LEDWEIGHT/WEIF,WEI,WEIN,ITEST,IRANPOS
       
C------------------------------------------------------------------
c      write(*,*)'in fsiw C1 C2 CN',C1,C2,CN 

C==> AC1,2 = "relativistic" Bohr radii for particle-nucleus systems

      C1N=C1*CN
      IF(C1N.NE.0.D0)AC1=137.036D0/(C1N*E1) !m1-->E1
      C2N=C2*CN
      IF(C2N.NE.0.D0)AC2=137.036D0/(C2N*E2) !m2-->E2
 
C-----------------------------------------------------------
      CALL FSIPN  !(WEIF) !weight due to particle-nucleus FSI
      JRAT=0
      CALL FSIWF  !(WEI)  !weight due to particle-particle-nucleus FSI
      WEIN=WEI
       IF(I3C*J.NE.0) THEN
      FF12=DCMPLX(1.D0,0.D0)
      FF21=DCMPLX(1.D0,0.D0)
      JRAT=1
      CALL VZ    !(WEIN) ! weight due to particle-particle FSI
       ENDIF
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
          
      SUBROUTINE VZ  !(WEI)
C==> Calculates the weight WEI due to FSI
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_JR/JRAT
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_SPIN/RHO(10)
      COMMON/FSI_ETA/ETA
      COMMON/FSI_AA/AA
      COMMON/FSI_FFF/F12,F21
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_FD/FD(10),RD(10)
      COMMON/FSI_RR/F(10)
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMMON/FSI_COULPH/EIDC
      COMMON/LEDWEIGHT/WEIF,WEI,WEIN,ITEST,IRANPOS
      COMPLEX*16 F,C,G,PSI12,PSI21
      COMPLEX*16 F12,F21
      COMPLEX*16 EIDC
      COMPLEX*8 Z8,CGAMMA
      COMMON/FSI_FFPN/FF12,FF21
      COMPLEX*16 FF12,FF21
c--mlv------------------------------------------------------
c      write(*,*)'------- WE are in VZ-----------'
c      write(*,*)'JRAT=', JRAT
c      write(*,*)'AK=', AK,'RP=',RP
c      write(*,*)'X Y Z', X,Y,Z
c       write(*,*)'PPX PPY PPZ',PPX,PPY,PPZ        
c      write(*,*)'ETA',ETA
c      write(*,*)'MSPIN',MSPIN
c      write(*,*)'ISI',ISI
c       write(*,*)'RP RPS',RP,RPS
c      write(*,*)'AA- BEFORE FIRT------',AA
c      write(*,*)'IQS',IQS
c      write(*,*)'ICH',ICH
c      write(*,*)'MSPIN',MSPIN
c      write(*,*)'C',(C(I),I=1,MSPIN)
c      write(*,*)'ACHR',ACHR
c      write(*,*)'F12',F12,'F21',F21,'FF12',FF12,'FF21',FF21

c--mlv------------------------------------------------------
      WEI=0.D0
      IF(JRAT.EQ.1)GOTO 11
      RHOS=AK*RP
      HS=X*PPX+Y*PPY+Z*PPZ
c      write(*,*)'rhos',rhos,'hs',hs
      IF(RHOS.LT.15.D0.AND.RHOS+DABS(HS).LT.20.D0)GOTO 2
C---Calc. EIDC=exp(i*Coul.Ph.);
C-- used in calc. of hypergeom. f-s in SEQA, FAS at k*R > 15, 20
c      write(*,*)'1111'
      Z8=CMPLX(1.,SNGL(ETA))
      Z8=CGAMMA(Z8) 
      EIDC=Z8/CABS(Z8)
c      write(*,*)'1111 z8 eidc',z8,eidc
      
 2    CALL FF(RHOS,HS)
 11   MSP=MSPIN
      
c      write(*,*)'before go to 4', wei
 
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
 
c      WRITE(*,*)'WEI111=',WEI

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
      RETURN
      
c      write(*,*)'before 4 wei',wei
 4    PSI12=FF12*F12
      IF(IQS.EQ.0)GOTO 50 ! the quantum statistics ON (OFF) if IQS=1(0)
      PSI21=FF21*F21
      JSIGN=-1
      DO 3 JJ=1,MSP
      JSIGN=-JSIGN
      G=PSI12+JSIGN*PSI21
c      write(*,*)'msp jsign psi12 psi21',msp,jsign,psi12,psi21
 3    WEI=WEI+RHO(JJ)*(DREAL(G)**2+DIMAG(G)**2)
c      write(*,*)'wei rho(1)G dreal(G) dimag(G)',wei,rho(1),G,
c     + dreal(G),dimag(G)
      GOTO 8
      
 50   WEI=DREAL(PSI12)**2+DIMAG(PSI12)**2
c      write(*,*)'WEI ---in VZ 50 before returne', WEI
      RETURN
 8    WEI=WEI/2
c      write(*,*)'WEI ---in VZ 8 before returne and END', WEI
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
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
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
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
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
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ETA/ETA
      COMMON/FSI_FFF/F12,F21
      COMPLEX*16 FF0,F12,F21
      C=DCOS(H)
      S=DSIN(H)
      F12=DCMPLX(C,-S)
      F21=DCMPLX(C,S)
      IF(ICH.EQ.0)RETURN
      RHOP=RHO+H
      RHOM=RHO-H
      F12=FF0(RHO,H)*F12
      F21=FF0(RHO,-H)*F21
      RETURN
      END
      
      FUNCTION FAS(RKS)
C-- FAS=F*ACHR
C---F is the confluent hypergeometric function at k*r >> 1
C-- ACHR=sqrt(ACH), where ACH is the Coulomb factor
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 FAS,EIDC,ZZ1
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
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
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ETA/ETA
      COMPLEX*16 ALF,ALF1,Z,S,A,FF0,FAS
      DATA ERR/1.D-5/
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
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
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
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
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
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,JJ,MSPIN
      COMMON/FSI_BP/B,P/FSI_DERIV/BPR,PPR/FSI_SHH/SH,CHH
      COMMON/FSI_DAK/DAK,HAC,IFUN
      HAC=.0D0
      IF(ICH.EQ.0)GOTO 1
      XH=AC*AK
      HCP=HC(XH)
      HPR=HCP+.1544313298D0
      ACH=ACP(XH)
      HAC=2*HCP/AC
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

c-------------#include "gen/pilot.h"
      FUNCTION DFRSIN(X)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      

      DIMENSION A(0:16),B(0:15),C1(0:25),C2(0:28)

      PARAMETER (Z1 = 1, R8 = Z1/8, R32 = Z1/32)

      DATA C0 /1.25331 41373 15500 3D0/

      DATA NA,NB,NC1,NC2 /16,15,25,28/

      DATA A( 0) / 0.76435 13866 41860 002D0/
      DATA A( 1) /-0.43135 54754 76601 793D0/
      DATA A( 2) / 0.43288 19997 97266 531D0/
      DATA A( 3) /-0.26973 31033 83871 110D0/
      DATA A( 4) / 0.08416 04532 08769 354D0/
      DATA A( 5) /-0.01546 52448 44613 820D0/
      DATA A( 6) / 0.00187 85542 34398 220D0/
      DATA A( 7) /-0.00016 26497 76188 875D0/
      DATA A( 8) / 0.00001 05739 76563 833D0/
      DATA A( 9) /-0.00000 05360 93398 892D0/
      DATA A(10) / 0.00000 00218 16584 549D0/
      DATA A(11) /-0.00000 00007 29016 212D0/
      DATA A(12) / 0.00000 00000 20373 325D0/
      DATA A(13) /-0.00000 00000 00483 440D0/
      DATA A(14) / 0.00000 00000 00009 865D0/
      DATA A(15) /-0.00000 00000 00000 175D0/
      DATA A(16) / 0.00000 00000 00000 003D0/

      DATA B( 0) / 0.63041 40431 45705 392D0/
      DATA B( 1) /-0.42344 51140 57053 335D0/
      DATA B( 2) / 0.37617 17264 33436 566D0/
      DATA B( 3) /-0.16249 48915 45095 674D0/
      DATA B( 4) / 0.03822 25577 86330 087D0/
      DATA B( 5) /-0.00564 56347 71321 909D0/
      DATA B( 6) / 0.00057 45495 19768 974D0/
      DATA B( 7) /-0.00004 28707 15321 020D0/
      DATA B( 8) / 0.00000 24512 07499 233D0/
      DATA B( 9) /-0.00000 01109 88418 409D0/
      DATA B(10) / 0.00000 00040 82497 317D0/
      DATA B(11) /-0.00000 00001 24498 302D0/
      DATA B(12) / 0.00000 00000 03200 484D0/
      DATA B(13) /-0.00000 00000 00070 324D0/
      DATA B(14) / 0.00000 00000 00001 336D0/
      DATA B(15) /-0.00000 00000 00000 022D0/

      DATA C1( 0) / 0.99056 04793 73497 549D0/
      DATA C1( 1) /-0.01218 35098 31478 997D0/
      DATA C1( 2) /-0.00248 27428 23113 060D0/
      DATA C1( 3) / 0.00026 60949 52647 247D0/
      DATA C1( 4) /-0.00000 10790 68987 406D0/
      DATA C1( 5) /-0.00000 48836 81753 933D0/
      DATA C1( 6) / 0.00000 09990 55266 368D0/
      DATA C1( 7) /-0.00000 00750 92717 372D0/
      DATA C1( 8) /-0.00000 00190 79487 573D0/
      DATA C1( 9) / 0.00000 00090 90797 293D0/
      DATA C1(10) /-0.00000 00019 66236 033D0/
      DATA C1(11) / 0.00000 00001 64772 911D0/
      DATA C1(12) / 0.00000 00000 63079 714D0/
      DATA C1(13) /-0.00000 00000 36432 219D0/
      DATA C1(14) / 0.00000 00000 10536 930D0/
      DATA C1(15) /-0.00000 00000 01716 438D0/
      DATA C1(16) /-0.00000 00000 00107 124D0/
      DATA C1(17) / 0.00000 00000 00204 099D0/
      DATA C1(18) /-0.00000 00000 00090 064D0/
      DATA C1(19) / 0.00000 00000 00025 506D0/
      DATA C1(20) /-0.00000 00000 00004 036D0/
      DATA C1(21) /-0.00000 00000 00000 570D0/
      DATA C1(22) / 0.00000 00000 00000 762D0/
      DATA C1(23) /-0.00000 00000 00000 363D0/
      DATA C1(24) / 0.00000 00000 00000 118D0/
      DATA C1(25) /-0.00000 00000 00000 025D0/

      DATA C2( 0) / 0.04655 77987 37516 4561D0/
      DATA C2( 1) / 0.04499 21302 01239 4140D0/
      DATA C2( 2) /-0.00175 42871 39651 4532D0/
      DATA C2( 3) /-0.00014 65340 02581 0678D0/
      DATA C2( 4) / 0.00003 91330 40863 0159D0/
      DATA C2( 5) /-0.00000 34932 28659 7731D0/
      DATA C2( 6) /-0.00000 03153 53003 2345D0/
      DATA C2( 7) / 0.00000 01876 58200 8529D0/
      DATA C2( 8) /-0.00000 00377 55280 4930D0/
      DATA C2( 9) / 0.00000 00026 65516 5010D0/
      DATA C2(10) / 0.00000 00010 88144 8122D0/
      DATA C2(11) /-0.00000 00005 35500 7671D0/
      DATA C2(12) / 0.00000 00001 31576 5447D0/
      DATA C2(13) /-0.00000 00000 15286 0881D0/
      DATA C2(14) /-0.00000 00000 03394 7646D0/
      DATA C2(15) / 0.00000 00000 02702 0267D0/
      DATA C2(16) /-0.00000 00000 00946 3142D0/
      DATA C2(17) / 0.00000 00000 00207 1565D0/
      DATA C2(18) /-0.00000 00000 00012 6931D0/
      DATA C2(19) /-0.00000 00000 00013 9756D0/
      DATA C2(20) / 0.00000 00000 00008 5929D0/
      DATA C2(21) /-0.00000 00000 00003 1070D0/
      DATA C2(22) / 0.00000 00000 00000 7515D0/
      DATA C2(23) /-0.00000 00000 00000 0648D0/
      DATA C2(24) /-0.00000 00000 00000 0522D0/
      DATA C2(25) / 0.00000 00000 00000 0386D0/
      DATA C2(26) /-0.00000 00000 00000 0165D0/
      DATA C2(27) / 0.00000 00000 00000 0050D0/
      DATA C2(28) /-0.00000 00000 00000 0009D0/

       V=ABS(X)
       IF(V .LT. 8) THEN
       Y=R8*V
       H=2*Y**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 4 I = NB,0,-1
       B0=B(I)+ALFA*B1-B2
       B2=B1
    4  B1=B0
       H=SQRT(V)*Y*(B0-B2)
      ELSE
       R=1/V
       H=10*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 5 I = NC1,0,-1
       B0=C1(I)+ALFA*B1-B2
       B2=B1
    5  B1=B0
       S=B0-H*B2
       B1=0
       B2=0
       DO 6 I = NC2,0,-1
       B0=C2(I)+ALFA*B1-B2
       B2=B1
    6  B1=B0
       H=C0-SQRT(R)*(S*COS(V)+(B0-H*B2)*SIN(V))
      END IF
      IF(X .LT. 0) H=-H
      GO TO 9

      ENTRY DFRCOS(X)

       V=ABS(X)
       IF(V .LT. 8) THEN
       H=R32*V**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = NA,0,-1
       B0=A(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=SQRT(V)*(B0-H*B2)
      ELSE
       R=1/V
       H=10*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 2 I = NC1,0,-1
       B0=C1(I)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       S=B0-H*B2
       B1=0
       B2=0
       DO 3 I = NC2,0,-1
       B0=C2(I)+ALFA*B1-B2
       B2=B1
    3  B1=B0
       H=C0-SQRT(R)*((B0-H*B2)*COS(V)-S*SIN(V))
      END IF
      IF(X .LT. 0) H=-H
    9 DFRSIN=H
      RETURN
      END

c--#include "gen/pilot.h"
      FUNCTION CGAMMA(Z)
 
      
      COMPLEX*8 CGAMMA
      COMPLEX*8 Z,U,V,F,H,S       
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT      
      PARAMETER (NAME = 'CGAMMA')

      DIMENSION C(0:15)

      PARAMETER (Z1 = 1, HF = Z1/2)

c--#if defined(CERNLIB_QF2C)
c--#include "gen/gcmpfun.inc"
c--#endif

      DATA PI /3.14159 26535 89793 24D0/
      DATA C1 /2.50662 82746 31000 50D0/

      DATA C( 0) / 41.62443 69164 39068D0/
      DATA C( 1) /-51.22424 10223 74774D0/
      DATA C( 2) / 11.33875 58134 88977D0/
      DATA C( 3) / -0.74773 26877 72388D0/
      DATA C( 4) /  0.00878 28774 93061D0/
      DATA C( 5) / -0.00000 18990 30264D0/
      DATA C( 6) /  0.00000 00019 46335D0/
      DATA C( 7) / -0.00000 00001 99345D0/
      DATA C( 8) /  0.00000 00000 08433D0/
      DATA C( 9) /  0.00000 00000 01486D0/
      DATA C(10) / -0.00000 00000 00806D0/
      DATA C(11) /  0.00000 00000 00293D0/
      DATA C(12) / -0.00000 00000 00102D0/
      DATA C(13) /  0.00000 00000 00037D0/
      DATA C(14) / -0.00000 00000 00014D0/
      DATA C(15) /  0.00000 00000 00006D0/

c----#if !defined(CERNLIB_QF2C)
c----#include "gen/gcmpfun.inc"
c----#endif

      COMPLEX*8 GREAL,GIMAG,XARG,YARG       

      COMPLEX*8  ZARG,GCONJG,GCMPLX                                                  
           GREAL( ZARG)=REAL( ZARG)                                                  
	   GIMAG( ZARG)=AIMAG(ZARG)                                                  
	   GCONJG(ZARG)=CONJG(ZARG)                                                  
c--           GCMPLX(XARG,YARG)= CMPLX(XARG,YARG)     

      U=Z
      X=U
      IF(GIMAG(U) .EQ. 0 .AND. -ABS(X) .EQ. INT(X)) THEN
       F=0
       H=0
       WRITE(ERRTXT,101) X
c-       CALL MTLPRT(NAME,'C305.1',ERRTXT)
      ELSE
       IF(X .GE. 1) THEN
        F=1
        V=U
       ELSEIF(X .GE. 0) THEN
        F=1/U
        V=1+U
       ELSE
        F=1
        V=1-U
       END IF
       H=1
       S=C(0)
       DO 1 K = 1,15
       H=((V-K)/(V+(K-1)))*H
    1  S=S+C(K)*H
       H=V+(4+HF)
       H=C1*EXP((V-HF)*LOG(H)-H)*S
       IF(X .LT. 0) H=PI/(SIN(PI*U)*H)
      ENDIF

c----#if !defined(CERNLIB_DOUBLE)
      CGAMMA=F*H
c----#endif

      RETURN
  101 FORMAT('ARGUMENT EQUALS NON-POSITIVE INTEGER = ',1P,E15.1)
      END


      SUBROUTINE FSIPN !(WEIF)
C  calculating particle-nucleus Coulomb Wave functions FFij
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  !part. momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, ! 4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
      COMMON/FSI_ICH1/ICH1
      COMMON/FSI_ETA/ETA
      COMMON/LEDWEIGHT/WEIF,WEI,WEIN,ITEST,IRANPOS
      COMMON/FSI_FFPN/FF12,FF21
      COMPLEX*16 FF1,FF12,FF21
      FF12=DCMPLX(1.D0,0.D0)
      FF21=DCMPLX(1.D0,0.D0)
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
 
C=======================================================
C
      SUBROUTINE FSIWF !(WEI)
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
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
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
      COMMON/LEDWEIGHT/WEIF,WEI,WEIN,ITEST,IRANPOS
      
      
cmlv
        IF(IRANPOS.EQ.0)THEN
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
c        WRITE(6,38)'RP ',RP,'X ',X,Y,Z,T
c        FORMAT(A7,E11.4,A7,4E11.4)
        ENDIF 

      CVK=(P12X*PPX+P12Y*PPY+P12Z*PPZ)/(P12*AK)
      V=P12/E12
 
      IF(ICH.EQ.0)GOTO 21
      XH=AC*AK
      ACH=ACP(XH)
      ACHR=DSQRT(ACH)
      ETA=0.D0
      IF(XH.NE.0.D0)ETA=1/XH
C---HCP, HPR needed (e.g. in GST) if ICH=1
      HCP=HC(XH)
      HPR=HCP+.1544313298D0
  21  CONTINUE
      MSP=MSPIN
      DO 30 JJ=1,MSP
      ISPIN=JJ
      IF(NS.NE.1)GOTO22
C---Calc. quantities for the square well potential;
C-- for LL > 5 the square well potential is not possible or available
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
      IF(LL.EQ.12.OR.LL.EQ.13)GOTO20             ! piN 
      IF(LL.EQ.8.OR.LL.EQ.9.OR.LL.EQ.18)GOTO20   ! Nd, dd
      IF(LL.EQ.14.OR.LL.EQ.17.OR.LL.EQ.23)GOTO27 ! K+K-, K-p, K0K0-b
      A1=RD(JJ)*FD(JJ)*AKS
      A2=1+.5D0*A1
      IF(ICH.EQ.1)A2=A2-2*HCP*FD(JJ)/AC
      AKF=AK*FD(JJ)
      IF(ICH.EQ.1)AKF=AKF*ACH
      C(JJ)=FD(JJ)/DCMPLX(A2,-AKF)
      GOTO30
 20   CONTINUE
C---Calc. scatt. ampl. C(JJ) for pipi, piN and Nd, dd
      JH=LL-7+2*JJ-2
      IF(LL.EQ.8.OR.LL.EQ.9)GPI2=GND(AKS,JH)
      IF(LL.EQ.18)GPI2=GDD(AKS,JJ)
      IF(LL.EQ.5.OR.LL.EQ.6.OR.LL.EQ.7)GPI2=GPIPI(AKS,2)
      IF(LL.EQ.12.OR.LL.EQ.13)GPI2=GPIN(AKS,2)
      C(JJ)=1.D0/DCMPLX(GPI2,-AK) !pi+pi+, nd, pd, pi+p, dd
      IF(LL.NE.5.AND.LL.NE.6.AND.LL.NE.13)GOTO27
      IF(LL.EQ.5.OR.LL.EQ.6)GPI1=GPIPI(AKS,1)
      IF(LL.EQ.13)GPI1=GPIN(AKS,1)
      IF(LL.EQ.5.OR.LL.EQ.13)
     c           C(JJ)=.6667D0/DCMPLX(GPI1,-AK)+.3333D0*C(JJ) !pi+pi-,pi-p
      IF(LL.EQ.6)C(JJ)=.3333D0/DCMPLX(GPI1,-AK)+.6667D0*C(JJ) !pi0pi0
 27   CONTINUE
C---Calc. K+K-, K0K0-b or K-p s-wave scatt. ampl.
      IF(LL.EQ.14.OR.LL.EQ.23)CALL CKKB
      IF(LL.EQ.17)C(JJ)=DCMPLX(3.29D0,3.55D0)
C---Calc. pi+pi-, pi+pi+, pd, pi+p, pi-p, K+K- or K-p s-wave scatt. ampl.
C-- divided by Coulomb penetration factor squared (if ICH=1)
      IF(ICH.EQ.0)GOTO 30
      AAK=ACH*AK
      HCP2=2*HCP/AC
      C(JJ)=1/(1/C(JJ)-HCP2+DCMPLX(0.D0,AK-AAK))
c      write(*,*)'c(jj)',c(jj)
 30   CONTINUE
C***********************************************************************
c      write(*,*)'before call vz in fsiwf wei',wei
      CALL VZ  !(WEI)
      RETURN
      END
