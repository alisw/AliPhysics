          
          SUBROUTINE fsiini
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
      COMMON/FSI_ACH/HPR,AC,ACH,ACHR,ISPIN,MSPIN
      COMMON/FSI_NS/LL,NS,ICH,ISI,IQS,I3C,I3S
      COMMON/FSI_FD/FD(10),RD(10)
      COMMON/FSI_C/C(10),AM,AMS,DM
      COMMON/FSI_CONS/PI,PI2,SPI,DR,W
           COMPLEX*16 C
      COMMON/FSI_AA/AA
      COMMON/FSI_AAPI/AAPI(20,2)/FSI_AAND/AAND(20,4)
      COMMON/FSI_AAPIN/AAPIN(20,2)
      COMMON/FSI_SW/RB(10),EB(10),BK(10),CDK(10),SDK(10),
     1              SBKRB(10),SDKK(10)
      COMMON/LEDWEIGHT/WEIF,WEI,WEIN,ITEST,IRANPOS 

c      Include 'common_fsi_poc.inc'
c      Include 'common_fsi_prf.inc'
c      Include 'common_fsi_spin.inc'
c      Include 'common_fsi_ach.inc'
c      Include 'common_fsi_ns.inc'
c      Include 'common_fsi_fd.inc'
c-mlv      Include 'common_fsi_c.inc'
c      Include 'common_fsi_cons.inc'
c      Include 'common_fsi_aa.inc'
c      Include 'common_fsi_aapi.inc'
c      Include 'common_fsi_aapin.inc'
c      Include 'common_fsi_sw.inc'
c      Include 'common_fsi_aand.inc'


      DIMENSION FDH(30,10),RDH(30,10),EBH(30,10),RBH(30,10)
      DIMENSION RHOH(30,10)
      DIMENSION AM1H(30),AM2H(30),C1H(30),C2H(30),MSPINH(30)
C============= declarations pour l'appel de READ_FILE()============
      CHARACTER*10 KEY
      CHARACTER*8  CH8
      INTEGER*4    INT4
      REAL*8       REAL8
      INTEGER*4    IERR
C
C--- mass of the first and second particle
      DATA AM1H/.93956563D0,.93827231D0,.93956563D0,3.72737978D0,
     C          .13957D0,.13498D0,.13957D0, .93956563D0, .93827231D0,
     C          4*.13957D0,4*.493677D0,
     C          2*1.87561339D0,2*2.80892165D0,2*.497672D0,
     C          1.87561339D0,3*.93827231D0,.93956563D0, 2*0.D0/
      DATA AM2H/.93956563D0,.93827231D0,.93827231D0,3.72737978D0,
     C          .13957D0,.13498D0,.13957D0, 2*1.87561339D0,
     C          2*.493677D0,2*.93827231D0,
     C          2*.493677D0,2*.93827231D0,
     C          1.87561339D0,3.72737978D0,2.80892165D0,3.72737978D0,
     C          2*.497672D0,2*2.80892165D0,3.72737978D0,
     C          2*1.115684D0,2*0.D0/
c--------|---------|---------|---------|---------|---------|---------|----------
C---  charge of the first and second particle
      DATA C1H/0.D0,1.D0,0.D0,2.D0, 1.D0,0.D0,1.D0,0.D0,1.D0,
     C         3*1.D0,-1.D0,3*1.D0,-1.D0,
     C         4*1.D0,2*0.D0,4*1.D0,0.D0, 2*0.D0/
      DATA C2H/0.D0,1.D0,1.D0,2.D0,-1.D0,0.D0,3*1.D0,
     C         -1.D0,3*1.D0,-1.D0,3*1.D0,
     C         1.D0,2.D0,1.D0,2.D0,2*0.D0,2*1.D0,2.D0,2*0.D0,2*0.D0/
C---MSPIN vs (LL)
      DATA MSPINH/3*2,4*1,2*2,8*1,3,1,2,1,2*1,2*2,1,2*2, 2*0/
C---Spin factors RHO vs (LL,ISPIN)
      DATA RHOH/3*.25D0, 4*1.D0, 2*.3333D0, 8*1.D0, 
     1          .1111D0,1.D0,.25D0,1.D0,2*1.D0,
     1          .3333D0,.25D0,1.D0,2*.25D0, 2*0.D0,
     2          3*.75D0, 4*0.D0, 2*.6667D0, 8*0.D0, 
     2          .3333D0,.0D0,.75D0,.0D0,2*0.D0,
     2          .6667D0,.75D0,0.D0,2*.75D0, 2*0.D0, 
     3          17*.0D0,.5556D0,3*0.D0, 7*0.D0,2*0.D0,210*0.D0/
C---Scattering length FD and effective radius RD in fm vs (LL,ISPIN)
      DATA FDH/17.0D0,7.8D0,23.7D0,2230.1218D0,.225D0,.081D0,-.063D0,
     1        -.65D0,-2.73D0,
     1        .137D0,-.071D0,-.148D0,.112D0,2*1.D-6,-.360D0,
     1        2*1.D-6,1.344D0,6*1.D-6,-5.628D0,2.18D0,2.40D0, 2*0.D0,
cc     2 -10.8D0,2*-5.4D0,4*0.D0,-6.35D0,-11.88D0,8*0.D0,9*0.D0,
     2        3*-5.4D0,4*0.D0,-6.35D0,-11.88D0,8*0.D0,9*0.D0,
     2        1.93D0,1.84D0,2*0.D0,
     3        240*0.D0/
c--------|---------|---------|---------|---------|---------|---------|----------     
      DATA RDH/2.7D0,2.8D0,2.7D0,1.12139906D0,-44.36D0,64.0D0,784.9D0,
     1  477.9D0, 2.27D0, 9*0.D0,-69.973D0, 6*0.D0,3.529D0,
     1  3.19D0,3.15D0, 2*0.D0,
     2  3*1.7D0,4*0.D0,2.0D0,2.63D0, 17*0.D0,3.35D0,3.37D0, 2*0.D0, 
     3  240*0.D0/
C---Corresponding square well parameters RB (width in fm) and
C-- EB =SQRT(-AM*U) (in GeV/c); U is the well height
      DATA RBH/2.545739D0,   2.779789D0, 2.585795D0, 5.023544D0,
     1 .124673D0, .3925180D0,.09D0, 2.D0, 4.058058D0, 17*0.D0, 
     1  2.252623D0, 2.278575D0, 2*0.D0,
     2         3*2.003144D0,
     2         4*0.D0, 2.D0, 4.132163D0, 17*0.D0, 
     2  2.272703D0, 2.256355D0, 2*0.D0, 
     3         240*0.D0/
      DATA EBH/.1149517D0,    .1046257D0,   .1148757D0, .1186010D0,
     1    .7947389D0,2.281208D0,8.7D0,.4D0,.1561219D0,17*0.D0,
     1    .1013293D0, .1020966D0, 2*0.D0,
     2         3*.1847221D0,
     2         4*0.D0, .4D0, .1150687D0, 17*0.D0, 
     2         .09736083D0, .09708310D0, 2*0.D0, 
     3         240*0.D0/
C=======< constants >========================
      W=1/.1973D0    ! from fm to 1/GeV
      PI=4*DATAN(1.D0)
      PI2=2*PI
      SPI=DSQRT(PI)
      DR=180.D0/PI   ! from radian to degree

c    WRITE(*,*)'from C++ to fortran W PI PI2 SPI DR',W,PI,PI2,SPI,DR

        AC1=1.D10
        AC2=1.D10
C=======< condition de calculs >=============
      NUNIT=11 ! for IBM or HP
C      NUNIT=4 ! for SUN in Prague
c-mlv      CALL readint4(NUNIT,'ITEST     ',ITEST)      
c-mlv      CALL readint4(NUNIT,'LL        ',LL)        ! Two-particle system
c-mlv      CALL readint4(NUNIT,'NS        ',NS) 
c      CALL READ_FILE(NUNIT,'ITEST     ',CHAR,ITEST,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'LL        ',CHAR,LL,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'NS        ',CHAR,NS,REAL8,IERR)


 
C---setting particle masses and charges
      AM1=AM1H(LL)
      AM2=AM2H(LL)
      C1=C1H(LL)
      C2=C2H(LL)
 
C-- Switches:
C   ISI=1(0)  the strong interaction between the two particles ON (OFF)
C   IQS=1(0)  the quantum statistics ON (OFF);
C             should be OFF for nonidentical particles
C   I3C=1(0)  the Coulomb interaction with the nucleus ON (OFF)
C   I3S=1(0)  the strong interaction with the nucleus ON (OFF)
C   ICH=1(0)  if C1*C2 is different from 0 (is equal to 0)
C-  to switch off the Coulomb force between the two particles
C   put ICH=0 and substitute the strong amplitude parameters by
C   the ones not affected by Coulomb interaction
 
       IF(ITEST.EQ.0)THEN
       ICH=0
       IF(C1*C2.NE.0.D0) ICH=1
       IQS=0
       IF(C1+AM1.EQ.C2+AM2) IQS=1
       I3S=0  ! only this option is available
       ISI=1
       I3C=0
       IF(CN*ICH.NE.0.D0) I3C=1
       ENDIF

       
c      IF(ITEST.EQ.1)THEN
c      NS=4 ! SPHER. WAVE
c      ICH=0
c      IQS=1
c      ISI=0
c      I3C=0
      
      
c-mlv      CALL readint4(NUNIT,'ICH       ',ICH)
c-mlv      CALL readint4(NUNIT,'IQS       ',IQS)
c-mlv      CALL readint4(NUNIT,'ISI       ',ISI)
c-mlv      CALL readint4(NUNIT,'I3C       ',I3C)
c      CALL READ_FILE(NUNIT,'ICH     ',CHAR,ICH,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'IQS     ',CHAR,IQS,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'ISI     ',CHAR,ISI,REAL8,IERR)
c      CALL READ_FILE(NUNIT,'I3C     ',CHAR,I3C,REAL8,IERR)
c      ENDIF
      
      write(*,*)'====itest ll ich iqs isi i3c===',itest,ll,ich, iqs, isi, i3c
      
C==================================================================
C---fm to 1/GeV
      DO 3 J1=1,30
      DO 3 J2=1,10
      FDH(J1,J2)=FDH(J1,J2)*W
      RDH(J1,J2)=RDH(J1,J2)*W
 3    RBH(J1,J2)=RBH(J1,J2)*W
C---calcul. twice the reduced mass (AM), the relative mass difference
C-- (DM) and the Bohr radius (AC)
      AM=2*AM1*AM2/(AM1+AM2)
      AMS=AM*AM
      DM=(AM1-AM2)/(AM1+AM2)
      AC=1.D10
      C12=C1*C2
      IF(C12.NE.0.D0)AC=2*137.036D0/(C12*AM)
C---Setting spin factors
      MSPIN=MSPINH(LL)
      MSP=MSPIN
      DO 91 ISPIN=1,10
  91  RHO(ISPIN)=RHOH(LL,ISPIN)
C---Integration limit AA in the spherical wave approximation
      AA=0.D0
cc      IF(NS.EQ.2.OR.NS.EQ.4)AA=.5D0 !!in 1/GeV --> 0.1 fm
      IF(NS.EQ.2.OR.NS.EQ.4)AA=6.D0 !!in 1/GeV --> 1.2 fm
C---Setting scatt. length (FD), eff. radius (RD) and, if possible,
C-- also the corresp. square well parameters (EB, RB)
      DO 55 JJ=1,MSP
      ISPIN=JJ
      FD(JJ)=FDH(LL,JJ)
      RD(JJ)=RDH(LL,JJ)
      EB(JJ)=EBH(LL,JJ)
      RB(JJ)=RBH(LL,JJ)
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
      IF(LL.EQ.12.OR.LL.EQ.13)THEN
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
      ENDIF
C---Calculation continues for any system (any LL)
 55   CONTINUE
      RETURN
      END


      
       FUNCTION GPIPI(X,J)                                                                    
C--- GPIPI = k*COTG(DELTA), X=k^2                                                            
C--  J=1(2) corresponds to isospin=0(2)                                                      
       IMPLICIT REAL*8 (A-H,O-Z)                                                              
c--       Include 'common_fsi_aapi.inc'                                                          
c--       Include 'common_fsi_c.inc'                                                             
       COMMON/FSI_AAPI/AAPI(20,2)                                                         
       COMMON/FSI_C/HELP(20),AM,AMS,DM                                                    
       OM=DSQRT(X+AMS)                                                                        
       XX=X/AMS                                                                               
       GPIPI=OM/AAPI(1,J)                                                                     
       GPIPI=GPIPI*(1+(AAPI(3,J)-AAPI(1,J)**2)*XX+AAPI(4,J)*XX*XX)                            
       GPIPI=GPIPI/(1+(AAPI(3,J)+AAPI(2,J)/AAPI(1,J))*XX)                                     
                                                                                       
       END   
      
       FUNCTION GPIN(X,J)                                                                    
C--- GPIN = k*COTG(DELTA), X=k^2                                                            
C--  J=1(2) corresponds to piN isospin=1/2(3/2)                                             
       IMPLICIT REAL*8 (A-H,O-Z)                                                             
c--       Include 'common_fsi_aapin.inc'                                                        
       COMMON/FSI_AAPIN/AAPIN(20,2)                                                      
       GPIN=1/AAPIN(1,J)+.5D0*AAPIN(2,J)*X                                                   
                                                                                       
       END    
