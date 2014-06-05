      SUBROUTINE CH3PISET(JJ)
C information on 3 pion sub-channel under construction obtained
C J3PI=1 means 3 prong 
C J3PI=2 means 1 prong
C to be initialized in routine DPHSAA of tauola.f
      COMMON /CHANOPT/ J3PI
      INTEGER          J3PI
      J3PI=JJ
      end
     
      SUBROUTINE CH3PIGET(JJ)
C information on 3 pion sub-channel under construction obtained
C J3PI=1 means 3 prong 
C J3PI=2 means 1 prong
C to be initialized in routine DPHSAA of tauola.f
      COMMON /CHANOPT/ J3PI
      INTEGER          J3PI
      IF (J3PI.EQ.1.OR.J3PI.EQ.2) THEN
        JJ=J3PI
      ELSE
       write(*,*) 'FROM value_parameter.f CH3PIGET, wrong J3PI=',J3PI
       stop
      ENDIF
      end

      SUBROUTINE GETFF2PIRHO(JJ)
      IMPLICIT NONE
      include '../parameter.inc'
      INTEGER JJ
      JJ = FF2PIRHO
      END

      SUBROUTINE OLACHNL(SIGN)
C provides sign of tau, to be used in CP dependent parts of current.
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      INTEGER            JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / IDFC  / IDFF
      INTEGER KTO
      REAL    SIGN
      IF     (KTOM.EQ.1.OR.KTOM.EQ.-1) THEN
        SIGN= IDFF/ABS(IDFF)
      ELSEIF (KTOM.EQ.2) THEN
        SIGN=-IDFF/ABS(IDFF)
      ELSE
        PRINT *, 'STOP IN OLACHNL: KTOM=',KTOM
        STOP
      ENDIF
      END

      FUNCTION COEFrr(I,J)
C clebsh gordan (or so ...)  coefs for 3 scalar final states
      implicit none
C TAUOLA RChL COEF(I,J) =  COEFr(I,J)
      REAL COEFr(1:5,0:7)
      REAL COEFrr
      DATA PI /3.141592653589793238462643/
      REAL PI
      DATA ICONT /0/
      INTEGER ICONT
      INTEGER I,J
      REAL FPIr

C initialization of FPI matrix defined in ...
C FPIc is to be used with cleo initialization
C FPIr is to be used with RChL initialization
C actual choice is made in ???


      DATA  FPIr /92.4E-3/


C initialization of COEF matrix defined in ...
C COEFc is to be used with cleo initialization
C COEFr is to be used with RChL initialization
      IF (ICONT.EQ.0) THEN
       ICONT=1
C
C********* COEFr(I,J) *******

       COEFr(1,0)= 1.
       COEFr(2,0)= -1.
       COEFr(3,0)= 0.
       COEFr(4,0)= 1.
       COEFr(5,0)= 0.

       COEFr(1,1)= 1.
       COEFr(2,1)= -1.
       COEFr(3,1)= 0.
       COEFr(4,1)= 1.
       COEFr(5,1)= 1.
C
       COEFr(1,2)=1.
       COEFr(2,2)= -1.
       COEFr(3,2)= 0.0
       COEFr(4,2)= 1.
       COEFr(5,2)=1.
C
       COEFr(1,3)= 0.
       COEFr(2,3)= 1.
       COEFr(3,3)= -1.
       COEFr(4,3)= 1.
       COEFr(5,3)= - 1.
C
       COEFr(1,4)= 1.0/SQRT(2.)/3.0
       COEFr(2,4)=-1.0/SQRT(2.)/3.0
       COEFr(3,4)= 0.0
       COEFr(4,4)= 0.0
       COEFr(5,4)= 0.0
C
       COEFr(1,5)=-SQRT(2.)/3.0
       COEFr(2,5)= SQRT(2.)/3.0
       COEFr(3,5)= 0.0
       COEFr(4,5)= 0.0
       COEFr(5,5)=-SQRT(2.)
C
       COEFr(1,6)= 1./3.
       COEFr(2,6)=-2./3.
       COEFr(3,6)= 2./3.
       COEFr(4,6)= 0.0
       COEFr(5,6)=-2.0
C
       COEFr(1,7)= 0.0
       COEFr(2,7)= 0.0
       COEFr(3,7)= 0.0
       COEFr(4,7)= 0.0
       COEFr(5,7)=-SQRT(2.0/3.0)
      ENDIF

      COEFrr=COEFr(I,J)
      END

      subroutine rchl_parameters(KAK)
      implicit none
C==============================================================================
C  Initialization, of '../parameter.inc' common block group 
C 
C  KAK may be equal to JAK of TAUOLA namespace, but it is not always the case 
C  Hard-coded  fit parameters:
C  rho, rhoprime, f2(1275), f0(1186), sigma(made up!)
C  The value of both the mass and width of resonances are taken 
C  from  fit to ALEPH data  (ref [1], Set 1) 
C  References: [1] arXiv: 0911.4436 [hep-ph] D.  Gomez Dumm et al
C              [2] arXiv: 0911.2640 [hep-ph] D.  Gomez Dumm et al. 
C              [3] P Roig, talk PhiPsi2011, Novosibirsk
C              [4] arXiv:0807.4883 [hep-ph] Diogo R. Boito et al.
C              [5] arXiv:0803.1786 [hep-ph] M. Jamin et al.
C  WARNING:    some of parameters require RERUN of da1wid_tot_rho1_gauss.f
C              pretabulating Q dependent a1 width,
C              directory RChL-currents/tabler/a1
C==============================================================================
      include '../parameter.inc'
      INTEGER KAK
      DATA IWARM/0/
      INTEGER IWARM
      INTEGER          J3PI
      COMMON /CHANOPT/ J3PI

         IF(KAK.EQ.4) THEN  
C  /MASS_RES/; resonances parameters initialization: 
C                           ! at present only for two pion mode non-default 
c                           ! values are used:
      mro = 0.77554d0
      mrho1 = 1.453d0
      grho1 = 0.50155D0
c   /PAR_RHOPRIME/; parameters of rho' and rho'' 
C                   used for 2 pion form factor, reference [3]
      COEF_GA =  0.14199D0    
      COEF_DE = -0.12623D0    
      phi_1   = -0.17377D0 
      phi_2   =  0.27632D0  
      grho2   =  0.41786D0  
      mrho2   =  1.8105d0
         ELSE IF(KAK.EQ.5) THEN 
      MRO   = 0.771849d0      !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MRHO1 = 1.35d0     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GRHO1 = 0.448379d0 !0.473287d0       !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f

      ELSE
      MRO   = 0.775     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MRHO1 = 1.465     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GRHO1 = 0.4       !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f

c   /PAR_RHOPRIME/; parameters of rho' and rho'' 
C                   used for 2 kaon form factor, reference [3]
c                   FOR THE MOMENT THEIR NUMERICAL VALUES COINCIDE WITH
c                   ONES FOR THE TWO PION MODE !!!!
      COEF_GA =  0.14199D0    
      COEF_DE = -0.12623D0    
      phi_1   = -0.17377D0 
      phi_2   =  0.27632D0  
      grho2   =  0.41786D0  
      mrho2   =  1.8105d0
         ENDIF


        IF(KAK.EQ.70)  THEN    ! non default values to be used
                               ! for KPI MODE NO FSR INTERACTION
c   /PAR_KPI/; parameters for Kpi mode, reference [4], table 4, row2
      MKST        = 0.943d0
      MKSTPR      = 1.374D0
      GAMMA_KST   = 0.06672d0  
      GAMMA_KSTPR = 0.240d0
      GAMMA_RCHT  =-0.039d0
        ELSE IF(KAK.EQ.71)  THEN     ! non default values to be used 
                                     ! for KPI MODE WITH FSR INTERACTION
c   parameters for Kpi mode, reference [5]
      MKST = 0.8953d0
      GAMMA_KST = 0.0475d0
      MKSTPR = 1.307d0
      GAMMA_KSTPR = 0.206d0
      GAMMA_RCHT = -0.043d0
        ELSE
C   /MASS_SCAL/; stable particles - final scalars 
      Mksp  = 0.89166d0   !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      Mks0  = 0.89610d0   !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MKST = (Mksp +Mks0)/2.
      MKSTPR = 1.374d0    !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GAMMA_KST   = 0.06672  
      GAMMA_KSTPR = 0.240
      GAMMA_RCHT  = -0.043 !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
        ENDIF

C  /RCHT_3PI/; model parameters; their value are from fit,
c             reference [1], set 1
C              CHANGE OF THEIR VALUES REQUIRES 
C              RERUN /tabler/a1/da1wid_tot_rho1_gauss.f
 
      IF(KAK.EQ.5) THEN
      FPI_RPT  = 0.091337d0
      FV_RPT   = 0.168652d0   
      FA_RPT   = 0.131425d0                 
      BETA_RHO = -0.318551d0 
      ELSE
      FPI_RPT  = 0.0924
      FV_RPT   = 0.18  
      FA_RPT   = 0.149                
      BETA_RHO =  -0.25  
      ENDIF
 
      FK_RPT = FPI_RPT*1.198d0  
      GV_RPT   = FPI_RPT*FPI_RPT/FV_RPT

c$$$c   It has to be used for a new parametrization of rho1 for 3pions, 
C$$$c   that is not checked yet
c$$$c        IF(KAK.EQ.5)  THEN  ! high energy behaviour imposes these relations   
c$$$c      GV_RPT   = 0.066  
c$$$c      FV1_RPT = 0.18D0
c$$$c      GV1_RPT = (FPI_RPT*FPI_RPT- FV_RPT*GV_RPT)/FV1_RPT
c$$$c        ELSE
c$$$c      GV_RPT   = FPI_RPT*FPI_RPT/FV_RPT
c$$$c      ENDIF

c   /SCAL_3PI/; parameters of sigma meson for 3 pion mode
C* Parameteres for the sigma contribution, using BW for sigma
      IF(KAK.EQ.5)  THEN
          IF (J3PI.EQ.1) THEN
      alpsig = -8.795938d0 
      betasig = 9.763701d0 
      gamsig =  1.264263d0 
      delsig =  0.656762d0 
      rsigma =  1.866913d0 
           ELSE IF (J3PI.EQ.2) THEN
      alpsig = 1.139486d0*0.63d0
      betasig = 1.139486d0*0.63d0
      gamsig =  0.889769d0*0.63d0
      delsig =  0.889769d0*0.63d0 
      rsigma =  0.000013d0
           ENDIF
       ENDIF

C  /MASS_RES/ 
      IF(KAK.EQ.5)  THEN
      MMA1  = 1.091865d0      !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
          IF (J3PI.EQ.1) THEN
      MSIG  = 0.487512d0 
      GSIG  = 0.70d0     
          ELSE IF(J3PI.EQ.2) THEN
      MSIG  = 0.55d0 
      GSIG  = 0.7d0 
          ENDIF
      ELSE
      MMA1 = 1.12
      MSIG = 0.475
      GSIG = 0.550
      ENDIF
        call rchl_REparam(0,IWARM,KAK)     
         IF (IWARM.EQ.1) RETURN  ! parameters below do not need 
         IWARM=1                 ! re-initialization

C  /MASS_RES/ 
      GRO   = 0.149d0     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MF2   = 1.275d0
      GF2   = 0.185d0
      MF0   = 1.186d0
      GF0   = 0.350d0
      MSG   = 0.860d0
      GSG   = 0.880d0
      MPHI  = 1.019d0     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GPHI  = 0.0042d0    !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MOM   = 0.781940d0  !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GOM   = 0.00843d0   !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f

C    /RES_MIXING_RCHT/; a parameter defines w-phi angle mixing
      THETA = 35.*PI/180. 

C    /FF0SCKPI/ a parameter normalized FFSC_KPI 
      F00 = 0.972


 

C   /MASS_SCAL/; stable particles - final scalars 
C              CHANGE OF THEIR VALUES (useful for some tests) REQUIRES, 
C              RERUN /tabler/a1/da1wid_tot_rho1_gauss.f

      MPIZ    = 0.1349766d0           !PKORB(1,7)    ! NEUTRAL PION MASS
      MPIC    = 0.13957018d0          !PKORB(1,8)    ! CHARGED PION MASS
      MMPI_AV = (MPIZ+2.*MPIC)/3.d0
      MKZ     = 0.497648d0            !PKORB(1,12)    ! NEUTRAL KAON MASS
      MKC     = 0.493677d0            !PKORB(1,11)    ! CHARGED KAON MASS
      MMK     = (MKC+MKZ)/2.d0
      MTAU    = 1.777
      MNUTA   = 0.001
      META    = 0.547d0
 


c   /PAR_KKPI/; parameters to describe KKpi modes, reference [2]
C              CHANGE OF THEIR VALUES REQUIRES 
C              RERUN /tabler/a1/da1wid_tot_rho1_gauss.f

      G2    =  mro/(192.*pi*pi*sqrt(2.)*FV_RPT)*3.
      G13   = -2.*g2 
      G4    = -0.72
      G5    = -0.6-2.*g4 
      C125  =  0.
      C1256 = -3/96./pi**2*FV_RPT*MRO/SQRT(2.)/FPI_RPT**2 
      C1235 =  0.
      C4    = -0.07
      D123  =  0.05 
      D3    = -MRO**2/(64.*PI*PI*FPI_RPT**2) 


c   /PAR_KPI/; parameters to describe Kpi mode, reference [4]
      Ht0 = -1.2400398216503017D-2
C     Ht0 = !!!!! TO ADD A FORMULAE FOR Ht0 (Jamin's email) !!!!!
      lap_KPI = 24.66e-3
      lapp_KPI = 11.99e-4
      c1_KPI = lap_KPI/mpic**2
      c2_KPI = (lapp_KPI - lap_kpi**2)/2.d0/mpic**4

c /KPISC_EM/; parameters for Kpi scalar FF from 
c  http://arxiv.org/pdf/1103.4855.pdf
      lnC = 0.20193d0
      lambda0 = 0.013139d0

c   /SCAL_3PI/; parameters of sigma meson for 3 pion mode
      a00_3piscal = 0.220
      b00_3piscal = 0.268/mmpi_av**2
      c00_3piscal = -0.0139/mmpi_av**4
      d00_3piscal = -0.00139/mmpi_av**6
      x00_3piscal = 36.77*mmpi_av**2
      a02_3piscal = -0.0444
      b02_3piscal = -0.0857/mmpi_av**2
      c02_3piscal = -0.00221/mmpi_av**4
      d02_3piscal = -0.000129/mmpi_av**6
      x02_3piscal = -21.62*mmpi_av**2
      MMF0 = 0.441

c /SCAL_3PI/; parameters for the scalar part 3 pion modes
c   Pablo private
      ALPHA0_3PI = 1.
      ALPHA1_3PI = 1.
      GAMMA0_3PI = 1.
      GAMMA1_3PI = 1.


C FFVEC: dipswitch for Final State interaction in two scalar modes
C     with FSI (default FFVEC =1) and
C     without FSI (FFVEC =0)

      FFVEC = 1

C  FFKPIVEC : parameter to choose the parametrization for 
C             vector Kpi form factor with FSI effects
C    FFKPIVEC = 0 parametrization Eqs.(17),(18) of [4]
C    FFKPIVEC = 1 parametrization Eq.(5) of [5]
C    FFKPIVEC = 2 parmetrization [4], total result
      FFKPIVEC = 2
C  FFKPISCAL : parameter to choose the parametrization for 
C             scalar Kpi form factor with FSI effects
C  FFKPISCAL = 0 no scalar contribution
C  FFKPISCAL = 1 parametrization of Mathias Jamin,adopted his private code
C  FFKPISCAL = 2 parametrization of Emilie Passerman,
C               adopted her private code []
      FFKPISCAL = 1 

C FFKKVEC: dipswitch for K0K- mode
C     with rho' and rho'' (FFKKVEC =1) and
C     without rho' and rho''  (default FFKKVEC =0)
      FFKKVEC = 0

C FF3PISCAL: dipswitch for the scalar contribution for 3 pion modes
C     with the scalar contribution ( default FF3PISCAL = 2) 
c     FF3PISCAL = 2 BW parametrization for sigma meson
c     FF3PISCAL = 1 simplified RCHT results
C     FF3PISCAL =0  no sigma contribution
      FF3PISCAL = 2

C  Implemetation of another parametrization rho1, not checked yet by tests 
C FF3PIRHOPR: dipswitch for the parametrization for rho' contribution 
C     For 3 pion modes
C     general parametrization  ( default FF3RHOPR =1) and
C     simplified  (FF3PIRHOPR =0)
      FF3PIRHOPR = 0


C FF2PIRHO: dipswitch for the two pion form factor (default FF2PIRHO = 1) 
C FF2PIRHO =1   RChL parametrization 
C FF2PIRHO = 2  Belle parametrization, 
C               all parameters par (1...11) of fit are free 
C FF2PIRHO = 3  Belle parametrization, 
C               parameters of fit are free 
C               except for fixed par(1)=F_pi(0)=1
      FF2PIRHO =2


C FCOUL: dipswitch for the Coulomb interaction 
C FCOUL = 1  with 
C FCOUL = 0  without
      FCOUL = 0

      call rchl_REparam(1,IWARM,KAK)
      return
      end

      subroutine rchl_REparam(IMODE,IWARM,KAK)
      include '../parameter.inc'
      common / PARAMS / P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,IUSE
      INTEGER                                                                  IUSE
      DOUBLE PRECISION  P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16
      DATA IUSE /0/

      IF(IUSE.EQ.0) RETURN
      
      IF (IMODE.EQ.-1) THEN
        IWARM=IWARM
      ELSE

C FF3PISCAL: dipswitch for the scalar contribution for 3 pion modes
C     with the scalar contribution ( default FF3PISCAL = 2) 
c     FF3PISCAL = 2 BW parametrization for sigma meson
c     FF3PISCAL = 1 simplified RCHT results
C     FF3PISCAL =0  no sigma contribution
c      FF3PISCAL = 2

C     CANDIDATES FOR PARAMETERS TO FIT with default values
C* Parameteres for the sigma contribution, using BW for sigma
      alpsig  = P1
      betasig = P2
      gamsig  = P3
      delsig  = P4
      rsigma  = P5


      MRO   = P6     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MRHO1 = P7     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GRHO1 = P8       !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
C  /MASS_RES/ 
      GRO   = P9     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MMA1  = P10    !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MSIG  = P11
      GSIG  = P12

C  /RCHT_3PI/; model parameters; their value are from fit,
c             reference [1], set 1
C              CHANGE OF THEIR VALUES REQUIRES 
C              RERUN /tabler/a1/da1wid_tot_rho1_gauss.f
      FPI_RPT  = P13
      FV_RPT   = P14
      FA_RPT   = P15
      BETA_RHO = P16
      FK_RPT = FPI_RPT*1.198d0  
      GV_RPT   = FPI_RPT*FPI_RPT/FV_RPT

      ENDIF
      return
      end
