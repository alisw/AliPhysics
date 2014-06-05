C==============================================================================
C  common blocks and constants 
C==============================================================================

      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT

      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG

      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META

      DOUBLE PRECISION        PI   
      DATA                    PI /3.141592653589793238462643D0/

      COMPLEX                 i
      DATA                    i/(0.,1.)/

      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3

      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
 
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2

      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI

      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL

      COMMON/FF0SCKPI/        F00
      REAL                    F00

      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0

      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal 
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma

      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal 
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma