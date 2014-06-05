C     JAK = 5  
C
      COMPLEX FUNCTION F3PI_RCHT(IFORM,QQ,SA,SB)
      IMPLICIT NONE
      INTEGER                    IFORM
      REAL                             QQ,SA,SB
C.......................................................................
C.
C.    F3PI - RchT version of the hadronic curent used in TAUOLA
C.    References: [1] arXiv:0911.4436 (hep-ph) P. Roig et al.
C.                eqs (5)-(9), (32) gives the main part of the current
C                 (the part without the sigma contribution) 
C.                [2] arXiv:1203.3955, 
C.                the manual for the 3pi current without the sigma contribution
C.                [3] http://annapurna.ifj.edu.pl/~wasm/RChL/RChLa1.pdf
C                     eq (3) the sigma meson contribution to the 3 pi current
C.    Inputs    : QQ,SA,SB - invariant masses**2  [GeV**2]
C.              : IFORM formfactor no.
C.    Outputs   : F3PI_RCHT formfactor value
C.
C.    COMMON    : RCHT_3PI content is defined in this routine
C.
C.    Calls     : functions from file ./funct_rpt.f,  ./fa1rchl.f
C.    Called    : from file ../formf.f
C.    Author    : O.S
C.    Created   :
C.    Modified  : 1. Feb 2011
C                 a part with scalars added on the 1st May 2012
C.......................................................................
       include '../parameter.inc'
       include '../funct_declar.inc'

      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS    
      INTEGER          J3PI  
      REAL wid_a1_fit 
*
      CHARACTER*(*) CRNAME
      PARAMETER(    CRNAME = 'F3PI' )
*
      INTEGER IDK
      DOUBLE PRECISION M1,M2,M3,M1SQ,M2SQ,M3SQ
      REAL GGMA1
      REAL S1,S2,S3,FACT_RPT
*
      REAL R_RPT,KAP_RPT,FACT_ADD_RPT
      COMPLEX ALP21_RPT,ALP22_RPT,ALP11_RPT,ALP12_RPT
     $       ,BETA1_RPT,BETA2_RPT,BETA1_RPT_RHO1,ALP11_RPT_RHO1
     $       ,ALP21_RPT_RHO1,ALP22_RPT_RHO1
      COMPLEX FA1RCHL
      INTEGER IFIRST
      DATA IFIRST/0/

C. GENERAL INITIALIZATION
C. ======================

      IF (IFIRST.EQ.0) THEN
        IFIRST = 1     
        PRINT *,' In F3pi: chpt + 1 resonance + 2 resonance (RchT)'
      END IF

C******************************************
C    Initilisation of the mass of the particles
C*****************************************
        call rchl_parameters(5)

C. We impose isospin symmetry requesting that charged and neutral pion mass
C. are equal. This may need to be changed MMPI_AV = (2.*MPIC+MPIZ)/3.

	M1 = MMPI_AV      
        M2 = MMPI_AV     
        M3 = MMPI_AV     


C. normalization factor to compensate on different  
C. convention for normalization constant  used in  TAUOLA and 
C. TAUOLA documentation on one side and paper [1] on other.
        FACT_ADD_RPT = 1./FPI_RPT



C. FUNCTION VALUE, GENERAL CASE
C. ============================

C     Function value set to 0
      F3PI_RCHT = CMPLX(0.,0.)

C.    First determine whether we are doing pi-2pi0 or 3pi.
C. we  CH3PI GET information (eg from phase space geneator of tauola.f)
C. whether it is 3 prong J3pi=1 or 1 prong J3pi=2 final state of 3 pion.
      CALL CH3PIGET(J3pi)

      IF (J3pi.EQ.2) THEN  ! it is pi 2pi0
          IDK     = 1
	  R_RPT   = -1. 
          KAP_RPT = 0.5 
      ELSE IF(J3pi.EQ.1) THEN              !  it is 3pi
          IDK     = 2
	  R_RPT   = 1.
	  KAP_RPT = 1.
      END IF

      M1SQ = M1*M1
      M2SQ = M2*M2
      M3SQ = M3*M3


C. Calculation, IFORM = 1 or 2.
C.   VECTOR 3 PION FORM FACTORS
C. ============================ 
      IF (IFORM.EQ.1.OR.IFORM.EQ.2) THEN
        S1 = SA     ! t variable in  [2] (!!! vec1 = v2 !!!)
        S2 = SB    ! s variable in [2]
        S3 = QQ-SA-SB+M1SQ+M2SQ+M3SQ

        IF (S3.LE.0..OR.S2.LE.0.) RETURN

C.    FUNCTIONS BETA_RPT, ALP1_RPT are coded in ./funct_rpt.f
C.    they are defined in Eq (6) of [2]
        BETA1_RPT = BETA_RPT(QQ,S1,S2,M1SQ,M2SQ,M3SQ,MRO,MRHO1,GRHO1)
	ALP11_RPT = ALP1_RPT(QQ,S1,S2,M1SQ,M2SQ,M3SQ,MRO,MRHO1,GRHO1)

        F3pi_RCHT = - 2.*SQRT(2.)/(3.*FPI_RPT) -
     $	   SQRT(2.)*FV_RPT*GV_RPT/(3.*FPI_RPT**3)*ALP11_RPT +
     $     4.*FA_RPT*GV_RPT/(3.*FPI_RPT**3)*BETA1_RPT*QQ
     $     *FA1RCHL(QQ)

C.    FA1RCHL(QQ) is the a1 propagator, is it coded in ./fa1rchl.f

     
        F3pi_RCHT =  F3pi_RCHT * R_RPT  
   
C.
C.      The contribution from the scalar (sigma) resonance
C.      
      IF(FF3PISCAL.EQ.1) THEN
C.      This parametrization does not fit
         IF(J3PI.EQ.1) THEN
        F3pi_RCHT = F3pi_RCHT + (
     &            sqrt(2.)*(R0scal_3pi(QQ,S1)+R0scal_3pi(QQ,S2)) +
     &            (R2scal_3pi(QQ,S1)+R2scal_3pi(QQ,S2))
     &            )
         ELSE IF(J3PI.EQ.2) THEN
        F3pi_RCHT = F3pi_RCHT - (
     &            (R0scal_3pi(QQ,S1)+R0scal_3pi(QQ,S2)) -
     &            sqrt(2.)*(R2scal_3pi(QQ,S1)+R2scal_3pi(QQ,S2))
     &            )
         ENDIF
       ELSE IF(FF3PISCAL.EQ.2) THEN
C.    A new parametrization 10.02.2013 for the scalar contribution
C.       analytical formulae in [3], eqs (3)-(6), 
C.       functions BWsig(Mm,Gg,Qx), FFsig(QQ,Qx) coded in ./funct_rpt.f
C.
         IF(J3PI.EQ.1) THEN
        F3pi_RCHT = F3pi_RCHT + 
     &              SQRT(2.)*FV_RPT*GV_RPT/(3.*FPI_RPT**3)*
     &             (alpsig*BWsig(MSIG,GSIG,S1)*FFsig(QQ,S1) +
     &              betasig*BWsig(MSIG,GSIG,S2)*FFsig(QQ,S2))
     &             +4.*FA_RPT*GV_RPT/(3.*FPI_RPT**3)*
     &             (gamsig*BWsig(MSIG,GSIG,S1)*FFsig(QQ,S1) +
     &              delsig*BWsig(MSIG,GSIG,S2)*FFsig(QQ,S2))
     $     *FA1RCHL(QQ) 
        ELSE IF(J3PI.EQ.2) THEN
        F3pi_RCHT = F3pi_RCHT - (
     &              SQRT(2.)*FV_RPT*GV_RPT/(3.*FPI_RPT**3)*
     &             alpsig*BWsig(MSIG,GSIG,S3)*FFsig(QQ,S3) 
     &             +4.*FA_RPT*GV_RPT/(3.*FPI_RPT**3)*
     &             gamsig*BWsig(MSIG,GSIG,S3)*FFsig(QQ,S3) 
     $     *FA1RCHL(QQ) )
         ENDIF

C.
C.     The Coulomb interaction effects for the final pions
C.          only for pi-pi+pi- should be included
C.     Functions fattcoul(m1,m2,s),frepcoul(m1,m2,s) are coded in ./funct_rpt.f
C.          
      IF (FCOUL.EQ.1.AND.(J3pi.EQ.1)) THEN !  it is 3pi
        F3PI_RCHT = F3PI_RCHT*sqrt(fattcoul(m2,m3,s1)
     &      *fattcoul(m1,m3,s2)*frepcoul(m1,m2,s3))
      END IF        
       ENDIF


C.
C. The factor 1/FACT_ADD_RPT = FPI_RPT comes to compensate an additional factor
C.  in the hadronic current F3pi_rcht  above compare with eq (6) in [2]


       F3pi_RCHT = F3pi_RCHT/FACT_ADD_RPT 

C. Calculation, for  IFORM = 3 is not needed
C. ======================= 
C.   F3PI_RCHT is set to zero in ../formf.f. 

C. Calculation, IFORM = 4.
C. PSEUDOSCALAR 3 PION FORM FACTOR
C. ======================= 
      ELSEIF (IFORM.EQ.4) THEN
        S1 = SA
        S2 = SB
        S3 = QQ-SA-SB+M1SQ+M2SQ+M3SQ

        IF (S3.LE.0..OR.S2.LE.0.) RETURN

C. Functions  ALP21_RPT,  ALP22_RPT are Eq(10) in [2]
C. GRHO_RCHT, GRHO1_RCHT s1 or s2 dependent widths of rho or rho1
C. coded in ./funct_rpt.f 
        ALP21_RPT = 3.*GV_RPT*S1*M1SQ*(S3-S2)/
     &	  (FV_RPT*QQ*(QQ-M1SQ)*(1.+BETA_RHO))*
     &    (1./(S1-MRO**2+i*MRO*GRHO_RCHT(S1,MRO))+
     &     BETA_RHO/(S1-MRHO1**2+i*MRHO1*GRHO1_RCHT(S1,MRHO1,GRHO1)))
        ALP22_RPT = 3.*GV_RPT*S2*M1SQ*(S3-S1)/
     &	  (FV_RPT*QQ*(QQ-M1SQ)*(1.+BETA_RHO))*
     &    (1./(S2-MRO**2+i*MRO*GRHO_RCHT(S2,MRO))+
     &     BETA_RHO/(S2-MRHO1**2+i*MRHO1*GRHO1_RCHT(S2,MRHO1,GRHO1)))

C.  PSEUDOSCALAR 3 PION FORM FACTOR. Eqs (9) in [2]
        F3PI_RCHT = SQRT(2.)/(3.*FPI_RPT*QQ*(QQ-M1SQ))*
     &      M1SQ*(3.*(S3-M3SQ)-QQ*(1.+2.*KAP_RPT*R_RPT))

        F3PI_RCHT = F3PI_RCHT -  sqrt(2.)*FV_RPT*GV_RPT/(3.*FPI_RPT**3)*
     &	  ( ALP21_RPT + ALP22_RPT)

C.
C. The factor 1/FACT_ADD_RPT = FPI_RPT comes to compensate an additional factor
C.  in the hadronic current F3pi_rcht  above compare with eq (9) in [2]
	F3PI_RCHT = R_RPT * F3PI_RCHT/FACT_ADD_RPT

C.
C.     The Coulomb interaction effects for the final pions
C.          only for pi-pi+pi- should be included
C.     Functions fattcoul(m1,m2,s),frepcoul(m1,m2,s) are coded in ./funct_rpt.f
C. 
      IF (FCOUL.EQ.1.AND.(J3pi.EQ.1)) THEN !  it is 3pi
        F3PI_RCHT = F3PI_RCHT*sqrt(fattcoul(m2,m3,s1)
     &      *fattcoul(m1,m3,s2)*frepcoul(m1,m2,s3))
      END IF

      END IF

      RETURN
      END

