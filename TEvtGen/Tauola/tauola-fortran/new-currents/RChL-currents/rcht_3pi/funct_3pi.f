C **********************************************************
C library of functions used in calculation of currents 
C references:
C [1]   arXiv:0911.4436 (hep-ph)  D. Gomez Dumm et al. (tau -> 3pi nu) 
C [2]   arXiv:0911.2640 (hep-ph) D. Gomezz Dumm et al. (tau -> KKpi nu)
C [3]   arXiv:0807.4883 (hep-ph) D. R. Boito et al. (tau -> Kpi nu)
C [4]   P. Roig, talk at (tau -> 2 pi nu)
C [5]   arXiv:0803.2039  (hep-ph) E. Arganda et al., Appendix B (tau -> KK nu)
C **********************************************************

      FUNCTION ALP1_RPT(Q,S,T,M1SQ,M2SQ,M3SQ,
     $                           MMM,MMM2,GGG)
      IMPLICIT NONE
      REAL                      Q,S,T 
      DOUBLE PRECISION                M1SQ,M2SQ,M3SQ,
     &                           MMM,MMM2,GGG
C **********************************************************
C     COMPLEX FUNCTION ALP1_RPT;  
C     one-resonance contribution to 
C     F3pi_rcht(iform=1,2) at f3pi_rcht.f;
C     corrected  formula (8) of REF [1] including rho1, 
C     the factor sqrt(2)*FV_RPT*GV_RPT/(3.*FPI_RPT**3) is included in
C     F3pi_rcht(iform1=1,2) at f3pi_rcht.f 
C **********************************************************
      REAL U
      include '../parameter.inc'
      include '../funct_declar.inc'

       U= Q-S-T+M1SQ+M2SQ+M3SQ
       ALP1_RPT = - 3.*S/(1.+BETA_RHO)*
     $     (1./(S-MMM**2+i*MMM*GRHO_RCHT(S,MMM))+
     $      BETA_RHO/(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG)))
     $      +(2.*GV_RPT/FV_RPT-1.)*
     $      ((2.*Q-2.*S-U)/(1.+BETA_RHO)*
     $      (1./(S-MMM**2+i*MMM*GRHO_RCHT(S,MMM))+
     $      BETA_RHO/(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG)))+
     $      (U-S)/(1.+BETA_RHO)*
     $      (1./(T-MMM**2+i*MMM*GRHO_RCHT(T,MMM))+
     $       BETA_RHO/(T-MMM2**2+i*MMM2*GRHO1_RCHT(T,MMM2,GGG))))

      RETURN

      END



      FUNCTION BETA_RPT(Q,S,T,M1SQ,M2SQ,M3SQ,MMM,MMM2,GGG)
      IMPLICIT NONE
      REAL                      Q,S,T
      DOUBLE PRECISION                M1SQ,M2SQ,M3SQ,MMM,MMM2,GGG
C **********************************************************
C     COMPLEX FUNCTION BETA_RPT ; two-resonance contribution to 
C     F3pi_rcht(iform=1,2) at f3pi_rcht.f;
C     corrected  formula (8) of REF [1] including rho1, 
C     the factor 4*FA_RPT*GV_RPT/(3.*FPI_RPT**3)*QQ/D_a1(QQ) 
C     is included in F3pi_rcht(iform1=1,2) at f3pi_rcht.f  
C **********************************************************
      include '../parameter.inc'
      include '../funct_declar.inc'

      REAL U,LAM0_RPT,LAM1_RPT,LAM2_RPT,FF1_RPT,FF2_RPT,FF_REL    

      FF_REL = FPI_RPT*FPI_RPT/(FV_RPT*FV_RPT)

       LAM1_RPT = FPI_RPT*FPI_RPT/(2.D0*SQRT(2.D0)*FA_RPT*GV_RPT)
       LAM2_RPT = -(1.-2.*FF_REL)*LAM1_RPT
       LAM0_RPT = (LAM1_RPT + LAM2_RPT)/4.

	U= Q-S-T+M1SQ+M2SQ+M3SQ
  
      FF1_RPT = -LAM0_RPT*M1SQ/Q +LAM1_RPT*S/Q+LAM2_RPT
      FF2_RPT = -LAM0_RPT*M1SQ/Q +LAM1_RPT*T/Q+LAM2_RPT

      BETA_RPT = -3.*(LAM1_RPT+LAM2_RPT)*S/(1.+BETA_RHO)*
     $     (1./(S-MMM**2+i*MMM*GRHO_RCHT(S,MMM))+
     $      BETA_RHO/(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG)))+
     $         FF1_RPT*(2.*Q+S-U)/(1.+BETA_RHO)*
     $     (1./(S-MMM**2+i*MMM*GRHO_RCHT(S,MMM))+
     $      BETA_RHO/(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG)))+
     $         FF2_RPT*(U-S)/(1.+BETA_RHO)*
     $     (1./(T-MMM**2+i*MMM*GRHO_RCHT(T,MMM))+
     $      BETA_RHO/(T-MMM2**2+i*MMM2*GRHO1_RCHT(T,MMM2,GGG)))

      RETURN

      END

      FUNCTION ALP1_RPT_RHO1(Q,S,T,M1SQ,M2SQ,M3SQ,
     $                           MMM,MMM2,GGG)
      IMPLICIT NONE
      REAL                      Q,S,T 
      DOUBLE PRECISION                M1SQ,M2SQ,M3SQ,
     &                           MMM,MMM2,GGG
C **********************************************************
C     COMPLEX FUNCTION ALP1_RPT_RHO1;  
C     one-resonance contribution to 
C     F3pi_rcht(iform=1,2) at f3pi_rcht.f;
C     corrected  formula (8) of REF [1] including rho1, 
C     the factor sqrt(2)*FV_RPT*GV_RPT/(3.*FPI_RPT**3) is included in
C     F3pi_rcht(iform1=1,2) at f3pi_rcht.f 
C **********************************************************
      REAL U
      include '../parameter.inc'
      include '../funct_declar.inc'

       U= Q-S-T+M1SQ+M2SQ+M3SQ
       ALP1_RPT_RHO1 = - 3.*S
     $      /(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG))
     $      +(2.*GV1_RPT/FV1_RPT-1.)*
     $      (  (2.*Q-2.*S-U)
     $      /(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG))+
     $      (U-S)/(T-MMM2**2+i*MMM2*GRHO1_RCHT(T,MMM2,GGG))  )

      RETURN

      END



      FUNCTION BETA_RPT_RHO1(Q,S,T,M1SQ,M2SQ,M3SQ,MMM,MMM2,GGG)
      IMPLICIT NONE
      REAL                      Q,S,T
      DOUBLE PRECISION                M1SQ,M2SQ,M3SQ,MMM,MMM2,GGG
C **********************************************************
C     COMPLEX FUNCTION BETA_RPT_RHO1 ; two-resonance contribution to 
C     F3pi_rcht(iform=1,2) at f3pi_rcht.f;
C     corrected  formula (8) of REF [1] including rho1, 
C     the factor 4*FA_RPT*GV_RPT/(3.*FPI_RPT**3)*QQ/D_a1(QQ) 
C     is included in F3pi_rcht(iform1=1,2) at f3pi_rcht.f  
C **********************************************************
      include '../parameter.inc'
      include '../funct_declar.inc'

      REAL U,LAM0_RPT,LAM1_RPT,LAM2_RPT,FF1_RPT,FF2_RPT,FF_REL    

      FF_REL = FPI_RPT*FPI_RPT/(FV1_RPT*FV1_RPT)

       LAM1_RPT = FPI_RPT*FPI_RPT/(2.D0*SQRT(2.D0)*FA_RPT*GV1_RPT)
       LAM2_RPT = -(1.-2.*FF_REL)*LAM1_RPT
       LAM0_RPT = (LAM1_RPT + LAM2_RPT)/4.

	U= Q-S-T+M1SQ+M2SQ+M3SQ
  
      FF1_RPT = -LAM0_RPT*M1SQ/Q +LAM1_RPT*S/Q+LAM2_RPT
      FF2_RPT = -LAM0_RPT*M1SQ/Q +LAM1_RPT*T/Q+LAM2_RPT

      BETA_RPT_RHO1 = -3.*(LAM1_RPT+LAM2_RPT)*S
     $     /(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG)) +
     $         FF1_RPT*(2.*Q+S-U)
     $     /(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG))+
     $         FF2_RPT*(U-S)
     $     /(T-MMM2**2+i*MMM2*GRHO1_RCHT(T,MMM2,GGG))

      RETURN

      END




