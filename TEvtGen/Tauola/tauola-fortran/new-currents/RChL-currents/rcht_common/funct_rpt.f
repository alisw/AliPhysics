
C **********************************************************
C library of functions used in calculation of currents 
C references:
C [1]   arXiv:0911.4436 (hep-ph)  D. Gomez Dumm et al. (tau -> 3pi nu) 
C [2]   arXiv:0911.2640 (hep-ph) D. Gomezz Dumm et al. (tau -> KKpi nu)
C [3]   arXiv:0807.4883 (hep-ph) D. R. Boito et al. (tau -> Kpi nu)
C [4]   P. Roig, talk at (tau -> 2 pi nu)
C [5]   arXiv:0803.2039  (hep-ph) E. Arganda et al., Appendix B (tau -> KK nu)
C **********************************************************
      FUNCTION GRHO_RCHT(XS,XMMM)
      IMPLICIT NONE
      REAL               XS
      DOUBLE PRECISION      XMMM
C **********************************************************
C     REAL FUNCTION Gamma Rho ; energy-dependent width of  rho meson
c     in SU(2) limit mpi=mpi0, mk=mk0;
C     formula (14) of REF [1]
C **********************************************************
       include '../parameter.inc'
      include '../funct_declar.inc'

       REAL MMPI_AV2,MMK_2

      MMPI_AV2 = MMPI_AV**2
      MMK_2 = MMK**2

      IF(XS.GE.(4.*MMK_2)) THEN
           GRHO_RCHT=XMMM*XS*((1.-4.*MMPI_AV2/XS)**1.5
     $              +0.5*(1.-4.*MMK_2/XS)**1.5)
     $                  /(96.*PI*FPI_RPT**2)
      ELSE IF((XS.GE.(4.*MMPI_AV2)).AND.(XS.LE.(4.*MMK_2))) THEN 
        GRHO_RCHT=XMMM*XS*(1.-4.*MMPI_AV2/XS)**1.5
     $            /(96.*PI*FPI_RPT**2)
      ELSE
        GRHO_RCHT = 0.
      ENDIF

      RETURN

      END


      FUNCTION GRHO1_RCHT(XS,XMMM,XGGG)
      IMPLICIT NONE
      REAL                XS
      DOUBLE PRECISION       XMMM,XGGG
C **********************************************************
C     REAL FUNCTION Gamma Rho1 ; energy-dependent width of  rho1 meson;
c     in SU(2) limit mpi=mpi0; only rho' -> 2pi loop;
C     formula (33) of REF [1]
C **********************************************************
      include '../parameter.inc'
      include '../funct_declar.inc'

       REAL MMPI_AV2
      MMPI_AV2 = MMPI_AV**2

      IF (XS.GE.(4.*MMPI_AV2)) THEN 
        GRHO1_RCHT=XGGG*SQRT(XMMM**2/XS)*  
     &      ((XS-4.*MMPI_AV2)/(XMMM**2-4.*MMPI_AV2))**1.5
      ELSE
        GRHO1_RCHT = 0.
      ENDIF

      RETURN
      END



      FUNCTION SIGP(SS)
      IMPLICIT NONE
      DOUBLE PRECISION               SS
C***********************************************************
C     DOUBLRE PRECISION FUNCTION 
c     of two body phase space threshold: equal mass scalars
C***********************************************************
      REAL TT
      include '../parameter.inc'
      include '../funct_declar.inc'
      
      TT = 1. - 4.*MMPI_AV**2/SS 

      IF (TT.GE.0) THEN 
      SIGP = SQRT(TT)
      ELSE
      SIGP = 0.
      ENDIF

      RETURN
      END




      FUNCTION LAMB_RCHT(X1,X2,X3)
      IMPLICIT NONE
      REAL               X1,X2,X3,ARG_RCHT      
C***********************************************************
C     REAL FUNCTION LAMDA of three body phase space
C***********************************************************
      include '../funct_declar.inc'

      ARG_RCHT = (X1-X2-X3)**2 - 4.*X3*X2
      IF(ARG_RCHT.GE.0.) THEN
      LAMB_RCHT = ARG_RCHT
      ELSE
      LAMB_RCHT = 0.
      ENDIF

      RETURN
      END



C******************** FUNCTIONS FOR 2 SCALAR MODES ***********************



C****************************************************************
      FUNCTION R0SCAL_3PI(QX,SX)
C****************************************************************
C     Complex Function: R_0 function (Pablo notes) 
C     for the scalar contribution for three pion modes
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      REAL            QX,SX,delta0_3piscal,xsx,xst
      DOUBLE PRECISION    XM1,XM2,DSX
      include '../funct_declar.inc'
       include '../parameter.inc'
c$$$      a00_3piscal = 0.220
c$$$      b00_3piscal = 0.268/mmpi_av**2
c$$$      c00_3piscal = -0.0139/mmpi_av**4
c$$$      d00_3piscal = -0.00139/mmpi_av**6
c$$$      x00_3piscal = 36.77*mmpi_av**2
c$$$c      MMF0 = 0.98
c$$$
      dsx = sx
      xsx = sx/4.*sigp(dsx)**2
      xst = sqrt(sx)

      if(sx.le.0.7) then
      delta0_3piscal = sigp(dsx)*(a00_3piscal + b00_3piscal*xsx +
     &               c00_3piscal*xsx**2 + d00_3piscal*xsx**3)

      delta0_3piscal = delta0_3piscal*
     &         (4*mmpi_av**2 - x00_3piscal)/(sx-x00_3piscal) 
      else if (xst.le.1.21) then
      delta0_3piscal = -10572.0+50658.0*xst-87903.0*xst**2+66886.0*xst**3
     & -18699.0*xst**4
      delta0_3piscal = delta0_3piscal*pi/180.0
      else
      delta0_3piscal = 255.0*pi/180.0
      endif
      

      delta0_3piscal = atan(delta0_3piscal)

      R0SCAL_3PI = ALPHA0_3PI/QX + 
     &          ALPHA1_3PI/QX**2*(SX - MMF0**2)

      R0SCAL_3PI = R0SCAL_3PI*(cos(delta0_3piscal) +
     &              i*sin(delta0_3piscal))

      RETURN
      END

C****************************************************************
      FUNCTION R2SCAL_3PI(QX,SX)
C****************************************************************
C     Complex Function: R_2 function (Pablo notes) 
C     for the scalar contribution for three pion modes
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      REAL            QX,SX,delta2_3piscal,xsx,xst
      DOUBLE PRECISION    XM1,XM2,DSX
      include '../funct_declar.inc'
       include '../parameter.inc'
c$$$      a02_3piscal = -0.0444
c$$$      b02_3piscal = -0.0857/mmpi_av**2
c$$$      c02_3piscal = -0.00221/mmpi_av**4
c$$$      d02_3piscal = -0.000129/mmpi_av**6
c$$$      x02_3piscal = -21.62*mmpi_av**2
c$$$c      MMF0 = 0.98

      dsx = sx
      xsx = sx/4.*sigp(dsx)**2
      xst = sqrt(sx)

      if(sx.le.0.7) then
      delta2_3piscal = sigp(dsx)*(a02_3piscal + b02_3piscal*xsx +
     &               c02_3piscal*xsx**2 + d02_3piscal*xsx**3)

      delta2_3piscal = delta2_3piscal*
     &         (4*mmpi_av**2 - x02_3piscal)/(sx-x02_3piscal) 
      else if(xst.le.1.21) then
      delta2_3piscal = 282.9-1314.9*xst+2153.4*xst**2-1574.5d0*xst**3+
     & 428.06d0*xst**4
      delta2_3piscal = delta2_3piscal*pi/180.0
      else
      delta2_3piscal = -27.0*pi/180.0
      endif

      delta2_3piscal = atan(delta2_3piscal)      

      R2SCAL_3PI = GAMMA0_3PI/QX + 
     &          GAMMA1_3PI/QX**2*(SX - MMF0**2)

      R2SCAL_3PI = R2SCAL_3PI*(cos(delta2_3piscal) +
     &              i*sin(delta2_3piscal))

      RETURN
      END




C*************************************************************************
C         Functions for sigma contributions
C*************************************************************************
      FUNCTION FFsig(QX,XX)
C **********************************************************
C     Complex Function:  
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      REAL            QX,XX,mm2
      DOUBLE PRECISION    XM1,XM2,xphi
      include '../funct_declar.inc'
      include '../parameter.inc'

      mm2 = MMPI_AV**2

      xphi = - rsigma**2* LAMB_RCHT(QX,XX,mm2)/(8.*QX)
   
      FFsig = dexp(xphi)


       RETURN
       END



C*************************************************************************
      FUNCTION BWsig(XM,XG,XQ)
C **********************************************************
C     Complex Function: S-wave Breit-Wigner  
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      REAL            XQ
      DOUBLE PRECISION    XM,XG,XM2,XXQ,GAMMA
      include '../funct_declar.inc'
      include '../parameter.inc'

      XXQ = XQ
      XM2 = XM**2
      GAMMA = XG*SIGP(XXQ)/SIGP(XM2)
      BWsig = XM*XM/CMPLX(XM*XM-XQ, -XM*GAMMA)


       RETURN
       END

C*************************************************************************
      FUNCTION DECOUL(mm1,mm3,ss2)
C **********************************************************
C     Real Function: Coulomb interaction effects of two particles
C                    with mm1 and mm2 
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      REAL            ss2
      DOUBLE PRECISION    mm1,mm3,betam1m3

      include '../funct_declar.inc'
      include '../parameter.inc'
C*******************************************
C   COMMON  block fixed in SUBROUTINE INIPHY
c*******************************************
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
      betam1m3 = 1.d0 - (mm1 +mm3)**2/ss2
      betam1m3 = dsqrt(betam1m3)

      if(ss2.gt.(mm1+mm3)**2)      
     &    decoul = pi/2.d0/betam1m3/ALFINV

       RETURN
       END

C*************************************************************************
      FUNCTION COUL3PART(mm1,mm2,mm3,ss1,ss2,ss3)
C **********************************************************
C     Real Function: Coulomb interaction effects of two particles
C                    with mm1 and mm2 
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      REAL            ss3,ss2,ss1
      DOUBLE PRECISION    mm1,mm2,mm3
      include '../funct_declar.inc'
      include '../parameter.inc'

      coul3part = decoul(mm1,mm3,ss2) + decoul(mm2,mm3,ss1) 
     &          - decoul(mm1,mm2,ss3)

       coul3part = exp(coul3part)

       RETURN
       END
C*************************************************************************
      FUNCTION fattcoul(mm1,mm3,ss2)
C **********************************************************
C     Real Function: Coulomb attraction of two particles
C                    with mm1 and mm3 
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      REAL            ss2
      DOUBLE PRECISION    mm1,mm3,betam1m3

      include '../funct_declar.inc'
      include '../parameter.inc'
C*******************************************
C   COMMON  block fixed in SUBROUTINE INIPHY
c*******************************************
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
      if(ss2.gt.(mm1+mm3)**2) then      
      betam1m3 = 2.*dsqrt(1.d0 - (mm1 +mm3)**2/ss2)
     &            /(1.+ (1.d0 - (mm1 +mm3)**2/ss2))


         fattcoul = 2.*pi/betam1m3/ALFINV
     &               /(1.-exp(-2.*pi/betam1m3/ALFINV))
       else
         fattcoul = 1
          endif

       RETURN
       END

C*************************************************************************
C*************************************************************************
      FUNCTION frepcoul(mm1,mm3,ss2)
C **********************************************************
C     Real Function: Coulomb repuslcion of two particles
C                    with mm1 and mm3 
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      REAL            ss2
      DOUBLE PRECISION    mm1,mm3,betam1m3

      include '../funct_declar.inc'
      include '../parameter.inc'
C*******************************************
C   COMMON  block fixed in SUBROUTINE INIPHY
c*******************************************
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
        if(ss2.gt.(mm1+mm3)**2) then      
      betam1m3 = 2.*dsqrt(1.d0 - (mm1 +mm3)**2/ss2)
     &           /(1.+ (1.d0 - (mm1 +mm3)**2/ss2))

         frepcoul = 2.*pi/betam1m3/ALFINV
     &               /(-1.+exp(2.*pi/betam1m3/ALFINV))
       else
         frepcoul = 1
          endif
       RETURN
       END

C*************************************************************************
