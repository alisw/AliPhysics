*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:41  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_NUCSTF_INC
#define CERNLIB_GEANT321_NUCSTF_INC
*
*
* nucstf.inc
*
*
*=== nucstf ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     NUClear STatement Functions:                                     *
*                                                                      *
*     Created on  14 january 1992  by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 04-feb-92     by    Alfredo Ferrari               *
*                                                                      *
*     Included in the following routines:                              *
*                                                                      *
*                        FPFRNC                                        *
*                        FRADNC                                        *
*                        FRHINC                                        *
*                        FRHONC                                        *
*                        SBCOMP                                        *
*                                                                      *
*----------------------------------------------------------------------*
*
*  NUClear STatement Functions: nucleus with a costant core and
*  linearly decreasing skin and halo
*  rho (r):
*     Statement functions:
      RHCORE (R) = RHOCEN
      RHSKIN (R) = RHOCEN * ( RADSKN - R ) / SKNEFF
      RHHALO (R) = RHOSKN * ( RADTOT - R ) / HALODP
*  other
*     RHCORE (R) = RHOCEN
*  Skin with a linear and a cubic term:
*     RHSKIN (R) = RHOCEN * ( ( RADSKN - R ) / SKNEFF + SK3PAR * ( R
*    &           - CPARWS )**3 )
*  Skin with a linear, a cubic and a quartic term:
*     RHSKIN (R) = RHOCEN * ( ( RADSKN - R ) / SKNEFF
*    &           + ( R - CPARWS )**3 * ( SK3PAR + SK4PAR
*    &           * ( R - CPARWS ) ) )
*     RHHALO (R) = RHOCEN * HABPAR * ( RADTOT - R )**1.5D+00
* / r      2
* | 4 pi r' rho (r') dr':
* / 0
*     Statement functions:
      RICORE (R) = 4.D+00 * PI / 3.D+00 * RHOCEN * R**3
      RISKIN (R) = RICORE (RADIU0) + 4.D+00 * PI / 3.D+00 * RHOCEN
     &           * ( RADSKN / SKNEFF * ( R**3 - RADIU0**3 )
     &           - 0.75D+00 * ( R**4 - RADIU0**4 ) / SKNEFF )
      RIHALO (R) = RISKIN (RADIU1) + 4.D+00 * PI / 3.D+00 * RHOSKN
     &           * ( RADTOT / HALODP * ( R**3 - RADIU1**3 )
     &           - 0.75D+00 * ( R**4 - RADIU1**4 ) / HALODP )
*  Other:
*     RICORE (R) = 4.D+00 * PI / 3.D+00 * RHOCEN * R**3
*  Skin with a linear and a cubic term:
*     RISKIN (R) = RICORE (RADIU0) + 4.D+00 * PI / 3.D+00 * RHOCEN
*    &           * ( ( RADSKN / SKNEFF - CPARWS**3 * SK3PAR )
*    &           * ( R**3 - RADIU0**3 ) - 0.75D+00 * ( 1.D+00 / SKNEFF
*    &           - 3.D+00 * SK3PAR * CPARWS**2 ) * ( R**4 - RADIU0**4 )
*    &           - 1.8D+00 * SK3PAR * CPARWS * ( R**5 - RADIU0**5 )
*    &           + 0.5D+00 * SK3PAR * ( R**6 - RADIU0**6 ) )
*  Skin with a linear, a cubic and a quartic term:
*     RISKIN (R) = RICORE (RADIU0) + 4.D+00 * PI / 3.D+00 * RHOCEN
*    &           * ( ( RADSKN / SKNEFF - CPARWS**3 * ( SK3PAR - CPARWS
*    &           * SK4PAR ) ) * ( R**3 - RADIU0**3 ) - 0.75D+00
*    &           * ( 1.D+00 / SKNEFF - CPARWS**2 * ( 3.D+00 * SK3PAR
*    &           - 4.D+00 * SK4PAR * CPARWS ) ) * ( R**4 - RADIU0**4 )
*    &           - 1.8D+00 * CPARWS * ( SK3PAR - 2.D+00 * SK4PAR
*    &           * CPARWS ) * ( R**5 - RADIU0**5 )
*    &           - ( 2.D+00 * SK4PAR * CPARWS - 0.5D+00 * SK3PAR )
*    &           * ( R**6 - RADIU0**6 ) + 3.D+00 / 7.D+00 * SK4PAR
*    &           * ( R**7 - RADIU0**7 ) )
*     RIHALO (R) = RISKIN (RADIU1) + 8.D+00 * PI * RHOCEN
*    &           * HABPAR * ( ( 0.2D+00 * RADTOT**2 - 2.D+00 / 7.D+00
*    &           * RADTOT * HALODP + HALODP**2 / 9.D+00 )
*    &           * ( SQRT (HALODP) )**5 - ( 0.2D+00 * RADTOT**2
*    &           - 2.D+00 / 7.D+00 * RADTOT * ( RADTOT - R )
*    &           + ( RADTOT - R )**2 / 9.D+00 )
*    &           * ( RADTOT - R )**2.5D+00 )
*  r (rho):
*     Statement functions:
      RACORE (RHO) = RADIU0
      RASKIN (RHO) = RADSKN - SKNEFF * RHO / RHOCEN
      RAHALO (RHO) = RADTOT - HALODP * RHO / RHOSKN
*  p_f (rho,i), i=1 proton, 2 neutron:
*     Statement functions:
      PFCORE (RHO,I) = PFRCEN (I)
      PFSKIN (RHO,I) = PFRCEN (I) * (RHO/RHOCEN)**0.3333333333333333D+00
      PFHALO (RHO,I) = PFRCEN (I) * (RHO/RHOCEN)**0.3333333333333333D+00
*  end NUClear STatement Functions
 

#endif
