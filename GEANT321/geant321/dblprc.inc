*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:38  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_DBLPRC_INC
#define CERNLIB_GEANT321_DBLPRC_INC
*
*
* dblprc.inc
*
#if defined(CERNLIB_DOUBLE)
*                                                                     *
*=== dblprc ==========================================================*
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
*      Dblprc: included in any routine                                *
*                                                                     *
*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  *
*  !!!! O N   M A C H I N E S   W H E R E   T H E   D O U B L E !!!!  *
*  !!!! P R E C I S I O N   I S   N O T   R E Q U I R E D  R E -!!!!  *
*  !!!! M O V E   T H E   D O U B L E   P R E C I S I O N       !!!!  *
*  !!!! S T A T E M E N T,  S E T   K A L G N M = 1   A N D     !!!!  *
*  !!!! C H A N G E   A L L   N U M E R I C A L   C O N S -     !!!!  *
*  !!!! T A N T S   T O   S I N G L E   P R E C I S I O N       !!!!  *
*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  *
*                                                                     *
*         Kalgnm = real address alignment, 2 for double precision,    *
*                  1 for single precision                             *
*         Anglgb = this parameter should be set equal to the machine  *
*                  "zero" with respect to unit                        *
*         Anglsq = this parameter should be set equal to the square   *
*                  of Anglgb                                          *
*         Axcssv = this parameter should be set equal to the number   *
*                  for which unity is negligible for the machine      *
*                  accuracy                                           *
*         Andrfl = "underflow" of the machine for floating point      *
*                  operation                                          *
*         Avrflw = "overflow"  of the machine for floating point      *
*                  operation                                          *
*         Ainfnt = code "infinite"                                    *
*         Azrzrz = code "zero"                                        *
*         Einfnt = natural logarithm of the code "infinite"           *
*         Ezrzrz = natural logarithm of the code "zero"               *
*         Onemns = 1- of the machine, it is 1 - 2 x Anglgb            *
*         Onepls = 1+ of the machine, it is 1 + 2 x Anglgb            *
*         Csnnrm = maximum tolerable error on cosine normalization,   *
*                  u**2+v**2+w**2: assuming a typical anglgb relative *
*                  error on each component we would get 2xanglgb: use *
*                  4xanglgb to avoid too many normalizations          *
*         Dmxtrn = "infinite" distance for transport (cm)             *
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( KALGNM = 2 )
      PARAMETER ( ANGLGB = 5.0D-16 )
      PARAMETER ( ANGLSQ = 2.5D-31 )
      PARAMETER ( AXCSSV = 0.2D+16 )
      PARAMETER ( ANDRFL = 1.0D-38 )
      PARAMETER ( AVRFLW = 1.0D+38 )
      PARAMETER ( AINFNT = 1.0D+30 )
      PARAMETER ( AZRZRZ = 1.0D-30 )
      PARAMETER ( EINFNT = +69.07755278982137 D+00 )
      PARAMETER ( EZRZRZ = -69.07755278982137 D+00 )
      PARAMETER ( ONEMNS = 0.999999999999999  D+00 )
      PARAMETER ( ONEPLS = 1.000000000000001  D+00 )
      PARAMETER ( CSNNRM = 2.0D-15 )
      PARAMETER ( DMXTRN = 1.0D+08 )
*
*======================================================================*
*======================================================================*
*=========                                                   ==========*
*=========    M A T H E M A T I C A L   C O N S T A N T S    ==========*
*=========                                                   ==========*
*======================================================================*
*======================================================================*
*                                                                      *
*   Numerical constants:                                               *
*                                                                      *
*         Zerzer = 0                                                   *
*         Oneone = 1                                                   *
*         Twotwo = 2                                                   *
*         Thrthr = 3                                                   *
*         Foufou = 4                                                   *
*         Fivfiv = 5                                                   *
*         Sixsix = 6                                                   *
*         Sevsev = 7                                                   *
*         Eigeig = 8                                                   *
*         Aninen = 9                                                   *
*         Tenten = 10                                                  *
*         Hlfhlf = 1/2                                                 *
*         Onethi = 1/3                                                 *
*         Twothi = 2/3                                                 *
*         Pipipi = Circumference / diameter                            *
*         Eneper = "e", base of natural logarithm                      *
*         Sqrent = square root of "e"                                  *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( ZERZER = 0.D+00 )
      PARAMETER ( ONEONE = 1.D+00 )
      PARAMETER ( TWOTWO = 2.D+00 )
      PARAMETER ( THRTHR = 3.D+00 )
      PARAMETER ( FOUFOU = 4.D+00 )
      PARAMETER ( FIVFIV = 5.D+00 )
      PARAMETER ( SIXSIX = 6.D+00 )
      PARAMETER ( SEVSEV = 7.D+00 )
      PARAMETER ( EIGEIG = 8.D+00 )
      PARAMETER ( ANINEN = 9.D+00 )
      PARAMETER ( TENTEN = 10.D+00 )
      PARAMETER ( HLFHLF = 0.5D+00 )
      PARAMETER ( ONETHI = ONEONE / THRTHR )
      PARAMETER ( TWOTHI = TWOTWO / THRTHR )
      PARAMETER ( PIPIPI = 3.1415926535897932270 D+00 )
      PARAMETER ( ENEPER = 2.7182818284590452354 D+00 )
      PARAMETER ( SQRENT = 1.6487212707001281468 D+00 )
*
*======================================================================*
*======================================================================*
*=========                                                   ==========*
*=========       P H Y S I C A L   C O N S T A N T S         ==========*
*=========                                                   ==========*
*======================================================================*
*======================================================================*
*                                                                      *
*   Primary constants:                                                 *
*                                                                      *
*         Clight = speed of light in cm s-1                            *
*         Avogad = Avogadro number                                     *
*         Amelgr = electron mass (g)                                   *
*         Plckbr = reduced Planck constant (erg s)                     *
*         Elccgs = elementary charge (CGS unit)                        *
*         Elcmks = elementary charge (MKS unit)                        *
*         Amugrm = Atomic mass unit (g)                                *
*         Ammumu = Muon mass (amu)                                     *
*                                                                      *
*   Derived constants:                                                 *
*                                                                      *
*         Alpfsc = Fine structure constant  = e^2/(hbar c)             *
*         Amelct = Electron mass (GeV) = 10^-16Amelgr Clight^2 / Elcmks*
*         Amugev = Atomic mass unit (GeV) = 10^-16Amelgr Clight^2      *
*                                           / Elcmks                   *
*         Ammuon = Muon mass (GeV) = Ammumu * Amugev                   *
*         Fscto2 = (Fine structure constant)^2                         *
*         Fscto3 = (Fine structure constant)^3                         *
*         Fscto4 = (Fine structure constant)^4                         *
*         Plabrc = Reduced Planck constant times the light velocity    *
*                  expressed in GeV fm                                 *
*         Rclsel = Classical electron radius (cm) = e^2 / (m_e c^2)    *
*   Conversion constants:                                              *
*         GeVMeV = from GeV to MeV                                     *
*         eMVGeV = from MeV to GeV                                     *
*         Raddeg = from radians to degrees                             *
*         Degrad = from degrees to radians                             *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( CLIGHT = 2.99792458         D+10 )
      PARAMETER ( AVOGAD = 6.0221367          D+23 )
      PARAMETER ( AMELGR = 9.1093897          D-28 )
      PARAMETER ( PLCKBR = 1.05457266         D-27 )
      PARAMETER ( ELCCGS = 4.8032068          D-10 )
      PARAMETER ( ELCMKS = 1.60217733         D-19 )
      PARAMETER ( AMUGRM = 1.6605402          D-24 )
      PARAMETER ( AMMUMU = 0.113428913        D+00 )
*     PARAMETER ( ALPFSC = 1.D+00 / 137.035989561D+00 )
*     PARAMETER ( FSCTO2 = ALPFSC * ALPFSC )
*     PARAMETER ( FSCTO3 = FSCTO2 * ALPFSC )
*     PARAMETER ( FSCTO4 = FSCTO3 * ALPFSC )
*    It is important to set the electron mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMELCT = 1.D-16 * AMELGR * CLIGHT * CLIGHT / ELCMKS )
*    It is important to set the amu mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMUGEV = 1.D-16 * AMUGRM * CLIGHT * CLIGHT / ELCMKS )
*    It is important to set the muon mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMMUON = AMMUMU * AMUGEV ELCMKS )
*     PARAMETER ( RCLSEL = ELCCGS * ELCCGS / CLIGHT / CLIGHT / AMELGR )
      PARAMETER ( ALPFSC = 7.2973530791728595 D-03 )
      PARAMETER ( FSCTO2 = 5.3251361962113614 D-05 )
      PARAMETER ( FSCTO3 = 3.8859399018437826 D-07 )
      PARAMETER ( FSCTO4 = 2.8357075508200407 D-09 )
      PARAMETER ( PLABRC = 0.197327053        D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMMUON = 0.105658389        D+00 )
      PARAMETER ( RCLSEL = 2.8179409183694872 D-13 )
      PARAMETER ( GEVMEV = 1.0                D+03 )
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( RADDEG = 180.D+00 / PIPIPI )
      PARAMETER ( DEGRAD = PIPIPI / 180.D+00 )
 
#endif
#if defined(CERNLIB_SINGLE)
      IMPLICIT REAL (A-H,O-Z)
      PARAMETER ( KALGNM = 1 )
      PARAMETER ( ANGLGB = 5.0e-16 )
      PARAMETER ( ANGLSQ = 2.5e-31 )
      PARAMETER ( AXCSSV = 0.2e+16 )
      PARAMETER ( ANDRFL = 1.0e-38 )
      PARAMETER ( AVRFLW = 1.0e+38 )
      PARAMETER ( AINFNT = 1.0e+30 )
      PARAMETER ( AZRZRZ = 1.0e-30 )
      PARAMETER ( EINFNT = +69.07755278982137 e+00 )
      PARAMETER ( EZRZRZ = -69.07755278982137 e+00 )
      PARAMETER ( ONEMNS = 0.999999999999999  e+00 )
      PARAMETER ( ONEPLS = 1.000000000000001  e+00 )
      PARAMETER ( CSNNRM = 2.0e-15 )
      PARAMETER ( DMXTRN = 1.0e+08 )
*
      PARAMETER ( ZERZER = 0.e+00 )
      PARAMETER ( ONEONE = 1.e+00 )
      PARAMETER ( TWOTWO = 2.e+00 )
      PARAMETER ( THRTHR = 3.e+00 )
      PARAMETER ( FOUFOU = 4.e+00 )
      PARAMETER ( FIVFIV = 5.e+00 )
      PARAMETER ( SIXSIX = 6.e+00 )
      PARAMETER ( SEVSEV = 7.e+00 )
      PARAMETER ( EIGEIG = 8.e+00 )
      PARAMETER ( ANINEN = 9.e+00 )
      PARAMETER ( TENTEN = 10.e+00 )
      PARAMETER ( HLFHLF = 0.5e+00 )
      PARAMETER ( ONETHI = ONEONE / THRTHR )
      PARAMETER ( TWOTHI = TWOTWO / THRTHR )
      PARAMETER ( PIPIPI = 3.1415926535897932270 e+00 )
      PARAMETER ( ENEPER = 2.7182818284590452354 e+00 )
      PARAMETER ( SQRENT = 1.6487212707001281468 e+00 )
*
      PARAMETER ( CLIGHT = 2.99792458         e+10 )
      PARAMETER ( AVOGAD = 6.0221367          e+23 )
      PARAMETER ( AMELGR = 9.1093897          e-28 )
      PARAMETER ( PLCKBR = 1.05457266         e-27 )
      PARAMETER ( ELCCGS = 4.8032068          e-10 )
      PARAMETER ( ELCMKS = 1.60217733         e-19 )
      PARAMETER ( AMUGRM = 1.6605402          e-24 )
      PARAMETER ( AMMUMU = 0.113428913        e+00 )
*     PARAMETER ( ALPFSC = 1.e+00 / 137.035989561e+00 )
*     PARAMETER ( FSCTO2 = ALPFSC * ALPFSC )
*     PARAMETER ( FSCTO3 = FSCTO2 * ALPFSC )
*     PARAMETER ( FSCTO4 = FSCTO3 * ALPFSC )
*    It is important to set the electron mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMELCT = 1.e-16 * AMELGR * CLIGHT * CLIGHT / ELCMKS )
*    It is important to set the amu mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMUGEV = 1.e-16 * AMUGRM * CLIGHT * CLIGHT / ELCMKS )
*    It is important to set the muon mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMMUON = AMMUMU * AMUGEV ELCMKS )
*     PARAMETER ( RCLSEL = ELCCGS * ELCCGS / CLIGHT / CLIGHT / AMELGR )
      PARAMETER ( ALPFSC = 7.2973530791728595 e-03 )
      PARAMETER ( FSCTO2 = 5.3251361962113614 e-05 )
      PARAMETER ( FSCTO3 = 3.8859399018437826 e-07 )
      PARAMETER ( FSCTO4 = 2.8357075508200407 e-09 )
      PARAMETER ( PLABRC = 0.197327053        e+00 )
      PARAMETER ( AMELCT = 0.51099906         e-03 )
      PARAMETER ( AMUGEV = 0.93149432         e+00 )
      PARAMETER ( AMMUON = 0.105658389        e+00 )
      PARAMETER ( RCLSEL = 2.8179409183694872 e-13 )
      PARAMETER ( GEVMEV = 1.0                e+03 )
      PARAMETER ( EMVGEV = 1.0                e-03 )
      PARAMETER ( RADDEG = 180.e+00 / PIPIPI )
      PARAMETER ( DEGRAD = PIPIPI / 180.e+00 )
 
#endif

#endif
