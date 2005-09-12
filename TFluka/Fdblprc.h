#ifndef FDBLPRC_H
#define FDBLPRC_H 1

#include "Rtypes.h"
#include "cfortran.h"
extern "C" {
//*$ create dblprc.add
//*copy dblprc
//*                                                                     *
//*=== dblprc ==========================================================*
//*                                                                     *
//*---------------------------------------------------------------------*
//*                                                                     *
//*      dblprc: included in any routine, machine, mathematical and     *
//*              physical constants plus global declarations            *
//*                                                                     *
//*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  *
//*  !!!! o n   m a c h i n e s   w h e r e   t h e   d o u b l e !!!!  *
//*  !!!! p r e c i s i o n   i s   n o t   r e q u i r e d  r e -!!!!  *
//*  !!!! m o v e   t h e   d o u b l e   p r e c i s i o n       !!!!  *
//*  !!!! s t a t e m e n t,  s e t   k a l g n m = 1   a n d     !!!!  *
//*  !!!! c h a n g e   a l l   n u m e r i c a l   c o n s -     !!!!  *
//*  !!!! t a n t s   t o   s i n g l e   p r e c i s i o n       !!!!  *
//*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  *
//*                                                                     *
//*         kalgnm = real address alignment, 2 for double precision,    *
//*                  1 for single precision                             *
//*         kalch8 = character*8 address alignment wrt the precision    *
//*                  defined by kalgnm (mostly 1 in all situations)     *
//*         i2algn = integer*2 address alignment wrt the normal integer *
//*                  precision (mostly 2, 4 for 64 bit integers)        *
//*         anglgb = this parameter should be set equal to the machine  *
//*                  "zero" with respect to unit                        *
//*         anglsq = this parameter should be set equal to the square   *
//*                  of anglgb                                          *
//*         axcssv = this parameter should be set equal to the number   *
//*                  for which unity is negligible for the machine      *
//*                  accuracy                                           *
//*         andrfl = "underflow" of the machine for floating point      *
//*                  operation                                          *
//*         avrflw = "overflow"  of the machine for floating point      *
//*                  operation                                          *
//*         ainfnt = code "infinite"                                    *
//*         azrzrz = code "zero"                                        *
//*         einfnt = natural logarithm of the code "infinite"           *
//*         ezrzrz = natural logarithm of the code "zero"               *
//*         excssv = natural logarithm of the code number for which     *
//*                  unit is negligible                                 *
//*         englgb = natural logarithm of the code "zero" with respect  *
//*                  to unit                                            *
//*         onemns = 1- of the machine, it is 1 - 2 x anglgb            *
//*         onepls = 1+ of the machine, it is 1 + 2 x anglgb            *
//*         csnnrm = maximum tolerable error on cosine normalization,   *
//*                  u**2+v**2+w**2: assuming a typical anglgb relative *
//*                  error on each component we would get 2xanglgb: use *
//*                  4xanglgb to avoid too many normalizations          *
//*         dmxtrn = "infinite" distance for transport (cm)             *
//*         rhflmn = minimal density for fluka (g/cm^3)                 *
//*                                                                     *
//*   "global" declarations:                                            *
//*         lfluka = set to true for a real (full) fluka run            *
//*         lgbias = set to true for a fully biased run                 *
//*         lgbana = set to true for a fully analogue run               *
//*         lflgeo = set to true when using the standard fluka geometry *
//*         loflts = set to true for special off-line testing of speci- *
//*                  fic routines                                       *
//*         lusrin = set to true if the user dependent initialization   *
//*                  routine usrini has been called at least once       *
//*         lnmgeo = set to true for a name-base geometry input         *
//*         lnminp = set to true for a name-base fluka input            *
//*         Lfrfmt = set to true for a free-format based Fluka input    *
//*         lfdrtr = set to true for going in/out feeder/flukam at each *
//*                  event                                              *
//*                                                                     *
//*---------------------------------------------------------------------*
//*                                                                     *
const Int_t kalgnm = 2;
const Int_t kalch8 = 1;
const Int_t i2algn = 2;
const Double_t anglgb = 5.0e-16;
const Double_t anglsq = 2.5e-31;
const Double_t axcssv = 0.2e+16;
const Double_t andrfl = 1.0e-38;
const Double_t avrflw = 1.0e+38;
const Double_t ainfnt = 1.0e+30;
const Double_t azrzrz = 1.0e-30;
const Double_t einfnt = +69.07755278982137e+00;
const Double_t ezrzrz = -69.07755278982137e+00;
const Double_t excssv = +35.23192357547063e+00;
const Double_t englgb = -35.23192357547063e+00;
const Double_t onemns = 0.999999999999999e+00;
const Double_t onepls = 1.000000000000001e+00;
const Double_t csnnrm = 2.0e-15;
const Double_t dmxtrn = 1.0e+08;
const Double_t rhflmn = 1.0e-10;
//*
//*======================================================================*
//*======================================================================*
//*=========                                                   ==========*
//*=========    m a t h e m a t i c a l   c o n s t a n t s    ==========*
//*=========                                                   ==========*
//*======================================================================*
//*======================================================================*
//*                                                                      *
//*   numerical constants (single precision):                            *
//*                                                                      *
//*         zersng = 0                                                   *
//*                                                                      *
//*   numerical constants (double precision):                            *
//*                                                                      *
//*         zerzer = 0                                                   *
//*         oneone = 1                                                   *
//*         twotwo = 2                                                   *
//*         thrthr = 3                                                   *
//*         foufou = 4                                                   *
//*         fivfiv = 5                                                   *
//*         sixsix = 6                                                   *
//*         sevsev = 7                                                   *
//*         eigeig = 8                                                   *
//*         aninen = 9                                                   *
//*         tenten = 10                                                  *
//*         eleven = 11                                                  *
//*         twelve = 12                                                  *
//*         fiften = 15                                                  *
//*         sixten = 16                                                  *
//*         hlfhlf = 1/2                                                 *
//*         onethi = 1/3                                                 *
//*         onefou = 1/4                                                 *
//*         onefiv = 1/5                                                 *
//*         onesix = 1/6                                                 *
//*         onesev = 1/7                                                 *
//*         oneeig = 1/8                                                 *
//*         twothi = 2/3                                                 *
//*         thrfou = 3/4                                                 *
//*         thrtwo = 3/2                                                 *
//*         pipipi = circumference / diameter                            *
//*         twopip = 2 x pipipi                                          *
//*         pip5o2 = 5/2 x pipipi                                        *
//*         pipisq = pipipi x pipipi                                     *
//*         pihalf = 1/2 x pipipi                                        *
//*         erfa00 = erf (oo) = 1/2 x square root of pi                  *
//*         sqtwpi = square root of 2xpi                                 *
//*         eulero = eulero's constant                                   *
//*         eulexp = exp ( eulero )                                      *
//*         e1m2eu = exp ( 1 - 2 eulero )                                *
//*         eneper = "e", base of natural logarithm                      *
//*         sqrent = square root of "e"                                  *
//*         sqrtwo = square root of  2                                   *
//*         sqrthr = square root of  3                                   *
//*         sqrfiv = square root of  5                                   *
//*         sqrsix = square root of  6                                   *
//*         sqrsev = square root of  7                                   *
//*         sqrt12 = square root of 12                                   *
//*         s2fwhm = 2 x square root of 2 x logarithm of 2               *
//*                                                                      *
//*----------------------------------------------------------------------*
//*
const Float_t  zersng = 0.e+00;
const Double_t zerzer = 0.e+00;
const Double_t oneone = 1.e+00;
const Double_t twotwo = 2.e+00;
const Double_t thrthr = 3.e+00;
const Double_t foufou = 4.e+00;
const Double_t fivfiv = 5.e+00;
const Double_t sixsix = 6.e+00;
const Double_t sevsev = 7.e+00;
const Double_t eigeig = 8.e+00;
const Double_t aninen = 9.e+00;
const Double_t tenten = 10.e+00;
const Double_t eleven = 11.e+00;
const Double_t twelve = 12.e+00;
const Double_t fiften = 15.e+00;
const Double_t sixten = 16.e+00;
const Double_t hlfhlf = 0.5e+00;
const Double_t onethi = oneone/thrthr;
const Double_t onefou = oneone/foufou;
const Double_t onefiv = oneone/fivfiv;
const Double_t onesix = oneone/sixsix;
const Double_t onesev = oneone/sevsev;
const Double_t oneeig = oneone/eigeig;
const Double_t twothi = twotwo/thrthr;
const Double_t thrfou = thrthr/foufou;
const Double_t thrtwo = thrthr/twotwo;
const Double_t fouthr = foufou/thrthr;    
const Double_t pipipi = 3.141592653589793238462643383279e+00;
const Double_t twopip = 6.283185307179586476925286766559e+00;
const Double_t pip5o2 = 7.853981633974483096156608458199e+00;
const Double_t pipisq = 9.869604401089358618834490999876e+00;
const Double_t pihalf = 1.570796326794896619231321691640e+00;
const Double_t erfa00 = 0.886226925452758013649083741671e+00;
const Double_t sqrtpi = 1.772453850905516027298167483341e+00;
const Double_t sqtwpi = 2.506628274631000502415765284811e+00;
const Double_t eulero = 0.577215664901532860606512e+00;
const Double_t eulexp = 1.781072417990197985236504e+00;
const Double_t eullog = -0.5495393129816448223376619e+00;
const Double_t e1m2eu = 0.8569023337737540831433017e+00;
const Double_t eneper = 2.718281828459045235360287471353e+00;
const Double_t sqrent = 1.648721270700128146848650787814e+00;
const Double_t sqrtwo = 1.414213562373095048801688724210e+00;
const Double_t sqrthr = 1.732050807568877293527446341506e+00;
const Double_t sqrfiv = 2.236067977499789696409173668731e+00;
const Double_t sqrsix = 2.449489742783178098197284074706e+00;
const Double_t sqrsev = 2.645751311064590590501615753639e+00;
const Double_t sqrt12 = 3.464101615137754587054892683012e+00;
const Double_t s2fwhm = 2.354820045030949e+00;
const Double_t twolog = 0.693147180559945309417232121458e+00;
//*
//*======================================================================*
//*======================================================================*
//*=========                                                   ==========*
//*=========       p h y s i c a l   c o n s t a n t s         ==========*
//*=========                                                   ==========*
//*======================================================================*
//*======================================================================*
//*                                                                      *
//*   primary constants:                                                 *
//*                                                                      *
//*         clight = speed of light in cm s-1                            *
//*         avogad = avogadro number                                     *
//*         boltzm = k boltzmann constant (j k-1)                        *
//*         amelgr = electron mass (g)                                   *
//*         plckbr = reduced planck constant (erg s)                     *
//*         elccgs = elementary charge (cgs unit)                        *
//*         elcmks = elementary charge (mks unit)                        *
//*         amugrm = atomic mass unit (g)                                *
//*         ammumu = muon    mass (amu)                                  *
//*         amprmu = proton  mass (amu)                                  *
//*         amnemu = neutron mass (amu)                                  *
//*                                                                      *
//*   derived constants:                                                 *
//*                                                                      *
//*         alpfsc = fine structure constant  = e^2/(hbar c) (cgs units) *
//*         amelct = electron mass (gev) = 10^-16amelgr clight^2 / elcmks*
//*         amugev = atomic mass unit (gev) = 10^-16amugrm clight^2      *
//*                                           / elcmks                   *
//*         ammuon = muon    mass (gev) = ammumu * amugev                *
//*         amprtn = proton  mass (gev) = amprmu * amugev                *
//*         amntrn = neutron mass (gev) = amnemu * amugev                *
//*         amdeut = deuteron mass (gev)                                 *
//*         amalph = alpha    mass (gev) (derived from the excess mass   *
//*                  and an (approximate) atomic binding not a really    *
//*                  measured constant)                                  *
//*         cougfm = e^2 (gev fm) = elccgs^2 / elcmks * 10^-7 * 10^-9    *
//*                * 10^13 (10^..=erg cm->joule cm->gev cm->gev fm       *
//*                it is equal to 0.00144 gev fm                         *
//*         fscto2 = (fine structure constant)^2                         *
//*         fscto3 = (fine structure constant)^3                         *
//*         fscto4 = (fine structure constant)^4                         *
//*         plabrc = reduced planck constant times the light velocity    *
//*                  expressed in gev fm                                 *
//*         rclsel = classical electron radius (cm) = e^2 / (m_e c^2)    *
//*         bltzmn = k boltzmann constant in gev k-1                     *
//*         a0bohr = bohr radius, hbar^2 / ( m_e e^2) (fm) = plabrc**2   *
//*                / amelct / cougfm, or equivalently,                   *
//*                plabrc / alpfsc / amelct                              *
//*         gfohb3 = fermi constant, g_f/(hbar c)^3, in gev^-2           *
//*         gfermi = fermi constant in gev fm^3                          *
//*         sin2tw = sin^2 theta_weinberg                                *
//*         prmgnm = proton  magnetic moment (magneton)                  *
//*         anmgnm = neutron magnetic moment (magneton)                  *
//*         s0thms = sigma_0 Thomson, 8/3 pi r_e^2 (mb)                  *
//*                                                                      *
//*   astronomical constants:                                            *
//*                                                                      *
//*         rearth = earth equatorial radius (cm)                        *
//*         auastu = astronomical unit       (cm)                        *
//*                                                                      *
//*   conversion constants:                                              *
//*                                                                      *
//*         gevmev = from gev to mev                                     *
//*         emvgev = from mev to gev                                     *
//*         gev2ev = from gev to  ev                                     *
//*         ev2gev = from ev  to gev                                     *
//*         algvmv = from gev to mev, log                                *
//*         raddeg = from radians to degrees                             *
//*         degrad = from degrees to radians                             *
//*         gevomg = from (photon) energy [gev] in 2pi x frequency [s^-1]*
//*         cmq2mb = from square centimetres to millibarns               *
//*                                                                      *
//*   useful constants:                                                  *
//*                                                                      *
//*         fertho = constant to be used in the fermi-thomas approxima-  *
//*                  ted expression for atomic binding energies          *
//*         expebn = exponent to be used in the fermi-thomas approxima-  *
//*                  ted expression for atomic binding energies          *
//*                    b_atomic (z) = fertho x z^expebn (gev)            *
//*         bexc12 = fermi-thomas approximated expression for 12-c ato-  *
//*                  mic binding energies (gev)                          *
//*         amunmu = difference between the atomic and nuclear mass units*
//*         amuc12 = "nuclear" mass unit = 1/12 m_nucl (12-c),           *
//*                  m_nucl (12-c) = m_atom (12-c) - 6 m_e + b_atom(12-c)*
//*                                                                      *
//*----------------------------------------------------------------------*
//*
const Double_t clight = 2.99792458e+10;
const Double_t avogad = 6.0221367e+23;
const Double_t boltzm = 1.380658e-23;
const Double_t amelgr = 9.1093897e-28;
const Double_t plckbr = 1.05457266e-27;
const Double_t elccgs = 4.8032068e-10;
const Double_t elcmks = 1.60217733e-19;
const Double_t amugrm = 1.6605402e-24;
const Double_t ammumu = 0.113428913e+00;
const Double_t amprmu = 1.007276470e+00;
const Double_t amnemu = 1.008664904e+00;
//* const Double_t alpfsc = 1.e+00 / 137.035989561e+00
//* const Double_t fscto2 = alpfsc * alpfsc
//* const Double_t fscto3 = fscto2 * alpfsc
//* const Double_t fscto4 = fscto3 * alpfsc
//*    it is important to set the electron mass exactly with the same
//*    rounding as in the mass tables, so use the explicit expression
//* const Double_t amelct = 1.e-16 * amelgr * clight * clight / elcmks
//*    it is important to set the amu mass exactly with the same
//*    rounding as in the mass tables, so use the explicit expression
//* const Double_t amugev = 1.e-16 * amugrm * clight * clight / elcmks
//*    it is important to set the muon,proton,neutron masses exactly with
//*    the same rounding as in the mass tables, so use the explicit
//*    expression
//* const Double_t ammuon = ammumu * amugev
//* const Double_t amprtn = amprmu * amugev
//* const Double_t amntrn = amnemu * amugev
//* const Double_t rclsel = elccgs * elccgs / clight / clight / amelgr
//* const Double_t bltzmn = boltzm / elcmks * 1.e-09
const Double_t alpfsc = 7.2973530791728595e-3;
const Double_t fscto2 = 5.3251361962113614e-5;
const Double_t fscto3 = 3.8859399018437826e-7;
const Double_t fscto4 = 2.8357075508200407e-9;
const Double_t plabrc = 0.197327053e+00;
const Double_t amelct = 0.51099906e-3;
const Double_t amugev = 0.93149432e+00;
const Double_t ammuon = 0.105658389e+00;
const Double_t amprtn = 0.93827231e+00;
const Double_t amntrn = 0.93956563e+00;
const Double_t amdeut = 1.87561339e+00;
const Double_t amalph = 3.72738025692891e+00;
const Double_t cougfm = elccgs*elccgs/elcmks*(1.e-7)*(1.e+13)*(1.e-9);
const Double_t rclsel = 2.8179409183694872e-13;
const Double_t alamb0 = twotwo * pipipi * rclsel / alpfsc;
const Double_t bltzmn = 8.617385e-14;
const Double_t a0bohr = plabrc/alpfsc/amelct;
const Double_t gfohb3 = 1.16639e-5;
const Double_t gfermi = gfohb3*plabrc*plabrc*plabrc;
const Double_t sin2tw = 0.2319e+00;
const Double_t prmgnm = 2.792847386e+00;
const Double_t anmgnm = -1.91304275e+00;
const Double_t rearth = 6.378140e+8;
const Double_t auastu = 1.4959787066e+13;
const Double_t gevmev = 1.0e+3;
const Double_t ev2gev = 1.0e-9;
const Double_t gev2ev = 1.0e+9;
const Double_t emvgev = 1.0e-3;
const Double_t cmq2mb = 1.0e+27;
const Double_t fmb2ba = 1.0e-3;
const Double_t bar2mb = 1.0e+3;
const Double_t fmb2fs = 1.0e-1;
const Double_t fms2mb = 1.0e+1;
const Double_t algvmv = 6.90775527898214e+00;
const Double_t raddeg = (180.e+00)/pipipi;
const Double_t degrad = pipipi/(180.e+00);
const Double_t gevomg = clight*(1.e+13)/plabrc;
const Double_t s0thms = eigeig / thrthr * pipipi * rclsel * rclsel * cmq2mb;
//*  old Fermi-Thomas parametrization of atomic binding energies:
//*     const Double_t fertho = 15.73       e-9
//*     const Double_t expebn = 7.e+00 / 3.e+00
//*     const Double_t bexc12 = fertho * 65.41634134195703e+00
//*  new Fermi-Thomas parametrization of atomic binding energies:
const Double_t fertho = 14.33e-9;
const Double_t expebn = 2.39e+00;
const Double_t bexc12 = fertho*72.40715579499394e+00;
const Double_t amunmu = hlfhlf*amelct-bexc12/12.e+00;
const Double_t amuc12 = amugev-amunmu;
//*  Old MeV units:
const Double_t amemev = gevmev * amelct;
//*

typedef struct {
   Int_t    lfluka;
   Int_t    lgbias;
   Int_t    lgbana;
   Int_t    lflgeo;
   Int_t    loflts;
   Int_t    lusrin;
   Int_t    lnmgeo;
   Int_t    lnminp;
   Int_t    lfrfmt;
   Int_t    lfdrtr;
   Int_t    kflgeo;
   Int_t    kfldnr;
} globalCommon;
#define GLOBAL COMMON_BLOCK(GLOBAL,global)
COMMON_BLOCK_DEF(globalCommon,GLOBAL);
}

#endif
