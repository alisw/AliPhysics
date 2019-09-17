// SigmaTotal.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains classes for cross section parametrizations.
// SigmaTotAux : base class for the different parametrizations.
// SigmaTotal  : top-level class, making use of the classes below.
// SigmaTotOwn : user-set values.
// SigmaSaSDL  : Schuler and Sjostrand, based on Donnachie and Landshoff.
// SigmaMBR    : Min Bias Rockefeller.
// SigmaABMST  : Appleby, Barlow, Molson, Serluca and Toader.
// SigmaRPP    : Review of Particle Physics 2014.

#ifndef Pythia8_SigmaTotal_H
#define Pythia8_SigmaTotal_H

#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/PythiaComplex.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// The SigmaTotAux is base class for the different parametrizations
// made accessible via SigmaTotal. Many variables are public, since
// they can only be accessed via the SigmaTotal methods anyway.

class SigmaTotAux {

public:

  // Constructor.
  SigmaTotAux() : isExpEl(), hasCou(), sigTot(), rhoOwn(), sigEl(), bEl(),
    sigTotCou(), sigElCou(), sigXB(), sigAX(), sigXX(), sigAXB(), idA(),
    idB(), tryCoulomb(), chgSgn(), tAbsMin(), lambda(), phaseCst(),
    particleDataPtr(), rndmPtr() {}

  // Destructor.
  virtual ~SigmaTotAux() {};

  // Store pointers and initialize data members.
  virtual void init( Info* , Settings&, ParticleData*, Rndm* ) {}

  // Calculate integrated total/elastic cross sections.
  // Usage: calcTotEl( idAin, idBin, sIn, mAin, mBin).
  virtual bool calcTotEl( int , int , double , double , double ) {return true;}

  // Store total and elastic cross section properties.
  bool   isExpEl, hasCou;
  double sigTot, rhoOwn, sigEl, bEl, sigTotCou, sigElCou;

  // Differential elastic cross section, d(sigma_el) / dt.
  virtual double dsigmaEl( double , bool = false, bool = true) {return 0.;}

  // Calculate integrated diffractive cross sections.
  // Usage: calcDiff(  idAin, idBin, sIn, mAin, mBin).
  virtual bool calcDiff(  int , int , double , double , double ) {return true;}

  // Store diffractive cross sections.
  double sigXB, sigAX, sigXX, sigAXB;

  // Differential single diffractive cross section,
  // xi * d(sigma_SD) / (dxi dt).
  virtual double dsigmaSD( double , double, bool = true, int = 0) {return 0.;}

  // Possibility to separate xi and t choices for diffraction.
  virtual bool splitDiff() {return false;}

  // Differential double diffractive cross section,
  // xi1 * xi2 * d(sigma_DD) / (dxi1 dxi2 dt).
  virtual double dsigmaDD( double , double , double , int = 0) {return 0.;}

  // Differential central diffractive cross section,
  // xi1 * xi2 * d(sigma_CD) / (dxi1 dxi2 dt1 dt2).
  virtual double dsigmaCD( double , double , double , double, int ) {
    return 0.;}

  // Minimal central diffractive mass.
  virtual double mMinCD() {return 1.;}

  // Standard methods to find t range of a 2 -> 2 process
  // and to check whether a given t value is in that range.
  pair<double,double> tRange( double sIn, double s1In, double s2In,
    double s3In,  double s4In) {
    double lambda12 = pow2( sIn - s1In - s2In) - 4. * s1In * s2In;
    double lambda34 = pow2( sIn - s3In - s4In) - 4. * s3In * s4In;
    if (lambda12 < 0. || lambda34 < 0.) return make_pair( 0., 0.);
    double tLow = -0.5 * (sIn - (s1In + s2In + s3In + s4In) + (s1In - s2In)
      * (s3In - s4In) / sIn + sqrtpos(lambda12 *  lambda34) / sIn);
    double tUpp = ( (s3In - s1In) * (s4In - s2In) + (s1In + s4In - s2In - s3In)
      * (s1In * s4In - s2In * s3In) / sIn ) / tLow;
    return make_pair( tLow, tUpp); }
  bool tInRange( double tIn, double sIn, double s1In, double s2In,
    double s3In,  double s4In) {
    pair<double, double> tRng = tRange( sIn, s1In, s2In, s3In, s4In);
    return (tIn > tRng.first && tIn < tRng.second); }

  // Commonly used proton form factor.
  double pFormFac(double tIn) {return (4. * SPROTON - 2.79 * tIn)
    / ((4. * SPROTON - tIn) * pow2(1. - tIn / 0.71)); }

protected:

  // Constants: could only be changed in the code itself.
  static const int    NPOINTS;
  static const double ALPHAEM, HBARC2, CONVERTEL, MPROTON, SPROTON, MPION,
                      SPION, GAMMAEUL, TABSREF, TABSMAX, MINSLOPEEL;

  // Initialization data, normally only set once.
  int    idA, idB;

  // Add Coulomb corrections to the elastic cross section.
  bool           tryCoulomb;
  double         chgSgn, tAbsMin, lambda, phaseCst;
  virtual bool   initCoulomb(Settings& settings,
    ParticleData* particleDataPtrIn);
  virtual bool   addCoulomb();
  virtual double dsigmaElCoulomb(double t);

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

};

//==========================================================================

// The SigmaTotal class contains parametrizations of total, elastic and
// diffractive cross sections, and of the respective slope parameter.

class SigmaTotal {

public:

  // Constructor.
  SigmaTotal() : isCalc(false), ispp(), modeTotEl(), modeTotElNow(),
    modeDiff(), modeDiffNow(), idAbsA(), idAbsB(), s(), sigND(),
    sigTotElPtr(NULL), sigDiffPtr(NULL), infoPtr(), settingsPtr(),
    particleDataPtr(), rndmPtr() {};

  // Destructor.
  virtual ~SigmaTotal() { if (sigTotElPtr) delete sigTotElPtr;
    if (sigDiffPtr) delete sigDiffPtr; }

  // Store pointers and initialize data members.
  void   init( Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn);

  // Calculate, or recalculate for new beams or new energy.
  bool   calc( int idA, int idB, double eCM);

  // Confirm that initialization worked.
  bool   hasSigmaTot()  {return isCalc;}

  // Total and integrated elastic cross sections.
  double sigmaTot()     {return sigTotElPtr->sigTotCou;}
  double rho()          {return sigTotElPtr->rhoOwn;}
  double sigmaEl()      {return sigTotElPtr->sigElCou;}
  bool   bElIsExp()     {return sigTotElPtr->isExpEl;}
  double bSlopeEl()     {return sigTotElPtr->bEl;}
  bool   hasCoulomb()   {return sigTotElPtr->hasCou;}

  // Total elastic cross section.
  bool calcTotEl( int idAin, int idBin, double sIn, double mAin, double mBin) {
    return sigTotElPtr->calcTotEl( idAin, idBin, sIn, mAin, mBin); }

  // Differential elastic cross section.
  double dsigmaEl( double t, bool useCoulomb = false,
    bool onlyPomerons = false) {
    return sigTotElPtr->dsigmaEl( t, useCoulomb, onlyPomerons); }

  // Integrated diffractive cross sections.
  double sigmaXB()  const {return sigDiffPtr->sigXB;}
  double sigmaAX()  const {return sigDiffPtr->sigAX;}
  double sigmaXX()  const {return sigDiffPtr->sigXX;}
  double sigmaAXB() const {return sigDiffPtr->sigAXB;}
  double sigmaND()  const {return sigND;}

  // Differential single diffractive cross section.
  double dsigmaSD( double xi, double t, bool isXB = true, int step = 0) {
    return sigDiffPtr->dsigmaSD( xi, t, isXB, step); }

  // Possibility to separate xi and t choices for diffraction.
  virtual bool splitDiff() {return sigDiffPtr->splitDiff();}

  // Differential double diffractive cross section.
  double dsigmaDD( double xi1, double xi2, double t, int step = 0) {
    return sigDiffPtr->dsigmaDD( xi1, xi2, t, step); }

  // Differential central diffractive cross section.
  double dsigmaCD( double xi1, double xi2, double t1, double t2, int step = 0)
    { return sigDiffPtr->dsigmaCD( xi1, xi2, t1, t2, step); }

  // Minimal central diffractive mass.
  double mMinCD() {return sigDiffPtr->mMinCD();}

  // Sample the VMD states for resolved photons.
  void chooseVMDstates(int idA, int idB, double eCM, int processCode);

  // Standard methods to find t range of a 2 -> 2 process
  // and to check whether a given t value is in that range.
  pair<double,double> tRange( double sIn, double s1In, double s2In,
    double s3In,  double s4In) {
    return sigDiffPtr->tRange( sIn, s1In, s2In, s3In, s4In); }
  bool tInRange( double tIn, double sIn, double s1In, double s2In,
    double s3In,  double s4In) {
    return sigDiffPtr->tInRange( tIn, sIn, s1In, s2In, s3In, s4In); }

private:

  // Constants: could only be changed in the code itself.
  static const double MMIN;

  // Initialization data, normally only set once.
  bool   isCalc, ispp;

  int    modeTotEl, modeTotElNow, modeDiff, modeDiffNow, idAbsA, idAbsB;
  double s, sigND;

  // Pointer to class that handles total and elastic cross sections.
  SigmaTotAux*  sigTotElPtr;

  // Pointer to class that handles diffractive cross sections.
  SigmaTotAux*  sigDiffPtr;

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the settings database.
  Settings*     settingsPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

};

//==========================================================================

// The SigmaTotOwn class parametrizes total, elastic and diffractive
// cross sections by user settings.

class SigmaTotOwn : public SigmaTotAux {

public:

  // Constructor.
  SigmaTotOwn() : dampenGap(), pomFlux(), s(), a0(), ap(), b0(), A1(), A2(),
    A3(), a1(), a2(), a3(), bMinDD(), ygap(), ypow(), expPygap(), mMinCDnow(),
    wtNow(), yNow(), yNow1(), yNow2(), b(), b1(), b2(), Q(), Q1(),
    Q2() {};

  // Store pointers and initialize data members.
  virtual void init( Info* , Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* );

  // Calculate integrated total/elastic cross sections.
  virtual bool calcTotEl( int idAin, int idBin, double , double , double);

  // Differential elastic cross section.
  virtual double dsigmaEl( double t, bool useCoulomb = false, bool = true);

  // Calculate integrated diffractive cross sections.
  virtual bool calcDiff(  int , int , double sIn, double , double ) {
    s = sIn; return true;}

  // Differential single diffractive cross section.
  virtual double dsigmaSD( double xi, double t, bool = true, int = 0);

  // Differential double diffractive cross section.
  virtual double dsigmaDD( double xi1, double xi2, double t, int = 0);

  // Differential central diffractive cross section.
  virtual double dsigmaCD( double xi1, double xi2, double t1, double t2,
    int = 0);

  // Minimal central diffractive mass.
  virtual double mMinCD() {return mMinCDnow;}

private:

  // Initialization data, normally only set once.
  bool   dampenGap;
  int    pomFlux;
  double s, a0, ap, b0, A1, A2, A3, a1, a2, a3, bMinDD, ygap, ypow, expPygap,
         mMinCDnow, wtNow, yNow, yNow1, yNow2, b, b1, b2, Q, Q1, Q2;

};

//==========================================================================

// The SigmaSaSDL class parametrizes total, elastic and diffractive
// cross sections according to Schuler and Sjostrand, starting from
// Donnachie and Landshoff.

class SigmaSaSDL : public SigmaTotAux {

public:

  // Constructor.
  SigmaSaSDL() : doDampen(), zeroAXB(), swapped(), sameSign(), idAbsA(),
    idAbsB(), iProc(), iHadA(), iHadB(), iHadAtmp(), iHadBtmp(), iProcVP(),
    iProcVV(), s(), mA(), mB(), bA(), bB(), maxXBOwn(), maxAXOwn(), maxXXOwn(),
    maxAXBOwn(), epsSaS(), sigmaPomP(), mPomP(), pPomP(), sigAXB2TeV(),
    mMin0(), cRes(), mRes0(), mMinCDnow(), alP2(), s0(), mMinXB(), mMinAX(),
    mMinAXB(), mResXB(), mResAX(), sResXB(), sResAX(), wtNow(), mAtmp(),
    mBtmp(), multVP(), multVV(), infoPtr() {};

  // Store pointers and initialize data members.
  virtual void init( Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* );

  // Calculate integrated total/elastic cross sections.
  virtual bool calcTotEl( int idAin, int idBin, double sIn, double mAin,
    double mBin);

  // Differential elastic cross section.
  virtual double dsigmaEl( double t, bool useCoulomb = false, bool = true);

  // Calculate integrated diffractive cross sections.
  virtual bool calcDiff(  int idAin, int idBin, double sIn, double mAin,
    double mBin);

  // Differential single diffractive cross section.
  virtual double dsigmaSD( double xi, double t, bool isXB, int = 0);

  // Differential double diffractive cross section.
  virtual double dsigmaDD( double xi1, double xi2, double t, int = 0);

  // Differential central diffractive cross section.
  virtual double dsigmaCD( double xi1, double xi2, double t1, double t2,
    int = 0);

  // Minimal central diffractive mass.
  virtual double mMinCD() {return mMinCDnow;}

private:

  // Constants: could only be changed in the code itself.
  static const int    IHADATABLE[], IHADBTABLE[], ISDTABLE[], IDDTABLE[], NVMD;
  static const double EPSILON, ETA, X[], Y[], BETA0[], BHAD[], ALPHAPRIME,
                      CONVERTSD, CONVERTDD, VMDMASS[4], GAMMAFAC[4],
                      CSD[10][8], CDD[10][9];

  // Initialization data, normally only set once, and result of calculation.
  bool   doDampen, zeroAXB, swapped, sameSign;
  int    idAbsA, idAbsB, iProc, iHadA, iHadB, iHadAtmp[4],
         iHadBtmp[4], iProcVP[4], iProcVV[4][4];
  double s, mA, mB, bA, bB, maxXBOwn, maxAXOwn, maxXXOwn, maxAXBOwn, epsSaS,
         sigmaPomP, mPomP, pPomP, sigAXB2TeV, mMin0,  cRes, mRes0, mMinCDnow,
         alP2, s0, mMinXB,  mMinAX, mMinAXB, mResXB, mResAX, sResXB,
         sResAX, wtNow, mAtmp[4], mBtmp[4], multVP[4], multVV[4][4];

  // Find combination of incoming beams.
  bool findBeamComb( int idAin, int idBin, double mAin, double mBin);

  // Pointer to various information on the generation.
  Info*         infoPtr;

};

//==========================================================================

// The SigmaMBR class parametrizes total and elastic cross sections
// according to the Minimum Bias Rockefeller (MBR) model.

class SigmaMBR : public SigmaTotAux {

public:

  // Constructor.
  SigmaMBR() : s(), sigSD(), sigDD(), sigCD(), eps(), alph(), beta0gev(),
    beta0mb(), sigma0mb(), sigma0gev(), m2min(), dyminSDflux(),
    dyminDDflux(), dyminCDflux(), dyminSD(), dyminDD(), dyminCD(),
    dyminSigSD(), dyminSigDD(), dyminSigCD(), a1(), a2(), b1(), b2(),
    sdpmax(), ddpmax(), dpepmax() {};

  // Initialize data members.
  virtual void init( Info* , Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* );

  // Calculate integrated total/elastic cross sections.
  virtual bool calcTotEl( int idAin, int idBin, double sIn, double , double );

  // Differential elastic cross section.
  virtual double dsigmaEl(double t, bool useCoulomb = false, bool = true);

  // Calculate integrated diffractive cross sections.
  virtual bool calcDiff(  int , int , double sIn, double , double );

  // In MBR choice of xi and t values are separated.
  virtual bool splitDiff() {return true;}

  // Differential single diffractive cross section.
  virtual double dsigmaSD( double xi, double t, bool = true, int step = 0);

  // Differential double diffractive cross section.
  virtual double dsigmaDD( double xi1, double xi2, double t, int step = 0);

  // Differential central diffractive cross section.
  virtual double dsigmaCD( double xi1, double xi2, double t1, double t2,
    int step = 0);

  // Minimal central diffractive mass.
  virtual double mMinCD() {return sqrt(m2min);}

private:

  // Integration of MBR cross sections and form factor approximation.
  static const int    NINTEG, NINTEG2;
  static const double FFA1, FFA2, FFB1, FFB2;

  // Parameters of MBR model.
  double s, sigSD, sigDD, sigCD, eps, alph, beta0gev,
         beta0mb, sigma0mb, sigma0gev, m2min, dyminSDflux, dyminDDflux,
         dyminCDflux, dyminSD, dyminDD, dyminCD, dyminSigSD, dyminSigDD,
         dyminSigCD, a1, a2, b1, b2, sdpmax, ddpmax, dpepmax;

  // The error function erf(x) should normally be in your math library,
  // but if not uncomment this simple parametrization by Sergei Winitzki.
  //double erf(double x) { double x2 = x * x; double kx2 = 0.147 * x2;
  //  double tmp = sqrt(1. - exp(-x2 * (4./M_PI + kx2) / (1. + kx2)));
  //  return ((x >= 0.) ? tmp : -tmp); }

};

//==========================================================================

// The SigmaABMST class parametrizes total and elastic cross sections
// according to Appleby, Barlow, Molson, Serluca and Toader (ABMST).

class SigmaABMST : public SigmaTotAux {

public:

  // Constructor.
  SigmaABMST() : ispp(), dampenGap(), useBMin(), modeSD(), modeDD(), modeCD(),
    s(), facEl(), m2minp(), m2minm(), alp0(), alpt(), s0(), c0(), ygap(),
    ypow(), expPygap(), multSD(), powSD(), multDD(), powDD(), multCD(),
    powCD(), mMinCDnow(), bMinSD(), bMinDD(), bMinCD() {};

  // Initialize data members.
  virtual void init( Info* , Settings& settings, ParticleData* ,
    Rndm* rndmPtrIn);

  // Calculate integrated total/elastic cross sections.
  virtual bool calcTotEl( int idAin, int idBin, double sIn, double , double );

  // Differential elastic cross section.
  virtual double dsigmaEl( double t, bool useCoulomb = false,
    bool onlyPomerons = false) {
    return facEl * pow2(abs(amplitude( t, useCoulomb, onlyPomerons)));}

  // Calculate integrated diffractive cross sections.
  virtual bool calcDiff(  int idAin, int idBin, double sIn, double , double );

  // Differential single diffractive cross section.
  virtual double dsigmaSD( double xi, double t, bool = true , int = 0);

  // Differential double diffractive cross section.
  virtual double dsigmaDD( double xi1, double xi2, double t, int = 0);

  // Differential central diffractive cross section.
  virtual double dsigmaCD( double xi1, double xi2, double t1, double t2,
    int = 0);

  // Minimal central diffractive mass.
  virtual double mMinCD() {return mMinCDnow;}

private:

  // Constants: could only be changed in the code itself.
  static const bool   MCINTDD;
  static const int    NPOINTSTSD, NPOINTSTDD, NPOINTMCDD, NPOINTMCCD;
  static const double EPSI[], ALPP[], NORM[], SLOPE[], FRACS[], TRIG[],
                      LAM2P, BAPPR[], LAM2FF, MRES[4], WRES[4], CRES[4],
                      AFAC[4], BFAC[4], CFAC[4], PPP[4], EPS[2], APR[2],
                      CPI[6], CNST[5], XIDIVSD, DXIRAWSD, DLNXIRAWSD,
                      XIDIVDD, DXIRAWDD, DLNXIRAWDD, BMCINTDD, BMCINTCD;

  // Initialization data, normally only set once.
  bool   ispp, dampenGap, useBMin;
  int    modeSD, modeDD, modeCD;
  double s, facEl, m2minp, m2minm, alp0[2], alpt[3], s0, c0,
         ygap, ypow, expPygap, multSD, powSD, multDD, powDD,
         multCD, powCD, mMinCDnow, bMinSD, bMinDD, bMinCD;

  // The scattering amplitude, from which cross sections are derived.
  complex amplitude( double t, bool useCoulomb = false,
    bool onlyPomerons = false) ;

  // Auxiliary method for repetitive part of amplitude.
  complex sModAlp( double sMod, double alpha) {
    return exp(complex( 0., -0.5 * M_PI * alpha)) * pow( sMod, alpha); }

  // Core differential single diffractive cross section.
  virtual double dsigmaSDcore(double xi, double t);

  // xi * d(sigma_SD) / (dxi dt) integrated over [t_min, t_max].
  double dsigmaSDintT( double xi, double tMinIn, double tMaxIn);

  // d(sigma_SD) / (dxi dt) integrated over [xi_min, xi_max] * [t_min, t_max].
  double dsigmaSDintXiT( double xiMinIn, double xiMaxIn, double tMinIn,
    double tMaxIn );

  // d(sigma_DD) / (dxi1 dxi2 dt) integrated with Monte Carlo sampling.
  double dsigmaDDintMC();

  // xi1 * xi2 * d(sigma_DD) / (dxi1 dxi2 dt) integrated over t.
  double dsigmaDDintT( double xi1, double xi2, double tMinIn, double tMaxIn);

  // xi1 * d(sigma_DD) / (dxi1 dxi2 dt) integrated over xi2 and t.
  double dsigmaDDintXi2T( double xi1, double xi2MinIn, double xi2MaxIn,
    double tMinIn, double tMaxIn);

  // d(sigma_DD) / (dxi1 dxi2 dt) integrated over xi1, xi2 and t.
  double dsigmaDDintXi1Xi2T( double xi1MinIn, double xi1MaxIn,
    double xi2MinIn, double xi2MaxIn, double tMinIn, double tMaxIn);

  // d(sigma_CD) / (dxi1 dxi2 dt1 dt2) integrated with Monte Carlo sampling.
  double dsigmaCDintMC();

};

//==========================================================================

// The SigmaRPP class parametrizes total and elastic cross sections
// according to the fit in Review of Particle Physics 2014.

class SigmaRPP : public SigmaTotAux {

public:

  // Constructor.
  SigmaRPP() : ispp(), s(), facEl() {};

  // Initialize data members.
  virtual void init( Info* , Settings& settings, ParticleData* , Rndm* ) {
    tryCoulomb = settings.flag("SigmaElastic:Coulomb");
    tAbsMin = settings.parm("SigmaElastic:tAbsMin"); }

  // Calculate integrated total/elastic cross sections.
  virtual bool calcTotEl( int idAin, int idBin, double sIn, double , double );

  // Differential elastic cross section.
  virtual double dsigmaEl( double t, bool useCoulomb = false, bool = true) {
    return facEl * pow2(abs(amplitude( t, useCoulomb))); }

private:

  // Constants: could only be changed in the code itself.
  static const double EPS1[], ALPP[], NORM[], BRPP[], KRPP[], LAM2FF;

  // Initialization data, normally only set once, and result of calculation.
  bool   ispp;
  double s, facEl;

  // The scattering amplitude, from which cross sections are derived.
  complex amplitude( double t, bool useCoulomb = false) ;

  // Auxiliary method for repetitive part of amplitude.
  complex sModAlp( double sMod, double alpha) {
    return exp(complex( 0., -0.5 * M_PI * alpha)) * pow( sMod, alpha); }

  // Auxiliary methods for Bessel J0 and J1 functions.
  complex besJ0( complex x);
  complex besJ1( complex x);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SigmaTotal_H
