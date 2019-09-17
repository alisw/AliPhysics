// SigmaTotal.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the SigmaTotal class,
// and the auxiliary classes SigmaTotAux, SigmaTotOwn, SigmaSaSDL, SigmaMBR,
// SigmaABMST and SigmaRPP.

#include "Pythia8/SigmaTotal.h"

namespace Pythia8 {

//==========================================================================

// The SigmaTotAux class.
// Base class for the individual implementations.

//--------------------------------------------------------------------------

// alpha_em(0).
const double SigmaTotAux::ALPHAEM    = 0.00729353;

// Conversion coefficients = 1/(16pi) * (mb <-> GeV^2).
const double SigmaTotAux::HBARC2     = 0.38938;
const double SigmaTotAux::CONVERTEL  = 0.0510925;

// Proton and pion masses, and their squares. Euler's constant.
const double SigmaTotAux::MPROTON    = 0.9382720;
const double SigmaTotAux::SPROTON    = 0.8803544;
const double SigmaTotAux::MPION      = 0.1349766;
const double SigmaTotAux::SPION      = 0.0182187;
const double SigmaTotAux::GAMMAEUL   = 0.577215665;

// Reference scale for nominal definition of t slope.
const double SigmaTotAux::TABSREF    = 2e-3;

// Numerical integration parameters.
const int    SigmaTotAux::NPOINTS    = 1000;
const double SigmaTotAux::TABSMAX    = 1.;
const double SigmaTotAux::MINSLOPEEL = 10.;

//--------------------------------------------------------------------------

// Initialize Coulomb correction parameters.

bool SigmaTotAux::initCoulomb(Settings& settings,
  ParticleData* particleDataPtrIn) {

  // Save pointer to particle database.
  particleDataPtr = particleDataPtrIn,

  // User-set values for total and elastic cross section.
  tryCoulomb = settings.flag("SigmaElastic:Coulomb");
  rhoOwn     = settings.parm("SigmaElastic:rho");
  tAbsMin    = settings.parm("SigmaElastic:tAbsMin");
  lambda     = settings.parm("SigmaElastic:lambda");
  phaseCst   = settings.parm("SigmaElastic:phaseConst");

  return true;

}

//--------------------------------------------------------------------------

// Add Coulomb corrections to the elastic and total cross sections.

bool SigmaTotAux::addCoulomb() {

  // Trivial case when there should be no Coulomb contribution.
  hasCou    = false;
  sigTotCou = sigTot;
  sigElCou  = sigEl;

  // Relative sign (or zero) for Coulomb term in elastic scattering.
  int iChA  = particleDataPtr->chargeType(idA);
  int iChB  = particleDataPtr->chargeType(idB);
  chgSgn    = 0.;
  if (iChA * iChB > 0) chgSgn =  1.;
  if (iChA * iChB < 0) chgSgn = -1.;

  // Done if no Coulomb corrections.
  if (!tryCoulomb || iChA * iChB == 0) return false;

  // Reduce hadronic part of elastic cross section by tMin cut.
  sigElCou    = sigEl * exp( - bEl * tAbsMin);
  if (tAbsMin < 0.9 * TABSMAX) {

    // Loop through t range according to dt/t^2.
    double sigCou = 0.;
    double sigInt = 0.;
    double xRel, tAbs, form2, phase;
    for (int i = 0; i < NPOINTS; ++i) {
      xRel    = (i + 0.5) / NPOINTS;
      tAbs    = tAbsMin * TABSMAX / (tAbsMin + xRel * (TABSMAX - tAbsMin));

      // Evaluate Coulomb and interference terms.
      form2   = pow4(lambda/(lambda + tAbs));
      sigCou += pow2(form2);
      phase   = chgSgn * ALPHAEM * (-phaseCst - log(0.5 * bEl * tAbs));
      sigInt += form2 * exp(-0.5 * bEl * tAbs) * tAbs
              * (rhoOwn * cos(phase) + sin(phase));
    }

    // Include common factors to give new elastic and total cross sections.
    sigCou   *= pow2(ALPHAEM) / (4. * CONVERTEL * tAbsMin) ;
    sigInt   *= - chgSgn * ALPHAEM * sigTot / tAbsMin;
    sigElCou += (sigCou + sigInt) / NPOINTS;
    hasCou    = true;
  }
  sigTotCou   = sigTot - sigEl + sigElCou;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Coulomb contribution to the differential elastic cross sections.

double SigmaTotAux::dsigmaElCoulomb( double t) {

  // Coulomb contribution and interference.
  double form2 = pow4(lambda/(lambda - t));
  double phase = chgSgn * ALPHAEM * (-phaseCst - log(-0.5 * bEl * t));
  return pow2(chgSgn * ALPHAEM * form2) / (4. * CONVERTEL * t*t)
         - chgSgn * ALPHAEM * form2 * sigTot * exp(0.5 * bEl * t)
         * (rhoOwn * cos(phase) + sin(phase)) / (-t);

}

//==========================================================================

// The SigmaTotal class.
// Parametrizes total, elastic and diffractive cross sections,
// making use of the other clases below.

//--------------------------------------------------------------------------

// Definitions of static variables.

// Minimum threshold below which no cross sections will be defined.
const double SigmaTotal::MMIN  = 2.;

//--------------------------------------------------------------------------

// Store pointer to Info and initialize data members.

void SigmaTotal::init( Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn) {

  // Store pointers.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  settingsPtr     = &settings;
  rndmPtr         = rndmPtrIn;

  // Choice of mode.
  modeTotEl  = settings.mode("SigmaTotal:mode");
  modeDiff   = settings.mode("SigmaDiffractive:mode");

}

//--------------------------------------------------------------------------

// Calculate, or recalculate for new beams or new energy.

bool SigmaTotal::calc(int idA, int idB, double eCM) {

  // Initial values.
  isCalc = ispp = false;
  s = eCM * eCM;

  // Find hadron masses and check that energy is enough.
  // For mesons use the corresponding vector meson masses.
  idAbsA = abs(idA);
  idAbsB = abs(idB);
  int idModA = (idAbsA < 100 || idAbsA > 1000) ? idAbsA : 10 * (idAbsA/10) + 3;
  int idModB = (idAbsB < 100 || idAbsB > 1000) ? idAbsB : 10 * (idAbsB/10) + 3;
  if (idAbsA == 22) idModA = 113;
  if (idAbsB == 22) idModB = 113;
  if (idAbsA == 990) idModA = idAbsA;
  if (idAbsB == 990) idModB = idAbsB;
  double mA  = particleDataPtr->m0(idModA);
  double mB  = particleDataPtr->m0(idModB);
  if (eCM < mA + mB + MMIN) {
    infoPtr->errorMsg("Error in SigmaTotal::calc: too low energy");
    return false;
  }

  // Most options only work for pp/ppbar, so may need to modify choice.
  // Treat a neutron like a proton (except no Coulomb term).
  modeTotElNow = modeTotEl;
  modeDiffNow  = modeDiff;
  if (idAbsA == 2112) idAbsA = 2212;
  if (idAbsB == 2112) idAbsB = 2212;
  if (idAbsA != 2212 || idAbsB != 2212) {
    modeTotElNow = min(1, modeTotEl);
    modeDiffNow  = min(1, modeDiff);
  }
  ispp = (idAbsA == 2212 && idAbsB == 2212 && idA * idB > 0);

  // Set up pointer to class that handles total and elastic cross sections.
  if (sigTotElPtr) delete sigTotElPtr;
  if      (modeTotElNow == 0) sigTotElPtr = new SigmaTotOwn();
  else if (modeTotElNow == 1) sigTotElPtr = new SigmaSaSDL();
  else if (modeTotElNow == 2) sigTotElPtr = new SigmaMBR();
  else if (modeTotElNow == 3) sigTotElPtr = new SigmaABMST();
  else                        sigTotElPtr = new SigmaRPP();

  // Initialize and calculate for selected total/elastic class.
  sigTotElPtr->init( infoPtr, *settingsPtr, particleDataPtr, rndmPtr);
  if ( !sigTotElPtr->calcTotEl( idA, idB, s, mA, mB) ) return false;

  // Set up pointer to class that handles diffractive cross sections.
  if (sigDiffPtr) delete sigDiffPtr;
  if      (modeDiffNow == 0) sigDiffPtr = new SigmaTotOwn();
  else if (modeDiffNow == 1) sigDiffPtr = new SigmaSaSDL();
  else if (modeDiffNow == 2) sigDiffPtr = new SigmaMBR();
  else                       sigDiffPtr = new SigmaABMST();

  // Initialize and calculate for selected diffractive class.
  if (sigDiffPtr != sigTotElPtr)
    sigDiffPtr->init( infoPtr, *settingsPtr, particleDataPtr, rndmPtr);
  if ( !sigDiffPtr->calcDiff( idA, idB, s, mA, mB) ) return false;

  // Inelastic nondiffractive by unitarity.
  double sigTotTmp = sigTotElPtr->sigTot;
  sigND = sigTotTmp - sigTotElPtr->sigEl - sigDiffPtr->sigXB
        - sigDiffPtr->sigAX - sigDiffPtr->sigXX - sigDiffPtr->sigAXB;
  if (sigND < 0.) {
    infoPtr->errorMsg("Error in SigmaTotal::init: sigND < 0");
    return false;
  } else if (sigND < 0.4 * sigTotTmp) infoPtr->errorMsg("Warning in "
    "SigmaTotal::init: sigND suspiciously low");

  // Done.
  isCalc = true;
  return true;

}

//--------------------------------------------------------------------------

// Sample the VMD states for resolved photons.

void SigmaTotal::chooseVMDstates(int idA, int idB, double eCM,
  int processCode) {

  // Constants and initial values.
  double gammaFac[4] = {2.2, 23.6, 18.4, 11.5};
  double alphaEM     = 0.00729353;
  double idVMD[4]    = {113, 223, 333, 443};
  double pVP[4]      = {0.};
  double pVV[4][4]   = {{0.}};
  double pSum        = 0.;
  int nVMD           = 4;

  // Values to start with.
  pair<int, int> idAB = make_pair(idA, idB);

  // gamma-gamma.
  if (idA == 22 && idB == 22) {
    for (int i = 0; i < nVMD; ++i)
    for (int j = 0; j < nVMD; ++j) {
      // Evaluate the cross sections individually.
      calc(idVMD[i], idVMD[j], eCM);
      pVV[i][j] = pow2(alphaEM) / (gammaFac[i] * gammaFac[j]);
      if      (processCode == 101) pVV[i][j] *= sigmaND();
      else if (processCode == 102) pVV[i][j] *= sigmaEl();
      else if (processCode == 103) pVV[i][j] *= sigmaXB();
      else if (processCode == 104) pVV[i][j] *= sigmaAX();
      else if (processCode == 105) pVV[i][j] *= sigmaXX();
      pSum     += pVV[i][j];
    }

    // Choose VMD states based on relative fractions.
    double pickMode = rndmPtr->flat() * pSum;
    bool pairFound  = false;
    for (int i = 0; i < nVMD; ++i) {
      for (int j = 0; j < nVMD; ++j) {
        pickMode -= pVV[i][j];
        if (pickMode < 0.) {
          idAB = make_pair(113 + 110 * i, 113 + 110 * j);
          pairFound = true;
          break;
        }
      }
      if (pairFound) break;
    }

  // gamma + p.
  } else if (idA == 22 && idB == 2212) {
    for (int i = 0; i < nVMD; ++i) {
      // Evaluate the cross sections individually.
      calc(idVMD[i], 2212, eCM);
      pVP[i] = alphaEM / gammaFac[i];
      if      (processCode == 101) pVP[i] *= sigmaND();
      else if (processCode == 102) pVP[i] *= sigmaEl();
      else if (processCode == 103) pVP[i] *= sigmaXB();
      else if (processCode == 104) pVP[i] *= sigmaAX();
      else if (processCode == 105) pVP[i] *= sigmaXX();
      pSum     += pVP[i];
    }

    // Choose VMD state based on relative fractions.
    double pickMode = rndmPtr->flat() * pSum;
    for (int i = 0; i < nVMD; ++i){
      pickMode -= pVP[i];
      if (pickMode < 0.){
         idAB = make_pair(113 + 110 * i, 2212);
         break;
      }
    }

  // p + gamma.
  } else if (idA == 2212 && idB == 22) {
    for (int i = 0; i < nVMD; ++i) {
      // Evaluate the cross sections individually.
      calc(2212, idVMD[i], eCM);
      pVP[i] = alphaEM / gammaFac[i];
      if      (processCode == 101) pVP[i] *= sigmaND();
      else if (processCode == 102) pVP[i] *= sigmaEl();
      else if (processCode == 103) pVP[i] *= sigmaXB();
      else if (processCode == 104) pVP[i] *= sigmaAX();
      else if (processCode == 105) pVP[i] *= sigmaXX();
      pSum     += pVP[i];
    }

    // Choose VMD state based on relative fractions.
    double pickMode = rndmPtr->flat() * pSum;
    for (int i = 0; i < nVMD; ++i){
      pickMode -= pVP[i];
      if (pickMode < 0.){
        idAB = make_pair(2212, 113 + 110 * i);
        break;
      }
    }
  }

  // Reset to the original cross section.
  calc(idA, idB, eCM);

  // Propagate the selected states to Info class for further usage.
  if (idAB.first == 113 || idAB.first == 223 || idAB.first == 333
      || idAB.first == 443) {
    double mA  = particleDataPtr->mSel(idAB.first);
    double scA = alphaEM / gammaFac[idAB.first/100 - 1];
    infoPtr->setVMDstateA(true, idAB.first, mA, scA);
  }
  if (idAB.second == 113 || idAB.second == 223 || idAB.second == 333
      || idAB.second == 443) {
    double mB  = particleDataPtr->mSel(idAB.second);
    double scB = alphaEM / gammaFac[idAB.second/100 - 1];
    infoPtr->setVMDstateB(true, idAB.second, mB, scB);
  }
}

//==========================================================================

// The SigmaTotOwn class.
// Parametrizes total, elastic and diffractive cross sections
// by user settings.

//--------------------------------------------------------------------------

// Initialize data members.

void SigmaTotOwn::init(Info* , Settings& settings,
  ParticleData* particleDataPtrIn, Rndm*) {

  // Main user-set values for total and elastic cross sections.
  sigTot  = settings.parm("SigmaTotal:sigmaTot");
  sigEl   = settings.parm("SigmaTotal:sigmaEl");
  bEl     = settings.parm("SigmaElastic:bSlope");

  // Initialize parameters for Coulomb corrections to elastic scattering.
  initCoulomb( settings, particleDataPtrIn);

  // User-set values for diffractive cross sections.
  sigXB   = settings.parm("SigmaTotal:sigmaXB");
  sigAX   = settings.parm("SigmaTotal:sigmaAX");
  sigXX   = settings.parm("SigmaTotal:sigmaXX");
  sigAXB  = settings.parm("SigmaTotal:sigmaAXB");

  // Set diffraction parameters.
  pomFlux = settings.mode("SigmaDiffractive:PomFlux");

  // Set up Pomeron flux constants, see HardDiffraction::init.
  a0      = 1. + settings.parm("SigmaDiffractive:PomFluxEpsilon");
  ap      = settings.parm("SigmaDiffractive:PomFluxAlphaPrime");

  // Schuler-Sjöstrand.
  if (pomFlux == 1) {
    b0    = 2.3;

  // Bruni-Ingelman.
  } else if (pomFlux == 2) {
    A1    = 6.38;
    A2    = 0.424;
    a1    = 8.;
    a2    = 3.;

  // Streng-Berger.
  } else if (pomFlux == 3) {
    a1    = 4.7;

  // Donnachie-Landshoff.
  } else if (pomFlux == 4) {
    A1    = 0.27;
    a1    = 8.38;
    A2    = 0.56;
    a2    = 3.78;
    A3    = 0.18;
    a3    = 1.36;

  // MBR.
  } else if (pomFlux == 5) {
    A1    = 0.9;
    a1    = 4.6;
    A2    = 0.1;
    a2    = 0.6;
    a0    = 1. + settings.parm("SigmaDiffractive:MBRepsilon");
    ap    = settings.parm("SigmaDiffractive:MBRalpha");

  // H1 Fit A, B.
  } else if (pomFlux == 6 || pomFlux == 7) {
    ap    = 0.06;
    b0    = 5.5;
    a0    = (pomFlux == 6) ? 1.1182 : 1.1110;
  }

  // b_min for double diffraction, suppression of small gaps, minimal CD mass.
  bMinDD    = settings.parm("SigmaDiffractive:OwnbMinDD");
  dampenGap = settings.flag("SigmaDiffractive:OwndampenGap");
  ygap      = settings.parm("SigmaDiffractive:Ownygap");
  ypow      = settings.parm("SigmaDiffractive:Ownypow");
  expPygap  = exp(ypow * ygap);
  mMinCDnow = settings.parm("SigmaDiffractive:OwnmMinCD");

}

//--------------------------------------------------------------------------

// With total and elastic cross section already set, only add Coulomb.

bool SigmaTotOwn::calcTotEl( int idAin, int idBin, double , double , double) {

  // Save some input.
  idA       = idAin;
  idB       = idBin;
  isExpEl   = true;

  // Possibly allow Coulomb correction + interference.
  addCoulomb();

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Differential elastic cross section.

double SigmaTotOwn::dsigmaEl( double t, bool useCoulomb, bool ) {

  // Hadronic contribution: simple exponential.
  double dsig = CONVERTEL * pow2(sigTot) * (1. + pow2(rhoOwn)) * exp(bEl * t);

  // Possibly add Coulomb contribution and interference.
  if (useCoulomb && hasCou) dsig += dsigmaElCoulomb(t);

  // Done.
  return dsig;

}

//--------------------------------------------------------------------------

// Differential single diffractive cross sections as provided by the
// x-weighted Pomeron flux, see HardDiffraction::xfPom for details.

double SigmaTotOwn::dsigmaSD(double xi, double t, bool , int ) {

  // Common setup.
  wtNow   = 1.;
  yNow    = -log(xi);

  // Schuler-Sjöstrand.
  if (pomFlux == 1) {
    b     = 2. * b0 + 2. * ap * yNow;
    wtNow = exp(b * t);

  // Bruni-Ingelman.
  } else if (pomFlux == 2) {
    wtNow = A1 * exp(a1 * t) + A2 * exp(a2 * t);

  // Streng-Berger.
  } else if (pomFlux == 3) {
    b     = a1 + 2. * ap * yNow;
    wtNow = pow( xi, 2. - 2. * a0) * exp(b * t);

  // Donnachie-Landshoff.
  } else if (pomFlux == 4) {
    Q     = 2. * ap * yNow;
    wtNow = pow( xi, 2. - 2. * a0) * ( A1 * exp((Q + a1) * t)
           + A2 * exp((Q + a2) * t) + A3 * exp((Q + a3) * t) );

  // MBR.
  } else if (pomFlux == 5) {
    Q     = 2. * ap * yNow;
    wtNow = pow( xi, 2. - 2. * a0)
          * (A1 * exp((Q + a1) * t) + A2 * exp((Q + a2) * t) );

  // H1 Fit A, B.
  } else if (pomFlux == 6 || pomFlux == 7) {
    b     = b0 + 2. * ap * yNow;
    wtNow = pow( xi, 2. - 2. * a0) * exp(b * t);
  }
  // Optionally dampen with 1 / (1 + exp( -p * (y - y_gap))).
  if (dampenGap) wtNow /= 1. + expPygap * pow( xi, ypow);

  // Done.
  return wtNow;

}

//--------------------------------------------------------------------------

// Differential double diffractive cross sections.
// Naive non-authorized generalization of the simple Pomeron fluxes above.

double SigmaTotOwn::dsigmaDD(double xi1, double xi2, double t, int ) {

  // Common setup.
  wtNow   = 1.;
  yNow    = - log(xi1 * xi2 * s / SPROTON);

  // Schuler-Sjöstrand.
  if (pomFlux == 1) {
    b     = max( bMinDD, 2. * ap * yNow);
    wtNow = exp(b * t);

  // Bruni-Ingelman.
  } else if (pomFlux == 2) {
    wtNow = A1 * exp(a1 * t) + A2 * exp(a2 * t);

  // Streng-Berger.
  } else if (pomFlux == 3) {
    b     = max( bMinDD, 2. * ap * yNow);
    wtNow = pow( xi1 * xi2, 2. - 2. * a0) * exp(b * t);

  // Donnachie-Landshoff.
  } else if (pomFlux == 4) {
    Q     = max( bMinDD, 2. * ap * yNow);
    wtNow = pow( xi1 * xi2, 2. - 2. * a0) * exp(Q * t);

  // MBR.
  } else if (pomFlux == 5) {
    Q     = max( bMinDD, 2. * ap * yNow);
    wtNow = pow( xi1 * xi2, 2. - 2. * a0) * exp(Q * t);

  // H1 Fit A, B.
  } else if (pomFlux == 6 || pomFlux == 7) {
    b     = max( bMinDD, 2. * ap * yNow);
    wtNow = pow( xi1 * xi2, 2. - 2. * a0) * exp(b * t);
  }

  // Optionally dampen with 1 / (1 + exp( -p * (y - y_gap))).
  if (dampenGap) wtNow /= 1. + expPygap * pow( xi1 * xi2 * s / SPROTON, ypow);

  // Done.
  return wtNow;

}

//--------------------------------------------------------------------------

// Differential central diffractive cross section.
// Naive non-authorized generalization of the simple Pomeron fluxes above.

double SigmaTotOwn::dsigmaCD( double xi1, double xi2, double t1, double t2,
  int ) {

  // Common setup.
  wtNow   = 1.;
  yNow1   = -log(xi1);
  yNow2   = -log(xi2);

  // Schuler-Sjöstrand.
  if (pomFlux == 1) {
    b1    = 2. * b0 + 2. * ap * yNow1;
    b2    = 2. * b0 + 2. * ap * yNow2;
    wtNow = exp(b1 * t1 + b2 * t2);

  // Bruni-Ingelman.
  } else if (pomFlux == 2) {
    wtNow = (A1 * exp(a1 * t1) + A2 * exp(a2 * t1))
          * (A1 * exp(a1 * t2) + A2 * exp(a2 * t2));

  // Streng-Berger.
  } else if (pomFlux == 3) {
    b1    = a1 + 2. * ap * yNow1;
    b2    = a1 + 2. * ap * yNow2;
    wtNow = pow( xi1 * xi2, 2. - 2. * a0) * exp(b1 * t1 + b2 * t2);

  // Donnachie-Landshoff.
  } else if (pomFlux == 4) {
    Q1    = 2. * ap * yNow1;
    Q2    = 2. * ap * yNow2;
    wtNow = pow( xi1 * xi2, 2. - 2. * a0) * ( A1 * exp((Q1 + a1) * t1)
          + A2 * exp((Q1 + a2) * t1) + A3 * exp((Q1 + a3) * t1) )
          * ( A1 * exp((Q2 + a1) * t2) + A2 * exp((Q2 + a2) * t2)
          + A3 * exp((Q2 + a3) * t2) );

  // MBR.
  } else if (pomFlux == 5) {
    Q1    = 2. * ap * yNow1;
    Q2    = 2. * ap * yNow2;
    wtNow = pow( xi1 * xi2, 2. - 2. * a0)
          * (A1 * exp((Q1 + a1) * t1) + A2 * exp((Q1 + a2) * t1) )
          * (A1 * exp((Q2 + a1) * t2) + A2 * exp((Q2 + a2) * t2) );

  // H1 Fit A, B.
  } else if (pomFlux == 6 || pomFlux == 7) {
    b1    = b0 + 2. * ap * yNow1;
    b2    = b0 + 2. * ap * yNow2;
    wtNow = pow( xi1 * xi2, 2. - 2. * a0) * exp(b1 * t1 + b2 * t2);
  }

  // Optionally dampen with 1 / (1 + exp( -p * (y - y_gap))) for both gaps.
  if (dampenGap) wtNow /= (1. + expPygap * pow( xi1, ypow))
                        * (1. + expPygap * pow( xi2, ypow));

  // Done.
  return wtNow;

}

//==========================================================================

// The SigmaSaSDL class.
// Formulae are taken from:
// G.A. Schuler and T. Sjostrand, Phys. Rev. D49 (1994) 2257,
//   Z. Phys. C73 (1997) 677
// which borrows some total cross sections from
// A. Donnachie and P.V. Landshoff, Phys. Lett. B296 (1992) 227.

// Implemented processes with their process number iProc:
// =  0 : p + p;     =  1 : pbar + p;
// =  2 : pi+ + p;   =  3 : pi- + p;     =  4 : pi0/rho0 + p;
// =  5 : phi + p;   =  6 : J/psi + p;
// =  7 : rho + rho; =  8 : rho + phi;   =  9 : rho + J/psi;
// = 10 : phi + phi; = 11 : phi + J/psi; = 12 : J/psi + J/psi.
// = 13 : gamma + p (preliminary); 14 : gamma + gamma (preliminary);
// = 15 : Pom + p (preliminary).
// For now a neutron is treated like a proton.

//--------------------------------------------------------------------------

// Definitions of static variables.
// Note that a lot of parameters are hardcoded as const here, rather
// than being interfaced for public change, since any changes would
// have to be done in a globally consistent manner. Which basically
// means a rewrite/replacement of the whole class.

// General constants in total cross section parametrization:
// sigmaTot = X * s^epsilon + Y * s^eta (pomeron + reggeon).
// Christine added Pom + p, gamma + gamma and gamma + p.
const double SigmaSaSDL::EPSILON = 0.0808;
const double SigmaSaSDL::ETA     = -0.4525;
const double SigmaSaSDL::X[] = { 21.70, 21.70, 13.63, 13.63, 13.63,
  10.01, 0.970, 8.56, 6.29, 0.609, 4.62, 0.447, 0.0434, 67.7e-3, 211e-6};
const double SigmaSaSDL::Y[] = { 56.08, 98.39, 27.56, 36.02, 31.79,
  1.51, -0.146, 13.08, -0.62, -0.060, 0.030, -0.0028, 0.00028, 129e-3,
  215e-6};

// Type of the two incoming hadrons as function of the process number:
// = 0 : p/n ; = 1 : pi/rho/omega; = 2 : phi; = 3 : J/psi.
const int SigmaSaSDL::IHADATABLE[] = { 0, 0, 1, 1, 1, 2, 3, 1, 1,
  1, 2, 2, 3};
const int SigmaSaSDL::IHADBTABLE[] = { 0, 0, 0, 0, 0, 0, 0, 1, 2,
  3, 2, 3, 3};

// Hadron-Pomeron coupling beta(t) = beta(0) * exp(b*t).
// 0 = Pom + p, 1 = Pom + pi/rho/omega, 2 = Pom + phi, 3 = Pom + J/psi.
const double SigmaSaSDL::BETA0[] = { 4.658, 2.926, 2.149, 0.208};
const double SigmaSaSDL::BHAD[]  = {   2.3,   1.4,   1.4,  0.23};

// Pomeron trajectory alpha(t) = 1 + epsilon + alpha' * t
const double SigmaSaSDL::ALPHAPRIME = 0.25;

// Conversion coefficients = 1/(16pi) * (mb <-> GeV^2) * (G_3P)^n,
// with n = 0 elastic, n = 1 single and n = 2 double diffractive.
const double SigmaSaSDL::CONVERTSD = 0.0336;
const double SigmaSaSDL::CONVERTDD = 0.0084;

// Factors related to VMD processes in gamma+gamma and gamma+p
// Gammafac = f_V^2 / (4 pi) [cf. Schuler-Sjostrand 1996)
const int SigmaSaSDL::NVMD = 4;
const double SigmaSaSDL::GAMMAFAC[4] = {2.2, 23.6, 18.4, 11.5};
const double SigmaSaSDL::VMDMASS[4]  = {0.77549, 0.78265, 1.01946,
  3.09692};

// Parameters and coefficients for single diffractive scattering.
const int SigmaSaSDL::ISDTABLE[] = { 0, 0, 1, 1, 1, 2, 3, 4, 5,
  6, 7, 8, 9};
const double SigmaSaSDL::CSD[10][8] = {
  { 0.213, 0.0, -0.47, 150., 0.213, 0.0, -0.47, 150., } ,
  { 0.213, 0.0, -0.47, 150., 0.267, 0.0, -0.47, 100., } ,
  { 0.213, 0.0, -0.47, 150., 0.232, 0.0, -0.47, 110., } ,
  { 0.213, 7.0, -0.55, 800., 0.115, 0.0, -0.47, 110., } ,
  { 0.267, 0.0, -0.46,  75., 0.267, 0.0, -0.46,  75., } ,
  { 0.232, 0.0, -0.46,  85., 0.267, 0.0, -0.48, 100., } ,
  { 0.115, 0.0, -0.50,  90., 0.267, 6.0, -0.56, 420., } ,
  { 0.232, 0.0, -0.48, 110., 0.232, 0.0, -0.48, 110., } ,
  { 0.115, 0.0, -0.52, 120., 0.232, 6.0, -0.56, 470., } ,
  { 0.115, 5.5, -0.58, 570., 0.115, 5.5, -0.58, 570.  } };

// Parameters and coefficients for double diffractive scattering.
const int SigmaSaSDL::IDDTABLE[] = { 0, 0, 1, 1, 1, 2, 3, 4, 5,
  6, 7, 8, 9};
const double SigmaSaSDL::CDD[10][9] = {
  { 3.11, -7.34,  9.71, 0.068, -0.42, 1.31, -1.37,  35.0,  118., } ,
  { 3.11, -7.10,  10.6, 0.073, -0.41, 1.17, -1.41,  31.6,   95., } ,
  { 3.12, -7.43,  9.21, 0.067, -0.44, 1.41, -1.35,  36.5,  132., } ,
  { 3.13, -8.18, -4.20, 0.056, -0.71, 3.12, -1.12,  55.2, 1298., } ,
  { 3.11, -6.90,  11.4, 0.078, -0.40, 1.05, -1.40,  28.4,   78., } ,
  { 3.11, -7.13,  10.0, 0.071, -0.41, 1.23, -1.34,  33.1,  105., } ,
  { 3.12, -7.90, -1.49, 0.054, -0.64, 2.72, -1.13,  53.1,  995., } ,
  { 3.11, -7.39,  8.22, 0.065, -0.44, 1.45, -1.36,  38.1,  148., } ,
  { 3.18, -8.95, -3.37, 0.057, -0.76, 3.32, -1.12,  55.6, 1472., } ,
  { 4.18, -29.2,  56.2, 0.074, -1.36, 6.67, -1.14, 116.2, 6532.  } };

//--------------------------------------------------------------------------

// Store pointer to Info and initialize data members.

void SigmaSaSDL::init( Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm*) {

  // Store pointer.
  infoPtr         = infoPtrIn;

  // Initialize parameters for Coulomb corrections to elastic scattering.
  initCoulomb( settings, particleDataPtrIn);

  // User-set values to dampen diffractive cross sections.
  doDampen   = settings.flag("SigmaDiffractive:dampen");
  maxXBOwn   = settings.parm("SigmaDiffractive:maxXB");
  maxAXOwn   = settings.parm("SigmaDiffractive:maxAX");
  maxXXOwn   = settings.parm("SigmaDiffractive:maxXX");
  maxAXBOwn  = settings.parm("SigmaDiffractive:maxAXB");

  // Parameters for diffractive systems.
  epsSaS     = settings.parm("SigmaDiffractive:SaSepsilon");
  sigmaPomP  = settings.parm("Diffraction:sigmaRefPomP");
  mPomP      = settings.parm("Diffraction:mRefPomP");
  pPomP      = settings.parm("Diffraction:mPowPomP");
  zeroAXB    = settings.flag("SigmaTotal:zeroAXB");
  sigAXB2TeV = settings.parm("SigmaTotal:sigmaAXB2TeV");

  // Diffractive mass spectrum starts at m + mMin0 and has a low-mass
  // enhancement, factor cRes, up to around m + mRes0. Special for CD.
  mMin0       = settings.parm("SigmaDiffractive:mMin");
  cRes        = settings.parm("SigmaDiffractive:lowMEnhance");
  mRes0       = settings.parm("SigmaDiffractive:mResMax");
  mMinCDnow   = settings.parm("SigmaDiffractive:mMinCD");

  // Derived quantities.
  alP2       = 2. * ALPHAPRIME;
  s0         = 1. / ALPHAPRIME;

}

//--------------------------------------------------------------------------

// Calculate total and (integrated) elastic cross sections.

bool SigmaSaSDL::calcTotEl( int idAin, int idBin, double sIn, double mAin,
  double mBin) {

  // Find appropriate combination of incoming beams.
  idA     = idAin;
  idB     = idBin;
  s       = sIn;
  isExpEl = true;
  if (!findBeamComb( idA, idB, mAin, mBin)) return false;
  double sEps = pow( s, EPSILON);
  double sEta = pow( s, ETA);

  // Ordinary hadron-hadron collisions. Total cross section.
  if (iProc < 13) {
    sigTot = X[iProc] * sEps + Y[iProc] * sEta;

    // Elastic slope parameter and cross section.
    bEl    = 2. * bA + 2. * bB + 4. * sEps - 4.2;
    sigEl  = CONVERTEL * pow2(sigTot) * (1. + pow2(rhoOwn)) / bEl;

  // gamma + p and p + gamma.
  } else if (iProc == 13){
    sigTot = X[iProc] * sEps + Y[iProc] * sEta;
    double sigElNow  = 0.;

    // Loop over VMD states on side A for elastic cross section.
    for (int iA = 0; iA < NVMD; ++iA){
      double bANow  = BHAD[iHadAtmp[iA]];
      double bBNow  = BHAD[iHadBtmp[iA]];
      double bElNow = 2. * bANow + 2. * bBNow + 4. * sEps - 4.2;
      double tmpTot = X[iProcVP[iA]] * sEps
                    + Y[iProcVP[iA]] * sEta;
      sigElNow     += multVP[iA] * CONVERTEL * pow2(tmpTot)
                   * (1. + pow2(rhoOwn)) / bElNow;
    }
    sigEl  = sigElNow;

  // gamma + gamma
  } else if (iProc == 14){
    sigTot = X[iProc] * sEps + Y[iProc] * sEta;
    double sigElNow  = 0.;

    // Loop over VMD states on side A and B for elastic cross section.
    for (int iA = 0; iA < NVMD; ++iA)
    for (int iB = 0; iB < NVMD; ++iB) {
      double bANow  = BHAD[iHadAtmp[iA]];
      double bBNow  = BHAD[iHadBtmp[iB]];
      double bElNow = 2. * bANow + 2. * bBNow + 4. * sEps - 4.2;
      double tmpTot = X[iProcVV[iA][iB]] * sEps
                    + Y[iProcVV[iA][iB]] * sEta;
      sigElNow     += multVV[iA][iB] * CONVERTEL * pow2(tmpTot)
                   * (1. + pow2(rhoOwn)) / bElNow;
    }
    sigEl  = sigElNow;

  // Primitive implementation of Pomeron + p. No elastic scattering.
  } else if (iProc == 15) {
    sigTot = sigmaPomP * pow( sqrt(s) / mPomP, pPomP);
    sigEl  = 0.;
  }

  // Possibly add Coulomb correction and interference.
  addCoulomb();

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Differential elastic cross section.

double SigmaSaSDL::dsigmaEl( double t, bool useCoulomb, bool ) {

  // Ordinary hadron-hadron collisions.
  double dsig = 0.;
  if (iProc < 13) {
    dsig = CONVERTEL * pow2(sigTot) * (1. + pow2(rhoOwn)) * exp(bEl * t);

  // gamma + p and p + gamma: loop over VMD states on side A.
  } else if (iProc == 13) {
    double sEps     = pow( s, EPSILON);
    double sEta     = pow( s, ETA);
    double dsigNow  = 0.;
    for (int iA = 0; iA < NVMD; ++iA){
      // Elastic slope parameter and cross section.
      double bANow  = BHAD[iHadAtmp[iA]];
      double bBNow  = BHAD[iHadBtmp[iA]];
      double bElNow = 2. * bANow + 2. * bBNow + 4. * sEps - 4.2;
      double tmpTot = X[iProcVP[iA]] * sEps
                    + Y[iProcVP[iA]] * sEta;
      // Hadronic contribution: simple exponential.
      dsigNow += multVP[iA] * CONVERTEL * pow2(tmpTot)
              * (1. + pow2(rhoOwn)) * exp(bElNow * t);
    }
    dsig = dsigNow;

  // gamma + gamma: loop over VMD states on side A and B.
  } else if (iProc == 14) {
    double sEps     = pow( s, EPSILON);
    double sEta     = pow( s, ETA);
    double dsigNow  = 0.;
    for (int iA = 0; iA < NVMD; ++iA)
    for (int iB = 0; iB < NVMD; ++iB){
      // Elastic slope parameter and cross section.
      double bANow  = BHAD[iHadAtmp[iA]];
      double bBNow  = BHAD[iHadBtmp[iB]];
      double bElNow = 2. * bANow + 2. * bBNow + 4. * sEps - 4.2;
      double tmpTot = X[iProcVV[iA][iB]] * sEps
                    + Y[iProcVV[iA][iB]] * sEta;
      // Hadronic contribution: simple exponential.
      dsigNow += multVV[iA][iB] * CONVERTEL * pow2(tmpTot)
              * (1. + pow2(rhoOwn)) * exp(bElNow * t);
    }
    dsig = dsigNow;

  // No elastic scattering for Pomeron + p.
  } else dsig = 0.;

  // Possibly add Coulomb contribution and interference.
  if (useCoulomb && hasCou) dsig += dsigmaElCoulomb(t);

  // Done.
  return dsig;

}

//--------------------------------------------------------------------------

// Diffractive cross sections.

bool SigmaSaSDL::calcDiff(  int idAin, int idBin, double sIn, double mAin,
  double mBin) {

  // Find appropriate combination of incoming beams.
  s   = sIn;
  idA = idAin;
  idB = idBin;
  mA  = mAin;
  mB  = mBin;
  if (!findBeamComb( idAin, idBin, mAin, mBin)) return false;
  sigXB = sigAX = sigXX = sigAXB = 0.;

  // Ordinary hadron-hadron collisions.
  if (iProc < 13) {
    // Lookup coefficients for single and double diffraction.
    int iSD = ISDTABLE[iProc];
    int iDD = IDDTABLE[iProc];
    double eCM = sqrt(s);
    double sum1, sum2, sum3, sum4;

    // Single diffractive scattering A + B -> X + B cross section.
    mMinXB          = mA + mMin0;
    double sMinXB   = pow2(mMinXB);
    mResXB          = mA + mRes0;
    sResXB          = pow2(mResXB);
    double sRMavgXB = mResXB * mMinXB;
    double sRMlogXB = log(1. + sResXB/sMinXB);
    double sMaxXB   = CSD[iSD][0] * s + CSD[iSD][1];
    double BcorrXB  = CSD[iSD][2] + CSD[iSD][3] / s;
    sum1  = log( (2.*bB + alP2 * log(s/sMinXB))
          / (2.*bB + alP2 * log(s/sMaxXB)) ) / alP2;
    sum2  = cRes * sRMlogXB / (2.*bB + alP2 * log(s/sRMavgXB) + BcorrXB) ;
    sigXB = CONVERTSD * X[iProc] * BETA0[iHadB] * max( 0., sum1 + sum2);

    // Single diffractive scattering A + B -> A + X cross section.
    mMinAX          = mB + mMin0;
    double sMinAX   = pow2(mMinAX);
    mResAX          = mB + mRes0;
    sResAX          = pow2(mResAX);
    double sRMavgAX = mResAX * mMinAX;
    double sRMlogAX = log(1. + sResAX/sMinAX);
    double sMaxAX   = CSD[iSD][4] * s + CSD[iSD][5];
    double BcorrAX  = CSD[iSD][6] + CSD[iSD][7] / s;
    sum1  = log( (2.*bA + alP2 * log(s/sMinAX))
          / (2.*bA + alP2 * log(s/sMaxAX)) ) / alP2;
    sum2  = cRes * sRMlogAX / (2.*bA + alP2 * log(s/sRMavgAX) + BcorrAX) ;
    sigAX = CONVERTSD * X[iProc] * BETA0[iHadA] * max( 0., sum1 + sum2);

    // Order single diffractive correctly.
    if (swapped) {
      swap( bB, bA);
      swap( sigXB, sigAX);
      swap( mMinXB, mMinAX);
      swap( mResXB, mResAX);
      swap( iHadA, iHadB);
    }

    // Double diffractive scattering A + B -> X1 + X2 cross section.
    double y0min = log( s * SPROTON / (sMinXB * sMinAX) ) ;
    double sLog  = log(s);
    double Delta0 = CDD[iDD][0] + CDD[iDD][1] / sLog
      + CDD[iDD][2] / pow2(sLog);
    sum1 = (y0min * (log( max( 1e-10, y0min/Delta0) ) - 1.) + Delta0)/ alP2;
    if (y0min < 0.) sum1 = 0.;
    double sMaxXX = s * ( CDD[iDD][3] + CDD[iDD][4] / sLog
      + CDD[iDD][5] / pow2(sLog) );
    double sLogUp = log( max( 1.1, s * s0 / (sMinXB * sRMavgAX) ));
    double sLogDn = log( max( 1.1, s * s0 / (sMaxXX * sRMavgAX) ));
    sum2   = cRes * log( sLogUp / sLogDn ) * sRMlogAX / alP2;
    sLogUp = log( max( 1.1, s * s0 / (sMinAX * sRMavgXB) ));
    sLogDn = log( max( 1.1, s * s0 / (sMaxXX * sRMavgXB) ));
    sum3   = cRes * log(sLogUp / sLogDn) * sRMlogXB / alP2;
    double BcorrXX =  CDD[iDD][6] + CDD[iDD][7] / eCM + CDD[iDD][8] / s;
    sum4   = pow2(cRes) * sRMlogAX * sRMlogXB
      / max( 0.1, alP2 * log( s * s0 / (sRMavgAX * sRMavgXB) ) + BcorrXX);
    sigXX  = CONVERTDD * X[iProc] * max( 0., sum1 + sum2 + sum3 + sum4);

    // Central diffractive scattering A + B -> A + X + B, only p and pbar.
    mMinAXB = 1.;
    if ( (idAbsA == 2212 || idAbsA == 2112)
      && (idAbsB == 2212 || idAbsB == 2112) && !zeroAXB) {
      double sMinAXB = pow2(mMinAXB);
      double sRefAXB = pow2(2000.);
      sigAXB = sigAXB2TeV * pow( log(0.06 * s / sMinAXB), 1.5 )
             / pow( log(0.06 * sRefAXB / sMinAXB), 1.5 );
    }

    // Option with user-requested damping of diffractive cross sections.
    if (doDampen) {
      sigXB  = sigXB  * maxXBOwn  / (sigXB  + maxXBOwn);
      sigAX  = sigAX  * maxAXOwn  / (sigAX  + maxAXOwn);
      sigXX  = sigXX  * maxXXOwn  / (sigXX  + maxXXOwn);
      sigAXB = (maxAXBOwn > 0.) ? sigAXB * maxAXBOwn / (sigAXB + maxAXBOwn)
             : 0.;
    }

  // gamma + p: loop over VMD states on side A.
  } else if (iProc == 13) {
    double sigAXNow  = 0.;
    double sigXBNow  = 0.;
    double sigXXNow  = 0.;

    for (int iA = 0; iA < NVMD; ++iA){

      // Lookup coefficients for single and double diffraction.
      int iSD = ISDTABLE[iProcVP[iA]];
      int iDD = IDDTABLE[iProcVP[iA]];
      double eCM = sqrt(s);
      double sum1, sum2, sum3, sum4;

      // Single diffractive scattering A + B -> X + B cross section.
      mMinXB          = mAtmp[iA] + mMin0;
      double sMinXB   = pow2(mMinXB);
      mResXB          = mAtmp[iA] + mRes0;
      sResXB          = pow2(mResXB);
      double sRMavgXB = mResXB * mMinXB;
      double sRMlogXB = log(1. + sResXB/sMinXB);
      double sMaxXB   = CSD[iSD][0] * s + CSD[iSD][1];
      double BcorrXB  = CSD[iSD][2] + CSD[iSD][3] / s;
      double bBNow    = BHAD[iHadBtmp[iA]];
      sum1      = log( (2.*bBNow + alP2 * log(s/sMinXB))
                / (2.*bBNow + alP2 * log(s/sMaxXB)) ) / alP2;
      sum2      = cRes * sRMlogXB
                / (2.*bBNow + alP2 * log(s/sRMavgXB) + BcorrXB) ;
      sigXBNow += multVP[iA] * CONVERTSD * X[iProcVP[iA]]
                * BETA0[iHadBtmp[iA]] * max( 0., sum1 + sum2);

      // Single diffractive scattering A + B -> A + X cross section.
      mMinAX          = mBtmp[iA] + mMin0;
      double sMinAX   = pow2(mMinAX);
      mResAX          = mBtmp[iA] + mRes0;
      sResAX          = pow2(mResAX);
      double sRMavgAX = mResAX * mMinAX;
      double sRMlogAX = log(1. + sResAX/sMinAX);
      double sMaxAX   = CSD[iSD][4] * s + CSD[iSD][5];
      double BcorrAX  = CSD[iSD][6] + CSD[iSD][7] / s;
      double bANow    = BHAD[iHadAtmp[iA]];
      sum1      = log( (2.*bANow + alP2 * log(s/sMinAX))
                / (2.*bANow + alP2 * log(s/sMaxAX)) ) / alP2;
      sum2      = cRes * sRMlogAX
                / (2.*bANow + alP2 * log(s/sRMavgAX) + BcorrAX) ;
      sigAXNow += multVP[iA] * CONVERTSD * X[iProcVP[iA]]
                * BETA0[iHadAtmp[iA]] * max( 0., sum1 + sum2);

      // Double diffractive scattering A + B -> X1 + X2 cross section.
      double y0min = log( s * SPROTON / (sMinXB * sMinAX) ) ;
      double sLog  = log(s);
      double Delta0 = CDD[iDD][0] + CDD[iDD][1] / sLog
                    + CDD[iDD][2] / pow2(sLog);
      sum1 = (y0min * (log( max( 1e-10, y0min/Delta0) ) - 1.) + Delta0)/ alP2;
      if (y0min < 0.) sum1 = 0.;
      double sMaxXX = s * ( CDD[iDD][3] + CDD[iDD][4] / sLog
                    + CDD[iDD][5] / pow2(sLog) );
      double sLogUp = log( max( 1.1, s * s0 / (sMinXB * sRMavgAX) ));
      double sLogDn = log( max( 1.1, s * s0 / (sMaxXX * sRMavgAX) ));
      sum2   = cRes * log( sLogUp / sLogDn ) * sRMlogAX / alP2;
      sLogUp = log( max( 1.1, s * s0 / (sMinAX * sRMavgXB) ));
      sLogDn = log( max( 1.1, s * s0 / (sMaxXX * sRMavgXB) ));
      sum3   = cRes * log(sLogUp / sLogDn) * sRMlogXB / alP2;
      double BcorrXX =  CDD[iDD][6] + CDD[iDD][7] / eCM + CDD[iDD][8] / s;
      sum4   = pow2(cRes) * sRMlogAX * sRMlogXB
        /max( 0.1, alP2 * log( s * s0 / (sRMavgAX * sRMavgXB) ) + BcorrXX);
      sigXXNow  += multVP[iA] * CONVERTDD * X[iProcVP[iA]]
                 * max( 0., sum1 + sum2 + sum3 + sum4);

    // End of VMD loop.
    }

    // Order single diffractive correctly.
    if (swapped) {
      swap( bB, bA);
      swap( sigXB, sigAX);
      swap( mMinXB, mMinAX);
      swap( mResXB, mResAX);
      swap( iHadA, iHadB);
      swap( sigXBNow, sigAXNow);
      for (int i = 0; i < NVMD; ++i){
        swap( iHadAtmp[i], iHadBtmp[i]);
        swap( mAtmp[i], mBtmp[i]);
      }
    }

    // Option with user-requested damping of diffractive cross sections.
    if (doDampen) {
      sigXBNow  = sigXBNow  * maxXBOwn  / (sigXBNow  + maxXBOwn);
      sigAXNow  = sigAXNow  * maxAXOwn  / (sigAXNow  + maxAXOwn);
      sigXXNow  = sigXXNow  * maxXXOwn  / (sigXXNow  + maxXXOwn);
    }
    sigAX  = sigAXNow;
    sigXB  = sigXBNow;
    sigXX  = sigXXNow;
    sigAXB = 0.;

  // gamma + gamma: loop over VMD states on sides A and B.
  } else if (iProc == 14) {
    double sigAXNow  = 0.;
    double sigXBNow  = 0.;
    double sigXXNow  = 0.;
    for (int iA = 0; iA < NVMD; ++iA)
    for (int iB = 0; iB < NVMD; ++iB){

      // Lookup coefficients for single and double diffraction.
      int iSD = ISDTABLE[iProcVV[iA][iB]];
      int iDD = IDDTABLE[iProcVV[iA][iB]];
      double eCM = sqrt(s);
      double sum1, sum2, sum3, sum4;

      // Single diffractive scattering A + B -> X + B cross section.
      mMinXB          = mAtmp[iA] + mMin0;
      double sMinXB   = pow2(mMinXB);
      mResXB          = mAtmp[iA] + mRes0;
      sResXB          = pow2(mResXB);
      double sRMavgXB = mResXB * mMinXB;
      double sRMlogXB = log(1. + sResXB/sMinXB);
      double sMaxXB   = CSD[iSD][0] * s + CSD[iSD][1];
      double BcorrXB  = CSD[iSD][2] + CSD[iSD][3] / s;
      double bBNow    = BHAD[iHadBtmp[iB]];
      sum1      = log( (2.*bBNow + alP2 * log(s/sMinXB))
                / (2.*bBNow + alP2 * log(s/sMaxXB)) ) / alP2;
      sum2      = cRes * sRMlogXB
                / (2.*bBNow + alP2 * log(s/sRMavgXB) + BcorrXB) ;
      sigXBNow += multVV[iA][iB] * CONVERTSD * X[iProcVV[iA][iB]]
                * BETA0[iHadBtmp[iB]] * max( 0., sum1 + sum2);

      // Single diffractive scattering A + B -> A + X cross section.
      mMinAX          = mBtmp[iB] + mMin0;
      double sMinAX   = pow2(mMinAX);
      mResAX          = mBtmp[iB] + mRes0;
      sResAX          = pow2(mResAX);
      double sRMavgAX = mResAX * mMinAX;
      double sRMlogAX = log(1. + sResAX/sMinAX);
      double sMaxAX   = CSD[iSD][4] * s + CSD[iSD][5];
      double BcorrAX  = CSD[iSD][6] + CSD[iSD][7] / s;
      double bANow    = BHAD[iHadAtmp[iA]];
      sum1      = log( (2.*bANow + alP2 * log(s/sMinAX))
                / (2.*bANow + alP2 * log(s/sMaxAX)) ) / alP2;
      sum2      = cRes * sRMlogAX
                / (2.*bANow + alP2 * log(s/sRMavgAX) + BcorrAX) ;
      sigAXNow += multVV[iA][iB] * CONVERTSD * X[iProcVV[iA][iB]]
                * BETA0[iHadAtmp[iA]] * max( 0., sum1 + sum2);

      // Double diffractive scattering A + B -> X1 + X2 cross section.
      double y0min = log( s * SPROTON / (sMinXB * sMinAX) ) ;
      double sLog  = log(s);
      double Delta0 = CDD[iDD][0] + CDD[iDD][1] / sLog
                    + CDD[iDD][2] / pow2(sLog);
      sum1 = (y0min * (log( max( 1e-10, y0min/Delta0) ) - 1.) + Delta0)/ alP2;
      if (y0min < 0.) sum1 = 0.;
      double sMaxXX = s * ( CDD[iDD][3] + CDD[iDD][4] / sLog
                    + CDD[iDD][5] / pow2(sLog) );
      double sLogUp = log( max( 1.1, s * s0 / (sMinXB * sRMavgAX) ));
      double sLogDn = log( max( 1.1, s * s0 / (sMaxXX * sRMavgAX) ));
      sum2   = cRes * log( sLogUp / sLogDn ) * sRMlogAX / alP2;
      sLogUp = log( max( 1.1, s * s0 / (sMinAX * sRMavgXB) ));
      sLogDn = log( max( 1.1, s * s0 / (sMaxXX * sRMavgXB) ));
      sum3   = cRes * log(sLogUp / sLogDn) * sRMlogXB / alP2;
      double BcorrXX =  CDD[iDD][6] + CDD[iDD][7] / eCM + CDD[iDD][8] / s;
      sum4   = pow2(cRes) * sRMlogAX * sRMlogXB
        /max( 0.1, alP2 * log( s * s0 / (sRMavgAX * sRMavgXB) ) + BcorrXX);
      sigXXNow  += multVV[iA][iB] * CONVERTDD * X[iProcVV[iA][iB]]
                 * max( 0., sum1 + sum2 + sum3 + sum4);

    // End of VMD loop.
    }

    // Option with user-requested damping of diffractive cross sections.
    if (doDampen) {
      sigXBNow  = sigXBNow  * maxXBOwn  / (sigXBNow  + maxXBOwn);
      sigAXNow  = sigAXNow  * maxAXOwn  / (sigAXNow  + maxAXOwn);
      sigXXNow  = sigXXNow  * maxXXOwn  / (sigXXNow  + maxXXOwn);
    }
    sigAX  = sigAXNow;
    sigXB  = sigXBNow;
    sigXX  = sigXXNow;
    sigAXB = 0.;

  // No diffractive scattering for Pomeron + p.
  } else return false;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Differential single diffractive cross sections.

double SigmaSaSDL::dsigmaSD(double xi, double t, bool isXB, int ) {

  // Calculate mass; later checked to be above threshold. Optional weight.
  double m2X   = xi * s;
  double mX    = sqrt(m2X);
  double epsWt = pow( m2X, -epsSaS);

  // Ordinary hadron-hadron collisions.
  if (iProc < 13) {

    // Separate XB and AX cases since A and B masses/couplings may differ.
    if (isXB) {
      if (mX < mMinXB || pow2(mX + mMinAX) > s) return 0.;

      // Return differential AB -> XB cross section value.
      double bXB = 2. * bB + alP2 * log(1. / xi) ;
      return CONVERTSD * X[iProc] * BETA0[iHadB] * exp(bXB * t) * (1. - xi)
        * ( 1. + cRes * sResXB / (sResXB + m2X) ) * epsWt;

    } else {
      if (mX < mMinAX || pow2(mX + mMinXB) > s) return 0.;

      // Return differential AB -> AX cross section value.
      double bAX = 2. * bA + alP2 * log(1. / xi) ;
      return CONVERTSD * X[iProc] * BETA0[iHadA] * exp(bAX * t) * (1. - xi)
        * ( 1. + cRes * sResAX / (sResAX + m2X) ) * epsWt;
    }

  // gamma + p: loop over VMD states on side A.
  } else if (iProc == 13) {
    double sigNow   = 0.;
    for (int iA = 0; iA < NVMD; ++iA){

      // Mass thresholds depend on VMD state.
      mMinXB = mAtmp[iA] + mMin0;
      mResXB = mAtmp[iA] + mRes0;
      sResXB = pow2(mResXB);
      mMinAX = mBtmp[iA] + mMin0;
      mResAX = mBtmp[iA] + mRes0;
      sResAX = pow2(mResAX);

      // Calculate and sum differential AB -> XB cross section value.
      if (isXB && mX > mMinXB && pow2(mX + mMinAX) < s) {
        double bBNow  = BHAD[iHadBtmp[iA]];
        double bXBNow = 2. * bBNow + alP2 * log(1. / xi) ;
        sigNow       += multVP[iA] * CONVERTSD * X[iProcVP[iA]]
                      * BETA0[iHadBtmp[iA]] * exp(bXBNow * t) * (1. - xi)
                      * ( 1. + cRes * sResXB / (sResXB + m2X) );

      // Calculate and sum differential AB -> AX cross section value.
      } else if (!isXB && mX > mMinAX && pow2(mX + mMinXB) < s) {
        double bANow  = BHAD[iHadAtmp[iA]];
        double bAXNow = 2. * bANow + alP2 * log(1. / xi) ;
        sigNow       += multVP[iA] *  CONVERTSD * X[iProcVP[iA]]
                      * BETA0[iHadAtmp[iA]] * exp(bAXNow * t) * (1. - xi)
                      * ( 1. + cRes * sResAX / (sResAX + m2X) );
      }
    }
    return sigNow * epsWt;

  // gamma + gamma: loop over VMD states on side A and B.
  } else if (iProc == 14) {
    double sigNow   = 0.;
    for (int iA = 0; iA < NVMD; ++iA)
    for (int iB = 0; iB < NVMD; ++iB){

      // Mass thresholds depend on VMD state.
      mMinXB = mAtmp[iA] + mMin0;
      mResXB = mAtmp[iA] + mRes0;
      sResXB = pow2(mResXB);
      mMinAX = mBtmp[iB] + mMin0;
      mResAX = mBtmp[iB] + mRes0;
      sResAX = pow2(mResAX);

      // Calculate and sum differential AB -> XB cross section value.
      if (isXB && mX > mMinXB && pow2(mX + mMinAX) < s) {
        double bBNow  = BHAD[iHadBtmp[iB]];
        double bXBNow = 2. * bBNow + alP2 * log(1. / xi) ;
        sigNow       += multVV[iA][iB] * CONVERTSD * X[iProcVV[iA][iB]]
                      * BETA0[iHadBtmp[iB]] * exp(bXBNow * t) * (1. - xi)
                      * ( 1. + cRes * sResXB / (sResXB + m2X) );

      // Calculate and sum differential AB -> AX cross section value.
      } else if (!isXB && mX > mMinAX && pow2(mX + mMinXB) < s) {
        double bANow  = BHAD[iHadAtmp[iA]];
        double bAXNow = 2. * bANow + alP2 * log(1. / xi) ;
        sigNow       += multVV[iA][iB] *  CONVERTSD * X[iProcVV[iA][iB]]
                      * BETA0[iHadAtmp[iA]] * exp(bAXNow * t) * (1. - xi)
                      * ( 1. + cRes * sResAX / (sResAX + m2X) );
      }
    }
    return sigNow * epsWt;

  // No diffractive scattering for Pomeron + p.
  } else return 0.;

}

//--------------------------------------------------------------------------

// Differential double diffractive cross section.

double SigmaSaSDL::dsigmaDD(double xi1, double xi2, double t, int ) {

  // Calculate masses and check that they are above threshold. Optional weight.
  double m2X1  = xi1 * s;
  double mX1   = sqrt(m2X1);
  double m2X2  = xi2 * s;
  double mX2   = sqrt(m2X2);
  double epsWt = pow( m2X1 * m2X2, -epsSaS);

  // Ordinary hadron-hadron collisions.
  if (iProc < 13) {
    if (mX1 < mMinXB || mX2 < mMinAX) return 0.;

    // Return differential AB -> X1X2 cross section value.
    double bXX = alP2 * log( exp(4.) + s * s0 / (m2X1 * m2X2) );
    return CONVERTDD * BETA0[iHadA] * BETA0[iHadB] * exp(bXX * t)
      * (1. - pow2(mX1 + mX2) / s)
      * (s * SPROTON / (s * SPROTON + m2X1 * m2X2))
      * ( 1. + cRes * sResXB / (sResXB + m2X1) )
      * ( 1. + cRes * sResAX / (sResAX + m2X2) ) * epsWt;

  // gamma + p: loop over VMD states on side A.
  } else if (iProc == 13){
    double sigNow   = 0.;
    for (int iA = 0; iA < NVMD; ++iA){

      // Mass thresholds depend on VMD state.
      mMinXB = mAtmp[iA] + mMin0;
      mResXB = mAtmp[iA] + mRes0;
      sResXB = pow2(mResXB);
      mMinAX = mBtmp[iA] + mMin0;
      mResAX = mBtmp[iA] + mRes0;
      sResAX = pow2(mResAX);

      // Calculate and sum differential AB -> X1X2 cross section value.
      if (mX1 > mMinXB && mX2 > mMinAX) {
        double bXX = alP2 * log( exp(4.) + s * s0 / (m2X1 * m2X2) );
        sigNow    += multVP[iA] * CONVERTDD * BETA0[iHadAtmp[iA]]
                   * BETA0[iHadBtmp[iA]] * exp(bXX * t)
                   * (1. - pow2(mX1 + mX2) / s)
                   * (s * SPROTON / (s * SPROTON + m2X1 * m2X2))
                   * ( 1. + cRes * sResXB / (sResXB + m2X1) )
                   * ( 1. + cRes * sResAX / (sResAX + m2X2) );
      }
    }
    return sigNow * epsWt;

  // gamma + gamma: loop over VMD states on side A and B.
  } else if (iProc == 14){
    double sigNow   = 0.;
    for (int iA = 0; iA < NVMD; ++iA)
    for (int iB = 0; iB < NVMD; ++iB){

      // Mass thresholds depend on VMD state.
      mMinXB = mAtmp[iA] + mMin0;
      mResXB = mAtmp[iA] + mRes0;
      sResXB = pow2(mResXB);
      mMinAX = mBtmp[iB] + mMin0;
      mResAX = mBtmp[iB] + mRes0;
      sResAX = pow2(mResAX);

      // Calculate and sum differential AB -> X1X2 cross section value.
      if (mX1 > mMinXB && mX2 > mMinAX) {
        double bXX = alP2 * log( exp(4.) + s * s0 / (m2X1 * m2X2) );
        sigNow    += multVV[iA][iB] * CONVERTDD * BETA0[iHadAtmp[iA]]
                   * BETA0[iHadBtmp[iB]] * exp(bXX * t)
                   * (1. - pow2(mX1 + mX2) / s)
                   * (s * SPROTON / (s * SPROTON + m2X1 * m2X2))
                   * ( 1. + cRes * sResXB / (sResXB + m2X1) )
                   * ( 1. + cRes * sResAX / (sResAX + m2X2) );
      }
    }
    return sigNow * epsWt;

  // No diffractive scattering for Pomeron + p.
  } else return 0.;

}

//--------------------------------------------------------------------------

// Differential central diffractive cross sections.
// Simple extension of single diffraction, but without normalization.

double SigmaSaSDL::dsigmaCD( double xi1, double xi2, double t1, double t2,
  int ) {

  // Ordinary hadron-hadron collisions.
  if (iProc < 13) {

    // Calculate mass; checked to be above threshold.
    double m2X = xi1 * xi2 * s;
    double mX  = sqrt(m2X);
    if (mX < mMinCDnow || pow2(mX + mA + mB) > s) return 0.;
    wtNow      = 1.;

    // Weight associated with differential AB -> AX cross section.
    double bAX = 2. * bA + alP2 * log(1. / xi1) ;
    wtNow     *= CONVERTSD * X[iProc] * BETA0[iHadA] * exp(bAX * t1)
               * (1. - xi1);

    // Weight associated with differential AB -> XB cross section.
    double bXB = 2. * bB + alP2 * log(1. / xi2) ;
    wtNow     *= CONVERTSD * X[iProc] * BETA0[iHadB] * exp(bXB * t2)
               * (1. - xi2);

    // Optional weight to tilt the mass spectrum.
    wtNow     *= pow( m2X, -epsSaS);

    // Done.
    return wtNow;

  // No central diffraction for gamma + p, gamma + gamma or Pomeron + p.
  } else return 0.;

}

//--------------------------------------------------------------------------

// Find beam combination.

bool SigmaSaSDL::findBeamComb( int idAin, int idBin, double mAin,
  double mBin) {

  // Order flavour of incoming hadrons: idAbsA < idAbsB (restore later).
  idAbsA = abs(idAin);
  idAbsB = abs(idBin);
  mA = mAin;
  mB = mBin;
  swapped = false;
  if (idAbsA > idAbsB) {
    swap( idAbsA, idAbsB);
    swap( mA, mB);
    swapped = true;
  }
  sameSign = (idAin * idBin > 0);

  // Find process number.
  iProc                                           = -1;
  if (idAbsA > 1000) {
    iProc                                         = (sameSign) ? 0 : 1;
  } else if (idAbsA > 100 && idAbsB > 1000) {
    iProc                                         = (sameSign) ? 2 : 3;
    if (idAbsA/10 == 11 || idAbsA/10 == 22) iProc = 4;
    if (idAbsA > 300) iProc                       = 5;
    if (idAbsA > 400) iProc                       = 6;
    if (idAbsA > 900) iProc                       = 15;
  } else if (idAbsA > 100) {
    iProc                                         = 7;
    if (idAbsB > 300) iProc                       = 8;
    if (idAbsB > 400) iProc                       = 9;
    if (idAbsA > 300) iProc                       = 10;
    if (idAbsA > 300 && idAbsB > 400) iProc       = 11;
    if (idAbsA > 400) iProc                       = 12;
  } else if (idAbsA == 22 || idAbsB == 22) {
    if (idAbsA == idAbsB) iProc                   = 14;
    if (idAbsA > 1000 || idAbsB > 1000) iProc     = 13;
  }
  if (iProc == -1) return false;

  // Set up global variables.
  iHadA = IHADATABLE[iProc];
  iHadB = IHADBTABLE[iProc];
  bA    = BHAD[iHadA];
  bB    = BHAD[iHadB];

  // Set up VMD global variables for gamma + p.
  if (iProc == 13){
    for (int i = 0; i < NVMD; ++i){
      // VMD always on side a
      mAtmp[i]    = VMDMASS[i];
      mBtmp[i]    = mB;
      iHadAtmp[i] = (i < 2) ? 1 : i;
      iHadBtmp[i] = 0;
      multVP[i]   = ALPHAEM / GAMMAFAC[i];
      if (i < 2) iProcVP[i]  = 4;
      else if (i == 2) iProcVP[i]  = 5;
      else if (i == 3) iProcVP[i]  = 6;
    }
  }

  // Set up VMD global variables for gamma + gamma.
  if (iProc == 14){
    for (int i = 0; i < NVMD; ++i){
      mAtmp[i]      = VMDMASS[i];
      mBtmp[i]      = VMDMASS[i];
      iHadAtmp[i]   = (i < 2) ? 1 : i;
      iHadBtmp[i]   = (i < 2) ? 1 : i;
      for (int j = 0; j < NVMD; ++j){
        multVV[i][j]  = pow2(ALPHAEM) / (GAMMAFAC[i] * GAMMAFAC[j]);
        if ( i < 2 ) {
          if ( j < 2)  iProcVV[i][j] = 7;
          else if ( j == 2) iProcVV[i][j] = 8;
          else if ( j == 3) iProcVV[i][j] = 9;
        } else if (i == 2) {
          if ( j < 2)  iProcVV[i][j] = 8;
          else if ( j == 2) iProcVV[i][j] = 10;
          else if ( j == 3) iProcVV[i][j] = 11;
        } else if ( i == 3) {
          if ( j < 2)  iProcVV[i][j] = 9;
          else if ( j == 2) iProcVV[i][j] = 11;
          else if ( j == 3) iProcVV[i][j] = 12;
        }
      }
    }
  }

  // Done.
  return true ;

}

//==========================================================================

// The SigmaMBR class.
// It parametrizes pp/ppbar total, elastic and diffractive cross sections
// according to the Minimum Bias Rockefeller (MBR) model.

// Integration of MBR cross sections.
const int    SigmaMBR::NINTEG  = 1000;
const int    SigmaMBR::NINTEG2 = 40;

// Proton form factor appoximation with two exponents, [FFB1,FFB2] = GeV^-2.
// Used for quick t-integration, while full expression used for t selection.
const double SigmaMBR::FFA1    = 0.9;
const double SigmaMBR::FFA2    = 0.1;
const double SigmaMBR::FFB1    = 4.6;
const double SigmaMBR::FFB2    = 0.6;

//--------------------------------------------------------------------------

// Initialize data members.

void SigmaMBR::init( Info*, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm*) {

  // Parameters for MBR model.
  eps         = settings.parm("SigmaDiffractive:MBRepsilon");
  alph        = settings.parm("SigmaDiffractive:MBRalpha");
  beta0gev    = settings.parm("SigmaDiffractive:MBRbeta0");
  beta0mb     = beta0gev * sqrt(HBARC2);
  sigma0mb    = settings.parm("SigmaDiffractive:MBRsigma0");
  sigma0gev   = sigma0mb/HBARC2;
  m2min       = settings.parm("SigmaDiffractive:MBRm2Min");
  dyminSDflux = settings.parm("SigmaDiffractive:MBRdyminSDflux");
  dyminDDflux = settings.parm("SigmaDiffractive:MBRdyminDDflux");
  dyminCDflux = settings.parm("SigmaDiffractive:MBRdyminCDflux");
  dyminSD     = settings.parm("SigmaDiffractive:MBRdyminSD");
  dyminDD     = settings.parm("SigmaDiffractive:MBRdyminDD");
  dyminCD     = settings.parm("SigmaDiffractive:MBRdyminCD") / 2.;
  dyminSigSD  = settings.parm("SigmaDiffractive:MBRdyminSigSD");
  dyminSigDD  = settings.parm("SigmaDiffractive:MBRdyminSigDD");
  dyminSigCD  = settings.parm("SigmaDiffractive:MBRdyminSigCD") / sqrt(2.);
  a1          = FFA1;
  a2          = FFA2;
  b1          = FFB1;
  b2          = FFB2;

  // Initialize parameters for Coulomb corrections to elastic scattering.
  initCoulomb( settings, particleDataPtrIn);

  // No rho parameter implemented.
  rhoOwn      = 0.;

}

//--------------------------------------------------------------------------

// Total and elastic cross section.

bool SigmaMBR::calcTotEl( int idAin, int idBin, double sIn, double , double) {

  // Save some input.
  idA       = idAin;
  idB       = idBin;
  s         = sIn;
  isExpEl   = true;

  // Total cross section and elastic/total parametrizations.
  double sCDF  = pow2(1800.);
  double ratio = 0.;
  if (s <= sCDF) {
    double sign = (idA * idB > 0) ? 1. : -1.;
    sigTot = 16.79 * pow(s, 0.104) + 60.81 * pow(s, -0.32)
           - sign * 31.68 * pow(s, -0.54);
    ratio  = 0.100 * pow(s, 0.06) + 0.421 * pow(s, -0.52)
           + sign * 0.160 * pow(s, -0.6);
  } else {
    double sigCDF = 80.03;
    double sF     = pow2(22.);
    sigTot = sigCDF + ( pow2( log(s / sF)) - pow2( log(sCDF / sF)) )
           * M_PI / (3.7 / HBARC2);
    ratio  = 0.066 + 0.0119 * log(s);
  }
  sigEl = sigTot * ratio;
  bEl   = CONVERTEL * pow2(sigTot) / sigEl;

  // Possibly add Coulomb correction and interference.
  addCoulomb();

   // Done.
  return true ;

}

//--------------------------------------------------------------------------

// Differential elastic cross section.

double SigmaMBR::dsigmaEl( double t, bool useCoulomb, bool ) {

  // Hadronic contribution: simple exponential.
  double dsig = sigEl * bEl * exp(bEl * t);

  // Possibly add Coulomb contribution and interference.
  if (useCoulomb && hasCou) dsig += dsigmaElCoulomb(t);

  // Done.
  return dsig;

}

//--------------------------------------------------------------------------

// Total diffractive cross sections, obtained from the ratio of
// two integrals: the Regge cross section and the renormalized flux.

bool SigmaMBR::calcDiff(  int , int , double sIn, double , double ) {

  // Common variables,
  s = sIn;
  double cflux, csig, c1, step, f;
  double dymin0 = 0.;

  // Calculate SD cross section.
  double dymaxSD = log(s / m2min);
  cflux          = pow2(beta0gev) / (16. * M_PI);
  csig           = cflux * sigma0mb;

  // SD flux.
  c1             = cflux;
  double fluxsd  = 0.;
  step           = (dymaxSD - dyminSDflux) / NINTEG;
  for (int i = 0; i < NINTEG; ++i) {
    double dy    = dyminSDflux + (i + 0.5) * step;
    f            = exp(2. * eps * dy) * ( (a1 / (b1 + 2. * alph * dy))
                 + (a2 / (b2 + 2. * alph * dy)) );
    f           *= 0.5 * (1. + erf( (dy - dyminSD) / dyminSigSD));
    fluxsd       = fluxsd + step * c1 * f;
  }
  if (fluxsd < 1.) fluxsd = 1.;

  // Regge SD cross section.
  c1             = csig * pow(s, eps);
  sigSD          = 0.;
  sdpmax         = 0.;
  step           = (dymaxSD - dymin0) / NINTEG;
  for (int i = 0; i < NINTEG; ++i) {
    double dy    = dymin0 + (i + 0.5) * step;
    f            = exp(eps * dy) * ( (a1 / (b1 + 2. * alph * dy))
                 + (a2 / (b2 + 2. * alph * dy)) );
    f           *= 0.5 * (1. + erf( (dy - dyminSD) / dyminSigSD));
    if (f > sdpmax) sdpmax = f;
    sigSD        = sigSD + step * c1 * f;
  }
  sdpmax        *= 1.01;
  sigSD         /= fluxsd;

  // Calculate DD cross section.
  // Note: dymaxDD = ln(s * s0 /mMin^4) with s0 = 1 GeV^2.
  double dymaxDD = log(s / pow2(m2min));
  cflux          = sigma0gev / (16. * M_PI);
  csig           = cflux * sigma0mb;

  // DD flux.
  c1             = cflux / (2. * alph);
  double fluxdd  = 0.;
  step           = (dymaxDD - dyminDDflux) / NINTEG;
  for (int i = 0; i < NINTEG; ++i) {
    double dy    = dyminDDflux + (i + 0.5) * step;
    f            = (dymaxDD - dy) * exp(2. * eps * dy)
                 * ( exp(-2. * alph * dy * exp(-dy))
                 - exp(-2. * alph * dy * exp(dy)) ) / dy;
    f           *= 0.5 * (1. + erf( (dy - dyminDD) / dyminSigDD));
    fluxdd       = fluxdd + step * c1 * f;
  }
  if (fluxdd < 1.) fluxdd = 1.;

  // Regge DD cross section.
  c1             = csig * pow(s, eps) / (2. * alph);
  ddpmax         = 0.;
  sigDD          = 0.;
  step           = (dymaxDD - dymin0) / NINTEG;
  for (int i = 0; i < NINTEG; ++i) {
    double dy    = dymin0 + (i + 0.5) * step;
    f            = (dymaxDD - dy) * exp(eps * dy)
                 * ( exp(-2. * alph * dy * exp(-dy))
                 - exp(-2. * alph * dy * exp(dy)) ) / dy;
    f           *= 0.5 * (1. + erf( (dy - dyminDD) / dyminSigDD));
    if (f > ddpmax) ddpmax = f;
    sigDD        = sigDD + step * c1 * f;
  }
  ddpmax        *= 1.01;
  sigDD         /= fluxdd;

  // Calculate DPE (CD) cross section.
  double dymaxCD = log(s / m2min);
  cflux          = pow4(beta0gev) / pow2(16. * M_PI);
  csig           = cflux * pow2(sigma0mb / beta0mb);
  double dy1, dy2, f1, f2, step2;

  // DPE flux.
  c1             = cflux;
  double fluxdpe = 0.;
  step           = (dymaxCD - dyminCDflux) / NINTEG;
  for (int i = 0; i < NINTEG; ++i) {
    double dy    = dyminCDflux + (i + 0.5) * step;
    f            = 0.;
    step2        = (dy - dyminCDflux) / NINTEG2;
    for (int j = 0; j < NINTEG2; ++j) {
      double yc  = -0.5 * (dy - dyminCDflux) + (j + 0.5) * step2;
      dy1        = 0.5 * dy - yc;
      dy2        = 0.5 * dy + yc;
      f1         = exp(2. * eps * dy1) * ( (a1 / (b1 + 2. * alph * dy1))
                 + (a2 / (b2 + 2. * alph * dy1)) );
      f2         = exp(2. * eps * dy2) * ( (a1 / (b1 + 2. * alph * dy2))
                 + (a2 / (b2 + 2. * alph * dy2)) );
      f1        *= 0.5 * (1. + erf( (dy1 - dyminCD) / dyminSigCD) );
      f2        *= 0.5 * (1. + erf( (dy2 - dyminCD) / dyminSigCD) );
      f         += f1 * f2 * step2;
    }
    fluxdpe     += step * c1 * f;
  }
  if (fluxdpe < 1.) fluxdpe = 1.;

  // Regge DPE (CD) cross section.
  c1             = csig * pow(s, eps);
  sigCD          = 0.;
  dpepmax        = 0;
  step           = (dymaxCD - dymin0) / NINTEG;
  for (int i = 0; i < NINTEG; ++i) {
    double dy    = dymin0 + (i + 0.5) * step;
    f            = 0.;
    step2        = (dy - dymin0) / NINTEG2;
    for (int j = 0; j < NINTEG2; ++j) {
      double yc  = -0.5 * (dy - dymin0) + (j + 0.5) * step2;
      dy1        = 0.5 * dy - yc;
      dy2        = 0.5 * dy + yc;
      f1         = exp(eps * dy1) * ( (a1 / (b1 + 2. * alph * dy1))
                 + (a2 / (b2 + 2. * alph * dy1)) );
      f2         = exp(eps * dy2) * ( (a1 / (b1 + 2. * alph * dy2))
                 + (a2 / (b2 + 2. * alph * dy2)) );
      f1        *= 0.5 * (1. + erf( (dy1 - dyminCD) / dyminSigCD) );
      f2        *= 0.5 * (1. + erf( (dy2 - dyminCD) / dyminSigCD) );
      f         += f1 * f2 * step2;
    }
    sigCD       += step * c1 * f;
    if ( f > dpepmax) dpepmax = f;
  }
  dpepmax       *= 1.01;
  sigCD        /= fluxdpe;

  // Output to standard names. Done.
  sigXB         = sigSD;
  sigAX         = sigSD;
  sigXX         = sigDD;
  sigAXB        = sigCD;
  return true;

}

//--------------------------------------------------------------------------

// Differential single diffractive cross sections, separated in xi and t steps.

double SigmaMBR::dsigmaSD(double xi, double t, bool , int step) {

  // Rapidity gap size.
  double dy = -log(xi);

  // Step 1: evaluate cross section in xi, integrated over t.
  if (step == 1) {
    if (xi * s < m2min) return 0.;
    return exp(eps * dy)
      * ( (a1 / (b1 + 2. * alph * dy)) + (a2 / (b2 + 2. * alph * dy)) )
      * 0.5 * (1. + erf( (dy - dyminSD) / dyminSigSD));

  // Step 2: evaluate cross section in t for fixed xi.
  } else if (step == 2) {
    return pow2(pFormFac(t)) * exp(2. * alph * dy * t);
  }

  // Done.
  return 0.;

}

//--------------------------------------------------------------------------

// Differential double diffractive cross sections, separated in xi and t steps.

double SigmaMBR::dsigmaDD(double xi1, double xi2, double t, int step) {

  // Rapidity gap size (with implicit scale s_0 = 1 GeV).
  double dy   = - log(xi1 * xi2 * s);

  // Step 1: evaluate cross section in xi1 and xi2, integrated over t.
  if (step == 1) {
    if (xi1 * s < m2min || xi2 * s < m2min || dy < 0.) return 0.;
    return exp(eps * dy) * ( exp(-2. * alph * dy * exp(-dy))
           - exp(-2. * alph * dy * exp(dy)) )  / dy
           * 0.5 * (1. + erf( (dy - dyminDD) / dyminSigDD));

  // Step 2: evaluate cross section in t for fixed xi1 and xi2.
  } else if (step == 2) {
    if ( t < -exp(dy) || t > -exp(-dy) ) return 0.;
    return exp(2. * alph * dy * t);
  }

  // Done.
  return 0.;

}

//--------------------------------------------------------------------------

// Differential central diffractive cross section, separated in xi and t steps.

double SigmaMBR::dsigmaCD(double xi1, double xi2, double t1, double t2,
  int step) {

  // Rapidity gap sizes.
  double dy1   = -log(xi1);
  double dy2   = -log(xi2);

  // Step 1: evaluate cross section in xi1 and xi2, integrated over t1 and t2.
  if (step == 1) {
    if (xi1 * xi2 * s < m2min) return 0.;
    double wt1 = exp(eps * dy1)
      * ( (a1 / (b1 + 2. * alph * dy1)) + (a2 / (b2 + 2. * alph * dy1)) )
      * 0.5 * (1. + erf( (dy1 - dyminCD) / dyminSigCD));
    double wt2 = exp(eps * dy2)
      * ( (a1 / (b1 + 2. * alph * dy2)) + (a2 / (b2 + 2. * alph * dy2)) )
      * 0.5 * (1. + erf( (dy2 - dyminCD) / dyminSigCD));
    return wt1 * wt2;

  // Step 2: evaluate cross section in t1 and t2 for fixed xi1 and xi2.
  } else if (step == 2) {
    return pow2(pFormFac(t1) * pFormFac(t2))
      * exp(2. * alph * (dy1 * t1 + dy2 * t2));
  }

  // Done.
  return 0.;

}

//==========================================================================

// The SigmaABMST class.
// It parametrizes pp/ppbar total, elastic and diffractive cross sections
// according to Appleby, Barlow, Molson, Serluca and Toader (ABMST).

//--------------------------------------------------------------------------

// Definitions of static variables.

// Parameters of parametrization: total and elastic cross sections.
const double SigmaABMST::EPSI[]  = {0.106231, 0.0972043, -0.510662, -0.302082};
const double SigmaABMST::ALPP[]  = { 0.0449029, 0.278037, 0.821595, 0.904556};
const double SigmaABMST::NORM[]  = { 228.359, 193.811, 518.686, 10.7843};
const double SigmaABMST::SLOPE[] = { 8.38, 3.78, 1.36};
const double SigmaABMST::FRACS[] = { 0.26, 0.56, 0.18};
const double SigmaABMST::TRIG[]  = { 0.3, 5.03};
const double SigmaABMST::LAM2P   = 0.521223;
const double SigmaABMST::BAPPR[] = { 8.5, 1.086};
const double SigmaABMST::LAM2FF  = 0.71;

// Parameters of parametrization: diffractive cross section.
const double SigmaABMST::MRES[4] = {1.44, 1.52, 1.68, 2.19};
const double SigmaABMST::WRES[4] = {0.325, 0.13, 0.14, 0.45};
const double SigmaABMST::CRES[4] = {3.07, 0.4149, 1.108, 0.9515};
const double SigmaABMST::AFAC[4] = {0.624529, 3.09088, 4.0, 177.217};
const double SigmaABMST::BFAC[4] = {2.5835, 4.51487, 3.03392, 5.86474};
const double SigmaABMST::CFAC[4] = {0., 0.186211, 10.0, 21.0029};
const double SigmaABMST::PPP[4]  = {0.4, 0.5, 0.4597, 5.7575};
const double SigmaABMST::EPS[2]  = {0.08, -0.4525};
const double SigmaABMST::APR[2]  = {0.25, 0.93};
const double SigmaABMST::CPI[6]  = {13.63, 0.0808, 31.79, -0.4525, 0.93, 14.4};
const double SigmaABMST::CNST[5] = {1., -0.05, -0.25, -1.15, 13.5};

// Parameters for integration over t and xi for SD, DD and CD.
const int    SigmaABMST::NPOINTSTSD = 200;
const double SigmaABMST::XIDIVSD    = 0.1;
const double SigmaABMST::DXIRAWSD   = 0.01;
const double SigmaABMST::DLNXIRAWSD = 0.1;
// For DD either Monte Carlo integration or (slower) grid one.
const bool   SigmaABMST::MCINTDD    = true;
const int    SigmaABMST::NPOINTSTDD = 20;
const double SigmaABMST::XIDIVDD    = 0.1;
const double SigmaABMST::DXIRAWDD   = 0.02;
const double SigmaABMST::DLNXIRAWDD = 0.1;
const double SigmaABMST::BMCINTDD   = 2.;
const int    SigmaABMST::NPOINTMCDD = 200000;
// For CD only Monte Carlo integration.
const double SigmaABMST::BMCINTCD   = 2.;
const int    SigmaABMST::NPOINTMCCD = 200000;

//--------------------------------------------------------------------------

// Initialize data members.

void SigmaABMST::init( Info* , Settings& settings, ParticleData* ,
   Rndm* rndmPtrIn) {

  // Save pointer.
  rndmPtr    = rndmPtrIn;

  // Common setup.
  m2minp     = pow2(MPROTON + MPION);
  m2minm     = pow2(MPROTON - MPION);

  // Allow Couplomb corrections for elastic scattering.
  tryCoulomb = settings.flag("SigmaElastic:Coulomb");
  tAbsMin    = settings.parm("SigmaElastic:tAbsMin");

  // Setup for single diffraction.
  modeSD     = settings.mode("SigmaDiffractive:ABMSTmodeSD");
  multSD     = settings.parm("SigmaDiffractive:ABMSTmultSD");
  powSD      = settings.parm("SigmaDiffractive:ABMSTpowSD");
  s0         = (modeSD % 2 == 0) ? 4000. : 100.;
  c0         = (modeSD % 2 == 0) ? 0.6 : 0.012;

  // Setup for double diffraction.
  modeDD     = settings.mode("SigmaDiffractive:ABMSTmodeDD");
  multDD     = settings.parm("SigmaDiffractive:ABMSTmultDD");
  powDD      = settings.parm("SigmaDiffractive:ABMSTpowDD");

  // Setup for central diffraction.
  modeCD     = settings.mode("SigmaDiffractive:ABMSTmodeCD");
  multCD     = settings.parm("SigmaDiffractive:ABMSTmultCD");
  powCD      = settings.parm("SigmaDiffractive:ABMSTpowCD");
  mMinCDnow  = settings.parm("SigmaDiffractive:ABMSTmMinCD");

  // Setup to dampen diffractive events with small rapidity gaps.
  dampenGap  = settings.flag("SigmaDiffractive:ABMSTdampenGap");
  ygap       = settings.parm("SigmaDiffractive:ABMSTygap");
  ypow       = settings.parm("SigmaDiffractive:ABMSTypow");
  expPygap   = exp(ypow * ygap);

  // Setup to force minimal t fall-off like exp(b_min * t).
  useBMin    = settings.flag("SigmaDiffractive:ABMSTuseBMin");
  bMinSD     = settings.parm("SigmaDiffractive:ABMSTbMinSD");
  bMinDD     = settings.parm("SigmaDiffractive:ABMSTbMinDD");
  bMinCD     = settings.parm("SigmaDiffractive:ABMSTbMinCD");

}

//--------------------------------------------------------------------------

// Calculate total and (integrated) elastic cross sections.

bool SigmaABMST::calcTotEl( int idAin, int idBin, double sIn, double ,
  double ) {

  // Find appropriate combination of incoming beams.
  idA     = idAin;
  idB     = idBin;
  ispp    = (idA * idB > 0);
  s       = sIn;
  facEl   = HBARC2 / (16. * M_PI);
  isExpEl = false;

  // Total cross section and the rho parameter using the full expression.
  complex amp = amplitude( 0., false, false);
  sigTot  = HBARC2 * imag(amp);
  rhoOwn  = real(amp) / imag(amp);

  // Total elastic cross section, by integration in exp( MINSLOPEEL * t).
  sigEl   = 0.;
  for (int i = 0; i < NPOINTS; ++i) {
    double y = (i + 0.5) / NPOINTS;
    double t = log(y) / MINSLOPEEL;
    sigEl   += dsigmaEl( t, false) / y;
  }
  sigEl  /= NPOINTS * MINSLOPEEL;

  // Approximate exponential slope.
  bEl = log( dsigmaEl( -TABSREF, false) / dsigmaEl( 0., false) ) / (-TABSREF);

  // Done if no Coulomb corrections.
  hasCou    = tryCoulomb;
  if (abs(idAin) == 2112 || abs(idBin) == 2112) hasCou = false;
  sigTotCou = sigTot;
  sigElCou  = sigEl;
  if (!hasCou) return true;

  // Reduce hadronic part of elastic cross section by tMin cut.
  sigElCou  = sigEl * exp( - bEl * tAbsMin);
  if (tAbsMin < 0.9 * TABSMAX) {

    // Loop through t range according to dt/t^2.
    double sumCou = 0.;
    double xRel, tAbs;
    for (int i = 0; i < NPOINTS; ++i) {
      xRel    = (i + 0.5) / NPOINTS;
      tAbs    = tAbsMin * TABSMAX / (tAbsMin + xRel * (TABSMAX - tAbsMin));

      // Evaluate cross section difference between with and without Coulomb.
      sumCou += pow2(tAbs) * (dsigmaEl( -tAbs, true)
                           -  dsigmaEl( -tAbs, false));
    }

    // Include common factors to give new elastic and total cross sections.
    sigElCou += sumCou * (TABSMAX - tAbsMin)/ (tAbsMin * TABSMAX * NPOINTS);
  }
  sigTotCou   = sigTot - sigEl + sigElCou;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// The scattering amplitude, from which cross sections are derived.

complex SigmaABMST::amplitude( double t, bool useCoulomb,
  bool onlyPomerons) {

  // Common values.
  double snu  = s - 2. * SPROTON + 0.5 * t;
  double ampt = FRACS[0] * exp(SLOPE[0] * t) + FRACS[1] * exp(SLOPE[1] * t)
              + FRACS[2] * exp(SLOPE[2] * t);
  complex amp[6], l2p[4], ll2p[4], d2p[4][3];

  // Two Pomeron and even and odd Reggeon exchange.
  for (int i = 0; i < 4; ++i)
    amp[i] = ((i < 3) ? complex(-NORM[i], 0.) : complex( 0., NORM[i]))
           * ampt * sModAlp( ALPP[i] * snu, 1. + EPSI[i] + ALPP[i] * t);

  // Two-pomeron exchange.
  amp[4] = complex(0., 0.);
  for (int i = 0; i < 4; ++i) {
    l2p[i]  = ALPP[i] * complex( log(ALPP[i] * snu), -0.5 * M_PI);
    ll2p[i] = (1. + EPSI[i]) * l2p[i] / ALPP[i];
    for (int k = 0; k < 3; ++k) d2p[i][k] = SLOPE[k] + l2p[i];
  }
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
  for (int k = 0; k < 3; ++k)
  for (int l = 0; l < 3; ++l) {
    complex part = NORM[i] * NORM[j] * exp( ll2p[i] + ll2p[j] )
                 * exp( t * d2p[i][k] * d2p[j][l] / (d2p[i][k] + d2p[j][l]) )
                 * FRACS[k] * FRACS[l] / (d2p[i][k] + d2p[j][l]);
    if (i == 3) part *= complex( 0., 1.);
    if (j == 3) part *= complex( 0., 1.);
    amp[4]      += part;
  }
  amp[4]        *= LAM2P * complex( 0., 1.) / (16. * M_PI * snu);

  // Triple-gluon exchange.
  amp[5] = sqrt(16. * M_PI / HBARC2) * TRIG[0] * ((t < -TRIG[1])
         ? 1. / pow4(t) :  exp(4. + 4. * t / TRIG[1]) / pow4(TRIG[1]));

  // Add up contributions.
  complex ampSum = 0.;
  if (onlyPomerons) ampSum = (amp[0] + amp[1]) / snu;
  else ampSum = (amp[0] + amp[1] + amp[2] + ((ispp) ? -amp[3] : amp[3])
              + amp[4]) / snu + ((ispp) ? amp[5] : -amp[5]);

  // Optional Coulomb term. Must not be used for t = 0.
  if (useCoulomb && t < 0.) {
    double bApp    = BAPPR[0] + BAPPR[1] * 0.5 * log(s);
    double phase   = (GAMMAEUL  + log( -0.5 * t * (bApp + 8. / LAM2FF)) +
                   - 4. * t / LAM2FF * log(- 4. * t / LAM2FF)
                   - 2. * t / LAM2FF) * (ispp ? 1. : -1.);
    complex ampCou = exp( complex( 0., -ALPHAEM * phase) ) * 8. * M_PI
                   * ALPHAEM * ampt / t;
    ampSum += (ispp) ? ampCou : -ampCou;
  }

  // Done.
  return ampSum;

}

//--------------------------------------------------------------------------

// Diffractive cross sections.

bool SigmaABMST::calcDiff( int idAin , int idBin, double sIn, double ,
  double ) {

  // Find appropriate combination of incoming beams.
  idA    = idAin;
  idB    = idBin;
  ispp   = (idA * idB > 0);
  s      = sIn;
  facEl  = HBARC2 / (16. * M_PI);

  // Total cross section needed for central diffraction.
  complex amp = amplitude( 0., false, true);
  sigTot = HBARC2 * imag(amp);

  // Single diffractive cross sections by grid integration.
  sigXB  = dsigmaSDintXiT( 0., 1., -100., 0.);
  sigAX  = sigXB;

  // Double diffractive cross section by Monte Carlo or grid integration.
  if (MCINTDD) sigXX = dsigmaDDintMC();
  else         sigXX = dsigmaDDintXi1Xi2T( 0., 1., 0., 1., -100., 0.);

  // Central diffractive cross section by Monte Carlo integration.
  sigAXB = dsigmaCDintMC();

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Variations on differential SD cross sections xi * dsigma / dxi dt.

double SigmaABMST::dsigmaSD(double xi, double t, bool, int) {

  // Calculate core SD cross section.
  double dSigSD = dsigmaSDcore( xi, t);

  // Optionally require falloff at least like exp(bMin * t);
  if (useBMin && bMinSD > 0.) {
    double dSigSDmx = dsigmaSDcore( xi, -SPION) * exp(bMinSD * t);
    if (dSigSD > dSigSDmx) dSigSD = dSigSDmx;
  }

  // Optionally dampen with 1 / (1 + exp( -p * (y - y_gap))).
  if (dampenGap) dSigSD /= 1. + expPygap * pow( xi, ypow);

  // Optionally multiply by s-dependent factor.
  if (modeSD > 1) dSigSD *= multSD * pow( s / SPROTON, powSD);

  // Done.
  return dSigSD;

}

//--------------------------------------------------------------------------

// Core differential single diffractive cross sections xi * dsigma / dxi dt.

double SigmaABMST::dsigmaSDcore(double xi, double t) {

  // Calculate mass and check it is above threshold. ABMST valid for |t| < 4.
  double m2X      = xi * s;
  if (m2X < m2minp) return 0.;
  double tAbs     = abs(t);
  if (modeSD % 2 == 0 && tAbs > 4.) return 0.;

  // Some basic parameters.
  double tmp      = 3. + c0 * pow2( log( s / s0) );
  double scaleFac = (s < s0) ? 1. : 3. / tmp;
  double mCut        = (s < s0) ? 3 : tmp;
  if (modeSD % 2 == 0) {
    scaleFac      = 1.;
    mCut          = (s < s0) ? 3. : 3. + c0 * log(s/s0);
  }
  double m2Cut    = pow2(mCut);
  double xiCut    = m2Cut / s;
  double xiThr    = m2minp / s;
  bool   isHighM  = (m2X > m2Cut);
  double xiNow    = (isHighM) ? xi : xiCut;
  double m2XNow   = xiNow * s;

  // Trajectories.
  for (int i = 0; i < 2; ++i) {
    alp0[i] = 1. + EPS[i];
    alpt[i] = alp0[i] + APR[i] * t;
  }
  alpt[2] = CPI[4] * (t - SPION);

  // Individual amplitudes: PPP remains unmodified.
  double ampPPP = pow(xiNow, alp0[0] - 2. * alpt[0]) * pow(s/CNST[0], EPS[0]);
  if (t > CNST[2]) ampPPP *= (PPP[0] + PPP[1] * t);
  else ampPPP *= (AFAC[0] * exp(BFAC[0] * t) + CFAC[0]) * t / (t + CNST[1]);
  if (t < CNST[3]) ampPPP *= (1. + PPP[2] * (tAbs + CNST[3])
                           + PPP[3] * pow2(tAbs + CNST[3]));
  double ampPPR = pow(xiNow, alp0[1] - 2. * alpt[0]) * pow(s/CNST[0], EPS[1]);
  double ampRRP = pow(xiNow, alp0[0] - 2. * alpt[1]) * pow(s/CNST[0], EPS[0]);
  double ampRRR = pow(xiNow, alp0[1] - 2. * alpt[1]) * pow(s/CNST[0], EPS[1]);

  // ABMST uses exponential slope plus constant term in t.
  if (modeSD % 2 == 0) {
    ampPPR *= (AFAC[1] * exp(BFAC[1] * t) + CFAC[1]);
    ampRRP *= (AFAC[2] * exp(BFAC[2] * t) + CFAC[2]);
    ampRRR *= (AFAC[3] * exp(BFAC[3] * t) + CFAC[3]);

  // Modified t form eliminates constant term and has broader exponential.
  } else {
    double modX[2], modX2[2], expX[2], sumX[2];
    double modY[3], modY2[3], expY[3], sumY[3];
    double den[3], modB[3], apC[3];
    for (int ia = 0; ia < 2; ++ia) {
      modX[ia]    = -2. * APR[ia] * log(xiNow);
      modX2[ia]   = pow2(modX[ia]);
      expX[ia]    = exp(-4. * modX[ia]);
      sumX[ia]    = 4. * modX[ia] + 1.;
    }
    for (int ib = 0; ib < 3; ++ib) {
      int ix = (ib == 0) ? 0 : 1;
      modY[ib]    = BFAC[ib+1] + modX[ix];
      modY2[ib]   = pow2(modY[ib]);
      expY[ib]    = exp(-4. * modY[ib]);
      sumY[ib]    = 4. * modY[ib] + 1.;
      den[ib]     = AFAC[ib+1] * modX2[ix] * (1. - expY[ib] * sumY[ib])
                  + CFAC[ib+1] * modY2[ib] * (1. - expX[ix] * sumX[ix]);
      modB[ib]    = (AFAC[ib+1] * modX2[ix] * modY[ib] * (1. - expY[ib])
                  + CFAC[ib+1] * modY2[ib] * modX[ix] * (1. - expX[ix]))
                  / den[ib] - modX[ix];
      apC[ib]     = pow2( AFAC[ib+1] * modX[ix] * (1. - expY[ib])
                  + CFAC[ib+1] * modY[ib] * (1. - expX[ix]) ) / den[ib];
    }
    ampPPR       *= apC[0] * exp(modB[0] * t);
    ampRRP       *= apC[1] * exp(modB[1] * t);
    ampRRR       *= apC[2] * exp(modB[2] * t);
  }

  // Form factor and pion amplitude, remains unmodified.
  double fFac     = pFormFac(t);
  double cnstPi   = CPI[5]/(4. * M_PI) * tAbs/pow2(t - SPION) * pow2(fFac);
  double sigPi    = CPI[0] * pow(m2XNow, CPI[1])
                  + CPI[2] * pow(m2XNow, CPI[3]);
  double ampPi    = cnstPi * sigPi * pow(xiNow, 1. - 2. * alpt[2]);

  // Total high-mass contribution. Done if at high masses.
  double ampHM    = scaleFac * (ampPPP + ampPPR + ampRRP + ampRRR + ampPi);
  if (isHighM) return xi * ampHM;

  // Begin handling of low-mass region.
  double ampBkg   = 0.;
  double ampRes   = 0.;
  double ampMatch = 0.;
  double ampLM    = 0.;

  // Resonance contribution.
  double qRef     = sqrt( (m2X - m2minp) * (m2X - m2minm) / (4. * m2X) );
  for (int i = 0; i < 4; ++i) {
    double m2Res  = pow2(MRES[i]);
    double qRes   = sqrt( (m2Res - m2minp) * (m2Res - m2minm) / (4. * m2Res) );
    double mGam   = MRES[i] * WRES[i] * pow( qRef / qRes, 2. * i + 3.)
                  * pow( (1. + 5. * qRes) / (1. + 5. * qRef), i + 1.);
    ampRes       += CRES[i] * mGam / (pow2(m2X - m2Res) + pow2(mGam));
    ampMatch     += CRES[i] * mGam / (pow2(m2Cut - m2Res) + pow2(mGam));
  }
  ampRes         *= exp( CNST[4] * (t - CNST[1]) ) / xi;
  ampMatch       *= exp( CNST[4] * (t - CNST[1]) ) / xiNow
                  * (xi - xiThr) / (xiNow - xiThr);

  // Background contribution.
  double dAmpPPP  = ampPPP * (alp0[0] - 2. * alpt[0]) / xiNow;
  double dAmpPPR  = ampPPR * (alp0[1] - 2. * alpt[0]) / xiNow;
  double dAmpRRP  = ampRRP * (alp0[0] - 2. * alpt[1]) / xiNow;
  double dAmpRRR  = ampRRR * (alp0[1] - 2. * alpt[1]) / xiNow;
  double dSigPi   = CPI[0] * CPI[1] * pow(m2XNow, CPI[1] - 1.)
                  + CPI[2] * CPI[3] * pow(m2XNow, CPI[3] - 1.);
  double dAmpPi   = cnstPi * (sigPi * (1. - 2. * alpt[2])
                  * pow(xiNow, -2. * alpt[2])
                  + dSigPi * pow(xiNow, 1. - 2. * alpt[2]) );
  double dAmpHM   = scaleFac * (dAmpPPP + dAmpPPR + dAmpRRP + dAmpRRR
                  + dAmpPi);

  // Original background which vanishes quadratically at threshold.
  if (modeSD % 2 == 0) {
    double coeff1 = (dAmpHM * (xiCut - xiThr) - ampHM) / pow2(xiCut - xiThr);
    double coeff2 = 2. * ampHM / (xiCut - xiThr) - dAmpHM;
    ampBkg        = coeff1 * pow2(xi - xiThr) + coeff2 * (xi - xiThr);

  // Modified form with combination of linear and quadratic description.
  } else {
    double xiAtM3 = 9. / s;
    double coeff1 = dAmpHM;
    double coeff2 = ampHM - coeff1 * (xiCut - xiThr);
    double coeff3 = - coeff2 / pow2(xiAtM3 - xiThr);
    double coeff4 = (2. * coeff1 * (xiAtM3 - xiThr) + 2. * coeff2)
                  / (xiAtM3 - xiThr) - coeff1;
    ampBkg        = (xi < xiAtM3) ? coeff3 * pow2(xi - xiThr)
                  + coeff4 * (xi - xiThr) : coeff1 * (xi - xiThr) + coeff2;
  }

  // Return total low-mass contribution.
  ampLM           = ampRes - ampMatch + ampBkg;
  return xi * ampLM;

}

//--------------------------------------------------------------------------

// Single diffractive cross sections integrated over [t_min, t_max].

double SigmaABMST::dsigmaSDintT(double xi, double tMinIn, double tMaxIn) {

  // Calculate t range. Check if range closed.
  double mu1   = SPROTON / s;
  double mu3   = xi;
  double rootv = (1. - 4. * mu1) * (pow2(1. - mu1 - mu3) - 4. * mu1 * mu3);
  if (rootv <= 0.) return 0.;
  double tMin  = -0.5 * s * (1. - 3. * mu1 - mu3 + sqrt(rootv));
  double tMax  = s * s * mu1 * pow2(mu3 - mu1) / tMin;
  tMin         = max( tMin, tMinIn);
  tMax         = min( tMax, tMaxIn);
  if (tMin >= tMax) return 0.;

  // Prepare integration.
  double slope = -0.5 * log(xi);
  double etMin = exp(slope * tMin);
  double etMax = exp(slope * tMax);

  // Do integration by uniform steps in exp(slope * t).
  double dsig  = 0.;
  double etNow, tNow;
  for (int i = 0; i < NPOINTSTSD; ++i) {
    etNow      = etMin + (i + 0.5) * (etMax - etMin) / NPOINTSTSD;
    tNow       = log(etNow) / slope;
    dsig      += dsigmaSD( xi, tNow, true, 0) / etNow;
  }

  // Normalize and done.
  dsig        *= (etMax - etMin) / (NPOINTSTSD * slope);
  return dsig;

}

//--------------------------------------------------------------------------

// Single diffractive cross section integrated over
// [xi_min, xi_max] * [t_min, t_max].

double SigmaABMST::dsigmaSDintXiT( double xiMinIn, double xiMaxIn,
  double tMinIn, double tMaxIn) {

  // Restrictions on xi range. Check it is not closed.
  double sig   = 0.;
  double xiMin = max(xiMinIn, m2minp / s);
  double xiMax = min(xiMaxIn, 1.);
  if (xiMin >= xiMax) return 0.;
  double xiNow;

  // Integration in xi: check size of affected range and adjust dxi.
  if (xiMax > XIDIVSD) {
    double xiMinRng = max( XIDIVSD, xiMin);
    int    nxi      = 2 + (xiMax - xiMinRng) / DXIRAWSD;
    double dxi      = (xiMax - xiMinRng) / nxi;
    for (int ixi = 0; ixi < nxi; ++ixi) {
      xiNow         =  xiMinRng + (ixi + 0.5) * dxi;
      sig          += dxi * dsigmaSDintT( xiNow, tMinIn, tMaxIn) / xiNow;
    }
  }

  // Integration in ln(xi): check size of affected range and adjust dlnxi.
  if (xiMin < XIDIVSD) {
    double xiMaxRng = min( XIDIVSD, xiMax);
    int    nlnxi    = 2 + log( xiMaxRng / xiMin) / DLNXIRAWSD;
    double dlnxi    = log( xiMaxRng / xiMin) / nlnxi;
    for (int ilnxi = 0; ilnxi < nlnxi; ++ilnxi) {
      xiNow         = xiMin * exp( dlnxi * (ilnxi + 0.5));
      sig          += dlnxi * dsigmaSDintT( xiNow, tMinIn, tMaxIn);
    }
  }

  // Done.
  return sig;

}

//--------------------------------------------------------------------------

// Differential double diffractive cross section.

double SigmaABMST::dsigmaDD(double xi1, double xi2, double t, int ) {

  // Calculate mass and check it is above threshold. ABMST valid for |t| < 4.
  double m2X1     = xi1 * s;
  double m2X2     = xi2 * s;
  if (m2X1 < m2minp || m2X2 < m2minp) return 0.;
  if (modeSD % 2 == 0 && abs(t) > 4.) return 0.;

  // dSigma_DD(x1, x2, t) = dSigma_SD(x1, t) * dSigma_SD(x2, t) / dSigma_el(t).
  // Elastic using only Pomerons.
  double dSigDD   = dsigmaSDcore( xi1, t) * dsigmaSDcore( xi2, t)
                  / dsigmaEl( t, false, true);

  // Minimal fall-off relative to t = 0 value.
  if (useBMin && bMinDD > 0.) {
    double dSigDDmx = dsigmaSDcore( xi1, -SPION) * dsigmaSDcore( xi2, -SPION)
                    * exp(bMinDD * t) / dsigmaEl( 0., false, true);
    if (dSigDD > dSigDDmx) dSigDD = dSigDDmx;
  }

  // Optionally dampen with 1 / (1 + exp( -p * (y - y_gap))).
  if (dampenGap) dSigDD /= 1. + expPygap * pow( xi1 * xi2 * s / SPROTON, ypow);

  // Optionally multiply by s-dependent factor.
  if (modeDD == 1) dSigDD *= multDD * pow( s / SPROTON, powDD);

  // Done.
  return dSigDD;

}

//--------------------------------------------------------------------------

// Double diffractive cross section integrated by Monte Carlo.

double SigmaABMST::dsigmaDDintMC() {

  // Set up parameters of integration.
  double sigSum = 0.;
  double xiMin  = m2minp / s;
  double mu1    = SPROTON / s;
  double xi1, xi2, t;

  // Integrate flat in dln(xi1) * dln(xi2) * exp(b_min t) dt.
  for (int iPoint = 0; iPoint < NPOINTMCDD; ++iPoint) {
    xi1   = pow( xiMin, rndmPtr->flat() );
    xi2   = pow( xiMin, rndmPtr->flat() );
    t     = log( rndmPtr->flat() ) / BMCINTDD;

    // Check that point is inside phase space.
    if (sqrt(xi1) + sqrt(xi2) > 1.) continue;
    if (!tInRange( t/s, 1., mu1, mu1, xi1, xi2)) continue;

    // Calculate and add cross section.
    sigSum += dsigmaDD( xi1, xi2, t) * exp(-BMCINTDD * t);
  }

  // Normalize and done.
  sigSum *= pow2(log(xiMin)) / (BMCINTDD * NPOINTMCDD);
  return sigSum;

}

//--------------------------------------------------------------------------

// Double diffractive cross section integrated over [t_min, t_max].

double SigmaABMST::dsigmaDDintT(double xi1, double xi2, double tMinIn,
  double tMaxIn) {

  // Calculate t range. Check if range closed.
  double mu1   = SPROTON / s;
  pair<double,double> tRng = tRange(1., mu1, mu1, xi1, xi2);
  double tMin  = max( s * tRng.first, tMinIn);
  double tMax  = min( s * tRng.second, tMaxIn);
  if (tMin >= tMax) return 0.;

  // Prepare integration.
  double slope = BMCINTDD;
  double etMin = exp(slope * tMin);
  double etMax = exp(slope * tMax);

  // Do integration by uniform steps in exp(slope * t).
  double dsig  = 0.;
  double etNow, tNow;
  for (int i = 0; i < NPOINTSTDD; ++i) {
    etNow      = etMin + (i + 0.5) * (etMax - etMin) / NPOINTSTDD;
    tNow       = log(etNow) / slope;
    dsig      += dsigmaDD( xi1, xi2, tNow) / etNow;
  }

  // Normalize and done.
  dsig        *= (etMax - etMin) / (NPOINTSTDD * slope);
  return dsig;

}

//--------------------------------------------------------------------------

// Double diffractive cross section integrated over
// [xi2_min, xi2_max] * [t_min, t_max].

double SigmaABMST::dsigmaDDintXi2T( double xi1, double xi2MinIn,
  double xi2MaxIn, double tMinIn, double tMaxIn) {

  // Restrictions on xi2 range. Check it is not closed.
  double dsig       = 0.;
  double xi2Min     = max( xi2MinIn, m2minp / s);
  double xi2Max     = min( xi2MaxIn, 1. + xi1 - 2. * sqrt(xi1));
  if (xi2Min >= xi2Max) return 0.;
  double xi2Now;

  // Integration in xi2: check size of affected range and adjust dxi2.
  if (xi2Max > XIDIVDD) {
    double xi2MinRng = max( XIDIVDD, xi2Min);
    int    nxi2     = 2 + (xi2Max - xi2MinRng) / DXIRAWDD;
    double dxi2     = (xi2Max - xi2MinRng) / nxi2;
    for (int ixi2 = 0; ixi2 < nxi2; ++ixi2) {
      xi2Now        =  xi2MinRng + (ixi2 + 0.5) * dxi2;
      dsig         += dxi2 * dsigmaDDintT( xi1, xi2Now, tMinIn, tMaxIn)
                    / xi2Now;
    }
  }

  // Integration in ln(xi2): check size of affected range and adjust dlnxi2.
  if (xi2Min < XIDIVDD) {
    double xi2MaxRng = min( XIDIVDD, xi2Max);
    int    nlnxi2   = 2 + log( xi2MaxRng / xi2Min) / DLNXIRAWDD;
    double dlnxi2   = log( xi2MaxRng / xi2Min) / nlnxi2;
    for (int ilnxi2 = 0; ilnxi2 < nlnxi2; ++ilnxi2) {
      xi2Now        = xi2Min * exp( dlnxi2 * (ilnxi2 + 0.5));
      dsig         += dlnxi2 * dsigmaDDintT( xi1, xi2Now, tMinIn, tMaxIn);
    }
  }

  // Done.
  return dsig;

}

//--------------------------------------------------------------------------

// Double diffractive cross section integrated over
// [xi1_min, xi1_max] * [xi2_min, xi2_max] * [t_min, t_max].

double SigmaABMST::dsigmaDDintXi1Xi2T( double xi1MinIn, double xi1MaxIn,
  double xi2MinIn, double xi2MaxIn, double tMinIn, double tMaxIn) {

  // Restrictions on xi1 range. Check it is not closed.
  double dsig       = 0.;
  double xi1Min     = max( xi1MinIn, m2minp / s);
  double xi1Max     = min( xi1MaxIn, 1.);
  if (xi1Min >= xi1Max) return 0.;
  double xi1Now;

  // Integration in xi1: check size of affected range and adjust dxi1.
  if (xi1Max > XIDIVDD) {
    double xi1MinRng = max( XIDIVDD, xi1Min);
    int    nxi1     = 2 + (xi1Max - xi1MinRng) / DXIRAWDD;
    double dxi1     = (xi1Max - xi1MinRng) / nxi1;
    for (int ixi1 = 0; ixi1 < nxi1; ++ixi1) {
      xi1Now        =  xi1MinRng + (ixi1 + 0.5) * dxi1;
      dsig         += dxi1 * dsigmaDDintXi2T( xi1Now, xi2MinIn, xi2MaxIn,
                      tMinIn, tMaxIn) / xi1Now;
    }
  }

  // Integration in ln(xi1): check size of affected range and adjust dlnxi1.
  if (xi1Min < XIDIVDD) {
    double xi1MaxRng = min( XIDIVDD, xi1Max);
    int    nlnxi1   = 2 + log( xi1MaxRng / xi1Min) / DLNXIRAWDD;
    double dlnxi1   = log( xi1MaxRng / xi1Min) / nlnxi1;
    for (int ilnxi1 = 0; ilnxi1 < nlnxi1; ++ilnxi1) {
      xi1Now        = xi1Min * exp( dlnxi1 * (ilnxi1 + 0.5));
      dsig         += dlnxi1 * dsigmaDDintXi2T( xi1Now, xi2MinIn, xi2MaxIn,
                      tMinIn, tMaxIn);
    }
  }

  // Done.
  return dsig;

}

//--------------------------------------------------------------------------

// Differential central diffractive cross section.

double SigmaABMST::dsigmaCD(double xi1, double xi2, double t1, double t2,
  int ) {

  // ABMST valid for |t| < 4.
  if (modeSD % 2 == 0 && max( abs(t1), abs(t2)) > 4.) return 0.;

  // dSigma_CD(x1, x2, t1, t2)
  // = dSigma_SD(x1, t1) * dSigma_SD(x2, t2) / sigma_tot.
  double dSigCD   = dsigmaSDcore( xi1, t1) * dsigmaSDcore( xi2, t2) / sigTot;

  // Minimal fall-off relative to t = 0 value.
  if (useBMin && bMinCD > 0.) {
    double dSigCDmx = dsigmaSDcore( xi1, -SPION) * dsigmaSDcore( xi2, -SPION)
                    * exp(bMinCD * (t1 + t2)) / sigTot;
    if (dSigCD > dSigCDmx) dSigCD = dSigCDmx;
  }

  // Optionally dampen with 1 / (1 + exp( -p * (y - y_gap))) for both gaps.
  if (dampenGap) dSigCD /= (1. + expPygap * pow( xi1, ypow))
                         * (1. + expPygap * pow( xi2, ypow));

  // Optionally multiply by s-dependent factor.
  if (modeCD == 1) dSigCD *= multCD * pow( s / SPROTON, powCD);

  // Done.
  return dSigCD;

}

//--------------------------------------------------------------------------

// Central diffractive cross section integrated by Monte Carlo.

double SigmaABMST::dsigmaCDintMC() {

  // Set up parameters of integration.
  double sigSum = 0.;
  double xiMin  = m2minp / s;
  double xi1, xi2, t1, t2;

  // Integrate flat in dln(xi1) * exp(b_min t1) dt1 * (same with xi2, t2).
  for (int iPoint = 0; iPoint < NPOINTMCCD; ++iPoint) {
    xi1   = pow( xiMin, rndmPtr->flat() );
    xi2   = pow( xiMin, rndmPtr->flat() );
    t1    = log( rndmPtr->flat() ) / BMCINTCD;
    t2    = log( rndmPtr->flat() ) / BMCINTCD;

    // Check that point is inside phase space.
    if (xi1 * xi2 < xiMin) continue;
    if (xi1 * xi2 + 2. * xiMin > 1.) continue;
    if (!tInRange( t1, s, SPROTON, SPROTON, SPROTON, SPROTON + xi1 * s))
      continue;
    if (!tInRange( t1, s, SPROTON, SPROTON, SPROTON, SPROTON + xi2 * s))
      continue;

    // Calculate and add cross section.
    sigSum += dsigmaCD( xi1, xi2, t1, t2) * exp(-BMCINTCD * (t1 + t2));
  }

  // Normalize and done.
  sigSum *= pow2(log(xiMin) / BMCINTCD) / NPOINTMCCD;
  return sigSum;

}

//==========================================================================

// The SigmaRPP class.
// It parametrizes pp/ppbar total and elastic cross sections according to
// the fit in Review of Particle Physics 2016.

//--------------------------------------------------------------------------

// Definitions of static variables.

// Parameters of parametrization.
const double SigmaRPP::EPS1[] = { 1., 0.614, 0.444, 1., 1., 1.};
const double SigmaRPP::ALPP[] = { 0.151, 0.8, 0.8, 0.947};
const double SigmaRPP::NORM[] = { 0.2478, 0.0078, 11.22, -0.150, 148.4, -26.6,
  -1.5, -0.0441, 0., 0.686, -3.82, -8.60, 64.1, 99.1, -58.0, 9.5};
const double SigmaRPP::BRPP[] = { 3.592, 0.622, 5.44, 0.205, 5.643, 1.92,
  0.41, 0., 0., 3.013, 2.572, 12.25, 2.611, 11.28, 1.27};
const double SigmaRPP::KRPP[] = { 0.3076, 0.0998, 1.678, 0.190, -26.1};
const double SigmaRPP::LAM2FF = 0.71;

//--------------------------------------------------------------------------

// Calculate total and (integrated) elastic cross sections.

bool SigmaRPP::calcTotEl( int idAin, int idBin, double sIn, double ,
  double ) {

  // Find appropriate combination of incoming beams.
  idA     = idAin;
  idB     = idBin;
  ispp    = (idA * idB > 0);
  s       = sIn;
  facEl   = CONVERTEL / (s * (s - 4. * SPROTON));
  isExpEl = false;

  // Total cross section and the rho parameter.
  complex amp = amplitude( 0., false);
  sigTot  = imag(amp) / sqrt(s * ( s - 4. * SPROTON));
  rhoOwn  = real(amp) / imag(amp);

  // Total elastic cross section, by integration in exp( MINSLOPEEL * t).
  sigEl   = 0.;
  for (int i = 0; i < NPOINTS; ++i) {
    double y = (i + 0.5) / NPOINTS;
    double t = log(y) / MINSLOPEEL;
    sigEl   += dsigmaEl( t, false) / y;
  }
  sigEl  /= NPOINTS * MINSLOPEEL;

  // Approximate exponential slope.
  bEl = log( dsigmaEl( -TABSREF, false) / dsigmaEl( 0., false) ) / (-TABSREF);

  // Done if no Coulomb corrections.
  hasCou    = tryCoulomb;
  if (abs(idAin) == 2112 || abs(idBin) == 2112) hasCou = false;
  sigTotCou = sigTot;
  sigElCou  = sigEl;
  if (!hasCou) return true;

  // Reduce hadronic part of elastic cross section by tMin cut.
  sigElCou  = sigEl * exp( - bEl * tAbsMin);
  if (tAbsMin < 0.9 * TABSMAX) {

    // Loop through t range according to dt/t^2.
    double sumCou = 0.;
    double xRel, tAbs;
    for (int i = 0; i < NPOINTS; ++i) {
      xRel    = (i + 0.5) / NPOINTS;
      tAbs    = tAbsMin * TABSMAX / (tAbsMin + xRel * (TABSMAX - tAbsMin));

      // Evaluate cross section difference between with and without Coulomb.
      sumCou += pow2(tAbs) * (dsigmaEl( -tAbs, true)
                           -  dsigmaEl( -tAbs, false));
    }

    // Include common factors to give new elastic and total cross sections.
    sigElCou += sumCou * (TABSMAX - tAbsMin)/ (tAbsMin * TABSMAX * NPOINTS);
  }
  sigTotCou   = sigTot - sigEl + sigElCou;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Amplitude.

complex SigmaRPP::amplitude( double t, bool useCoulomb) {

  // Modified s-related values.
  double  shat   = s - 2. * SPROTON + 0.5 * t;
  complex stlog  = complex( log(shat), -0.5 * M_PI);
  complex taut   = sqrt(abs(t)) * stlog;

  // Trajectories.
  double aP      = EPS1[0] + ALPP[0] * t;
  double aRpos   = EPS1[1] + ALPP[1] * t;
  double aRneg   = EPS1[2] + ALPP[2] * t;
  double aO      = EPS1[3] + ALPP[3] * t;
  double aOP     = EPS1[4] + ALPP[0] * ALPP[3] * t / (ALPP[0] + ALPP[3]);
  double aPP     = EPS1[5] + 0.5 * ALPP[0] * t;
  double aRPpos  = EPS1[1] + ALPP[0] * ALPP[1] * t / (ALPP[0] + ALPP[1]);
  double aRPneg  = EPS1[2] + ALPP[0] * ALPP[2] * t / (ALPP[0] + ALPP[2]);

  // Even terms.
  complex besArg = KRPP[0] * taut;
  complex besJ0n = besJ0(besArg);
  complex besJ1n = besJ1(besArg);
  complex besRat = (abs(besArg) < 0.01) ? 1. : 2. * besJ1n / besArg;
  complex fPosH  = complex( 0., shat) * (NORM[0] * besRat
                 * exp(BRPP[0] * t) * stlog * stlog
                 + NORM[1] * besJ0n * exp(BRPP[1] * t) * stlog
                 + NORM[2] * (besJ0n - besArg * besJ1n) * exp(BRPP[2] * t));
  complex fPosP  = -NORM[3] * exp(BRPP[3] * t) * sModAlp( shat, aP);
  complex fPosPP = (-NORM[4] / stlog) * exp(BRPP[4] * t) * sModAlp( shat, aPP);
  complex fPosR  = -NORM[5] * exp(BRPP[5] * t) * sModAlp( shat, aRpos);
  complex fPosRP = t * (NORM[6] / stlog) * exp(BRPP[6] * t)
                 * sModAlp( shat, aRPpos);
  complex nPos   = complex( 0., -shat) * NORM[7] * stlog * t
                 * pow( 1. - t / KRPP[2], -5.);
  complex fPos   = fPosH + fPosP + fPosPP + fPosR + fPosRP + nPos;

  // Odd terms.
  complex fNegMO = shat * (NORM[9] * cos(KRPP[1] * taut) * exp(BRPP[9] * t)
                 * stlog + NORM[10] * exp(BRPP[10] * t));
  complex fNegO  = complex( 0., NORM[11]) * exp(BRPP[11] * t)
                 * sModAlp( shat, aO) * (1. + KRPP[4] * t);
  complex fNegOP = (complex( 0., NORM[12]) / stlog) * exp(BRPP[12] * t)
                 * sModAlp( shat, aOP);
  complex fNegR  = complex( 0., -NORM[13]) * exp(BRPP[13] * t)
                 * sModAlp( shat, aRneg);
  complex fNegRP = t * (complex( 0., -NORM[14]) / stlog) * exp(BRPP[14] * t)
                 * sModAlp( shat, aRPneg);
  complex nNeg   = -shat * NORM[15] * stlog * t * pow( 1. - t / KRPP[3], -5.);
  complex fNeg   = fNegMO + fNegO + fNegOP + fNegR + fNegRP + nNeg;

  // Combine nuclear part.
  complex ampSum = ispp ? fPos + fNeg : fPos - fNeg;

  // Optional Coulomb term. Must not be used for t = 0.
  complex ampCou = 0.;
  if (useCoulomb && t < 0.) {
    double bAppr = imag(ampSum) / ( sqrt(s * ( s - 4. * SPROTON))
                 * 4. * M_PI * HBARC2 );
    double phase = (log( -0.5 * t * (bAppr + 8. / LAM2FF)) + GAMMAEUL
                 - 4. * t / LAM2FF * log(- 4. * t / LAM2FF)
                 - 2. * t / LAM2FF) * (ispp ? -1. : 1.);
    ampCou       = exp( complex( 0., ALPHAEM * phase) ) * 8. * M_PI * HBARC2
                 * ALPHAEM * s / t * pow(1 - t / LAM2FF, -4.);
  }

  // Combine and return.
  return ispp ? ampSum + ampCou : ampSum - ampCou;

}

//--------------------------------------------------------------------------

// Complex Bessel functions J0 and J1.

complex SigmaRPP::besJ0( complex x) {
  int mMax    = 5. + 5. * abs(x);
  complex z   = 0.25 * x * x;
  complex term = 1.;
  complex sum  = term;
  for (int m = 1; m < mMax; ++m) {
    term *= - z / double(m * m);
    sum  += term;
  }
  return sum;
}

complex SigmaRPP::besJ1( complex x) {
  int mMax    = 5. + 5. * abs(x);
  complex z   = 0.25 * x * x;
  complex term = 0.5 * x;
  complex sum  = term;
  for (int m = 1; m < mMax; ++m) {
    term *= - z / double(m * (m+1));
    sum  += term;
  }
  return sum;
}

//==========================================================================

} // end namespace Pythia8
