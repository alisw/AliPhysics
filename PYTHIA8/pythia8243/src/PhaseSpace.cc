// PhaseSpace.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// PhaseSpace and PhaseSpace2to2tauyz classes.

#include "Pythia8/PhaseSpace.h"

namespace Pythia8 {

//==========================================================================

// The PhaseSpace class.
// Base class for phase space generators.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of trial maxima around which maximum search is performed.
const int    PhaseSpace::NMAXTRY        = 2;

// Number of three-body trials in phase space optimization.
const int    PhaseSpace::NTRY3BODY      = 20;

// Maximum cross section increase, just in case true maximum not found.
const double PhaseSpace::SAFETYMARGIN   = 1.05;

// Small number to avoid division by zero.
const double PhaseSpace::TINY           = 1e-20;

// Fraction of total weight that is shared evenly between all shapes.
const double PhaseSpace::EVENFRAC       = 0.4;

// Two cross sections with a small relative error are assumed same.
const double PhaseSpace::SAMESIGMA      = 1e-6;

// Do not allow resonance to have mass below 2 m_e.
const double PhaseSpace::MRESMINABS     = 0.001;

// Do not include resonances peaked too far outside allowed mass region.
const double PhaseSpace::WIDTHMARGIN    = 20.;

// Special optimization treatment when two resonances at almost same mass.
const double PhaseSpace::SAMEMASS       = 0.01;

// Minimum phase space left when kinematics constraints are combined.
const double PhaseSpace::MASSMARGIN     = 0.01;

// When using Breit-Wigners in 2 -> 2 raise maximum weight estimate.
const double PhaseSpace::EXTRABWWTMAX   = 1.25;

// Size of Breit-Wigner threshold region, for mass selection biasing.
const double PhaseSpace::THRESHOLDSIZE  = 3.;

// Step size in optimal-mass search, for mass selection biasing.
const double PhaseSpace::THRESHOLDSTEP  = 0.2;

// Minimal rapidity range for allowed open range (in 2 -> 3).
const double PhaseSpace::YRANGEMARGIN  = 1e-6;

// Cutoff for f_e^e at x < 1 - 10^{-10} to be used in phase space selection.
// Note: the ...MIN quantities come from 1 - x_max or 1 - tau_max.
const double PhaseSpace::LEPTONXMIN     = 1e-10;
const double PhaseSpace::LEPTONXMAX     = 1. - 1e-10;
const double PhaseSpace::LEPTONXLOGMIN  = log(1e-10);
const double PhaseSpace::LEPTONXLOGMAX  = log(1. - 1e-10);
const double PhaseSpace::LEPTONTAUMIN   = 2e-10;

// Safety to avoid division with unreasonably small value for z selection.
const double PhaseSpace::SHATMINZ       = 1.;

// Regularization for small pT2min in z = cos(theta) selection.
const double PhaseSpace::PT2RATMINZ     = 0.0001;

// These numbers are hardwired empirical parameters,
// intended to speed up the M-generator.
const double PhaseSpace::WTCORRECTION[11] = { 1., 1., 1.,
  2., 5., 15., 60., 250., 1250., 7000., 50000. };

//--------------------------------------------------------------------------

// Perform simple initialization and store pointers.

void PhaseSpace::init(bool isFirst, SigmaProcess* sigmaProcessPtrIn,
  Info* infoPtrIn, Settings* settingsPtrIn, ParticleData* particleDataPtrIn,
  Rndm* rndmPtrIn, BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  Couplings* couplingsPtrIn, SigmaTotal* sigmaTotPtrIn,
  UserHooks* userHooksPtrIn) {

  // Store input pointers for future use.
  sigmaProcessPtr = sigmaProcessPtrIn;
  infoPtr         = infoPtrIn;
  settingsPtr     = settingsPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  beamAPtr        = beamAPtrIn;
  beamBPtr        = beamBPtrIn;
  couplingsPtr    = couplingsPtrIn;
  sigmaTotPtr     = sigmaTotPtrIn;
  userHooksPtr    = userHooksPtrIn;

  // Some commonly used beam information.
  idA             = beamAPtr->id();
  idB             = beamBPtr->id();
  mA              = beamAPtr->m();
  mB              = beamBPtr->m();
  eCM             = infoPtr->eCM();
  s               = eCM * eCM;

  // Flag if lepton beams, and if non-resolved ones.
  hasLeptonBeamA      = beamAPtr->isLepton();
  hasLeptonBeamB      = beamBPtr->isLepton();
  hasTwoLeptonBeams   = hasLeptonBeamA && hasLeptonBeamB;
  hasOneLeptonBeam = (hasLeptonBeamA || hasLeptonBeamB) && !hasTwoLeptonBeams;
  bool hasPointLepton = (hasLeptonBeamA && beamAPtr->isUnresolved())
                     || (hasLeptonBeamB && beamBPtr->isUnresolved());

  // Flags also for unresolved photons.
  hasPointGammaA       = beamAPtr->isGamma() && beamAPtr->isUnresolved();
  hasPointGammaB       = beamBPtr->isGamma() && beamBPtr->isUnresolved();
  hasOnePointParticle  = (hasOneLeptonBeam  && hasPointLepton)
    || ( hasPointGammaA && !hasPointGammaB)
    || (!hasPointGammaA &&  hasPointGammaB);
  hasTwoPointParticles = (hasTwoLeptonBeams && hasPointLepton)
    || ( hasPointGammaA && hasPointGammaB);

  // Flag if photons from leptons.
  bool beamHasResGamma = beamAPtr->hasResGamma() && beamBPtr->hasResGamma();

  // Set flags for (un)resolved photons according to gammaModes.
  if ( beamAPtr->isGamma() && beamBPtr->isGamma() ) {

    int beamAGammaMode = beamAPtr->getGammaMode();
    int beamBGammaMode = beamBPtr->getGammaMode();

    if ( beamAGammaMode == 2 && beamBGammaMode != 2 ) {
      hasOnePointParticle = true;
      hasPointGammaA = true;
    }
    if ( beamBGammaMode == 2 && beamAGammaMode != 2 ) {
      hasOnePointParticle = true;
      hasPointGammaB = true;
    }
    if ( beamAGammaMode == 2 && beamBGammaMode == 2 ) {
      hasTwoPointParticles = true;
      hasPointGammaA = true;
      hasPointGammaB = true;
    }
  }

  // Standard phase space cuts.
  if (isFirst || settingsPtr->flag("PhaseSpace:sameForSecond")) {
    mHatGlobalMin      = settingsPtr->parm("PhaseSpace:mHatMin");
    mHatGlobalMax      = settingsPtr->parm("PhaseSpace:mHatMax");
    pTHatGlobalMin     = settingsPtr->parm("PhaseSpace:pTHatMin");
    pTHatGlobalMax     = settingsPtr->parm("PhaseSpace:pTHatMax");

  // Optionally separate phase space cuts for second hard process.
  } else {
    mHatGlobalMin      = settingsPtr->parm("PhaseSpace:mHatMinSecond");
    mHatGlobalMax      = settingsPtr->parm("PhaseSpace:mHatMaxSecond");
    pTHatGlobalMin     = settingsPtr->parm("PhaseSpace:pTHatMinSecond");
    pTHatGlobalMax     = settingsPtr->parm("PhaseSpace:pTHatMaxSecond");
  }

  // Cutoff against divergences at pT -> 0.
  pTHatMinDiverge      = settingsPtr->parm("PhaseSpace:pTHatMinDiverge");

  // Special cut on DIS Q2 = -tHat.
  Q2GlobalMin          = settingsPtr->parm("PhaseSpace:Q2Min");
  hasQ2Min             = ( Q2GlobalMin >= pow2(pTHatMinDiverge) );

  // For photons from lepton beams match the cuts to gm+gm system cuts.
  if ( beamHasResGamma ) {
    double Wmax         = settingsPtr->parm("Photon:Wmax");
    if ( (mHatGlobalMax > Wmax) || mHatGlobalMax < 0.) mHatGlobalMax = Wmax;
  }

  // When to use Breit-Wigners.
  useBreitWigners      = settingsPtr->flag("PhaseSpace:useBreitWigners");
  minWidthBreitWigners = settingsPtr->parm("PhaseSpace:minWidthBreitWigners");
  minWidthNarrowBW     = settingsPtr->parm("PhaseSpace:minWidthNarrowBW");

  // Whether generation is with variable energy.
  doEnergySpread       = settingsPtr->flag("Beams:allowMomentumSpread")
                      || settingsPtr->flag("Beams:allowVariableEnergy");

  // Flags for maximization information and violation handling.
  showSearch           = settingsPtr->flag("PhaseSpace:showSearch");
  showViolation        = settingsPtr->flag("PhaseSpace:showViolation");
  increaseMaximum      = settingsPtr->flag("PhaseSpace:increaseMaximum");

  // Know whether a Z0 is pure Z0 or admixed with gamma*.
  gmZmodeGlobal        = settingsPtr->mode("WeakZ0:gmZmode");

  // Flags if user should be allowed to reweight cross section.
  canModifySigma   = (userHooksPtr != 0)
                   ? userHooksPtr->canModifySigma() : false;
  canBiasSelection = (userHooksPtr != 0)
                   ? userHooksPtr->canBiasSelection() : false;

  // Parameters for simplified reweighting of 2 -> 2 processes.
  canBias2Sel      = settingsPtr->flag("PhaseSpace:bias2Selection");
  bias2SelPow      = settingsPtr->parm("PhaseSpace:bias2SelectionPow");
  bias2SelRef      = settingsPtr->parm("PhaseSpace:bias2SelectionRef");
  if (canBias2Sel) pTHatGlobalMin = max( pTHatGlobalMin, pTHatMinDiverge);

  // Default event-specific kinematics properties.
  x1H             = 1.;
  x2H             = 1.;
  m3              = 0.;
  m4              = 0.;
  m5              = 0.;
  s3              = m3 * m3;
  s4              = m4 * m4;
  s5              = m5 * m5;
  mHat            = eCM;
  sH              = s;
  tH              = 0.;
  uH              = 0.;
  pTH             = 0.;
  theta           = 0.;
  phi             = 0.;
  runBW3H         = 1.;
  runBW4H         = 1.;
  runBW5H         = 1.;

  // Default cross section information.
  sigmaNw         = 0.;
  sigmaMx         = 0.;
  sigmaPos        = 0.;
  sigmaNeg        = 0.;
  newSigmaMx      = false;
  biasWt          = 1.;

}

//--------------------------------------------------------------------------

// Allow for nonisotropic decays when ME's available.

void PhaseSpace::decayKinematics( Event& process) {

  // Identify sets of sister partons.
  int iResEnd = 4;
  for (int iResBeg = 5; iResBeg < process.size(); ++iResBeg) {
    if (iResBeg <= iResEnd) continue;
    iResEnd = iResBeg;
    while ( iResEnd < process.size() - 1
      && process[iResEnd + 1].mother1() == process[iResBeg].mother1()
      && process[iResEnd + 1].mother2() == process[iResBeg].mother2() )
      ++iResEnd;

    // Check that at least one of them is a resonance.
    bool hasRes = false;
    for (int iRes = iResBeg; iRes <= iResEnd; ++iRes)
      if ( !process[iRes].isFinal() ) hasRes = true;
    if ( !hasRes ) continue;

    // Evaluate matrix element and decide whether to keep kinematics.
    double decWt = sigmaProcessPtr->weightDecay( process, iResBeg, iResEnd);
    if (decWt < 0.) infoPtr->errorMsg("Warning in PhaseSpace::decay"
      "Kinematics: negative angular weight");
    if (decWt > 1.) infoPtr->errorMsg("Warning in PhaseSpace::decay"
      "Kinematics: angular weight above unity");
    while (decWt < rndmPtr->flat() ) {

      // Find resonances for which to redo decay angles.
      for (int iRes = iResBeg; iRes < process.size(); ++iRes) {
        if ( process[iRes].isFinal() ) continue;
        int iResMother = iRes;
        while (iResMother > iResEnd)
          iResMother = process[iResMother].mother1();
        if (iResMother < iResBeg) continue;

        // Do decay of this mother isotropically in phase space.
        decayKinematicsStep( process, iRes);

      // End loop over resonance decay chains.
      }

      // Ready to allow new test of matrix element.
      decWt = sigmaProcessPtr->weightDecay( process, iResBeg, iResEnd);
      if (decWt < 0.) infoPtr->errorMsg("Warning in PhaseSpace::decay"
        "Kinematics: negative angular weight");
      if (decWt > 1.) infoPtr->errorMsg("Warning in PhaseSpace::decay"
        "Kinematics: angular weight above unity");
    }

  // End loop over sets of sister resonances/partons.
  }

}

//--------------------------------------------------------------------------

// Reselect decay products momenta isotropically in phase space.
// Does not redo secondary vertex position!

void PhaseSpace::decayKinematicsStep( Event& process, int iRes) {

   // Multiplicity and mother mass and four-momentum.
   int    i1   = process[iRes].daughter1();
   int    mult = process[iRes].daughter2() + 1 - i1;
   double m0   = process[iRes].m();
   Vec4   pRes = process[iRes].p();

  // Description of two-body decays as simple special case.
  if (mult == 2) {

    // Products and product masses.
    int    i2   = i1 + 1;
    double m1t  = process[i1].m();
    double m2t  = process[i2].m();

    // Energies and absolute momentum in the rest frame.
    double e1   = 0.5 * (m0*m0 + m1t*m1t - m2t*m2t) / m0;
    double e2   = 0.5 * (m0*m0 + m2t*m2t - m1t*m1t) / m0;
    double p12  = 0.5 * sqrtpos( (m0 - m1t - m2t) * (m0 + m1t + m2t)
      * (m0 + m1t - m2t) * (m0 - m1t + m2t) ) / m0;

    // Pick isotropic angles to give three-momentum.
    double cosTheta = 2. * rndmPtr->flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi12    = 2. * M_PI * rndmPtr->flat();
    double pX       = p12 * sinTheta * cos(phi12);
    double pY       = p12 * sinTheta * sin(phi12);
    double pZ       = p12 * cosTheta;

    // Fill four-momenta in mother rest frame and then boost to lab frame.
    Vec4 p1(  pX,  pY,  pZ, e1);
    Vec4 p2( -pX, -pY, -pZ, e2);
    p1.bst( pRes );
    p2.bst( pRes );

    // Done for two-body decay.
    process[i1].p( p1 );
    process[i2].p( p2 );
    return;
  }

  // Description of three-body decays as semi-simple special case.
  if (mult == 3) {

    // Products and product masses.
    int    i2      = i1 + 1;
    int    i3      = i2 + 1;
    double m1t     = process[i1].m();
    double m2t     = process[i2].m();
    double m3t     = process[i3].m();
    double mDiff   = m0 - (m1t + m2t + m3t);

    // Kinematical limits for 2+3 mass. Maximum phase-space weight.
    double m23Min  = m2t + m3t;
    double m23Max  = m0 - m1t;
    double p1Max   = 0.5 * sqrtpos( (m0 - m1t - m23Min)
      * (m0 + m1t + m23Min) * (m0 + m1t - m23Min)
      * (m0 - m1t + m23Min) ) / m0;
    double p23Max  = 0.5 * sqrtpos( (m23Max - m2t - m3t)
      * (m23Max + m2t + m3t) * (m23Max + m2t - m3t)
      * (m23Max - m2t + m3t) ) / m23Max;
    double wtPSmax = 0.5 * p1Max * p23Max;

    // Pick an intermediate mass m23 flat in the allowed range.
    double wtPS, m23, p1Abs, p23Abs;
    do {
      m23 = m23Min + rndmPtr->flat() * mDiff;

      // Translate into relative momenta and find phase-space weight.
      p1Abs  = 0.5 * sqrtpos( (m0 - m1t - m23) * (m0 + m1t + m23)
        * (m0 + m1t - m23) * (m0 - m1t + m23) ) / m0;
      p23Abs = 0.5 * sqrtpos( (m23 - m2t - m3t) * (m23 + m2t + m3t)
        * (m23 + m2t - m3t) * (m23 - m2t + m3t) ) / m23;
      wtPS   = p1Abs * p23Abs;

    // If rejected, try again with new invariant masses.
    } while ( wtPS < rndmPtr->flat() * wtPSmax );

    // Set up m23 -> m2 + m3 isotropic in its rest frame.
    double cosTheta = 2. * rndmPtr->flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi23    = 2. * M_PI * rndmPtr->flat();
    double pX       = p23Abs * sinTheta * cos(phi23);
    double pY       = p23Abs * sinTheta * sin(phi23);
    double pZ       = p23Abs * cosTheta;
    double e2       = sqrt( m2t*m2t + p23Abs*p23Abs);
    double e3       = sqrt( m3t*m3t + p23Abs*p23Abs);
    Vec4 p2(  pX,  pY,  pZ, e2);
    Vec4 p3( -pX, -pY, -pZ, e3);

    // Set up 0 -> 1 + 23 isotropic in its rest frame.
    cosTheta        = 2. * rndmPtr->flat() - 1.;
    sinTheta        = sqrt(1. - cosTheta*cosTheta);
    phi23           = 2. * M_PI * rndmPtr->flat();
    pX              = p1Abs * sinTheta * cos(phi23);
    pY              = p1Abs * sinTheta * sin(phi23);
    pZ              = p1Abs * cosTheta;
    double e1       = sqrt( m1t*m1t + p1Abs*p1Abs);
    double e23      = sqrt( m23*m23 + p1Abs*p1Abs);
    Vec4 p1( pX, pY, pZ, e1);

    // Boost 2 + 3 to the 0 rest frame and then boost to lab frame.
    Vec4 p23( -pX, -pY, -pZ, e23);
    p2.bst( p23 );
    p3.bst( p23 );
    p1.bst( pRes );
    p2.bst( pRes );
    p3.bst( pRes );

    // Done for three-body decay.
    process[i1].p( p1 );
    process[i2].p( p2 );
    process[i3].p( p3 );
    return;
  }

  // Do a multibody decay using the M-generator algorithm.

  // Set up masses and four-momenta in a vector, with mother in slot 0.
  vector<double> mProd;
  mProd.push_back( m0);
  for (int i = i1; i <= process[iRes].daughter2(); ++i)
    mProd.push_back( process[i].m() );
  vector<Vec4> pProd;
  pProd.push_back( pRes);

  // Sum of daughter masses.
  double mSum    = mProd[1];
  for (int i = 2; i <= mult; ++i) mSum += mProd[i];
  double mDiff   = m0 - mSum;

  // Begin setup of intermediate invariant masses.
  vector<double> mInv;
  for (int i = 0; i <= mult; ++i) mInv.push_back( mProd[i]);

  // Calculate the maximum weight in the decay.
  double wtPSmax = 1. / WTCORRECTION[mult];
  double mMaxWT  = mDiff + mProd[mult];
  double mMinWT  = 0.;
  for (int i = mult - 1; i > 0; --i) {
    mMaxWT      += mProd[i];
    mMinWT      += mProd[i+1];
    double mNow  = mProd[i];
    wtPSmax *= 0.5 * sqrtpos( (mMaxWT - mMinWT - mNow)
      * (mMaxWT + mMinWT + mNow) * (mMaxWT + mMinWT - mNow)
      * (mMaxWT - mMinWT + mNow) ) / mMaxWT;
  }

  // Begin loop to find the set of intermediate invariant masses.
  vector<double> rndmOrd;
  double wtPS;
  do {
    wtPS  = 1.;

    // Find and order random numbers in descending order.
    rndmOrd.resize(0);
    rndmOrd.push_back(1.);
    for (int i = 1; i < mult - 1; ++i) {
      double rndm = rndmPtr->flat();
      rndmOrd.push_back(rndm);
      for (int j = i - 1; j > 0; --j) {
        if (rndm > rndmOrd[j]) swap( rndmOrd[j], rndmOrd[j+1] );
        else break;
      }
    }
    rndmOrd.push_back(0.);

    // Translate into intermediate masses and find weight.
    for (int i = mult - 1; i > 0; --i) {
      mInv[i] = mInv[i+1] + mProd[i] + (rndmOrd[i-1] - rndmOrd[i]) * mDiff;
      wtPS   *= 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - mProd[i])
        * (mInv[i] + mInv[i+1] + mProd[i]) * (mInv[i] + mInv[i+1] - mProd[i])
        * (mInv[i] - mInv[i+1] + mProd[i]) ) / mInv[i];
    }

  // If rejected, try again with new invariant masses.
  } while ( wtPS < rndmPtr->flat() * wtPSmax );

  // Perform two-particle decays in the respective rest frame.
  vector<Vec4> pInv;
  pInv.resize(mult + 1);
  for (int i = 1; i < mult; ++i) {
    double p12 = 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - mProd[i])
      * (mInv[i] + mInv[i+1] + mProd[i]) * (mInv[i] + mInv[i+1] - mProd[i])
      * (mInv[i] - mInv[i+1] + mProd[i]) ) / mInv[i];

    // Isotropic angles give three-momentum.
    double cosTheta = 2. * rndmPtr->flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phiLoc   = 2. * M_PI * rndmPtr->flat();
    double pX       = p12 * sinTheta * cos(phiLoc);
    double pY       = p12 * sinTheta * sin(phiLoc);
    double pZ       = p12 * cosTheta;

    // Calculate energies, fill four-momenta.
    double eHad     = sqrt( mProd[i]*mProd[i] + p12*p12);
    double eInv     = sqrt( mInv[i+1]*mInv[i+1] + p12*p12);
    pProd.push_back( Vec4( pX, pY, pZ, eHad) );
    pInv[i+1].p( -pX, -pY, -pZ, eInv);
  }
  pProd.push_back( pInv[mult] );

  // Boost decay products to the mother rest frame and on to lab frame.
  pInv[1] = pProd[0];
  for (int iFrame = mult - 1; iFrame > 0; --iFrame)
    for (int i = iFrame; i <= mult; ++i) pProd[i].bst(pInv[iFrame]);

  // Done for multibody decay.
  for (int i = 1; i <= mult; ++i)
    process[i1 + i - 1].p( pProd[i] );
  return;

}

//--------------------------------------------------------------------------

// Determine how 3-body phase space should be sampled.

void PhaseSpace::setup3Body() {

  // Check for massive t-channel propagator particles.
  int idTchan1    = abs( sigmaProcessPtr->idTchan1() );
  int idTchan2    = abs( sigmaProcessPtr->idTchan2() );
  mTchan1         = (idTchan1 == 0) ? pTHatMinDiverge
                                    : particleDataPtr->m0(idTchan1);
  mTchan2         = (idTchan2 == 0) ? pTHatMinDiverge
                                    : particleDataPtr->m0(idTchan2);
  sTchan1         = mTchan1 * mTchan1;
  sTchan2         = mTchan2 * mTchan2;

  // Find coefficients of different pT2 selection terms. Mirror choice.
  frac3Pow1       = sigmaProcessPtr->tChanFracPow1();
  frac3Pow2       = sigmaProcessPtr->tChanFracPow2();
  frac3Flat       = 1. - frac3Pow1 - frac3Pow2;
  useMirrorWeight = sigmaProcessPtr->useMirrorWeight();

}

//--------------------------------------------------------------------------

// Determine how phase space should be sampled.

bool PhaseSpace::setupSampling123(bool is2, bool is3) {

  // Optional printout.
  if (showSearch) cout <<  "\n PYTHIA Optimization printout for "
    << sigmaProcessPtr->name() << "\n \n" << scientific << setprecision(3);

  // Check that open range in tau (+ set tauMin, tauMax).
  if (!limitTau(is2, is3)) return false;

  // Reset coefficients and matrices of equation system to solve.
  int binTau[8], binY[8], binZ[8];
  double vecTau[8], matTau[8][8], vecY[8], matY[8][8], vecZ[8], matZ[8][8];
  for (int i = 0; i < 8; ++i) {
    tauCoef[i] = 0.;
    yCoef[i]   = 0.;
    zCoef[i]   = 0.;
    binTau[i]  = 0;
    binY[i]    = 0;
    binZ[i]    = 0;
    vecTau[i]  = 0.;
    vecY[i]    = 0.;
    vecZ[i]    = 0.;
    for (int j = 0; j < 8; ++j) {
      matTau[i][j] = 0.;
      matY[i][j]   = 0.;
      matZ[i][j]   = 0.;
    }
  }
  sigmaMx  = 0.;
  sigmaNeg = 0.;

  // Number of used coefficients/points for each dimension: tau, y, c.
  nTau = (hasTwoPointParticles) ? 1 : 2;
  nY   = (hasOnePointParticle || hasTwoPointParticles) ? 1 : 5;
  nZ   = (is2) ? 5 : 1;

  // Identify if any resonances contribute in s-channel.
  idResA = sigmaProcessPtr->resonanceA();
  if (idResA != 0) {
     mResA = particleDataPtr->m0(idResA);
     GammaResA = particleDataPtr->mWidth(idResA);
     if (mHatMin > mResA + WIDTHMARGIN * GammaResA || (mHatMax > 0.
       && mHatMax < mResA - WIDTHMARGIN * GammaResA) ) idResA = 0;
  }
  idResB = sigmaProcessPtr->resonanceB();
  if (idResB != 0) {
     mResB = particleDataPtr->m0(idResB);
     GammaResB = particleDataPtr->mWidth(idResB);
     if (mHatMin > mResB + WIDTHMARGIN * GammaResB || (mHatMax > 0.
       && mHatMax < mResB - WIDTHMARGIN * GammaResB) ) idResB = 0;
  }
  if (idResA == 0 && idResB != 0) {
    idResA = idResB;
    mResA = mResB;
    GammaResA = GammaResB;
    idResB = 0;
  }

  // More sampling in tau if resonances in s-channel.
  if (idResA !=0 && !hasTwoPointParticles) {
    nTau += 2;
    tauResA = mResA * mResA / s;
    widResA = mResA * GammaResA / s;
  }
  if (idResB != 0 && !hasTwoPointParticles) {
    nTau += 2;
    tauResB = mResB * mResB / s;
    widResB = mResB * GammaResB / s;
  }

  // More sampling in tau (and different in y) if incoming lepton beams.
  if (hasTwoLeptonBeams && !hasTwoPointParticles) ++nTau;

  // Special case when both resonances have same mass.
  sameResMass = false;
  if (idResB != 0 && abs(mResA - mResB) < SAMEMASS * (GammaResA + GammaResB))
    sameResMass = true;

  // Default z value and weight required for 2 -> 1. Number of dimensions.
  z = 0.;
  wtZ = 1.;
  int nVar = (is2) ? 3 : 2;

  // Initial values, to be modified later.
  tauCoef[0] = 1.;
  yCoef[1]   = 0.5;
  yCoef[2]   = 0.5;
  zCoef[0]   = 1.;

  // Step through grid in tau. Set limits on y and z generation.
  for (int iTau = 0; iTau < nTau; ++iTau) {
    double posTau = 0.5;
    if (sameResMass && iTau > 1 && iTau < 6) posTau = (iTau < 4) ? 0.4 : 0.6;
    selectTau( iTau, posTau, is2);
    if (!limitY()) continue;
    if (is2 && !limitZ()) continue;

    // Step through grids in y and z.
    for (int iY = 0; iY < nY; ++iY) {
      selectY( iY, 0.5);
      for (int iZ = 0; iZ < nZ; ++iZ) {
        if (is2) selectZ( iZ, 0.5);
        double sigmaTmp = 0.;

        // 2 -> 1: calculate cross section, weighted by phase-space volume.
        if (!is2 && !is3) {
          sigmaProcessPtr->set1Kin( x1H, x2H, sH);
          sigmaTmp = sigmaProcessPtr->sigmaPDF(true);
          sigmaTmp *= wtTau * wtY;

        // 2 -> 2: calculate cross section, weighted by phase-space volume
        // and Breit-Wigners for masses
        } else if (is2) {
          sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4,
            runBW3H, runBW4H);
          sigmaTmp = sigmaProcessPtr->sigmaPDF(true);
          sigmaTmp *= wtTau * wtY * wtZ * wtBW;

        // 2 -> 3: repeat internal 3-body phase space several times and
        // keep maximal cross section, weighted by phase-space volume
        // and Breit-Wigners for masses
        } else if (is3) {
          for (int iTry3 = 0; iTry3 < NTRY3BODY; ++iTry3) {
            if (!select3Body()) continue;
            sigmaProcessPtr->set3Kin( x1H, x2H, sH, p3cm, p4cm, p5cm,
              m3, m4, m5, runBW3H, runBW4H, runBW5H);
            double sigmaTry = sigmaProcessPtr->sigmaPDF(true);
            sigmaTry *= wtTau * wtY * wt3Body * wtBW;
            if (sigmaTry > sigmaTmp) sigmaTmp = sigmaTry;
          }
        }

        // Allow possibility for user to modify cross section. (3body??)
        if (canModifySigma) sigmaTmp
           *= userHooksPtr->multiplySigmaBy( sigmaProcessPtr, this, false);
        if (canBiasSelection) sigmaTmp
           *= userHooksPtr->biasSelectionBy( sigmaProcessPtr, this, false);
        if (canBias2Sel) sigmaTmp *= pow( pTH / bias2SelRef, bias2SelPow);

        // Check if current maximum exceeded.
        if (sigmaTmp > sigmaMx) sigmaMx = sigmaTmp;

        // Optional printout. Protect against negative cross sections.
        if (showSearch) cout << " tau =" << setw(11) << tau << "  y ="
          << setw(11) << y << "  z =" << setw(11) << z
          << "  sigma =" << setw(11) << sigmaTmp << "\n";
        if (sigmaTmp < 0.) sigmaTmp = 0.;

        // Sum up tau cross-section pieces in points used.
        if (!hasTwoPointParticles) {
          binTau[iTau]      += 1;
          vecTau[iTau]      += sigmaTmp;
          matTau[iTau][0]   += 1. / intTau0;
          matTau[iTau][1]   += (1. / intTau1) / tau;
          if (idResA != 0) {
            matTau[iTau][2] += (1. / intTau2) / (tau + tauResA);
            matTau[iTau][3] += (1. / intTau3)
              * tau / ( pow2(tau - tauResA) + pow2(widResA) );
          }
          if (idResB != 0) {
            matTau[iTau][4] += (1. / intTau4) / (tau + tauResB);
            matTau[iTau][5] += (1. / intTau5)
              * tau / ( pow2(tau - tauResB) + pow2(widResB) );
          }
          if (hasTwoLeptonBeams) matTau[iTau][nTau - 1] += (1. / intTau6)
              * tau / max( LEPTONTAUMIN, 1. - tau);
        }

        // Sum up y cross-section pieces in points used.
        if (!hasOnePointParticle && !hasTwoPointParticles) {
          binY[iY]      += 1;
          vecY[iY]      += sigmaTmp;
          matY[iY][0]   += (yMax / intY0) / cosh(y);
          matY[iY][1]   += (yMax / intY12) * (y + yMax);
          matY[iY][2]   += (yMax / intY12) * (yMax - y);
          if (!hasTwoLeptonBeams) {
            matY[iY][3] += (yMax / intY34) * exp(y);
            matY[iY][4] += (yMax / intY34) * exp(-y);
          } else {
            matY[iY][3] += (yMax / intY56)
              / max( LEPTONXMIN, 1. - exp( y - yMax) );
            matY[iY][4] += (yMax / intY56)
              / max( LEPTONXMIN, 1. - exp(-y - yMax) );
          }
        }

        // Integrals over z expressions at tauMax, to be used below.
        if (is2) {
          double p2AbsMax   = 0.25 * (pow2(tauMax * s - s3 - s4)
            - 4. * s3 * s4) / (tauMax * s);
          double zMaxMax    = sqrtpos( 1. - pT2HatMin / p2AbsMax );
          double zPosMaxMax = max(ratio34, unity34 + zMaxMax);
          double zNegMaxMax = max(ratio34, unity34 - zMaxMax);
          double intZ0Max   = 2. * zMaxMax;
          double intZ12Max  = log( zPosMaxMax / zNegMaxMax);
          double intZ34Max  = 1. / zNegMaxMax - 1. / zPosMaxMax;

          // Sum up z cross-section pieces in points used.
          binZ[iZ]    += 1;
          vecZ[iZ]    += sigmaTmp;
          matZ[iZ][0] += 1.;
          matZ[iZ][1] += (intZ0Max / intZ12Max) / zNeg;
          matZ[iZ][2] += (intZ0Max / intZ12Max) / zPos;
          matZ[iZ][3] += (intZ0Max / intZ34Max) / pow2(zNeg);
          matZ[iZ][4] += (intZ0Max / intZ34Max) / pow2(zPos);
        }

      // End of loops over phase space points.
      }
    }
  }

  // Fail if no non-vanishing cross sections.
  if (sigmaMx <= 0.) {
    sigmaMx = 0.;
    return false;
  }


  // Solve respective equation system for better phase space coefficients.
  if (!hasTwoPointParticles) solveSys( nTau, binTau, vecTau, matTau, tauCoef);
  if (!hasOnePointParticle && !hasTwoPointParticles)
    solveSys( nY, binY, vecY, matY, yCoef);
  if (is2) solveSys( nZ, binZ, vecZ, matZ, zCoef);
  if (showSearch) cout << "\n";

  // Provide cumulative sum of coefficients.
  tauCoefSum[0] = tauCoef[0];
    yCoefSum[0] =   yCoef[0];
    zCoefSum[0] =   zCoef[0];
  for (int i = 1; i < 8; ++ i) {
    tauCoefSum[i] = tauCoefSum[i - 1] + tauCoef[i];
      yCoefSum[i] =   yCoefSum[i - 1] +   yCoef[i];
      zCoefSum[i] =   zCoefSum[i - 1] +   zCoef[i];
  }
  // The last element should be > 1 to be on safe side in selection below.
  tauCoefSum[nTau - 1] = 2.;
    yCoefSum[nY   - 1] = 2.;
    zCoefSum[nZ   - 1] = 2.;


  // Begin find two most promising maxima among same points as before.
  int iMaxTau[NMAXTRY + 2], iMaxY[NMAXTRY + 2], iMaxZ[NMAXTRY + 2];
  double sigMax[NMAXTRY + 2];
  int nMax = 0;

  // Scan same grid as before in tau, y, z.
  for (int iTau = 0; iTau < nTau; ++iTau) {
    double posTau = 0.5;
    if (sameResMass && iTau > 1 && iTau < 6) posTau = (iTau < 4) ? 0.4 : 0.6;
    selectTau( iTau, posTau, is2);
    if (!limitY()) continue;
    if (is2 && !limitZ()) continue;
    for (int iY = 0; iY < nY; ++iY) {
      selectY( iY, 0.5);
      for (int iZ = 0; iZ < nZ; ++iZ) {
        if (is2) selectZ( iZ, 0.5);
        double sigmaTmp = 0.;

        // 2 -> 1: calculate cross section, weighted by phase-space volume.
        if (!is2 && !is3) {
          sigmaProcessPtr->set1Kin( x1H, x2H, sH);
          sigmaTmp = sigmaProcessPtr->sigmaPDF(true);
          sigmaTmp *= wtTau * wtY;

        // 2 -> 2: calculate cross section, weighted by phase-space volume
        // and Breit-Wigners for masses
        } else if (is2) {
          sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4,
            runBW3H, runBW4H);
          sigmaTmp = sigmaProcessPtr->sigmaPDF(true);
          sigmaTmp *= wtTau * wtY * wtZ * wtBW;

        // 2 -> 3: repeat internal 3-body phase space several times and
        // keep maximal cross section, weighted by phase-space volume
        // and Breit-Wigners for masses
        } else if (is3) {
          for (int iTry3 = 0; iTry3 < NTRY3BODY; ++iTry3) {
            if (!select3Body()) continue;
            sigmaProcessPtr->set3Kin( x1H, x2H, sH, p3cm, p4cm, p5cm,
              m3, m4, m5, runBW3H, runBW4H, runBW5H);
            double sigmaTry = sigmaProcessPtr->sigmaPDF(true);
            sigmaTry *= wtTau * wtY * wt3Body * wtBW;
            if (sigmaTry > sigmaTmp) sigmaTmp = sigmaTry;
          }
        }

        // Allow possibility for user to modify cross section. (3body??)
        if (canModifySigma) sigmaTmp
           *= userHooksPtr->multiplySigmaBy( sigmaProcessPtr, this, false);
        if (canBiasSelection) sigmaTmp
           *= userHooksPtr->biasSelectionBy( sigmaProcessPtr, this, false);
        if (canBias2Sel) sigmaTmp *= pow( pTH / bias2SelRef, bias2SelPow);

        // Optional printout. Protect against negative cross section.
        if (showSearch) cout << " tau =" << setw(11) << tau << "  y ="
          << setw(11) << y << "  z =" << setw(11) << z
          << "  sigma =" << setw(11) << sigmaTmp << "\n";
        if (sigmaTmp < 0.) sigmaTmp = 0.;

        // Check that point is not simply mirror of already found one.
        bool mirrorPoint = false;
        for (int iMove = 0; iMove < nMax; ++iMove)
          if (abs(sigmaTmp - sigMax[iMove]) < SAMESIGMA
            * (sigmaTmp + sigMax[iMove])) mirrorPoint = true;

        // Add to or insert in maximum list. Only first two count.
        if (!mirrorPoint) {
          int iInsert = 0;
          for (int iMove = nMax - 1; iMove >= -1; --iMove) {
            iInsert = iMove + 1;
            if (iInsert == 0 || sigmaTmp < sigMax[iMove]) break;
            iMaxTau[iMove + 1] = iMaxTau[iMove];
            iMaxY[iMove + 1] = iMaxY[iMove];
            iMaxZ[iMove + 1] = iMaxZ[iMove];
            sigMax[iMove + 1] = sigMax[iMove];
          }
          iMaxTau[iInsert] = iTau;
          iMaxY[iInsert] = iY;
          iMaxZ[iInsert] = iZ;
          sigMax[iInsert] = sigmaTmp;
          if (nMax < NMAXTRY) ++nMax;
        }

      // Found two most promising maxima.
      }
    }
  }
  if (showSearch) cout << "\n";

  // Read out starting position for search.
  sigmaMx = sigMax[0];
  int beginVar = (hasTwoPointParticles) ? 2 : 0;
  if (hasOnePointParticle) beginVar = 1;
  for (int iMax = 0; iMax < nMax; ++iMax) {
    int iTau = iMaxTau[iMax];
    int iY = iMaxY[iMax];
    int iZ = iMaxZ[iMax];
    double tauVal = 0.5;
    double yVal = 0.5;
    double zVal = 0.5;
    int iGrid;
    double varVal, varNew, deltaVar, marginVar, sigGrid[3];

    // Starting point and step size in parameter space.
    for (int iRepeat = 0; iRepeat < 2; ++iRepeat) {
      // Run through (possibly a subset of) tau, y and z.
      for (int iVar = beginVar; iVar < nVar; ++iVar) {
        bool isTauVar = iVar == 0 || (beginVar == 1 && iVar == 1);
        if (isTauVar) varVal = tauVal;
        else if (iVar == 1) varVal = yVal;
        else varVal = zVal;
        deltaVar = (iRepeat == 0) ? 0.1
          : max( 0.01, min( 0.05, min( varVal - 0.02, 0.98 - varVal) ) );
        marginVar = (iRepeat == 0) ? 0.02 : 0.002;
        int moveStart = (iRepeat == 0 && isTauVar) ? 0 : 1;
        for (int move = moveStart; move < 9; ++move) {

          // Define new parameter-space point by step in one dimension.
          if (move == 0) {
            iGrid = 1;
            varNew = varVal;
          } else if (move == 1) {
            iGrid = 2;
            varNew = varVal + deltaVar;
          } else if (move == 2) {
            iGrid = 0;
            varNew = varVal - deltaVar;
          } else if (sigGrid[2] >= max( sigGrid[0], sigGrid[1])
            && varVal + 2. * deltaVar < 1. - marginVar) {
            varVal += deltaVar;
            sigGrid[0] = sigGrid[1];
            sigGrid[1] = sigGrid[2];
            iGrid = 2;
            varNew = varVal + deltaVar;
          } else if (sigGrid[0] >= max( sigGrid[1], sigGrid[2])
            && varVal - 2. * deltaVar > marginVar) {
            varVal -= deltaVar;
            sigGrid[2] = sigGrid[1];
            sigGrid[1] = sigGrid[0];
            iGrid = 0;
            varNew = varVal - deltaVar;
          } else if (sigGrid[2] >= sigGrid[0]) {
            deltaVar *= 0.5;
            varVal += deltaVar;
            sigGrid[0] = sigGrid[1];
            iGrid = 1;
            varNew = varVal;
          } else {
            deltaVar *= 0.5;
            varVal -= deltaVar;
            sigGrid[2] = sigGrid[1];
            iGrid = 1;
            varNew = varVal;
          }

          // Convert to relevant variables and find derived new limits.
          bool insideLimits = true;
          if (isTauVar) {
            tauVal = varNew;
            selectTau( iTau, tauVal, is2);
            if (!limitY()) insideLimits = false;
            if (is2 && !limitZ()) insideLimits = false;
            if (insideLimits) {
              selectY( iY, yVal);
              if (is2) selectZ( iZ, zVal);
            }
          } else if (iVar == 1) {
            yVal = varNew;
            selectY( iY, yVal);
          } else if (iVar == 2) {
            zVal = varNew;
            selectZ( iZ, zVal);
          }

          // Evaluate cross-section.
          double sigmaTmp = 0.;
          if (insideLimits) {

            // 2 -> 1: calculate cross section, weighted by phase-space volume.
            if (!is2 && !is3) {
              sigmaProcessPtr->set1Kin( x1H, x2H, sH);
              sigmaTmp = sigmaProcessPtr->sigmaPDF(true);
              sigmaTmp *= wtTau * wtY;

            // 2 -> 2: calculate cross section, weighted by phase-space volume
            // and Breit-Wigners for masses
            } else if (is2) {
              sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4,
                runBW3H, runBW4H);
              sigmaTmp = sigmaProcessPtr->sigmaPDF(true);
              sigmaTmp *= wtTau * wtY * wtZ * wtBW;

            // 2 -> 3: repeat internal 3-body phase space several times and
            // keep maximal cross section, weighted by phase-space volume
            // and Breit-Wigners for masses
            } else if (is3) {
              for (int iTry3 = 0; iTry3 < NTRY3BODY; ++iTry3) {
                if (!select3Body()) continue;
                sigmaProcessPtr->set3Kin( x1H, x2H, sH, p3cm, p4cm, p5cm,
                  m3, m4, m5, runBW3H, runBW4H, runBW5H);
                double sigmaTry = sigmaProcessPtr->sigmaPDF(true);
                sigmaTry *= wtTau * wtY * wt3Body * wtBW;
                if (sigmaTry > sigmaTmp) sigmaTmp = sigmaTry;
              }
            }

            // Allow possibility for user to modify cross section.
            if (canModifySigma) sigmaTmp
              *= userHooksPtr->multiplySigmaBy( sigmaProcessPtr, this, false);
            if (canBiasSelection) sigmaTmp
              *= userHooksPtr->biasSelectionBy( sigmaProcessPtr, this, false);
            if (canBias2Sel) sigmaTmp *= pow( pTH / bias2SelRef, bias2SelPow);

            // Optional printout. Protect against negative cross section.
            if (showSearch) cout << " tau =" << setw(11) << tau << "  y ="
              << setw(11) << y << "  z =" << setw(11) << z
              << "  sigma =" << setw(11) << sigmaTmp << "\n";
            if (sigmaTmp < 0.) sigmaTmp = 0.;
          }

          // Save new maximum. Final maximum.
          sigGrid[iGrid] = sigmaTmp;
          if (sigmaTmp > sigmaMx) sigmaMx = sigmaTmp;
        }
      }
    }
  }
  sigmaMx *= SAFETYMARGIN;
  sigmaPos = sigmaMx;

  // Optional printout.
  if (showSearch) cout << "\n Final maximum = "  << setw(11) << sigmaMx
    << endl;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Select a trial kinematics phase space point.
// Note: by In is meant the integral over the quantity multiplying
// coefficient cn. The sum of cn is normalized to unity.

bool PhaseSpace::trialKin123(bool is2, bool is3, bool inEvent) {

  // Allow for possibility that energy varies from event to event.
  if (doEnergySpread) {
    eCM       = infoPtr->eCM();
    s         = eCM * eCM;

    // Find shifted tauRes values.
    if (idResA !=0 && !hasTwoPointParticles) {
      tauResA = mResA * mResA / s;
      widResA = mResA * GammaResA / s;
    }
    if (idResB != 0 && !hasTwoPointParticles) {
      tauResB = mResB * mResB / s;
      widResB = mResB * GammaResB / s;
    }
  }

  // Choose tau according to h1(tau)/tau, where
  // h1(tau) = c0/I0 + (c1/I1) * 1/tau
  // + (c2/I2) / (tau + tauResA)
  // + (c3/I3) * tau / ((tau - tauResA)^2 + widResA^2)
  // + (c4/I4) / (tau + tauResB)
  // + (c5/I5) * tau / ((tau - tauResB)^2 + widResB^2)
  // + (c6/I6) * tau / (1 - tau).
  if (!limitTau(is2, is3)) return false;
  int iTau = 0;
  if (!hasTwoPointParticles) {
    double rTau = rndmPtr->flat();
    while (rTau > tauCoefSum[iTau]) ++iTau;
  }
  selectTau( iTau, rndmPtr->flat(), is2);

  // Choose y according to h2(y), where
  // h2(y) = (c0/I0) * 1/cosh(y)
  // + (c1/I1) * (y-ymin) + (c2/I2) * (ymax-y)
  // + (c3/I3) * exp(y) + (c4/i4) * exp(-y) (for hadron; for lepton instead)
  // + (c5/I5) * 1 / (1 - exp(y-ymax)) + (c6/I6) * 1 / (1 - exp(ymin-y)).
  if (!limitY()) return false;
  int iY = 0;
  if (!hasOnePointParticle && !hasTwoPointParticles) {
    double rY = rndmPtr->flat();
    while (rY > yCoefSum[iY]) ++iY;
  }
  selectY( iY, rndmPtr->flat());

  // Choose z = cos(thetaHat) according to h3(z), where
  // h3(z) = c0/I0 + (c1/I1) * 1/(A - z) + (c2/I2) * 1/(A + z)
  // + (c3/I3) * 1/(A - z)^2 + (c4/I4) * 1/(A + z)^2,
  // where A = 1 + 2*(m3*m4/sH)^2 (= 1 for massless products).
  if (is2) {
    if (!limitZ()) return false;
    int iZ = 0;
    double rZ = rndmPtr->flat();
    while (rZ > zCoefSum[iZ]) ++iZ;
    selectZ( iZ, rndmPtr->flat());
  }

  // 2 -> 1: calculate cross section, weighted by phase-space volume.
  if (!is2 && !is3) {
    sigmaProcessPtr->set1Kin( x1H, x2H, sH);
    sigmaNw  = sigmaProcessPtr->sigmaPDF();
    sigmaNw *= wtTau * wtY;

  // 2 -> 2: calculate cross section, weighted by phase-space volume
  // and Breit-Wigners for masses
  } else if (is2) {
    sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4, runBW3H, runBW4H);
    sigmaNw  = sigmaProcessPtr->sigmaPDF();
    sigmaNw *= wtTau * wtY * wtZ * wtBW;

  // 2 -> 3: also sample internal 3-body phase, weighted by
  // 2 -> 1 phase-space volume and Breit-Wigners for masses
  } else if (is3) {
    if (!select3Body()) sigmaNw = 0.;
    else {
      sigmaProcessPtr->set3Kin( x1H, x2H, sH, p3cm, p4cm, p5cm,
         m3, m4, m5, runBW3H, runBW4H, runBW5H);
      sigmaNw  = sigmaProcessPtr->sigmaPDF();
      sigmaNw *= wtTau * wtY * wt3Body * wtBW;
    }
  }

  // Allow possibility for user to modify cross section.
  if (canModifySigma) sigmaNw
    *= userHooksPtr->multiplySigmaBy( sigmaProcessPtr, this, inEvent);
  if (canBiasSelection) sigmaNw
    *= userHooksPtr->biasSelectionBy( sigmaProcessPtr, this, inEvent);
  if (canBias2Sel) sigmaNw *= pow( pTH / bias2SelRef, bias2SelPow);

  // Check if maximum violated.
  newSigmaMx = false;
  if (sigmaNw > sigmaMx) {
    infoPtr->errorMsg("Warning in PhaseSpace2to2tauyz::trialKin: "
      "maximum for cross section violated");

    // Violation strategy 1: increase maximum (always during initialization).
    if (increaseMaximum || !inEvent) {
      double violFact = SAFETYMARGIN * sigmaNw / sigmaMx;
      sigmaMx = SAFETYMARGIN * sigmaNw;
      newSigmaMx = true;
      if (showViolation) {
        if (violFact < 9.99) cout << fixed;
        else                 cout << scientific;
        cout << " PYTHIA Maximum for " << sigmaProcessPtr->name()
             << " increased by factor " << setprecision(3) << violFact
             << " to " << scientific << sigmaMx << endl;
      }

    // Violation strategy 2: weight event (done in ProcessContainer).
    } else if (showViolation && sigmaNw > sigmaPos) {
      double violFact = sigmaNw / sigmaMx;
      if (violFact < 9.99) cout << fixed;
      else                 cout << scientific;
      cout << " PYTHIA Maximum for " << sigmaProcessPtr->name()
           << " exceeded by factor " << setprecision(3) << violFact << endl;
      sigmaPos = sigmaNw;
    }
  }

  // Check if negative cross section.
  if (sigmaNw < sigmaNeg) {
    infoPtr->errorMsg("Warning in PhaseSpace2to2tauyz::trialKin:"
      " negative cross section set 0", "for " +  sigmaProcessPtr->name() );
    sigmaNeg = sigmaNw;

    // Optional printout of (all) violations.
    if (showViolation) cout << " PYTHIA Negative minimum for "
      << sigmaProcessPtr->name() << " changed to " << scientific
      << setprecision(3) << sigmaNeg << endl;
  }
  if (sigmaNw < 0.) sigmaNw = 0.;

  // Set event weight, where relevant.
  biasWt = (canBiasSelection) ? userHooksPtr->biasedSelectionWeight() : 1.;
  if (canBias2Sel) biasWt /= pow( pTH / bias2SelRef, bias2SelPow);

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Find range of allowed tau values.

bool PhaseSpace::limitTau(bool is2, bool is3) {

  // Trivial reply for unresolved lepton beams.
  if (hasTwoPointParticles) {
    tauMin = 1.;
    tauMax = 1.;
    return true;
  }

  // Requirements from allowed mHat range and allowed Q2Min.
  tauMin = sHatMin / s;
  if (is2 && hasQ2Min && Q2GlobalMin + s3 + s4 > sHatMin)
    tauMin = (Q2GlobalMin + s3 + s4) / s;
  tauMax = (mHatMax < mHatMin) ? 1. : min( 1., sHatMax / s);

  // Requirements from allowed pT range and masses.
  if (is2 || is3) {
    double mT3Min = sqrt(s3 + pT2HatMin);
    double mT4Min = sqrt(s4 + pT2HatMin);
    double mT5Min = (is3) ? sqrt(s5 + pT2HatMin) : 0.;
    tauMin = max( tauMin, pow2(mT3Min + mT4Min + mT5Min) / s);
  }

  // Check that there is an open range.
  return (tauMax > tauMin);
}

//--------------------------------------------------------------------------

// Find range of allowed y values.

bool PhaseSpace::limitY() {

  // Trivial reply for unresolved lepton beams.
  if (hasTwoPointParticles) {
    yMax = 1.;
    return true;
  }

  // Requirements from selected tau value. Trivial for one unresolved beam.
  yMax = -0.5 * log(tau);
  if (hasOnePointParticle) return true;

  // For lepton beams requirements from cutoff for f_e^e.
  double yMaxMargin = (hasTwoLeptonBeams) ? yMax + LEPTONXLOGMAX : yMax;

  // Check that there is an open range.
  return (yMaxMargin > 0.);
}

//--------------------------------------------------------------------------

// Find range of allowed z = cos(theta) values.

bool PhaseSpace::limitZ() {

  // Default limits.
  zMin = 0.;
  zMax = 1.;

  // Requirements from pTHat limits.
  zMax = sqrtpos( 1. - pT2HatMin / p2Abs );
  if (pTHatMax > pTHatMin) zMin = sqrtpos( 1. - pT2HatMax / p2Abs );

  // Check that there is an open range so far.
  hasNegZ = false;
  hasPosZ = false;
  if (zMax < zMin) return false;

  // Define two individual ranges.
  hasNegZ = true;
  hasPosZ = true;
  zNegMin = -zMax;
  zNegMax = -zMin;
  zPosMin =  zMin;
  zPosMax =  zMax;

  // Optionally introduce Q2 = -tHat cut.
  if (hasQ2Min) {
    double zMaxQ2 = (sH - s3 - s4 - 2. * Q2GlobalMin) / (2. * pAbs * mHat);
    if (zMaxQ2 > zPosMin) {
      if (zMaxQ2 < zPosMax) zPosMax = zMaxQ2;
    } else {
      hasPosZ = false;
      zPosMax = zPosMin;
      if (zMaxQ2 > zNegMin) {
        if (zMaxQ2 < zNegMax) zNegMax = zMaxQ2;
      } else {
        hasNegZ = false;
        zNegMin = zNegMax;
      }
    }
  }

  // Check that there is an open range.
  return hasNegZ;
}

//--------------------------------------------------------------------------

// Select tau according to a choice of shapes.

void PhaseSpace::selectTau(int iTau, double tauVal, bool is2) {

  // Trivial reply for unresolved lepton beams.
  if (hasTwoPointParticles) {
    tau = 1.;
    wtTau = 1.;
    sH = s;
    mHat = sqrt(sH);
    if (is2) {
      p2Abs = 0.25 * (pow2(sH - s3 - s4) - 4. * s3 * s4) / sH;
      pAbs = sqrtpos( p2Abs );
    }
    return;
  }

  // Contributions from s-channel resonances.
  double tRatA = 0.;
  double aLowA = 0.;
  double aUppA = 0.;
  if (idResA !=0) {
    tRatA = ((tauResA + tauMax) / (tauResA + tauMin)) * (tauMin / tauMax);
    aLowA = atan( (tauMin - tauResA) / widResA);
    aUppA = atan( (tauMax - tauResA) / widResA);
  }
  double tRatB = 0.;
  double aLowB = 0.;
  double aUppB = 0.;
  if (idResB != 0) {
    tRatB = ((tauResB + tauMax) / (tauResB + tauMin)) * (tauMin / tauMax);
    aLowB = atan( (tauMin - tauResB) / widResB);
    aUppB = atan( (tauMax - tauResB) / widResB);
  }

  // Contributions from 1 / (1 - tau)  for lepton beams.
  double aLowT = 0.;
  double aUppT = 0.;
  if (hasTwoLeptonBeams) {
    aLowT = log( max( LEPTONTAUMIN, 1. - tauMin) );
    aUppT = log( max( LEPTONTAUMIN, 1. - tauMax) );
    intTau6 = aLowT - aUppT;
  }

  // Select according to 1/tau or 1/tau^2.
  if (iTau == 0) tau = tauMin * pow( tauMax / tauMin, tauVal);
  else if (iTau == 1) tau = tauMax * tauMin
    / (tauMin + (tauMax - tauMin) * tauVal);

  // Select according to 1 / (1 - tau) for lepton beams.
  else if (hasTwoLeptonBeams && iTau == nTau - 1)
    tau = 1. - exp( aUppT + intTau6 * tauVal );

  // Select according to 1 / (tau * (tau + tauRes)) or
  // 1 / ((tau - tauRes)^2 + widRes^2) for resonances A and B.
  else if (iTau == 2) tau = tauResA * tauMin
    / ((tauResA + tauMin) * pow( tRatA, tauVal) - tauMin);
  else if (iTau == 3) tau = tauResA + widResA
    * tan( aLowA + (aUppA - aLowA) * tauVal);
  else if (iTau == 4) tau = tauResB * tauMin
    / ((tauResB + tauMin) * pow( tRatB, tauVal) - tauMin);
  else if (iTau == 5) tau = tauResB + widResB
    * tan( aLowB + (aUppB - aLowB) * tauVal);

  // Phase space weight in tau.
  intTau0 = log( tauMax / tauMin);
  intTau1 = (tauMax - tauMin) / (tauMax * tauMin);
  double invWtTau = (tauCoef[0] / intTau0) + (tauCoef[1] / intTau1) / tau;
  if (idResA != 0) {
    intTau2 = -log(tRatA) / tauResA;
    intTau3 = (aUppA - aLowA) / widResA;
    invWtTau += (tauCoef[2] / intTau2) / (tau + tauResA)
      + (tauCoef[3] / intTau3) * tau / ( pow2(tau - tauResA) + pow2(widResA) );
  }
  if (idResB != 0) {
    intTau4 = -log(tRatB) / tauResB;
    intTau5 = (aUppB - aLowB) / widResB;
    invWtTau += (tauCoef[4] / intTau4) / (tau + tauResB)
      + (tauCoef[5] / intTau5) * tau / ( pow2(tau - tauResB) + pow2(widResB) );
  }
  if (hasTwoLeptonBeams)
    invWtTau += (tauCoef[nTau - 1] / intTau6)
      * tau / max( LEPTONTAUMIN, 1. - tau);
  wtTau = 1. / invWtTau;

  // Calculate sHat and absolute momentum of outgoing partons.
  sH = tau * s;
  mHat = sqrt(sH);
  if (is2) {
    p2Abs = 0.25 * (pow2(sH - s3 - s4) - 4. * s3 * s4) / sH;
    pAbs = sqrtpos( p2Abs );
  }

}

//--------------------------------------------------------------------------

// Select y according to a choice of shapes.

void PhaseSpace::selectY(int iY, double yVal) {

  // Trivial reply for two unresolved lepton beams.
  if (hasTwoPointParticles) {
    y = 0.;
    wtY = 1.;
    x1H = 1.;
    x2H = 1.;
    return;
  }

  // Trivial replies for one unresolved lepton beam.
  if (hasOnePointParticle) {
    if (hasLeptonBeamA || hasPointGammaA) {
      y   = yMax;
      x1H = 1.;
      x2H = tau;
    } else {
      y   = -yMax;
      x1H = tau;
      x2H = 1.;
    }
    wtY = 1.;
    return;
  }

  // For lepton beams skip options 3&4 and go straight to 5&6.
  if (hasTwoLeptonBeams && iY > 2) iY += 2;

  // Standard expressions used below.
  double expYMax = exp( yMax );
  double expYMin = exp(-yMax );
  double atanMax = atan( expYMax );
  double atanMin = atan( expYMin );
  double aUppY = (hasTwoLeptonBeams)
    ? log( max( LEPTONXMIN, LEPTONXMAX / tau - 1. ) ) : 0.;
  double aLowY = LEPTONXLOGMIN;

  // 1 / cosh(y).
  if (iY == 0) y = log( tan( atanMin + (atanMax - atanMin) * yVal ) );

  // y - y_min or mirrored y_max - y.
  else if (iY <= 2) y = yMax * (2. * sqrt(yVal) - 1.);

  // exp(y) or mirrored exp(-y).
  else if (iY <= 4) y = log( expYMin + (expYMax - expYMin) * yVal );

  // 1 / (1 - exp(y - y_max)) or mirrored 1 / (1 - exp(y_min - y)).
  else y = yMax - log( 1. + exp(aLowY + (aUppY - aLowY) * yVal) );

  // Mirror two cases.
  if (iY == 2 || iY == 4 || iY == 6) y = -y;

  // Phase space integral in y.
  intY0  = 2. * (atanMax - atanMin);
  intY12 = 0.5 * pow2(2. * yMax);
  intY34 = expYMax - expYMin;
  intY56 = aUppY - aLowY;
  double invWtY = (yCoef[0] / intY0) / cosh(y)
     + (yCoef[1] / intY12) * (y + yMax) + (yCoef[2] / intY12) * (yMax - y);
  if (!hasTwoLeptonBeams) invWtY
    += (yCoef[3] / intY34) * exp(y)     + (yCoef[4] / intY34) * exp(-y);
  else invWtY
    += (yCoef[3] / intY56) / max( LEPTONXMIN, 1. - exp( y - yMax) )
    +  (yCoef[4] / intY56) / max( LEPTONXMIN, 1. - exp(-y - yMax) );
  wtY = 1. / invWtY;

  // Calculate x1 and x2.
  x1H = sqrt(tau) * exp(y);
  x2H = sqrt(tau) * exp(-y);
}

//--------------------------------------------------------------------------

// Select z = cos(theta) according to a choice of shapes.
// The selection is split in the positive- and negative-z regions,
// since a pTmax cut can remove the region around z = 0.
// Furthermore, a Q2 (= -tHat) cut can make the two regions asymmetric.

void PhaseSpace::selectZ(int iZ, double zVal) {

  // Mass-dependent dampening of pT -> 0 limit.
  ratio34 = max(TINY, 2. * s3 * s4 / pow2(sH));
  unity34 = 1. + ratio34;
  double ratiopT2 = 2. * pT2HatMin / max( SHATMINZ, sH);
  if (ratiopT2 < PT2RATMINZ) ratio34 = max( ratio34, ratiopT2);

  // Common expressions of unity - z and unity + z limits, protected from 0.
  double zNegMinM = max(ratio34, unity34 - zNegMin);
  double zNegMaxM = max(ratio34, unity34 - zNegMax);
  double zPosMinM = max(ratio34, unity34 - zPosMin);
  double zPosMaxM = max(ratio34, unity34 - zPosMax);
  double zNegMinP = max(ratio34, unity34 + zNegMin);
  double zNegMaxP = max(ratio34, unity34 + zNegMax);
  double zPosMinP = max(ratio34, unity34 + zPosMin);
  double zPosMaxP = max(ratio34, unity34 + zPosMax);

  // Evaluate integrals over negative and positive z ranges.
  // Flat in z.
  double area0Neg = zNegMax - zNegMin;
  double area0Pos = zPosMax - zPosMin;
  double area0    = area0Neg + area0Pos;
  // 1 / (unity34 - z).
  double area1Neg = log(zNegMinM / zNegMaxM);
  double area1Pos = log(zPosMinM / zPosMaxM);
  double area1    = area1Neg + area1Pos;
  // 1 / (unity34 + z).
  double area2Neg = log(zNegMaxP / zNegMinP);
  double area2Pos = log(zPosMaxP / zPosMinP);
  double area2    = area2Neg + area2Pos;
  // 1 / (unity34 - z)^2.
  double area3Neg = 1. / zNegMaxM - 1. / zNegMinM;
  double area3Pos = 1. / zPosMaxM - 1. / zPosMinM;
  double area3    = area3Neg + area3Pos;
  // 1 / (unity34 + z)^2.
  double area4Neg = 1. / zNegMinP - 1. / zNegMaxP;
  double area4Pos = 1. / zPosMinP - 1. / zPosMaxP;
  double area4    = area4Neg + area4Pos;

  // Pick z value according to alternatives.
  // Flat in z.
  if (iZ == 0) {
    if (!hasPosZ || zVal * area0 < area0Neg) {
      double zValMod = zVal * area0 / area0Neg;
      z = zNegMin + zValMod * area0Neg;
    } else {
      double zValMod = (zVal * area0 - area0Neg) / area0Pos;
      z = zPosMin + zValMod * area0Pos;
    }

  // 1 / (unity34 - z).
  } else if (iZ == 1) {
    if (!hasPosZ || zVal * area1 < area1Neg) {
      double zValMod = zVal * area1 / area1Neg;
      z = unity34 - zNegMinM * pow(zNegMaxM / zNegMinM, zValMod);
    } else {
      double zValMod = (zVal * area1 - area1Neg)/ area1Pos;
      z = unity34 - zPosMinM * pow(zPosMaxM / zPosMinM, zValMod);
    }

  // 1 / (unity34 + z).
  } else if (iZ == 2) {
    if (!hasPosZ || zVal * area2 < area2Neg) {
      double zValMod = zVal * area2 / area2Neg;
      z = zNegMinP * pow(zNegMaxP / zNegMinP, zValMod) - unity34;
    } else {
      double zValMod = (zVal * area2 - area2Neg)/ area2Pos;
      z = zPosMinP * pow(zPosMaxP / zPosMinP, zValMod) - unity34;
    }

  // 1 / (unity34 - z)^2.
  } else if (iZ == 3) {
    if (!hasPosZ || zVal * area3 < area3Neg) {
      double zValMod = zVal * area3 / area3Neg;
      z = unity34 - 1. / (1./zNegMinM + area3Neg * zValMod);
    } else {
      double zValMod = (zVal * area3 - area3Neg)/ area3Pos;
      z = unity34 - 1. / (1./zPosMinM + area3Pos * zValMod);
    }

  // 1 / (unity34 + z)^2.
  } else if (iZ == 4) {
    if (!hasPosZ || zVal * area4 < area4Neg) {
      double zValMod = zVal * area4 / area4Neg;
      z = 1. / (1./zNegMinP - area4Neg * zValMod) - unity34;
    } else {
      double zValMod = (zVal * area4 - area4Neg)/ area4Pos;
      z = 1. / (1./zPosMinP - area4Pos * zValMod) - unity34;
    }
  }

  // Safety check for roundoff errors. Combinations with z.
  if (z < 0.) z = min( zNegMax, max( zNegMin, z));
  else        z = min( zPosMax, max( zPosMin, z) );
  zNeg = max(ratio34, unity34 - z);
  zPos = max(ratio34, unity34 + z);

  // Phase space integral in z.
  wtZ = mHat * pAbs / ( (zCoef[0] / area0) + (zCoef[1] / area1) / zNeg
    + (zCoef[2] / area2) / zPos + (zCoef[3] / area3) / pow2(zNeg)
    + (zCoef[4] / area4) / pow2(zPos) );

  // Calculate tHat and uHat. Also gives pTHat.
  double sH34 = -0.5 * (sH - s3 - s4);
  double tHuH = pow2(sH34) * (1. - z) * (1. + z) + s3 * s4 * pow2(z);
  if (z < 0.) {
    tH = sH34 + mHat * pAbs * z;
    uH = tHuH / tH;
  } else {
    uH = sH34 - mHat * pAbs * z;
    tH = tHuH / uH;
  }
  pTH = sqrtpos( (tH * uH - s3 * s4) / sH);

}

//--------------------------------------------------------------------------

// Select three-body phase space according to a cylindrically based form
// that can be chosen to favour low pT based on the form of propagators.

bool PhaseSpace::select3Body() {

  // Upper and lower limits of pT choice for 4 and 5.
  double m35S = pow2(m3 + m5);
  double pT4Smax = 0.25 * ( pow2(sH - s4 - m35S) - 4. * s4 * m35S ) / sH;
  if (pTHatMax > pTHatMin) pT4Smax = min( pT2HatMax, pT4Smax);
  double pT4Smin = pT2HatMin;
  double m34S = pow2(m3 + m4);
  double pT5Smax = 0.25 * ( pow2(sH - s5 - m34S) - 4. * s5 * m34S ) / sH;
  if (pTHatMax > pTHatMin) pT5Smax = min( pT2HatMax, pT5Smax);
  double pT5Smin = pT2HatMin;

  // Check that pT ranges not closed.
  if ( pT4Smax < pow2(pTHatMin + MASSMARGIN) ) return false;
  if ( pT5Smax < pow2(pTHatMin + MASSMARGIN) ) return false;

  // Select pT4S according to c0 + c1/(M^2 + pT^2) + c2/(M^2 + pT^2)^2.
  double pTSmaxProp = pT4Smax + sTchan1;
  double pTSminProp = pT4Smin + sTchan1;
  double pTSratProp = pTSmaxProp / pTSminProp;
  double pTSdiff    = pT4Smax - pT4Smin;
  double rShape     = rndmPtr->flat();
  double pT4S       = 0.;
  if (rShape < frac3Flat) pT4S = pT4Smin + rndmPtr->flat() * pTSdiff;
  else if (rShape < frac3Flat + frac3Pow1) pT4S = max( pT2HatMin,
    pTSminProp * pow( pTSratProp, rndmPtr->flat() ) - sTchan1 );
  else pT4S = max( pT2HatMin, pTSminProp * pTSmaxProp
    / (pTSminProp + rndmPtr->flat()* pTSdiff) - sTchan1 );
  double wt4 = pTSdiff / ( frac3Flat
    + frac3Pow1 * pTSdiff / (log(pTSratProp) * (pT4S + sTchan1))
    + frac3Pow2 * pTSminProp * pTSmaxProp / pow2(pT4S + sTchan1) );

  // Select pT5S according to c0 + c1/(M^2 + pT^2) + c2/(M^2 + pT^2)^2.
  pTSmaxProp  = pT5Smax + sTchan2;
  pTSminProp  = pT5Smin + sTchan2;
  pTSratProp  = pTSmaxProp / pTSminProp;
  pTSdiff     = pT5Smax - pT5Smin;
  rShape      = rndmPtr->flat();
  double pT5S = 0.;
  if (rShape < frac3Flat) pT5S = pT5Smin + rndmPtr->flat() * pTSdiff;
  else if (rShape < frac3Flat + frac3Pow1) pT5S = max( pT2HatMin,
    pTSminProp * pow( pTSratProp, rndmPtr->flat() ) - sTchan2 );
  else pT5S = max( pT2HatMin, pTSminProp * pTSmaxProp
    / (pTSminProp + rndmPtr->flat()* pTSdiff) - sTchan2 );
  double wt5 = pTSdiff / ( frac3Flat
    + frac3Pow1 * pTSdiff / (log(pTSratProp) * (pT5S + sTchan2))
    + frac3Pow2 * pTSminProp * pTSmaxProp / pow2(pT5S + sTchan2) );

  // Select azimuthal angles and check that third pT in range.
  double phi4 = 2. * M_PI * rndmPtr->flat();
  double phi5 = 2. * M_PI * rndmPtr->flat();
  double pT3S = max( 0., pT4S + pT5S + 2. * sqrt(pT4S * pT5S)
    * cos(phi4 - phi5) );
  if ( pT3S < pT2HatMin || (pTHatMax > pTHatMin && pT3S > pT2HatMax) )
    return false;

  // Calculate transverse masses and check that phase space not closed.
  double sT3 = s3 + pT3S;
  double sT4 = s4 + pT4S;
  double sT5 = s5 + pT5S;
  double mT3 = sqrt(sT3);
  double mT4 = sqrt(sT4);
  double mT5 = sqrt(sT5);
  if ( mT3 + mT4 + mT5 + MASSMARGIN > mHat ) return false;

  // Select rapidity for particle 3 and check that phase space not closed.
  double m45S = pow2(mT4 + mT5);
  double y3max = log( ( sH + sT3 - m45S + sqrtpos( pow2(sH - sT3 - m45S)
    - 4 * sT3 * m45S ) ) / (2. * mHat * mT3) );
  if (y3max < YRANGEMARGIN) return false;
  double y3    = (2. * rndmPtr->flat() - 1.) * (1. -  YRANGEMARGIN) * y3max;
  double pz3   = mT3 * sinh(y3);
  double e3    = mT3 * cosh(y3);

  // Find momentum transfers in the two mirror solutions (in 4-5 frame).
  double pz45  = -pz3;
  double e45   = mHat - e3;
  double sT45  = e45 * e45 - pz45 * pz45;
  double lam45 = sqrtpos( pow2(sT45 - sT4 - sT5) - 4. * sT4 * sT5 );
  if (lam45 < YRANGEMARGIN * sH) return false;
  double lam4e = sT45 + sT4 - sT5;
  double lam5e = sT45 + sT5 - sT4;
  double tFac  = -0.5 * mHat / sT45;
  double t1Pos = tFac * (e45 - pz45) * (lam4e - lam45);
  double t1Neg = tFac * (e45 - pz45) * (lam4e + lam45);
  double t2Pos = tFac * (e45 + pz45) * (lam5e - lam45);
  double t2Neg = tFac * (e45 + pz45) * (lam5e + lam45);

  // Construct relative mirror weights and make choice.
  double wtPosUnnorm = 1.;
  double wtNegUnnorm = 1.;
  if (useMirrorWeight) {
    wtPosUnnorm  = 1./ pow2( (t1Pos - sTchan1) * (t2Pos - sTchan2) );
    wtNegUnnorm  = 1./ pow2( (t1Neg - sTchan1) * (t2Neg - sTchan2) );
  }
  double wtPos   = wtPosUnnorm / (wtPosUnnorm + wtNegUnnorm);
  double wtNeg   = wtNegUnnorm / (wtPosUnnorm + wtNegUnnorm);
  double epsilon = (rndmPtr->flat() < wtPos) ? 1. : -1.;

  // Construct four-vectors in rest frame of subprocess.
  double px4 = sqrt(pT4S) * cos(phi4);
  double py4 = sqrt(pT4S) * sin(phi4);
  double px5 = sqrt(pT5S) * cos(phi5);
  double py5 = sqrt(pT5S) * sin(phi5);
  double pz4 = 0.5 * (pz45 * lam4e + epsilon * e45 * lam45) / sT45;
  double pz5 = pz45 - pz4;
  double e4  = sqrt(sT4 + pz4 * pz4);
  double e5  = sqrt(sT5 + pz5 * pz5);
  p3cm       = Vec4( -(px4 + px5), -(py4 + py5), pz3, e3);
  p4cm       = Vec4( px4, py4, pz4, e4);
  p5cm       = Vec4( px5, py5, pz5, e5);

  // Total weight to associate with kinematics choice.
  wt3Body    = wt4 * wt5 * (2. * y3max) / (128. * pow3(M_PI) * lam45);
  wt3Body   *= (epsilon > 0.) ? 1. / wtPos : 1. / wtNeg;

  // Cross section of form |M|^2/(2 sHat) dPS_3 so need 1/(2 sHat).
  wt3Body   /= (2. * sH);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Solve linear equation system for better phase space coefficients.

void PhaseSpace::solveSys( int n, int bin[8], double vec[8],
  double mat[8][8], double coef[8]) {

  // Optional printout.
  if (showSearch) {
    cout << "\n Equation system: " << setw(5) << bin[0];
    for (int j = 0; j < n; ++j) cout << setw(12) << mat[0][j];
    cout << setw(12) << vec[0] << "\n";
    for (int i = 1; i < n; ++i) {
      cout << "                  " << setw(5) << bin[i];
      for (int j = 0; j < n; ++j) cout << setw(12) << mat[i][j];
      cout << setw(12) << vec[i] << "\n";
    }
  }

  // Local variables.
  double vecNor[8], coefTmp[8];
  for (int i = 0; i < n; ++i) coefTmp[i] = 0.;

  // Check if equation system solvable.
  bool canSolve = true;
  for (int i = 0; i < n; ++i) if (bin[i] == 0) canSolve = false;
  double vecSum = 0.;
  for (int i = 0; i < n; ++i) vecSum += vec[i];
  if (abs(vecSum) < TINY) canSolve = false;

  // Solve to find relative importance of cross-section pieces.
  if (canSolve) {
    for (int i = 0; i < n; ++i) vecNor[i] = max( 0.1, vec[i] / vecSum);
    for (int k = 0; k < n - 1; ++k) {
      for (int i = k + 1; i < n; ++i) {
        if (abs(mat[k][k]) < TINY) {canSolve = false; break;}
        double ratio = mat[i][k] / mat[k][k];
        vec[i] -= ratio * vec[k];
        for (int j = k; j < n; ++j) mat[i][j] -= ratio * mat[k][j];
      }
      if (!canSolve) break;
    }
    if (canSolve) {
      for (int k = n - 1; k >= 0; --k) {
        for (int j = k + 1; j < n; ++j) vec[k] -= mat[k][j] * coefTmp[j];
        coefTmp[k] = vec[k] / mat[k][k];
      }
    }
  }

  // Share evenly if failure.
  if (!canSolve) for (int i = 0; i < n; ++i) {
    coefTmp[i] = 1.;
    vecNor[i] = 0.1;
    if (vecSum > TINY) vecNor[i] = max(0.1, vec[i] / vecSum);
  }

  // Normalize coefficients, with piece shared democratically.
  double coefSum = 0.;
  vecSum = 0.;
  for (int i = 0; i < n; ++i) {
    coefTmp[i] = max( 0., coefTmp[i]);
    coefSum += coefTmp[i];
    vecSum += vecNor[i];
  }
  if (coefSum > 0.) for (int i = 0; i < n; ++i) coef[i] = EVENFRAC / n
    + (1. - EVENFRAC) * 0.5 * (coefTmp[i] / coefSum + vecNor[i] / vecSum);
  else for (int i = 0; i < n; ++i) coef[i] = 1. / n;

  // Optional printout.
  if (showSearch) {
    cout << " Solution:             ";
    for (int i = 0; i < n; ++i) cout << setw(12) << coef[i];
    cout << "\n";
  }
}

//--------------------------------------------------------------------------

// Setup mass selection for one resonance at a time - part 1.

void PhaseSpace::setupMass1(int iM) {

  // Identity for mass seletion; is 0 also for light quarks (not yet selected).
  if (iM == 3) idMass[iM] = abs(sigmaProcessPtr->id3Mass());
  if (iM == 4) idMass[iM] = abs(sigmaProcessPtr->id4Mass());
  if (iM == 5) idMass[iM] = abs(sigmaProcessPtr->id5Mass());

  // Masses and widths of resonances.
  if (idMass[iM] == 0) {
    mPeak[iM]  = 0.;
    mWidth[iM] = 0.;
    mMin[iM]   = 0.;
    mMax[iM]   = 0.;
  } else {
    mPeak[iM]  = particleDataPtr->m0(idMass[iM]);
    mWidth[iM] = particleDataPtr->mWidth(idMass[iM]);
    mMin[iM]   = max( MRESMINABS, particleDataPtr->mMin(idMass[iM]) );
    mMax[iM]   = particleDataPtr->mMax(idMass[iM]);
    // gmZmode == 1 means pure photon propagator; set at lower mass limit.
    if (idMass[iM] == 23 && gmZmode == 1) mPeak[iM] = mMin[iM];
  }

  // Mass and width combinations for Breit-Wigners.
  sPeak[iM]       = mPeak[iM] * mPeak[iM];
  useBW[iM]       = useBreitWigners && (mWidth[iM] > minWidthBreitWigners);
  useNarrowBW[iM] = useBreitWigners && !useBW[iM]
                  && (mWidth[iM] > minWidthNarrowBW);
  if (!useBW[iM] && !useNarrowBW[iM]) mWidth[iM] = 0.;
  mw[iM]          = mPeak[iM] * mWidth[iM];
  wmRat[iM]       = (idMass[iM] == 0 || mPeak[iM] == 0.)
                  ? 0. : mWidth[iM] / mPeak[iM];

  // Simple Breit-Wigner range, upper edge to be corrected subsequently.
  if (useBW[iM]) {
    mLower[iM] = mMin[iM];
    mUpper[iM] = mHatMax;
  }

}

//--------------------------------------------------------------------------

// Setup mass selection for one resonance at a time - part 2.

void PhaseSpace::setupMass2(int iM, double distToThresh) {

  // Store reduced Breit-Wigner range.
  if (mMax[iM] > mMin[iM]) mUpper[iM] = min( mUpper[iM], mMax[iM]);
  sLower[iM]     = mLower[iM] * mLower[iM];
  sUpper[iM]     = mUpper[iM] * mUpper[iM];

  // Prepare to select m3 by BW + flat + 1/s_3.
  // Determine relative coefficients by allowed mass range.
  if (distToThresh > THRESHOLDSIZE) {
    fracFlatS[iM] = 0.1;
    fracFlatM[iM] = 0.1;
    fracInv[iM]   = 0.1;
  } else if (distToThresh > - THRESHOLDSIZE) {
    fracFlatS[iM] = 0.25 - 0.15 * distToThresh / THRESHOLDSIZE;
    fracInv [iM]  = 0.15 - 0.05 * distToThresh / THRESHOLDSIZE;
  } else {
   fracFlatS[iM]  = 0.3;
   fracFlatM[iM]  = 0.1;
   fracInv[iM]    = 0.2;
  }

  // For gamma*/Z0: increase 1/s_i part and introduce 1/s_i^2 part.
  fracInv2[iM]   = 0.;
  if (idMass[iM] == 23 && gmZmode == 0) {
    fracFlatS[iM] *= 0.5;
    fracFlatM[iM] *= 0.5;
    fracInv[iM]    = 0.5 * fracInv[iM] + 0.25;
    fracInv2[iM]   = 0.25;
  } else if (idMass[iM] == 23 && gmZmode == 1) {
    fracFlatS[iM] = 0.1;
    fracFlatM[iM] = 0.1;
    fracInv[iM]   = 0.35;
    fracInv2[iM]  = 0.35;
  }

  // Normalization integrals for the respective contribution.
  atanLower[iM]  = atan( (sLower[iM] - sPeak[iM])/ mw[iM] );
  atanUpper[iM]  = atan( (sUpper[iM] - sPeak[iM])/ mw[iM] );
  intBW[iM]      = atanUpper[iM] - atanLower[iM];
  intFlatS[iM]    = sUpper[iM] - sLower[iM];
  intFlatM[iM]    = mUpper[iM] - mLower[iM];
  intInv[iM]     = log( sUpper[iM] / sLower[iM] );
  intInv2[iM]    = 1./sLower[iM] - 1./sUpper[iM];

}

//--------------------------------------------------------------------------

// Select Breit-Wigner-distributed or fixed masses.

void PhaseSpace::trialMass(int iM) {

  // References to masses to be set.
  double& mSet = (iM == 3) ? m3 : ( (iM == 4) ? m4 : m5 );
  double& sSet = (iM == 3) ? s3 : ( (iM == 4) ? s4 : s5 );

  // Distribution for m_i is BW + flat(s) + 1/sqrt(s_i) + 1/s_i + 1/s_i^2.
  if (useBW[iM]) {
    double pickForm = rndmPtr->flat();
    if (pickForm > fracFlatS[iM] + fracFlatM[iM] + fracInv[iM] + fracInv2[iM])
      sSet = sPeak[iM] + mw[iM] * tan( atanLower[iM]
           + rndmPtr->flat() * intBW[iM] );
    else if (pickForm > fracFlatM[iM] + fracInv[iM] + fracInv2[iM])
      sSet = sLower[iM] + rndmPtr->flat() * (sUpper[iM] - sLower[iM]);
    else if (pickForm > fracInv[iM] + fracInv2[iM])
      sSet = pow2(mLower[iM] + rndmPtr->flat() * (mUpper[iM] - mLower[iM]));
    else if (pickForm > fracInv2[iM])
      sSet = sLower[iM] * pow( sUpper[iM] / sLower[iM], rndmPtr->flat() );
    else sSet = sLower[iM] * sUpper[iM]
      / (sLower[iM] + rndmPtr->flat() * (sUpper[iM] - sLower[iM]));
    mSet = sqrt(sSet);

  // Distribution for m_i is simple BW.
  } else if (useNarrowBW[iM]) {
    mSet = particleDataPtr->mSel(idMass[iM]);
    sSet = mSet * mSet;

  // Else m_i is fixed at peak value.
  } else {
    mSet = mPeak[iM];
    sSet = sPeak[iM];
  }

}

//--------------------------------------------------------------------------

// Naively a fixed-width Breit-Wigner is used to pick the mass.
// Here come the correction factors for
// (i) preselection from BW + flat in s_i + 1/sqrt(s_i) + 1/s_i  + 1/s_i^2,
// (ii) reduced allowed mass range,
// (iii) running width, i.e. m0*Gamma0 -> s*Gamma0/m0.
// In the end, the weighted distribution is a running-width BW.

double PhaseSpace::weightMass(int iM) {

  // Reference to mass and to Breit-Wigner weight to be set.
  double& mSet   = (iM == 3) ? m3 : ( (iM == 4) ? m4 : m5 );
  double& sSet   = (iM == 3) ? s3 : ( (iM == 4) ? s4 : s5 );
  double& runBWH = (iM == 3) ? runBW3H : ( (iM == 4) ? runBW4H : runBW5H );

  // Default weight if no Breit-Wigner.
  runBWH = 1.;
  if (!useBW[iM]) return 1.;

  // Weight of generated distribution.
  double genBW
    = (1. - fracFlatS[iM] - fracFlatM[iM] - fracInv[iM] - fracInv2[iM])
      * mw[iM] / ( (pow2(sSet - sPeak[iM]) + pow2(mw[iM])) * intBW[iM])
    + fracFlatS[iM] / intFlatS[iM]
    + fracFlatM[iM] / (2. * mSet * intFlatM[iM])
    + fracInv[iM] / (sSet * intInv[iM])
    + fracInv2[iM] / (sSet*sSet * intInv2[iM]);

  // Weight of distribution with running width in Breit-Wigner.
  double mwRun = sSet * wmRat[iM];
  //?? Alternative recipe, taking into account that decay channels close
  // at different mass thresholds Needs refining, e.g. no doublecouting
  // with openFrac and difference numerator/denominator.
  //double mwRun = mSet * particleDataPtr->resWidthOpen(idMass[iM], mSet);
  runBWH = mwRun / (pow2(sSet - sPeak[iM]) + pow2(mwRun)) / M_PI;

  // Done.
  return (runBWH / genBW);

}

//==========================================================================

// PhaseSpace2to1tauy class.
// 2 -> 1 kinematics for normal subprocesses.

//--------------------------------------------------------------------------

// Set limits for resonance mass selection.

bool PhaseSpace2to1tauy::setupMass() {

  // Treat Z0 as such or as gamma*/Z0
  gmZmode         = gmZmodeGlobal;
  int gmZmodeProc = sigmaProcessPtr->gmZmode();
  if (gmZmodeProc >= 0) gmZmode = gmZmodeProc;

  // Mass limits for current resonance.
  int idRes = abs(sigmaProcessPtr->resonanceA());
  int idTmp = abs(sigmaProcessPtr->resonanceB());
  if (idTmp > 0) idRes = idTmp;
  double mResMin = (idRes == 0) ? 0. : particleDataPtr->mMin(idRes);
  double mResMax = (idRes == 0) ? 0. : particleDataPtr->mMax(idRes);

  // Compare with global mass limits and pick tighter of them.
  mHatMin = max( mResMin, mHatGlobalMin);
  sHatMin = mHatMin*mHatMin;
  mHatMax = eCM;
  if (mResMax > mResMin) mHatMax = min( mHatMax, mResMax);
  if (mHatGlobalMax > mHatGlobalMin) mHatMax = min( mHatMax, mHatGlobalMax);
  sHatMax = mHatMax*mHatMax;

  // Default Breit-Wigner weight.
  wtBW = 1.;

  // Fail if mass window (almost) closed.
  return (mHatMax > mHatMin + MASSMARGIN);

}

//--------------------------------------------------------------------------

// Construct the four-vector kinematics from the trial values.

bool PhaseSpace2to1tauy::finalKin() {

  // Particle masses; incoming always on mass shell.
  mH[1] = 0.;
  mH[2] = 0.;
  mH[3] = mHat;

  // Incoming partons along beam axes. Outgoing has sum of momenta.
  pH[1] = Vec4( 0., 0.,  0.5 * eCM * x1H, 0.5 * eCM * x1H);
  pH[2] = Vec4( 0., 0., -0.5 * eCM * x2H, 0.5 * eCM * x2H);
  pH[3] = pH[1] + pH[2];

  // Done.
  return true;
}

//==========================================================================

// PhaseSpace2to2tauyz class.
// 2 -> 2 kinematics for normal subprocesses.

//--------------------------------------------------------------------------

// Set up for fixed or Breit-Wigner mass selection.

bool PhaseSpace2to2tauyz::setupMasses() {

  // Treat Z0 as such or as gamma*/Z0
  gmZmode         = gmZmodeGlobal;
  int gmZmodeProc = sigmaProcessPtr->gmZmode();
  if (gmZmodeProc >= 0) gmZmode = gmZmodeProc;

  // Set sHat limits - based on global limits only.
  mHatMin = mHatGlobalMin;
  sHatMin = mHatMin*mHatMin;
  mHatMax = eCM;
  if (mHatGlobalMax > mHatGlobalMin) mHatMax = min( eCM, mHatGlobalMax);
  sHatMax = mHatMax*mHatMax;

  // Masses and widths of resonances.
  setupMass1(3);
  setupMass1(4);

  // Reduced mass range when two massive particles.
  if (useBW[3]) mUpper[3] -= (useBW[4]) ? mMin[4] : mPeak[4];
  if (useBW[4]) mUpper[4] -= (useBW[3]) ? mMin[3] : mPeak[3];

  // If closed phase space then unallowed process.
  bool physical = true;
  if (useBW[3] && mUpper[3] < mLower[3] + MASSMARGIN) physical = false;
  if (useBW[4] && mUpper[4] < mLower[4] + MASSMARGIN) physical = false;
  if (!useBW[3] && !useBW[4] && mHatMax < mPeak[3] + mPeak[4] + MASSMARGIN)
    physical = false;
  if (!physical) return false;

  // If either particle is massless then need extra pTHat cut.
  pTHatMin   = pTHatGlobalMin;
  if (mPeak[3] < pTHatMinDiverge || mPeak[4] < pTHatMinDiverge)
    pTHatMin = max( pTHatMin, pTHatMinDiverge);
  pT2HatMin  = pTHatMin * pTHatMin;
  pTHatMax   = pTHatGlobalMax;
  pT2HatMax  = pTHatMax * pTHatMax;

  // Prepare to select m3 by BW + flat + 1/s_3.
  if (useBW[3]) {
    double distToThreshA = (mHatMax - mPeak[3] - mPeak[4]) * mWidth[3]
      / (pow2(mWidth[3]) + pow2(mWidth[4]));
    double distToThreshB = (mHatMax - mPeak[3] - mMin[4]) / mWidth[3];
    double distToThresh = min( distToThreshA, distToThreshB);
    setupMass2(3, distToThresh);
  }

  // Prepare to select m4 by BW + flat + 1/s_4.
  if (useBW[4]) {
    double distToThreshA = (mHatMax - mPeak[3] - mPeak[4]) * mWidth[4]
      / (pow2(mWidth[3]) + pow2(mWidth[4]));
    double distToThreshB = (mHatMax - mMin[3] - mPeak[4]) / mWidth[4];
    double distToThresh = min( distToThreshA, distToThreshB);
    setupMass2(4, distToThresh);
  }

  // Initialization masses. Special cases when constrained phase space.
  m3 = (useBW[3]) ? min(mPeak[3], mUpper[3]) : mPeak[3];
  m4 = (useBW[4]) ? min(mPeak[4], mUpper[4]) : mPeak[4];
  if (m3 + m4 + THRESHOLDSIZE * (mWidth[3] + mWidth[4]) + MASSMARGIN
    > mHatMax) {
    if (useBW[3] && useBW[4]) physical = constrainedM3M4();
    else if (useBW[3]) physical = constrainedM3();
    else if (useBW[4]) physical = constrainedM4();
  }
  s3 = m3*m3;
  s4 = m4*m4;

  // Correct selected mass-spectrum to running-width Breit-Wigner.
  // Extra safety margin for maximum search.
  wtBW = 1.;
  if (useBW[3]) wtBW *= weightMass(3) * EXTRABWWTMAX;
  if (useBW[4]) wtBW *= weightMass(4) * EXTRABWWTMAX;

  // Done.
  return physical;

}


//--------------------------------------------------------------------------

// Select Breit-Wigner-distributed or fixed masses.

bool PhaseSpace2to2tauyz::trialMasses() {

  // By default vanishing cross section.
  sigmaNw = 0.;
  wtBW = 1.;

  // Pick m3 and m4 independently.
  trialMass(3);
  trialMass(4);

  // If outside phase space then reject event.
  if (m3 + m4 + MASSMARGIN > mHatMax) return false;

  // Correct selected mass-spectrum to running-width Breit-Wigner.
  if (useBW[3]) wtBW *= weightMass(3);
  if (useBW[4]) wtBW *= weightMass(4);

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Construct the four-vector kinematics from the trial values.

bool PhaseSpace2to2tauyz::finalKin() {

  // Assign masses to particles assumed massless in matrix elements.
  int id3 = sigmaProcessPtr->id(3);
  int id4 = sigmaProcessPtr->id(4);
  if (idMass[3] == 0) { m3 = particleDataPtr->m0(id3); s3 = m3*m3; }
  if (idMass[4] == 0) { m4 = particleDataPtr->m0(id4); s4 = m4*m4; }

  // Sometimes swap tHat <-> uHat to reflect chosen final-state order.
  if (sigmaProcessPtr->swappedTU()) {
    swap(tH, uH);
    z = -z;
  }

  // Check that phase space still open after new mass assignment.
  if (m3 + m4 + MASSMARGIN > mHat) {
    infoPtr->errorMsg("Warning in PhaseSpace2to2tauyz::finalKin: "
      "failed after mass assignment");
    return false;
  }
  p2Abs = 0.25 * (pow2(sH - s3 - s4) - 4. * s3 * s4) / sH;
  pAbs = sqrtpos( p2Abs );

  // Particle masses; incoming always on mass shell.
  mH[1] = 0.;
  mH[2] = 0.;
  mH[3] = m3;
  mH[4] = m4;

  // Special kinematics for direct photon+hadron (massless+massive) to fulfill
  // s = x1 * x2 * sHat and to retain the momentum of the massless photon beam.
  if ( hasPointGammaA && beamBPtr->isHadron() ) {
    double eCM1 = 0.5 * ( s + pow2(mA) - pow2(mB) ) / eCM;
    double eCM2 = 0.25 * x2H * s / eCM1;
    pH[1] = Vec4( 0., 0.,  eCM1, eCM1);
    pH[2] = Vec4( 0., 0., -eCM2, eCM2);
  } else if ( hasPointGammaB && beamAPtr->isHadron() ) {
    double eCM2 = 0.5 * ( s - pow2(mA) + pow2(mB) ) / eCM;
    double eCM1 = 0.25 * x1H * s / eCM2;
    pH[1] = Vec4( 0., 0.,  eCM1, eCM1);
    pH[2] = Vec4( 0., 0., -eCM2, eCM2);

  // Special kinematics for DIS to preserve lepton mass.
  } else if ( ( (beamAPtr->isLepton() && beamBPtr->isHadron())
             || (beamBPtr->isLepton() && beamAPtr->isHadron()) )
             && !settingsPtr->flag("PDF:lepton2gamma") ) {
    mH[1] = mA;
    mH[2] = mB;
    double pzAcm = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
      * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
    double eAcm  = sqrt( mH[1]*mH[1] + pzAcm*pzAcm);
    double pzBcm = -pzAcm;
    double eBcm  = sqrt( mH[2]*mH[2] + pzBcm*pzBcm);
    pH[1] = Vec4( 0., 0., pzAcm * x1H, eAcm * x1H);
    pH[2] = Vec4( 0., 0., pzBcm * x2H, eBcm * x2H);

  // Default kinematics with incoming partons along beam axes.
  } else {
    pH[1] = Vec4( 0., 0., 0.5 * eCM * x1H, 0.5 * eCM * x1H);
    pH[2] = Vec4( 0., 0., -0.5 * eCM * x2H, 0.5 * eCM * x2H);
  }

  // Outgoing partons initially in collision CM frame along beam axes.
  pH[3] = Vec4( 0., 0.,  pAbs, 0.5 * (sH + s3 - s4) / mHat);
  pH[4] = Vec4( 0., 0., -pAbs, 0.5 * (sH + s4 - s3) / mHat);

  // Then rotate and boost them to overall CM frame.
  theta = acos(z);
  phi   = 2. * M_PI * rndmPtr->flat();
  betaZ = (x1H - x2H)/(x1H + x2H);
  pH[3].rot( theta, phi);
  pH[4].rot( theta, phi);
  pH[3].bst( 0., 0., betaZ);
  pH[4].bst( 0., 0., betaZ);
  pTH = pAbs * sin(theta);

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Special choice of m3 and m4 when mHatMax push them off mass shell.
// Vary x in expression m3 + m4 = mHatMax - x * (Gamma3 + Gamma4).
// For each x try to put either 3 or 4 as close to mass shell as possible.
// Maximize BW_3 * BW_4 * beta_34, where latter approximate phase space.

bool PhaseSpace2to2tauyz::constrainedM3M4() {

  // Initial values.
  bool foundNonZero = false;
  double wtMassMax = 0.;
  double m3WtMax = 0.;
  double m4WtMax = 0.;
  double xMax = (mHatMax - mLower[3] - mLower[4]) / (mWidth[3] + mWidth[4]);
  double xStep = THRESHOLDSTEP * min(1., xMax);
  double xNow = 0.;
  double wtMassXbin, wtMassMaxOld, m34, mT34Min, wtMassNow,
    wtBW3Now, wtBW4Now, beta34Now;

  // Step through increasing x values.
  do {
    xNow += xStep;
    wtMassXbin = 0.;
    wtMassMaxOld = wtMassMax;
    m34 = mHatMax - xNow * (mWidth[3] + mWidth[4]);

    // Study point where m3 as close as possible to on-shell.
    m3 = min( mUpper[3], m34 - mLower[4]);
    if (m3 > mPeak[3]) m3 = max( mLower[3], mPeak[3]);
    m4 = m34 - m3;
    if (m4 < mLower[4]) {m4 = mLower[4]; m3 = m34 - m4;}

    // Check that inside phase space limit set by pTmin.
    mT34Min = sqrt(m3*m3 + pT2HatMin) + sqrt(m4*m4 + pT2HatMin);
    if (mT34Min < mHatMax) {

      // Breit-Wigners and beta factor give total weight.
      wtMassNow = 0.;
      if (m3 > mLower[3] && m3 < mUpper[3] && m4 > mLower[4]
        && m4 < mUpper[4]) {
        wtBW3Now = mw[3] / ( pow2(m3*m3 - sPeak[3]) + pow2(mw[3]) );
        wtBW4Now = mw[4] / ( pow2(m4*m4 - sPeak[4]) + pow2(mw[4]) );
        beta34Now = sqrt( pow2(mHatMax*mHatMax - m3*m3 - m4*m4)
          - pow2(2. * m3 * m4) ) / (mHatMax*mHatMax);
        wtMassNow = wtBW3Now * wtBW4Now * beta34Now;
      }

      // Store new maximum, if any.
      if (wtMassNow > wtMassXbin) wtMassXbin = wtMassNow;
      if (wtMassNow > wtMassMax) {
        foundNonZero = true;
        wtMassMax = wtMassNow;
        m3WtMax = m3;
        m4WtMax = m4;
      }
    }

    // Study point where m4 as close as possible to on-shell.
    m4 = min( mUpper[4], m34 - mLower[3]);
    if (m4 > mPeak[4]) m4 = max( mLower[4], mPeak[4]);
    m3 = m34 - m4;
    if (m3 < mLower[3]) {m3 = mLower[3]; m4 = m34 - m3;}

    // Check that inside phase space limit set by pTmin.
    mT34Min = sqrt(m3*m3 + pT2HatMin) + sqrt(m4*m4 + pT2HatMin);
    if (mT34Min < mHatMax) {

      // Breit-Wigners and beta factor give total weight.
      wtMassNow = 0.;
      if (m3 > mLower[3] && m3 < mUpper[3] && m4 > mLower[4]
        && m4 < mUpper[4]) {
        wtBW3Now = mw[3] / ( pow2(m3*m3 - sPeak[3]) + pow2(mw[3]) );
        wtBW4Now = mw[4] / ( pow2(m4*m4 - sPeak[4]) + pow2(mw[4]) );
        beta34Now = sqrt( pow2(mHatMax*mHatMax - m3*m3 - m4*m4)
          - pow2(2. * m3 * m4) ) / (mHatMax*mHatMax);
        wtMassNow = wtBW3Now * wtBW4Now * beta34Now;
      }

      // Store new maximum, if any.
      if (wtMassNow > wtMassXbin) wtMassXbin = wtMassNow;
      if (wtMassNow > wtMassMax) {
        foundNonZero = true;
        wtMassMax = wtMassNow;
        m3WtMax = m3;
        m4WtMax = m4;
      }
    }

  // Continue stepping if increasing trend and more x range available.
  } while ( (!foundNonZero || wtMassXbin > wtMassMaxOld)
    && xNow < xMax - xStep);

  // Restore best values for subsequent maximization. Return.
  m3 = m3WtMax;
  m4 = m4WtMax;
  return foundNonZero;

}

//--------------------------------------------------------------------------

// Special choice of m3 when mHatMax pushes it off mass shell.
// Vary x in expression m3 = mHatMax - m4 - x * Gamma3.
// Maximize BW_3 * beta_34, where latter approximate phase space.

bool PhaseSpace2to2tauyz::constrainedM3() {

  // Initial values.
  bool foundNonZero = false;
  double wtMassMax = 0.;
  double m3WtMax = 0.;
  double mT4Min = sqrt(m4*m4 + pT2HatMin);
  double xMax = (mHatMax - mLower[3] - m4) / mWidth[3];
  double xStep = THRESHOLDSTEP * min(1., xMax);
  double xNow = 0.;
  double wtMassNow, mT34Min, wtBW3Now, beta34Now;

  // Step through increasing x values; gives m3 unambiguously.
  do {
    xNow += xStep;
    wtMassNow = 0.;
    m3 = mHatMax - m4 - xNow * mWidth[3];

    // Check that inside phase space limit set by pTmin.
    mT34Min = sqrt(m3*m3 + pT2HatMin) + mT4Min;
    if (mT34Min < mHatMax) {

      // Breit-Wigner and beta factor give total weight.
      wtBW3Now = mw[3] / ( pow2(m3*m3 - sPeak[3]) + pow2(mw[3]) );
      beta34Now = sqrt( pow2(mHatMax*mHatMax - m3*m3 - m4*m4)
        - pow2(2. * m3 * m4) ) / (mHatMax*mHatMax);
      wtMassNow = wtBW3Now * beta34Now;

      // Store new maximum, if any.
      if (wtMassNow > wtMassMax) {
        foundNonZero = true;
        wtMassMax = wtMassNow;
        m3WtMax = m3;
      }
    }

  // Continue stepping if increasing trend and more x range available.
  } while ( (!foundNonZero || wtMassNow > wtMassMax)
    && xNow < xMax - xStep);

  // Restore best value for subsequent maximization. Return.
  m3 = m3WtMax;
  return foundNonZero;

}

//--------------------------------------------------------------------------

// Special choice of m4 when mHatMax pushes it off mass shell.
// Vary x in expression m4 = mHatMax - m3 - x * Gamma4.
// Maximize BW_4 * beta_34, where latter approximate phase space.

bool PhaseSpace2to2tauyz::constrainedM4() {

  // Initial values.
  bool foundNonZero = false;
  double wtMassMax = 0.;
  double m4WtMax = 0.;
  double mT3Min = sqrt(m3*m3 + pT2HatMin);
  double xMax = (mHatMax - mLower[4] - m3) / mWidth[4];
  double xStep = THRESHOLDSTEP * min(1., xMax);
  double xNow = 0.;
  double wtMassNow, mT34Min, wtBW4Now, beta34Now;

  // Step through increasing x values; gives m4 unambiguously.
  do {
    xNow += xStep;
    wtMassNow = 0.;
    m4 = mHatMax - m3 - xNow * mWidth[4];

    // Check that inside phase space limit set by pTmin.
    mT34Min = mT3Min + sqrt(m4*m4 + pT2HatMin);
    if (mT34Min < mHatMax) {

      // Breit-Wigner and beta factor give total weight.
      wtBW4Now = mw[4] / ( pow2(m4*m4 - sPeak[4]) + pow2(mw[4]) );
      beta34Now = sqrt( pow2(mHatMax*mHatMax - m3*m3 - m4*m4)
        - pow2(2. * m3 * m4) ) / (mHatMax*mHatMax);
      wtMassNow = wtBW4Now * beta34Now;

      // Store new maximum, if any.
      if (wtMassNow > wtMassMax) {
        foundNonZero = true;
        wtMassMax = wtMassNow;
        m4WtMax = m4;
      }
    }

  // Continue stepping if increasing trend and more x range available.
  } while ( (!foundNonZero || wtMassNow > wtMassMax)
    && xNow < xMax - xStep);

  // Restore best value for subsequent maximization.
  m4 = m4WtMax;
  return foundNonZero;

}

//--------------------------------------------------------------------------

// Calculate the cross section with rescaled sHat.

void PhaseSpace2to2tauyz::rescaleSigma(double sHatNew){

  // With massless matrix element derive tHat without masses.
  if ( idMass[3] == 0 ) s3 = 0.;
  if ( idMass[4] == 0 ) s4 = 0.;

  // Update variables according to new sHat.
  sH    = sHatNew;
  double sH34 = -0.5 * (sH - s3 - s4);
  p2Abs = 0.25 * (pow2(sH - s3 - s4) - 4. * s3 * s4) / sH;
  pAbs  = sqrtpos( p2Abs );
  mHat  = sqrt(sH);
  tH    = sH34 + mHat * pAbs * z;
  uH    = sH34 - mHat * pAbs * z;
  pTH   = sqrtpos( (tH * uH - s3 * s4) / sH);

  // Calculate the cross section for the process with rescaled kinematics
  // if original cross section non-zero.
  if (sigmaNw > TINY) {
    sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4, runBW3H, runBW4H);
    sigmaNw  = sigmaProcessPtr->sigmaPDF(false, true);
    sigmaNw *= wtTau * wtY * wtZ * wtBW;
    if (canBias2Sel) sigmaNw *= pow( pTH / bias2SelRef, bias2SelPow);
  }

}

//--------------------------------------------------------------------------

// Rescales the momenta of incoming and outgoing partons according to
// new sHat sampled in GammaKinematics.

void PhaseSpace2to2tauyz::rescaleMomenta( double sHatNew){

  // Loop over initial and final partons.
  for (int i = 0; i <= 1; ++i){

    // Either final or initial partons.
    int iPartonA = (i == 0) ? 1 : 3;
    int iPartonB = (i == 0) ? 2 : 4;

    // Original momenta of partons.
    Vec4 pA = p( iPartonA);
    Vec4 pB = p( iPartonB);

    // Calculate new momenta in CM-frame.
    double m2A = pow2( m( iPartonA) );
    double m2B = pow2( m( iPartonB) );
    double eA  = 0.5 * ( sHatNew + m2A - m2B) / sqrt( sHatNew);
    double eB  = 0.5 * ( sHatNew + m2B - m2A) / sqrt( sHatNew);
    double pz  = 0.5 * sqrtpos( pow2(sHatNew - m2A - m2B) - 4.0 * m2A * m2B )
               / sqrt( sHatNew);
    Vec4 pANew( 0, 0,  pz, eA );
    Vec4 pBNew( 0, 0, -pz, eB );

    // Find boost to original frame.
    RotBstMatrix MtoCMinc;
    MtoCMinc.toCMframe( pA, pB);
    MtoCMinc.invert();

    // Boost outgoing partons to original frame and replace the momenta.
    pANew.rotbst( MtoCMinc);
    pBNew.rotbst( MtoCMinc);
    setP( iPartonA, pANew);
    setP( iPartonB, pBNew);
  }
}

//--------------------------------------------------------------------------

// Calculates the cross-section weight for overestimated photon flux.

double PhaseSpace2to2tauyz::weightGammaPDFApprox(){

  // No need for reweighting if only direct photons.
  if (beamAPtr->getGammaMode() == 2 && beamBPtr->getGammaMode() == 2)
    return 1.;
  if ( (beamAPtr->getGammaMode() == 2 && beamBPtr->isHadron())
       || (beamBPtr->getGammaMode() == 2 && beamAPtr->isHadron()) )
    return 1.;

  // Get the combined x and x_gamma values and derive x'.
  double x1GammaHadr = beamAPtr->xGammaHadr();
  double x2GammaHadr = beamBPtr->xGammaHadr();
  double x1Gamma     = beamAPtr->xGamma();
  double x2Gamma     = beamBPtr->xGamma();
  double x1Hadr      = x1GammaHadr / x1Gamma;
  double x2Hadr      = x2GammaHadr / x2Gamma;

  // For photon-hadron case do not reweight the hadron side.
  if ( beamAPtr->isHadron() || beamAPtr->getGammaMode() == 2 ) {
    x1GammaHadr = -1.;
    x1Gamma     = -1.;
  }
  if ( beamBPtr->isHadron() || beamBPtr->getGammaMode() == 2 ) {
    x2GammaHadr = -1.;
    x2Gamma     = -1.;
  }

  // Calculate the over-estimated PDF convolution and the current one.
  double sigmaOver = sigmaProcessPtr->sigmaPDF(false, false, true,
                                               x1GammaHadr, x2GammaHadr);
  double sigmaCorr = sigmaProcessPtr->sigmaPDF(false, false, true,
                                               x1Hadr, x2Hadr);

  // Make sure that the overestimate is finite.
  if (sigmaOver < TINY) return 0.;

  // Return weight.
  return sigmaCorr / sigmaOver;

}

//==========================================================================

// PhaseSpace2to2elastic class.
// 2 -> 2 kinematics set up for elastic scattering.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Max number of tries to find acceptable t.
const int    PhaseSpace2to2elastic::NTRY     = 1000;

// Conversion coefficient (mb <-> GeV^2).
const double PhaseSpace2to2elastic::HBARC2   = 0.38938;

// Width and relative importance of two exponentials.
const double PhaseSpace2to2elastic::BNARROW  = 10.;
const double PhaseSpace2to2elastic::BWIDE    = 1.;
const double PhaseSpace2to2elastic::WIDEFRAC = 0.1;
const double PhaseSpace2to2elastic::TOFFSET  = -0.2;

//--------------------------------------------------------------------------

// Form of phase space sampling already fixed, so no optimization.
// However, need to read out relevant parameters from SigmaTotal.

bool PhaseSpace2to2elastic::setupSampling() {

  // Flag if a photon inside lepton beam.
  hasGamma = settingsPtr->flag("PDF:lepton2gamma");

  // Flag if photon has a VMD state.
  hasVMD = infoPtr->isVMDstateA() || infoPtr->isVMDstateB();

  // If not photoproduction, calculate the cross-section estimates directly.
  if (!hasGamma) {

    // Find maximum = value of cross section.
    sigmaNw    = sigmaProcessPtr->sigmaHatWrap();

  // For photoproduction calculate the estimates for photon-hadron system.
  } else {

    // Total cross section using a photon instead of the actual beam.
    idAgm   = gammaKinPtr->idInA();
    idBgm   = gammaKinPtr->idInB();
    sigmaTotPtr->calc( idAgm, idBgm, eCM);
    sigmaProcessPtr->setIdInDiff(idAgm, idBgm);

    // Zero mass for photons from lepton beam.
    if (idAgm == 22) mA = 0.;
    if (idBgm == 22) mB = 0.;

    // Use the elastic cross section for overestimate.
    sigmaMxGm = sigmaTotPtr->sigmaEl();
    sigmaNw   = gammaKinPtr->setupSoftPhaseSpaceSampling(sigmaMxGm);
  }

  // Save the maximum value.
  sigmaMx    = sigmaNw;

  // Character of elastic generation.
  isOneExp   = sigmaTotPtr->bElIsExp();
  useCoulomb = sigmaTotPtr->hasCoulomb();
  alphaEM0   = settingsPtr->parm("StandardModel:alphaEM0");

  // Squared and outgoing masses of particles.
  // Recalculated later if photon fluctuates into different vector mesons.
  s1         = mA * mA;
  s2         = mB * mB;
  m3         = mA;
  m4         = mB;

  // Determine maximum possible t range.
  lambda12S  = pow2(s - s1 - s2) - 4. * s1 * s2 ;
  tLow       = - lambda12S / s;
  tUpp       = (useCoulomb) ? -settingsPtr->parm("SigmaElastic:tAbsMin") : 0.;

  // Upper estimate as sum of two exponentials and a Coulomb.
  // VMD: Start with bNarrow but recalculate when VM state sampled in trialKin.
  bSlope1    = (isOneExp && !hasVMD) ? sigmaTotPtr->bSlopeEl() : BNARROW;
  bSlope2    = BWIDE;
  sigRef1    = sigmaTotPtr->dsigmaEl( tUpp, false);
  if (isOneExp) {
    sigNorm1 = sigRef1 / bSlope1;
    if (useCoulomb) sigNorm1 *= 2.;
    sigNorm2 = 0.;
  } else {
    sigRef2  = sigmaTotPtr->dsigmaEl( tUpp + TOFFSET, false);
    sigRef   = (sigRef1 > 2. * sigRef2) ? 2. * sigRef1 : 5. * sigRef2;
    rel2     = exp((bSlope2 - bSlope1) * tUpp) * WIDEFRAC / (1. - WIDEFRAC);
    sigNorm1 = sigRef / (bSlope1 + rel2 * bSlope2);
    sigNorm2 = sigNorm1 * rel2;
  }
  sigNorm3   = (useCoulomb) ? -2. * HBARC2 * 4. * M_PI * pow2(alphaEM0)
               / tUpp : 0.;
  sigNormSum = sigNorm1 + sigNorm2 + sigNorm3;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select a trial kinematics phase space point. Perform full
// Monte Carlo acceptance/rejection at this stage.

bool PhaseSpace2to2elastic::trialKin( bool, bool ) {

  // Allow for possibility that energy varies from event to event.
  if (doEnergySpread) {
    eCM       = infoPtr->eCM();
    s         = eCM * eCM;
    if (!hasVMD) {
      lambda12S = pow2(s - s1 - s2) - 4. * s1 * s2 ;
      tLow      = - lambda12S / s;
    }
  }

  // Sample kinematics for gamma+gamma(hadron) sub-event and reject
  // to account for over sampling.
  if (hasGamma) {

    // Current weight.
    double wt = 1.0;

    // Sample gamma kinematics.
    if (!gammaKinPtr->trialKinSoftPhaseSpaceSampling() ) return false;

    // Calculate the cross sections with invariant mass of the sub-system.
    double eCMsub = gammaKinPtr->eCMsub();
    sigmaTotPtr->calc( idAgm, idBgm, eCMsub );

    // Correct for the over-estimated sigma.
    double sigmaW  = sigmaTotPtr->sigmaEl();
    double wtSigma = sigmaW/sigmaMxGm;

    // Calculate the total weight and warn if unphysical weight.
    wt *= wtSigma * gammaKinPtr->weight();
    if ( wt > 1. ) infoPtr->errorMsg("Warning in PhaseSpace2to2elastic::"
      "trialKin: weight above unity");

    // Correct for over-estimated cross section and x_gamma limits.
    if ( wt < rndmPtr->flat() ) return false;

    // For accepted kinematics use the sub-collision energy.
    eCM       = eCMsub;
    s         = eCM * eCM;
    if (!hasVMD) {
      lambda12S = pow2( s - s1 - s2) - 4. * s1 * s2 ;
      tLow      = - lambda12S / s;
    }
  }

  // Elastically scattered photon in vector-meson state.
  if (hasVMD) {

    // Sample VMD states and calculate kinematics with vector-meson mass(es).
    if (hasGamma) sigmaTotPtr->chooseVMDstates(idAgm, idBgm, eCM, 102);
    else          sigmaTotPtr->chooseVMDstates(idA, idB, eCM, 102);
    m3 = (infoPtr->isVMDstateA()) ? infoPtr->mVMDA() : mA;
    m4 = (infoPtr->isVMDstateB()) ? infoPtr->mVMDB() : mB;
    s3 = m3 * m3;
    s4 = m4 * m4;

    // New kinematics with sampled VM masses and possibly new energy.
    lambda12 = sqrtpos( pow2( s - s1 - s2) - 4. * s1 * s2 );
    lambda34 = sqrtpos( pow2( s - s3 - s4) - 4. * s3 * s4 );
    tempA    = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4) / s;
    tempB    = lambda12 * lambda34 / s;
    tempC    = (s3 - s1) * (s4 - s2) + (s1 + s4 - s2 - s3)
             * (s1 * s4 - s2 * s3) / s;
    tLow     = -0.5 * (tempA + tempB);
    tUpp     = (useCoulomb) ? -settingsPtr->parm("SigmaElastic:tAbsMin")
             : tempC / tLow;

    // Recalculate the elastic cross section with the sampled state
    // to save bEl.
    int idAVMD = infoPtr->isVMDstateA() ? infoPtr->idVMDA() : idA;
    int idBVMD = infoPtr->isVMDstateB() ? infoPtr->idVMDB() : idB;
    sigmaTotPtr->calcTotEl(idAVMD, idBVMD, s, m3, m4);

    // Find the b-slope for the selected VM state and derive new sigmaNorm.
    bSlope1    = (isOneExp) ? sigmaTotPtr->bSlopeEl() : BNARROW;
    bSlope2    = BWIDE;
    sigRef1    = sigmaTotPtr->dsigmaEl( tUpp, false);

    // Need to recalculate the cross section estimates for the given VM state.
    if (isOneExp) {
      sigNorm1 = sigRef1 / bSlope1;
      if (useCoulomb) sigNorm1 *= 2.;
      sigNorm2 = 0.;
    } else {
      sigRef2  = sigmaTotPtr->dsigmaEl( tUpp + TOFFSET, false);
      sigRef   = (sigRef1 > 2. * sigRef2) ? 2. * sigRef1 : 5. * sigRef2;
      rel2     = exp((bSlope2 - bSlope1) * tUpp) * WIDEFRAC / (1. - WIDEFRAC);
      sigNorm1 = sigRef / (bSlope1 + rel2 * bSlope2);
      sigNorm2 = sigNorm1 * rel2;
    }
    sigNorm3   = (useCoulomb) ? -2. * HBARC2 * 4. * M_PI * pow2(alphaEM0)
               / tUpp : 0.;
    sigNormSum = sigNorm1 + sigNorm2 + sigNorm3;
  }

  // Repeated tries until accepted.
  double rNow, bNow, sigNow, sigEst;
  int loop = 0;
  do {
    ++loop;
    if (loop == NTRY) {
      infoPtr->errorMsg("Error in PhaseSpace2to2elastic::trialKin: "
        " quit after repeated tries");
      return false;
    }
    rNow = rndmPtr->flat() * sigNormSum;
    if (useCoulomb && rNow > sigNorm1 + sigNorm2) tH = tUpp / rndmPtr->flat();
    else {
      bNow = (rNow < sigNorm1) ? bSlope1 : bSlope2;
      tH   = tUpp + log( rndmPtr->flat() ) / bNow;
    }
    sigNow = sigmaTotPtr->dsigmaEl( tH, useCoulomb);
    sigEst = sigNorm1 * bSlope1 * exp( bSlope1 * (tH - tUpp))
           + sigNorm2 * bSlope2 * exp( bSlope2 * (tH - tUpp));
    if (useCoulomb) sigEst += sigNorm3 * (-tUpp) / pow2(tH);
  } while (tH < tLow || sigNow < sigEst * rndmPtr->flat());
  if (sigNow > 1.01 * sigEst) infoPtr->errorMsg("Warning in "
    "PhaseSpace2to2elastic::trialKin: cross section maximum violated");

  if (!hasVMD){
    // Careful reconstruction of scattering angle.
    double tRat     = s * tH / lambda12S;
    double cosTheta = min(1., max(-1., 1. + 2. * tRat ) );
    double sinTheta = 2. * sqrtpos( -tRat * (1. + tRat) );
    theta           = asin( min(1., sinTheta));
    if (cosTheta < 0.) theta = M_PI - theta;
  } else {
    // Careful reconstruction of scattering angle.
    double cosTheta = min(1., max(-1., (tempA + 2. * tH) / tempB));
    double sinTheta = 2. * sqrtpos( -(tempC + tempA * tH + tH * tH) ) / tempB;
    theta = asin( min(1., sinTheta));
    if (cosTheta < 0.) theta = M_PI - theta;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Construct the four-vector kinematics from the trial values.

bool PhaseSpace2to2elastic::finalKin() {

  // Particle masses.
  mH[1] = mA;
  mH[2] = mB;
  mH[3] = m3;
  mH[4] = m4;

  // Use smapled vector-meson masses when relevant.
  if (!hasVMD) {

    // Incoming particles along beam axes.
    pAbs = 0.5 * sqrtpos(lambda12S) / eCM;
    pH[1] = Vec4( 0., 0.,  pAbs, 0.5 * (s + s1 - s2) / eCM);
    pH[2] = Vec4( 0., 0., -pAbs, 0.5 * (s + s2 - s1) / eCM);

    // Outgoing particles initially along beam axes.
    pH[3] = Vec4( 0., 0.,  pAbs, 0.5 * (s + s1 - s2) / eCM);
    pH[4] = Vec4( 0., 0., -pAbs, 0.5 * (s + s2 - s1) / eCM);

  } else {

    // Incoming particles along beam axes.
    pAbs = 0.5 * lambda12 / eCM;
    pH[1] = Vec4( 0., 0.,  pAbs, 0.5 * (s + s1 - s2) / eCM);
    pH[2] = Vec4( 0., 0., -pAbs, 0.5 * (s + s2 - s1) / eCM);

    // Outgoing particles initially along beam axes.
    pAbs = 0.5 * lambda34 / eCM;
    pH[3] = Vec4( 0., 0.,  pAbs, 0.5 * (s + s3 - s4) / eCM);
    pH[4] = Vec4( 0., 0., -pAbs, 0.5 * (s + s4 - s3) / eCM);

  }

  // Then rotate them
  phi = 2. * M_PI * rndmPtr->flat();
  pH[3].rot( theta, phi);
  pH[4].rot( theta, phi);

  // Set some further info for completeness.
  x1H = 1.;
  x2H = 1.;
  sH = s;
  uH = 2. * (s1 + s2) - sH - tH;
  mHat = eCM;
  p2Abs = pAbs * pAbs;
  betaZ = 0.;
  pTH = pAbs * sin(theta);

  // Save the sampled photon kinematics.
  if (hasGamma) gammaKinPtr->finalize();

  // Done.
  return true;

}

//==========================================================================

// PhaseSpace2to2diffractive class.
// 2 -> 2 kinematics set up for diffractive scattering.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of tries to find acceptable (m^2, t) set.
const int PhaseSpace2to2diffractive::NTRY          = 2500;

// t is sampled according to sum of four distributions.
const double PhaseSpace2to2diffractive::BWID1      = 8.;
const double PhaseSpace2to2diffractive::BWID2      = 2.;
const double PhaseSpace2to2diffractive::BWID3      = 0.5;
const double PhaseSpace2to2diffractive::BWID4      = 0.2;
const double PhaseSpace2to2diffractive::FWID1SD    = 1.;
const double PhaseSpace2to2diffractive::FWID2SD    = 0.2;
const double PhaseSpace2to2diffractive::FWID3SD    = 0.1;
const double PhaseSpace2to2diffractive::FWID4SD    = 0.1;
const double PhaseSpace2to2diffractive::FWID1DD    = 0.1;
const double PhaseSpace2to2diffractive::FWID2DD    = 1.;
const double PhaseSpace2to2diffractive::FWID3DD    = 0.5;
const double PhaseSpace2to2diffractive::FWID4DD    = 0.2;

// Safety margins for upper estimate of cross section.
const double PhaseSpace2to2diffractive::MAXFUDGESD = 2.;
const double PhaseSpace2to2diffractive::MAXFUDGEDD = 2.;
const double PhaseSpace2to2diffractive::MAXFUDGET  = 4.;

// Safety margin so sum of masses not too close to eCM.
const double PhaseSpace2to2diffractive::DIFFMASSMARGIN = 0.2;

// Squared proton mass.
const double PhaseSpace2to2diffractive::SPROTON = 0.8803544;

//--------------------------------------------------------------------------

// Form of phase space sampling already fixed, so no optimization.
// However, need to find upper estimate at t = 0.

bool PhaseSpace2to2diffractive::setupSampling() {

  // Flag if a photon inside lepton beam.
  hasGamma = settingsPtr->flag("PDF:lepton2gamma");

  // Flag if photon has a VMD state.
  hasVMD = infoPtr->isVMDstateA() || infoPtr->isVMDstateB();

  // If not photoproduction, calculate the cross-section estimates directly.
  if (!hasGamma) {

    // Find maximum = value of cross section.
    sigmaNw    = sigmaProcessPtr->sigmaHatWrap();

  // For photoproduction calculate the estimates for photon-hadron system.
  } else {

    // Total cross section using a photon instead of the actual beam.
    idAgm   = gammaKinPtr->idInA();
    idBgm   = gammaKinPtr->idInB();
    sigmaTotPtr->calc( idAgm, idBgm, eCM);
    sigmaProcessPtr->setIdInDiff(idAgm, idBgm);

    // Zero mass for photons from lepton beam.
    if (idAgm == 22) mA = 0.;
    if (idBgm == 22) mB = 0.;

    // Use the correct diffractive cross section for overestimate.
    sigmaMxGm = 0.;
    if      (isDiffA && isSD)    sigmaMxGm = sigmaTotPtr->sigmaAX();
    else if (isDiffB && isSD)    sigmaMxGm = sigmaTotPtr->sigmaXB();
    else if (isDiffB && isDiffA) sigmaMxGm = sigmaTotPtr->sigmaXX();
    sigmaNw   = gammaKinPtr->setupSoftPhaseSpaceSampling(sigmaMxGm);
  }

  // Save the maximum value.
  sigmaMx    = sigmaNw;

  // Masses of particles and minimal masses of diffractive states.
  // COR: Take VMD states into account already here, because of maximal cross
  // section calculation below. Minimal VMD mass is the rho mass.
  double mPi   = particleDataPtr->m0(211);
  double mRho  = particleDataPtr->m0(113);
  double mAtmp = (infoPtr->isVMDstateA()) ? mRho : mA;
  double mBtmp = (infoPtr->isVMDstateB()) ? mRho : mB;
  m3ElDiff     = (isDiffA) ? mAtmp + mPi : mAtmp;
  m4ElDiff     = (isDiffB) ? mBtmp + mPi : mBtmp;
  s1           = mA * mA;
  s2           = mB * mB;
  s3           = pow2( m3ElDiff);
  s4           = pow2( m4ElDiff);

  // Initial kinematics value.
  lambda12 = sqrtpos( pow2( s - s1 - s2) - 4. * s1 * s2 );

  // Scenarios with separate handling of xi and t (currently only MBR).
  // Step 0 = both xi and t, 1 = xi only, 2 = t only.
  splitxit = sigmaTotPtr->splitDiff();
  int step = (splitxit) ? 1 : 0;

  // Find maximal cross section xi * dsigma / (dxi dt) at t = 0.
  sigMax     = 0.;
  if (isSD) {
    xiMin      = (isDiffA) ? s3 / s : s4 / s;
    for (int i = 0; i < 100; ++i) {
      xiNow  = pow( xiMin, 0.01 * i + 0.005);
      sigNow = sigmaTotPtr->dsigmaSD( xiNow, 0., isDiffA, step);
      if (sigNow > sigMax) sigMax = sigNow;
    }

  // Find maximal cross section xi1 * xi2 * dsigma / (dxi1 dxi2 dt) at t = 0.
  } else {
    xiMin      = max( s3, s4) / s;
    xiMax      = sqrt( SPROTON / s);
    for (int i = 0; i < 100; ++i) {
      xiNow  = xiMin * pow( xiMax / xiMin, 0.01 * i + 0.005);
      sigNow = sigmaTotPtr->dsigmaDD( xiNow, xiNow, 0., step);
      if (sigNow > sigMax) sigMax = sigNow;
    }
  }
  sigMax *= (isSD ? MAXFUDGESD : MAXFUDGEDD);

  // Combinations of t sampling parameters.
  fWid1    = (isSD ? FWID1SD : FWID1DD);
  fWid2    = (isSD ? FWID2SD : FWID2DD);
  fWid3    = (isSD ? FWID3SD : FWID3DD);
  fWid4    = (isSD ? FWID4SD : FWID4DD);
  fbWid1   = fWid1 * BWID1;
  fbWid2   = fWid2 * BWID2;
  fbWid3   = fWid3 * BWID3;
  fbWid4   = fWid4 * BWID4;
  fbWid1234 = fbWid1 + fbWid2 + fbWid3 + fbWid4;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select a trial kinematics phase space point. Perform full
// Monte Carlo acceptance/rejection at this stage.

bool PhaseSpace2to2diffractive::trialKin( bool, bool ) {

  // Allow for possibility that energy varies from event to event.
  if (doEnergySpread) {
    eCM       = infoPtr->eCM();
    s         = eCM * eCM;
    lambda12 = sqrtpos( pow2( s - s1 - s2) - 4. * s1 * s2 );
  }

  // Sample kinematics for gamma+gamma(hadron) sub-event and reject
  // to account for over sampling.
  if (hasGamma) {

    // Current weight.
    double wt = 1.0;

    // Sample gamma kinematics.
    if (!gammaKinPtr->trialKinSoftPhaseSpaceSampling() ) return false;

    // Calculate the cross sections with invariant mass of the sub-system.
    double eCMsub = gammaKinPtr->eCMsub();
    sigmaTotPtr->calc( idAgm, idBgm, eCMsub );

    // Correct for the over-estimated sigmaZZ.
    double sigmaW = 0.;
    if      (isDiffA && isSD)    sigmaW = sigmaTotPtr->sigmaAX();
    else if (isDiffB && isSD)    sigmaW = sigmaTotPtr->sigmaXB();
    else if (isDiffB && isDiffA) sigmaW = sigmaTotPtr->sigmaXX();
    double wtSigma = sigmaW/sigmaMxGm;

    // Calculate the total weight and warn if unphysical weight.
    wt *= wtSigma * gammaKinPtr->weight();
    if ( wt > 1. ) infoPtr->errorMsg("Warning in PhaseSpace2to2diffractive::"
      "trialKin: weight above unity");

    // Correct for over-estimated cross section and x_gamma limits.
    if ( wt < rndmPtr->flat() ) return false;

    // For accepted kinematics use the sub-collision energy.
    eCM       = eCMsub;
    s         = eCM * eCM;
    lambda12  = sqrtpos( pow2( s - s1 - s2) - 4. * s1 * s2 );
  }

  // Sample VMD states with possibly different CM energy.
  double mAtmp = mA;
  double mBtmp = mB;
  if (hasVMD) {

    // Find the diffractive process.
    int processCode = 101;
    if      (isDiffA && isSD)    processCode = 104;
    else if (isDiffB && isSD)    processCode = 103;
    else if (isDiffB && isDiffA) processCode = 105;

    // Sampled VM states.
    if (hasGamma) sigmaTotPtr->chooseVMDstates(idAgm, idBgm, eCM, processCode);
    else          sigmaTotPtr->chooseVMDstates(idA, idB, eCM, processCode);

    // Now choose proper VMD mass. Special handling for minimal
    // diffractive mass for J/Psi as we require at least two D-mesons to
    // be produced by string breaking.
    double mPi   = particleDataPtr->m0(211);
    double mD    = particleDataPtr->m0(411);
    mAtmp        = (infoPtr->isVMDstateA()) ? infoPtr->mVMDA() : mA;
    mBtmp        = (infoPtr->isVMDstateB()) ? infoPtr->mVMDB() : mB;
    m3ElDiff     = (isDiffA) ? mAtmp + mPi : mAtmp;
    m4ElDiff     = (isDiffB) ? mBtmp + mPi : mBtmp;
    if (isDiffA && infoPtr->idVMDA() == 443) m3ElDiff = 2. * mD;
    if (isDiffB && infoPtr->idVMDB() == 443) m4ElDiff = 2. * mD;
    s3           = pow2( m3ElDiff);
    s4           = pow2( m4ElDiff);

  }

  // Normally xi and t in one step, but possible to split into two.
  int nStep = (splitxit) ? 2 : 1;
  for (int iStep = 0; iStep < nStep; ++iStep) {
    int step = (splitxit) ? iStep + 1 : 0;

    // Loop over attempts to set up masses and t consistently.
    for (int loop = 0; ; ++loop) {
      if (loop == NTRY) {
        infoPtr->errorMsg("Error in PhaseSpace2to2diffractive::trialKin: "
          " quit after repeated tries");
        return false;
      }

      // Select diffractive mass(es) according to dm^2/m^2. COR change
      // from mA -> mAtmp?
      if (iStep == 0) {
        m3 = (isDiffA) ? m3ElDiff * pow( max(mAtmp, eCM - m4ElDiff) / m3ElDiff,
          rndmPtr->flat()) : m3ElDiff;
        m4 = (isDiffB) ? m4ElDiff * pow( max(mBtmp, eCM - m3ElDiff) / m4ElDiff,
          rndmPtr->flat()) : m4ElDiff;
        if (m3 + m4 + DIFFMASSMARGIN >= eCM) continue;
        s3 = m3 * m3;
        s4 = m4 * m4;
      }

      // Select t according to exp(b*t), b picked among four options.
      if (step != 1) {
        double pickb = rndmPtr->flat() * (fWid1 + fWid2 + fWid3 + fWid4);
        bNow = (pickb < fWid1) ? BWID1
           : ( (pickb < fWid1 + fWid2) ? BWID2
           : ( (pickb < fWid1 + fWid2 + fWid3) ? BWID3 : BWID4 ) );
        tH   = log(rndmPtr->flat()) / bNow;

        // Check whether m^2 and t choices are consistent.
        lambda34 = sqrtpos( pow2( s - s3 - s4) - 4. * s3 * s4 );
        tempA    = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4) / s;
        tempB    = lambda12 *  lambda34 / s;
        tempC    = (s3 - s1) * (s4 - s2) + (s1 + s4 - s2 - s3)
                 * (s1 * s4 - s2 * s3) / s;
        tLow     = -0.5 * (tempA + tempB);
        tUpp     = tempC / tLow;
        if (tH < tLow || tH > tUpp) continue;
      }

      // Evaluate single or double diffractive cross section.
      if (isSD) {
        xiNow     = (isDiffA) ? s3 / s : s4 / s;
        sigNow    = sigmaTotPtr->dsigmaSD( xiNow, tH, isDiffA, step);
      } else {
        sigNow    = sigmaTotPtr->dsigmaDD( s3 / s, s4 / s, tH, step);
      }

      // Maximum weight based on sampling strategy.
      tWeight   = ( fbWid1 * exp( BWID1 * tH) + fbWid2 * exp(BWID2 * tH)
        + fbWid3 * exp(BWID3 * tH) + fbWid4 * exp(BWID4 * tH) ) / fbWid1234;
      sigMaxNow = (step == 0) ? sigMax * tWeight
                : ( (step == 1) ? sigMax : MAXFUDGET * tWeight );

      // Check for maximum violations. Possibly break out of the loop.
      if (sigNow > sigMaxNow) infoPtr->errorMsg("Error in PhaseSpace2to2"
        "diffractive::trialKin: maximum cross section violated");
      if (sigNow > rndmPtr->flat() * sigMaxNow) break;

    // End of loops over tries and steps.
    }
  }

  // Careful reconstruction of scattering angle.
  double cosTheta = min(1., max(-1., (tempA + 2. * tH) / tempB));
  double sinTheta = 2. * sqrtpos( -(tempC + tempA * tH + tH * tH) ) / tempB;
  theta = asin( min(1., sinTheta));
  if (cosTheta < 0.) theta = M_PI - theta;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Construct the four-vector kinematics from the trial values.

bool PhaseSpace2to2diffractive::finalKin() {

  // Particle masses; incoming always on mass shell.
  mH[1] = mA;
  mH[2] = mB;
  mH[3] = m3;
  mH[4] = m4;

  // Incoming particles along beam axes.
  pAbs = 0.5 * lambda12 / eCM;
  pH[1] = Vec4( 0., 0.,  pAbs, 0.5 * (s + s1 - s2) / eCM);
  pH[2] = Vec4( 0., 0., -pAbs, 0.5 * (s + s2 - s1) / eCM);

  // Outgoing particles initially along beam axes.
  pAbs = 0.5 * lambda34 / eCM;
  pH[3] = Vec4( 0., 0.,  pAbs, 0.5 * (s + s3 - s4) / eCM);
  pH[4] = Vec4( 0., 0., -pAbs, 0.5 * (s + s4 - s3) / eCM);

  // Then rotate them
  phi = 2. * M_PI * rndmPtr->flat();
  pH[3].rot( theta, phi);
  pH[4].rot( theta, phi);

  // Set some further info for completeness (but Info can be for subprocess).
  x1H   = 1.;
  x2H   = 1.;
  sH    = s;
  uH    = s1 + s2 + s3 + s4 - sH - tH;
  mHat  = eCM;
  p2Abs = pAbs * pAbs;
  betaZ = 0.;
  pTH   = pAbs * sin(theta);

  // Save the sampled photon kinematics.
  if (hasGamma) gammaKinPtr->finalize();

  // Done.
  return true;

}

//==========================================================================

// PhaseSpace2to3diffractive class.
// 2 -> 3 kinematics set up for central diffractive scattering.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of tries to find acceptable (xi1, xi2, t1, t2) set.
const int PhaseSpace2to3diffractive::NTRY = 2500;

// t is sampled according to sum of three distributions.
const double PhaseSpace2to3diffractive::BWID1    = 8.;
const double PhaseSpace2to3diffractive::BWID2    = 4.;
const double PhaseSpace2to3diffractive::BWID3    = 1.;
const double PhaseSpace2to3diffractive::FWID1    = 1.;
const double PhaseSpace2to3diffractive::FWID2    = 0.4;
const double PhaseSpace2to3diffractive::FWID3    = 0.1;

// Safety margins for upper estimate of cross section.
const double PhaseSpace2to3diffractive::MAXFUDGECD = 2.5;
const double PhaseSpace2to3diffractive::MAXFUDGET  = 10.0;

// Safety margin so sum of masses not too close to eCM.
const double PhaseSpace2to3diffractive::DIFFMASSMARGIN = 0.2;

//--------------------------------------------------------------------------

// Set up for phase space sampling.

bool PhaseSpace2to3diffractive::setupSampling() {

  // Find maximum = value of cross section.
  sigmaNw  = sigmaProcessPtr->sigmaHatWrap();
  sigmaMx  = sigmaNw;

  // Squared masses of particles and minimal mass of diffractive states.
  s1       = mA * mA;
  s2       = mB * mB;
  s3       = s1;
  s4       = s2;
  m5min    = sigmaTotPtr->mMinCD();
  s5min    = m5min * m5min;

  // Scenarios with separate handling of xi and t (currently only MBR).
  // Step 0 = both xi and t, 1 = xi only, 2 = t only.
  splitxit = sigmaTotPtr->splitDiff();
  int step = (splitxit) ? 1 : 0;

  // Find maximal cross section xi1 * xi2 * dsigma / (dxi1 dxi2 dt1 dt2)
  // at t1 = t2 = 0 and grid in (xi1, xi2).
  sigMax   = 0.;
  xiMin    = s5min / s;
  for (int i = 0; i < 100; ++i)
  for (int j = 0; j <= i; ++j) {
    xi1    = pow( xiMin, 0.01 * i + 0.005);
    xi2    = pow( xiMin, 0.01 * j + 0.005);
    if (xi1 * xi2 > xiMin) {
      sigNow = sigmaTotPtr->dsigmaCD( xi1, xi2, 0., 0., step);
      if (sigNow > sigMax) sigMax = sigNow;
    }
  }
  sigMax *= MAXFUDGECD;

  // Combinations of t sampling parameters.
  fWid1    = FWID1;
  fWid2    = FWID2;
  fWid3    = FWID3;
  fbWid1   = fWid1 * BWID1;
  fbWid2   = fWid2 * BWID2;
  fbWid3   = fWid3 * BWID3;
  fbWid123 = fbWid1 + fbWid2 + fbWid3;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select a trial kinematics phase space point. Perform full
// Monte Carlo acceptance/rejection at this stage.

bool PhaseSpace2to3diffractive::trialKin( bool, bool ) {

  // Allow for possibility that energy varies from event to event.
  if (doEnergySpread) {
    eCM = infoPtr->eCM();
    s   = eCM * eCM;
  }

  // Trivial kinematics of incoming hadrons.
  pAbs = 0.5 * sqrtpos( pow2( s - s1 - s2) - 4. * s1 * s2 ) / eCM;
  p1.p( 0., 0.,  pAbs, 0.5 * (s + s1 - s2) / eCM);
  p2.p( 0., 0., -pAbs, 0.5 * (s + s2 - s1) / eCM);

  // Temporary variables during generation.
  double pickb, bNow, tNow, sx1, sx2, sx3, sx4, t1, t2, tWeight1, tWeight2,
    lambda12, lambda34, tempA, tempB, tempC, cosTheta, sinTheta, pz, pT,
    deltaE, p2Esum, facP;
  xi1 = xi2 = t1 = t2 = 0.;

  // Normally xi and t in one step, but possible to split into two.
  int nStep = (splitxit) ? 2 : 1;
  for (int iStep = 0; iStep < nStep; ++iStep) {
    int step = (splitxit) ? iStep + 1 : 0;

    // Loop over attempts to set up xi1, xi2, t1, t2 consistently.
    for (int loop = 0; ; ++loop) {
      if (loop == NTRY) {
        infoPtr->errorMsg("Error in PhaseSpace2to3diffractive::trialKin: "
        " quit after repeated tries");
        return false;
      }

      // Select mass^2 = xi1 * xi2 * s according to dxi_1/xi_1 * dxi_2/xi_2.
      if (iStep == 0) {
        do {
          xi1 = pow( s5min / s, rndmPtr->flat());
          xi2 = pow( s5min / s, rndmPtr->flat());
          s5  = xi1 * xi2 * s;
          m5  = sqrt(s5);
        } while (m5 < m5min || mA + mB + m5 + DIFFMASSMARGIN > eCM);
      }

      // Select t1, t2 according to exp(b*t), b picked among three options.
      if (step != 1) {
        bool tryAgain = false;
        for (int i = 0; i < 2; ++i) {
          pickb = rndmPtr->flat() * (fWid1 + fWid2 + fWid3);
          bNow  = (pickb < fWid1) ? BWID1
                : ( (pickb < fWid1 + fWid2) ? BWID2 : BWID3 );
          tNow  = log(rndmPtr->flat()) / bNow;

          // Check whether xi and t choices are consistent on each side.
          sx1   = (i == 0) ? s1 : s2;
          sx2   = (i == 0) ? s2 : s1;
          sx3   = sx1;
          sx4   = (i == 0) ? s2 + xi1 * s : s1 + xi2 * s;
          if (sqrt(sx3) + sqrt(sx4) + DIFFMASSMARGIN > eCM) tryAgain = true;
          if (!tInRange(tNow, s, sx1, sx2, sx3, sx4)) tryAgain = true;
          if (tryAgain) break;
          if (i == 0) t1 = tNow;
          else        t2 = tNow;
        }
        if (tryAgain) continue;
      }

      // Evaluate central diffractive cross section.
      sigNow    = sigmaTotPtr->dsigmaCD( xi1, xi2, t1, t2, step);

      // Maximum weight based on sampling strategy.
      tWeight1  = ( fbWid1 * exp( BWID1 * t1) + fbWid2 * exp(BWID2 * t1)
                + fbWid3 * exp(BWID3 * t1) ) / fbWid123;
      tWeight2  = ( fbWid1 * exp( BWID1 * t2) + fbWid2 * exp(BWID2 * t2)
                + fbWid3 * exp(BWID3 * t2) ) / fbWid123;
      sigMaxNow = (step == 0) ? sigMax * tWeight1 * tWeight2
                : ( (step == 1) ? sigMax : MAXFUDGET * tWeight1 * tWeight2 );

      // Check for maximum violations. Possibly break out of the loop.
      if (sigNow > sigMaxNow) infoPtr->errorMsg("Error in PhaseSpace2to3"
        "diffractive::trialKin: maximum cross section violated");
      if (sigNow > rndmPtr->flat() * sigMaxNow) break;

    // End of loops over tries and steps.
    }
  }

  // Careful reconstruction of scattering angles.
  for (int i = 0; i < 2; ++i) {
    tNow     = (i == 0) ? t1 : t2;
    sx1      = (i == 0) ? s1 : s2;
    sx2      = (i == 0) ? s2 : s1;
    sx3      = sx1;
    sx4      = (i == 0) ? s2 + xi1 * s : s1 + xi2 * s;
    lambda12 = sqrtpos( pow2( s - sx1 - sx2) - 4. * sx1 * sx2 );
    lambda34 = sqrtpos( pow2( s - sx3 - sx4) - 4. * sx3 * sx4 );
    tempA    = s - (sx1 + sx2 + sx3 + sx4) + (sx1 - sx2) * (sx3 - sx4) / s;
    tempB    = lambda12 *  lambda34 / s;
    tempC    = (sx3 - sx1) * (sx4 - sx2) + (sx1 + sx4 - sx2 - sx3)
             * (sx1 * sx4 - sx2 * sx3) / s;
    cosTheta = min(1., max(-1., (tempA + 2. * tNow) / tempB));
    sinTheta = 2. * sqrtpos( -(tempC + tempA * tNow + tNow * tNow) ) / tempB;
    theta    = asin( min(1., sinTheta));
    if (cosTheta < 0.) theta = M_PI - theta;

    // Pick random phi angles. Save outgoing four-vectors.
    pAbs     = 0.5 * lambda34 / eCM;
    pz       = (i == 0) ? pAbs * cos(theta) : -pAbs * cos(theta);
    pT       = pAbs * sin(theta);
    phi      = 2. * M_PI * rndmPtr->flat();
    Vec4& pNow = (i == 0) ? p3 : p4;
    pNow.p( pT * cos(phi), pT * sin(phi), pz, sqrt(pAbs * pAbs + sx1) );
  }

  // Force wanted diffractive mass and rescale three-momenta to fix energies.
  p5     = (p1 - p3) + (p2 - p4);
  p5.e( sqrt(s5 + p5.pAbs2()) );
  for (int iter = 0; iter < 5; ++iter) {
    deltaE = eCM - p3.e() - p4.e() - p5.e();
    if (abs(deltaE) < 1e-10 * eCM) break;
    p2Esum = p3.pAbs2() / p3.e() + p4.pAbs2() / p4.e() + p5.pAbs2() / p5.e();
    facP   = 1. + deltaE / p2Esum;
    p3.rescale3(facP);
    p4.rescale3(facP);
    p5.rescale3(facP);
    p3.e( sqrt(s1 + p3.pAbs2()) );
    p4.e( sqrt(s2 + p4.pAbs2()) );
    p5.e( sqrt(s5 + p5.pAbs2()) );
  }

  return true;
}

//--------------------------------------------------------------------------

// Construct the four-vector kinematics from the trial values.

bool PhaseSpace2to3diffractive::finalKin() {

  // Particle four-momenta and masses.
  pH[1] = p1;
  pH[2] = p2;
  pH[3] = p3;
  pH[4] = p4;
  pH[5] = p5;
  mH[1] = mA;
  mH[2] = mB;
  mH[3] = mA;
  mH[4] = mB;
  mH[5] = m5;

  // Set some further info for completeness (but Info can be for subprocess).
  x1H   = 1.;
  x2H   = 1.;
  sH    = s;
  tH    = (p1 - p3).m2Calc();
  uH    = (p2 - p4).m2Calc();
  mHat  = eCM;
  p2Abs = pAbs * pAbs;
  betaZ = 0.;
  // Store average pT of three final particles for documentation.
  pTH   = (p3.pT() + p4.pT() + p5.pT()) / 3.;

  // Done.
  return true;

}

//==========================================================================

// PhaseSpace2to2nondiffractive class.
// 2 -> 2 kinematics for non-diffractive events.

//--------------------------------------------------------------------------

// Set up for phase space sampling. Trivial if not photoproduction.

bool PhaseSpace2to2nondiffractive::setupSampling(){

  // Flag if a photon inside lepton beam.
  hasGamma = settingsPtr->flag("PDF:lepton2gamma");

  // Default behaviour with usual hadron beams.
  if (!hasGamma) {
    sigmaNw = sigmaProcessPtr->sigmaHat();
    sigmaMx = sigmaNw;

  // Derive overestimate for sigmaND for photons in leptons.
  } else {
    idAgm     = gammaKinPtr->idInA();
    idBgm     = gammaKinPtr->idInB();
    sigmaTotPtr->calc( idAgm, idBgm, eCM);
    sigmaMxGm = sigmaTotPtr->sigmaND();
    sigmaNw   = gammaKinPtr->setupSoftPhaseSpaceSampling(sigmaMxGm);
    sigmaMx   = sigmaNw;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select a trial kinematics phase space point. Trivial if not photoproduction.

bool PhaseSpace2to2nondiffractive::trialKin( bool, bool) {

  // Sample kinematics for gamma+gamma(hadron) sub-event and reject
  // to account for over sampling.
  if (hasGamma) {

    // Current weight.
    double wt = 1.0;

    // Sample gamma kinematics.
    if (!gammaKinPtr->trialKinSoftPhaseSpaceSampling() ) return false;

    // Correct for the estimated sigmaND.
    sigmaTotPtr->calc( idAgm, idBgm, gammaKinPtr->eCMsub() );
    double wtSigma = sigmaTotPtr->sigmaND()/sigmaMxGm;

    // Calculate the total weight and warn if unphysical weight.
    wt *= wtSigma * gammaKinPtr->weight();
    if ( wt > 1. ) infoPtr->errorMsg("Warning in PhaseSpace2to2nondiffractive"
      "::trialKin: weight above unity");

    // Correct for over-estimated cross section and x_gamma limits.
    if ( wt < rndmPtr->flat() ) return false;

  }

  // Done.
  return true;
}

//==========================================================================

// PhaseSpace2to3tauycyl class.
// 2 -> 3 kinematics for normal subprocesses.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of Newton-Raphson iterations of kinematics when masses introduced.
const int PhaseSpace2to3tauycyl::NITERNR = 5;

//--------------------------------------------------------------------------

// Set up for fixed or Breit-Wigner mass selection.

bool PhaseSpace2to3tauycyl::setupMasses() {

  // Treat Z0 as such or as gamma*/Z0
  gmZmode         = gmZmodeGlobal;
  int gmZmodeProc = sigmaProcessPtr->gmZmode();
  if (gmZmodeProc >= 0) gmZmode = gmZmodeProc;

  // Set sHat limits - based on global limits only.
  mHatMin   = mHatGlobalMin;
  sHatMin   = mHatMin*mHatMin;
  mHatMax   = eCM;
  if (mHatGlobalMax > mHatGlobalMin) mHatMax = min( eCM, mHatGlobalMax);
  sHatMax   = mHatMax*mHatMax;

  // Masses and widths of resonances.
  setupMass1(3);
  setupMass1(4);
  setupMass1(5);

  // Reduced mass range - do not make it as fancy as in two-body case.
  if (useBW[3]) mUpper[3] -= (mPeak[4] + mPeak[5]);
  if (useBW[4]) mUpper[4] -= (mPeak[3] + mPeak[5]);
  if (useBW[5]) mUpper[5] -= (mPeak[3] + mPeak[4]);

  // If closed phase space then unallowed process.
  bool physical = true;
  if (useBW[3] && mUpper[3] < mLower[3] + MASSMARGIN) physical = false;
  if (useBW[4] && mUpper[4] < mLower[4] + MASSMARGIN) physical = false;
  if (useBW[5] && mUpper[5] < mLower[5] + MASSMARGIN) physical = false;
  if (!useBW[3] && !useBW[4] && !useBW[5] && mHatMax < mPeak[3]
    + mPeak[4] + mPeak[5] + MASSMARGIN) physical = false;
  if (!physical) return false;

  // No extra pT precautions in massless limit - assumed fixed by ME's.
  pTHatMin  = pTHatGlobalMin;
  pT2HatMin = pTHatMin * pTHatMin;
  pTHatMax  = pTHatGlobalMax;
  pT2HatMax = pTHatMax * pTHatMax;

  // Prepare to select m3 by BW + flat + 1/s_3.
  if (useBW[3]) {
    double distToThreshA = (mHatMax - mPeak[3] - mPeak[4] - mPeak[5])
      * mWidth[3] / (pow2(mWidth[3]) + pow2(mWidth[4]) + pow2(mWidth[5]));
    double distToThreshB = (mHatMax - mPeak[3] - mMin[4] - mMin[5])
      / mWidth[3];
    double distToThresh = min( distToThreshA, distToThreshB);
    setupMass2(3, distToThresh);
  }

  // Prepare to select m4 by BW + flat + 1/s_3.
  if (useBW[4]) {
    double distToThreshA = (mHatMax - mPeak[3] - mPeak[4] - mPeak[5])
      * mWidth[4] / (pow2(mWidth[3]) + pow2(mWidth[4]) + pow2(mWidth[5]));
    double distToThreshB = (mHatMax - mPeak[4] - mMin[3] - mMin[5])
      / mWidth[4];
    double distToThresh = min( distToThreshA, distToThreshB);
    setupMass2(4, distToThresh);
  }

  // Prepare to select m5 by BW + flat + 1/s_3.
  if (useBW[5]) {
    double distToThreshA = (mHatMax - mPeak[3] - mPeak[4] - mPeak[5])
      * mWidth[5] / (pow2(mWidth[3]) + pow2(mWidth[4]) + pow2(mWidth[5]));
    double distToThreshB = (mHatMax - mPeak[5] - mMin[3] - mMin[4])
      / mWidth[5];
    double distToThresh = min( distToThreshA, distToThreshB);
    setupMass2(5, distToThresh);
  }

  // Initialization masses. For now give up when constrained phase space.
  m3 = (useBW[3]) ? min(mPeak[3], mUpper[3]) : mPeak[3];
  m4 = (useBW[4]) ? min(mPeak[4], mUpper[4]) : mPeak[4];
  m5 = (useBW[5]) ? min(mPeak[5], mUpper[5]) : mPeak[5];
  if (m3 + m4 + m5 + MASSMARGIN > mHatMax) physical = false;
  s3 = m3*m3;
  s4 = m4*m4;
  s5 = m5*m5;

  // Correct selected mass-spectrum to running-width Breit-Wigner.
  // Extra safety margin for maximum search.
  wtBW = 1.;
  if (useBW[3]) wtBW *= weightMass(3) * EXTRABWWTMAX;
  if (useBW[4]) wtBW *= weightMass(4) * EXTRABWWTMAX;
  if (useBW[5]) wtBW *= weightMass(5) * EXTRABWWTMAX;

  // Done.
  return physical;

}

//--------------------------------------------------------------------------

// Select Breit-Wigner-distributed or fixed masses.

bool PhaseSpace2to3tauycyl::trialMasses() {

  // By default vanishing cross section.
  sigmaNw = 0.;
  wtBW = 1.;

  // Pick m3, m4 and m5 independently.
  trialMass(3);
  trialMass(4);
  trialMass(5);

  // If outside phase space then reject event.
  if (m3 + m4 + m5 + MASSMARGIN > mHatMax) return false;

  // Correct selected mass-spectrum to running-width Breit-Wigner.
  if (useBW[3]) wtBW *= weightMass(3);
  if (useBW[4]) wtBW *= weightMass(4);
  if (useBW[5]) wtBW *= weightMass(5);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Construct the four-vector kinematics from the trial values.

bool PhaseSpace2to3tauycyl::finalKin() {

  // Assign masses to particles assumed massless in matrix elements.
  int id3 = sigmaProcessPtr->id(3);
  int id4 = sigmaProcessPtr->id(4);
  int id5 = sigmaProcessPtr->id(5);
  if (idMass[3] == 0) { m3 = particleDataPtr->m0(id3); s3 = m3*m3; }
  if (idMass[4] == 0) { m4 = particleDataPtr->m0(id4); s4 = m4*m4; }
  if (idMass[5] == 0) { m5 = particleDataPtr->m0(id5); s5 = m5*m5; }

  // Check that phase space still open after new mass assignment.
  if (m3 + m4 + m5 + MASSMARGIN > mHat) {
    infoPtr->errorMsg("Warning in PhaseSpace2to3tauycyl::finalKin: "
      "failed after mass assignment");
    return false;
  }

  // Particle masses; incoming always on mass shell.
  mH[1] = 0.;
  mH[2] = 0.;
  mH[3] = m3;
  mH[4] = m4;
  mH[5] = m5;

  // Incoming partons along beam axes.
  pH[1] = Vec4( 0., 0.,  0.5 * eCM * x1H, 0.5 * eCM * x1H);
  pH[2] = Vec4( 0., 0., -0.5 * eCM * x2H, 0.5 * eCM * x2H);

  // Begin three-momentum rescaling to compensate for masses.
  if (idMass[3] == 0 || idMass[4] == 0 || idMass[5] == 0) {
    double p3S = p3cm.pAbs2();
    double p4S = p4cm.pAbs2();
    double p5S = p5cm.pAbs2();
    double fac = 1.;
    double e3, e4, e5, value, deriv;

    // Iterate rescaling solution five times, using Newton-Raphson.
    for (int i = 0; i < NITERNR; ++i) {
      e3    = sqrt(s3 + fac * p3S);
      e4    = sqrt(s4 + fac * p4S);
      e5    = sqrt(s5 + fac * p5S);
      value = e3 + e4 + e5 - mHat;
      deriv = 0.5 * (p3S / e3 + p4S / e4 + p5S / e5);
      fac  -= value / deriv;
    }

    // Rescale momenta appropriately.
    double facRoot = sqrt(fac);
    p3cm.rescale3( facRoot );
    p4cm.rescale3( facRoot );
    p5cm.rescale3( facRoot );
    p3cm.e( sqrt(s3 + fac * p3S) );
    p4cm.e( sqrt(s4 + fac * p4S) );
    p5cm.e( sqrt(s5 + fac * p5S) );
  }

  // Outgoing partons initially in collision CM frame along beam axes.
  pH[3] = p3cm;
  pH[4] = p4cm;
  pH[5] = p5cm;

  // Then boost them to overall CM frame
  betaZ = (x1H - x2H)/(x1H + x2H);
  pH[3].rot( theta, phi);
  pH[4].rot( theta, phi);
  pH[3].bst( 0., 0., betaZ);
  pH[4].bst( 0., 0., betaZ);
  pH[5].bst( 0., 0., betaZ);

  // Store average pT of three final particles for documentation.
  pTH = (p3cm.pT() + p4cm.pT() + p5cm.pT()) / 3.;

  // Done.
  return true;
}

//==========================================================================

// The PhaseSpace2to3yyycyl class.
// Phase space for 2 -> 3 QCD processes, 1 + 2 -> 3 + 4 + 5 set up in
// y3, y4, y5, pT2_3, pT2_5, phi_3 and phi_5, and with R separation cut.

//--------------------------------------------------------------------------

//  Sample the phase space of the process.

bool PhaseSpace2to3yyycyl::setupSampling() {

  // Phase space cuts specifically for 2 -> 3 QCD processes.
  pTHat3Min            = settingsPtr->parm("PhaseSpace:pTHat3Min");
  pTHat3Max            = settingsPtr->parm("PhaseSpace:pTHat3Max");
  pTHat5Min            = settingsPtr->parm("PhaseSpace:pTHat5Min");
  pTHat5Max            = settingsPtr->parm("PhaseSpace:pTHat5Max");
  RsepMin              = settingsPtr->parm("PhaseSpace:RsepMin");
  R2sepMin             = pow2(RsepMin);

  // If both beams are baryons then softer PDF's than for mesons/Pomerons.
  hasBaryonBeams = ( beamAPtr->isBaryon() && beamBPtr->isBaryon() );

  // Work with massless partons.
  for (int i = 0; i < 6; ++i) mH[i] = 0.;

  // Constrain to possible cuts at current CM energy and check consistency.
  pT3Min = pTHat3Min;
  pT3Max = pTHat3Max;
  if (pT3Max < pT3Min) pT3Max = 0.5 * eCM;
  pT5Min = pTHat5Min;
  pT5Max = pTHat5Max;
  if (pT5Max < pT5Min) pT5Max = 0.5 * eCM;
  if (pT5Max > pT3Max || pT5Min > pT3Min || pT3Min + 2. * pT5Min > eCM) {
    infoPtr->errorMsg("Error in PhaseSpace2to3yyycyl::setupSampling: "
    "inconsistent pT limits in 3-body phase space");
    return false;
  }

  // Loop over some configurations where cross section could be maximal.
  // In all cases all sum p_z = 0, for maximal PDF weights.
  // Also pT3 and R45 are minimal, while pT5 may vary.
  sigmaMx = 0.;
  double pT5EffMax = min( pT5Max, 0.5 * pT3Min / cos(0.5 * RsepMin) );
  double pT3EffMin = max( pT3Min, 2.0 * pT5Min * cos(0.5 * RsepMin) ) ;
  double sinhR = sinh(0.5 * RsepMin);
  double coshR = cosh(0.5 * RsepMin);
  for (int iStep = 0; iStep < 120; ++iStep) {

    // First kind: |phi4 - phi5| = R, all p_z = 0, i.e. separation in phi.
    if (iStep < 10) {
      pT3   = pT3EffMin;
      pT5   = pT5Min * pow( pT5EffMax / pT5Min, iStep / 9.);
      double pTRat    = pT5 / pT3;
      double sin2Rsep = pow2( sin(RsepMin) );
      double cosPhi35 = - cos(RsepMin) * sqrtpos(1. - sin2Rsep
        * pow2(pTRat)) - sin2Rsep * pTRat;
      cosPhi35 = max( cosPhi35, cos(M_PI - 0.5 * RsepMin) );
      double sinPhi35 = sqrt(1. - pow2(cosPhi35));
      pT4   = sqrt( pow2(pT3) + pow2(pT5) + 2. * pT3 * pT5 * cosPhi35);
      p3cm  = pT3 * Vec4( 1., 0., 0., 1.);
      p4cm  = Vec4(-pT3 - pT5 * cosPhi35, -pT5 * sinPhi35, 0., pT4);
      p5cm  = pT5 * Vec4( cosPhi35, sinPhi35, 0., 1.);

    // Second kind: |y4 - y5| = R, phi4 = phi5, i.e. separation in y.
    } else {
      pT5   = pT5Min * pow( pT5Max / pT5Min, iStep%10 / 9. );
      pT3   = max( pT3Min, 2. * pT5);
      pT4   = pT3 - pT5;
      p4cm  = pT4 * Vec4( -1., 0.,  sinhR, coshR );
      p5cm  = pT5 * Vec4( -1., 0., -sinhR, coshR );
      y3    = -1.2 + 0.2 * (iStep/10);
      p3cm  = pT3 * Vec4( 1., 0., sinh(y3), cosh(y3));
      betaZ = (p3cm.pz() + p4cm.pz() + p5cm.pz())
            / (p3cm.e()  + p4cm.e()  + p5cm.e());
      p3cm.bst( 0., 0., -betaZ);
      p4cm.bst( 0., 0., -betaZ);
      p5cm.bst( 0., 0., -betaZ);
    }

    // Find cross section in chosen phase space point.
    pInSum = p3cm + p4cm + p5cm;
    x1H   = pInSum.e() /  eCM;
    x2H   = x1H;
    sH    = pInSum.m2Calc();
    sigmaProcessPtr->set3Kin( x1H, x2H, sH, p3cm, p4cm, p5cm,
      0., 0., 0., 1., 1., 1.);
    sigmaNw = sigmaProcessPtr->sigmaPDF();

    // Multiply by Jacobian.
    double flux  = 1. /(8. * pow2(sH) * pow5(2. * M_PI));
    double pTRng = pow2(M_PI)
      * pow4(pT3) * (1./pow2(pT3Min) - 1./pow2(pT3Max))
      * pow2(pT5) * 2.* log(pT5Max/pT5Min);
    double yRng  = 8. * log(eCM / pT3) * log(eCM / pT4) * log(eCM / pT5);
    sigmaNw *= SAFETYMARGIN * flux * pTRng * yRng;

    // Update to largest maximum.
    if (showSearch && sigmaNw > sigmaMx) cout << "\n New sigmamax is "
      << scientific << setprecision(3) << sigmaNw << " for x1 = " << x1H
      << " x2 = " << x2H << " sH = " << sH << endl << " p3 = " << p3cm
      << " p4 = " << p4cm << " p5 = " << p5cm;
    if (sigmaNw > sigmaMx) sigmaMx = sigmaNw;
  }
  sigmaPos = sigmaMx;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

//  Sample the phase space of the process.

bool PhaseSpace2to3yyycyl::trialKin(bool inEvent, bool) {

  // Allow for possibility that energy varies from event to event.
  if (doEnergySpread) {
    eCM   = infoPtr->eCM();
    s     = eCM * eCM;
  }
  sigmaNw = 0.;

  // Constrain to possible cuts at current CM energy and check consistency.
  pT3Min = pTHat3Min;
  pT3Max = pTHat3Max;
  if (pT3Max < pT3Min) pT3Max = 0.5 * eCM;
  pT5Min = pTHat5Min;
  pT5Max = pTHat5Max;
  if (pT5Max < pT5Min) pT5Max = 0.5 * eCM;
  if (pT5Max > pT3Max || pT5Min > pT3Min || pT3Min + 2. * pT5Min > eCM) {
    infoPtr->errorMsg("Error in PhaseSpace2to3yyycyl::trialKin: "
    "inconsistent pT limits in 3-body phase space");
    return false;
  }

  // Pick pT3 according to d^2(pT3)/pT3^4 and pT5 to d^2(pT5)/pT5^2.
  pT3    = pT3Min * pT3Max / sqrt( pow2(pT3Min) +
    rndmPtr->flat() * (pow2(pT3Max) - pow2(pT3Min)) );
  pT5Max = min(pT5Max, pT3);
  if (pT5Max < pT5Min) return false;
  pT5    = pT5Min * pow( pT5Max / pT5Min, rndmPtr->flat() );

  // Pick azimuthal angles flat and reconstruct pT4, between pT3 and pT5.
  phi3   = 2. * M_PI * rndmPtr->flat();
  phi5   = 2. * M_PI * rndmPtr->flat();
  pT4    = sqrt( pow2(pT3) + pow2(pT5) + 2. * pT3 * pT5 * cos(phi3 - phi5) );
  if (pT4 > pT3 || pT4 < pT5) return false;
  phi4   = atan2( -(pT3 * sin(phi3) + pT5 * sin(phi5)),
                  -(pT3 * cos(phi3) + pT5 * cos(phi5)) );

  // Pick rapidities flat in allowed ranges.
  y3Max  = log(eCM / pT3);
  y4Max  = log(eCM / pT4);
  y5Max  = log(eCM / pT5);
  y3     = y3Max * (2. * rndmPtr->flat() - 1.);
  y4     = y4Max * (2. * rndmPtr->flat() - 1.);
  y5     = y5Max * (2. * rndmPtr->flat() - 1.);

  // Reject some events at large rapidities to improve efficiency.
  // (Works for baryons, not pions or Pomerons if they have hard PDF's.)
  double WTy = (hasBaryonBeams) ? (1. - pow2(y3/y3Max))
             * (1. - pow2(y4/y4Max)) * (1. - pow2(y5/y5Max)) : 1.;
  if (WTy < rndmPtr->flat()) return false;

  // Check that any pair separated more then RsepMin in (y, phi) space.
  dphi   = abs(phi3 - phi4);
  if (dphi > M_PI) dphi = 2. * M_PI - dphi;
  if (pow2(y3 - y4) + pow2(dphi) < R2sepMin) return false;
  dphi   = abs(phi3 - phi5);
  if (dphi > M_PI) dphi = 2. * M_PI - dphi;
  if (pow2(y3 - y5) + pow2(dphi) < R2sepMin) return false;
  dphi   = abs(phi4 - phi5);
  if (dphi > M_PI) dphi = 2. * M_PI - dphi;
  if (pow2(y4 - y5) + pow2(dphi) < R2sepMin) return false;

  // Reconstruct all four-vectors.
  pH[3]  = pT3 * Vec4( cos(phi3), sin(phi3), sinh(y3), cosh(y3) );
  pH[4]  = pT4 * Vec4( cos(phi4), sin(phi4), sinh(y4), cosh(y4) );
  pH[5]  = pT5 * Vec4( cos(phi5), sin(phi5), sinh(y5), cosh(y5) );
  pInSum = pH[3] + pH[4] + pH[5];

  // Check that x values physical and sHat in allowed range.
  x1H    = (pInSum.e() + pInSum.pz()) /  eCM;
  x2H    = (pInSum.e() - pInSum.pz()) /  eCM;
  if (x1H >= 1. || x2H >= 1.) return false;
  sH     = pInSum.m2Calc();
  if ( sH < pow2(mHatGlobalMin) ||
    (mHatGlobalMax > mHatGlobalMin && sH > pow2(mHatGlobalMax)) )
    return false;

  // Boost four-vectors to rest frame of collision.
  betaZ  = (x1H - x2H)/(x1H + x2H);
  p3cm   = pH[3];    p3cm.bst( 0., 0., -betaZ);
  p4cm   = pH[4];    p4cm.bst( 0., 0., -betaZ);
  p5cm   = pH[5];    p5cm.bst( 0., 0., -betaZ);

  // Find cross section in chosen phase space point.
  sigmaProcessPtr->set3Kin( x1H, x2H, sH, p3cm, p4cm, p5cm,
    0., 0., 0., 1., 1., 1.);
  sigmaNw = sigmaProcessPtr->sigmaPDF();

  // Multiply by Jacobian. Correct for rejection of large rapidities.
  double flux  = 1. /(8. * pow2(sH) * pow5(2. * M_PI));
  double yRng  = 8. * y3Max * y4Max * y5Max;
  double pTRng = pow2(M_PI)
    * pow4(pT3) * (1./pow2(pT3Min) - 1./pow2(pT3Max))
    * pow2(pT5) * 2.* log(pT5Max/pT5Min);
  sigmaNw *= flux * yRng * pTRng / WTy;

  // Allow possibility for user to modify cross section.
  if (canModifySigma) sigmaNw
    *= userHooksPtr->multiplySigmaBy( sigmaProcessPtr, this, inEvent);
  if (canBiasSelection) sigmaNw
    *= userHooksPtr->biasSelectionBy( sigmaProcessPtr, this, inEvent);
  if (canBias2Sel) sigmaNw *= pow( pTH / bias2SelRef, bias2SelPow);

  // Check if maximum violated.
  newSigmaMx = false;
  if (sigmaNw > sigmaMx) {
    infoPtr->errorMsg("Warning in PhaseSpace2to3yyycyl::trialKin: "
      "maximum for cross section violated");

    // Violation strategy 1: increase maximum (always during initialization).
    if (increaseMaximum || !inEvent) {
      double violFact = SAFETYMARGIN * sigmaNw / sigmaMx;
      sigmaMx = SAFETYMARGIN * sigmaNw;
      newSigmaMx = true;
      if (showViolation) {
        if (violFact < 9.99) cout << fixed;
        else                 cout << scientific;
        cout << " PYTHIA Maximum for " << sigmaProcessPtr->name()
             << " increased by factor " << setprecision(3) << violFact
             << " to " << scientific << sigmaMx << endl;
      }

    // Violation strategy 2: weight event (done in ProcessContainer).
    } else if (showViolation && sigmaNw > sigmaPos) {
      double violFact = sigmaNw / sigmaMx;
      if (violFact < 9.99) cout << fixed;
      else                 cout << scientific;
      cout << " PYTHIA Maximum for " << sigmaProcessPtr->name()
           << " exceeded by factor " << setprecision(3) << violFact << endl;
      sigmaPos = sigmaNw;
    }
  }

  // Check if negative cross section.
  if (sigmaNw < sigmaNeg) {
    infoPtr->errorMsg("Warning in PhaseSpace2to3yyycyl::trialKin:"
      " negative cross section set 0", "for " +  sigmaProcessPtr->name() );
    sigmaNeg = sigmaNw;

    // Optional printout of (all) violations.
    if (showViolation) cout << " PYTHIA Negative minimum for "
      << sigmaProcessPtr->name() << " changed to " << scientific
      << setprecision(3) << sigmaNeg << endl;
  }
  if (sigmaNw < 0.) sigmaNw = 0.;


  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Construct the final kinematics of the process: not much left

bool PhaseSpace2to3yyycyl::finalKin() {

  // Work with massless partons.
  for (int i = 0; i < 6; ++i) mH[i] = 0.;

  // Ibncoming partons to collision.
  pH[1] = 0.5 * (pInSum.e() + pInSum.pz()) * Vec4( 0., 0.,  1., 1.);
  pH[2] = 0.5 * (pInSum.e() - pInSum.pz()) * Vec4( 0., 0., -1., 1.);

  // Some quantities meaningless for 2 -> 3. pT defined as average value.
  tH    = 0.;
  uH    = 0.;
  pTH = (pH[3].pT() + pH[4].pT() + pH[5].pT()) / 3.;
  theta = 0.;
  phi   = 0.;

  return true;
}


//==========================================================================

// The PhaseSpaceLHA class.
// A derived class for Les Houches events.
// Note: arbitrary subdivision into PhaseSpaceLHA and SigmaLHAProcess tasks.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// LHA convention with cross section in pb forces conversion to mb.
const double PhaseSpaceLHA::CONVERTPB2MB  = 1e-9;

//--------------------------------------------------------------------------

// Find maximal cross section for comparison with internal processes.

bool PhaseSpaceLHA::setupSampling() {

  // Find which strategy Les Houches events are produced with.
  strategy = lhaUpPtr->strategy();
  stratAbs = abs(strategy);
  if (strategy == 0 || stratAbs > 4) {
    ostringstream stratCode;
    stratCode << strategy;
    infoPtr->errorMsg("Error in PhaseSpaceLHA::setupSampling: unknown "
      "Les Houches Accord weighting stategy", stratCode.str());
    return false;
  }

  // Number of contributing processes.
  nProc = lhaUpPtr->sizeProc();

  // Loop over all processes. Read out maximum and cross section.
  xMaxAbsSum = 0.;
  xSecSgnSum = 0.;
  int    idPr;
  double xMax, xSec, xMaxAbs;
  for (int iProc = 0 ; iProc < nProc; ++iProc) {
    idPr = lhaUpPtr->idProcess(iProc);
    xMax = lhaUpPtr->xMax(iProc);
    xSec = lhaUpPtr->xSec(iProc);

    // Check for inconsistencies between strategy and stored values.
    if ( (strategy == 1 || strategy == 2) && xMax < 0.) {
      infoPtr->errorMsg("Error in PhaseSpaceLHA::setupSampling: "
        "negative maximum not allowed");
      return false;
    }
    if ( ( strategy == 2 || strategy == 3) && xSec < 0.) {
      infoPtr->errorMsg("Error in PhaseSpaceLHA::setupSampling: "
        "negative cross section not allowed");
      return false;
    }

    // Store maximal cross sections for later choice.
    if      (stratAbs == 1) xMaxAbs = abs(xMax);
    else if (stratAbs  < 4) xMaxAbs = abs(xSec);
    else                    xMaxAbs = 1.;
    idProc.push_back( idPr );
    xMaxAbsProc.push_back( xMaxAbs );

    // Find sum and convert to mb.
    xMaxAbsSum += xMaxAbs;
    xSecSgnSum += xSec;
  }
  sigmaMx  = xMaxAbsSum * CONVERTPB2MB;
  sigmaSgn = xSecSgnSum * CONVERTPB2MB;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Construct the next process, by interface to Les Houches class.

bool PhaseSpaceLHA::trialKin( bool, bool repeatSame ) {

  // Must select process type in some cases.
  int idProcNow = 0;
  if (repeatSame) idProcNow = idProcSave;
  else if (stratAbs <= 2) {
    double xMaxAbsRndm = xMaxAbsSum * rndmPtr->flat();
    int iProc = -1;
    do    xMaxAbsRndm -= xMaxAbsProc[++iProc];
    while (xMaxAbsRndm > 0. && iProc < nProc - 1);
    idProcNow = idProc[iProc];
  }

  // Generate Les Houches event. Return if fail (= end of file).
  bool physical = lhaUpPtr->setEvent(idProcNow);
  if (!physical) return false;

  // Find which process was generated.
  int    idPr = lhaUpPtr->idProcess();
  int    iProc = 0;
  for (int iP = 0; iP < int(idProc.size()); ++iP)
    if (idProc[iP] == idPr) iProc = iP;
  idProcSave = idPr;

  // Extract cross section and rescale according to strategy.
  double wtPr = lhaUpPtr->weight();
  if      (stratAbs ==  1) sigmaNw = wtPr * CONVERTPB2MB
    * xMaxAbsSum / xMaxAbsProc[iProc];
  else if (stratAbs ==  2) sigmaNw = (wtPr / abs(lhaUpPtr->xMax(iProc)))
    * sigmaMx;
  else if (strategy ==  3) sigmaNw = sigmaMx;
  else if (strategy == -3 && wtPr > 0.) sigmaNw =  sigmaMx;
  else if (strategy == -3)              sigmaNw = -sigmaMx;
  else if (stratAbs ==  4) sigmaNw = wtPr * CONVERTPB2MB;

  // Set x scales.
  x1H = lhaUpPtr->x1();
  x2H = lhaUpPtr->x2();

  // Done.
  return true;

}

//==========================================================================

} // end namespace Pythia8
