// PhaseSpace.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
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

// MBR: form factor appoximation with two exponents, [FFB1,FFB2] = GeV-2.
const double PhaseSpace::FFA1 = 0.9;
const double PhaseSpace::FFA2 = 0.1;
const double PhaseSpace::FFB1 = 4.6;
const double PhaseSpace::FFB2 = 0.6;

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
  hasOnePointLepton   = hasOneLeptonBeam  && hasPointLepton;
  hasTwoPointLeptons  = hasTwoLeptonBeams && hasPointLepton;

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

  // When to use Breit-Wigners.
  useBreitWigners      = settingsPtr->flag("PhaseSpace:useBreitWigners");
  minWidthBreitWigners = settingsPtr->parm("PhaseSpace:minWidthBreitWigners");

  // Whether generation is with variable energy.
  doEnergySpread       = settingsPtr->flag("Beams:allowMomentumSpread");

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
  for (int i = 1; i < mult; ++i)
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

bool PhaseSpace::setupSampling123(bool is2, bool is3, ostream& os) {

  // Optional printout.
  if (showSearch) os <<  "\n PYTHIA Optimization printout for "
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
  nTau = (hasTwoPointLeptons) ? 1 : 2;
  nY   = (hasOnePointLepton || hasTwoPointLeptons) ? 1 : 5;
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
  if (idResA !=0 && !hasTwoPointLeptons) {
    nTau += 2;
    tauResA = mResA * mResA / s;
    widResA = mResA * GammaResA / s;
  }
  if (idResB != 0 && !hasTwoPointLeptons) {
    nTau += 2;
    tauResB = mResB * mResB / s;
    widResB = mResB * GammaResB / s;
  }

  // More sampling in tau (and different in y) if incoming lepton beams.
  if (hasTwoLeptonBeams && !hasTwoPointLeptons) ++nTau;

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
          sigmaTmp = sigmaProcessPtr->sigmaPDF();
          sigmaTmp *= wtTau * wtY;

        // 2 -> 2: calculate cross section, weighted by phase-space volume
        // and Breit-Wigners for masses
        } else if (is2) {
          sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4,
            runBW3H, runBW4H);
          sigmaTmp = sigmaProcessPtr->sigmaPDF();
          sigmaTmp *= wtTau * wtY * wtZ * wtBW;

        // 2 -> 3: repeat internal 3-body phase space several times and
        // keep maximal cross section, weighted by phase-space volume
        // and Breit-Wigners for masses
        } else if (is3) {
          for (int iTry3 = 0; iTry3 < NTRY3BODY; ++iTry3) {
            if (!select3Body()) continue;
            sigmaProcessPtr->set3Kin( x1H, x2H, sH, p3cm, p4cm, p5cm,
              m3, m4, m5, runBW3H, runBW4H, runBW5H);
            double sigmaTry = sigmaProcessPtr->sigmaPDF();
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
        if (showSearch) os << " tau =" << setw(11) << tau << "  y ="
          << setw(11) << y << "  z =" << setw(11) << z
          << "  sigma =" << setw(11) << sigmaTmp << "\n";
        if (sigmaTmp < 0.) sigmaTmp = 0.;

        // Sum up tau cross-section pieces in points used.
        if (!hasTwoPointLeptons) {
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
        if (!hasOnePointLepton && !hasTwoPointLeptons) {
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
  if (!hasTwoPointLeptons) solveSys( nTau, binTau, vecTau, matTau, tauCoef);
  if (!hasOnePointLepton && !hasTwoPointLeptons)
    solveSys( nY, binY, vecY, matY, yCoef);
  if (is2) solveSys( nZ, binZ, vecZ, matZ, zCoef);
  if (showSearch) os << "\n";

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
          sigmaTmp = sigmaProcessPtr->sigmaPDF();
          sigmaTmp *= wtTau * wtY;

        // 2 -> 2: calculate cross section, weighted by phase-space volume
        // and Breit-Wigners for masses
        } else if (is2) {
          sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4,
            runBW3H, runBW4H);
          sigmaTmp = sigmaProcessPtr->sigmaPDF();
          sigmaTmp *= wtTau * wtY * wtZ * wtBW;

        // 2 -> 3: repeat internal 3-body phase space several times and
        // keep maximal cross section, weighted by phase-space volume
        // and Breit-Wigners for masses
        } else if (is3) {
          for (int iTry3 = 0; iTry3 < NTRY3BODY; ++iTry3) {
            if (!select3Body()) continue;
            sigmaProcessPtr->set3Kin( x1H, x2H, sH, p3cm, p4cm, p5cm,
              m3, m4, m5, runBW3H, runBW4H, runBW5H);
            double sigmaTry = sigmaProcessPtr->sigmaPDF();
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
        if (showSearch) os << " tau =" << setw(11) << tau << "  y ="
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
  if (showSearch) os << "\n";

  // Read out starting position for search.
  sigmaMx = sigMax[0];
  int beginVar = (hasTwoPointLeptons) ? 2 : 0;
  if (hasOnePointLepton) beginVar = 1;
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
              sigmaTmp = sigmaProcessPtr->sigmaPDF();
              sigmaTmp *= wtTau * wtY;

            // 2 -> 2: calculate cross section, weighted by phase-space volume
            // and Breit-Wigners for masses
            } else if (is2) {
              sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4,
                runBW3H, runBW4H);
              sigmaTmp = sigmaProcessPtr->sigmaPDF();
              sigmaTmp *= wtTau * wtY * wtZ * wtBW;

            // 2 -> 3: repeat internal 3-body phase space several times and
            // keep maximal cross section, weighted by phase-space volume
            // and Breit-Wigners for masses
            } else if (is3) {
              for (int iTry3 = 0; iTry3 < NTRY3BODY; ++iTry3) {
                if (!select3Body()) continue;
                sigmaProcessPtr->set3Kin( x1H, x2H, sH, p3cm, p4cm, p5cm,
                  m3, m4, m5, runBW3H, runBW4H, runBW5H);
                double sigmaTry = sigmaProcessPtr->sigmaPDF();
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
            if (showSearch) os << " tau =" << setw(11) << tau << "  y ="
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
  if (showSearch) os << "\n Final maximum = "  << setw(11) << sigmaMx << endl;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Select a trial kinematics phase space point.
// Note: by In is meant the integral over the quantity multiplying
// coefficient cn. The sum of cn is normalized to unity.

bool PhaseSpace::trialKin123(bool is2, bool is3, bool inEvent, ostream& os) {

  // Allow for possibility that energy varies from event to event.
  if (doEnergySpread) {
    eCM       = infoPtr->eCM();
    s         = eCM * eCM;

    // Find shifted tauRes values.
    if (idResA !=0 && !hasTwoPointLeptons) {
      tauResA = mResA * mResA / s;
      widResA = mResA * GammaResA / s;
    }
    if (idResB != 0 && !hasTwoPointLeptons) {
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
  if (!hasTwoPointLeptons) {
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
  if (!hasOnePointLepton && !hasTwoPointLeptons) {
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
        if (violFact < 9.99) os << fixed;
        else                 os << scientific;
        os << " PYTHIA Maximum for " << sigmaProcessPtr->name()
           << " increased by factor " << setprecision(3) << violFact
           << " to " << scientific << sigmaMx << endl;
      }

    // Violation strategy 2: weight event (done in ProcessContainer).
    } else if (showViolation && sigmaNw > sigmaPos) {
      double violFact = sigmaNw / sigmaMx;
      if (violFact < 9.99) os << fixed;
      else                 os << scientific;
      os << " PYTHIA Maximum for " << sigmaProcessPtr->name()
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
    if (showViolation) os << " PYTHIA Negative minimum for "
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
  if (hasTwoPointLeptons) {
    tauMin = 1.;
    tauMax = 1.;
    return true;
  }

  // Requirements from allowed mHat range.
  tauMin = sHatMin / s;
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
  if (hasTwoPointLeptons) {
    yMax = 1.;
    return true;
  }

  // Requirements from selected tau value. Trivial for one unresolved beam.
  yMax = -0.5 * log(tau);
  if (hasOnePointLepton) return true;

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

  // Check that there is an open range.
  return (zMax > zMin);
}

//--------------------------------------------------------------------------

// Select tau according to a choice of shapes.

void PhaseSpace::selectTau(int iTau, double tauVal, bool is2) {

  // Trivial reply for unresolved lepton beams.
  if (hasTwoPointLeptons) {
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
  if (hasTwoPointLeptons) {
    y = 0.;
    wtY = 1.;
    x1H = 1.;
    x2H = 1.;
    return;
  }

  // Trivial replies for one unresolved lepton beam.
  if (hasOnePointLepton) {
    if (hasLeptonBeamA) {
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

void PhaseSpace::selectZ(int iZ, double zVal) {

  // Mass-dependent dampening of pT -> 0 limit.
  ratio34 = max(TINY, 2. * s3 * s4 / pow2(sH));
  unity34 = 1. + ratio34;
  double ratiopT2 = 2. * pT2HatMin / max( SHATMINZ, sH);
  if (ratiopT2 < PT2RATMINZ) ratio34 = max( ratio34, ratiopT2);

  // Common expressions in z limits.
  double zPosMax = max(ratio34, unity34 + zMax);
  double zNegMax = max(ratio34, unity34 - zMax);
  double zPosMin = max(ratio34, unity34 + zMin);
  double zNegMin = max(ratio34, unity34 - zMin);

  // Flat in z.
  if (iZ == 0) {
    if (zVal < 0.5) z = -(zMax + (zMin - zMax) * 2. * zVal);
    else z = zMin + (zMax - zMin) * (2. * zVal - 1.);

  // 1 / (unity34 - z).
  } else if (iZ == 1) {
    double areaNeg = log(zPosMax / zPosMin);
    double areaPos = log(zNegMin / zNegMax);
    double area = areaNeg + areaPos;
    if (zVal * area < areaNeg) {
      double zValMod = zVal * area / areaNeg;
      z = unity34 - zPosMax * pow(zPosMin / zPosMax, zValMod);
    } else {
      double zValMod = (zVal * area - areaNeg)/ areaPos;
      z = unity34 - zNegMin * pow(zNegMax / zNegMin, zValMod);
    }

  // 1 / (unity34 + z).
  } else if (iZ == 2) {
    double areaNeg = log(zNegMin / zNegMax);
    double areaPos = log(zPosMax / zPosMin);
    double area = areaNeg + areaPos;
    if (zVal * area < areaNeg) {
      double zValMod = zVal * area / areaNeg;
      z = zNegMax * pow(zNegMin / zNegMax, zValMod) - unity34;
    } else {
      double zValMod = (zVal * area - areaNeg)/ areaPos;
      z = zPosMin * pow(zPosMax / zPosMin, zValMod) - unity34;
    }

  // 1 / (unity34 - z)^2.
  } else if (iZ == 3) {
    double areaNeg = 1. / zPosMin - 1. / zPosMax;
    double areaPos = 1. / zNegMax - 1. / zNegMin;
    double area = areaNeg + areaPos;
    if (zVal * area < areaNeg) {
      double zValMod = zVal * area / areaNeg;
      z = unity34 - 1. / (1./zPosMax + areaNeg * zValMod);
    } else {
      double zValMod = (zVal * area - areaNeg)/ areaPos;
      z = unity34 - 1. / (1./zNegMin + areaPos * zValMod);
    }

  // 1 / (unity34 + z)^2.
  } else if (iZ == 4) {
    double areaNeg = 1. / zNegMax - 1. / zNegMin;
    double areaPos = 1. / zPosMin - 1. / zPosMax;
    double area = areaNeg + areaPos;
    if (zVal * area < areaNeg) {
      double zValMod = zVal * area / areaNeg;
      z = 1. / (1./zNegMax - areaNeg * zValMod) - unity34;
    } else {
      double zValMod = (zVal * area - areaNeg)/ areaPos;
      z = 1. / (1./zPosMin - areaPos * zValMod) - unity34;
    }
  }

  // Safety check for roundoff errors. Combinations with z.
  if (z < 0.) z = min(-zMin, max(-zMax, z));
  else z = min(zMax, max(zMin, z));
  zNeg = max(ratio34, unity34 - z);
  zPos = max(ratio34, unity34 + z);

  // Phase space integral in z.
  double intZ0 = 2. * (zMax - zMin);
  double intZ12 = log( (zPosMax * zNegMin) / (zPosMin * zNegMax) );
  double intZ34 = 1. / zPosMin - 1. / zPosMax + 1. / zNegMax
    - 1. / zNegMin;
  wtZ = mHat * pAbs / ( (zCoef[0] / intZ0) + (zCoef[1] / intZ12) / zNeg
    + (zCoef[2] / intZ12) / zPos + (zCoef[3] / intZ34) / pow2(zNeg)
    + (zCoef[4] / intZ34) / pow2(zPos) );

  // Calculate tHat and uHat. Also gives pTHat.
  double sH34 = -0.5 * (sH - s3 - s4);
  tH  = sH34 + mHat * pAbs * z;
  uH  = sH34 - mHat * pAbs * z;
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
  double mat[8][8], double coef[8], ostream& os) {

  // Optional printout.
  if (showSearch) {
    os << "\n Equation system: " << setw(5) << bin[0];
    for (int j = 0; j < n; ++j) os << setw(12) << mat[0][j];
    os << setw(12) << vec[0] << "\n";
    for (int i = 1; i < n; ++i) {
      os << "                  " << setw(5) << bin[i];
      for (int j = 0; j < n; ++j) os << setw(12) << mat[i][j];
      os << setw(12) << vec[i] << "\n";
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
    os << " Solution:             ";
    for (int i = 0; i < n; ++i) os << setw(12) << coef[i];
    os << "\n";
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
    mMin[iM]   = particleDataPtr->mMin(idMass[iM]);
    mMax[iM]   = particleDataPtr->mMax(idMass[iM]);
    // gmZmode == 1 means pure photon propagator; set at lower mass limit.
    if (idMass[iM] == 23 && gmZmode == 1) mPeak[iM] = mMin[iM];
  }

  // Mass and width combinations for Breit-Wigners.
  sPeak[iM]    = mPeak[iM] * mPeak[iM];
  useBW[iM]    = useBreitWigners && (mWidth[iM] > minWidthBreitWigners);
  if (!useBW[iM]) mWidth[iM] = 0.;
  mw[iM]       = mPeak[iM] * mWidth[iM];
  wmRat[iM]    = (idMass[iM] == 0 || mPeak[iM] == 0.)
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
    fracFlat[iM] = 0.1;
    fracInv[iM]  = 0.1;
  } else if (distToThresh > - THRESHOLDSIZE) {
    fracFlat[iM] = 0.25 - 0.15 * distToThresh / THRESHOLDSIZE;
    fracInv [iM] = 0.15 - 0.05 * distToThresh / THRESHOLDSIZE;
  } else {
   fracFlat[iM]  = 0.4;
   fracInv[iM]   = 0.2;
  }

  // For gamma*/Z0: increase 1/s_i part and introduce 1/s_i^2 part.
  fracInv2[iM]   = 0.;
  if (idMass[iM] == 23 && gmZmode == 0) {
    fracFlat[iM] *= 0.5;
    fracInv[iM]  = 0.5 * fracInv[iM] + 0.25;
    fracInv2[iM] = 0.25;
  } else if (idMass[iM] == 23 && gmZmode == 1) {
    fracFlat[iM] = 0.1;
    fracInv[iM]  = 0.4;
    fracInv2[iM] = 0.4;
  }

  // Normalization integrals for the respective contribution.
  atanLower[iM]  = atan( (sLower[iM] - sPeak[iM])/ mw[iM] );
  atanUpper[iM]  = atan( (sUpper[iM] - sPeak[iM])/ mw[iM] );
  intBW[iM]      = atanUpper[iM] - atanLower[iM];
  intFlat[iM]    = sUpper[iM] - sLower[iM];
  intInv[iM]     = log( sUpper[iM] / sLower[iM] );
  intInv2[iM]    = 1./sLower[iM] - 1./sUpper[iM];

}

//--------------------------------------------------------------------------

// Select Breit-Wigner-distributed or fixed masses.

void PhaseSpace::trialMass(int iM) {

  // References to masses to be set.
  double& mSet = (iM == 3) ? m3 : ( (iM == 4) ? m4 : m5 );
  double& sSet = (iM == 3) ? s3 : ( (iM == 4) ? s4 : s5 );

  // Distribution for m_i is BW + flat + 1/s_i + 1/s_i^2.
  if (useBW[iM]) {
    double pickForm = rndmPtr->flat();
    if (pickForm > fracFlat[iM] + fracInv[iM] + fracInv2[iM])
      sSet = sPeak[iM] + mw[iM] * tan( atanLower[iM]
           + rndmPtr->flat() * intBW[iM] );
    else if (pickForm > fracInv[iM] + fracInv2[iM])
      sSet = sLower[iM] + rndmPtr->flat() * (sUpper[iM] - sLower[iM]);
    else if (pickForm > fracInv2[iM])
      sSet = sLower[iM] * pow( sUpper[iM] / sLower[iM], rndmPtr->flat() );
    else sSet = sLower[iM] * sUpper[iM]
      / (sLower[iM] + rndmPtr->flat() * (sUpper[iM] - sLower[iM]));
    mSet = sqrt(sSet);

  // Else m_i is fixed at peak value.
  } else {
    mSet = mPeak[iM];
    sSet = sPeak[iM];
  }

}

//--------------------------------------------------------------------------

// Naively a fixed-width Breit-Wigner is used to pick the mass.
// Here come the correction factors for
// (i) preselection according to BW + flat in s_i + 1/s_i + 1/s_i^2,
// (ii) reduced allowed mass range,
// (iii) running width, i.e. m0*Gamma0 -> s*Gamma0/m0.
// In the end, the weighted distribution is a running-width BW.

double PhaseSpace::weightMass(int iM) {

  // Reference to mass and to Breit-Wigner weight to be set.
  double& sSet   = (iM == 3) ? s3 : ( (iM == 4) ? s4 : s5 );
  double& runBWH = (iM == 3) ? runBW3H : ( (iM == 4) ? runBW4H : runBW5H );

  // Default weight if no Breit-Wigner.
  runBWH = 1.;
  if (!useBW[iM]) return 1.;

  // Weight of generated distribution.
  double genBW  = (1. - fracFlat[iM] - fracInv[iM] - fracInv2[iM])
      * mw[iM] / ( (pow2(sSet - sPeak[iM]) + pow2(mw[iM])) * intBW[iM])
      + fracFlat[iM] / intFlat[iM] + fracInv[iM] / (sSet * intInv[iM])
      + fracInv2[iM] / (sSet*sSet * intInv2[iM]);

  // Weight of distribution with running width in Breit-Wigner.
  double mwRun = sSet * wmRat[iM];
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
  pH[1] = Vec4( 0., 0., 0.5 * eCM * x1H, 0.5 * eCM * x1H);
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

  // Incoming partons along beam axes.
  pH[1] = Vec4( 0., 0., 0.5 * eCM * x1H, 0.5 * eCM * x1H);
  pH[2] = Vec4( 0., 0., -0.5 * eCM * x2H, 0.5 * eCM * x2H);

  // Outgoing partons initially in collision CM frame along beam axes.
  pH[3] = Vec4( 0., 0.,  pAbs, 0.5 * (sH + s3 - s4) / mHat);
  pH[4] = Vec4( 0., 0., -pAbs, 0.5 * (sH + s4 - s3) / mHat);

  // Then rotate and boost them to overall CM frame.
  theta = acos(z);
  phi = 2. * M_PI * rndmPtr->flat();
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

//==========================================================================

// PhaseSpace2to2elastic class.
// 2 -> 2 kinematics set up for elastic scattering.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum positive/negative argument for exponentiation.
const double PhaseSpace2to2elastic::EXPMAX = 50.;

// Conversion coefficients = 1/(16pi) * (mb <-> GeV^2).
const double PhaseSpace2to2elastic::CONVERTEL = 0.0510925;

//--------------------------------------------------------------------------

// Form of phase space sampling already fixed, so no optimization.
// However, need to read out relevant parameters from SigmaTotal.

bool PhaseSpace2to2elastic::setupSampling() {

  // Find maximum = value of cross section.
  sigmaNw    = sigmaProcessPtr->sigmaHatWrap();
  sigmaMx    = sigmaNw;

  // Squared and outgoing masses of particles.
  s1         = mA * mA;
  s2         = mB * mB;
  m3         = mA;
  m4         = mB;

  // Elastic slope.
  bSlope     = sigmaTotPtr->bSlopeEl();

  // Determine maximum possible t range.
  lambda12S  = pow2(s - s1 - s2) - 4. * s1 * s2 ;
  tLow       = - lambda12S / s;
  tUpp       = 0.;

  // Production model with Coulomb corrections need more parameters.
  useCoulomb =  settingsPtr->flag("SigmaTotal:setOwn")
             && settingsPtr->flag("SigmaElastic:setOwn");
  if (useCoulomb) {
    sigmaTot = sigmaTotPtr->sigmaTot();
    rho      = settingsPtr->parm("SigmaElastic:rho");
    lambda   = settingsPtr->parm("SigmaElastic:lambda");
    tAbsMin  = settingsPtr->parm("SigmaElastic:tAbsMin");
    phaseCst = settingsPtr->parm("SigmaElastic:phaseConst");
    alphaEM0 = settingsPtr->parm("StandardModel:alphaEM0");

    // Relative rate of nuclear and Coulombic parts in trials.
    tUpp     = -tAbsMin;
    sigmaNuc = CONVERTEL * pow2(sigmaTot) * (1. + rho*rho) / bSlope
             * exp(-bSlope * tAbsMin);
    sigmaCou = (useCoulomb) ?
               pow2(alphaEM0) / (4. * CONVERTEL * tAbsMin) : 0.;
    signCou  = (idA == idB) ? 1. : -1.;

  // Dummy values.
  } else {
    sigmaNuc = sigmaNw;
    sigmaCou = 0.;
  }

  // Calculate coefficient of generation.
  tAux       = exp( max(-EXPMAX, bSlope * (tLow - tUpp)) ) - 1.;

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
    lambda12S = pow2(s - s1 - s2) - 4. * s1 * s2 ;
    tLow      = - lambda12S / s;
    tAux      = exp( max(-EXPMAX, bSlope * (tLow - tUpp)) ) - 1.;
  }

  // Select t according to exp(bSlope*t).
  if (!useCoulomb || sigmaNuc > rndmPtr->flat() * (sigmaNuc + sigmaCou))
   tH = tUpp + log(1. + tAux * rndmPtr->flat()) / bSlope;

 // Select t according to 1/t^2.
  else tH = tLow * tUpp / (tUpp + rndmPtr->flat() * (tLow - tUpp));

  // Correction factor for ratio full/simulated.
  if (useCoulomb) {
    double sigmaN   = CONVERTEL * pow2(sigmaTot) * (1. + rho*rho)
                    * exp(bSlope * tH);
    double alpEM    = couplingsPtr->alphaEM(-tH);
    double sigmaC   = pow2(alpEM) / (4. * CONVERTEL * tH*tH);
    double sigmaGen = 2. * (sigmaN + sigmaC);
    double form2    = pow4(lambda/(lambda - tH));
    double phase    = signCou * alpEM
                    * (-phaseCst - log(-0.5 * bSlope * tH));
    double sigmaCor = sigmaN + pow2(form2) * sigmaC
      - signCou * alpEM * sigmaTot * (form2 / (-tH))
      *  exp(0.5 * bSlope * tH) * (rho * cos(phase) + sin(phase));
    sigmaNw         = sigmaMx * sigmaCor / sigmaGen;
  }

  // Careful reconstruction of scattering angle.
  double tRat       = s * tH / lambda12S;
  double cosTheta   = min(1., max(-1., 1. + 2. * tRat ) );
  double sinTheta   = 2. * sqrtpos( -tRat * (1. + tRat) );
  theta             = asin( min(1., sinTheta));
  if (cosTheta < 0.) theta = M_PI - theta;

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

  // Incoming particles along beam axes.
  pAbs = 0.5 * sqrtpos(lambda12S) / eCM;
  pH[1] = Vec4( 0., 0.,  pAbs, 0.5 * (s + s1 - s2) / eCM);
  pH[2] = Vec4( 0., 0., -pAbs, 0.5 * (s + s2 - s1) / eCM);

  // Outgoing particles initially along beam axes.
  pH[3] = Vec4( 0., 0.,  pAbs, 0.5 * (s + s1 - s2) / eCM);
  pH[4] = Vec4( 0., 0., -pAbs, 0.5 * (s + s2 - s1) / eCM);

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
const int PhaseSpace2to2diffractive::NTRY = 500;

// Maximum positive/negative argument for exponentiation.
const double PhaseSpace2to2diffractive::EXPMAX = 50.;

// Safety margin so sum of diffractive masses not too close to eCM.
const double PhaseSpace2to2diffractive::DIFFMASSMARGIN = 0.2;

//--------------------------------------------------------------------------

// Form of phase space sampling already fixed, so no optimization.
// However, need to read out relevant parameters from SigmaTotal.

bool PhaseSpace2to2diffractive::setupSampling() {

  // Pomeron flux parametrization, and parameters of some options.
  PomFlux      = settingsPtr->mode("Diffraction:PomFlux");
  epsilonPF    = settingsPtr->parm("Diffraction:PomFluxEpsilon");
  alphaPrimePF = settingsPtr->parm("Diffraction:PomFluxAlphaPrime");

  // Find maximum = value of cross section.
  sigmaNw = sigmaProcessPtr->sigmaHatWrap();
  sigmaMx = sigmaNw;

  // Masses of particles and minimal masses of diffractive states.
  m3ElDiff = (isDiffA) ? sigmaTotPtr->mMinXB()  : mA;
  m4ElDiff = (isDiffB) ? sigmaTotPtr->mMinAX()  : mB;
  s1 = mA * mA;
  s2 = mB * mB;
  s3 = pow2( m3ElDiff);
  s4 = pow2( m4ElDiff);

  // Determine maximum possible t range and coefficient of generation.
  lambda12 = sqrtpos( pow2( s - s1 - s2) - 4. * s1 * s2 );
  lambda34 = sqrtpos( pow2( s - s3 - s4) - 4. * s3 * s4 );
  double tempA = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4) / s;
  double tempB = lambda12 *  lambda34 / s;
  double tempC = (s3 - s1) * (s4 - s2) + (s1 + s4 - s2 - s3)
    * (s1 * s4 - s2 * s3) / s;
  tLow  = -0.5 * (tempA + tempB);
  tUpp  = tempC / tLow;

  // Default for all parametrization-specific parameters.
  cRes = sResXB = sResAX = sProton = bMin = bSlope = bSlope1 = bSlope2
       = probSlope1 = xIntPF = xtCorPF = mp24DL = coefDL = tAux
       = tAux1 = tAux2 = 0.;

  // Schuler&Sjostrand: parameters of low-mass-resonance enhancement.
  if (PomFlux == 1) {
    cRes = sigmaTotPtr->cRes();
    sResXB = pow2( sigmaTotPtr->mResXB());
    sResAX = pow2( sigmaTotPtr->mResAX());
    sProton = sigmaTotPtr->sProton();

    // Schuler&Sjostrand: lower limit diffractive slope.
    if      (!isDiffB) bMin = sigmaTotPtr->bMinSlopeXB();
    else if (!isDiffA) bMin = sigmaTotPtr->bMinSlopeAX();
    else               bMin = sigmaTotPtr->bMinSlopeXX();
    tAux = exp( max(-EXPMAX, bMin * (tLow - tUpp)) ) - 1.;

  // Bruni&Ingelman: relative weight of two diffractive slopes.
  } else if (PomFlux == 2) {
    bSlope1     = 8.0;
    probSlope1  = 6.38 * ( exp(max(-EXPMAX, bSlope1 * tUpp))
                -  exp(max(-EXPMAX, bSlope1 * tLow)) ) / bSlope1;
    bSlope2     = 3.0;
    double pS2  = 0.424 * ( exp(max(-EXPMAX, bSlope2 * tUpp))
                -  exp(max(-EXPMAX, bSlope2 * tLow)) ) / bSlope2;
    probSlope1 /= probSlope1 + pS2;
    tAux1 = exp( max(-EXPMAX, bSlope1 * (tLow - tUpp)) ) - 1.;
    tAux2 = exp( max(-EXPMAX, bSlope2 * (tLow - tUpp)) ) - 1.;

  // Streng&Berger (RapGap): diffractive slope, power of mass spectrum.
  } else if (PomFlux == 3) {
    bSlope        = 4.7;
    double xPowPF = 1. - 2. * (1. + epsilonPF);
    xIntPF        = 2. * (1. + xPowPF);
    xtCorPF       = 2. * alphaPrimePF;
    tAux          = exp( max(-EXPMAX, bSlope  * (tLow - tUpp)) ) - 1.;

  // Donnachie&Landshoff (RapGap):  power of mass spectrum.
  } else if (PomFlux == 4) {
    mp24DL        = 4. * pow2(particleDataPtr->m0(2212));
    double xPowPF = 1. - 2. * (1. + epsilonPF);
    xIntPF        = 2. * (1. + xPowPF);
    xtCorPF       = 2. * alphaPrimePF;
    // Upper estimate of t dependence, for preliminary choice.
    coefDL               = 0.85;
    tAux1                = 1. / pow3(1. - coefDL * tLow);
    tAux2                = 1. / pow3(1. - coefDL * tUpp);

  // MBR model.
  } else if (PomFlux == 5) {
    eps        = settingsPtr->parm("Diffraction:MBRepsilon");
    alph       = settingsPtr->parm("Diffraction:MBRalpha");
    alph2      = alph * alph;
    m2min      = settingsPtr->parm("Diffraction:MBRm2Min");
    dyminSD    = settingsPtr->parm("Diffraction:MBRdyminSD");
    dyminDD    = settingsPtr->parm("Diffraction:MBRdyminDD");
    dyminSigSD = settingsPtr->parm("Diffraction:MBRdyminSigSD");
    dyminSigDD = settingsPtr->parm("Diffraction:MBRdyminSigDD");

    // Max f(dy) for Von Neumann method, from SigmaTot.
    sdpmax= sigmaTotPtr->sdpMax();
    ddpmax= sigmaTotPtr->ddpMax();

  // H1 Fit A/B.
  } else if (PomFlux == 6 || PomFlux == 7) {
    bSlope        = 5.5;
    epsilonPF     =  (PomFlux == 6) ? 0.1182 : 0.1110;
    alphaPrimePF  = 0.06;
    double xPowPF = 1. - 2. * (1. + epsilonPF);
    xIntPF        = 2. * (1. + xPowPF);
    xtCorPF       = 2. * alphaPrimePF;
    tAux          = exp( max(-EXPMAX, bSlope  * (tLow - tUpp)) ) - 1.;
  }

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
    lambda34 = sqrtpos( pow2( s - s3 - s4) - 4. * s3 * s4 );
    double tempA = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4) / s;
    double tempB = lambda12 *  lambda34 / s;
    double tempC = (s3 - s1) * (s4 - s2) + (s1 + s4 - s2 - s3)
      * (s1 * s4 - s2 * s3) / s;
    tLow  = -0.5 * (tempA + tempB);
    tUpp  = tempC / tLow;
    if (PomFlux == 1) {
      tAux = exp( max(-EXPMAX, bMin * (tLow - tUpp)) ) - 1.;
    } else if (PomFlux == 2) {
      tAux1 = exp( max(-EXPMAX, bSlope1 * (tLow - tUpp)) ) - 1.;
      tAux2 = exp( max(-EXPMAX, bSlope2 * (tLow - tUpp)) ) - 1.;
    } else if (PomFlux == 3 || PomFlux == 6 || PomFlux == 7) {
      tAux          = exp( max(-EXPMAX, bSlope  * (tLow - tUpp)) ) - 1.;
    } else if (PomFlux == 4) {
      tAux1                = 1. / pow3(1. - coefDL * tLow);
      tAux2                = 1. / pow3(1. - coefDL * tUpp);
    }
  }

  // Loop over attempts to set up masses and t consistently.
  for (int loop = 0; ; ++loop) {
    if (loop == NTRY) {
      infoPtr->errorMsg("Error in PhaseSpace2to2diffractive::trialKin: "
        " quit after repeated tries");
      return false;
    }

    // Schuler and Sjostrand:
    if (PomFlux == 1) {

      // Select diffractive mass(es) according to dm^2/m^2.
      m3 = (isDiffA) ? m3ElDiff * pow( max(mA, eCM - m4ElDiff) / m3ElDiff,
        rndmPtr->flat()) : m3ElDiff;
      m4 = (isDiffB) ? m4ElDiff * pow( max(mB, eCM - m3ElDiff) / m4ElDiff,
        rndmPtr->flat()) : m4ElDiff;
      if (m3 + m4 + DIFFMASSMARGIN >= eCM) continue;
      s3 = m3 * m3;
      s4 = m4 * m4;

      // Additional mass factors, including resonance enhancement.
      if (isDiffA && !isDiffB) {
        double facXB = (1. - s3 / s)
          * (1. + cRes * sResXB / (sResXB + s3));
        if (facXB < rndmPtr->flat() * (1. + cRes)) continue;
      } else if (isDiffB && !isDiffA) {
        double facAX = (1. - s4 / s)
          * (1. + cRes * sResAX / (sResAX + s4));
        if (facAX < rndmPtr->flat() * (1. + cRes)) continue;
      } else {
        double facXX = (1. - pow2(m3 + m4) / s)
          * (s * sProton / (s * sProton + s3 * s4))
          * (1. + cRes * sResXB / (sResXB + s3))
          * (1. + cRes * sResAX / (sResAX + s4));
        if (facXX < rndmPtr->flat() * pow2(1. + cRes)) continue;
      }

      // Select t according to exp(bMin*t) and correct to right slope.
      tH = tUpp + log(1. + tAux * rndmPtr->flat()) / bMin;
      double bDiff = 0.;
      if (isDiffA && !isDiffB) bDiff = sigmaTotPtr->bSlopeXB(s3) - bMin;
      else if (!isDiffA) bDiff = sigmaTotPtr->bSlopeAX(s4) - bMin;
      else bDiff = sigmaTotPtr->bSlopeXX(s3, s4) - bMin;
      bDiff = max(0., bDiff);
      if (exp( max(-EXPMAX, bDiff * (tH - tUpp)) ) < rndmPtr->flat()) continue;

    // Bruni and Ingelman:
    } else if (PomFlux == 2) {

      // Select diffractive mass(es) according to dm^2/m^2.
      m3 = (isDiffA) ? m3ElDiff * pow( max(mA, eCM - m4ElDiff) / m3ElDiff,
        rndmPtr->flat()) : m3ElDiff;
      m4 = (isDiffB) ? m4ElDiff * pow( max(mB, eCM - m3ElDiff) / m4ElDiff,
        rndmPtr->flat()) : m4ElDiff;
      if (m3 + m4 + DIFFMASSMARGIN >= eCM) continue;
      s3 = m3 * m3;
      s4 = m4 * m4;

      // Select t according to exp(bSlope*t) with two possible slopes.
      tH = (rndmPtr->flat() < probSlope1)
         ? tUpp + log(1. + tAux1 * rndmPtr->flat()) / bSlope1
         : tUpp + log(1. + tAux2 * rndmPtr->flat()) / bSlope2;

    // Streng and Berger et al. (RapGap) & H1 Fit A/B:
    } else if (PomFlux == 3 || PomFlux == 6 || PomFlux == 7) {

      // Select diffractive mass(es) according to dm^2/(m^2)^(1 + 2 epsilon).
      m3 = m3ElDiff;
      m4 = m4ElDiff;
      if (isDiffA) {
        double s3MinPow = pow( m3ElDiff, xIntPF );
        double s3MaxPow = pow( max(mA, eCM - m4ElDiff), xIntPF );
        m3 = pow( s3MinPow + rndmPtr->flat() * (s3MaxPow - s3MinPow),
                  1. / xIntPF );
      }
      if (isDiffB) {
        double s4MinPow = pow( m4ElDiff, xIntPF );
        double s4MaxPow = pow( max(mB, eCM - m3ElDiff), xIntPF );
        m4 = pow( s4MinPow + rndmPtr->flat() * (s4MaxPow - s4MinPow),
                  1. / xIntPF );
      }
      if (m3 + m4 + DIFFMASSMARGIN >= eCM) continue;
      s3 = m3 * m3;
      s4 = m4 * m4;

      // Select t according to exponential and weigh by x_P^(2 alpha' |t|).
      tH = tUpp + log(1. + tAux * rndmPtr->flat()) / bSlope;
      if ( isDiffA && pow( s3 / s, xtCorPF * abs(tH) ) < rndmPtr->flat() )
        continue;
      if ( isDiffB && pow( s4 / s, xtCorPF * abs(tH) ) < rndmPtr->flat() )
        continue;

    // Donnachie and Landshoff (RapGap):
    } else if (PomFlux == 4) {

      // Select diffractive mass(es) according to dm^2/(m^2)^(1 + 2 epsilon).
      m3 = m3ElDiff;
      m4 = m4ElDiff;
      if (isDiffA) {
        double s3MinPow = pow( m3ElDiff, xIntPF );
        double s3MaxPow = pow( max(mA, eCM - m4ElDiff), xIntPF );
        m3 = pow( s3MinPow + rndmPtr->flat() * (s3MaxPow - s3MinPow),
                  1. / xIntPF );
      }
      if (isDiffB) {
        double s4MinPow = pow( m4ElDiff, xIntPF );
        double s4MaxPow = pow( max(mB, eCM - m3ElDiff), xIntPF );
        m4 = pow( s4MinPow + rndmPtr->flat() * (s4MaxPow - s4MinPow),
                  1. / xIntPF );
      }
      if (m3 + m4 + DIFFMASSMARGIN >= eCM) continue;
      s3 = m3 * m3;
      s4 = m4 * m4;

      // Select t according to power and weigh by x_P^(2 alpha' |t|).
      tH = - (1. / pow( tAux1 + rndmPtr->flat() * (tAux2 - tAux1), 1./3.)
         - 1.) / coefDL;
      double wDL = pow2( (mp24DL - 2.8 * tH) / (mp24DL - tH) )
                 / pow4( 1. - tH / 0.7);
      double wMX = 1. / pow4( 1. - coefDL * tH);
      if (wDL < rndmPtr->flat() * wMX) continue;
      if ( isDiffA && pow( s3 / s, xtCorPF * abs(tH) ) < rndmPtr->flat() )
        continue;
      if ( isDiffB && pow( s4 / s, xtCorPF * abs(tH) ) < rndmPtr->flat() )
        continue;

    // MBR model:
    } else if (PomFlux == 5) {
      m3 = mA;
      m4 = mB;
      double xi, P, yRnd, dy;

      // MBR double diffractive.
      if (isDiffA && isDiffB) {
        dymin0 = 0.;
        dymax  = log(s/pow2(m2min));

        // Von Neumann method to generate dy, uses ddpmax from SigmaTot.
        do {
          dy = dymin0 + (dymax - dymin0) * rndmPtr->flat();
          P  = (dymax - dy) * exp(eps*dy) * ( exp(-2. * alph * dy * exp(-dy))
             - exp(-2. * alph * dy * exp(dy)) ) / dy;
          // Suppress smaller gap, smooth transition to non-diffractive.
          P *= 0.5 * (1 + erf( ( dy - dyminDD) / dyminSigDD ) );
          if (P > ddpmax) {
            ostringstream osWarn;
            osWarn << "ddpmax = " << scientific << setprecision(3)
                   << ddpmax << " " << P << " " << dy << endl;
            infoPtr->errorMsg("Warning in PhaseSpace2to2diffractive::"
              "trialKin for double diffraction:", osWarn.str());
          }
          yRnd = ddpmax * rndmPtr->flat();
        } while (yRnd > P);

        double y0max = (dymax - dy)/2.;
        double y0min = -y0max;
        double y0    = y0min + (y0max - y0min) * rndmPtr->flat();
        am1          = sqrt( eCM * exp( -y0 - dy/2. ) );
        am2          = sqrt( eCM * exp(  y0 - dy/2. ) );

        // Generate 4-momentum transfer, t from exp.
        double b = 2. * alph * dy;
        tUpp     = -exp( -dy );
        tLow     = -exp( dy );
        tAux     = exp( b * (tLow - tUpp) ) - 1.;
        t        = tUpp + log(1. + tAux * rndmPtr->flat()) / b;
        m3       = am1;
        m4       = am2;
        tH       = t;

      // MBR single diffractive.
      } else if (isDiffA || isDiffB) {
        dymin0 = 0.;
        dymax  = log(s/m2min);

        // Von Neumann method to generate dy, uses sdpmax from SigmaTot.
        do {
          dy = dymin0 + (dymax - dymin0) * rndmPtr->flat();
          P  = exp(eps * dy) * ( (FFA1 / (FFB1 + 2. * alph * dy) )
             + (FFA2 / (FFB2 + 2. * alph * dy) ) );
          // Suppress smaller gap.
          P *= 0.5 * (1. + erf( (dy - dyminSD) / dyminSigSD) );
          if (P > sdpmax) {
            ostringstream osWarn;
            osWarn << "sdpmax = " << scientific << setprecision(3)
                   << sdpmax << " " << P << " " << dy << endl;
            infoPtr->errorMsg("Warning in PhaseSpace2to2diffractive::"
              "trialKin for single diffraction:", osWarn.str());
          }
          yRnd = sdpmax * rndmPtr->flat();
        } while (yRnd > P);
        xi  = exp( -dy );
        amx = sqrt( xi * s );

        // Generate 4-momentum transfer, t. First exponent, then FF*exp.
        double tmin = -s1 * xi * xi / (1 - xi);
        do {
          t          = tmin + log(1. - rndmPtr->flat());
          double pFF = (4. * s1 - 2.8 * t) / ( (4. * s1 - t)
                     * pow2(1. - t / 0.71) );
          P          = pow2(pFF) * exp(2. * alph * dy * t);
          yRnd       = exp(t) * rndmPtr->flat();
        } while (yRnd > P);
        if(isDiffA) m3 = amx;
        if(isDiffB) m4 = amx;
        tH = t;
      }

      // End of MBR model code.
      s3 = m3 * m3;
      s4 = m4 * m4;
    }

    // Check whether m^2 and t choices are consistent.
    lambda34 = sqrtpos( pow2( s - s3 - s4) - 4. * s3 * s4 );
    double tempA = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4) / s;
    double tempB = lambda12 *  lambda34 / s;
    double tempC = (s3 - s1) * (s4 - s2) + (s1 + s4 - s2 - s3)
      * (s1 * s4 - s2 * s3) / s;
    double tLowNow = -0.5 * (tempA + tempB);
    double tUppNow = tempC / tLowNow;
    if (tH < tLowNow || tH > tUppNow) continue;

    // Careful reconstruction of scattering angle.
    double cosTheta = min(1., max(-1., (tempA + 2. * tH) / tempB));
    double sinTheta = 2. * sqrtpos( -(tempC + tempA * tH + tH * tH) )
      / tempB;
    theta = asin( min(1., sinTheta));
    if (cosTheta < 0.) theta = M_PI - theta;

    // Found acceptable kinematics, so no more looping. Done
    break;
  }
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

  // Done.
  return true;

}

//==========================================================================

// PhaseSpace2to3diffractive class.
// 2 -> 3 kinematics set up for central diffractive scattering.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of tries to find acceptable (m^2, t1, t2) set.
const int PhaseSpace2to3diffractive::NTRY = 500;
const int PhaseSpace2to3diffractive::NINTEG2 = 40;

// Maximum positive/negative argument for exponentiation.
const double PhaseSpace2to3diffractive::EXPMAX = 50.;

// Minimal mass of central diffractive system.
const double PhaseSpace2to3diffractive::DIFFMASSMIN = 0.8;

// Safety margin so sum of diffractive masses not too close to eCM.
const double PhaseSpace2to3diffractive::DIFFMASSMARGIN = 0.2;

//--------------------------------------------------------------------------

// Set up for phase space sampling.

bool PhaseSpace2to3diffractive::setupSampling() {

  // Pomeron flux parametrization, and parameters of some options.
  PomFlux      = settingsPtr->mode("Diffraction:PomFlux");
  epsilonPF    = settingsPtr->parm("Diffraction:PomFluxEpsilon");
  alphaPrimePF = settingsPtr->parm("Diffraction:PomFluxAlphaPrime");

  // Find maximum = value of cross section.
  sigmaNw      = sigmaProcessPtr->sigmaHatWrap();
  sigmaMx      = sigmaNw;

  // Squared masses of particles and minimal mass of diffractive states.
  s1           = mA * mA;
  s2           = mB * mB;
  m5min        = sigmaTotPtr->mMinAXB();
  s5min        = m5min * m5min;

  // Loop over two cases: s4 = (X + B)^2 and s3 = (A + X)^2.
  for (int i = 0; i < 2; ++i) {
    s3 = (i == 0) ? s1 : pow2(mA + m5min);
    s4 = (i == 0) ? pow2(mB + m5min) : s2;

    // Determine maximum possible t range and coefficient of generation.
    double lambda12 = sqrtpos( pow2( s - s1 - s2) - 4. * s1 * s2 );
    double lambda34 = sqrtpos( pow2( s - s3 - s4) - 4. * s3 * s4 );
    double tempA    = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4) / s;
    double tempB    = lambda12 *  lambda34 / s;
    double tempC    = (s3 - s1) * (s4 - s2) + (s1 + s4 - s2 - s3)
                    * (s1 * s4 - s2 * s3) / s;
    tLow[i]         = -0.5 * (tempA + tempB);
    tUpp[i]         = tempC / tLow[i];
  }
  s3 = s1;
  s4 = s2;

  // Default for all parametrization-specific parameters.
  bSlope1 = bSlope2 = bSlope = xIntPF = xIntInvPF = xtCorPF = mp24DL
    = coefDL = 0.;
  for (int i = 0; i < 2; ++i)
    bMin[i] = tAux[i] = probSlope1[i] = tAux1[i] = tAux2[i] = 0.;

  // Schuler&Sjostrand: lower limit diffractive slope.
  if (PomFlux == 1) {
    bMin[0] = sigmaTotPtr->bMinSlopeAX();
    tAux[0] = exp( max(-EXPMAX, bMin[0] * (tLow[0] - tUpp[0])) ) - 1.;
    bMin[1] = sigmaTotPtr->bMinSlopeXB();
    tAux[1] = exp( max(-EXPMAX, bMin[1] * (tLow[1] - tUpp[1])) ) - 1.;

  // Bruni&Ingelman: relative weight of two diffractive slopes.
  } else if (PomFlux == 2) {
    bSlope1     = 8.0;
    bSlope2     = 3.0;
    for (int i = 0; i < 2; ++i) {
      probSlope1[i]  = 6.38 * ( exp(max(-EXPMAX, bSlope1 * tUpp[i]))
                     -  exp(max(-EXPMAX, bSlope1 * tLow[i])) ) / bSlope1;
      double pS2     = 0.424 * ( exp(max(-EXPMAX, bSlope2 * tUpp[i]))
                     -  exp(max(-EXPMAX, bSlope2 * tLow[i])) ) / bSlope2;
      probSlope1[i] /= probSlope1[i] + pS2;
      tAux1[i] = exp( max(-EXPMAX, bSlope1 * (tLow[i] - tUpp[i])) ) - 1.;
      tAux2[i] = exp( max(-EXPMAX, bSlope2 * (tLow[i] - tUpp[i])) ) - 1.;
    }

  // Streng&Berger (RapGap): diffractive slope, power of mass spectrum.
  } else if (PomFlux == 3) {
    bSlope        = 4.7;
    double xPowPF = 1. - 2. * (1. + epsilonPF);
    xIntPF        = 1. + xPowPF;
    xIntInvPF     = 1. / xIntPF;
    xtCorPF       = 2. * alphaPrimePF;
    tAux[0]       = exp( max(-EXPMAX, bSlope  * (tLow[0] - tUpp[0])) ) - 1.;
    tAux[1]       = exp( max(-EXPMAX, bSlope  * (tLow[1] - tUpp[1])) ) - 1.;

  // Donnachie&Landshoff (RapGap):  power of mass spectrum.
  } else if (PomFlux == 4) {
    mp24DL        = 4. * pow2(particleDataPtr->m0(2212));
    double xPowPF = 1. - 2. * (1. + epsilonPF);
    xIntPF        = 1. + xPowPF;
    xIntInvPF     = 1. / xIntPF;
    xtCorPF       = 2. * alphaPrimePF;
    // Upper estimate of t dependence, for preliminary choice.
    coefDL        = 0.85;
    tAux1[0]      = 1. / pow3(1. - coefDL * tLow[0]);
    tAux2[0]      = 1. / pow3(1. - coefDL * tUpp[0]);
    tAux1[1]      = 1. / pow3(1. - coefDL * tLow[1]);
    tAux2[1]      = 1. / pow3(1. - coefDL * tUpp[1]);

  // Setup for the MBR model.
  } else if (PomFlux == 5) {
    epsMBR        = settingsPtr->parm("Diffraction:MBRepsilon");
    alphMBR       = settingsPtr->parm("Diffraction:MBRalpha");
    m2minMBR      = settingsPtr->parm("Diffraction:MBRm2Min");
    dyminMBR      = settingsPtr->parm("Diffraction:MBRdyminCD");
    dyminSigMBR   = settingsPtr->parm("Diffraction:MBRdyminSigCD");
    dyminInvMBR   = sqrt(2.) / dyminSigMBR;
    // Max f(dy) for Von Neumann method, dpepmax from SigmaTot.
    dpepmax       = sigmaTotPtr->dpepMax();

  // H1 Fit A/B.
  } else if (PomFlux == 6 || PomFlux == 7) {
    bSlope        = 5.5;
    epsilonPF     = (PomFlux == 6) ? 0.1182 : 0.1110;
    alphaPrimePF  = 0.06;
    double xPowPF = 1. - 2. * (1. + epsilonPF);
    xIntPF        = 1. + xPowPF;
    xIntInvPF     = 1. / xIntPF;
    xtCorPF       = 2. * alphaPrimePF;
    tAux[0]       = exp( max(-EXPMAX, bSlope  * (tLow[0] - tUpp[0])) ) - 1.;
    tAux[1]       = exp( max(-EXPMAX, bSlope  * (tLow[1] - tUpp[1])) ) - 1.;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select a trial kinematics phase space point. Perform full
// Monte Carlo acceptance/rejection at this stage.

bool PhaseSpace2to3diffractive::trialKin( bool, bool ) {

  // Allow for possibility that energy varies from event to event.
  if (doEnergySpread) {
    eCM       = infoPtr->eCM();
    s         = eCM * eCM;
    for (int i = 0; i < 2; ++i) {
      s3 = (i == 0) ? s1 : pow2(mA + m5min);
      s4 = (i == 0) ? pow2(mB + m5min) : s2;
      double lambda12 = sqrtpos( pow2( s - s1 - s2) - 4. * s1 * s2 );
      double lambda34 = sqrtpos( pow2( s - s3 - s4) - 4. * s3 * s4 );
      double tempA    = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4) / s;
      double tempB    = lambda12 *  lambda34 / s;
      double tempC    = (s3 - s1) * (s4 - s2) + (s1 + s4 - s2 - s3)
                      * (s1 * s4 - s2 * s3) / s;
      tLow[i]         = -0.5 * (tempA + tempB);
      tUpp[i]         = tempC / tLow[i];
    }
    s3 = s1;
    s4 = s2;
    if (PomFlux == 1) {
      tAux[0] = exp( max(-EXPMAX, bMin[0] * (tLow[0] - tUpp[0])) ) - 1.;
      tAux[1] = exp( max(-EXPMAX, bMin[1] * (tLow[1] - tUpp[1])) ) - 1.;
    } else if (PomFlux == 2) {
      for (int i = 0; i < 2; ++i) {
        tAux1[i] = exp( max(-EXPMAX, bSlope1 * (tLow[i] - tUpp[i])) ) - 1.;
        tAux2[i] = exp( max(-EXPMAX, bSlope2 * (tLow[i] - tUpp[i])) ) - 1.;
      }
    } else if (PomFlux == 3 || PomFlux == 6 || PomFlux == 7) {
      tAux[0]       = exp( max(-EXPMAX, bSlope  * (tLow[0] - tUpp[0])) ) - 1.;
      tAux[1]       = exp( max(-EXPMAX, bSlope  * (tLow[1] - tUpp[1])) ) - 1.;
    } else if (PomFlux == 4) {
      tAux1[0]      = 1. / pow3(1. - coefDL * tLow[0]);
      tAux2[0]      = 1. / pow3(1. - coefDL * tUpp[0]);
      tAux1[1]      = 1. / pow3(1. - coefDL * tLow[1]);
      tAux2[1]      = 1. / pow3(1. - coefDL * tUpp[1]);
    }
  }

  // Trivial kinematics of incoming hadrons.
  double lambda12 = sqrtpos( pow2( s - s1 - s2) - 4. * s1 * s2 );
  pAbs            = 0.5 * lambda12 / eCM;
  p1.p( 0., 0.,  pAbs, 0.5 * (s + s1 - s2) / eCM);
  p2.p( 0., 0., -pAbs, 0.5 * (s + s2 - s1) / eCM);

  // Loop over attempts to set up mass, t1, t2 consistently.
  for (int loop = 0; ; ++loop) {
    if (loop == NTRY) {
      infoPtr->errorMsg("Error in PhaseSpace2to3diffractive::trialKin: "
      " quit after repeated tries");
      return false;
    }
    double xi1 = 0.;
    double xi2 = 0.;
    double tVal[2] = { 0., 0.};

    // Schuler and Sjostrand:
    if (PomFlux == 1) {

      // Select mass according to dxi_1/xi_1 * dxi_2/xi_2 * (1 - m^2/s).
      do {
        xi1 = pow( s5min / s, rndmPtr->flat());
        xi2 = pow( s5min / s, rndmPtr->flat());
        s5 = xi1 * xi2 * s;
      } while (s5 < s5min || xi1 * xi2 > rndmPtr->flat());
      if (mA + mB + sqrt(s5) + DIFFMASSMARGIN >= eCM) continue;

      // Select t according to exp(bMin*t) and correct to right slope.
      bool tryAgain = false;
      for (int i = 0; i < 2; ++i) {
        tVal[i] = tUpp[i] + log(1. + tAux[i] * rndmPtr->flat()) / bMin[i];
        double bDiff = (i == 0) ? sigmaTotPtr->bSlopeAX(s2 + xi1 * s)
                                : sigmaTotPtr->bSlopeXB(s1 + xi2 * s);
        bDiff = max(0., bDiff - bMin[i]);
        if (exp( max(-EXPMAX, bDiff * (tVal[i] - tUpp[i])) )
          < rndmPtr->flat()) tryAgain = true;
      }
      if (tryAgain) continue;

    // Bruni and Ingelman:
    } else if (PomFlux == 2) {

      // Select mass according to dxi_1/xi_1 * dxi_2/xi_2.
      do {
        xi1 = pow( s5min / s, rndmPtr->flat());
        xi2 = pow( s5min / s, rndmPtr->flat());
        s5 = xi1 * xi2 * s;
      } while (s5 < s5min);
      if (mA + mB + sqrt(s5) + DIFFMASSMARGIN >= eCM) continue;

      // Select t according to exp(bSlope*t) with two possible slopes.
      for (int i = 0; i < 2; ++i)
        tVal[i] = (rndmPtr->flat() < probSlope1[i])
                ? tUpp[i] + log(1. + tAux1[i] * rndmPtr->flat()) / bSlope1
                : tUpp[i] + log(1. + tAux2[i] * rndmPtr->flat()) / bSlope2;

    // Streng and Berger et al. (RapGap) and H1 Fit A/B:
    } else if (PomFlux == 3 || PomFlux == 6 || PomFlux == 7) {

      // Select mass by dxi_1 * dxi_2 / (xi_1 * xi_2)^(1 + 2 epsilon).
      double sMinPow = pow( s5min / s, xIntPF);
      do {
        xi1 = pow( sMinPow + rndmPtr->flat() * (1. - sMinPow), xIntInvPF );
        xi2 = pow( sMinPow + rndmPtr->flat() * (1. - sMinPow), xIntInvPF );
        s5 = xi1 * xi2 * s;
      } while (s5 < s5min);
      if (mA + mB + sqrt(s5) + DIFFMASSMARGIN >= eCM) continue;

      // Select t according to exponential and weigh by x_P^(2 alpha' |t|).
      bool tryAgain = false;
      for (int i = 0; i < 2; ++i) {
        tVal[i] = tUpp[i] + log(1. + tAux[i] * rndmPtr->flat()) / bSlope;
        double xi = (i == 0) ? xi1 : xi2;
        if ( pow( xi, xtCorPF * abs(tVal[i]) ) < rndmPtr->flat() )
          tryAgain = true;
      }
      if (tryAgain) continue;

    // Donnachie and Landshoff (RapGap):
    } else if (PomFlux == 4) {

      // Select mass by dxi_1 * dxi_2 / (xi_1 * xi_2)^(1 + 2 epsilon).
      double sMinPow = pow( s5min / s, xIntPF);
      do {
        xi1 = pow( sMinPow + rndmPtr->flat() * (1. - sMinPow), xIntInvPF );
        xi2 = pow( sMinPow + rndmPtr->flat() * (1. - sMinPow), xIntInvPF );
        s5 = xi1 * xi2 * s;
      } while (s5 < s5min);
      if (mA + mB + sqrt(s5) + DIFFMASSMARGIN >= eCM) continue;

      // Select t according to power and weigh by x_P^(2 alpha' |t|).
      bool tryAgain = false;
      for (int i = 0; i < 2; ++i) {
        tVal[i] = - (1. / pow( tAux1[i] + rndmPtr->flat()
                * (tAux2[i] - tAux1[i]), 1./3.) - 1.) / coefDL;
        double wDL = pow2( (mp24DL - 2.8 * tVal[i]) / (mp24DL - tVal[i]) )
                   / pow4( 1. - tVal[i] / 0.7);
        double wMX = 1. / pow4( 1. - coefDL * tVal[i]);
        if (wDL < rndmPtr->flat() * wMX) tryAgain = true;
        double xi = (i == 0) ? xi1 : xi2;
        if ( pow( xi, xtCorPF * abs(tVal[i]) ) < rndmPtr->flat() )
          tryAgain = true;
      }
      if (tryAgain) continue;

    // The MBR model (PomFlux == 5).
    } else if (PomFlux == 5) {
      double dymin0 = 0.;
      double dymax  = log(s/m2minMBR);
      double f1, f2, step2, dy, yc, ycmin, ycmax, dy1, dy2,
             P, P1, P2, yRnd, yRnd1, yRnd2;

      // Von Neumann method to generate dy, uses dpepmax from SigmaTot.
      do {
        dy    = dymin0 + (dymax - dymin0) * rndmPtr->flat();
        P     = 0.;
        step2 = (dy - dymin0) / NINTEG2;
        for (int j = 0; j < NINTEG2 ; ++j) {
          yc  = -(dy - dymin0) / 2. + (j + 0.5) * step2;
          dy1 = 0.5 * dy - yc;
          dy2 = 0.5 * dy + yc;
          f1  = exp(epsMBR * dy1) * ( (FFA1 / (FFB1 + 2. * alphMBR * dy1) )
              + (FFA2 / (FFB2 + 2. * alphMBR * dy1) ) );
          f2  = exp(epsMBR * dy2) * ( (FFA1 / (FFB1 + 2. * alphMBR * dy2) )
              + (FFA2 / (FFB2 + 2. * alphMBR * dy2) ) );
          f1 *= 0.5 * (1. + erf( (dy1 - 0.5 * dyminMBR) * dyminInvMBR ));
          f2 *= 0.5 * (1. + erf( (dy2 - 0.5 * dyminMBR) * dyminInvMBR ));
          P  += f1 * f2 * step2;
        }
        if (P > dpepmax) {
          ostringstream osWarn;
          osWarn << "dpepmax = " << scientific << setprecision(3)
                 << dpepmax << " " << P << " " << dy << endl;
          infoPtr->errorMsg("Warning in PhaseSpace2to2diffractive::"
            "trialKin for central diffraction:", osWarn.str());
        }
        yRnd = dpepmax * rndmPtr->flat();

        // Generate dyc.
        ycmax = (dy - dymin0) / 2.;
        ycmin = -ycmax;
        yc    = ycmin + (ycmax - ycmin) * rndmPtr->flat();

        // xi1, xi2 from dy and dy0.
        dy1 = 0.5 * dy + yc;
        dy2 = 0.5 * dy - yc;
        P1  = 0.5 * (1. + erf( (dy1 - 0.5 * dyminMBR) * dyminInvMBR ));
        P2  = 0.5 * (1. + erf( (dy2 - 0.5 * dyminMBR) * dyminInvMBR ));
        yRnd1 = rndmPtr->flat();
        yRnd2 = rndmPtr->flat();
      } while( !(yRnd < P && yRnd1 < P1 && yRnd2 < P2) );
      xi1 = exp( -dy1 );
      xi2 = exp( -dy2 );

      // Generate t1 at vertex1. First exponent, then FF*exp.
      double tmin  = -s1 * xi1 * xi1 / (1. - xi1);
      do {
        t1         = tmin + log(1. - rndmPtr->flat());
        double pFF = (4. * s1 - 2.8 * t1) / ( (4. * s1 - t1)
                   * pow2(1. - t1 / 0.71));
        P          = pow2(pFF) * exp(2. * alphMBR * dy1 * t1);
        yRnd       = exp(t1) * rndmPtr->flat();
      } while (yRnd > P);

      // Generate t2 at vertex2. First exponent, then FF*exp.
      tmin         = -s2 * xi2 * xi2 / (1. - xi2);
      do {
        t2         = tmin + log(1. - rndmPtr->flat());
        double pFF = (4. * s2 - 2.8 * t2) / ((4. * s2 - t2)
                   * pow2(1. - t2 / 0.71));
        P          = pow2(pFF) * exp(2. * alphMBR * dy2 * t2);
        yRnd       = exp(t2) * rndmPtr->flat();
      } while (yRnd > P);

    }

    // Checks and kinematics construction four first options.
    double pz3 = 0.;
    double pz4 = 0.;
    double pT3 = 0.;
    double pT4 = 0.;
    if (PomFlux != 5) {

      // Check whether m^2 (i.e. xi) and t choices are consistent.
      bool tryAgain   = false;
      for (int i = 0; i < 2; ++i) {
        double sx1 = (i == 0) ? s1 : s2;
        double sx2 = (i == 0) ? s2 : s1;
        double sx3 = sx1;
        double sx4 = (i == 0) ? s2 + xi1 * s : s1 + xi2 * s;
        if (sqrt(sx3) + sqrt(sx4) + DIFFMASSMARGIN > eCM) tryAgain = true;
        double lambda34 = sqrtpos( pow2( s - sx3 - sx4) - 4. * sx3 * sx4 );
        double tempA    = s - (sx1 + sx2 + sx3 + sx4)
                        + (sx1 - sx2) * (sx3 - sx4) / s;
        double tempB    = lambda12 * lambda34 / s;
        double tempC    = (sx3 - sx1) * (sx4 - sx2) + (sx1 + sx4 - sx2 - sx3)
                        * (sx1 * sx4 - sx2 * sx3) / s;
        double tLowNow  = -0.5 * (tempA + tempB);
        double tUppNow  = tempC / tLowNow;
        if (tVal[i] < tLowNow || tVal[i] > tUppNow) tryAgain = true;
        if (tryAgain) break;

        // Careful reconstruction of scattering angle.
        double cosTheta = min(1., max(-1., (tempA + 2. * tVal[i]) / tempB));
        double sinTheta = 2. * sqrtpos( -(tempC + tempA * tVal[i]
                        + tVal[i] * tVal[i]) ) / tempB;
        theta           = asin( min(1., sinTheta));
        if (cosTheta < 0.) theta = M_PI - theta;
        double pAbs34   = 0.5 * lambda34 / eCM;
        if (i == 0) {
          pz3   =  pAbs34 * cos(theta);
          pT3   =  pAbs34 * sin(theta);
        } else {
          pz4   = -pAbs34 * cos(theta);
          pT4   =  pAbs34 * sin(theta);
        }
      }
      if (tryAgain) continue;
      t1        = tVal[0];
      t2        = tVal[1];

    // Kinematics construction in the MBR model.
    } else {
      pz3       =  pAbs * (1. - xi1);
      pz4       = -pAbs * (1. - xi2);
      pT3       =  sqrt( (1. - xi1) * abs(t1) - s1 * pow2(xi1) );
      pT4       =  sqrt( (1. - xi2) * abs(t2) - s2 * pow2(xi2) );
    }

    // Common final steps of kinematics.
    double phi3 = 2. * M_PI * rndmPtr->flat();
    double phi4 = 2. * M_PI * rndmPtr->flat();
    p3.p( pT3 * cos(phi3), pT3 * sin(phi3), pz3,
          sqrt(pz3 * pz3 + pT3 * pT3 + s1) );
    p4.p( pT4 * cos(phi4), pT4 * sin(phi4), pz4,
          sqrt(pz4 * pz4 + pT4 * pT4 + s2) );

    // Central dissociated system, from Pomeron-Pomeron 4 vectors.
    p5   = (p1 - p3) + (p2 - p4);
    mHat = p5.mCalc();

    // If acceptable diffractive mass then no more looping.
    if (mHat > DIFFMASSMIN) break;
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
  mH[5] = mHat;

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
  pH[1] = Vec4( 0., 0., 0.5 * eCM * x1H, 0.5 * eCM * x1H);
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
// Note: here cout is used for output, not os. Change??

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
