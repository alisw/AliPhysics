// MultipartonInteractions.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// SigmaMultiparton and MultipartonInteractions classes.

#include "Pythia8/MultipartonInteractions.h"

// Internal headers for special processes.
#include "Pythia8/SigmaQCD.h"
#include "Pythia8/SigmaEW.h"
#include "Pythia8/SigmaOnia.h"

namespace Pythia8 {

//==========================================================================

// The SigmaMultiparton class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// The sum of outgoing masses must not be too close to the cm energy.
const double SigmaMultiparton::MASSMARGIN = 0.1;

// Fraction of time not the dominant "gluon t-channel" process is picked.
const double SigmaMultiparton::OTHERFRAC  = 0.2;

//--------------------------------------------------------------------------

// Initialize the generation process for given beams.

bool SigmaMultiparton::init(int inState, int processLevel, Info* infoPtr,
    Settings* settingsPtr, ParticleData* particleDataPtr, Rndm* rndmPtrIn,
    BeamParticle* beamAPtr, BeamParticle* beamBPtr, Couplings* couplingsPtr) {

  // Store input pointer for future use.
  rndmPtr          = rndmPtrIn;

  // Reset vector sizes (necessary in case of re-initialization).
  if (sigmaT.size() > 0) {
    for (int i = 0; i < int(sigmaT.size()); ++i) delete sigmaT[i];
    sigmaT.resize(0);
  }
  if (sigmaU.size() > 0) {
    for (int i = 0; i < int(sigmaU.size()); ++i) delete sigmaU[i];
    sigmaU.resize(0);
  }

  // Always store mimimal set of processes:QCD 2 -> 2 t-channel.

  // Gluon-gluon instate.
  if (inState == 0) {
    sigmaT.push_back( new Sigma2gg2gg() );
    sigmaU.push_back( new Sigma2gg2gg() );

  // Quark-gluon instate.
  } else if (inState == 1) {
    sigmaT.push_back( new Sigma2qg2qg() );
    sigmaU.push_back( new Sigma2qg2qg() );

  // Quark-(anti)quark instate.
  } else {
    sigmaT.push_back( new Sigma2qq2qq() );
    sigmaU.push_back( new Sigma2qq2qq() );
  }

  // Normally store QCD processes to new flavour.
  if (processLevel > 0) {
    if (inState == 0) {
      sigmaT.push_back( new Sigma2gg2qqbar() );
      sigmaU.push_back( new Sigma2gg2qqbar() );
      sigmaT.push_back( new Sigma2gg2QQbar(4, 121) );
      sigmaU.push_back( new Sigma2gg2QQbar(4, 121) );
      sigmaT.push_back( new Sigma2gg2QQbar(5, 123) );
      sigmaU.push_back( new Sigma2gg2QQbar(5, 123) );
    } else if (inState == 2) {
      sigmaT.push_back( new Sigma2qqbar2gg() );
      sigmaU.push_back( new Sigma2qqbar2gg() );
      sigmaT.push_back( new Sigma2qqbar2qqbarNew() );
      sigmaU.push_back( new Sigma2qqbar2qqbarNew() );
      sigmaT.push_back( new Sigma2qqbar2QQbar(4, 122) );
      sigmaU.push_back( new Sigma2qqbar2QQbar(4, 122) );
      sigmaT.push_back( new Sigma2qqbar2QQbar(5, 124) );
      sigmaU.push_back( new Sigma2qqbar2QQbar(5, 124) );
    }
  }

  // Optionally store electroweak processes, mainly photon production.
  if (processLevel > 1) {
    if (inState == 0) {
      sigmaT.push_back( new Sigma2gg2ggamma() );
      sigmaU.push_back( new Sigma2gg2ggamma() );
      sigmaT.push_back( new Sigma2gg2gammagamma() );
      sigmaU.push_back( new Sigma2gg2gammagamma() );
    } else if (inState == 1) {
      sigmaT.push_back( new Sigma2qg2qgamma() );
      sigmaU.push_back( new Sigma2qg2qgamma() );
    } else if (inState == 2) {
      sigmaT.push_back( new Sigma2qqbar2ggamma() );
      sigmaU.push_back( new Sigma2qqbar2ggamma() );
      sigmaT.push_back( new Sigma2ffbar2gammagamma() );
      sigmaU.push_back( new Sigma2ffbar2gammagamma() );
      sigmaT.push_back( new Sigma2ffbar2ffbarsgm() );
      sigmaU.push_back( new Sigma2ffbar2ffbarsgm() );
    }
    if (inState >= 2) {
      sigmaT.push_back( new Sigma2ff2fftgmZ() );
      sigmaU.push_back( new Sigma2ff2fftgmZ() );
      sigmaT.push_back( new Sigma2ff2fftW() );
      sigmaU.push_back( new Sigma2ff2fftW() );
    }
  }

  // Optionally store charmonium and bottomonium production.
  if (processLevel > 2) {
    SigmaOniaSetup charmonium(infoPtr, settingsPtr, particleDataPtr, 4);
    SigmaOniaSetup bottomonium(infoPtr, settingsPtr, particleDataPtr, 5);
    if (inState == 0) {
      charmonium.setupSigma2gg(sigmaT, true);
      charmonium.setupSigma2gg(sigmaU, true);
      bottomonium.setupSigma2gg(sigmaT, true);
      bottomonium.setupSigma2gg(sigmaU, true);
    } else if (inState == 1) {
      charmonium.setupSigma2qg(sigmaT, true);
      charmonium.setupSigma2qg(sigmaU, true);
      bottomonium.setupSigma2qg(sigmaT, true);
      bottomonium.setupSigma2qg(sigmaU, true);
    } else if (inState == 2) {
      charmonium.setupSigma2qq(sigmaT, true);
      charmonium.setupSigma2qq(sigmaU, true);
      bottomonium.setupSigma2qq(sigmaT, true);
      bottomonium.setupSigma2qq(sigmaU, true);
    }
  }

  // Resize arrays to match sizes above.
  nChan = sigmaT.size();
  needMasses.resize(nChan);
  m3Fix.resize(nChan);
  m4Fix.resize(nChan);
  sHatMin.resize(nChan);
  sigmaTval.resize(nChan);
  sigmaUval.resize(nChan);

  // Initialize the processes.
  for (int i = 0; i < nChan; ++i) {
    sigmaT[i]->init( infoPtr, settingsPtr, particleDataPtr, rndmPtr,
      beamAPtr, beamBPtr, couplingsPtr);
    sigmaT[i]->initProc();
    sigmaU[i]->init( infoPtr, settingsPtr, particleDataPtr, rndmPtr,
      beamAPtr, beamBPtr, couplingsPtr);
    sigmaU[i]->initProc();

    // Prepare for massive kinematics (but fixed masses!) where required.
    needMasses[i] = false;
    int id3Mass =  sigmaT[i]->id3Mass();
    int id4Mass =  sigmaT[i]->id4Mass();
    m3Fix[i] = 0.;
    m4Fix[i] = 0.;
    if (id3Mass > 0 || id4Mass > 0) {
      needMasses[i] = true;
      m3Fix[i] =  particleDataPtr->m0(id3Mass);
      m4Fix[i] =  particleDataPtr->m0(id4Mass);
    }
    sHatMin[i] = pow2( m3Fix[i] + m4Fix[i] + MASSMARGIN);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Calculate cross section summed over possibilities.

double SigmaMultiparton::sigma( int id1, int id2, double x1, double x2,
  double sHat, double tHat, double uHat, double alpS, double alpEM,
  bool restore, bool pickOtherIn) {

  // Choose either the dominant process (in slot 0) or the rest of them.
  if (restore) pickOther = pickOtherIn;
  else         pickOther = (rndmPtr->flat() < OTHERFRAC);

  // Iterate over all subprocesses.
  sigmaTsum = 0.;
  sigmaUsum = 0.;
  for (int i = 0; i < nChan; ++i) {
    sigmaTval[i] = 0.;
    sigmaUval[i] = 0.;

    // Skip the not chosen processes.
    if (i == 0 && pickOther) continue;
    if (i > 0 && !pickOther) continue;

    // t-channel-sampling contribution.
    if (sHat > sHatMin[i]) {
      sigmaT[i]->set2KinMPI( x1, x2, sHat, tHat, uHat,
        alpS, alpEM, needMasses[i], m3Fix[i], m4Fix[i]);
      sigmaTval[i] = sigmaT[i]->sigmaHatWrap(id1, id2);
      sigmaT[i]->pickInState(id1, id2);
      // Correction factor for tHat rescaling in massive kinematics.
      if (needMasses[i]) sigmaTval[i] *= sigmaT[i]->sHBetaMPI() / sHat;
      sigmaTsum += sigmaTval[i];
    }

    // u-channel-sampling contribution.
    if (sHat > sHatMin[i]) {
      sigmaU[i]->set2KinMPI( x1, x2, sHat, uHat, tHat,
        alpS, alpEM, needMasses[i], m3Fix[i], m4Fix[i]);
      sigmaUval[i] = sigmaU[i]->sigmaHatWrap( id1, id2);
      sigmaU[i]->pickInState(id1, id2);
      // Correction factor for tHat rescaling in massive kinematics.
      if (needMasses[i]) sigmaUval[i] *= sigmaU[i]->sHBetaMPI() / sHat;
      sigmaUsum += sigmaUval[i];
    }

  // Average of t- and u-channel sampling; corrected for not selected channels.
  }
  double sigmaAvg = 0.5 * (sigmaTsum + sigmaUsum);
  if (pickOther) sigmaAvg /= OTHERFRAC;
  if (!pickOther) sigmaAvg /= (1. - OTHERFRAC);
  return sigmaAvg;

}

//--------------------------------------------------------------------------

// Return one subprocess, picked according to relative cross sections.

SigmaProcess* SigmaMultiparton::sigmaSel() {

  // Decide between t- and u-channel-sampled kinematics.
  pickedU = (rndmPtr->flat() * (sigmaTsum + sigmaUsum) < sigmaUsum);

  // Pick one of t-channel-sampled processes.
  if (!pickedU) {
    double sigmaRndm = sigmaTsum * rndmPtr->flat();
    int    iPick = -1;
    do     sigmaRndm -= sigmaTval[++iPick];
    while  (sigmaRndm > 0.);
    return sigmaT[iPick];

  // Pick one of u-channel-sampled processes.
  } else {
    double sigmaRndm = sigmaUsum * rndmPtr->flat();
    int    iPick = -1;
    do     sigmaRndm -= sigmaUval[++iPick];
    while  (sigmaRndm > 0.);
    return sigmaU[iPick];
  }

}

//==========================================================================

// The MultipartonInteractions class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Factorization scale pT2 by default, but could be shifted to pT2 + pT02.
// (A priori possible, but problems for flavour threshold interpretation.)
const bool   MultipartonInteractions::SHIFTFACSCALE = false;

// Pick one parton to represent rescattering cross section to speed up.
const bool   MultipartonInteractions::PREPICKRESCATTER = true;

// Naive upper estimate of cross section too pessimistic, so reduce by this.
const double MultipartonInteractions::SIGMAFUDGE    = 0.8;

// The r value above, picked to allow a flatter correct/trial cross section.
const double MultipartonInteractions::RPT20         = 0.25;

// Reduce pT0 by factor pT0STEP if sigmaInt < SIGMASTEP * sigmaND.
const double MultipartonInteractions::PT0STEP       = 0.9;
const double MultipartonInteractions::SIGMASTEP     = 1.1;

// Stop if pT0 or pTmin fall below this, or alpha_s blows up.
const double MultipartonInteractions::PT0MIN        = 0.2;

// Refuse too low expPow in impact parameter profile.
const double MultipartonInteractions::EXPPOWMIN     = 0.4;

// Define low-b region by interaction rate above given number.
const double MultipartonInteractions::PROBATLOWB    = 0.6;

// Basic step size for b integration; sometimes modified.
const double MultipartonInteractions::BSTEP         = 0.01;

// Stop b integration when integrand dropped enough.
const double MultipartonInteractions::BMAX          = 1e-8;

// Do not allow too large argument to exp function.
const double MultipartonInteractions::EXPMAX        = 50.;

// Convergence criterion for k iteration.
const double MultipartonInteractions::KCONVERGE     = 1e-7;

// Conversion of GeV^{-2} to mb for cross section.
const double MultipartonInteractions::CONVERT2MB    = 0.389380;

// Stay away from division by zero in Jacobian for tHat -> pT2.
const double MultipartonInteractions::ROOTMIN       = 0.01;

// No need to reinitialize parameters if energy close to previous.
const double MultipartonInteractions::ECMDEV        = 0.01;

// Settings for x-dependent matter profile:
// Number of bins in b (with these settings, no bStep increase and
// reintegration needed with a1 ~ 0.20 up to ECM ~ 40TeV).
const int    MultipartonInteractions::XDEP_BBIN     = 500;
// Initial value of a0.
const double MultipartonInteractions::XDEP_A0       = 1.0;
// Width of form ( XDEP_A1 + a1 * log(1 / x) ).
const double MultipartonInteractions::XDEP_A1       = 1.0;
// Initial step size in b and increment.
const double MultipartonInteractions::XDEP_BSTEP    = 0.02;
const double MultipartonInteractions::XDEP_BSTEPINC = 0.01;
// Overlap-weighted cross section in last bin of b must be beneath
// XDEP_CUTOFF * sigmaInt.
const double MultipartonInteractions::XDEP_CUTOFF   = 1e-4;
// a0 is calculated in units of sqrt(mb), so convert to fermi.
const double MultipartonInteractions::XDEP_SMB2FM   = sqrt(0.1);

// Only write warning when weight clearly above unity.
const double MultipartonInteractions::WTACCWARN     = 1.1;

//--------------------------------------------------------------------------

// Initialize the generation process for given beams.

bool MultipartonInteractions::init( bool doMPIinit, int iDiffSysIn,
  Info* infoPtrIn, Settings& settings, ParticleData* particleDataPtr,
  Rndm* rndmPtrIn, BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  Couplings* couplingsPtrIn, PartonSystems* partonSystemsPtrIn,
  SigmaTotal* sigmaTotPtrIn, UserHooks* userHooksPtrIn, ostream& os) {

  // Store input pointers for future use. Done if no initialization.
  iDiffSys         = iDiffSysIn;
  infoPtr          = infoPtrIn;
  rndmPtr          = rndmPtrIn;
  beamAPtr         = beamAPtrIn;
  beamBPtr         = beamBPtrIn;
  couplingsPtr     = couplingsPtrIn;
  partonSystemsPtr = partonSystemsPtrIn;
  sigmaTotPtr      = sigmaTotPtrIn;
  userHooksPtr     = userHooksPtrIn;
  if (!doMPIinit) return false;

  // If both beams are baryons then softer PDF's than for mesons/Pomerons.
  hasBaryonBeams = ( beamAPtr->isBaryon() && beamBPtr->isBaryon() );

  // Matching in pT of hard interaction to further interactions.
  pTmaxMatch     = settings.mode("MultipartonInteractions:pTmaxMatch");

  //  Parameters of alphaStrong generation.
  alphaSvalue    = settings.parm("MultipartonInteractions:alphaSvalue");
  alphaSorder    = settings.mode("MultipartonInteractions:alphaSorder");
  alphaSnfmax    = settings.mode("StandardModel:alphaSnfmax");

  // Parameters of alphaEM generation.
  alphaEMorder   = settings.mode("MultipartonInteractions:alphaEMorder");

  //  Parameters of cross section generation.
  Kfactor        = settings.parm("MultipartonInteractions:Kfactor");

  // Regularization of QCD evolution for pT -> 0.
  pT0Ref         = settings.parm("MultipartonInteractions:pT0Ref");
  ecmRef         = settings.parm("MultipartonInteractions:ecmRef");
  ecmPow         = settings.parm("MultipartonInteractions:ecmPow");
  pTmin          = settings.parm("MultipartonInteractions:pTmin");

  // Impact parameter profile: nondiffractive topologies.
  if (iDiffSys == 0) {
    bProfile     = settings.mode("MultipartonInteractions:bProfile");
    coreRadius   = settings.parm("MultipartonInteractions:coreRadius");
    coreFraction = settings.parm("MultipartonInteractions:coreFraction");
    expPow       = settings.parm("MultipartonInteractions:expPow");
    expPow       = max(EXPPOWMIN, expPow);
    // x-dependent impact parameter profile.
    a1           = settings.parm("MultipartonInteractions:a1");

  // Impact parameter profile: diffractive topologies.
  } else {
    bProfile     = settings.mode("Diffraction:bProfile");
    coreRadius   = settings.parm("Diffraction:coreRadius");
    coreFraction = settings.parm("Diffraction:coreFraction");
    expPow       = settings.parm("Diffraction:expPow");
    expPow       = max(EXPPOWMIN, expPow);
  }

  // Common choice of "pT" scale for determining impact parameter.
  bSelScale      = settings.mode("MultipartonInteractions:bSelScale");

  // Process sets to include in machinery.
  processLevel   = settings.mode("MultipartonInteractions:processLevel");

  // Parameters of rescattering description.
  allowRescatter = settings.flag("MultipartonInteractions:allowRescatter");
  allowDoubleRes
    = settings.flag("MultipartonInteractions:allowDoubleRescatter");
  rescatterMode  = settings.mode("MultipartonInteractions:rescatterMode");
  ySepResc       = settings.parm("MultipartonInteractions:ySepRescatter");
  deltaYResc     = settings.parm("MultipartonInteractions:deltaYRescatter");

  // Rescattering not yet implemented for x-dependent impact profile.
  if (bProfile == 4) allowRescatter = false;

  // A global recoil FSR stategy restricts rescattering.
  globalRecoilFSR     = settings.flag("TimeShower:globalRecoil");
  nMaxGlobalRecoilFSR = settings.mode("TimeShower:nMaxGlobalRecoil");

  // Various other parameters.
  nQuarkIn       = settings.mode("MultipartonInteractions:nQuarkIn");
  nSample        = settings.mode("MultipartonInteractions:nSample");

  // Optional dampening at small pT's when large multiplicities.
  enhanceScreening = settings.mode("MultipartonInteractions:enhanceScreening");

  // Parameters for diffractive systems.
  sigmaPomP      = settings.parm("Diffraction:sigmaRefPomP");
  mPomP          = settings.parm("Diffraction:mRefPomP");
  pPomP          = settings.parm("Diffraction:mPowPomP");
  mMinPertDiff   = settings.parm("Diffraction:mMinPert");

  // Possibility to allow user veto of MPI
  canVetoMPI = (userHooksPtr != 0) ? userHooksPtr->canVetoMPIEmission()
             : false;

  // Some common combinations for double Gaussian, as shorthand.
  if (bProfile == 2) {
    fracA        = pow2(1. - coreFraction);
    fracB        = 2. * coreFraction * (1. - coreFraction);
    fracC        = pow2(coreFraction);
    radius2B     = 0.5 * (1. + pow2(coreRadius));
    radius2C     = pow2(coreRadius);

  // Some common combinations for exp(b^pow), as shorthand.
  } else if (bProfile == 3) {
    hasLowPow    = (expPow < 2.);
    expRev       = 2. / expPow - 1.;
  }

  // Initialize alpha_strong generation.
  alphaS.init( alphaSvalue, alphaSorder, alphaSnfmax, false);
  double Lambda3 = alphaS.Lambda3();

  // Initialize alphaEM generation.
  alphaEM.init( alphaEMorder, &settings);

  // Attach matrix-element calculation objects.
  sigma2gg.init( 0, processLevel, infoPtr, &settings, particleDataPtr,
    rndmPtr, beamAPtr, beamBPtr, couplingsPtr);
  sigma2qg.init( 1, processLevel, infoPtr, &settings, particleDataPtr,
    rndmPtr, beamAPtr, beamBPtr, couplingsPtr);
  sigma2qqbarSame.init( 2, processLevel, infoPtr, &settings, particleDataPtr,
    rndmPtr, beamAPtr, beamBPtr, couplingsPtr);
  sigma2qq.init( 3, processLevel, infoPtr, &settings, particleDataPtr,
    rndmPtr,  beamAPtr, beamBPtr, couplingsPtr);

  // Calculate invariant mass of system.
  eCM          = infoPtr->eCM();
  sCM          = eCM * eCM;
  mMaxPertDiff = eCM;
  eCMsave      = eCM;

  // Get the total inelastic and nondiffractive cross section.
  if (!sigmaTotPtr->hasSigmaTot()) return false;
  bool isNonDiff = (iDiffSys == 0);
  sigmaND = sigmaTotPtr->sigmaND();
  double sigmaMaxViol = 0.;

  // Output initialization info - first part.
  bool showMPI = settings.flag("Init:showMultipartonInteractions");
  if (showMPI) {
    os << "\n *-------  PYTHIA Multiparton Interactions Initialization  "
       << "---------* \n"
       << " |                                                        "
       << "          | \n";
    if (isNonDiff)
      os << " |                   sigmaNonDiffractive = " << fixed
         << setprecision(2) << setw(7) << sigmaND << " mb               | \n";
    else if (iDiffSys == 1)
      os << " |                          diffraction XB                "
         << "          | \n";
    else if (iDiffSys == 2)
      os << " |                          diffraction AX                "
         << "          | \n";
    else if (iDiffSys == 3)
      os << " |                          diffraction AXB               "
         << "          | \n";
    os << " |                                                        "
       << "          | \n";
  }

  // For diffraction need to cover range of diffractive masses.
  nStep     = (iDiffSys == 0) ? 1 : 5;
  eStepSize = (nStep < 2) ? 1.
            : log(mMaxPertDiff / mMinPertDiff) / (nStep - 1.);
  for (int iStep = 0; iStep < nStep; ++iStep) {

    // Update and output current diffractive mass and
    // fictitious Pomeron-proton cross section for normalization.
    if (nStep > 1) {
      eCM = mMinPertDiff * pow( mMaxPertDiff / mMinPertDiff,
            iStep / (nStep - 1.) );
      sCM = eCM * eCM;
      sigmaND = sigmaPomP * pow( eCM / mPomP, pPomP);
      if (showMPI) os << " |    diffractive mass = " << scientific
        << setprecision(3) << setw(9) << eCM << " GeV and sigmaNorm = "
        << fixed << setw(6) << sigmaND << " mb    | \n";
    }

    // Set current pT0 scale.
    pT0 = pT0Ref * pow(eCM / ecmRef, ecmPow);

    // The pT0 value may need to be decreased, if sigmaInt < sigmaND.
    double pT4dSigmaMaxBeg = 0.;
    for ( ; ; ) {

      // Derived pT kinematics combinations.
      pT20         = pT0*pT0;
      pT2min       = pTmin*pTmin;
      pTmax        = 0.5*eCM;
      pT2max       = pTmax*pTmax;
      pT20R        = RPT20 * pT20;
      pT20minR     = pT2min + pT20R;
      pT20maxR     = pT2max + pT20R;
      pT20min0maxR = pT20minR * pT20maxR;
      pT2maxmin    = pT2max - pT2min;

      // Provide upper estimate of interaction rate d(Prob)/d(pT2).
      upperEnvelope();

      // Setup binning in b for x-dependent matter profile.
      if (bProfile == 4) {
        sigmaIntWgt.resize(XDEP_BBIN);
        sigmaSumWgt.resize(XDEP_BBIN);
        bstepNow = XDEP_BSTEP;
      }

      // Integrate the parton-parton interaction cross section.
      pT4dSigmaMaxBeg = pT4dSigmaMax;
      jetCrossSection();

      // If the overlap-weighted cross section has not fallen below
      // cutoff, then increase bin size in b and reintegrate.
      while (bProfile == 4
        && sigmaIntWgt[XDEP_BBIN - 1] > XDEP_CUTOFF * sigmaInt) {
        bstepNow += XDEP_BSTEPINC;
        jetCrossSection();
      }

      // Sufficiently big SigmaInt or reduce pT0; maybe also pTmin.
      if (sigmaInt > SIGMASTEP * sigmaND) break;
      if (showMPI) os << fixed << setprecision(2) << " |    pT0 = "
        << setw(5) << pT0 << " gives sigmaInteraction = " << setw(7)
        << sigmaInt << " mb: rejected     | \n";
      if (pTmin > pT0) pTmin *= PT0STEP;
      pT0 *= PT0STEP;

      // Give up if pT0 and pTmin fall too low.
      if ( max(pT0, pTmin) < max(PT0MIN, Lambda3) ) {
        infoPtr->errorMsg("Error in MultipartonInteractions::init:"
          " failed to find acceptable pT0 and pTmin");
        infoPtr->setTooLowPTmin(true);
        return false;
      }
    }

    // Output for accepted pT0.
    if (showMPI) os << fixed << setprecision(2) << " |    pT0 = "
      << setw(5) << pT0 << " gives sigmaInteraction = "<< setw(7)
      << sigmaInt << " mb: accepted     | \n";

    // Calculate factor relating matter overlap and interaction rate.
    overlapInit();

    // Maximum violation relative to first estimate.
    sigmaMaxViol = max( sigmaMaxViol, pT4dSigmaMax / pT4dSigmaMaxBeg);

    // Save values calculated.
    if (nStep > 1) {
      pT0Save[iStep]          = pT0;
      pT4dSigmaMaxSave[iStep] = pT4dSigmaMax;
      pT4dProbMaxSave[iStep]  = pT4dProbMax;
      sigmaIntSave[iStep]     = sigmaInt;
      for (int j = 0; j <= 100; ++j) sudExpPTSave[iStep][j] = sudExpPT[j];
      zeroIntCorrSave[iStep]  = zeroIntCorr;
      normOverlapSave[iStep]  = normOverlap;
      kNowSave[iStep]         = kNow;
      bAvgSave[iStep]         = bAvg;
      bDivSave[iStep]         = bDiv;
      probLowBSave[iStep]     = probLowB;
      fracAhighSave[iStep]    = fracAhigh;
      fracBhighSave[iStep]    = fracBhigh;
      fracChighSave[iStep]    = fracBhigh;
      fracABChighSave[iStep]  = fracABChigh;
      cDivSave[iStep]         = cDiv;
      cMaxSave[iStep]         = cMax;
   }

  // End of loop over diffractive masses.
  }

  // Output details for x-dependent matter profile.
  if (bProfile == 4 && showMPI)
    os << " |                                              "
       << "                    | \n"
       << fixed << setprecision(2)
       << " |  x-dependent matter profile: a1 = " << a1 << ", "
       << "a0 = " << a0now * XDEP_SMB2FM << ", bStep = "
       << bstepNow << "  | \n";

  // End initialization printout.
  if (showMPI) os << " |                                              "
     << "                    | \n"
     << " *-------  End PYTHIA Multiparton Interactions Initialization"
     << "  -----* " << endl;

  // Amount of violation from upperEnvelope to jetCrossSection.
  if (sigmaMaxViol > 1.) {
    ostringstream osWarn;
    osWarn << "by factor " << fixed << setprecision(3) << sigmaMaxViol;
    infoPtr->errorMsg("Warning in MultipartonInteractions::init:"
      " maximum increased", osWarn.str());
  }

  // Reset statistics.
  SigmaMultiparton* dSigma;
  for (int i = 0; i < 4; ++i) {
    if      (i == 0) dSigma = &sigma2gg;
    else if (i == 1) dSigma = &sigma2qg;
    else if (i == 2) dSigma = &sigma2qqbarSame;
    else             dSigma = &sigma2qq;
    int nProc = dSigma->nProc();
    for (int iProc = 0; iProc < nProc; ++iProc)
      nGen[ dSigma->codeProc(iProc) ] = 0;
  }

  // Additional setup for x-dependent matter profile.
  if (bProfile == 4) {
    sigmaIntWgt.clear();
    sigmaSumWgt.clear();
  }
  // No preselection of sea/valence content and initialise a0.
  vsc1 = 0;
  vsc2 = 0;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Reset impact parameter choice and update the CM energy.
// For diffraction also interpolate parameters to current CM energy.

void MultipartonInteractions::reset( ) {

  // Reset impact parameter choice.
  bIsSet      = false;
  bSetInFirst = false;

  // Update CM energy. Done if not diffraction and not new energy.
  eCM = infoPtr->eCM();
  sCM = eCM * eCM;
  if (nStep == 1 || abs( eCM / eCMsave - 1.) < ECMDEV) return;

  // Set fictitious Pomeron-proton cross section for diffractive system.
  sigmaND = sigmaPomP * pow( eCM / mPomP, pPomP);

  // Current interpolation point.
  eCMsave   = eCM;
  eStepSave = log(eCM / mMinPertDiff) / eStepSize;
  iStepFrom = max( 0, min( nStep - 2, int( eStepSave) ) );
  iStepTo   = iStepFrom + 1;
  eStepTo   = max( 0., min( 1., eStepSave - iStepFrom) );
  eStepFrom = 1. - eStepTo;

  // Update pT0 and combinations derived from it.
  pT0           = eStepFrom * pT0Save[iStepFrom]
                + eStepTo   * pT0Save[iStepTo];
  pT20          = pT0*pT0;
  pT2min        = pTmin*pTmin;
  pTmax         = 0.5*eCM;
  pT2max        = pTmax*pTmax;
  pT20R         = RPT20 * pT20;
  pT20minR      = pT2min + pT20R;
  pT20maxR      = pT2max + pT20R;
  pT20min0maxR  = pT20minR * pT20maxR;
  pT2maxmin     = pT2max - pT2min;

  // Update other parameters used in pT choice.
  pT4dSigmaMax  = eStepFrom * pT4dSigmaMaxSave[iStepFrom]
                + eStepTo   * pT4dSigmaMaxSave[iStepTo];
  pT4dProbMax   = eStepFrom * pT4dProbMaxSave[iStepFrom]
                + eStepTo   * pT4dProbMaxSave[iStepTo];
  sigmaInt      = eStepFrom * sigmaIntSave[iStepFrom]
                + eStepTo   * sigmaIntSave[iStepTo];
  for (int j = 0; j <= 100; ++j)
    sudExpPT[j] = eStepFrom * sudExpPTSave[iStepFrom][j]
                + eStepTo   * sudExpPTSave[iStepTo][j];

  // Update parameters related to the impact-parameter picture.
  zeroIntCorr   = eStepFrom * zeroIntCorrSave[iStepFrom]
                + eStepTo   * zeroIntCorrSave[iStepTo];
  normOverlap   = eStepFrom * normOverlapSave[iStepFrom]
                + eStepTo   * normOverlapSave[iStepTo];
  kNow          = eStepFrom * kNowSave[iStepFrom]
                + eStepTo   * kNowSave[iStepTo];
  bAvg          = eStepFrom * bAvgSave[iStepFrom]
                + eStepTo   * bAvgSave[iStepTo];
  bDiv          = eStepFrom * bDivSave[iStepFrom]
                + eStepTo   * bDivSave[iStepTo];
  probLowB      = eStepFrom * probLowBSave[iStepFrom]
                + eStepTo   * probLowBSave[iStepTo];
  fracAhigh     = eStepFrom * fracAhighSave[iStepFrom]
                + eStepTo   * fracAhighSave[iStepTo];
  fracBhigh     = eStepFrom * fracBhighSave[iStepFrom]
                + eStepTo   * fracBhighSave[iStepTo];
  fracChigh     = eStepFrom * fracChighSave[iStepFrom]
                + eStepTo   * fracChighSave[iStepTo];
  fracABChigh   = eStepFrom * fracABChighSave[iStepFrom]
                + eStepTo   * fracABChighSave[iStepTo];
  cDiv          = eStepFrom * cDivSave[iStepFrom]
                + eStepTo   * cDivSave[iStepTo];
  cMax          = eStepFrom * cMaxSave[iStepFrom]
                + eStepTo   * cMaxSave[iStepTo];

}

//--------------------------------------------------------------------------

// Select first = hardest pT in nondiffractive process.
// Requires separate treatment at low and high b values.

void MultipartonInteractions::pTfirst() {

  // Pick impact parameter and thereby interaction rate enhancement.
  // This is not used for the x-dependent matter profile, which
  // instead uses trial interactions.
  if (bProfile == 4) isAtLowB = false;
  else               overlapFirst();
  bSetInFirst = true;
  double WTacc;

  // At low b values evolve downwards with Sudakov.
  if (isAtLowB) {
    pT2 = pT2max;
    do {

      // Pick a pT using a quick-and-dirty cross section estimate.
      pT2 = fastPT2(pT2);

      // If fallen below lower cutoff then need to restart at top.
      if (pT2 < pT2min) {
        pT2 = pT2max;
        WTacc = 0.;

      // Else pick complete kinematics and evaluate cross-section correction.
      } else {
        WTacc = sigmaPT2scatter(true) / dSigmaApprox;
        if (WTacc > WTACCWARN) infoPtr->errorMsg("Warning in "
            "MultipartonInteractions::pTfirst: weight above unity");
      }

    // Loop until acceptable pT and acceptable kinematics.
    } while (WTacc < rndmPtr->flat() || !dSigmaDtSel->final2KinMPI());

  // At high b values make preliminary pT choice without Sudakov factor.
  } else {

    while (true) {
      do {
        pT2 = pT20min0maxR / (pT20minR + rndmPtr->flat() * pT2maxmin) - pT20R;

        // Evaluate upper estimate of cross section for this pT2 choice.
        dSigmaApprox = pT4dSigmaMax / pow2(pT2 + pT20R);

        // Pick complete kinematics and evaluate cross-section correction.
        WTacc = sigmaPT2scatter(true) / dSigmaApprox;

        // Evaluate and include Sudakov factor.
        if (bProfile != 4) WTacc *= sudakov( pT2, enhanceB);

        // Warn for weight above unity
        if (WTacc > WTACCWARN) infoPtr->errorMsg("Warning in "
            "MultipartonInteractions::pTfirst: weight above unity");

      // Loop until acceptable pT and acceptable kinematics.
      } while (WTacc < rndmPtr->flat() || !dSigmaDtSel->final2KinMPI());

      // For x-dependent matter profile, use trial interactions to
      // generate Sudakov, otherwise done.
      if (bProfile != 4) break;
      else {
        // Save details of the original hard interaction.
        pT2Save      = pT2;
        id1Save      = id1;
        id2Save      = id2;
        x1Save       = x1;
        x2Save       = x2;
        sHatSave     = sHat;
        tHatSave     = tHat;
        uHatSave     = uHat;
        alpSsave     = alpS;
        alpEMsave    = alpEM;
        pT2FacSave   = pT2Fac;
        pT2RenSave   = pT2Ren;
        xPDF1nowSave = xPDF1now;
        xPDF2nowSave = xPDF2now;
        // Save accepted kinematics and pointer to SigmaProcess.
        dSigmaDtSel->saveKin();
        dSigmaDtSelSave = dSigmaDtSel;

        // Put x1, x2 information into beam pointers to get correct
        // PDF rescaling in trial interaction (for hard process, can
        // be sea or valence, not companion).
        beamAPtr->append( 0, id1, x1);
        beamAPtr->xfISR( 0, id1, x1, pT2Fac * pT2Fac);
        vsc1 = beamAPtr->pickValSeaComp();
        beamBPtr->append( 0, id2, x2);
        beamBPtr->xfISR( 0, id2, x2, pT2Fac * pT2Fac);
        vsc2 = beamBPtr->pickValSeaComp();

        // Pick b according to O(b, x1, x2).
        double w1    = XDEP_A1 + a1 * log(1. / x1);
        double w2    = XDEP_A1 + a1 * log(1. / x2);
        double fac   = a02now * (w1 * w1 + w2 * w2);
        double expb2 = rndmPtr->flat();
        b2now  = - fac * log(expb2);
        bNow   = sqrt(b2now);

        // Enhancement factor for the hard process and overestimate
        // for fastPT2. Note that existing framework has a (1. / sigmaND)
        // present.
        enhanceB    = sigmaND / M_PI / fac * expb2;
        enhanceBmax = sigmaND / 2. / M_PI / a02now *
                      exp( -b2now / 2. / a2max );

        // Trial interaction with dummy event.
        Event evDummy;
        double pTtrial = pTnext(pTmax, pTmin, evDummy);

        // Restore beams.
        beamAPtr->clear();
        beamBPtr->clear();

        // Accept if fallen beneath factorisation scale.
        if (pTtrial < sqrt(pT2FacSave)) {
          // Restore previous values and original sigma.
          swap(pT2,      pT2Save);
          swap(pT2Fac,   pT2FacSave);
          swap(pT2Ren,   pT2RenSave);
          swap(id1,      id1Save);
          swap(id2,      id2Save);
          swap(x1,       x1Save);
          swap(x2,       x2Save);
          swap(sHat,     sHatSave);
          swap(tHat,     tHatSave);
          swap(uHat,     uHatSave);
          swap(alpS,     alpSsave);
          swap(alpEM,    alpEMsave);
          swap(xPDF1now, xPDF1nowSave);
          swap(xPDF2now, xPDF2nowSave);
          if (dSigmaDtSel == dSigmaDtSelSave) dSigmaDtSel->swapKin();
          else swap(dSigmaDtSel, dSigmaDtSelSave);

          // Accept.
          bIsSet = true;
          break;
        }
      } // if (bProfile == 4)
    } // while (true)

  // End handling for high b.
  }

}

//--------------------------------------------------------------------------

// Set up kinematics for first = hardest pT in nondiffractive process.

void MultipartonInteractions::setupFirstSys( Event& process) {

  // Last beam-status particles. Offset relative to normal beam locations.
  int sizeProc = process.size();
  int nBeams   = 3;
  for (int i = 3; i < sizeProc; ++i)
    if (process[i].statusAbs() < 20) nBeams = i + 1;
  int nOffset  = nBeams - 3;

  // Remove any partons of previous failed interactions.
  if (sizeProc > nBeams) {
    process.popBack( sizeProc - nBeams);
    process.initColTag();
  }

  // Entries 3 and 4, now to be added, come from 1 and 2.
  process[1 + nOffset].daughter1(3 + nOffset);
  process[2 + nOffset].daughter1(4 + nOffset);

  // Negate beam status, if not already done. (Case with offset beams.)
  process[1 + nOffset].statusNeg();
  process[2 + nOffset].statusNeg();

  // Loop over four partons and offset info relative to subprocess itself.
  int colOffset = process.lastColTag();
  for (int i = 1; i <= 4; ++i) {
    Particle parton = dSigmaDtSel->getParton(i);
    if (i <= 2) parton.status(-21);
    else        parton.status(23);
    if (i <= 2) parton.mothers( i + nOffset, 0);
    else        parton.mothers( 3 + nOffset, 4 + nOffset);
    if (i <= 2) parton.daughters( 5 + nOffset, 6 + nOffset);
    else        parton.daughters( 0, 0);
    int col = parton.col();
    if (col > 0) parton.col( col + colOffset);
    int acol = parton.acol();
    if (acol > 0) parton.acol( acol + colOffset);

    // Put the partons into the event record.
    process.append(parton);
  }

  // Set scale from which to begin evolution.
  process.scale(  sqrt(pT2Fac) );

  // Info on subprocess - specific to mimimum-bias events.
  string nameSub = dSigmaDtSel->name();
  int codeSub    = dSigmaDtSel->code();
  int nFinalSub  = dSigmaDtSel->nFinal();
  double pTMPI   = dSigmaDtSel->pTMPIFin();
  infoPtr->setSubType( iDiffSys, nameSub, codeSub, nFinalSub);
  if (iDiffSys == 0) infoPtr->setTypeMPI( codeSub, pTMPI, 0, 0,
    enhanceB / zeroIntCorr);

  // Further standard info on process.
  infoPtr->setPDFalpha( iDiffSys, id1, id2, x1, x2, xPDF1now, xPDF2now,
    pT2Fac, alpEM, alpS, pT2Ren, 0.);
  double m3    = dSigmaDtSel->m(3);
  double m4    = dSigmaDtSel->m(4);
  double theta = dSigmaDtSel->thetaMPI();
  double phi   = dSigmaDtSel->phiMPI();
  infoPtr->setKin( iDiffSys, id1, id2, x1, x2, sHat, tHat, uHat, sqrt(pT2),
    m3, m4, theta, phi);

}

//--------------------------------------------------------------------------

// Find whether to limit maximum scale of emissions.

bool MultipartonInteractions::limitPTmax( Event& event) {

  // User-set cases.
  if (pTmaxMatch == 1) return true;
  if (pTmaxMatch == 2) return false;

  // Always restrict SoftQCD processes.
  if (infoPtr->isNonDiffractive() || infoPtr->isDiffractiveA()
    || infoPtr->isDiffractiveB() || infoPtr->isDiffractiveC() ) return true;

  // Look if only quarks (u, d, s, c, b), gluons and photons in final state.
  bool onlyQGP1 = true;
  bool onlyQGP2 = true;
  int  n21      = 0;
  for (int i = 5; i < event.size(); ++i) {
    if (event[i].status() == -21) ++n21;
    else if (n21 == 0) {
      int idAbs = event[i].idAbs();
      if (idAbs > 5 && idAbs != 21 && idAbs != 22) onlyQGP1 = false;
    } else if (n21 == 2) {
      int idAbs = event[i].idAbs();
      if (idAbs > 5 && idAbs != 21 && idAbs != 22) onlyQGP2 = false;
    }
  }

  // If two hard interactions then limit if one only contains q/g/gamma.
  bool onlyQGP = (n21 == 2) ? (onlyQGP1 || onlyQGP2) : onlyQGP1;
  return (onlyQGP);

}

//--------------------------------------------------------------------------

// Select next pT in downwards evolution.

double MultipartonInteractions::pTnext( double pTbegAll, double pTendAll,
  Event& event) {

  // Initial values.
  bool   pickRescatter, acceptKin;
  double dSigmaScatter, dSigmaRescatter, WTacc;
  double pT2end = pow2( max(pTmin, pTendAll) );

  // With the x-dependent matter profile, and minimum bias or diffraction,
  // it is possible to reuse values that have been stored during the trial
  // interactions for a slight speedup.
  // bIsSet is false during trial interactions, counter 21 in case partonLevel
  // is retried and counter 22 for the first pass through partonLevel.
  if (bProfile == 4 && bIsSet && bSetInFirst && infoPtr->getCounter(21) == 1
    && infoPtr->getCounter(22) == 1) {
    if (pT2Save < pT2end) return 0.;
    pT2      = pT2Save;
    pT2Fac   = pT2FacSave;
    pT2Ren   = pT2RenSave;
    id1      = id1Save;
    id2      = id2Save;
    x1       = x1Save;
    x2       = x2Save;
    sHat     = sHatSave;
    tHat     = tHatSave;
    uHat     = uHatSave;
    alpS     = alpSsave;
    alpEM    = alpEMsave;
    xPDF1now = xPDF1nowSave;
    xPDF2now = xPDF2nowSave;
    if (dSigmaDtSel == dSigmaDtSelSave) dSigmaDtSel->swapKin();
    else dSigmaDtSel = dSigmaDtSelSave;
    return sqrt(pT2);
  }

  // Do not allow rescattering while still FSR with global recoil.
  bool allowRescatterNow = allowRescatter;
  if (globalRecoilFSR && partonSystemsPtr->sizeOut(0) <= nMaxGlobalRecoilFSR)
    allowRescatterNow = false;

  // Initial pT2 value.
  pT2 = pow2(pTbegAll);

  // Find the set of already scattered partons on the two sides.
  if (allowRescatterNow) findScatteredPartons( event);

  // Pick a pT2 using a quick-and-dirty cross section estimate.
  do {
    do {
      pT2 = fastPT2(pT2);
      if (pT2 < pT2end) return 0.;

      // Initial values: no rescattering.
      i1Sel           = 0;
      i2Sel           = 0;
      dSigmaSum       = 0.;
      pickRescatter   = false;

      // Pick complete kinematics and evaluate interaction cross-section.
      dSigmaScatter   = sigmaPT2scatter(false);

      // Also cross section from rescattering if allowed.
      dSigmaRescatter = (allowRescatterNow) ? sigmaPT2rescatter( event) : 0.;

      // Normalize to dSigmaApprox, which was set in fastPT2 above.
      WTacc = (dSigmaScatter + dSigmaRescatter) / dSigmaApprox;
      if (WTacc > WTACCWARN) infoPtr->errorMsg("Warning in "
        "MultipartonInteractions::pTnext: weight above unity");

      // Idea suggested by Gosta Gustafson: increased screening in events
      // with large activity can be simulated by pT0_eff = sqrt(n) * pT0.
      if (enhanceScreening > 0) {
        int nSysNow     = infoPtr->nMPI() + 1;
        if (enhanceScreening == 2) nSysNow += infoPtr->nISR();
        double WTscreen = pow2( (pT2 + pT20) / (pT2 + nSysNow * pT20) );
        WTacc          *= WTscreen;
      }

      // x-dependent matter profile overlap weighting.
      if (bProfile == 4) {
        double w1   = XDEP_A1 + a1 * log(1. / x1);
        double w2   = XDEP_A1 + a1 * log(1. / x2);
        double fac  = a02now * (w1 * w1 + w2 * w2);
        // Correct enhancement factor and weight
        enhanceBnow = sigmaND / M_PI / fac * exp( - b2now / fac);
        double oWgt = enhanceBnow / enhanceBmax;
        if (oWgt > 1.0000000001) infoPtr->errorMsg("Warning in Multiparton"
          "Interactions::pTnext: overlap weight above unity");
        WTacc *= oWgt;
      }

      // Decide whether to keep the event based on weight.
    } while (WTacc < rndmPtr->flat());

    // When rescattering possible: new interaction or rescattering?
    if (allowRescatterNow) {
      pickRescatter = (i1Sel > 0 || i2Sel > 0);

      // Restore kinematics for selected scattering/rescattering.
      id1      = id1Sel;
      id2      = id2Sel;
      x1       = x1Sel;
      x2       = x2Sel;
      sHat     = sHatSel;
      tHat     = tHatSel;
      uHat     = uHatSel;
      sigma2Sel->sigma( id1, id2, x1, x2, sHat, tHat, uHat, alpS, alpEM,
        true, pickOtherSel);
    }

    // Pick one of the possible channels summed above.
    dSigmaDtSel = sigma2Sel->sigmaSel();
    if (sigma2Sel->swapTU()) swap( tHat, uHat);

    // Decide to keep event based on kinematics (normally always OK).
    // Rescattering: need to provide incoming four-vectors and masses.
    if (pickRescatter) {
      Vec4 p1Res = (i1Sel == 0) ? 0.5 * eCM * x1Sel * Vec4( 0., 0.,  1., 1.)
                                 : event[i1Sel].p();
      Vec4 p2Res = (i2Sel == 0) ? 0.5 * eCM * x2Sel * Vec4( 0., 0., -1., 1.)
                                   : event[i2Sel].p();
      double m1Res = (i1Sel == 0) ? 0. :  event[i1Sel].m();
      double m2Res = (i2Sel == 0) ? 0. :  event[i2Sel].m();
      acceptKin = dSigmaDtSel->final2KinMPI( i1Sel, i2Sel, p1Res, p2Res,
        m1Res, m2Res);
    // New interaction: already stored values enough.
    } else acceptKin = dSigmaDtSel->final2KinMPI();
  } while (!acceptKin);

  // Done.
  return sqrt(pT2);

}

//--------------------------------------------------------------------------

// Set up the kinematics of the 2 -> 2 scattering process,
// and store the scattering in the event record.

bool MultipartonInteractions::scatter( Event& event) {

  // Last beam-status particles. Offset relative to normal beam locations.
  int sizeProc = event.size();
  int nBeams   = 3;
  for (int i = 3; i < sizeProc; ++i)
    if (event[i].statusAbs() < 20) nBeams = i + 1;
  int nOffset  = nBeams - 3;

  // Loop over four partons and offset info relative to subprocess itself.
  int motherOffset = event.size();
  int colOffset = event.lastColTag();
  for (int i = 1; i <= 4; ++i) {
    Particle parton = dSigmaDtSel->getParton(i);
    if (i <= 2 ) parton.mothers( i + nOffset, 0);
    else parton.mothers( motherOffset, motherOffset + 1);
    if (i <= 2 ) parton.daughters( motherOffset + 2, motherOffset + 3);
    else parton.daughters( 0, 0);
    int col = parton.col();
    if (col > 0) parton.col( col + colOffset);
    int acol = parton.acol();
    if (acol > 0) parton.acol( acol + colOffset);

    // Put the partons into the event record.
    event.append(parton);
  }

  // Allow veto of MPI. If so restore event record to before scatter.
  if (canVetoMPI && userHooksPtr->doVetoMPIEmission(sizeProc, event)) {
    event.popBack(event.size() - sizeProc);
    return false;
  }

  // Store participating partons as a new set in list of all systems.
  int iSys = partonSystemsPtr->addSys();
  partonSystemsPtr->setInA(iSys, motherOffset);
  partonSystemsPtr->setInB(iSys, motherOffset + 1);
  partonSystemsPtr->addOut(iSys, motherOffset + 2);
  partonSystemsPtr->addOut(iSys, motherOffset + 3);
  partonSystemsPtr->setSHat(iSys, sHat);

  // Tag double rescattering graphs that annihilate one initial colour.
  bool annihil1 = false;
  bool annihil2 = false;
  if (i1Sel > 0 && i2Sel > 0) {
    if (event[motherOffset].col() == event[motherOffset + 1].acol()
      && event[motherOffset].col() > 0) annihil1 = true;
    if (event[motherOffset].acol() == event[motherOffset + 1].col()
      && event[motherOffset].acol() > 0) annihil2 = true;
  }

  // Beam remnant A: add scattered partons to list.
  BeamParticle& beamA = *beamAPtr;
  int iA = beamA.append( motherOffset, id1, x1);

  // Find whether incoming partons are valence or sea, so prepared for ISR.
  if (i1Sel == 0) {
    beamA.xfISR( iA, id1, x1, pT2);
    beamA.pickValSeaComp();

  // Remove rescattered parton from final state and change history.
  // Propagate existing colour labels throught graph.
  } else {
    beamA[iA].companion(-10);
    event[i1Sel].statusNeg();
    event[i1Sel].daughters( motherOffset, motherOffset);
    event[motherOffset].mothers( i1Sel, i1Sel);
    int colOld = event[i1Sel].col();
    if (colOld > 0) {
      int colNew = event[motherOffset].col();
      for (int i = motherOffset; i < motherOffset + 4; ++i) {
        if (event[i].col() == colNew) event[i].col( colOld);
        if (event[i].acol() == colNew) event[i].acol( colOld);
      }
    }
    int acolOld = event[i1Sel].acol();
    if (acolOld > 0) {
      int acolNew = event[motherOffset].acol();
      for (int i = motherOffset; i < motherOffset + 4; ++i) {
        if (event[i].col() == acolNew) event[i].col( acolOld);
        if (event[i].acol() == acolNew) event[i].acol( acolOld);
      }
    }
  }

  // Beam remnant B: add scattered partons to list.
  BeamParticle& beamB = *beamBPtr;
  int iB = beamB.append( motherOffset + 1, id2, x2);

  // Find whether incoming partons are valence or sea, so prepared for ISR.
  if (i2Sel == 0) {
    beamB.xfISR( iB, id2, x2, pT2);
    beamB.pickValSeaComp();

  // Remove rescattered parton from final state and change history.
  // Propagate existing colour labels throught graph.
  } else {
    beamB[iB].companion(-10);
    event[i2Sel].statusNeg();
    event[i2Sel].daughters( motherOffset + 1, motherOffset + 1);
    event[motherOffset + 1].mothers( i2Sel, i2Sel);
    int colOld = event[i2Sel].col();
    if (colOld > 0) {
      int colNew = event[motherOffset + 1].col();
      for (int i = motherOffset; i < motherOffset + 4; ++i) {
        if (event[i].col() == colNew) event[i].col( colOld);
        if (event[i].acol() == colNew) event[i].acol( colOld);
      }
    }
    int acolOld = event[i2Sel].acol();
    if (acolOld > 0) {
      int acolNew = event[motherOffset + 1].acol();
      for (int i = motherOffset; i < motherOffset + 4; ++i) {
        if (event[i].col() == acolNew) event[i].col( acolOld);
        if (event[i].acol() == acolNew) event[i].acol( acolOld);
      }
    }
  }

  // Annihilating colour in double rescattering requires relabelling
  // of one colour into the other in the whole preceding event.
  if (annihil1 || annihil2) {
    int colLeft = (annihil1) ? event[motherOffset].col()
                : event[motherOffset].acol();
    int mother1 = event[motherOffset].mother1();
    int mother2 = event[motherOffset + 1].mother1();
    int colLost = (annihil1)
                ? event[mother1].col() + event[mother2].acol() - colLeft
                : event[mother1].acol() + event[mother2].col() - colLeft;
    for (int i = 0; i < motherOffset; ++i) {
      if (event[i].col()  == colLost) event[i].col(  colLeft );
      if (event[i].acol() == colLost) event[i].acol( colLeft );
    }
  }

  // Store info on subprocess code and rescattered partons.
  int    codeMPI = dSigmaDtSel->code();
  double pTMPI   = dSigmaDtSel->pTMPIFin();
  if (iDiffSys == 0) infoPtr->setTypeMPI( codeMPI, pTMPI, i1Sel, i2Sel,
    enhanceBnow / zeroIntCorr);
  partonSystemsPtr->setPTHat(iSys, pTMPI);

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Determine constant in d(Prob)/d(pT2) < const / (pT2 + r * pT20)^2.

void MultipartonInteractions::upperEnvelope() {

  // Initially determine constant in jet cross section upper estimate
  // d(sigma_approx)/d(pT2) < const / (pT2 + r * pT20)^2.
  pT4dSigmaMax = 0.;

  // Loop thorough allowed pT range logarithmically evenly.
  for (int iPT = 0; iPT < 100; ++iPT) {
    double pT = pTmin * pow( pTmax/pTmin, 0.01 * (iPT + 0.5) );
    pT2       = pT*pT;
    pT2shift  = pT2 + pT20;
    pT2Ren    = pT2shift;
    pT2Fac    = (SHIFTFACSCALE) ? pT2shift : pT2;
    xT        = 2. * pT / eCM;

    // Evaluate parton density sums at x1 = x2 = xT.
    double xPDF1sumMax = (9./4.) * beamAPtr->xf(21, xT, pT2Fac);
    for (int id = 1; id <= nQuarkIn; ++id)
      xPDF1sumMax += beamAPtr->xf( id, xT, pT2Fac)
                   + beamAPtr->xf(-id, xT, pT2Fac);
    double xPDF2sumMax = (9./4.) * beamBPtr->xf(21, xT, pT2Fac);
    for (int id = 1; id <= nQuarkIn; ++id)
      xPDF2sumMax += beamBPtr->xf( id, xT, pT2Fac)
                   + beamBPtr->xf(-id, xT, pT2Fac);

    // Evaluate alpha_strong and _EM, matrix element and phase space volume.
    alpS  = alphaS.alphaS(pT2Ren);
    alpEM = alphaEM.alphaEM(pT2Ren);
    double dSigmaPartonApprox = CONVERT2MB * Kfactor * 0.5 * M_PI
      * pow2(alpS / pT2shift);
    double yMax = log(1./xT + sqrt(1./(xT*xT) - 1.));
    double volumePhSp = pow2(2. * yMax);

    // Final comparison to determine upper estimate.
    double dSigmaApproxNow = SIGMAFUDGE * xPDF1sumMax * xPDF2sumMax
      * dSigmaPartonApprox * volumePhSp;
    double pT4dSigmaNow = pow2(pT2 + pT20R) * dSigmaApproxNow;
    if ( pT4dSigmaNow > pT4dSigmaMax) pT4dSigmaMax = pT4dSigmaNow;
  }

  // Get wanted constant by dividing by the nondiffractive cross section.
  pT4dProbMax = pT4dSigmaMax / sigmaND;

}

//--------------------------------------------------------------------------

// Integrate the parton-parton interaction cross section,
// using stratified Monte Carlo sampling.
// Store result in pT bins for use as Sudakov form factors.

void MultipartonInteractions::jetCrossSection() {

  // Common factor from bin size in dpT2 / (pT2 + r * pT20)^2 and statistics.
  double sigmaFactor = (1. / pT20minR - 1. / pT20maxR) / (100. * nSample);

  // Reset overlap-weighted cross section for x-dependent matter profile.
  if (bProfile == 4) for (int bBin = 0; bBin < XDEP_BBIN; bBin++)
    sigmaIntWgt[bBin] = 0.;

  // Loop through allowed pT range evenly in dpT2 / (pT2 + r * pT20)^2.
  sigmaInt         = 0.;
  double dSigmaMax = 0.;
  sudExpPT[100]  = 0.;

  for (int iPT = 99; iPT >= 0; --iPT) {
    double sigmaSum = 0.;

    // Reset pT-binned overlap-weigted integration.
    if (bProfile == 4) for (int bBin = 0; bBin < XDEP_BBIN; bBin++)
      sigmaSumWgt[bBin] = 0.;

    // In each pT bin sample a number of random pT values.
    for (int iSample = 0; iSample < nSample; ++iSample) {
      double mappedPT2 = 1. - 0.01 * (iPT + rndmPtr->flat());
      pT2 = pT20min0maxR / (pT20minR + mappedPT2 * pT2maxmin) - pT20R;

      // Evaluate cross section dSigma/dpT2 in phase space point.
      double dSigma = sigmaPT2scatter(true);

      // Multiply by (pT2 + r * pT20)^2 to compensate for pT sampling. Sum.
      dSigma   *= pow2(pT2 + pT20R);
      sigmaSum += dSigma;
      if (dSigma > dSigmaMax) dSigmaMax = dSigma;

      // Overlap-weighted cross section for x-dependent matter profile.
      // Note that dSigma can be 0. when points are rejected.
      if (bProfile == 4 && dSigma > 0.) {
        double w1  = XDEP_A1 + a1 * log(1. / x1);
        double w2  = XDEP_A1 + a1 * log(1. / x2);
        double fac = XDEP_A0 * XDEP_A0 * (w1 * w1 + w2 * w2);
        double b   = 0.5 * bstepNow;
        for (int bBin = 0; bBin < XDEP_BBIN; bBin++) {
          double wgt = exp( - b * b / fac ) / fac / M_PI;
          sigmaSumWgt[bBin] += dSigma * wgt;
          b += bstepNow;
        }
      }
    }

    // Store total cross section and exponent of Sudakov.
    sigmaSum *= sigmaFactor;
    sigmaInt += sigmaSum;
    sudExpPT[iPT] = sudExpPT[iPT + 1] + sigmaSum / sigmaND;

    // Sum overlap-weighted cross section
    if (bProfile == 4) for (int bBin = 0; bBin < XDEP_BBIN; bBin++) {
      sigmaSumWgt[bBin] *= sigmaFactor;
      sigmaIntWgt[bBin] += sigmaSumWgt[bBin];
    }

  // End of loop over pT values.
  }

  // Update upper estimate of differential cross section. Done.
  if (dSigmaMax  > pT4dSigmaMax) {
    pT4dSigmaMax = dSigmaMax;
    pT4dProbMax  = dSigmaMax / sigmaND;
  }

}

//--------------------------------------------------------------------------

// Evaluate "Sudakov form factor" for not having a harder interaction
// at the selected b value, given the pT scale of the event.

double MultipartonInteractions::sudakov(double pT2sud, double enhance) {

  // Find bin the pT2 scale falls in.
  double xBin = (pT2sud - pT2min) * pT20maxR
    / (pT2maxmin * (pT2sud + pT20R));
  xBin = max(1e-6, min(100. - 1e-6, 100. * xBin) );
  int iBin = int(xBin);

  // Interpolate inside bin. Optionally include enhancement factor.
  double sudExp = sudExpPT[iBin] + (xBin - iBin)
    * (sudExpPT[iBin + 1] - sudExpPT[iBin]);
  return exp( -enhance * sudExp);

}

//--------------------------------------------------------------------------

// Pick a trial next pT, based on a simple upper estimate of the
// d(sigma)/d(pT2) spectrum.

double MultipartonInteractions::fastPT2( double pT2beg) {

  // Use d(Prob)/d(pT2) < pT4dProbMax / (pT2 + r * pT20)^2.
  double pT20begR       = pT2beg + pT20R;
  double pT4dProbMaxNow = pT4dProbMax * enhanceBmax;
  double pT2try         = pT4dProbMaxNow * pT20begR
    / (pT4dProbMaxNow - pT20begR * log(rndmPtr->flat())) - pT20R;

  // Save cross section associated with ansatz above. Done.
  dSigmaApprox = pT4dSigmaMax / pow2(pT2try + pT20R);
  return pT2try;

}

//--------------------------------------------------------------------------

// Calculate the actual cross section to decide whether fast choice is OK.
// Select flavours and kinematics for interaction at given pT.
// Slightly different treatment for first interaction and subsequent ones.

double MultipartonInteractions::sigmaPT2scatter(bool isFirst) {

  // Derive recormalization and factorization scales, amd alpha_strong/em.
  pT2shift = pT2 + pT20;
  pT2Ren   = pT2shift;
  pT2Fac   = (SHIFTFACSCALE) ? pT2shift : pT2;
  alpS  = alphaS.alphaS(pT2Ren);
  alpEM = alphaEM.alphaEM(pT2Ren);

  // Derive rapidity limits from chosen pT2.
  xT       = 2. * sqrt(pT2) / eCM;
  if (xT >= 1.) return 0.;
  xT2      = xT*xT;
  double yMax = log(1./xT + sqrt(1./xT2 - 1.));

  // Select rapidities y3 and y4 of the two produced partons.
  double y3 = yMax * (2. * rndmPtr->flat() - 1.);
  double y4 = yMax * (2. * rndmPtr->flat() - 1.);
  y = 0.5 * (y3 + y4);

  // Failure if x1 or x2 exceed what is left in respective beam.
  x1 = 0.5 * xT * (exp(y3) + exp(y4));
  x2 = 0.5 * xT * (exp(-y3) + exp(-y4));
  if (isFirst && iDiffSys == 0) {
    if (x1 > 1. || x2 > 1.) return 0.;
  } else {
    if (x1 > beamAPtr->xMax() || x2 > beamBPtr->xMax()) return 0.;
  }
  tau = x1 * x2;

  // Begin evaluate parton densities at actual x1 and x2.
  double xPDF1[21];
  double xPDF1sum = 0.;
  double xPDF2[21];
  double xPDF2sum = 0.;

  // For first interaction use normal densities.
  if (isFirst) {
    for (int id = -nQuarkIn; id <= nQuarkIn; ++id) {
      if (id == 0) xPDF1[10] = (9./4.) * beamAPtr->xf(21, x1, pT2Fac);
      else xPDF1[id+10] = beamAPtr->xf(id, x1, pT2Fac);
      xPDF1sum += xPDF1[id+10];
    }
    for (int id = -nQuarkIn; id <= nQuarkIn; ++id) {
      if (id == 0) xPDF2[10] = (9./4.) * beamBPtr->xf(21, x2, pT2Fac);
      else xPDF2[id+10] = beamBPtr->xf(id, x2, pT2Fac);
      xPDF2sum += xPDF2[id+10];
    }

  // For subsequent interactions use rescaled densities.
  } else {
    for (int id = -nQuarkIn; id <= nQuarkIn; ++id) {
      if (id == 0) xPDF1[10] = (9./4.) * beamAPtr->xfMPI(21, x1, pT2Fac);
      else xPDF1[id+10] = beamAPtr->xfMPI(id, x1, pT2Fac);
      xPDF1sum += xPDF1[id+10];
    }
    for (int id = -nQuarkIn; id <= nQuarkIn; ++id) {
      if (id == 0) xPDF2[10] = (9./4.) * beamBPtr->xfMPI(21, x2, pT2Fac);
      else xPDF2[id+10] = beamBPtr->xfMPI(id, x2, pT2Fac);
      xPDF2sum += xPDF2[id+10];
    }
  }

  // Select incoming flavours according to actual PDF's.
  id1 = -nQuarkIn - 1;
  double temp = xPDF1sum * rndmPtr->flat();
  do { xPDF1now = xPDF1[(++id1) + 10]; temp -= xPDF1now; }
  while (temp > 0. && id1 < nQuarkIn);
  if (id1 == 0) id1 = 21;
  id2 = -nQuarkIn-1;
  temp = xPDF2sum * rndmPtr->flat();
  do { xPDF2now = xPDF2[(++id2) + 10]; temp -= xPDF2now;}
  while (temp > 0. && id2 < nQuarkIn);
  if (id2 == 0) id2 = 21;

  // Assign pointers to processes relevant for incoming flavour choice:
  // g + g, q + g, q + qbar (same flavour), q + q(bar) (the rest).
  // Factor 4./9. per incoming gluon to compensate for preweighting.
  SigmaMultiparton* sigma2Tmp;
  double gluFac = 1.;
  if (id1 == 21 && id2 == 21) {
    sigma2Tmp = &sigma2gg;
    gluFac = 16. / 81.;
  } else if (id1 == 21 || id2 == 21) {
    sigma2Tmp = &sigma2qg;
    gluFac = 4. / 9.;
  } else if (id1 == -id2) sigma2Tmp = &sigma2qqbarSame;
  else sigma2Tmp = &sigma2qq;

  // Prepare to generate differential cross sections.
  sHat        = tau * sCM;
  double root = sqrtpos(1. - xT2 / tau);
  tHat        = -0.5 * sHat * (1. - root);
  uHat        = -0.5 * sHat * (1. + root);

  // Evaluate cross sections, include possibility of K factor.
  // Use kinematics for two symmetrical configurations (tHat <-> uHat).
  // (Not necessary, but makes for better MC efficiency.)
  double dSigmaPartonCorr = Kfactor * gluFac
    * sigma2Tmp->sigma( id1, id2, x1, x2, sHat, tHat, uHat, alpS, alpEM);

  // Combine cross section, pdf's and phase space integral.
  double volumePhSp = pow2(2. * yMax);
  double dSigmaScat = dSigmaPartonCorr * xPDF1sum * xPDF2sum * volumePhSp;

  // Dampen cross section at small pT values; part of formalism.
  // Use pT2 corrected for massive kinematics at this step.??
  // double pT2massive = dSigmaDtSel->pT2MPI();
  double pT2massive = pT2;
  dSigmaScat *= pow2( pT2massive / (pT2massive + pT20) );

  // Sum up total contribution for all scatterings and rescatterings.
  dSigmaSum  += dSigmaScat;

  // Save values for comparison with rescattering processes.
  i1Sel        = 0;
  i2Sel        = 0;
  id1Sel       = id1;
  id2Sel       = id2;
  x1Sel        = x1;
  x2Sel        = x2;
  sHatSel      = sHat;
  tHatSel      = tHat;
  uHatSel      = uHat;
  sigma2Sel    = sigma2Tmp;
  pickOtherSel = sigma2Tmp->pickedOther();

  // For first interaction: pick one of the possible channels summed above.
  if (isFirst) {
    dSigmaDtSel = sigma2Tmp->sigmaSel();
    if (sigma2Tmp->swapTU()) swap( tHat, uHat);
  }

  // Done.
  return dSigmaScat;
}

//--------------------------------------------------------------------------

// Find the partons that are allowed to rescatter.

void MultipartonInteractions::findScatteredPartons( Event& event) {

  // Reset arrays.
  scatteredA.resize(0);
  scatteredB.resize(0);
  double yTmp, probA, probB;

  // Loop though the event record and catch "final" partons.
  for (int i = 0; i < event.size(); ++i)
  if (event[i].isFinal() && (event[i].idAbs() <= nQuarkIn
  || event[i].id() == 21)) {
    yTmp = event[i].y();

    // Different strategies to determine which partons may rescatter.
    switch(rescatterMode) {

    // Case 0: step function at origin
    case 0:
      if ( yTmp > 0.) scatteredA.push_back( i);
      if (-yTmp > 0.) scatteredB.push_back( i);
      break;

    // Case 1: step function as ySepResc.
    case 1:
      if ( yTmp > ySepResc) scatteredA.push_back( i);
      if (-yTmp > ySepResc) scatteredB.push_back( i);
      break;

    // Case 2: linear rise from ySep - deltaY to ySep + deltaY.
    case 2:
      probA = 0.5 * (1. + ( yTmp - ySepResc) / deltaYResc);
      if (probA > rndmPtr->flat()) scatteredA.push_back( i);
      probB = 0.5 * (1. + (-yTmp - ySepResc) / deltaYResc);
      if (probB > rndmPtr->flat()) scatteredB.push_back( i);
      break;

    // Case 3: rise like (1/2) * ( 1 + tanh((y - ySep) / deltaY) ).
    // Use that (1/2) (1 + tanh(x)) = 1 / (1 + exp(-2x)).
    case 3:
      probA = 1. / (1. + exp(-2. * ( yTmp - ySepResc) / deltaYResc));
      if (probA > rndmPtr->flat()) scatteredA.push_back( i);
      probB = 1. / (1. + exp(-2. * (-yTmp - ySepResc) / deltaYResc));
      if (probB > rndmPtr->flat()) scatteredB.push_back( i);
      break;

    // Case 4 and undefined values: all partons can rescatter.
    default:
      scatteredA.push_back( i);
      scatteredB.push_back( i);
      break;

    // End of loop over partons. Done.
    }
  }

}

//--------------------------------------------------------------------------

// Rescattering contribution summed over all already scattered partons.
// Calculate the actual cross section to decide whether fast choice is OK.
// Select flavours and kinematics for interaction at given pT.

double MultipartonInteractions::sigmaPT2rescatter( Event& event) {

  // Derive xT scale from chosen pT2.
  xT       = 2. * sqrt(pT2) / eCM;
  if (xT >= 1.) return 0.;
  xT2      = xT*xT;

  //  Pointer to cross section and total rescatter contribution.
  SigmaMultiparton* sigma2Tmp;
  double dSigmaResc = 0.;

  // Normally save time by picking one random scattered parton from side A.
  int nScatA = scatteredA.size();
  int iScatA = -1;
  if (PREPICKRESCATTER && nScatA > 0) iScatA = min( nScatA - 1,
    int( rndmPtr->flat() * double(nScatA) ) );

  // Loop over all already scattered partons from side A.
  for (int iScat = 0; iScat < nScatA; ++iScat) {
    if (PREPICKRESCATTER && iScat != iScatA) continue;
    int iA         = scatteredA[iScat];
    int id1Tmp     = event[iA].id();

    // Calculate maximum possible sHat and check whether big enough.
    double x1Tmp   = event[iA].pPos() / eCM;
    double sHatMax = x1Tmp * beamBPtr->xMax() * sCM;
    if (sHatMax < 4. * pT2) continue;

    // Select rapidity separation between two produced partons.
    double dyMax   = log( sqrt(0.25 * sHatMax / pT2)
                   * (1. + sqrtpos(1. - 4. * pT2 / sHatMax)) );
    double dy      = dyMax * (2. * rndmPtr->flat() - 1.);

    // Reconstruct kinematical variables, especially x2.
    // Incoming c/b masses handled better if tau != x1 * x2.
    double m2Tmp   = event[iA].m2();
    double sHatTmp = m2Tmp + 4. * pT2 * pow2(cosh(dy));
    double x2Tmp   = (sHatTmp - m2Tmp) /(x1Tmp * sCM);
    double tauTmp  = sHatTmp / sCM;
    double root    = sqrtpos(1. - xT2 / tauTmp);
    double tHatTmp = -0.5 * sHatTmp * (1. - root);
    double uHatTmp = -0.5 * sHatTmp * (1. + root);

    // Begin evaluate parton densities at actual x2.
    double xPDF2[21];
    double xPDF2sum = 0.;

    // Use rescaled densities, with preweighting 9/4 for gluons.
    for (int id = -nQuarkIn; id <= nQuarkIn; ++id) {
      if (id == 0) xPDF2[10] = (9./4.) * beamBPtr->xfMPI(21, x2Tmp, pT2Fac);
      else xPDF2[id+10] = beamBPtr->xfMPI(id, x2Tmp, pT2Fac);
      xPDF2sum += xPDF2[id+10];
    }

    // Select incoming flavour according to actual PDF's.
    int id2Tmp = -nQuarkIn - 1;
    double temp = xPDF2sum * rndmPtr->flat();
    do { xPDF2now = xPDF2[(++id2Tmp) + 10]; temp -= xPDF2now;}
    while (temp > 0. && id2Tmp < nQuarkIn);
    if (id2Tmp == 0) id2Tmp = 21;

    // Assign pointers to processes relevant for incoming flavour choice:
    // g + g, q + g, q + qbar (same flavour), q + q(bar) (the rest).
    // Factor 4./9. for incoming gluon 2 to compensate for preweighting.
    if      (id1Tmp == 21 && id2Tmp == 21) sigma2Tmp = &sigma2gg;
    else if (id1Tmp == 21 || id2Tmp == 21) sigma2Tmp = &sigma2qg;
    else if (id1Tmp == -id2Tmp)            sigma2Tmp = &sigma2qqbarSame;
    else                                   sigma2Tmp = &sigma2qq;
    double gluFac = (id2Tmp == 21) ? 4. / 9. : 1.;

    // Evaluate cross sections, include possibility of K factor.
    // Use kinematics for two symmetrical configurations (tHat <-> uHat).
    // (Not necessary, but makes for better MC efficiency.)
    double dSigmaPartonCorr = Kfactor * gluFac
      * sigma2Tmp->sigma( id1Tmp, id2Tmp, x1Tmp, x2Tmp, sHatTmp, tHatTmp,
        uHatTmp, alpS, alpEM);

    // Combine cross section, pdf's and phase space integral.
    double volumePhSp = 4. *dyMax;
    double dSigmaCorr = dSigmaPartonCorr * xPDF2sum * volumePhSp;

    // Dampen cross section at small pT values; part of formalism.
    // Use pT2 corrected for massive kinematics at this step.
    //?? double pT2massive = dSigmaDtSel->pT2MPI();
    double pT2massive = pT2;
    dSigmaCorr *= pow2( pT2massive / (pT2massive + pT20) );

    // Compensate for prepicked rescattering if appropriate.
    if (PREPICKRESCATTER) dSigmaCorr *= nScatA;

    // Sum up total contribution for all scatterings or rescattering only.
    dSigmaSum  += dSigmaCorr;
    dSigmaResc += dSigmaCorr;

    // Determine whether current rescattering should be selected.
    if (dSigmaCorr > rndmPtr->flat() * dSigmaSum) {
      i1Sel        = iA;
      i2Sel        = 0;
      id1Sel       = id1Tmp;
      id2Sel       = id2Tmp;
      x1Sel        = x1Tmp;
      x2Sel        = x2Tmp;
      sHatSel      = sHatTmp;
      tHatSel      = tHatTmp;
      uHatSel      = uHatTmp;
      sigma2Sel    = sigma2Tmp;
      pickOtherSel = sigma2Tmp->pickedOther();
    }
  }

  // Normally save time by picking one random scattered parton from side B.
  int nScatB = scatteredB.size();
  int iScatB = -1;
  if (PREPICKRESCATTER && nScatB > 0) iScatB = min( nScatB - 1,
    int( rndmPtr->flat() * double(nScatB) ) );

  // Loop over all already scattered partons from side B.
  for (int iScat = 0; iScat < nScatB; ++iScat) {
    if (PREPICKRESCATTER && iScat != iScatB) continue;
    int iB         = scatteredB[iScat];
    int id2Tmp     = event[iB].id();

    // Calculate maximum possible sHat and check whether big enough.
    double x2Tmp   = event[iB].pNeg() / eCM;
    double sHatMax = beamAPtr->xMax() * x2Tmp * sCM;
    if (sHatMax < 4. * pT2) continue;

    // Select rapidity separation between two produced partons.
    double dyMax   = log( sqrt(0.25 * sHatMax / pT2)
                   * (1. + sqrtpos(1. - 4. * pT2 / sHatMax)) );
    double dy      = dyMax * (2. * rndmPtr->flat() - 1.);

    // Reconstruct kinematical variables, especially x1.
    // Incoming c/b masses handled better if tau != x1 * x2.
    double m2Tmp   = event[iB].m2();
    double sHatTmp = m2Tmp + 4. * pT2 * pow2(cosh(dy));
    double x1Tmp   = (sHatTmp - m2Tmp) /(x2Tmp * sCM);
    double tauTmp  = sHatTmp / sCM;
    double root    = sqrtpos(1. - xT2 / tauTmp);
    double tHatTmp = -0.5 * sHatTmp * (1. - root);
    double uHatTmp = -0.5 * sHatTmp * (1. + root);

    // Begin evaluate parton densities at actual x1.
    double xPDF1[21];
    double xPDF1sum = 0.;

    // Use rescaled densities, with preweighting 9/4 for gluons.
    for (int id = -nQuarkIn; id <= nQuarkIn; ++id) {
      if (id == 0) xPDF1[10] = (9./4.) * beamAPtr->xfMPI(21, x1Tmp, pT2Fac);
      else xPDF1[id+10] = beamAPtr->xfMPI(id, x1Tmp, pT2Fac);
      xPDF1sum += xPDF1[id+10];
    }

    // Select incoming flavour according to actual PDF's.
    int id1Tmp = -nQuarkIn - 1;
    double temp = xPDF1sum * rndmPtr->flat();
    do { xPDF1now = xPDF1[(++id1Tmp) + 10]; temp -= xPDF1now;}
    while (temp > 0. && id1Tmp < nQuarkIn);
    if (id1Tmp == 0) id1Tmp = 21;

    // Assign pointers to processes relevant for incoming flavour choice:
    // g + g, q + g, q + qbar (same flavour), q + q(bar) (the rest).
    // Factor 4./9. for incoming gluon 2 to compensate for preweighting.
    if      (id1Tmp == 21 && id2Tmp == 21) sigma2Tmp = &sigma2gg;
    else if (id1Tmp == 21 || id2Tmp == 21) sigma2Tmp = &sigma2qg;
    else if (id1Tmp == -id2Tmp)            sigma2Tmp = &sigma2qqbarSame;
    else                                   sigma2Tmp = &sigma2qq;
    double gluFac = (id1Tmp == 21) ? 4. / 9. : 1.;

    // Evaluate cross sections, include possibility of K factor.
    // Use kinematics for two symmetrical configurations (tHat <-> uHat).
    // (Not necessary, but makes for better MC efficiency.)
    double dSigmaPartonCorr = Kfactor * gluFac
      * sigma2Tmp->sigma( id1Tmp, id2Tmp, x1Tmp, x2Tmp, sHatTmp, tHatTmp,
        uHatTmp, alpS, alpEM);

    // Combine cross section, pdf's and phase space integral.
    double volumePhSp = 4. *dyMax;
    double dSigmaCorr = dSigmaPartonCorr * xPDF1sum * volumePhSp;

    // Dampen cross section at small pT values; part of formalism.
    // Use pT2 corrected for massive kinematics at this step.
    //?? double pT2massive = dSigmaDtSel->pT2MPI();
    double pT2massive = pT2;
    dSigmaCorr *= pow2( pT2massive / (pT2massive + pT20) );

    // Compensate for prepicked rescattering if appropriate.
    if (PREPICKRESCATTER) dSigmaCorr *= nScatB;

    // Sum up total contribution for all scatterings or rescattering only.
    dSigmaSum  += dSigmaCorr;
    dSigmaResc += dSigmaCorr;

    // Determine whether current rescattering should be selected.
    if (dSigmaCorr > rndmPtr->flat() * dSigmaSum) {
      i1Sel        = 0;
      i2Sel        = iB;
      id1Sel       = id1Tmp;
      id2Sel       = id2Tmp;
      x1Sel        = x1Tmp;
      x2Sel        = x2Tmp;
      sHatSel      = sHatTmp;
      tHatSel      = tHatTmp;
      uHatSel      = uHatTmp;
      sigma2Sel    = sigma2Tmp;
      pickOtherSel = sigma2Tmp->pickedOther();
    }
  }

  // Double rescatter: loop over already scattered partons from both sides.
  // As before, allow to use only one parton per side to speed up.
  if (allowDoubleRes) {
    for (int iScat1 = 0; iScat1 < nScatA; ++iScat1) {
      if (PREPICKRESCATTER && iScat1 != iScatA) continue;
      int iA           = scatteredA[iScat1];
      int id1Tmp       = event[iA].id();
      for (int iScat2 = 0; iScat2 < nScatB; ++iScat2) {
        if (PREPICKRESCATTER && iScat2 != iScatB) continue;
        int iB         = scatteredB[iScat2];
        int id2Tmp     = event[iB].id();

        // Calculate current sHat and check whether big enough.
        double sHatTmp = (event[iA].p() + event[iB].p()).m2Calc();
        if (sHatTmp < 4. * pT2) continue;

        // Reconstruct other kinematical variables. (x values not needed.)
        double x1Tmp   = event[iA].pPos() / eCM;
        double x2Tmp   = event[iB].pNeg() / eCM;
        double tauTmp  = sHatTmp / sCM;
        double root    = sqrtpos(1. - xT2 / tauTmp);
        double tHatTmp = -0.5 * sHatTmp * (1. - root);
        double uHatTmp = -0.5 * sHatTmp * (1. + root);

        // Assign pointers to processes relevant for incoming flavour choice:
        // g + g, q + g, q + qbar (same flavour), q + q(bar) (the rest).
        if      (id1Tmp == 21 && id2Tmp == 21) sigma2Tmp = &sigma2gg;
        else if (id1Tmp == 21 || id2Tmp == 21) sigma2Tmp = &sigma2qg;
        else if (id1Tmp == -id2Tmp)            sigma2Tmp = &sigma2qqbarSame;
        else                                   sigma2Tmp = &sigma2qq;

        // Evaluate cross sections, include possibility of K factor.
        // Use kinematics for two symmetrical configurations (tHat <-> uHat).
        // (Not necessary, but makes for better MC efficiency.)
        double dSigmaPartonCorr = Kfactor
          * sigma2Tmp->sigma( id1Tmp, id2Tmp, x1Tmp, x2Tmp, sHatTmp, tHatTmp,
            uHatTmp, alpS, alpEM);

        // Combine cross section and Jacobian tHat -> pT2
        // (with safety minvalue).
        double dSigmaCorr = dSigmaPartonCorr / max(ROOTMIN, root);

        // Dampen cross section at small pT values; part of formalism.
        // Use pT2 corrected for massive kinematics at this step.
        //?? double pT2massive = dSigmaDtSel->pT2MPI();
        double pT2massive = pT2;
        dSigmaCorr *= pow2( pT2massive / (pT2massive + pT20) );

        // Compensate for prepicked rescattering if appropriate.
        if (PREPICKRESCATTER) dSigmaCorr *= nScatA * nScatB;

        // Sum up total contribution for all scatterings or rescattering only.
        dSigmaSum  += dSigmaCorr;
        dSigmaResc += dSigmaCorr;

        // Determine whether current rescattering should be selected.
        if (dSigmaCorr > rndmPtr->flat() * dSigmaSum) {
          i1Sel        = iA;
          i2Sel        = iB;
          id1Sel       = id1Tmp;
          id2Sel       = id2Tmp;
          x1Sel        = x1Tmp;
          x2Sel        = x2Tmp;
          sHatSel      = sHatTmp;
          tHatSel      = tHatTmp;
          uHatSel      = uHatTmp;
          sigma2Sel    = sigma2Tmp;
          pickOtherSel = sigma2Tmp->pickedOther();
        }
      }
    }
  }

  // Done.
  return dSigmaResc;
}


//--------------------------------------------------------------------------

// Calculate factor relating matter overlap and interaction rate,
// i.e. k in <n_interaction(b)> = k * overlap(b) (neglecting
// n_int = 0 corrections and energy-momentum conservation effects).

void MultipartonInteractions::overlapInit() {

  // Initial values for iteration. Step size of b integration.
  nAvg = sigmaInt / sigmaND;
  kNow = 0.5;
  int stepDir = 1;
  double deltaB = BSTEP;
  if (bProfile == 2) deltaB *= min( 0.5, 2.5 * coreRadius);
  if (bProfile == 3) deltaB *= max(1., pow(2. / expPow, 1. / expPow));

  // Further variables, with dummy initial values.
  double nNow           = 0.;
  double kLow           = 0.;
  double nLow           = 0.;
  double kHigh          = 0.;
  double nHigh          = 0.;
  double overlapNow     = 0.;
  double probNow        = 0.;
  double overlapInt     = 0.5;
  double probInt        = 0.;
  double probOverlapInt = 0.;
  double bProbInt       = 0.;
  normPi                = 1. / (2. * M_PI);

  // Subdivision into low-b and high-b region by interaction rate.
  bool pastBDiv = false;
  double overlapHighB = 0.;

  // For x-dependent matter profile, try to find a0 rather than k.
  // Existing framework is used, but with substitutions:
  //   a0 tuned according to Int( Pint(b), d^2b ) = sigmaND,
  //   nAvg = sigmaND, kNow = a0now, kLow = a0low, kHigh = a0high,
  //   nNow = probInt, nLow = probIntLow, nHigh = probIntHigh.
  double rescale2 = 1.;
  if (bProfile == 4) {
    nAvg = sigmaND;
    kNow = XDEP_A0 / 2.0;
  }

  // First close k into an interval by binary steps,
  // then find k by successive interpolation.
  do {
    if (stepDir == 1) kNow *= 2.;
    else if (stepDir == -1) kNow *= 0.5;
    else kNow = kLow + (nAvg - nLow) * (kHigh - kLow) / (nHigh - nLow);

    // Overlap trivial if no impact parameter dependence.
    if (bProfile <= 0 || bProfile > 4) {
      probInt        = 0.5 * M_PI * (1. - exp(-kNow));
      probOverlapInt = probInt / M_PI;
      bProbInt       = probInt;

      // Ratio of b-integrated k * overlap / (1 - exp( - k * overlap)).
      nNow = M_PI * kNow * overlapInt / probInt;

    // Else integrate overlap over impact parameter.
    } else if (bProfile < 4) {

      // Reset integrals.
      overlapInt     = (bProfile == 3) ? 0. : 0.5;
      probInt        = 0.;
      probOverlapInt = 0.;
      bProbInt       = 0.;
      pastBDiv       = false;
      overlapHighB   = 0.;

      // Step in b space.
      double b       = -0.5 * deltaB;
      double bArea   = 0.;
      do {
        b           += deltaB;
        bArea        = 2. * M_PI * b * deltaB;

        // Evaluate overlap at current b value.
        if (bProfile == 1) {
          overlapNow = normPi * exp( -b*b);
        } else if (bProfile == 2) {
          overlapNow = normPi * ( fracA * exp( -min(EXPMAX, b*b))
            + fracB * exp( -min(EXPMAX, b*b / radius2B)) / radius2B
            + fracC * exp( -min(EXPMAX, b*b / radius2C)) / radius2C );
        } else {
          overlapNow  = normPi * exp( -pow( b, expPow));
          overlapInt += bArea * overlapNow;
        }
        if (pastBDiv) overlapHighB += bArea * overlapNow;

        // Calculate interaction probability and integrate.
        probNow         = 1. - exp( -min(EXPMAX, M_PI * kNow * overlapNow));
        probInt        += bArea * probNow;
        probOverlapInt += bArea * overlapNow * probNow;
        bProbInt       += b * bArea * probNow;

        // Check when interaction probability has dropped sufficiently.
        if (!pastBDiv && probNow < PROBATLOWB) {
          bDiv     = b + 0.5 * deltaB;
          pastBDiv = true;
        }

      // Continue out in b until overlap too small.
      } while (b < 1. || b * probNow > BMAX);

      // Ratio of b-integrated k * overlap / (1 - exp( - k * overlap)).
      nNow = M_PI * kNow * overlapInt / probInt;

    // x-dependent matter profile.
    } else if (bProfile == 4) {
      rescale2  = pow2(kNow / XDEP_A0);
      probInt   = 0.;
      double b  = 0.5 * bstepNow;
      for (int bBin = 0; bBin < XDEP_BBIN; bBin++) {
        double bArea   = 2. * M_PI * b * bstepNow;
        double pIntNow = 1 - exp( -min(EXPMAX, sigmaIntWgt[bBin] / rescale2) );
        probInt += bArea * rescale2 * pIntNow;
        b += bstepNow;
      }
      nNow = probInt;
    }

    // Replace lower or upper limit of k.
    if (nNow < nAvg) {
      kLow = kNow;
      nLow = nNow;
      if (stepDir == -1) stepDir = 0;
    } else {
      kHigh = kNow;
      nHigh = nNow;
      if (stepDir == 1) stepDir = -1;
    }

  // Continue iteration until convergence.
  } while (abs(nNow - nAvg) > KCONVERGE * nAvg);

  // Save relevant final numbers for overlap values.
  if (bProfile >= 0 && bProfile < 4) {
    double avgOverlap  = probOverlapInt / probInt;
    zeroIntCorr = probOverlapInt / overlapInt;
    normOverlap = normPi * zeroIntCorr / avgOverlap;
    bAvg = bProbInt / probInt;

  // Values for x-dependent matter profile.
  } else if (bProfile == 4) {
    // bAvg        = Int(b * Pint(b), d2b)      / sigmaND.
    // zeroIntCorr = Int(<n(b)> * Pint(b), d2b) / sigmaInt.
    bAvg        = 0.;
    zeroIntCorr = 0.;
    double b    = 0.5 * bstepNow;
    for (int bBin = 0; bBin < XDEP_BBIN; bBin++) {
      double bArea   = 2. * M_PI * b * bstepNow;
      double pIntNow = 1. - exp( -min(EXPMAX, sigmaIntWgt[bBin] / rescale2) );
      bAvg          += sqrt(rescale2) * b * bArea * rescale2 * pIntNow;
      zeroIntCorr   += bArea * sigmaIntWgt[bBin] * pIntNow;
      b             += bstepNow;
    }
    bAvg        /= nNow;
    zeroIntCorr /= sigmaInt;

    // Other required values.
    a0now  = kNow;
    infoPtr->seta0MPI(a0now * XDEP_SMB2FM);
    a02now = a0now * a0now;
    double xMin = 2. * pTmin / infoPtr->eCM();
    a2max  = a0now * (XDEP_A1 + a1 * log(1. / xMin));
    a2max *= a2max;
  }

  // Relative rates for preselection of low-b and high-b region.
  // Other useful combinations for subsequent selection.
  if (bProfile > 0 && bProfile <= 3) {
    probLowB = M_PI * bDiv*bDiv;
    double probHighB = M_PI * kNow * overlapHighB;
    if (bProfile == 1) probHighB = M_PI * kNow * 0.5 * exp( -bDiv*bDiv);
    else if (bProfile == 2) {
      fracAhigh   = fracA * exp( -bDiv*bDiv);
      fracBhigh   = fracB * exp( -bDiv*bDiv / radius2B);
      fracChigh   = fracC * exp( -bDiv*bDiv / radius2C);
      fracABChigh = fracAhigh + fracBhigh + fracChigh;
      probHighB   = M_PI * kNow * 0.5 * fracABChigh;
    } else {
      cDiv = pow( bDiv, expPow);
      cMax = max(2. * expRev, cDiv);
    }
    probLowB /= (probLowB + probHighB);
  }

}

//--------------------------------------------------------------------------

// Pick impact parameter and interaction rate enhancement beforehand,
// i.e. before even the hardest interaction for minimum-bias events.

void MultipartonInteractions::overlapFirst() {

  // Trivial values if no impact parameter dependence.
  if (bProfile <= 0 || bProfile > 4) {
    bNow     = bAvg;
    enhanceB = zeroIntCorr;
    bIsSet   = true;
    isAtLowB = true;
    return;
  }

  // Preliminary choice between and inside low-b and high-b regions.
  double overlapNow = 0.;
  double probAccept = 0.;
  do {

    // Treatment in low-b region: pick b flat in area.
    if (rndmPtr->flat() < probLowB) {
      isAtLowB = true;
      bNow = bDiv * sqrt(rndmPtr->flat());

      // Evaluate overlap and from that acceptance probability.
      if (bProfile == 1) overlapNow = normPi * exp( -bNow*bNow);
      else if (bProfile == 2) overlapNow = normPi *
        ( fracA * exp( -bNow*bNow)
        + fracB * exp( -bNow*bNow / radius2B) / radius2B
        + fracC * exp( -bNow*bNow / radius2C) / radius2C );
      else overlapNow = normPi * exp( -pow( bNow, expPow));
      probAccept = 1. - exp( -min(EXPMAX, M_PI * kNow * overlapNow));

    // Treatment in high-b region: pick b according to overlap.
    } else {
      isAtLowB = false;

      // For simple and double Gaussian pick b according to exp(-b^2 / r^2).
      if (bProfile == 1) {
        bNow = sqrt(bDiv*bDiv - log(rndmPtr->flat()));
        overlapNow = normPi * exp( -min(EXPMAX, bNow*bNow));
      } else if (bProfile == 2) {
        double pickFrac = rndmPtr->flat() * fracABChigh;
        if (pickFrac < fracAhigh)
          bNow = sqrt(bDiv*bDiv - log(rndmPtr->flat()));
        else if (pickFrac < fracAhigh + fracBhigh)
          bNow = sqrt(bDiv*bDiv - radius2B * log(rndmPtr->flat()));
        else bNow = sqrt(bDiv*bDiv - radius2C * log(rndmPtr->flat()));
        overlapNow = normPi * ( fracA * exp( -min(EXPMAX, bNow*bNow))
          + fracB * exp( -min(EXPMAX, bNow*bNow / radius2B)) / radius2B
          + fracC * exp( -min(EXPMAX, bNow*bNow / radius2C)) / radius2C );

      // For exp( - b^expPow) transform to variable c = b^expPow so that
      // f(b) = b * exp( - b^expPow) -> f(c) = c^r * exp(-c) with r = expRev.
      // case hasLowPow: expPow < 2 <=> r > 0: preselect according to
      // f(c) < N exp(-c/2) and then accept with N' * c^r * exp(-c/2).
      } else if (hasLowPow) {
        double cNow, acceptC;
        do {
          cNow = cDiv - 2. * log(rndmPtr->flat());
          acceptC = pow(cNow / cMax, expRev) * exp( -0.5 * (cNow - cMax));
        } while (acceptC < rndmPtr->flat());
        bNow = pow( cNow, 1. / expPow);
        overlapNow = normPi * exp( -cNow);

      // case !hasLowPow: expPow >= 2 <=> - 1 < r < 0: preselect according to
      // f(c) < N exp(-c) and then accept with N' * c^r.
      } else {
        double cNow, acceptC;
        do {
          cNow = cDiv - log(rndmPtr->flat());
          acceptC = pow(cNow / cDiv, expRev);
        } while (acceptC < rndmPtr->flat());
        bNow = pow( cNow, 1. / expPow);
        overlapNow = normPi * exp( -cNow);
      }
      double temp = M_PI * kNow *overlapNow;
      probAccept = (1. - exp( -min(EXPMAX, temp))) / temp;
    }

  // Confirm choice of b value. Derive enhancement factor.
  } while (probAccept < rndmPtr->flat());

  // Same enhancement for hardest process and all subsequent MPI
  enhanceB = enhanceBmax = enhanceBnow = (normOverlap / normPi) * overlapNow;

  // Done.
  bIsSet = true;

}

//--------------------------------------------------------------------------

// Pick impact parameter and interaction rate enhancement afterwards,
// i.e. after a hard interaction is known but before rest of MPI treatment.

void MultipartonInteractions::overlapNext(Event& event, double pTscale) {

  // Default, valid for bProfile = 0. Also initial Sudakov.
  enhanceB = zeroIntCorr;
  if (bProfile <= 0 || bProfile > 4) return;

  // Alternative choices of event scale for Sudakov in (pT, b) space.
  if (bSelScale == 1) {
    vector<double> mmT;
    for (int i = 5; i < event.size(); ++i) if (event[i].isFinal()) {
      mmT.push_back( event[i].m() + event[i].mT() );
      for (int j = int(mmT.size()) - 1; j > 0; --j)
        if (mmT[j] > mmT[j - 1]) swap( mmT[j], mmT[j - 1] );
    }
    pTscale = 0.5 * mmT[0];
    for (int j = 1; j < int(mmT.size()); ++j) pTscale += mmT[j] / (j + 1.);
  } else if (bSelScale == 2) pTscale = event.scale();
  double pT2scale = pTscale*pTscale;

  // Use trial interaction for x-dependent matter profile.
  if (bProfile == 4) {
    double pTtrial = 0.;
    do {
      // Pick b according to wanted O(b, x1, x2).
      double expb2 = rndmPtr->flat();
      double w1    = XDEP_A1 + a1 * log(1. / infoPtr->x1());
      double w2    = XDEP_A1 + a1 * log(1. / infoPtr->x2());
      double fac   = a02now * (w1 * w1 + w2 * w2);
      b2now  = - fac * log(expb2);
      bNow   = sqrt(b2now);

      // Enhancement factor for the hard process and overestimate
      // for fastPT2. Note that existing framework has a (1. / sigmaND)
      // present.
      enhanceB    = sigmaND / M_PI / fac * expb2;
      enhanceBmax = sigmaND / 2. / M_PI / a02now
                  * exp( -b2now / 2. / a2max );

      // Trial interaction. Keep going until pTtrial < pTscale.
      pTtrial = pTnext(pTmax, pTmin, event);
    } while (pTtrial > pTscale);
    bIsSet = true;
    return;
  }

  // Begin loop over pT-dependent rejection of b value.
  do {

    // Flat enhancement distribution for simple Gaussian.
    if (bProfile == 1) {
      double expb2 = rndmPtr->flat();
      // Same enhancement for hardest process and all subsequent MPI.
      enhanceB = enhanceBmax = enhanceBnow = normOverlap * expb2;
      bNow = sqrt( -log(expb2));

    // For double Gaussian go the way via b, according to each Gaussian.
    } else if (bProfile == 2) {
      double bType = rndmPtr->flat();
      double b2 = -log( rndmPtr->flat() );
      if (bType < fracA) ;
      else if (bType < fracA + fracB) b2 *= radius2B;
      else b2 *= radius2C;
      // Same enhancement for hardest process and all subsequent MPI.
      enhanceB = enhanceBmax = enhanceBnow = normOverlap *
        ( fracA * exp( -min(EXPMAX, b2))
        + fracB * exp( -min(EXPMAX, b2 / radius2B)) / radius2B
        + fracC * exp( -min(EXPMAX, b2 / radius2C)) / radius2C );
      bNow = sqrt(b2);

    // For exp( - b^expPow) transform to variable c = b^expPow so that
    // f(b) = b * exp( - b^expPow) -> f(c) = c^r * exp(-c) with r = expRev.
    // case hasLowPow: expPow < 2 <=> r > 0:
    // f(c) < r^r exp(-r) for c < 2r, < (2r)^r exp(-r-c/2) for c > 2r.
    } else if (bProfile == 3 && hasLowPow) {
      double cNow, acceptC;
      double probLowC = expRev / (expRev + pow(2., expRev) * exp( - expRev));
      do {
        if (rndmPtr->flat() < probLowC) {
          cNow = 2. * expRev * rndmPtr->flat();
          acceptC = pow( cNow / expRev, expRev) * exp(expRev - cNow);
        } else {
          cNow = 2. * (expRev - log( rndmPtr->flat() ));
          acceptC = pow(0.5 * cNow / expRev, expRev)
                  * exp(expRev - 0.5 * cNow);
        }
      } while (acceptC < rndmPtr->flat());
      // Same enhancement for hardest process and all subsequent MPI.
      enhanceB = enhanceBmax = enhanceBnow = normOverlap *exp(-cNow);
      bNow = pow( cNow, 1. / expPow);

    // case !hasLowPow: expPow >= 2 <=> - 1 < r < 0:
    // f(c) < c^r for c < 1,  f(c) < exp(-c) for c > 1.
    } else if (bProfile == 3 && !hasLowPow) {
      double cNow, acceptC;
      double probLowC = expPow / (2. * exp(-1.) + expPow);
      do {
        if (rndmPtr->flat() < probLowC) {
          cNow = pow( rndmPtr->flat(), 0.5 * expPow);
          acceptC = exp(-cNow);
        } else {
          cNow = 1. - log( rndmPtr->flat() );
          acceptC = pow( cNow, expRev);
        }
      } while (acceptC < rndmPtr->flat());
      // Same enhancement for hardest process and all subsequent MPI.
      enhanceB = enhanceBmax = enhanceBnow = normOverlap * exp(-cNow);
      bNow = pow( cNow, 1. / expPow);
    }

  // Evaluate "Sudakov form factor" for not having a harder interaction.
  } while (sudakov(pT2scale, enhanceB) < rndmPtr->flat());

  // Done.
  bIsSet = true;
}

//--------------------------------------------------------------------------

// Printe statistics on number of multiparton-interactions processes.

void MultipartonInteractions::statistics(bool resetStat, ostream& os) {

  // Header.
  os << "\n *-------  PYTHIA Multiparton Interactions Statistics  -----"
     << "---*\n"
     << " |                                                            "
     << " |\n"
     << " |  Note: excludes hardest subprocess if already listed above "
     << " |\n"
     << " |                                                            "
     << " |\n"
     << " | Subprocess                               Code |       Times"
     << " |\n"
     << " |                                               |            "
     << " |\n"
     << " |------------------------------------------------------------"
     << "-|\n"
     << " |                                               |            "
     << " |\n";

  // Loop over existing processes. Sum of all subprocesses.
  int numberSum = 0;
  for ( map<int, int>::iterator iter = nGen.begin(); iter != nGen.end();
    ++iter) {
    int code = iter -> first;
    int number = iter->second;
    numberSum += number;

    // Find process name that matches code.
    string name = " ";
    bool foundName = false;
    SigmaMultiparton* dSigma;
    for (int i = 0; i < 4; ++i) {
      if      (i == 0) dSigma = &sigma2gg;
      else if (i == 1) dSigma = &sigma2qg;
      else if (i == 2) dSigma = &sigma2qqbarSame;
      else             dSigma = &sigma2qq;
      int nProc = dSigma->nProc();
      for (int iProc = 0; iProc < nProc; ++iProc)
      if (dSigma->codeProc(iProc) == code) {
        name = dSigma->nameProc(iProc);
        foundName = true;
      }
      if (foundName) break;
    }

    // Print individual process info.
    os << " | " << left << setw(40) << name << right << setw(5) << code
       << " | " << setw(11) << number << " |\n";
  }

  // Print summed process info.
  os << " |                                                            "
     << " |\n"
     << " | " << left << setw(45) << "sum" << right << " | " << setw(11)
       << numberSum  << " |\n";

    // Listing finished.
  os << " |                                               |            "
     << " |\n"
     << " *-------  End PYTHIA Multiparton Interactions Statistics ----"
     << "-*" << endl;

  // Optionally reset statistics contents.
  if (resetStat) resetStatistics();

}

//==========================================================================

} // end namespace Pythia8
