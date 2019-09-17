// SimpleTimeShower.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// SimpleTimeShower class.

#include "Pythia8/SimpleTimeShower.h"

namespace Pythia8 {

//==========================================================================

// The SimpleTimeShower class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Minimal allowed c and b quark masses, for flavour thresholds.
const double SimpleTimeShower::MCMIN        = 1.2;
const double SimpleTimeShower::MBMIN        = 4.0;

// For small x approximate 1 - sqrt(1 - x) by x/2.
const double SimpleTimeShower::SIMPLIFYROOT = 1e-8;

// Do not allow x too close to 0 or 1 in matrix element expressions.
// Warning: cuts into phase space for E_CM > 2 * pTmin * sqrt(1/XMARGIN),
// i.e. will become problem roughly for E_CM > 10^6 GeV.
const double SimpleTimeShower::XMARGIN      = 1e-12;
const double SimpleTimeShower::XMARGINCOMB  = 1e-4;

// Lower limit on PDF value in order to avoid division by zero.
const double SimpleTimeShower::TINYPDF      = 1e-10;

// Big starting value in search for smallest invariant-mass pair.
const double SimpleTimeShower::LARGEM2      = 1e20;

// In g -> q qbar or gamma -> f fbar require m2_pair > this * m2_q/f.
const double SimpleTimeShower::THRESHM2      = 4.004;

// Never pick pT so low that alphaS is evaluated too close to Lambda_3.
const double SimpleTimeShower::LAMBDA3MARGIN = 1.1;

// Rescatter: rescattering + ISR + FSR + primordial kT can lead to
//            systems not locally conserving momentum.
// Fix up momentum in intermediate systems with rescattering
const bool   SimpleTimeShower::FIXRESCATTER          = true;
// Veto negative energies when using FIXRESCATTER option.
const bool   SimpleTimeShower::VETONEGENERGY         = false;
// Do not allow too large time- or spacelike virtualities in fixing-up.
const double SimpleTimeShower::MAXVIRTUALITYFRACTION = 0.5;
// Do not allow too large negative spacelike energy in system rest frame.
const double SimpleTimeShower::MAXNEGENERGYFRACTION  = 0.7;

// Fudge extra weight for overestimation of weak shower t-channel correction.
const double SimpleTimeShower::WEAKPSWEIGHT = 5.;

// Extra overestimate of g -> q qbar branching rate for DGLAP comparison.
const double SimpleTimeShower::WG2QEXTRA = 20.;

// Limit on size of number of rejections for uncertainty variations.
const double SimpleTimeShower::REJECTFACTOR = 0.1;

// Limit on probability for uncertainty variations.
const double SimpleTimeShower::PROBLIMIT = 0.99;

//--------------------------------------------------------------------------

// Initialize alphaStrong, alphaEM and related pTmin parameters.

void SimpleTimeShower::init( BeamParticle* beamAPtrIn,
  BeamParticle* beamBPtrIn) {

  // Store input pointers for future use.
  beamAPtr           = beamAPtrIn;
  beamBPtr           = beamBPtrIn;

  // Main flags.
  doQCDshower        = settingsPtr->flag("TimeShower:QCDshower");
  doQEDshowerByQ     = settingsPtr->flag("TimeShower:QEDshowerByQ");
  doQEDshowerByL     = settingsPtr->flag("TimeShower:QEDshowerByL");
  doQEDshowerByOther = settingsPtr->flag("TimeShower:QEDshowerByOther");
  doQEDshowerByGamma = settingsPtr->flag("TimeShower:QEDshowerByGamma");
  doWeakShower       = settingsPtr->flag("TimeShower:weakShower");
  doMEcorrections    = settingsPtr->flag("TimeShower:MEcorrections");
  doMEextended       = settingsPtr->flag("TimeShower:MEextended");
  if (!doMEcorrections) doMEextended = false;
  doMEafterFirst     = settingsPtr->flag("TimeShower:MEafterFirst");
  doPhiPolAsym       = settingsPtr->flag("TimeShower:phiPolAsym");
  doPhiPolAsymHard   = settingsPtr->flag("TimeShower:phiPolAsymHard");
  doInterleave       = settingsPtr->flag("TimeShower:interleave");
  allowBeamRecoil    = settingsPtr->flag("TimeShower:allowBeamRecoil");
  dampenBeamRecoil   = settingsPtr->flag("TimeShower:dampenBeamRecoil");
  recoilToColoured   = settingsPtr->flag("TimeShower:recoilToColoured");
  allowMPIdipole     = settingsPtr->flag("TimeShower:allowMPIdipole");

  // If SimpleSpaceShower does dipole recoil then SimpleTimeShower must adjust.
  doDipoleRecoil     = settingsPtr->flag("SpaceShower:dipoleRecoil");
  if (doDipoleRecoil) allowBeamRecoil  = true;
  if (doDipoleRecoil) dampenBeamRecoil = false;

  // Matching in pT of hard interaction or MPI to shower evolution.
  pTmaxMatch         = settingsPtr->mode("TimeShower:pTmaxMatch");
  pTdampMatch        = settingsPtr->mode("TimeShower:pTdampMatch");
  pTmaxFudge         = settingsPtr->parm("TimeShower:pTmaxFudge");
  pTmaxFudgeMPI      = settingsPtr->parm("TimeShower:pTmaxFudgeMPI");
  pTdampFudge        = settingsPtr->parm("TimeShower:pTdampFudge");

  // Charm and bottom mass thresholds.
  mc                 = max( MCMIN, particleDataPtr->m0(4));
  mb                 = max( MBMIN, particleDataPtr->m0(5));
  m2c                = mc * mc;
  m2b                = mb * mb;

  // Parameters of scale choices.
  renormMultFac     = settingsPtr->parm("TimeShower:renormMultFac");
  factorMultFac     = settingsPtr->parm("TimeShower:factorMultFac");
  useFixedFacScale  = settingsPtr->flag("TimeShower:useFixedFacScale");
  fixedFacScale2    = pow2(settingsPtr->parm("TimeShower:fixedFacScale"));

  // Parameters of alphaStrong generation.
  alphaSvalue        = settingsPtr->parm("TimeShower:alphaSvalue");
  alphaSorder        = settingsPtr->mode("TimeShower:alphaSorder");
  alphaSnfmax        = settingsPtr->mode("StandardModel:alphaSnfmax");
  alphaSuseCMW       = settingsPtr->flag("TimeShower:alphaSuseCMW");
  alphaS2pi          = 0.5 * alphaSvalue / M_PI;

  // Initialize alphaStrong generation.
  alphaS.init( alphaSvalue, alphaSorder, alphaSnfmax, alphaSuseCMW);

  // Lambda for 5, 4 and 3 flavours.
  Lambda3flav        = alphaS.Lambda3();
  Lambda4flav        = alphaS.Lambda4();
  Lambda5flav        = alphaS.Lambda5();
  Lambda5flav2       = pow2(Lambda5flav);
  Lambda4flav2       = pow2(Lambda4flav);
  Lambda3flav2       = pow2(Lambda3flav);

  // Parameters of QCD evolution. Warn if pTmin must be raised.
  nGluonToQuark      = settingsPtr->mode("TimeShower:nGluonToQuark");
  weightGluonToQuark = settingsPtr->mode("TimeShower:weightGluonToQuark");
  scaleGluonToQuark  = settingsPtr->parm("TimeShower:scaleGluonToQuark");
  extraGluonToQuark  = (weightGluonToQuark%4 == 3) ? WG2QEXTRA : 1.;
  recoilDeadCone     = settingsPtr->flag("TimeShower:recoilDeadCone");
  pTcolCutMin        = settingsPtr->parm("TimeShower:pTmin");
  if (pTcolCutMin > LAMBDA3MARGIN * Lambda3flav / sqrt(renormMultFac))
    pTcolCut         = pTcolCutMin;
  else {
    pTcolCut         = LAMBDA3MARGIN * Lambda3flav / sqrt(renormMultFac);
    ostringstream newPTcolCut;
    newPTcolCut << fixed << setprecision(3) << pTcolCut;
    infoPtr->errorMsg("Warning in TimeShower::init: pTmin too low",
                      ", raised to " + newPTcolCut.str() );
    infoPtr->setTooLowPTmin(true);
  }
  pT2colCut          = pow2(pTcolCut);

  // Parameters of alphaEM generation.
  alphaEMorder       = settingsPtr->mode("TimeShower:alphaEMorder");

  // Initialize alphaEM generation.
  alphaEM.init( alphaEMorder, settingsPtr);

  // Parameters of QED evolution.
  nGammaToQuark      = settingsPtr->mode("TimeShower:nGammaToQuark");
  nGammaToLepton     = settingsPtr->mode("TimeShower:nGammaToLepton");
  pTchgQCut          = settingsPtr->parm("TimeShower:pTminChgQ");
  pT2chgQCut         = pow2(pTchgQCut);
  pTchgLCut          = settingsPtr->parm("TimeShower:pTminChgL");
  pT2chgLCut         = pow2(pTchgLCut);
  mMaxGamma          = settingsPtr->parm("TimeShower:mMaxGamma");
  m2MaxGamma         = pow2(mMaxGamma);

  // Parameters of weak evolution.
  weakMode           = settingsPtr->mode("TimeShower:weakShowerMode");
  pTweakCut          = settingsPtr->parm("TimeShower:pTminWeak");
  pT2weakCut         = pow2(pTweakCut);
  weakEnhancement    = settingsPtr->parm("WeakShower:enhancement");
  singleWeakEmission = settingsPtr->flag("WeakShower:singleEmission");
  vetoWeakJets       = settingsPtr->flag("WeakShower:vetoWeakJets");
  vetoWeakDeltaR2    = pow2(settingsPtr->parm("WeakShower:vetoWeakDeltaR"));
  weakExternal       = settingsPtr->flag("WeakShower:externalSetup");

  // Consisteny check for gamma -> f fbar variables.
  if (nGammaToQuark <= 0 && nGammaToLepton <= 0) doQEDshowerByGamma = false;

  // Possibility of a global recoil stategy, e.g. for MC@NLO.
  globalRecoil       = settingsPtr->flag("TimeShower:globalRecoil");
  nMaxGlobalRecoil   = settingsPtr->mode("TimeShower:nMaxGlobalRecoil");
  // Number of splittings produced with global recoil.
  nMaxGlobalBranch   = settingsPtr->mode("TimeShower:nMaxGlobalBranch");
  // Number of partons in Born-like events, to distinguish between S and H.
  nFinalBorn         = settingsPtr->mode("TimeShower:nPartonsInBorn");
  // Flag to allow to start from a scale smaller than scalup.
  globalRecoilMode   = settingsPtr->mode("TimeShower:globalRecoilMode");
  // Flag to allow to start from a scale smaller than scalup.
  limitMUQ           = settingsPtr->flag("TimeShower:limitPTmaxGlobal");

  // Fraction and colour factor of gluon emission off onium octat state.
  octetOniumFraction = settingsPtr->parm("TimeShower:octetOniumFraction");
  octetOniumColFac   = settingsPtr->parm("TimeShower:octetOniumColFac");

  // Z0 and W+- properties needed for gamma/Z0 mixing and weak showers.
  mZ                 = particleDataPtr->m0(23);
  gammaZ             = particleDataPtr->mWidth(23);
  thetaWRat          = 1. / (16. * coupSMPtr->sin2thetaW()
                       * coupSMPtr->cos2thetaW());
  mW                 = particleDataPtr->m0(24);
  gammaW             = particleDataPtr->mWidth(24);

  // May have to fix up recoils related to rescattering.
  allowRescatter     = settingsPtr->flag("PartonLevel:MPI")
    && settingsPtr->flag("MultipartonInteractions:allowRescatter");

  // Hidden Valley scenario with further shower activity.
  doHVshower         = settingsPtr->flag("HiddenValley:FSR");
  nCHV               = settingsPtr->mode("HiddenValley:Ngauge");
  alphaHVfix         = settingsPtr->parm("HiddenValley:alphaFSR");
  alphaHVorder       = (nCHV > 1 )
                     ? settingsPtr->mode("HiddenValley:alphaOrder") : 0;
  nFlavHV            = settingsPtr->mode("HiddenValley:nFlav");
  LambdaHV           = settingsPtr->parm("HiddenValley:Lambda");
  pThvCut            = settingsPtr->parm("HiddenValley:pTminFSR");
  CFHV               = (nCHV == 1) ? 1. : (nCHV * nCHV - 1.)/(2. * nCHV);
  idHV               = (nCHV == 1) ? 4900022 : 4900021;
  mHV                = particleDataPtr->m0(idHV);
  brokenHVsym        = (nCHV == 1 && mHV > 0.);
  if (pThvCut < LambdaHV) {
    pThvCut         = LAMBDA3MARGIN * LambdaHV;
    ostringstream newPTcolCut;
    newPTcolCut << fixed << setprecision(3) << pThvCut;
    infoPtr->errorMsg("Warning in SimpleTimeShower::init: Hidden Valley ",
                      "pTmin too low, raised to " + newPTcolCut.str() );
  }
  pT2hvCut           = pThvCut * pThvCut;

  // Possibility of two predetermined hard emissions in event.
  doSecondHard      = settingsPtr->flag("SecondHard:generate");
  twoHard            = doSecondHard;

  // Possibility to allow user veto of emission step.
  hasUserHooks       = (userHooksPtr != 0);
  canVetoEmission    = hasUserHooks && userHooksPtr->canVetoFSREmission();

  // Set initial value, just in case.
  dopTdamp           = false;
  pT2damp            = 0.;

  // Default values for the weak shower.
  hasWeaklyRadiated  = false;

  // Disallow simultaneous splitting and trial emission enhancements.
  canEnhanceEmission = hasUserHooks && userHooksPtr->canEnhanceEmission();
  canEnhanceTrial    = hasUserHooks && userHooksPtr->canEnhanceTrial();
  if (canEnhanceEmission && canEnhanceTrial) {
    infoPtr->errorMsg("Error in SimpleTimeShower::init: Enhance for both "
    "actual and trial emissions not possible. Both switched off.");
    canEnhanceEmission = false;
    canEnhanceTrial    = false;
  }

  // Initialize variables set in pTnext but not in showerQED.
  doTrialNow = false;
  canEnhanceET = false;

  // Properties for enhanced emissions.
  splittingNameSel   = "";
  splittingNameNow   = "";
  enhanceFactors.clear();

  // Enable automated uncertainty variations.
  nVarQCD            = 0;
  doUncertainties    = settingsPtr->flag("UncertaintyBands:doVariations")
                    && initUncertainties();
  doUncertaintiesNow = doUncertainties;
  uVarNflavQ         = settingsPtr->mode("UncertaintyBands:nFlavQ");
  uVarMPIshowers     = settingsPtr->flag("UncertaintyBands:MPIshowers");
  cNSpTmin           = settingsPtr->parm("UncertaintyBands:cNSpTmin");
  uVarpTmin2         = pT2colCut;
  uVarpTmin2        *= settingsPtr->parm("UncertaintyBands:FSRpTmin2Fac");
  int varType        = settingsPtr->mode("UncertaintyBands:type");
  noResVariations    = (varType == 1) ? true: false;
  noProcVariations   = (varType == 2) ? true: false;
  overFactor         = settingsPtr->parm("UncertaintyBands:overSampleFSR");

  // Possibility to set parton vertex information.
  doPartonVertex     = settingsPtr->flag("PartonVertex:setVertex")
                    && (partonVertexPtr != 0);

}

//--------------------------------------------------------------------------

// Find whether to limit maximum scale of emissions.
// Also allow for dampening at factorization or renormalization scale.

bool SimpleTimeShower::limitPTmax( Event& event, double Q2Fac, double Q2Ren) {

  // Find whether to limit pT. Begin by user-set cases.
  twoHard = doSecondHard;
  bool dopTlimit = false;
  dopTlimit1 = dopTlimit2 = false;
  int nHeavyCol = 0;
  if      (pTmaxMatch == 1) dopTlimit = dopTlimit1 = dopTlimit2 = true;
  else if (pTmaxMatch == 2) dopTlimit = dopTlimit1 = dopTlimit2 = false;

  // Always restrict SoftQCD processes.
  else if (infoPtr->isNonDiffractive() || infoPtr->isDiffractiveA()
    || infoPtr->isDiffractiveB() || infoPtr->isDiffractiveC() )
    dopTlimit = dopTlimit1 = dopTlimit2 = true;

  // Look if any quark (u, d, s, c, b), gluon or photon in final state.
  // Also count number of heavy coloured particles, like top.
  else {
    int n21 = 0;
    int iBegin = 5 + beamOffset;
    for (int i = iBegin; i < event.size(); ++i) {
      if (event[i].status() == -21) ++n21;
      else if (n21 == 0) {
        int idAbs = event[i].idAbs();
        if (idAbs <= 5 || idAbs == 21 || idAbs == 22) dopTlimit1 = true;
        if ( (event[i].col() != 0 || event[i].acol() != 0)
          && idAbs > 5 && idAbs != 21 ) ++nHeavyCol;
      } else if (n21 == 2) {
        int idAbs = event[i].idAbs();
        if (idAbs <= 5 || idAbs == 21 || idAbs == 22) dopTlimit2 = true;
      }
    }
    twoHard = (n21 == 2);
    dopTlimit = (twoHard) ? (dopTlimit1 && dopTlimit2) : dopTlimit1;
  }

  // Dampening at factorization or renormalization scale; only for hardest.
  dopTdamp   = false;
  pT2damp    = 0.;
  if (!dopTlimit1 && (pTdampMatch == 1 || pTdampMatch == 2)) {
    dopTdamp = true;
    pT2damp  = pow2(pTdampFudge) * ((pTdampMatch == 1) ? Q2Fac : Q2Ren);
  }
  if (!dopTlimit1 && nHeavyCol > 1 && (pTdampMatch == 3 || pTdampMatch == 4)) {
    dopTdamp = true;
    pT2damp  = pow2(pTdampFudge) * ((pTdampMatch == 3) ? Q2Fac : Q2Ren);
  }

  // Done.
  return dopTlimit;

}

//--------------------------------------------------------------------------

// Top-level routine to do a full time-like shower in resonance decay.

int SimpleTimeShower::shower( int iBeg, int iEnd, Event& event, double pTmax,
  int nBranchMax) {

  // Add new system, automatically with two empty beam slots.
  int iSys = partonSystemsPtr->addSys();

  // Loop over allowed range to find all final-state particles.
  // Check if they all have same single mother => resonance decay.
  Vec4 pSum;
  bool isResDec = true;
  int  iRes     = -1;
  for (int i = iBeg; i <= iEnd; ++i) {
    if (event[i].isFinal()) {
      partonSystemsPtr->addOut( iSys, i);
      pSum += event[i].p();
      // Look for common resonance-decay mother.
      if ( event[i].mother2() != 0 && event[i].mother1() != event[i].mother2())
        isResDec = false;
      else if ( iRes != -1 && event[i].mother1() != iRes)
        isResDec = false;
      else
        iRes = event[i].mother1();
    }
  }
  partonSystemsPtr->setSHat( iSys, pSum.m2Calc() );
  if (isResDec) partonSystemsPtr->setInRes( iSys, iRes);

  // Let prepare routine do the setup.
  dopTlimit1        = true;
  dopTlimit2        = true;
  dopTdamp          = false;
  hasWeaklyRadiated = false;
  prepare( iSys, event, true);

  // Begin evolution down in pT from hard pT scale.
  int nBranch  = 0;
  pTLastBranch = 0.;
  do {
    double pTtimes = pTnext( event, pTmax, 0.);

    // Do a final-state emission (if allowed).
    if (pTtimes > 0.) {
      if (branch( event)) {
        ++nBranch;
        pTLastBranch = pTtimes;
      }
      pTmax = pTtimes;
    }

    // Keep on evolving until nothing is left to be done.
    else pTmax = 0.;
  } while (pTmax > 0. && (nBranchMax <= 0 || nBranch < nBranchMax));

  // Return number of emissions that were performed.
  return nBranch;

}

//--------------------------------------------------------------------------

// Top-level routine for QED radiation in hadronic decay to two leptons.
// Intentionally only does photon radiation, i.e. no photon branchings.

int SimpleTimeShower::showerQED( int i1, int i2, Event& event, double pTmax) {

  // Add new system, automatically with two empty beam slots.
  int iSys = partonSystemsPtr->addSys();
  partonSystemsPtr->addOut( iSys, i1);
  partonSystemsPtr->addOut( iSys, i2);
  partonSystemsPtr->setSHat( iSys, m2(event[i1], event[i2]) );
  // Add incoming (decaying) particle; assumed = event[i1].mother1()
  partonSystemsPtr->setInRes( iSys, event[i1].mother1() );

  // Charge type of two leptons tells whether MEtype is gamma*/Z0 or W+-.
  int iChg1  = event[i1].chargeType();
  int iChg2  = event[i2].chargeType();
  int MEtype = (iChg1 + iChg2 == 0) ? 102 : 101;

  // Fill dipole-ends list.
  dipEnd.resize(0);
  if (iChg1 != 0) dipEnd.push_back( TimeDipoleEnd(i1, i2, pTmax,
      0, iChg1, 0, 0, 0, iSys, MEtype, i2) );
  if (iChg2 != 0) dipEnd.push_back( TimeDipoleEnd(i2, i1, pTmax,
      0, iChg2, 0, 0, 0, iSys, MEtype, i1) );

  // Begin evolution down in pT from hard pT scale.
  int nBranch  = 0;
  pTLastBranch = 0.;
  do {

    // Begin loop over all possible radiating dipole ends.
    dipSel  = 0;
    iDipSel = -1;
    double pT2sel = 0.;
    for (int iDip = 0; iDip < int(dipEnd.size()); ++iDip) {
      TimeDipoleEnd& dip = dipEnd[iDip];

      // Dipole properties.
      dip.mRad  = event[dip.iRadiator].m();
      dip.mRec  = event[dip.iRecoiler].m();
      dip.mDip  = m( event[dip.iRadiator], event[dip.iRecoiler] );
      dip.m2Rad = pow2(dip.mRad);
      dip.m2Rec = pow2(dip.mRec);
      dip.m2Dip = pow2(dip.mDip);

      // Find maximum evolution scale for dipole.
      dip.m2DipCorr    = pow2(dip.mDip - dip.mRec) - dip.m2Rad;
      double pTbegDip  = min( pTmax, dip.pTmax );
      double pT2begDip = min( pow2(pTbegDip), 0.25 * dip.m2DipCorr);

      // Do QED evolution where relevant.
      dip.pT2 = 0.;
      if (pT2begDip > pT2sel) {
        pT2nextQED( pT2begDip, pT2sel, dip, event);

        // Update if found larger pT than current maximum. End dipole loop.
        if (dip.pT2 > pT2sel) {
          pT2sel  = dip.pT2;
          dipSel  = &dip;
          iDipSel = iDip;
        }
      }
    }
    double pTsel = (dipSel == 0) ? 0. : sqrt(pT2sel);

    // Do a final-state emission (if allowed).
    if (pTsel > 0.) {

      // Find initial radiator and recoiler particles in dipole branching.
      int iRadBef      = dipSel->iRadiator;
      int iRecBef      = dipSel->iRecoiler;
      Particle& radBef = event[iRadBef];
      Particle& recBef = event[iRecBef];
      Vec4 pRadBef     = event[iRadBef].p();
      Vec4 pRecBef     = event[iRecBef].p();

      // Construct kinematics in dipole rest frame; massless emitter.
      double pTorig       = sqrt( dipSel->pT2);
      double eRadPlusEmt  = 0.5 * (dipSel->m2Dip + dipSel->m2 - dipSel->m2Rec)
        / dipSel->mDip;
      double e2RadPlusEmt = pow2(eRadPlusEmt);
      double pzRadPlusEmt = 0.5 * sqrtpos( pow2(dipSel->m2Dip - dipSel->m2
        - dipSel->m2Rec) - 4. * dipSel->m2 * dipSel->m2Rec ) / dipSel->mDip;
      double pT2corr = dipSel->m2 * (e2RadPlusEmt * dipSel->z
        * (1. - dipSel->z) - 0.25 * dipSel->m2) / pow2(pzRadPlusEmt);
      double pTcorr       = sqrtpos( pT2corr );
      double pzRad        = (e2RadPlusEmt * dipSel->z - 0.5 * dipSel->m2)
        / pzRadPlusEmt;
      double pzEmt        = (e2RadPlusEmt * (1. - dipSel->z)
        - 0.5 * dipSel->m2) / pzRadPlusEmt;
      double mRad         = dipSel->mRad;
      double mEmt         = 0.;

      // Kinematics reduction for radiator mass.
      double m2Ratio    = dipSel->m2Rad / dipSel->m2;
      pTorig           *= 1. - m2Ratio;
      pTcorr           *= 1. - m2Ratio;
      pzRad            += pzEmt * m2Ratio;
      pzEmt            *= 1. - m2Ratio;

      // Store kinematics of branching in dipole rest frame.
      double phi = 2. * M_PI * rndmPtr->flat();
      Vec4 pRad = Vec4( pTcorr * cos(phi), pTcorr * sin(phi), pzRad,
        sqrt( pow2(pTcorr) + pow2(pzRad) + pow2(mRad) ) );
      Vec4 pEmt = Vec4( -pRad.px(), -pRad.py(), pzEmt,
        sqrt( pow2(pTcorr) + pow2(pzEmt) + pow2(mEmt) ) );
      Vec4 pRec = Vec4( 0., 0., -pzRadPlusEmt,
        sqrt( pow2(pzRadPlusEmt) + dipSel->m2Rec ) );

      // Rotate and boost dipole products to the event frame.
      RotBstMatrix M;
      M.fromCMframe(pRadBef, pRecBef);
      pRad.rotbst(M);
      pEmt.rotbst(M);
      pRec.rotbst(M);

      // Define new particles from dipole branching.
      Particle rad = Particle(radBef.id(), 51, iRadBef, 0, 0, 0,
        radBef.col(), radBef.acol(), pRad, mRad, pTsel);
      Particle emt = Particle(22, 51, iRadBef, 0, 0, 0,
        0, 0, pEmt, mEmt, pTsel);
      Particle rec = Particle(recBef.id(),  52, iRecBef, iRecBef, 0, 0,
        recBef.col(), recBef.acol(), pRec, dipSel->mRec, pTsel);

      // ME corrections can lead to branching being rejected.
      if (dipSel->MEtype == 0
        || findMEcorr( dipSel, rad, rec, emt, false) > rndmPtr->flat() ) {

        // Shower may occur at a displaced vertex, or for unstable particle.
        if (radBef.hasVertex()) {
          rad.vProd( radBef.vProd() );
          emt.vProd( radBef.vProd() );
        }
        if (recBef.hasVertex()) rec.vProd( recBef.vProd() );
        rad.tau( event[iRadBef].tau() );
        rec.tau( event[iRecBef].tau() );

        // Put new particles into the event record.
        int iRad = event.append(rad);
        int iEmt = event.append(emt);
        event[iRadBef].statusNeg();
        event[iRadBef].daughters( iRad, iEmt);
        int iRec = event.append(rec);
        event[iRecBef].statusNeg();
        event[iRecBef].daughters( iRec, iRec);

        // Update to new dipole ends.
        dipSel->iRadiator = iRad;
        dipSel->iRecoiler = iRec;
        dipSel->pTmax = pTsel;

        // Update other dipoles that also involved the radiator or recoiler.
        for (int i = 0; i < int(dipEnd.size()); ++i) if (i != iDipSel) {
          if (dipEnd[i].iRadiator  == iRadBef) dipEnd[i].iRadiator  = iRad;
          if (dipEnd[i].iRecoiler  == iRadBef) dipEnd[i].iRecoiler  = iRad;
          if (dipEnd[i].iMEpartner == iRadBef) dipEnd[i].iMEpartner = iRad;
          if (dipEnd[i].iRadiator  == iRecBef) dipEnd[i].iRadiator  = iRec;
          if (dipEnd[i].iRecoiler  == iRecBef) dipEnd[i].iRecoiler  = iRec;
          if (dipEnd[i].iMEpartner == iRecBef) dipEnd[i].iMEpartner = iRec;
        }

        // Done with branching
        ++nBranch;
        pTLastBranch = pTsel;
      }
      pTmax = pTsel;
    }

    // Keep on evolving until nothing is left to be done.
    else pTmax = 0.;
  } while (pTmax > 0.);

  // Return number of emissions that were performed.
  return nBranch;

}

//--------------------------------------------------------------------------

// Global recoil: reset counters and store locations of outgoing partons.

void SimpleTimeShower::prepareGlobal( Event& event) {

  // Global recoils: reset some counters.
  nGlobal    = 0;
  nHard      = 0;
  nProposed.clear();
  hardPartons.resize(0);
  nFinalBorn = settingsPtr->mode("TimeShower:nPartonsInBorn");

  // Global recoils: store positions of hard outgoing partons.
  // No global recoil for H events.
  int nHeavyCol = 0;
  if (globalRecoil) {
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].isFinal() && event[i].colType() != 0)
        hardPartons.push_back(i);
      if ( event[i].isFinal() && event[i].idAbs() > 5 && event[i].idAbs() != 21
          && (event[i].col() != 0 || event[i].acol() != 0))
        ++nHeavyCol;
    }
    nHard = hardPartons.size();
    if (nFinalBorn > 0 && nHard > nFinalBorn) {
      hardPartons.resize(0);
      nHard = 0;
    }
  }

  // Reset nFinalBorn on an event-by-event basis.
  string nNow = infoPtr->getEventAttribute("npNLO",true);
  if (nNow != "" && nFinalBorn == -1){
    nFinalBorn = max(0, atoi((char*)nNow.c_str()));
    // Add number of heavy coloured objects in lowest multiplicity state.
    nFinalBorn += nHeavyCol;
  }

}

//--------------------------------------------------------------------------

// Prepare system for evolution; identify ME.

void SimpleTimeShower::prepare( int iSys, Event& event, bool limitPTmaxIn) {

  // Reset W/Z radiation flag at first call for new event.
  if (iSys == 0) hasWeaklyRadiated = false;

  // Reset dipole-ends list for first interaction and for resonance decays.
  int iInA = partonSystemsPtr->getInA(iSys);
  int iInB = partonSystemsPtr->getInB(iSys);
  if (iSys == 0 || iInA == 0) dipEnd.resize(0);
  int dipEndSizeBeg = dipEnd.size();

  // No dipoles for 2 -> 1 processes.
  if (partonSystemsPtr->sizeOut(iSys) < 2) return;

  // In case of DPS overwrite limitPTmaxIn by saved value.
  if (twoHard && iSys == 0) limitPTmaxIn = dopTlimit1;
  if (twoHard && iSys == 1) limitPTmaxIn = dopTlimit2;

  // Reset number of proposed splittings. Used for global recoil.
  // First check if this system belongs to the hard scattering.
  bool isHard = false;
  for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
    int ii = partonSystemsPtr->getOut( iSys, i);
    for (int iHard = 0; iHard < int(hardPartons.size()); ++iHard) {
      if ( event[ii].isAncestor(hardPartons[iHard])
        || ii == hardPartons[iHard]){
        isHard = true;
        break;
      }
    }
    if (isHard) break;
  }
  // If the system belongs to the hard scattering, initialise
  // counter of proposed emissions.
  if (isHard &&  nProposed.find(iSys) == nProposed.end() )
    nProposed.insert(make_pair(iSys,0));
  partonSystemsPtr->setHard(iSys, isHard);

  // Loop through final state of system to find possible dipole ends.
  for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
    int iRad = partonSystemsPtr->getOut( iSys, i);

    if (event[iRad].isFinal() && event[iRad].scale() > 0.) {

      // Identify colour octet onium state. Check whether QCD shower allowed.
      int idRad    = event[iRad].id();
      int idRadAbs = abs(idRad);
      bool isOctetOnium = particleDataPtr->isOctetHadron(idRad);
      bool doQCD = doQCDshower;
      if (doQCD && isOctetOnium)
        doQCD = (rndmPtr->flat() < octetOniumFraction);

      // Find dipole end formed by colour index.
      int colTag = event[iRad].col();
      if (doQCD && colTag > 0) setupQCDdip( iSys, i,  colTag,  1, event,
        isOctetOnium, limitPTmaxIn);

      // Find dipole end formed by anticolour index.
      int acolTag = event[iRad].acol();
      if (doQCD && acolTag > 0) setupQCDdip( iSys, i, acolTag, -1, event,
        isOctetOnium, limitPTmaxIn);

      // Find "charge-dipole" and "photon-dipole" ends.
      int  chgType  = event[iRad].chargeType();
      bool doChgDip = (chgType != 0)
                   && ( ( doQEDshowerByQ && event[iRad].isQuark() )
                     || ( doQEDshowerByL && event[iRad].isLepton() )
                     || ( doQEDshowerByOther && event[iRad].isResonance() ) );
      int  gamType  = (idRad == 22) ? 1 : 0;
      bool doGamDip = (gamType == 1) && doQEDshowerByGamma;
      if (doChgDip || doGamDip) setupQEDdip( iSys, i, chgType, gamType,
         event, limitPTmaxIn);

      // Find weak diple ends.
      if (doWeakShower && (iSys == 0 || !partonSystemsPtr->hasInAB(iSys))
        && (event[iRad].isQuark()  || event[iRad].isLepton())
          && (!weakExternal || iSys != 0)  ) {
        if (weakMode == 0 || weakMode == 1)
          setupWeakdip( iSys, i, 1, event, limitPTmaxIn);
        if (weakMode == 0 || weakMode == 2)
          setupWeakdip( iSys, i, 2, event, limitPTmaxIn);
      }

      // Find Hidden Valley dipole ends.
      bool isHVrad =  (idRadAbs > 4900000 && idRadAbs < 4900007)
                   || (idRadAbs > 4900010 && idRadAbs < 4900017)
                   || (idRadAbs > 4900100 && idRadAbs < 4900109);
      if (doHVshower && isHVrad) setupHVdip( iSys, i, event, limitPTmaxIn);

    // End loop over system final state. Have now found the dipole ends.
    }
  }

  // Special setup for weak dipoles if they are setup externally.
  if (doWeakShower && weakExternal && iSys == 0)
    setupWeakdipExternal(event, limitPTmaxIn);

  // Loop through dipole ends to find matrix element corrections.
  for (int iDip = dipEndSizeBeg; iDip < int(dipEnd.size()); ++iDip)
    findMEtype( event, dipEnd[iDip]);

  // Update dipole list after a multiparton interactions rescattering.
  if (iSys > 0 && ( (iInA > 0 && event[iInA].status() == -34)
    || (iInB > 0 && event[iInB].status() == -34) ) )
    rescatterUpdate( iSys, event);

}

//--------------------------------------------------------------------------

// Update dipole list after a multiparton interactions rescattering.

void SimpleTimeShower::rescatterUpdate( int iSys, Event& event) {

  // Loop over two incoming partons in system; find their rescattering mother.
  // (iOut is outgoing from old system = incoming iIn of rescattering system.)
  for (int iResc = 0; iResc < 2; ++iResc) {
    int iIn = (iResc == 0) ? partonSystemsPtr->getInA(iSys)
                           : partonSystemsPtr->getInB(iSys);
    if (iIn == 0 || event[iIn].status() != -34) continue;
    int iOut = event[iIn].mother1();

    // Loop over all dipoles.
    int dipEndSize = dipEnd.size();
    for (int iDip = 0; iDip < dipEndSize; ++iDip) {
      TimeDipoleEnd& dipNow = dipEnd[iDip];

      // Kill dipoles where rescattered parton is radiator.
      if (dipNow.iRadiator == iOut) {
        dipNow.colType = 0;
        dipNow.chgType = 0;
        dipNow.gamType = 0;
        continue;
      }
      // No matrix element for dipoles between scatterings.
      if (dipNow.iMEpartner == iOut) {
        dipNow.MEtype     =  0;
        dipNow.iMEpartner = -1;
      }

      // Update dipoles where outgoing rescattered parton is recoiler.
      if (dipNow.iRecoiler == iOut) {
        int iRad = dipNow.iRadiator;

        // Colour dipole: recoil in final state, initial state or new.
        if (dipNow.colType > 0) {
          int colTag = event[iRad].col();
          bool done  = false;
          for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
            int iRecNow = partonSystemsPtr->getOut( iSys, i);
            if (event[iRecNow].acol() == colTag) {
              dipNow.iRecoiler = iRecNow;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              done             = true;
              break;
            }
          }
          if (!done) {
            int iIn2 = (iResc == 0) ? partonSystemsPtr->getInB(iSys)
                                    : partonSystemsPtr->getInA(iSys);
            if (event[iIn2].col() == colTag) {
              dipNow.iRecoiler = iIn2;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              int isrType      = event[iIn2].mother1();
              // This line in case mother is a rescattered parton.
              while (isrType > 2 + beamOffset)
                isrType = event[isrType].mother1();
              if (isrType > 2) isrType -= beamOffset;
              dipNow.isrType   = isrType;
              done             = true;
            }
          }
          // If above options failed, then create new dipole.
          if (!done) {
            int iRadNow = partonSystemsPtr->getIndexOfOut(dipNow.system, iRad);
            if (iRadNow != -1)
              setupQCDdip(dipNow.system, iRadNow, event[iRad].col(), 1,
                          event, dipNow.isOctetOnium, true);
            else
              infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterUpdate:"
              " failed to locate radiator in system");

            dipNow.colType = 0;
            dipNow.chgType = 0;
            dipNow.gamType = 0;

            infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterUpdate: "
            "failed to locate new recoiling colour partner");
          }

        // Anticolour dipole: recoil in final state, initial state or new.
        } else if (dipNow.colType < 0) {
          int  acolTag = event[iRad].acol();
          bool done    = false;
          for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
            int iRecNow = partonSystemsPtr->getOut( iSys, i);
            if (event[iRecNow].col() == acolTag) {
              dipNow.iRecoiler = iRecNow;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              done             = true;
              break;
            }
          }
          if (!done) {
            int iIn2 = (iResc == 0) ? partonSystemsPtr->getInB(iSys)
                                    : partonSystemsPtr->getInA(iSys);
            if (event[iIn2].acol() == acolTag) {
              dipNow.iRecoiler = iIn2;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              int isrType      = event[iIn2].mother1();
              // This line in case mother is a rescattered parton.
              while (isrType > 2 + beamOffset)
                isrType = event[isrType].mother1();
              if (isrType > 2) isrType -= beamOffset;
              dipNow.isrType   = isrType;
              done             = true;
            }
          }
          // If above options failed, then create new dipole.
          if (!done) {
            int iRadNow = partonSystemsPtr->getIndexOfOut(dipNow.system, iRad);
            if (iRadNow != -1)
              setupQCDdip(dipNow.system, iRadNow, event[iRad].acol(), -1,
                          event, dipNow.isOctetOnium, true);
            else
              infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterUpdate:"
              " failed to locate radiator in system");

            dipNow.colType = 0;
            dipNow.chgType = 0;
            dipNow.gamType = 0;

            infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterUpdate:"
            " failed to locate new recoiling colour partner");
          }

        // Charge or photon dipoles: same flavour in final or initial state.
        } else if (dipNow.chgType != 0 || dipNow.gamType != 0) {
          int  idTag = event[dipNow.iRecoiler].id();
          bool done  = false;
          for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
            int iRecNow = partonSystemsPtr->getOut( iSys, i);
            if (event[iRecNow].id() == idTag) {
              dipNow.iRecoiler = iRecNow;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              done             = true;
              break;
            }
          }
          if (!done) {
            int iIn2 = (iResc == 0) ? partonSystemsPtr->getInB(iSys)
                                    : partonSystemsPtr->getInA(iSys);
            if (event[iIn2].id() == -idTag) {
              dipNow.iRecoiler = iIn2;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              int isrType      = event[iIn2].mother1();
              // This line in case mother is a rescattered parton.
              while (isrType > 2 + beamOffset)
                isrType = event[isrType].mother1();
              if (isrType > 2) isrType -= beamOffset;
              dipNow.isrType   = isrType;
              done             = true;
            }
          }
          // If above options failed, then create new dipole
          if (!done) {
            int iRadNow = partonSystemsPtr->getIndexOfOut(dipNow.system, iRad);
            if (iRadNow != -1)
              setupQEDdip(dipNow.system, iRadNow, dipNow.chgType,
                          dipNow.gamType, event, true);
            else
              infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterUpdate:"
              " failed to locate radiator in system");

            dipNow.colType = 0;
            dipNow.chgType = 0;
            dipNow.gamType = 0;

            infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterUpdate:"
            " failed to locate new recoiling charge partner");
          }
        }
      }

    // End of loop over dipoles and two incoming sides.
    }
  }

}

//--------------------------------------------------------------------------

// Update dipole list after each ISR emission (so not used for resonances).

  void SimpleTimeShower::update( int iSys, Event& event, bool hasWeakRad) {

  // Start list of rescatterers that gave further changed systems in ISR.
  vector<int> iRescatterer;

  // Find new and old positions of incoming partons in the system.
  vector<int> iNew, iOld;
  iNew.push_back( partonSystemsPtr->getInA(iSys) );
  iOld.push_back( event[iNew[0]].daughter2() );
  iNew.push_back( partonSystemsPtr->getInB(iSys) );
  iOld.push_back( event[iNew[1]].daughter2() );

  // Ditto for outgoing partons, except the newly created one.
  int sizeOut = partonSystemsPtr->sizeOut(iSys) - 1;
  for (int i = 0; i < sizeOut; ++i) {
    int iNow = partonSystemsPtr->getOut(iSys, i);
    iNew.push_back( iNow );
    iOld.push_back( event[iNow].mother1() );
    // Add non-final to list of rescatterers.
    if (!event[iNow].isFinal()) iRescatterer.push_back( iNow );
  }
  int iNewNew = partonSystemsPtr->getOut(iSys, sizeOut);

  // Swap beams to let 0 be side on which branching occured.
  if (event[iNew[0]].status() != -41) {
    swap( iNew[0], iNew[1]);
    swap( iOld[0], iOld[1]);
  }

  // Loop over all dipole ends belonging to the system
  // or to the recoil system, if different.
  for (int iDip = 0; iDip < int(dipEnd.size()); ++iDip)
  if (dipEnd[iDip].system == iSys || dipEnd[iDip].systemRec == iSys) {
    TimeDipoleEnd& dipNow = dipEnd[iDip];

    // Replace radiator (always in final state so simple).
    for (int i = 2; i < 2 + sizeOut; ++i)
    if (dipNow.iRadiator == iOld[i]) {
      dipNow.iRadiator = iNew[i];
      break;
    }

    // Replace ME partner (always in final state, if exists, so simple).
    for (int i = 2; i < 2 + sizeOut; ++i)
    if (dipNow.iMEpartner == iOld[i]) {
      dipNow.iMEpartner = iNew[i];
      break;
    }

    // Recoiler: by default pick old one, only moved. Note excluded beam.
    int iRec    = 0;
    int colRad  = event[dipNow.iRadiator].col();
    int acolRad = event[dipNow.iRadiator].acol();
    if (dipNow.systemRec == iSys) {
      for (int i = 1; i < 2 + sizeOut; ++i)
      if (dipNow.iRecoiler == iOld[i]) {
        iRec = iNew[i];
        break;
      }

      // QCD recoiler: check if colour hooks up with new final parton.
      if (dipNow.colType > 0 && event[iNewNew].acol() == colRad) {
        iRec = iNewNew;
        dipNow.isrType = 0;
      }
      if (dipNow.colType < 0 && event[iNewNew].col() == acolRad) {
        iRec = iNewNew;
        dipNow.isrType = 0;
      }

      // QCD recoiler: check if colour hooks up with new beam parton.
      if (iRec == 0 && dipNow.colType > 0 && event[iNew[0]].col() == colRad)
        iRec = iNew[0];
      if (iRec == 0 && dipNow.colType < 0 && event[iNew[0]].acol() == acolRad)
        iRec = iNew[0];

      // QCD recoiler: emergency catch of mismatches e.g. when gluinos.
      if (iRec == 0 && dipNow.colType > 0 ) {
        for (int i = 2; i < 2 + sizeOut; ++i)
        if (event[iNew[i]].acol() == colRad) {
          iRec = iNew[i];
          dipNow.isrType = 0;
          break;
        }
      }
      if (iRec == 0 && dipNow.colType < 0 ) {
        for (int i = 2; i < 2 + sizeOut; ++i)
        if (event[iNew[i]].col() == acolRad) {
          iRec = iNew[i];
          dipNow.isrType = 0;
          break;
        }
      }

      // QED/photon recoiler: either to new particle or remains to beam.
      if ( iRec == 0 && (dipNow.chgType != 0 || dipNow.gamType != 0) ) {
        if ( event[iNew[0]].chargeType() == 0 ) {
          iRec = iNewNew;
          dipNow.isrType = 0;
        } else {
          iRec = iNew[0];
        }
      }

    // Recoiler in another system: keep it as is.
    } else iRec = dipNow.iRecoiler;

    // Done. Kill dipole if failed to find new recoiler.
    dipNow.iRecoiler = iRec;
    if ( iRec == 0 && (dipNow.colType != 0 || dipNow.chgType != 0
      || dipNow.gamType != 0) ) {
      dipNow.colType = 0;
      dipNow.chgType = 0;
      dipNow.gamType = 0;
      infoPtr->errorMsg("Error in SimpleTimeShower::update: "
      "failed to locate new recoiling partner");
    }

    // Kill weak dipoles if ISR emitted W/Z
    // and only a single weak emission is allowed.
    if (hasWeakRad && singleWeakEmission && dipNow.weakType != 0)
      dipNow.weakType = 0;
  }

  // Set the weak radiated variable to true if already radiated.
  if (hasWeakRad) hasWeaklyRadiated = true;

  // Find new dipole end formed by colour index.
  int colTag = event[iNewNew].col();
  if (doQCDshower && colTag > 0)
    setupQCDdip( iSys, sizeOut, colTag, 1, event, false, true);

  // Find new dipole end formed by anticolour index.
  int acolTag = event[iNewNew].acol();
  if (doQCDshower && acolTag > 0)
    setupQCDdip( iSys, sizeOut, acolTag, -1, event, false, true);


  // Find new "charge-dipole" and "photon-dipole" ends.
  int  chgType  = event[iNewNew].chargeType();
  bool doChgDip = (chgType != 0)
                  && ( ( doQEDshowerByQ && event[iNewNew].isQuark()  )
                    || ( doQEDshowerByL && event[iNewNew].isLepton() ) );
  int  gamType  = (event[iNewNew].id() == 22) ? 1 : 0;
  bool doGamDip = (gamType == 1) && doQEDshowerByGamma;
  if (doChgDip || doGamDip)
    setupQEDdip( iSys, sizeOut, chgType, gamType, event, true);

  // Find new weak dipole.
  // Uses the size of dipEnd to tell whether a new dipole is added.
  unsigned int nDips = dipEnd.size();
  if (doWeakShower && (event[iNewNew].isQuark() || event[iNewNew].isLepton())
      && !(hasWeaklyRadiated && singleWeakEmission)
      && (iSys == 0 || !partonSystemsPtr->hasInAB(iSys))) {

    if (weakMode == 0 || weakMode == 1)
      setupWeakdip( iSys, sizeOut, 1, event, true);
    // If added new dipole update the ME correction and me partner.
    if (nDips != dipEnd.size()) {
      nDips = dipEnd.size();
      dipEnd.back().MEtype = 200;
      dipEnd.back().iMEpartner = dipEnd.back().iRecoiler;
    }

    if (weakMode == 0 || weakMode == 2)
      setupWeakdip( iSys, sizeOut, 2, event, true);
    // If added new dipole, update the ME correction and me partner.
    if (nDips != dipEnd.size()) {
      nDips = dipEnd.size();
      dipEnd.back().MEtype = 205;
      dipEnd.back().iMEpartner = dipEnd.back().iRecoiler;
    }
  }

  // Start iterate over list of rescatterers - may be empty.
  int iRescNow = -1;
  while (++iRescNow < int(iRescatterer.size())) {

    // Identify systems that rescatterers belong to.
    int iOutNew    = iRescatterer[iRescNow];
    int iInNew     = event[iOutNew].daughter1();
    int iSysResc   = partonSystemsPtr->getSystemOf(iInNew, true);

    // Find new and old positions of incoming partons in the system.
    iNew.resize(0);
    iOld.resize(0);
    iNew.push_back( partonSystemsPtr->getInA(iSysResc) );
    iOld.push_back( event[iNew[0]].daughter1() );
    iNew.push_back( partonSystemsPtr->getInB(iSysResc) );
    iOld.push_back( event[iNew[1]].daughter1() );

    // Ditto for outgoing partons.
    sizeOut = partonSystemsPtr->sizeOut(iSysResc);
    for (int i = 0; i < sizeOut; ++i) {
      int iNow = partonSystemsPtr->getOut(iSysResc, i);
      iNew.push_back( iNow );
      iOld.push_back( event[iNow].mother1() );
      // Add non-final to list of rescatterers.
      if (!event[iNow].isFinal()) iRescatterer.push_back( iNow );
    }

    // Loop over all dipole ends belonging to the system
    // or to the recoil system, if different.
    for (int iDip = 0; iDip < int(dipEnd.size()); ++iDip)
    if (dipEnd[iDip].system == iSysResc
      || dipEnd[iDip].systemRec == iSysResc) {
      TimeDipoleEnd& dipNow = dipEnd[iDip];

      // Replace radiator (always in final state so simple).
      for (int i = 2; i < 2 + sizeOut; ++i)
      if (dipNow.iRadiator == iOld[i]) {
        dipNow.iRadiator = iNew[i];
        break;
      }

      // Replace ME partner (always in final state, if exists, so simple).
      for (int i = 2; i < 2 + sizeOut; ++i)
      if (dipNow.iMEpartner == iOld[i]) {
        dipNow.iMEpartner = iNew[i];
        break;
      }

      // Replace recoiler.
      for (int i = 0; i < 2 + sizeOut; ++i)
      if (dipNow.iRecoiler == iOld[i]) {
        dipNow.iRecoiler = iNew[i];
        break;
      }
    }

  // End iterate over list of rescatterers.
  }

}

//--------------------------------------------------------------------------

// Setup a dipole end for a QCD colour charge.

void SimpleTimeShower::setupQCDdip( int iSys, int i, int colTag, int colSign,
  Event& event, bool isOctetOnium, bool limitPTmaxIn) {

  // Initial values.
  int iRad     = partonSystemsPtr->getOut(iSys, i);
  int iRec     = 0;
  int sizeAll  = partonSystemsPtr->sizeAll(iSys);
  int sizeOut  = partonSystemsPtr->sizeOut(iSys);
  // Number of potential recoilers; decide if beams included or not.
  int sizeRec  = ( allowBeamRecoil && partonSystemsPtr->hasInAB(iSys) ) ?
    sizeAll : sizeOut;
  int sizeInRec    = sizeRec - sizeOut;
  int sizeInNonRec = sizeAll - sizeOut - sizeInRec;
  int iOffset  = i + sizeAll - sizeOut;
  bool otherSystemRec = false;
  bool allowInitial   = partonSystemsPtr->hasInAB(iSys);
  // PS dec 2010: possibility to allow for several recoilers and each with
  // flexible normalization
  bool   isFlexible   = false;
  double flexFactor   = 1.0;
  vector<int> iRecVec(0);

  // Colour: other end by same index in beam or opposite in final state.
  // Exclude rescattered incoming and not final outgoing.
  if (colSign > 0) {
    for (int j = 0; j < sizeRec; ++j) {
      if (j + sizeInNonRec != iOffset) {
        int iRecNow = partonSystemsPtr->getAll(iSys, j + sizeInNonRec);
        if ( ( j <  sizeInRec && event[iRecNow].col()  == colTag
               && !event[iRecNow].isRescatteredIncoming() )
             || ( j >= sizeInRec && event[iRecNow].acol() == colTag
                  && event[iRecNow].isFinal() ) ) {
          iRec = iRecNow;
          break;
        }
      }
    }
  }

  // Anticolour: other end by same index in beam or opposite in final state.
  // Exclude rescattered incoming and not final outgoing.
  if (colSign < 0) {
    for (int j = 0; j < sizeRec; ++j) {
      if (j + sizeInNonRec != iOffset) {
        int iRecNow = partonSystemsPtr->getAll(iSys, j + sizeInNonRec);
        if ( ( j <  sizeInRec && event[iRecNow].acol()  == colTag
               && !event[iRecNow].isRescatteredIncoming() )
             || ( j >= sizeInRec && event[iRecNow].col() == colTag
                  && event[iRecNow].isFinal() ) ) {
          iRec = iRecNow;
          break;
        }
      }
    }
  }

  // Resonance decays:
  // other end to nearest recoiler in same system final state,
  // by (p_i + p_j)^2 - (m_i + m_j)^2 = 2 (p_i p_j - m_i m_j).
  // (junction colours more involved, so keep track if junction colour)
  bool hasJunction = false;
  if (iRec == 0 && !allowInitial) {
    for (int iJun = 0; iJun < event.sizeJunction(); ++ iJun) {
      // For types 1&2, all legs in final state
      // For types 3&4, two legs in final state
      // For types 5&6, one leg in final state
      int iBeg = (event.kindJunction(iJun)-1)/2;
      for (int iLeg = iBeg; iLeg < 3; ++iLeg)
        if (event.endColJunction( iJun, iLeg) == colTag) hasJunction  = true;
    }
    double ppMin = LARGEM2;
    for (int j = 0; j < sizeOut; ++j) {
      if (j != i) {
        int iRecNow  = partonSystemsPtr->getOut(iSys, j);
        if (!event[iRecNow].isFinal()) continue;
        double ppNow = event[iRecNow].p() * event[iRad].p()
          - event[iRecNow].m() * event[iRad].m();
        if (ppNow < ppMin) {
          iRec  = iRecNow;
          ppMin = ppNow;
        }
      }
    }
  }

  // If no success then look for matching (anti)colour anywhere in final state.
  if ( iRec == 0 || (!doInterleave && allowMPIdipole
    && !event[iRec].isFinal()) ) {
    for (int j = 0; j < event.size(); ++j) {
      if (event[j].isFinal()) {
        if ( (colSign > 0 && event[j].acol() == colTag)
             || (colSign < 0 && event[j].col()  == colTag) ) {
          iRec = j;
          otherSystemRec = true;
          break;
        }
      }
    }

    // If no success then look for match to non-rescattered in initial state.
    if (iRec == 0 && allowInitial) {
      for (int iSysR = 0; iSysR < partonSystemsPtr->sizeSys(); ++iSysR)
      if (iSysR != iSys) {
        int j = partonSystemsPtr->getInA(iSysR);
        if (j > 0 && event[j].isRescatteredIncoming()) j = 0;
        if (j > 0 && ( (colSign > 0 && event[j].col() == colTag)
          || (colSign < 0 && event[j].acol()  == colTag) ) ) {
          iRec = j;
          otherSystemRec = true;
          break;
        }
        j = partonSystemsPtr->getInB(iSysR);
        if (j > 0 && event[j].isRescatteredIncoming()) j = 0;
        if (j > 0 && ( (colSign > 0 && event[j].col() == colTag)
          || (colSign < 0 && event[j].acol()  == colTag) ) ) {
          iRec = j;
          otherSystemRec = true;
          break;
        }
      }
    }
  }

  // Junctions (PS&ND dec 2010)
  // For types 1&2: all legs in final state
  //                half-strength dipoles between all legs
  // For types 3&4, two legs in final state
  //                full-strength dipole between final-state legs
  // For types 5&6, one leg in final state
  //                no final-state dipole end

  if (hasJunction) {
    for (int iJun = 0; iJun < event.sizeJunction(); ++ iJun) {
      int kindJun = event.kindJunction(iJun);
      int iBeg = (kindJun-1)/2;
      for (int iLeg = iBeg; iLeg < 3; ++iLeg) {
        if (event.endColJunction( iJun, iLeg) == colTag) {
          // For types 5&6, no other leg to recoil against. Switch off if
          // no other particles at all, since radiation then handled by ISR.
          // Example: qq -> ~t* : no radiation off ~t*
          // Allow radiation + recoil if unconnected partners available
          // Example: qq -> ~t* -> tbar ~chi0 : allow radiation off tbar,
          //                                    with ~chi0 as recoiler
          if (kindJun >= 5) {
            if (sizeOut == 1) return;
            else break;
          }
          // For junction types 3 & 4, span one full-strength dipole
          // (only look inside same decay system)
          else if (kindJun >= 3) {
            int iLegRec = 3-iLeg;
            int colTagRec = event.endColJunction( iJun, iLegRec);
            for (int j = 0; j < sizeOut; ++j) if (j != i) {
                int iRecNow  = partonSystemsPtr->getOut(iSys, j);
                if (!event[iRecNow].isFinal()) continue;
                if ( (colSign > 0 && event[iRecNow].col()  == colTagRec)
                  || (colSign < 0 && event[iRecNow].acol() == colTagRec) ) {
                  // Only accept if staying inside same system
                  iRec = iRecNow;
                  break;
                }
              }
          }
          // For junction types 1 & 2, span two half-strength dipoles
          // (only look inside same decay system)
          else {
            // Loop over two half-strength dipole connections
            for (int jLeg = 1; jLeg <= 2; jLeg++) {
              int iLegRec = (iLeg + jLeg) % 3;
              int colTagRec = event.endColJunction( iJun, iLegRec);
              for (int j = 0; j < sizeOut; ++j) if (j != i) {
                  int iRecNow  = partonSystemsPtr->getOut(iSys, j);
                  if (!event[iRecNow].isFinal()) continue;
                  if ( (colSign > 0 && event[iRecNow].col()  == colTagRec)
                    || (colSign < 0 && event[iRecNow].acol() == colTagRec) ) {
                    // Store recoilers in temporary array
                    iRecVec.push_back(iRecNow);
                    // Set iRec != 0 for checks below
                    iRec = iRecNow;
                  }
                }
            }

          }     // End if-then-else of junction kinds

        }       // End if leg has right colour tag
      }         // End of loop over junction legs
    }           // End loop over junctions

  }             // End main junction if

  // If fail, then other end to nearest recoiler in same system final state,
  // by (p_i + p_j)^2 - (m_i + m_j)^2 = 2 (p_i p_j - m_i m_j).
  if (iRec == 0) {
    double ppMin = LARGEM2;
    for (int j = 0; j < sizeOut; ++j) if (j != i) {
      int iRecNow  = partonSystemsPtr->getOut(iSys, j);
      if (!event[iRecNow].isFinal()) continue;
      double ppNow = event[iRecNow].p() * event[iRad].p()
                   - event[iRecNow].m() * event[iRad].m();
      if (ppNow < ppMin) {
        iRec  = iRecNow;
        ppMin = ppNow;
      }
    }
  }

  // If fail, then other end to nearest recoiler in any system final state,
  // by (p_i + p_j)^2 - (m_i + m_j)^2 = 2 (p_i p_j - m_i m_j).
  if (iRec == 0) {
    double ppMin = LARGEM2;
    for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow)
    if (iRecNow != iRad && event[iRecNow].isFinal()) {
      double ppNow = event[iRecNow].p() * event[iRad].p()
                   - event[iRecNow].m() * event[iRad].m();
      if (ppNow < ppMin) {
        iRec  = iRecNow;
        otherSystemRec = true;
        ppMin = ppNow;
      }
    }
  }

  // PS dec 2010: make sure iRec is stored in iRecVec
  if (iRecVec.size() == 0 && iRec != 0) iRecVec.push_back(iRec);

  // Remove any zero recoilers from normalization
  int nRec = iRecVec.size();
  for (unsigned int mRec = 0; mRec < iRecVec.size(); ++mRec)
    if (iRecVec[mRec] <= 0) nRec--;
  if (nRec >= 2) {
    isFlexible = true;
    flexFactor = 1.0/nRec;
  }

  // Check for failure to locate any recoiler
  if ( nRec <= 0 ) {
    infoPtr->errorMsg("Error in SimpleTimeShower::setupQCDdip: "
                      "failed to locate any recoiling partner");
    return;
  }

  // Store dipole colour end(s).
  for (unsigned int mRec = 0; mRec < iRecVec.size(); ++mRec) {
    iRec = iRecVec[mRec];
    if (iRec <= 0) continue;
    // Max scale either by parton scale or by half dipole mass.
    double pTmax = event[iRad].scale();
    if (limitPTmaxIn) {
      if (iSys == 0 || (iSys == 1 && twoHard)) pTmax *= pTmaxFudge;
      else if (sizeInRec > 0) pTmax *= pTmaxFudgeMPI;
    } else pTmax = 0.5 * m( event[iRad], event[iRec]);
    int colType  = (event[iRad].id() == 21) ? 2 * colSign : colSign;
    int isrType  = (event[iRec].isFinal()) ? 0 : event[iRec].mother1();
    // This line in case mother is a rescattered parton.
    while (isrType > 2 + beamOffset) isrType = event[isrType].mother1();
    if (isrType > 2) isrType -= beamOffset;
    dipEnd.push_back( TimeDipoleEnd( iRad, iRec, pTmax,
      colType, 0, 0, 0, isrType, iSys, -1, -1, 0, isOctetOnium) );

    // If hooked up with other system then find which.
    if (otherSystemRec) {
      int systemRec = partonSystemsPtr->getSystemOf(iRec, true);
      if (systemRec >= 0) dipEnd.back().systemRec = systemRec;
      dipEnd.back().MEtype = 0;
    }

    // PS dec 2010
    // If non-unity (flexible) normalization, set normalization factor
    if (isFlexible) {
      dipEnd.back().isFlexible = true;
      dipEnd.back().flexFactor = flexFactor;
    }
  }

}

//--------------------------------------------------------------------------

// Setup a dipole end for a QED colour charge or a photon.
// No failsafe choice of recoiler, so gradually widen search.

void SimpleTimeShower::setupQEDdip( int iSys, int i, int chgType, int gamType,
  Event& event, bool limitPTmaxIn) {

  // Initial values. Find if allowed to hook up beams.
  int iRad     = partonSystemsPtr->getOut(iSys, i);
  int idRad    = event[iRad].id();
  int iRec     = 0;
  int sizeAll  = partonSystemsPtr->sizeAll(iSys);
  int sizeOut  = partonSystemsPtr->sizeOut(iSys);
  // Number of potential recoilers; decide if beams included or not.
  int sizeRec  = ( allowBeamRecoil && partonSystemsPtr->hasInAB(iSys) )
    ? sizeAll : sizeOut;
  int sizeInRec    = sizeRec - sizeOut;
  int sizeInNonRec = sizeAll - sizeOut - sizeInRec;
  int iOffset  = i + sizeAll - sizeOut;
  double ppMin = LARGEM2;
  bool hasRescattered = false;
  bool otherSystemRec = false;

  // Find nearest same- (opposide-) flavour recoiler in initial (final)
  // state of same system, excluding rescattered (in or out) partons.
  // Also find if system is involved in rescattering.
  // Note: (p_i + p_j)2 - (m_i + m_j)2 = 2 (p_i p_j - m_i m_j).
  for (int j = 0; j < sizeRec; ++j) {
    if (j + sizeInNonRec != iOffset) {
      int iRecNow  = partonSystemsPtr->getAll(iSys, j + sizeInNonRec);
      if ( (j <  sizeInRec && !event[iRecNow].isRescatteredIncoming())
           || (j >= sizeInRec && event[iRecNow].isFinal()) ) {
        if ( (j <  sizeInRec && event[iRecNow].id() ==  idRad)
             || (j >= sizeInRec && event[iRecNow].id() == -idRad) ) {
          double ppNow = event[iRecNow].p() * event[iRad].p()
            - event[iRecNow].m() * event[iRad].m();
          if (ppNow < ppMin) {
            iRec  = iRecNow;
            ppMin = ppNow;
          }
        }
      } else hasRescattered = true;
    }
  }

  // If rescattering then find nearest opposite-flavour recoiler
  // anywhere in final state.
  if (iRec == 0 && hasRescattered) {
    for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow)
    if (event[iRecNow].id() == -idRad && event[iRecNow].isFinal()) {
      double ppNow = event[iRecNow].p() * event[iRad].p()
                   - event[iRecNow].m() * event[iRad].m();
      if (ppNow < ppMin) {
        iRec  = iRecNow;
        ppMin = ppNow;
        otherSystemRec = true;
      }
    }
  }

  // Find nearest recoiler in same system, charge-squared-weighted,
  // including initial state, but excluding rescatterer.
  if (iRec == 0) {
    for (int j = 0; j < sizeRec; ++j) {
      if (j + sizeInNonRec != iOffset) {
        int iRecNow = partonSystemsPtr->getAll(iSys, j + sizeInNonRec);
        int chgTypeRecNow = event[iRecNow].chargeType();
        if (chgTypeRecNow == 0) continue;
        if ( (j <  sizeInRec && !event[iRecNow].isRescatteredIncoming())
             || (j >= sizeInRec && event[iRecNow].isFinal()) ) {
          double ppNow = (event[iRecNow].p() * event[iRad].p()
                          -  event[iRecNow].m() * event[iRad].m())
            / pow2(chgTypeRecNow);
          if (ppNow < ppMin) {
            iRec  = iRecNow;
            ppMin = ppNow;
          }
        }
      }
    }
  }

  // If rescattering then find nearest recoiler in the final state,
  // charge-squared-weighted.
  if (iRec == 0 && hasRescattered) {
    for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow)
    if (iRecNow != iRad && event[iRecNow].isFinal()) {
      int chgTypeRecNow = event[iRecNow].chargeType();
      if (chgTypeRecNow != 0 && event[iRecNow].isFinal()) {
        double ppNow = (event[iRecNow].p() * event[iRad].p()
                     -  event[iRecNow].m() * event[iRad].m())
                     / pow2(chgTypeRecNow);
        if (ppNow < ppMin) {
          iRec  = iRecNow;
          ppMin = ppNow;
          otherSystemRec = true;
        }
      }
    }
  }

  // Find any nearest recoiler in final state of same system.
  if (iRec == 0) {
    for (int j = 0; j < sizeOut; ++j) {
      if (j != i) {
        int iRecNow  = partonSystemsPtr->getOut(iSys, j);
        double ppNow = event[iRecNow].p() * event[iRad].p()
          - event[iRecNow].m() * event[iRad].m();
        if (ppNow < ppMin) {
          iRec  = iRecNow;
          ppMin = ppNow;
        }
      }
    }
  }

  // Find any nearest recoiler in final state.
  if (iRec == 0) {
    for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow) {
      if (iRecNow != iRad && event[iRecNow].isFinal()) {
        double ppNow = event[iRecNow].p() * event[iRad].p()
          - event[iRecNow].m() * event[iRad].m();
        if (ppNow < ppMin) {
          iRec  = iRecNow;
          ppMin = ppNow;
          otherSystemRec = true;
        }
      }
    }
  }

  // Fill charge-dipole or photon-dipole end.
  if (iRec > 0) {
    // Max scale either by parton scale or by half dipole mass.
    double pTmax = event[iRad].scale();
    if (limitPTmaxIn) {
      if (iSys == 0 || (iSys == 1 && twoHard)) pTmax *= pTmaxFudge;
      else if (sizeInRec > 0) pTmax *= pTmaxFudgeMPI;
    } else pTmax = 0.5 * m( event[iRad], event[iRec]);
    int isrType = (event[iRec].isFinal()) ? 0 : event[iRec].mother1();
    // This line in case mother is a rescattered parton.
    while (isrType > 2 + beamOffset) isrType = event[isrType].mother1();
    if (isrType > 2) isrType -= beamOffset;
    dipEnd.push_back( TimeDipoleEnd(iRad, iRec, pTmax,
      0, chgType, gamType, 0, isrType, iSys, -1) );

    // If hooked up with other system then find which.
    if (otherSystemRec) {
      int systemRec = partonSystemsPtr->getSystemOf(iRec);
      if (systemRec >= 0) dipEnd.back().systemRec = systemRec;
      dipEnd.back().MEtype = 0;
    }

  // Failure to find other end of dipole.
  } else {
    infoPtr->errorMsg("Error in SimpleTimeShower::setupQEDdip: "
      "failed to locate any recoiling partner");
  }

}

//--------------------------------------------------------------------------

 // Setup a dipole end for weak W or Z emission.

void SimpleTimeShower::setupWeakdip( int iSys, int i, int weakType,
  Event& event, bool limitPTmaxIn) {

  // Initial values. Find if allowed to hook up beams.
  int iRad     = partonSystemsPtr->getOut(iSys, i);
  int idRad    = event[iRad].id();
  int iRec     = 0;
  int sizeAll  = partonSystemsPtr->sizeAll(iSys);
  int sizeOut  = partonSystemsPtr->sizeOut(iSys);
  // Only allow weak dipoles to take outgoing particles as recoiler.
  int sizeRec  = sizeOut;
  int sizeInRec     = sizeRec - sizeOut;
  int sizeInNonRec  = sizeAll - sizeInRec - sizeOut;
  int iOffset  = i + sizeAll - sizeOut;
  double ppMin = LARGEM2;
  bool hasRescattered = false;
  bool otherSystemRec = false;

  // Find nearest same- (opposide-) flavour recoiler in initial (final)
  // state of same system, excluding rescattered (in or out) partons.
  // Also find if system is involved in rescattering.
  // Note: (p_i + p_j)2 - (m_i + m_j)2 = 2 (p_i p_j - m_i m_j).
  for (int j = 0; j < sizeRec; ++j)
    if (j + sizeInNonRec != iOffset) {
      int iRecNow  = partonSystemsPtr->getAll(iSys, j + sizeInNonRec);
      if ( (j <  sizeInRec && !event[iRecNow].isRescatteredIncoming())
        || (j >= sizeInRec && event[iRecNow].isFinal()) ) {
        if ( (j <  sizeInRec && event[iRecNow].id() ==  idRad)
          || (j >= sizeInRec && event[iRecNow].id() == -idRad) ) {
          double ppNow = event[iRecNow].p() * event[iRad].p()
                       - event[iRecNow].m() * event[iRad].m();
          if (ppNow < ppMin) {
          iRec  = iRecNow;
          ppMin = ppNow;
          }
        }
      } else hasRescattered = true;
    }

  // If rescattering then find nearest opposite-flavour recoiler
  // anywhere in final state.
  if (iRec == 0 && hasRescattered) {
    for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow)
    if (event[iRecNow].id() == -idRad && event[iRecNow].isFinal()) {
      double ppNow = event[iRecNow].p() * event[iRad].p()
                   - event[iRecNow].m() * event[iRad].m();
      if (ppNow < ppMin) {
        iRec  = iRecNow;
        ppMin = ppNow;
        otherSystemRec = true;
      }
    }
  }

  // Find nearest recoiler in same system, weak-charge-squared-weighted,
  // including initial state, but excluding rescatterer.
  if (iRec == 0)
  for (int j = 0; j < sizeRec; ++j) if (j + sizeInNonRec != iOffset) {
    int iRecNow = partonSystemsPtr->getAll(iSys, j + sizeInNonRec);
    if (abs(event[iRecNow].id()) >= 20 || weakType < 1
      || weakType > 2) continue;
    double weakCoupNow = 1.;
    if (weakType == 2) weakCoupNow = coupSMPtr->vf2(event[iRecNow].idAbs())
      + coupSMPtr->af2(event[iRecNow].idAbs());
    if ( (j <  sizeInRec && !event[iRecNow].isRescatteredIncoming())
      || (j >= sizeInRec && event[iRecNow].isFinal()) ) {
      double ppNow = (event[iRecNow].p() * event[iRad].p()
                   -  event[iRecNow].m() * event[iRad].m()) / weakCoupNow;
      if (ppNow < ppMin) {
        iRec  = iRecNow;
        ppMin = ppNow;
      }
    }
  }

  // If rescattering then find nearest recoiler in the final state,
  // weak-charge-squared-weighted.
  if (iRec == 0 && hasRescattered) {
    for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow)
    if (iRecNow != iRad && event[iRecNow].isFinal()) {
      if (abs(event[iRecNow].id()) >= 20 || weakType < 1
        || weakType > 2) continue;
      double weakCoupNow = 1.;
      if (weakType == 2) weakCoupNow = coupSMPtr->vf2(event[iRecNow].idAbs())
        + coupSMPtr->af2(event[iRecNow].idAbs());
      double ppNow = (event[iRecNow].p() * event[iRad].p()
                   - event[iRecNow].m() * event[iRad].m()) / weakCoupNow;
      if (ppNow < ppMin) {
        iRec  = iRecNow;
        ppMin = ppNow;
        otherSystemRec = true;
      }
    }
  }

  // Find any nearest recoiler in final state of same system.
  if (iRec == 0)
  for (int j = 0; j < sizeOut; ++j) if (j != i) {
    int iRecNow  = partonSystemsPtr->getOut(iSys, j);
    double ppNow = event[iRecNow].p() * event[iRad].p()
                 - event[iRecNow].m() * event[iRad].m();
    if (ppNow < ppMin) {
      iRec  = iRecNow;
      ppMin = ppNow;
    }
  }

  // Find any nearest recoiler in final state.
  if (iRec == 0)
  for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow)
  if (iRecNow != iRad && event[iRecNow].isFinal()) {
    double ppNow = event[iRecNow].p() * event[iRad].p()
                 - event[iRecNow].m() * event[iRad].m();
    if (ppNow < ppMin) {
      iRec  = iRecNow;
      ppMin = ppNow;
      otherSystemRec = true;
    }
  }

  // Fill in weak dipole-end.
  if (iRec > 0) {

    // Calculate 2 -> 2 kinematics, needed for finding ISR fermion line.
    Vec4 p3weak = event[3].p();
    Vec4 p4weak = event[4].p();
    double tHat = (event[iRad].p() - p3weak).m2Calc();
    double uHat = (event[iRad].p() - p4weak).m2Calc();

    // Find correct helicity.
    int weakPol = (rndmPtr->flat() > 0.5) ? -1 : 1;
    // Check if particle has already gotten a helicity.
    if (event[iRad].intPol() == 1 || event[iRad].intPol() == -1)
      weakPol = event[iRad].intPol();
    // If particle come from ISR radiation.
    else if (event[iRad].statusAbs() > 40) {
      if (event[event[iRad].mother1()].idAbs() < 20)
        weakPol = event[event[iRad].mother1()].intPol();
      else if (int(event[iRad].sisterList(true).size()) != 0)
        weakPol = event[event[iRad].sisterList(true)[0]].intPol();
    }
    // If it is not a 2 to 2 process, always use recoiler.
    else if (infoPtr->nFinal() != 2) {
      if (event[iRec].intPol() == 1 || event[iRec].intPol() == -1)
        weakPol = event[iRec].intPol();
    }
    // If s-channel, choose same spin as recoiler.
    else if (idRad == - event[iRec].id()) {
      if (event[iRec].intPol() == 1 || event[iRec].intPol() == -1)
        weakPol = event[iRec].intPol();
    }
    // if W-decay, choose always left handed.
    else if (event[event[iRad].mother1()].idAbs() == 24) weakPol = -1;
    // If four particles of the same type.
    else if (idRad == event[iRec].id()) {
      if (uHat*uHat/(tHat*tHat + uHat*uHat) > 0.5) weakPol = event[3].intPol();
      else weakPol = event[4].intPol();
    }
    // For different particle types, choose correct fermion line.
    else if (event[3].id() == idRad) weakPol = event[3].intPol();
    else if (event[4].id() == idRad) weakPol = event[4].intPol();
    // If weak ISR is turned off, this would try to use polarization
    // that is not set as expected. In this case use random polarization.
    if (weakPol > 1) weakPol = (rndmPtr->flat() > 0.5) ? -1 : 1;
    event[iRad].pol(weakPol);

    // Max scale either by parton scale or by half dipole mass.
    double pTmax = event[iRad].scale();
    if (limitPTmaxIn) {
      if (iSys == 0) pTmax *= pTmaxFudge;
      if (iSys > 0 && sizeInRec > 0) pTmax *= pTmaxFudgeMPI;
    } else pTmax = 0.5 * m( event[iRad], event[iRec]);
    int isrType = (event[iRec].isFinal()) ? 0 : event[iRec].mother1();
    // This line in case mother is a rescattered parton.
    while (isrType > 2 + beamOffset) isrType = event[isrType].mother1();
    if (isrType > 2) isrType -= beamOffset;
    // No right-handed W emission.

    if (weakType == 1 && weakPol == 1) return;
    dipEnd.push_back( TimeDipoleEnd(iRad, iRec, pTmax,
      0, 0, 0, weakType, isrType, iSys, -1, -1, weakPol) );

    // If hooked up with other system then find which.
    if (otherSystemRec) {
      int systemRec = partonSystemsPtr->getSystemOf(iRec);
      if (systemRec >= 0) dipEnd.back().systemRec = systemRec;
      dipEnd.back().MEtype = 0;
    }

  // Failure to find other end of dipole.
  } else {
    infoPtr->errorMsg("Error in SimpleTimeShower::setupWeakdip: "
      "failed to locate any recoiling partner");
  }
}

//--------------------------------------------------------------------------

// Special setup for weak dipoles if already specified in info ptr.
void SimpleTimeShower::setupWeakdipExternal(Event& event, bool limitPTmaxIn) {

  // Get information.
  vector<pair<int,int> > weakDipoles = infoPtr->getWeakDipoles();
  vector<int> weakModes = infoPtr->getWeakModes();
  weakMomenta = infoPtr->getWeakMomenta();
  weak2to2lines = infoPtr->getWeak2to2lines();
  weakHardSize = int(weakModes.size());

  // Loop over dipoles.
  for (int i = 0; i < int(weakDipoles.size()); ++i) {
    // Only consider FSR dipoles.
    if (event[weakDipoles[i].first].status() > 0) {
      // Find ME.
      int iRad = weakDipoles[i].first;
      int iRec = weakDipoles[i].second;

      // Find MEtype.
      int MEtypeWeak = 0;
      if (weakModes[weakDipoles[i].first] == 1) MEtypeWeak = 200;
      else if (weakModes[weakDipoles[i].first] == 2) MEtypeWeak = 201;
      else if (weakModes[weakDipoles[i].first] == 3) MEtypeWeak = 202;
      else MEtypeWeak = 203;

      // Find correct polarization, if it is already set use it.
      // Otherwise pick randomly.
      int weakPol = (rndmPtr->flat() > 0.5) ? -1 : 1;
      if (event[weakDipoles[i].first].intPol() != 9)
        weakPol = event[weakDipoles[i].first].intPol();
      else if (event[weakDipoles[i].second].intPol() != 9) {
        if (event[weakDipoles[i].second].status() < 0)
          weakPol = event[weakDipoles[i].second].intPol();
        else
          weakPol = -event[weakDipoles[i].second].intPol();
      }
      event[weakDipoles[i].first].pol(weakPol);

      // Max scale either by parton scale or by half dipole mass.
      double pTmax = event[iRad].scale();

      if (limitPTmaxIn) {
        pTmax *= pTmaxFudge;
      } else pTmax = 0.5 * m( event[iRad], event[iRec]);

      // Recoiler is always final state.
      int isrType = 0;

      // No right-handed W emission.
      // Add the dipoles.
      if ( (weakMode == 0 || weakMode == 1) && weakPol == -1)
        dipEnd.push_back( TimeDipoleEnd(iRad, iRec, pTmax,
          0, 0, 0, 1, isrType, 0, MEtypeWeak, -1, weakPol) );

      if (weakMode == 0 || weakMode == 2)
         dipEnd.push_back( TimeDipoleEnd(iRad, iRec, pTmax,
          0, 0, 0, 2, isrType, 0, MEtypeWeak +5, -1, weakPol) );

    }
  }

  for (int i = 0;i < int(dipEnd.size()); ++i) {
    Vec4 p3weak, p4weak;
    if (dipEnd[i].MEtype > 200) {
      int i2to2Mother = dipEnd[i].iRadiator;
      while (i2to2Mother >= weakHardSize)
        i2to2Mother = event[i2to2Mother].mother1();
      if (weak2to2lines[2] == i2to2Mother) {
        p3weak = weakMomenta[0];
        p4weak = weakMomenta[1];
      } else {
        p3weak = weakMomenta[1];
        p4weak = weakMomenta[0];
      }
    }
  }

}

//--------------------------------------------------------------------------

// Setup a dipole end for a Hidden Valley colour charge.

void SimpleTimeShower::setupHVdip( int iSys, int i, Event& event,
  bool limitPTmaxIn) {

  // Initial values.
  int iRad    = partonSystemsPtr->getOut(iSys, i);
  int iRec    = 0;
  int idRad   = event[iRad].id();
  int sizeOut = partonSystemsPtr->sizeOut(iSys);

  // Hidden Valley colour positive for positive id, and vice versa.
  // Find opposte HV colour in final state of same system.
  for (int j = 0; j < sizeOut; ++j) if (j != i) {
    int iRecNow = partonSystemsPtr->getOut(iSys, j);
    int idRec   = event[iRecNow].id();
    if ( (abs(idRec) > 4900000 && abs(idRec) < 4900017)
      && idRad * idRec < 0) {
      iRec = iRecNow;
      break;
    }
  }

  // Else find heaviest other final-state in same system.
  // (Intended for decays; should mainly be two-body so unique.)
  double mMax = -sqrt(LARGEM2);
   if (iRec == 0)
  for (int j = 0; j < sizeOut; ++j) if (j != i) {
    int iRecNow = partonSystemsPtr->getOut(iSys, j);
    if (event[iRecNow].m() > mMax) {
      iRec = iRecNow;
      mMax = event[iRecNow].m();
    }
  }

  // Set up dipole end, or report failure.
  if (iRec > 0) {
    // Max scale either by parton scale or by half dipole mass.
    double pTmax = event[iRad].scale();
    if (limitPTmaxIn) {
      if (iSys == 0 || (iSys == 1 && twoHard)) pTmax *= pTmaxFudge;
    } else pTmax = 0.5 * m( event[iRad], event[iRec]);
    int colvType  = (event[iRad].id() > 0) ? 1 : -1;
    dipEnd.push_back( TimeDipoleEnd( iRad, iRec, pTmax, 0, 0, 0, 0, 0,
      iSys, -1, -1, 0, false, true, colvType) );
  } else infoPtr->errorMsg("Error in SimpleTimeShower::setupHVdip: "
      "failed to locate any recoiling partner");

}

//--------------------------------------------------------------------------

// Select next pT in downwards evolution of the existing dipoles.

double SimpleTimeShower::pTnext( Event& event, double pTbegAll,
  double pTendAll, bool isFirstTrial, bool doTrialIn) {

  // Begin loop over all possible radiating dipole ends.
  dipSel  = 0;
  iDipSel = -1;
  double pT2sel = pTendAll * pTendAll;

  // Check if enhanced emissions should be applied.
  doTrialNow    = doTrialIn;
  canEnhanceET  = (!doTrialNow && canEnhanceEmission)
               || ( doTrialNow && canEnhanceTrial);

  // Starting values for enhanced emissions.
  splittingNameSel = "";
  splittingNameNow = "";
  enhanceFactors.clear();
  if (hasUserHooks) userHooksPtr->setEnhancedTrial(0., 1.);

  for (int iDip = 0; iDip < int(dipEnd.size()); ++iDip) {
    TimeDipoleEnd& dip = dipEnd[iDip];
    dip.pAccept        = 1.0;

    // Check if this system is part of the hard scattering
    // (including resonance decay products).
    bool hardSystem = partonSystemsPtr->getHard(dip.system);
    bool isQCD = event[dip.iRadiator].colType() != 0;

    // Check if global recoil should be used.
    useLocalRecoilNow = !(globalRecoil && hardSystem
      && partonSystemsPtr->sizeOut(dip.system) <= nMaxGlobalRecoil);

    // Do not use global recoil if the radiator line has already branched.
    if (globalRecoilMode == 1 && isQCD) {
      if (globalRecoil && hardSystem) useLocalRecoilNow = true;
      for (int iHard = 0; iHard < int(hardPartons.size()); ++iHard)
        if ( event[dip.iRadiator].isAncestor(hardPartons[iHard]) )
          useLocalRecoilNow = false;
      // Check if global recoil should be used.
      if ( !globalRecoil || nGlobal >= nMaxGlobalBranch )
        useLocalRecoilNow = true;
    // Switch off global recoil after first trial emission.
    } else if (globalRecoilMode == 2 && isQCD) {
      useLocalRecoilNow = !(globalRecoil && hardSystem
        && nProposed.find(dip.system) != nProposed.end()
        && nProposed[dip.system]-infoPtr->getCounter(40) == 0);
      int nFinal = 0;
      for (int k = 0; k < int(event.size()); ++k)
        if ( event[k].isFinal() && event[k].colType() != 0) nFinal++;
      bool isFirst = (nHard == nFinal);

      // Switch off global recoil after first emission
      if ( globalRecoil && doInterleave && !isFirst )
        useLocalRecoilNow = true;
      // No global recoil for H-events.
      if ( nFinalBorn > 0 && nHard > nFinalBorn )
        useLocalRecoilNow = true;
    }

    // Dipole properties; normal local recoil.
    dip.mRad   = event[dip.iRadiator].m();
    if (useLocalRecoilNow) {
      dip.mRec = event[dip.iRecoiler].m();
      dip.mDip = m( event[dip.iRadiator], event[dip.iRecoiler] );

    // Dipole properties, alternative global recoil. Squares.
    } else {
      Vec4 pSumGlobal;
      // Include all particles in all hard systems (hard production system,
      // systems of resonance decay products) in the global recoil momentum.
      for (int iS = 0; iS < partonSystemsPtr->sizeSys(); ++iS) {
        for (int i = 0; i < partonSystemsPtr->sizeOut(iS); ++i) {
          int ii = partonSystemsPtr->getOut( iS, i);
          bool hasHardAncestor = event[ii].statusAbs() < 23;
          for (int iHard = 0; iHard < int(hardPartons.size()); ++iHard) {
            if ( event[ii].isAncestor(hardPartons[iHard])
              || ii == hardPartons[iHard]
              || (event[ii].status() == 23 && event[ii].colType() == 0) )
              hasHardAncestor = true;
          }
          if (hasHardAncestor && ii !=  dip.iRadiator && event[ii].isFinal() )
            pSumGlobal += event[ii].p();
        }
      }
      dip.mRec = pSumGlobal.mCalc();
      dip.mDip = m( event[dip.iRadiator].p(), pSumGlobal);
    }
    dip.m2Rad  = pow2(dip.mRad);
    dip.m2Rec  = pow2(dip.mRec);
    dip.m2Dip  = pow2(dip.mDip);

    // Find maximum evolution scale for dipole.
    dip.m2DipCorr    = pow2(dip.mDip - dip.mRec) - dip.m2Rad;
    double pTbegDip = min( pTbegAll, dip.pTmax );
    double pT2begDip = min( pow2(pTbegDip), 0.25 * dip.m2DipCorr);

    // For global recoil, always set the starting scale for first emission.
    bool isFirstWimpy = !useLocalRecoilNow && (pTmaxMatch == 1)
                      && nProposed.find(dip.system) != nProposed.end()
                      && (nProposed[dip.system] - infoPtr->getCounter(40) == 0
                      || isFirstTrial);
    double muQ        = (infoPtr->scalup() > 0.) ? infoPtr->scalup()
                      : infoPtr->QFac();
    if (isFirstWimpy && !limitMUQ) pT2begDip = pow2(muQ);
    else if (isFirstWimpy && limitMUQ) {
      // Find mass of colour dipole.
      double mS   = event[dip.iRecoiler].m();
      double mD   = m( event[dip.iRadiator], event[dip.iRecoiler] );
      double m2DC = pow2(mD - mS) - pow2(dip.mRad);
      // Choose minimal scale.
      pT2begDip = min( pow2(muQ), min(pow2(pTbegDip), 0.25 * m2DC) );
    }

    // Do not try splitting if the corrected dipole mass is negative.
    dip.pT2 = 0.;
    if (dip.m2DipCorr < 0.) {
      infoPtr->errorMsg("Warning in SimpleTimeShower::pTnext: "
      "negative dipole mass.");
      continue;
    }

    // Do QCD, QED, weak or HV evolution if it makes sense.
    if (pT2begDip > pT2sel) {
      if      (dip.colType != 0)
        pT2nextQCD(pT2begDip, pT2sel, dip, event);
      else if (dip.chgType != 0 || dip.gamType != 0)
        pT2nextQED(pT2begDip, pT2sel, dip, event);
      else if (dip.weakType != 0)
        pT2nextWeak(pT2begDip, pT2sel, dip, event);
      else if (dip.colvType != 0)
        pT2nextHV(pT2begDip, pT2sel, dip, event);

      // Update if found larger pT than current maximum.
      if (dip.pT2 > pT2sel) {
        pT2sel  = dip.pT2;
        dipSel  = &dip;
        iDipSel = iDip;
        splittingNameSel = splittingNameNow;
      }
    }
  }

  // Update the number of proposed timelike emissions.
  if (dipSel != 0 && nProposed.find(dipSel->system) != nProposed.end())
    ++nProposed[dipSel->system];

  // Return nonvanishing value if found pT bigger than already found.
  return (dipSel == 0) ? 0. : sqrt(pT2sel);

}

//--------------------------------------------------------------------------

// Evolve a QCD dipole end.

void SimpleTimeShower::pT2nextQCD(double pT2begDip, double pT2sel,
  TimeDipoleEnd& dip, Event& event) {

  // Lower cut for evolution. Return if no evolution range.
  double pT2endDip = max( pT2sel, pT2colCut );
  if (pT2begDip < pT2endDip) return;

  // For dipole recoil: no emission if the radiator is a quark,
  // since then a unified description is in SpaceShower.
  int    colTypeAbs = abs(dip.colType);
  if (doDipoleRecoil && dip.isrType != 0 && colTypeAbs == 1) return;

  // Upper estimate for matrix element weighting and colour factor.
  // Special cases for triplet recoiling against gluino and octet onia.
  // Note that g -> g g and g -> q qbar are split on two sides.
  double wtPSglue   = 2.;
  double colFac     = (colTypeAbs == 1) ? 4./3. : 3./2.;
  if (dip.MEgluinoRec)  colFac  = 3.;
  if (dip.isOctetOnium) colFac *= 0.5 * octetOniumColFac;
  // PS dec 2010. Include possibility for flexible normalization,
  // e.g., for dipoles stretched to junctions or to switch off radiation.
  if (dip.isFlexible)   colFac *= dip.flexFactor;
  double wtPSqqbar  = (colTypeAbs == 2)
    ? 0.25 * nGluonToQuark * extraGluonToQuark : 0.;

  // Variables used inside evolution loop. (Mainly dummy start values.)
  dip.pT2              = pT2begDip;
  int    nFlavour      = 3;
  double zMinAbs       = 0.5;
  double pT2min        = pT2endDip;
  double b0            = 4.5;
  double Lambda2       = Lambda3flav2;
  double emitCoefGlue  = 0.;
  double emitCoefQqbar = 0.;
  double emitCoefTot   = 0.;
  double wt            = 0.;
  bool   mustFindRange = true;

  // Add more headRoom if doing uncertainty variations
  // (to ensure at least a minimal number of failed branchings).
  doUncertaintiesNow   = doUncertainties;
  if (!uVarMPIshowers && dip.system != 0
    && partonSystemsPtr->hasInAB(dip.system)) doUncertaintiesNow = false;
  double overFac       = doUncertaintiesNow ? overFactor : 1.0;

  // Set default values for enhanced emissions.
  bool isEnhancedQ2QG, isEnhancedG2QQ, isEnhancedG2GG;
  isEnhancedQ2QG = isEnhancedG2QQ = isEnhancedG2GG = false;
  double enhanceNow = 1.;
  string nameNow = "";

  // Begin evolution loop towards smaller pT values.
  do {

    // Default values for current tentative emission.
    isEnhancedQ2QG = isEnhancedG2QQ = isEnhancedG2GG = false;
    enhanceNow = 1.;
    nameNow = "";

    // Initialize evolution coefficients at the beginning and
    // reinitialize when crossing c and b flavour thresholds.
    if (mustFindRange) {

      // Determine overestimated z range; switch at c and b masses.
      if (dip.pT2 > m2b) {
        nFlavour = 5;
        pT2min   = max( m2b, pT2endDip);
        b0       = 23./6.;
        Lambda2  = Lambda5flav2;
      } else if (dip.pT2 > m2c) {
        nFlavour = 4;
        pT2min   = max( m2c, pT2endDip);
        b0       = 25./6.;
        Lambda2  = Lambda4flav2;
      } else {
        nFlavour = 3;
        pT2min   = pT2endDip;
        b0       = 27./6.;
        Lambda2  = Lambda3flav2;
      }
      // A change of renormalization scale expressed by a change of Lambda.
      Lambda2 /= renormMultFac;

      // Calculate allowed z range; fail if it is too tiny.
      zMinAbs = 0.5 - sqrtpos( 0.25 - pT2min / dip.m2DipCorr );
      if (zMinAbs < SIMPLIFYROOT) zMinAbs = pT2min / dip.m2DipCorr;
      if (zMinAbs > 0.499) { dip.pT2 = 0.; return; }

      // Find emission coefficient for X -> X g.
      emitCoefGlue = overFac * wtPSglue * colFac * log(1. / zMinAbs - 1.);
      // Optionally enhanced branching rate.
      if (canEnhanceET && colTypeAbs == 2)
        emitCoefGlue *= userHooksPtr->enhanceFactor("fsr:G2GG");
      if (canEnhanceET && colTypeAbs == 1)
        emitCoefGlue *= userHooksPtr->enhanceFactor("fsr:Q2QG");

      // For dipole recoil: no g -> g g branching, since in SpaceShower.
      if (doDipoleRecoil && dip.isrType != 0 && colTypeAbs == 2)
        emitCoefGlue = 0.;

      // Find emission coefficient for g -> q qbar.
      emitCoefTot  = emitCoefGlue;
      if (colTypeAbs == 2 && event[dip.iRadiator].id() == 21) {
        emitCoefQqbar = overFac * wtPSqqbar * (1. - 2. * zMinAbs);
        // Optionally enhanced branching rate.
        if (canEnhanceET)
          emitCoefQqbar *= userHooksPtr->enhanceFactor("fsr:G2QQ");
        emitCoefTot  += emitCoefQqbar;
      }

      // Initialization done for current range.
      mustFindRange = false;
    }

    // Pick pT2 (in overestimated z range) for fixed alpha_strong.
    if (alphaSorder == 0) {
      dip.pT2 = dip.pT2 * pow( rndmPtr->flat(),
        1. / (alphaS2pi * emitCoefTot) );

    // Ditto for first-order alpha_strong.
    } else if (alphaSorder == 1) {
      dip.pT2 = Lambda2 * pow( dip.pT2 / Lambda2,
        pow( rndmPtr->flat(), b0 / emitCoefTot) );

      // For second order reject by second term in alpha_strong expression.
    } else {
      do dip.pT2 = Lambda2 * pow( dip.pT2 / Lambda2,
        pow( rndmPtr->flat(), b0 / emitCoefTot) );
      while (alphaS.alphaS2OrdCorr(renormMultFac * dip.pT2) < rndmPtr->flat()
        && dip.pT2 > pT2min);
    }
    wt = 0.;

    // If crossed c or b thresholds: continue evolution from threshold.
    if (nFlavour == 5 && dip.pT2 < m2b) {
      mustFindRange = true;
      dip.pT2       = m2b;
    } else if ( nFlavour == 4 && dip.pT2 < m2c) {
      mustFindRange = true;
      dip.pT2       = m2c;

    // Abort evolution if below cutoff scale, or below another branching.
    } else {
      if ( dip.pT2 < pT2endDip) { dip.pT2 = 0.; return; }

      // Pick kind of branching: X -> X g or g -> q qbar.
      dip.flavour  = 21;
      dip.mFlavour = 0.;
      if (colTypeAbs == 2 && emitCoefQqbar > rndmPtr->flat()
        * emitCoefTot) dip.flavour = 0;

      // Pick z: either dz/(1-z) or flat dz.
      if (dip.flavour == 21) {
        dip.z = 1. - zMinAbs * pow( 1. / zMinAbs - 1., rndmPtr->flat() );
      } else {
        dip.z = zMinAbs + (1. - 2. * zMinAbs) * rndmPtr->flat();
      }

      // Do not accept branching if outside allowed z range.
      double zMin = 0.5 - sqrtpos( 0.25 - dip.pT2 / dip.m2DipCorr );
      if (zMin < SIMPLIFYROOT) zMin = dip.pT2 / dip.m2DipCorr;
      dip.m2 = dip.m2Rad + dip.pT2 / (dip.z * (1. - dip.z));
      if (dip.z > zMin && dip.z < 1. - zMin
        && dip.m2 * dip.m2Dip < dip.z * (1. - dip.z)
        * pow2(dip.m2Dip + dip.m2 - dip.m2Rec) ) {

        // Flavour choice for g -> q qbar.
        if (dip.flavour == 0) {
          dip.flavour  = min(5, 1 + int(nGluonToQuark * rndmPtr->flat()));
          dip.mFlavour = particleDataPtr->m0(dip.flavour);
        }

        if (dip.flavour == 21
          && (colTypeAbs == 1 || colTypeAbs == 3) ) {
          nameNow = "fsr:Q2QG";
          // Optionally enhanced branching rate.
          if (canEnhanceET) {
            double enhance = userHooksPtr->enhanceFactor(nameNow);
            if (enhance != 1.) {
              enhanceNow = enhance;
              isEnhancedQ2QG = true;
            }
          }
        } else if (dip.flavour == 21) {
          nameNow = "fsr:G2GG";
          // Optionally enhanced branching rate.
          if (canEnhanceET) {
            double enhance = userHooksPtr->enhanceFactor(nameNow);
            if (enhance != 1.) {
              enhanceNow = enhance;
              isEnhancedG2GG = true;
            }
          }
        } else {
          if      (dip.flavour <  4) nameNow = "fsr:G2QQ";
          else if (dip.flavour == 4) nameNow = "fsr:G2QQ:cc";
          else                       nameNow = "fsr:G2QQ:bb";
          // Optionally enhanced branching rate.
          if (canEnhanceET) {
            double enhance = userHooksPtr->enhanceFactor(nameNow);
            if (enhance != 1.) {
              enhanceNow = enhance;
              isEnhancedG2QQ = true;
            }
          }
        }

        // No z weight, except threshold, if to do ME corrections later on.
        if (dip.MEtype > 0) {
          wt = 1.;
          if (dip.flavour < 10 && dip.m2 < THRESHM2 * pow2(dip.mFlavour))
            wt = 0.;

        // z weight for X -> X g.
        } else if (dip.flavour == 21
          && (colTypeAbs == 1 || colTypeAbs == 3) ) {
          wt = (1. + pow2(dip.z)) / wtPSglue;

        // z weight for g -> g g; optional suppression for massive recoiler.
        } else if (dip.flavour == 21) {
          wt = (1. + pow3(dip.z)) / wtPSglue;
          if (recoilDeadCone && dip.mRec > 0.) {
            double r2G = dip.m2Rec / dip.m2Dip;
            double x1G = (1. - r2G + dip.m2 / dip.m2Dip) * dip.z;
            double x2G =  1. + r2G - dip.m2 / dip.m2Dip;
            wt *= 1. - (r2G / max(XMARGIN, x1G + x2G - 1. - r2G))
              * (max(XMARGIN, 1. + r2G - x2G) / max(XMARGIN,1. - r2G - x1G));
          }

        // z weight for g -> q qbar: different options.
        } else {
          double ratioQ = pow2(dip.mFlavour) / dip.m2;
          double betaQ  = sqrtpos( 1. - 4. * ratioQ );
          if (weightGluonToQuark%4 == 1) {
            wt = betaQ * ( pow2(dip.z) + pow2(1. - dip.z) );
          } else if (weightGluonToQuark%4 == 2) {
            wt = betaQ * ( pow2(dip.z) + pow2(1. - dip.z)
               + 8. * ratioQ * dip.z * (1. - dip.z) );
          } else {
            double m2Rat = dip.m2 / dip.m2DipCorr;
            double zCosThe = ((1. + m2Rat) * dip.z - m2Rat) / (1. - m2Rat);
            wt = betaQ * ( pow2(zCosThe) + pow2(1. - zCosThe)
               + 8. * ratioQ * zCosThe * (1. - zCosThe) )
               * (1. + m2Rat) / ((1. - m2Rat) * extraGluonToQuark) ;
            if (weightGluonToQuark%4 == 0) wt *= pow3(1. - m2Rat);
          }
          if (weightGluonToQuark > 4 && alphaSorder > 0)
            wt *= log(dip.pT2 / Lambda2)
                / log(scaleGluonToQuark * dip.m2 / Lambda2);
        }

        // Cancel out extra uncertainty-band headroom factors.
        wt /= overFac;

        // Suppression factors for dipole to beam remnant.
        if (dip.isrType != 0 && useLocalRecoilNow) {
          BeamParticle& beam = (dip.isrType == 1) ? *beamAPtr : *beamBPtr;
          int iSysRec = dip.systemRec;
          double xOld = beam[iSysRec].x();
          double xNew = xOld * (1. + (dip.m2 - dip.m2Rad) /
            (dip.m2Dip - dip.m2Rad));
          double xMaxAbs = beam.xMax(iSysRec);
          if (xMaxAbs < 0.) {
            infoPtr->errorMsg("Warning in SimpleTimeShower::pT2nextQCD: "
            "xMaxAbs negative");
            return;
          }

          // New: Ensure that no x-value larger than unity is picked. Only
          // necessary for imprecise LHE input.
          if (xNew > 1.) wt = 0.;

          // Firstly reduce by PDF ratio.
          if (xNew > xMaxAbs) wt = 0.;
          else {
            int idRec     = event[dip.iRecoiler].id();
            pdfScale2 = (useFixedFacScale) ? fixedFacScale2
              : factorMultFac * dip.pT2;
            double pdfOld = max ( TINYPDF,
              beam.xfISR( iSysRec, idRec, xOld, pdfScale2) );
            double pdfNew =
              beam.xfISR( iSysRec, idRec, xNew, pdfScale2);
            wt *= min( 1., pdfNew / pdfOld);
          }

          // Secondly optionally reduce by 4 pT2_hard / (4 pT2_hard + m2).
          if (dampenBeamRecoil) {
            double pTpT = sqrt(event[dip.iRadiator].pT2() * dip.pT2);
            wt *= pTpT / (pTpT + dip.m2);
          }
        }

        // Optional dampening of large pT values in hard system.
        if (dopTdamp && dip.system == 0 && dip.MEtype == 0)
          wt *= pT2damp / (dip.pT2 + pT2damp);
      }
    }

    // If doing uncertainty variations, postpone accept/reject to branch().
    if (wt > 0. && dip.pT2 > pT2min && doUncertaintiesNow) {
      dip.pAccept = wt;
      wt          = 1.0;
    }

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < rndmPtr->flat());

  // Store outcome of enhanced branching rate analysis.
  splittingNameNow = nameNow;
  if (canEnhanceET) {
    if (isEnhancedQ2QG) storeEnhanceFactor(dip.pT2,"fsr:Q2QG", enhanceNow);
    if (isEnhancedG2QQ) storeEnhanceFactor(dip.pT2,"fsr:G2QQ", enhanceNow);
    if (isEnhancedG2GG) storeEnhanceFactor(dip.pT2,"fsr:G2GG", enhanceNow);
  }

}

//--------------------------------------------------------------------------

// Evolve a QED dipole end, either charged or photon.

void SimpleTimeShower::pT2nextQED(double pT2begDip, double pT2sel,
  TimeDipoleEnd& dip, Event& event) {

  // Lower cut for evolution. Return if no evolution range.
  double pT2chgCut = (dip.chgType != 0 && abs(dip.chgType) != 3)
    ? pT2chgQCut : pT2chgLCut;
  double pT2endDip = max( pT2sel, pT2chgCut );
  if (pT2begDip < pT2endDip) return;

  // Emission of photon or photon branching.
  bool hasCharge = (dip.chgType != 0);

  // Default values.
  double wtPSgam     = 0.;
  double chg2Sum     = 0.;
  double chg2SumL    = 0.;
  double chg2SumQ    = 0.;
  double zMinAbs     = 0.;
  double emitCoefTot = 0.;

  // alpha_em at maximum scale provides upper estimate.
  double alphaEMmax  = alphaEM.alphaEM(renormMultFac * dip.m2DipCorr);
  double alphaEM2pi  = alphaEMmax / (2. * M_PI);

  // Set default values for enhanced emissions.
  bool isEnhancedQ2QA, isEnhancedA2LL, isEnhancedA2QQ;
  isEnhancedQ2QA = isEnhancedA2LL = isEnhancedA2QQ = false;
  double enhanceNow = 1.;
  string nameNow = "";

  // Emission: upper estimate for matrix element weighting; charge factor.
  if (hasCharge) {
    wtPSgam     = 2.;
    double chg2 = pow2(dip.chgType / 3.);

    // Determine overestimated z range. Find evolution coefficient.
    zMinAbs = 0.5 - sqrtpos( 0.25 - pT2endDip / dip.m2DipCorr );
    if (zMinAbs < SIMPLIFYROOT) zMinAbs = pT2endDip / dip.m2DipCorr;
    emitCoefTot = alphaEM2pi * chg2 * wtPSgam * log(1. / zMinAbs - 1.);
    // Optionally enhanced branching rate.
    if (canEnhanceET) emitCoefTot *= userHooksPtr->enhanceFactor("fsr:Q2QA");

  // Branching: sum of squared charge factors for lepton and quark daughters.
  } else {
    chg2SumL = max(0, min(3, nGammaToLepton));
    if      (nGammaToQuark > 4) chg2SumQ = 11. / 9.;
    else if (nGammaToQuark > 3) chg2SumQ = 10. / 9.;
    else if (nGammaToQuark > 2) chg2SumQ =  6. / 9.;
    else if (nGammaToQuark > 1) chg2SumQ =  5. / 9.;
    else if (nGammaToQuark > 0) chg2SumQ =  1. / 9.;

    // Optionally enhanced branching rate.
    if (canEnhanceET) chg2SumL *= userHooksPtr->enhanceFactor("fsr:A2LL");
    if (canEnhanceET) chg2SumQ *= userHooksPtr->enhanceFactor("fsr:A2QQ");

    // Total sum of squared charge factors. Find evolution coefficient.
    chg2Sum     = chg2SumL + 3. * chg2SumQ;
    emitCoefTot = alphaEM2pi * chg2Sum * extraGluonToQuark;
  }

  // Variables used inside evolution loop.
  dip.pT2 = pT2begDip;
  double wt;

  // Begin evolution loop towards smaller pT values.
  do {


    // Default values for current tentative emission.
    isEnhancedQ2QA = isEnhancedA2LL = isEnhancedA2QQ = false;
    enhanceNow = 1.;
    nameNow = "";

    // Pick pT2 (in overestimated z range).
    dip.pT2 = dip.pT2 * pow(rndmPtr->flat(), 1. / emitCoefTot);
    wt = 0.;

    // Abort evolution if below cutoff scale, or below another branching.
    if ( dip.pT2 < pT2endDip) { dip.pT2 = 0.; return; }

    // Pick z according to dz/(1-z) or flat.
    if (hasCharge) dip.z = 1. - zMinAbs
      * pow( 1. / zMinAbs - 1., rndmPtr->flat() );
    else           dip.z = rndmPtr->flat();

    // Do not accept branching if outside allowed z range.
    double zMin = 0.5 - sqrtpos( 0.25 - dip.pT2 / dip.m2DipCorr );
    if (zMin < SIMPLIFYROOT) zMin = dip.pT2 / dip.m2DipCorr;
    if (dip.z <= zMin || dip.z >= 1. - zMin) continue;
    dip.m2 = dip.m2Rad + dip.pT2 / (dip.z * (1. - dip.z));
    if (dip.m2 * dip.m2Dip < dip.z * (1. - dip.z)
        * pow2(dip.m2Dip + dip.m2 - dip.m2Rec)
      // For gamma -> f fbar also impose maximum mass.
      && (hasCharge || dip.m2 < m2MaxGamma) ) {

      // Photon emission: unique flavour choice.
      if (hasCharge) {
        dip.flavour = 22;
        dip.mFlavour = 0.;

      // Photon branching: either lepton or quark flavour choice.
      } else {
        if (rndmPtr->flat() * chg2Sum < chg2SumL)
          dip.flavour  = 9 + 2 * min(3, 1 + int(chg2SumL * rndmPtr->flat()));
        else {
          double rndmQ = 9. * chg2SumQ * rndmPtr->flat();
          if      (rndmQ <  1.) dip.flavour = 1;
          else if (rndmQ <  5.) dip.flavour = 2;
          else if (rndmQ <  6.) dip.flavour = 3;
          else if (rndmQ < 10.) dip.flavour = 4;
          else                  dip.flavour = 5;
        }
        dip.mFlavour = particleDataPtr->m0(dip.flavour);
      }


      if (hasCharge) {
        nameNow = "fsr:Q2QA";
        // Optionally enhanced branching rate.
        if (canEnhanceET) {
          double enhance = userHooksPtr->enhanceFactor(nameNow);
          if (enhance != 1.) {
            enhanceNow = enhance;
            isEnhancedQ2QA = true;
          }
        }
      } else if (dip.flavour > 10) {
        nameNow = "fsr:A2LL";
        // Optionally enhanced branching rate.
        if (canEnhanceET) {
          double enhance = userHooksPtr->enhanceFactor(nameNow);
          if (enhance != 1.) {
            enhanceNow = enhance;
            isEnhancedA2LL = true;
          }
        }
      } else {
        nameNow = "fsr:A2QQ";
        // Optionally enhanced branching rate.
        if (canEnhanceET) {
          double enhance = userHooksPtr->enhanceFactor(nameNow);
          if (enhance != 1.) {
            enhanceNow = enhance;
            isEnhancedA2QQ = true;
          }
        }
      }

      // No z weight, except threshold, if to do ME corrections later on.
      if (dip.MEtype > 0) {
        wt = 1.;
        if (dip.flavour < 20 && dip.m2 < THRESHM2 * pow2(dip.mFlavour))
          wt = 0.;

      // z weight for X -> X gamma.
      } else if (hasCharge) {
        wt = (1. + pow2(dip.z)) / wtPSgam;

      // z weight for gamma -> f fbar; different options.
      } else {
        double ratioF = pow2(dip.mFlavour) / dip.m2;
        double betaF  = sqrtpos( 1. - 4. * ratioF );
        if (weightGluonToQuark%4 == 1) {
          wt = betaF * ( pow2(dip.z) + pow2(1. - dip.z) );
        } else if (weightGluonToQuark%4 == 2) {
          wt = betaF * ( pow2(dip.z) + pow2(1. - dip.z)
             + 8. * ratioF * dip.z * (1. - dip.z) );
        } else {
          double m2Rat = dip.m2 / dip.m2DipCorr;
          double zCosThe = ((1. + m2Rat) * dip.z - m2Rat) / (1. - m2Rat);
          wt = betaF * ( pow2(zCosThe) + pow2(1. - zCosThe)
             + 8. * ratioF * zCosThe * (1. - zCosThe) )
             * (1. + m2Rat) / ((1. - m2Rat) * extraGluonToQuark) ;
          if (weightGluonToQuark%4 == 0) wt *= pow3(1. - m2Rat);
        }
      }

      // Correct to current value of alpha_EM.
      double aEMscale = dip.pT2;
      if (dip.flavour < 20 && weightGluonToQuark > 4)
        aEMscale = scaleGluonToQuark * dip.m2;
      double alphaEMnow = alphaEM.alphaEM(renormMultFac * aEMscale);
      wt *= (alphaEMnow / alphaEMmax);

      // Suppression factors for dipole to beam remnant.
      if (dip.isrType != 0 && useLocalRecoilNow) {
        BeamParticle& beam = (dip.isrType == 1) ? *beamAPtr : *beamBPtr;
        int iSys    = dip.system;
        double xOld = beam[iSys].x();
        double xNew = xOld * (1. + (dip.m2 - dip.m2Rad) /
          (dip.m2Dip - dip.m2Rad));
        double xMaxAbs = beam.xMax(iSys);
        if (xMaxAbs < 0.) {
          infoPtr->errorMsg("Warning in SimpleTimeShower::pT2nextQED: "
          "xMaxAbs negative");
          return;
        }

        // New: Ensure that no x-value larger than unity is picked. Only
        // necessary for imprecise LHE input.
        if (xNew > 1.) wt = 0.;

        // Firstly reduce by PDF ratio.
        if (xNew > xMaxAbs) wt = 0.;
        else {
          int idRec     = event[dip.iRecoiler].id();
          pdfScale2 = (useFixedFacScale) ? fixedFacScale2
            : factorMultFac * dip.pT2;
          double pdfOld = max ( TINYPDF,
            beam.xfISR( iSys, idRec, xOld, pdfScale2) );
          double pdfNew =
            beam.xfISR( iSys, idRec, xNew, pdfScale2);
          wt *= min( 1., pdfNew / pdfOld);
        }

        // Secondly optionally reduce by 4 pT2_hard / (4 pT2_hard + m2).
        if (dampenBeamRecoil) {
          double pT24 = 4. * event[dip.iRadiator].pT2();
          wt *= pT24 / (pT24 + dip.m2);
        }
      }

      // Optional dampening of large pT values in hard system.
      if (dopTdamp && dip.system == 0 && dip.MEtype == 0)
        wt *= pT2damp / (dip.pT2 + pT2damp);
    }

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < rndmPtr->flat());

  // Store outcome of enhanced branching rate analysis.
  splittingNameNow = nameNow;
  if (canEnhanceET) {
    if (isEnhancedQ2QA) storeEnhanceFactor(dip.pT2,"fsr:Q2QA", enhanceNow);
    if (isEnhancedA2LL) storeEnhanceFactor(dip.pT2,"fsr:A2LL", enhanceNow);
    if (isEnhancedA2QQ) storeEnhanceFactor(dip.pT2,"fsr:A2QQ", enhanceNow);
  }

}

//--------------------------------------------------------------------------

// Evolve a weak-emission dipole end.

void SimpleTimeShower::pT2nextWeak(double pT2begDip, double pT2sel,
  TimeDipoleEnd& dip, Event& event) {

  // Lower cut for evolution. Return if no evolution range.
  double pT2endDip = max( pT2sel, pT2weakCut );
  if (pT2begDip < pT2endDip) return;

  // Default values.
  double wtPSgam     = 0.;
  double zMinAbs     = 0.;
  double emitCoefTot = 0.;

  // alpha_em at maximum scale provides upper estimate.
  double alphaEMmax  = alphaEM.alphaEM(renormMultFac * pT2begDip);
  double alphaEM2pi  = alphaEMmax / (2. * M_PI);

  // Emission: upper estimate for matrix element weighting; charge factor.
  wtPSgam     = 8.;

  // Determine overestimated z range. Find evolution coefficient.
  zMinAbs = 0.5 - sqrtpos( 0.25 - pT2endDip / dip.m2DipCorr );
  if (zMinAbs < SIMPLIFYROOT) zMinAbs = pT2endDip / dip.m2DipCorr;

  // Determine weak coupling.
  double weakCoupling = 0.;
  // W-radiation, with additional factor of two, from only having
  // left-handed fermions.
  if (dip.weakType == 1)
    weakCoupling = 2. * alphaEM2pi / (4. * coupSMPtr->sin2thetaW());
  // Z-radiation, split between left and right fermions.
  else if (dip.weakType == 2 && dip.weakPol == -1)
    weakCoupling = alphaEM2pi * thetaWRat
      * pow2(2. * coupSMPtr->lf( event[dip.iRadiator].idAbs() ));
  else
    weakCoupling = alphaEM2pi * thetaWRat
      * pow2(2. * coupSMPtr->rf( event[dip.iRadiator].idAbs() ));

  // Set default values for enhanced emissions.
  bool isEnhancedQ2QW;
  isEnhancedQ2QW = false;
  double enhanceNow = 1.;
  string nameNow = "";

  // Variables used inside evolution loop.
  emitCoefTot = weakEnhancement * weakCoupling
    * wtPSgam * log(1. / zMinAbs - 1.);
  // Fudge factor to correct for weak PS not being an overestimate of ME.
  if ( dip.MEtype == 201 || dip.MEtype == 202 || dip.MEtype == 203
    || dip.MEtype == 206 || dip.MEtype == 207 || dip.MEtype == 208 )
    emitCoefTot *= WEAKPSWEIGHT;
  dip.pT2 = pT2begDip;
  double wt;

  // Optionally enhanced branching rate.
  if (canEnhanceET) emitCoefTot *= userHooksPtr->enhanceFactor("fsr:Q2QW");

  // Begin evolution loop towards smaller pT values.
  do {

    // Default values for current tentative emission.
    isEnhancedQ2QW = false;
    enhanceNow = 1.;
    nameNow = "";

    // Pick pT2 (in overestimated z range).
    dip.pT2 = dip.pT2 * pow(rndmPtr->flat(), 1. / emitCoefTot);
    wt = 0.;

    // Abort evolution if below cutoff scale, or below another branching.
    if ( dip.pT2 < pT2endDip) {dip.pT2 = 0.; return; }

    // Pick z according to dz/(1-z) or flat.
    dip.z = 1. - zMinAbs * pow( 1. / zMinAbs - 1., rndmPtr->flat() );

    // Do not accept branching if outside allowed z range.
    double zMin = 0.5 - sqrtpos( 0.25 - dip.pT2 / dip.m2DipCorr );
    if (zMin < SIMPLIFYROOT) zMin = dip.pT2 / dip.m2DipCorr;
    dip.m2 = dip.m2Rad + dip.pT2 / (dip.z * (1. - dip.z));
    if (dip.z > zMin && dip.z < 1. - zMin
      && dip.m2 * dip.m2Dip < dip.z * (1. - dip.z)
        * pow2(dip.m2Dip + dip.m2 - dip.m2Rec) ) {

      // Check whether emission of W+ or W-, or else Z0.
      if (dip.weakType == 1) {
        dip.flavour = (event[dip.iRadiator].id() > 0) ? 24 : -24;
        if (event[dip.iRadiator].idAbs() % 2 == 1) dip.flavour = -dip.flavour;
      } else if (dip.weakType == 2) dip.flavour = 23;

      // Set mass of emitted particle, with Breit-Wigner distribution.
      dip.mFlavour = particleDataPtr->mSel( dip.flavour);

      // No z weight, except threshold, if to do ME corrections later on.
      // Here no pure shower mode exists, always needs ME corrections.
      if (dip.MEtype > 0) wt = 1.;

      // Correct to current value of alpha_EM.
      double alphaEMnow = alphaEM.alphaEM(renormMultFac * dip.pT2);
      wt *= (alphaEMnow / alphaEMmax);

      nameNow = "fsr:Q2QW";
      // Optionally enhanced branching rate.
      if (canEnhanceET) {
        double enhance = userHooksPtr->enhanceFactor(nameNow);
        if (enhance != 1.) {
          enhanceNow = enhance;
          isEnhancedQ2QW = true;
        }
      }

      // Suppression factors for dipole to beam remnant.
      if (dip.isrType != 0 && useLocalRecoilNow) {
        BeamParticle& beam = (dip.isrType == 1) ? *beamAPtr : *beamBPtr;
        int iSys    = dip.system;
        double xOld = beam[iSys].x();
        double xNew = xOld * (1. + (dip.m2 - dip.m2Rad) /
          (dip.m2Dip - dip.m2Rad));
        double xMaxAbs = beam.xMax(iSys);
        if (xMaxAbs < 0.) {
          infoPtr->errorMsg("Warning in SimpleTimeShower::pT2nextWeak: "
          "xMaxAbs negative");
          return;
        }

        // New: Ensure that no x-value larger than unity is picked. Only
        // necessary for imprecise LHE input.
        if (xNew > 1.) wt = 0.;

        // Firstly reduce by PDF ratio.
        if (xNew > xMaxAbs) wt = 0.;
        else {
          int idRec     = event[dip.iRecoiler].id();
          pdfScale2 = (useFixedFacScale) ? fixedFacScale2
            : factorMultFac * dip.pT2;
          double pdfOld = max ( TINYPDF,
            beam.xfISR( iSys, idRec, xOld, pdfScale2) );
          double pdfNew =
            beam.xfISR( iSys, idRec, xNew, pdfScale2);
          wt *= min( 1., pdfNew / pdfOld);
        }
        // Secondly optionally reduce by 4 pT2_hard / (4 pT2_hard + m2).
        if (dampenBeamRecoil) {
          double pT24 = 4. * event[dip.iRadiator].pT2();
          wt *= pT24 / (pT24 + dip.m2);
        }
      }
    }

    // Optional dampening of large pT values in hard system.
    if (dopTdamp && dip.system == 0) wt *= pT2damp / (dip.pT2 + pT2damp);

    // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < rndmPtr->flat());

  // Store outcome of enhanced branching rate analysis.
  splittingNameNow = nameNow;
  if (canEnhanceET && isEnhancedQ2QW)
    storeEnhanceFactor(dip.pT2,"fsr:Q2QW", enhanceNow);

}

//--------------------------------------------------------------------------

// Evolve a Hidden Valley dipole end.

void SimpleTimeShower::pT2nextHV(double pT2begDip, double pT2sel,
  TimeDipoleEnd& dip, Event& ) {

  // Lower cut for evolution. Return if no evolution range.
  double pT2endDip = max( pT2sel, pT2hvCut );
  if (pT2begDip < pT2endDip) return;

  // C_F * alpha_HV/2 pi.
  int    colvTypeAbs = abs(dip.colvType);
  double colvFac     = (colvTypeAbs == 1) ? CFHV : 0.5 * nCHV;
  double alphaHV2pi  = alphaHVfix / (2. * M_PI);
  double b0HV        = (11. /6. * nCHV - 2. / 6. * nFlavHV);

  // Determine overestimated z range. Find evolution coefficient.
  double zMinAbs = 0.5 - sqrtpos( 0.25 - pT2endDip / dip.m2DipCorr );
  if (zMinAbs < SIMPLIFYROOT) zMinAbs = pT2endDip / dip.m2DipCorr;
  double emitCoefTot = colvFac * 2. * log(1. / zMinAbs - 1.);
  double LambdaHV2 = pow2(LambdaHV);

  // Variables used inside evolution loop.
  dip.pT2 = pT2begDip;
  double wt;

  // Set default values for enhanced emissions.
  bool isEnhancedQ2QHV;
  isEnhancedQ2QHV = false;
  double enhanceNow = 1.;
  string nameNow = "";

  // Optionally enhanced branching rate.
  if (canEnhanceET) emitCoefTot *= userHooksPtr->enhanceFactor("fsr:Q2QHV");

  // Begin evolution loop towards smaller pT values.
  do {

    // Default values for current tentative emission.
    isEnhancedQ2QHV = false;
    enhanceNow = 1.;
    nameNow = "";

    // Pick pT2 (in overestimated z range), fixed or first-order alpha_strong.
    if (alphaHVorder == 0) {
      dip.pT2 = dip.pT2 * pow( rndmPtr->flat(),
        1. / (alphaHV2pi * emitCoefTot) );
    } else if (alphaHVorder == 1) {
      dip.pT2 = LambdaHV2 * pow( dip.pT2 / LambdaHV2,
        pow( rndmPtr->flat(), b0HV / emitCoefTot) );
    }
    wt = 0.;

    // Abort evolution if below cutoff scale, or below another branching.
    if ( dip.pT2 < pT2endDip) { dip.pT2 = 0.; return; }

    // Pick z according to dz/(1-z).
    dip.z = 1. - zMinAbs * pow( 1. / zMinAbs - 1., rndmPtr->flat() );

    // Do not accept branching if outside allowed z range.
    double zMin = 0.5 - sqrtpos( 0.25 - dip.pT2 / dip.m2DipCorr );
    if (zMin < SIMPLIFYROOT) zMin = dip.pT2 / dip.m2DipCorr;
    dip.m2 = dip.m2Rad + dip.pT2 / (dip.z * (1. - dip.z));
    if (dip.z > zMin && dip.z < 1. - zMin
      && dip.m2 * dip.m2Dip < dip.z * (1. - dip.z)
        * pow2(dip.m2Dip + dip.m2 - dip.m2Rec) ) {

      // HV gamma or gluon emission: unique flavour choice.
      dip.flavour  = idHV;
      dip.mFlavour = mHV;

      // No z weight, except threshold, if to do ME corrections later on.
      if (dip.MEtype > 0) wt = 1.;

      // z weight for X -> X g_HV.
      else if (colvTypeAbs == 1) wt = (1. + pow2(dip.z)) / 2.;
      else wt = (1. + pow3(dip.z)) / 2.;

      nameNow = "fsr:Q2QHV";
      // Optionally enhanced branching rate.
      if (canEnhanceET) {
        double enhance = userHooksPtr->enhanceFactor(nameNow);
        if (enhance != 1.) {
          enhanceNow = enhance;
          isEnhancedQ2QHV = true;
        }
      }

    }

    // Optional dampening of large pT values in hard system.
    if (dopTdamp && dip.system == 0 && dip.MEtype == 0)
      wt *= pT2damp / (dip.pT2 + pT2damp);

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < rndmPtr->flat());

  // Store outcome of enhanced branching rate analysis.
  splittingNameNow = nameNow;
  if (canEnhanceET && isEnhancedQ2QHV)
    storeEnhanceFactor(dip.pT2,"fsr:Q2QHV", enhanceNow);

}

//--------------------------------------------------------------------------

// ME corrections and kinematics that may give failure.
// Notation: radBef, recBef = radiator, recoiler before emission,
//           rad, rec, emt = radiator, recoiler, emitted efter emission.
//           (rad, emt distinguished by colour flow for g -> q qbar.)

bool SimpleTimeShower::branch( Event& event, bool isInterleaved) {

  // Check if this system is part of the hard scattering
  // (including resonance decay products).
  bool hardSystem = partonSystemsPtr->getHard(dipSel->system);
  bool isQCD = event[dipSel->iRadiator].colType() != 0;

  // Check if global recoil should be used in resonance showers.
  useLocalRecoilNow = !(globalRecoil && hardSystem
    && partonSystemsPtr->sizeOut(dipSel->system) <= nMaxGlobalRecoil);

  // Do not use global recoil if the radiator line has already branched.
  if (globalRecoilMode == 1 && isQCD) {
    if ( globalRecoil && hardSystem) useLocalRecoilNow = true;
    for (int iHard = 0; iHard < int(hardPartons.size()); ++iHard)
      if ( event[dipSel->iRadiator].isAncestor(hardPartons[iHard]) )
        useLocalRecoilNow = false;
    // Check if global recoil should be used.
    if ( !globalRecoil || nGlobal >= nMaxGlobalBranch )
      useLocalRecoilNow = true;

  // Switch off global recoil after first trial emission
  } else if (globalRecoilMode == 2 && isQCD) {
    useLocalRecoilNow = !(globalRecoil
      && nProposed.find(dipSel->system) != nProposed.end()
      && nProposed[dipSel->system] - infoPtr->getCounter(40) == 1);
    // Check if global recoil should be used.
    int nFinal = 0;
    for (int i = 0; i < int(event.size()); ++i)
      if ( event[i].isFinal() && event[i].colType() != 0) nFinal++;
    bool isFirst = (nHard == nFinal);
    if ( globalRecoil && doInterleave && !isFirst )
      useLocalRecoilNow = true;
    // No global recoil for H-events.
    if ( nFinalBorn > 0 && nHard > nFinalBorn )
      useLocalRecoilNow = true;
  }

  // Check if the first emission should be studied for removal.
  bool canMergeFirst = (mergingHooksPtr != 0)
                     ? mergingHooksPtr->canVetoEmission() : false;

  int npartons = 0, nfinal = 0, nw = 0, nz = 0;
  for ( int i = 0; i < event.size(); ++i) {if(event[i].isFinal() ) {
      nfinal++;
      if (event[i].colType() != 0) npartons++;
      if (event[i].id() == 23) nz++;
      if (event[i].idAbs() == 24) nw++;
    }
  }

  // Find initial radiator and recoiler particles in dipole branching.
  int iRadBef      = dipSel->iRadiator;
  int iRecBef      = dipSel->iRecoiler;
  Particle& radBef = event[iRadBef];
  Particle& recBef = event[iRecBef];

  // Find their momenta, with special sum for global recoil.
  Vec4 pRadBef     = event[iRadBef].p();
  Vec4 pRecBef;
  vector<int> iGRecBef, iGRec;
  if (useLocalRecoilNow) pRecBef =  event[iRecBef].p();
  else {
    // Include all particles in all hard systems (hard production system,
    // systems of resonance decay products) in the global recoil momentum.
    for (int iS = 0; iS < partonSystemsPtr->sizeSys(); ++iS) {
      for (int i = 0; i < partonSystemsPtr->sizeOut(iS); ++i) {
        int iG = partonSystemsPtr->getOut( iS, i);
        bool hasHardAncestor = event[iG].statusAbs() < 23;
        for (int iHard = 0; iHard < int(hardPartons.size()); ++iHard)
          if ( event[iG].isAncestor(hardPartons[iHard])
            || iG == hardPartons[iHard]
            || (event[iG].status() == 23 && event[iG].colType() == 0))
            hasHardAncestor = true;
        if (hasHardAncestor && iG != dipSel->iRadiator
          && event[iG].isFinal() ) {
          iGRecBef.push_back(iG);
          pRecBef += event[iG].p();
        }
      }
    }
  }

  // Find old incoming momenta for weak shower t-channel ME correction.
  Vec4 p3weak, p4weak;
  if (dipSel->MEtype >= 200 && dipSel->MEtype <= 210) {
    p3weak = event[3].p();
    p4weak = event[4].p();
  }
  if ( dipSel->MEtype == 201 || dipSel->MEtype == 202
    || dipSel->MEtype == 203 || dipSel->MEtype == 206
    || dipSel->MEtype == 207 || dipSel->MEtype == 208) {
    if (!weakExternal) {
      // Trace back to original mother. MPI not allowed to radiate weakly.
      int i2to2Mother = iRadBef;
      while (i2to2Mother != 5 && i2to2Mother != 6 && i2to2Mother != 0)
        i2to2Mother = event[i2to2Mother].mother1();
      if (i2to2Mother == 0) return false;

      // u d -> u d  && u g -> u g.
      if (event[3].id() != event[4].id()) {
        if (event[3].id() == event[i2to2Mother].id());
        else if (event[4].id() == event[i2to2Mother].id())
          swap(p3weak, p4weak);
        // In case of no match, assign random combination.
        else if (rndmPtr->flat() < 0.5) swap(p3weak, p4weak);
      }
      // u u -> u u, assign random combination.
      else if (rndmPtr->flat() < 0.5) swap(p3weak, p4weak);
    } else {
      int i2to2Mother = iRadBef;
      while (i2to2Mother >= weakHardSize)
        i2to2Mother = event[i2to2Mother].mother1();
      if (weak2to2lines[2] == i2to2Mother) {
        p3weak = weakMomenta[0];
        p4weak = weakMomenta[1];
      } else {
        p3weak = weakMomenta[1];
        p4weak = weakMomenta[0];
      }
    }
  }

  // Default flavours and colour tags for new particles in dipole branching.
  int idRad        = radBef.id();
  int idEmt        = dipSel->flavour;
  int colRad       = radBef.col();
  int acolRad      = radBef.acol();
  int colEmt       = 0;
  int acolEmt      = 0;
  iSysSel          = dipSel->system;
  int iSysSelRec   = dipSel->systemRec;

  // Sometimes need to patch up colType in junction systems.
  int colTypeTmp   = dipSel->colType;
  int colTypeRec   = particleDataPtr->colType( recBef.id() );
  // Negate colour type if recoiler is initial-state quark.
  if (!recBef.isFinal()) colTypeRec = -colTypeRec;
  int colTypeRad   = particleDataPtr->colType( idRad );
  // Perform junction tests for all colour (anti)triplets.
  if (colTypeRec ==  1 && colTypeTmp > 0) colTypeTmp = -colTypeTmp;
  if (colTypeRec == -1 && colTypeTmp < 0) colTypeTmp = -colTypeTmp;
  if (colTypeRad ==  1 && colTypeTmp < 0) colTypeTmp = -colTypeTmp;
  if (colTypeRad == -1 && colTypeTmp > 0) colTypeTmp = -colTypeTmp;

  // Default OK for photon, photon_HV or gluon_HV emission.
  if (dipSel->flavour == 22 || dipSel->flavour == idHV) {
  // New colour tag required for gluon emission.
  } else if (dipSel->flavour == 21 && colTypeTmp > 0) {
    colEmt  = colRad;
    colRad  = event.nextColTag();
    acolEmt = colRad;
  } else if (dipSel->flavour == 21) {
    acolEmt = acolRad;
    acolRad = event.nextColTag();
    colEmt  = acolRad;
  // New flavours for g -> q qbar; split colours.
  } else if (colTypeTmp > 0) {
    idEmt   = dipSel->flavour ;
    idRad   = -idEmt;
    colEmt  = colRad;
    colRad  = 0;
  } else if (colTypeTmp < 0) {
    idEmt   = -dipSel->flavour ;
    idRad   = -idEmt;
    acolEmt = acolRad;
    acolRad = 0;
  // New flavours for gamma -> f fbar, and maybe also colours.
  } else if (dipSel->gamType == 1 && rndmPtr->flat() > 0.5) {
    idEmt   = -dipSel->flavour ;
    idRad   = -idEmt;
    if (idRad < 10) colRad = event.nextColTag();
    acolEmt = colRad;
  } else if (dipSel->gamType == 1) {
    idEmt   = dipSel->flavour ;
    idRad   = -idEmt;
    if (idEmt < 10) colEmt = event.nextColTag();
    acolRad = colEmt;
  }

  // Change fermion flavour by W emissions.
  int idRadSv = idRad;
  if (abs(idEmt) == 24) {
    if (rndmPtr->flat() > coupSMPtr->V2CKMsum(idRad)) return false;
    idRad = coupSMPtr->V2CKMpick(idRad);
  }

  // Construct kinematics in dipole rest frame:
  // begin simple (like g -> g g).
  double pTorig       = sqrt( dipSel->pT2);
  double eRadPlusEmt  = 0.5 * (dipSel->m2Dip + dipSel->m2 - dipSel->m2Rec)
    / dipSel->mDip;
  double e2RadPlusEmt = pow2(eRadPlusEmt);
  double pzRadPlusEmt = 0.5 * sqrtpos( pow2(dipSel->m2Dip - dipSel->m2
    - dipSel->m2Rec) - 4. * dipSel->m2 * dipSel->m2Rec ) / dipSel->mDip;
  double pT2corr = dipSel->m2 * (e2RadPlusEmt * dipSel->z * (1. - dipSel->z)
                      - 0.25 * dipSel->m2) / pow2(pzRadPlusEmt);
  double pTcorr       = sqrtpos( pT2corr );
  double pzRad        = (e2RadPlusEmt * dipSel->z - 0.5 * dipSel->m2)
                      / pzRadPlusEmt;
  double pzEmt        = (e2RadPlusEmt * (1. - dipSel->z) - 0.5 * dipSel->m2)
                      / pzRadPlusEmt;
  // Radiator flavour changed if W emission, so find new mass.
  double mRad         = (idRad == idRadSv) ? dipSel->mRad
                      : particleDataPtr->m0(idRad);
  double m2Rad        = pow2(mRad);
  double mEmt         = 0.;

  // Kinematics reduction for f -> f W/Z when m_f > 0 (and m_W/Z > 0)
  // or q -> q gamma_v when m_q > 0 and m_gamma_v > 0.
  if ( dipSel->weakType != 0
    || (abs(dipSel->colvType) == 1 && dipSel->mFlavour > 0.) ) {
    mEmt              = dipSel->mFlavour;
    if (pow2(mRad + mEmt) > dipSel->m2) return false;
    double m2Emt      = pow2(mEmt);
    double lambda     = sqrtpos( pow2(dipSel->m2 - m2Rad - m2Emt)
                      - 4. * m2Rad * m2Emt );
    kRad              = 0.5 * (dipSel->m2 - lambda + m2Emt - m2Rad)
                      / dipSel->m2;
    kEmt              = 0.5 * (dipSel->m2 - lambda + m2Rad - m2Emt)
                      / dipSel->m2;
    pTorig           *= 1. - kRad - kEmt;
    pTcorr           *= 1. - kRad - kEmt;
    double pzMove     = kRad * pzRad - kEmt * pzEmt;
    pzRad            -= pzMove;
    pzEmt            += pzMove;

  // Kinematics reduction for q -> q g/gamma/g_HV when m_q > 0.
  } else if (abs(dipSel->colType) == 1 || dipSel->chgType != 0
    || abs(dipSel->colvType) == 1) {
    pTorig           *= 1. - dipSel->m2Rad / dipSel->m2;
    pTcorr           *= 1. - dipSel->m2Rad / dipSel->m2;
    pzRad            += pzEmt * dipSel->m2Rad / dipSel->m2;
    pzEmt            *= 1. - dipSel->m2Rad / dipSel->m2;

  // Kinematics reduction for g -> q qbar or gamma -> f fbar when m_f > 0;
  } else if (abs(dipSel->flavour) < 20) {
    mEmt              = dipSel->mFlavour;
    mRad              = mEmt;
    double beta       = sqrtpos( 1. - 4. * pow2(mEmt) / dipSel->m2 );
    pTorig           *= beta;
    pTcorr           *= beta;
    pzRad             = 0.5 * ( (1. + beta) * pzRad + (1. - beta) * pzEmt );
    pzEmt             = pzRadPlusEmt - pzRad;
  }

  // Reject g emission where mass effects have reduced pT below cutoff.
  if (idEmt == 21 && pTorig < pTcolCut) return false;

  // Find rest frame and angles of original dipole.
  RotBstMatrix M;
  M.fromCMframe(pRadBef, pRecBef);
  RotBstMatrix M1;
  M1.fromCMframe(pRadBef, pRecBef);
  RotBstMatrix M2;
  M2.toCMframe(pRadBef, pRecBef);

  // Evaluate coefficient of azimuthal asymmetry from gluon polarization.
  findAsymPol( event, dipSel);

  // Begin construction of new dipole kinematics: pick azimuthal angle.
  Vec4 pRad, pEmt, pRec;
  double wtPhi = 1.;
  do {
    double phi = 2. * M_PI * rndmPtr->flat();

    // Define kinematics of branching in dipole rest frame.
    pRad = Vec4( pTcorr * cos(phi), pTcorr * sin(phi), pzRad,
      sqrt( pow2(pTcorr) + pow2(pzRad) + pow2(mRad) ) );
    pEmt = Vec4( -pRad.px(), -pRad.py(), pzEmt,
      sqrt( pow2(pTcorr) + pow2(pzEmt) + pow2(mEmt) ) );
    pRec = Vec4( 0., 0., -pzRadPlusEmt, sqrt( pow2(pzRadPlusEmt)
      + dipSel->m2Rec ) );

    // Rotate and boost dipole products to the event frame.
    pRad.rotbst(M);
    pEmt.rotbst(M);
    pRec.rotbst(M);

    // New: To avoid instabilities for violent boosts, ensure that an incoming
    // recoiler always has zero px and py.
    if (dipSel->isrType != 0) {
      if (abs(pRec.px()) > 0.) {
        double phixx = pRec.phi();
        RotBstMatrix rot_by_pphi;
        rot_by_pphi.rot(0.,-phixx);
        pRec.rotbst( rot_by_pphi);
        double thetaxx = pRec.theta();
        if ( pRec.px() < 0. ) thetaxx *= -1.;
        if ( pRec.pz() < 0.) thetaxx += M_PI;
        RotBstMatrix rot_by_ptheta;
        rot_by_ptheta.rot(-thetaxx, 0.);
        pRec.rotbst( rot_by_ptheta );
      }
    }

    // Azimuthal phi weighting: loop to new phi value if required.
    if (dipSel->asymPol != 0.) {
      Vec4 pAunt = event[dipSel->iAunt].p();
      double cosPhi = cosphi( pRad, pAunt, pRadBef );
      wtPhi = ( 1. + dipSel->asymPol * (2. * pow2(cosPhi) - 1.) )
        / ( 1. + abs(dipSel->asymPol) );
    }
  } while (wtPhi < rndmPtr->flat()) ;

  // Kinematics when recoiler is initial-state parton.
  int isrTypeNow  = dipSel->isrType;
  int isrTypeSave = isrTypeNow;
  if (!useLocalRecoilNow) isrTypeNow = 0;
  if (isrTypeNow != 0) pRec = 2. * recBef.p() - pRec;

  // New: Return if the x-value for the incoming recoiler is nonsense.
  if ( isrTypeNow != 0 && 2.*pRec.e()/event[0].m() > 1. ) {
    infoPtr->errorMsg("Error in SimpleTimeShower::branch: "
            "Larger than unity Bjorken x value");
    return false;
  }

  // PS dec 2010: check if radiator has flexible normalization
  bool isFlexible = dipSel->isFlexible;

  // Define new particles from dipole branching.
  double pTsel = sqrt(dipSel->pT2);
  Particle rad = Particle(idRad, 51, iRadBef, 0, 0, 0,
    colRad, acolRad, pRad, mRad, pTsel);
  Particle emt = Particle(idEmt, 51, iRadBef, 0, 0, 0,
    colEmt, acolEmt, pEmt, mEmt, pTsel);

  // Recoiler either in final or in initial state
  Particle rec = (isrTypeNow == 0)
    ? Particle(recBef.id(),  52, iRecBef, iRecBef, 0, 0,
      recBef.col(), recBef.acol(), pRec, dipSel->mRec, pTsel)
    : Particle(recBef.id(), -53, 0, 0, iRecBef, iRecBef,
      recBef.col(), recBef.acol(), pRec, 0., 0.);

  // Special checks to set weak particles status equal to 56.
  // This is needed for decaying the particles. Also set polarisation.
  if (emt.idAbs() == 23 || emt.idAbs() == 24) {
    emt.status(56);
    event[iRadBef].pol( dipSel->weakPol );
    rad.pol( dipSel->weakPol );
  }

  // Recover delayed shower-accept probability for uncertainty variations.
  double pAccept = dipSel->pAccept;

  // ME corrections can lead to branching being rejected.
  if (dipSel->MEtype > 0) {
    Particle& partner = (dipSel->iMEpartner == iRecBef)
      ? rec : event[dipSel->iMEpartner];
    double pMEC = findMEcorr( dipSel, rad, partner, emt);
    if (dipSel->MEtype >= 200 && dipSel->MEtype <= 210)
      pMEC *= findMEcorrWeak( dipSel, rad.p(), partner.p(), emt.p(),
        p3weak, p4weak, event[iRadBef].p(), event[iRecBef].p());
    pAccept *= pMEC;
  }

  // Decide if we are going to accept or reject this branching.
  // (Without wasting time generating random numbers if pAccept = 1.)
  bool acceptEvent = true;
  if (pAccept < 1.0) acceptEvent = (rndmPtr->flat() < pAccept);

  // Determine if this FSR is part of process or resonance showering
  bool inResonance = !partonSystemsPtr->hasInAB(iSysSel);

  // If doing uncertainty variations, calculate accept/reject reweightings.
  doUncertaintiesNow = doUncertainties;
  // Check if variations are allowed in MPIs.
  if (!uVarMPIshowers && iSysSel != 0 && !inResonance)
    doUncertaintiesNow = false;

  // Check if to allow variations in resonance decays.
  if (noResVariations && inResonance) doUncertaintiesNow = false;

  // Check if to allow variations in process.
  if (noProcVariations && iSysSel==0 && !inResonance)
    doUncertaintiesNow = false;

  // Check if below cutoff for calculating variations
  if ( dipSel->pT2 < uVarpTmin2 ) doUncertaintiesNow = false;

  // Early return if allowed.
  if (!doUncertaintiesNow && !acceptEvent) return false;

  // Rescatter: if the recoiling partner is not in the same system
  //            as the radiator, fix up intermediate systems (can lead
  //            to emissions being vetoed)
  if (allowRescatter && FIXRESCATTER && isInterleaved
    && iSysSel != iSysSelRec) {
    Vec4 pNew = rad.p() + emt.p();
    if (!rescatterPropagateRecoil(event, pNew)) return false;
  }

  // For photon-beam recoiler check that room for remnants after branching.
  if ( isrTypeNow != 0 ) {
    BeamParticle& beamRec = (isrTypeNow == 1) ? *beamAPtr : *beamBPtr;
    if ( beamRec.isGamma() ) {
      // If recoiler kinematics fixed by ISR can't act as recoiler.
      if ( !beamRec.resolvedGamma() ) return false;
      BeamParticle& beamOther = (isrTypeNow == 1) ? *beamBPtr : *beamAPtr;
      bool physical   = true;
      double xRec     = 2. * pRec.e() / (beamRec.e() + beamOther.e());
      double sCM      = m2( beamRec.p(), beamOther.p());
      double eCM      = sqrt(sCM);
      // One-remnant system.
      if ( !beamOther.resolvedGamma() ) {
        physical = beamRec.roomFor1Remnant(beamRec[0].id(), xRec, eCM);
      // Two-remnants systems.
      } else {
        physical = beamOther.roomFor2Remnants(beamRec[0].id(), xRec, eCM);
      }

      if (!physical) return false;
    }
  }

  // Save properties to be restored in case of user-hook veto of emission.
  int eventSizeOld = event.size();
  int iRadStatusV  = event[iRadBef].status();
  int iRadDau1V    = event[iRadBef].daughter1();
  int iRadDau2V    = event[iRadBef].daughter2();
  int iRecStatusV  = event[iRecBef].status();
  int iRecMot1V    = event[iRecBef].mother1();
  int iRecMot2V    = event[iRecBef].mother2();
  int iRecDau1V    = event[iRecBef].daughter1();
  int iRecDau2V    = event[iRecBef].daughter2();
  int beamOff1     = 1 + beamOffset;
  int beamOff2     = 2 + beamOffset;
  int ev1Dau1V     = event[beamOff1].daughter1();
  int ev2Dau1V     = event[beamOff2].daughter1();

  // Shower may occur at a displaced vertex.
  if (radBef.hasVertex()) {
    rad.vProd( radBef.vProd() );
    emt.vProd( radBef.vProd() );
  }
  if (recBef.hasVertex()) rec.vProd( recBef.vProd() );

  // Put new particles into the event record.
  int iRad = event.append(rad);
  int iEmt = event.append(emt);

  // Allow setting of new parton production vertex.
  if (doPartonVertex) partonVertexPtr->vertexFSR( iEmt, event);

  // Mark original dipole partons as branched and set daughters/mothers.
  event[iRadBef].statusNeg();
  event[iRadBef].daughters( iRad, iEmt);
  int iRec = 0;
  if (useLocalRecoilNow) {
    iRec = event.append(rec);
    if (isrTypeNow == 0) {
      event[iRecBef].statusNeg();
      event[iRecBef].daughters( iRec, iRec);
    } else {
      event[iRecBef].mothers( iRec, iRec);
      event[iRec].mothers( iRecMot1V, iRecMot2V);
      if (iRecMot1V == beamOff1) event[beamOff1].daughter1( iRec);
      if (iRecMot1V == beamOff2) event[beamOff2].daughter1( iRec);
    }

  // Global recoil: need to find relevant rotation+boost for recoilers:
  // boost+rotate to rest frame, boost along z axis, rotate+boost back.
  } else {
    RotBstMatrix MG = M;
    MG.invert();
    double pzRecBef = -0.5 * sqrtpos( pow2(dipSel->m2Dip - dipSel->m2Rad
      - dipSel->m2Rec) - 4. * dipSel->m2Rad * dipSel->m2Rec ) / dipSel->mDip;
    double eRecBef  = sqrt( pow2(pzRecBef) + dipSel->m2Rec);
    double pzRecAft = -pzRadPlusEmt;
    double eRecAft  = sqrt( pow2(pzRecAft) + dipSel->m2Rec);
    MG.bst( Vec4(0., 0., pzRecBef, eRecBef), Vec4(0., 0., pzRecAft, eRecAft) );
    MG.rotbst( M);

    // Global recoil: copy particles, and rotate+boost momenta (not vertices).
    for (int iG = 0; iG < int(iGRecBef.size()); ++iG) {
      iRec = event.copy( iGRecBef[iG], 52);
      iGRec.push_back( iRec);
      Vec4 pGRec = event[iRec].p();
      pGRec.rotbst( MG);
      event[iRec].p( pGRec);
    }
  }

  // Allow veto of branching. If so restore event record to before emission.
  if ( (canVetoEmission && userHooksPtr->doVetoFSREmission( eventSizeOld,
    event, iSysSel, inResonance))
    || (canMergeFirst && mergingHooksPtr->doVetoEmission( event )) ) {
    event.popBack( event.size() - eventSizeOld);
    event[iRadBef].status( iRadStatusV);
    event[iRadBef].daughters( iRadDau1V, iRadDau2V);
    if (useLocalRecoilNow && isrTypeNow == 0) {
      event[iRecBef].status( iRecStatusV);
      event[iRecBef].daughters( iRecDau1V, iRecDau2V);
    } else if (useLocalRecoilNow) {
      event[iRecBef].mothers( iRecMot1V, iRecMot2V);
      if (iRecMot1V == beamOff1) event[beamOff1].daughter1( ev1Dau1V);
      if (iRecMot1V == beamOff2) event[beamOff2].daughter1( ev2Dau1V);
    } else {
      for (int iG = 0; iG < int(iGRecBef.size()); ++iG) {
        event[iGRecBef[iG]].statusPos();
        event[iGRecBef[iG]].daughters( 0, 0);
      }
    }
    return false;
  }
  // Default settings for uncertainty calculations.
  double weight = 1.;
  double vp = 0.;
  bool vetoedEnhancedEmission = false;

  // Calculate event weight for enhanced emission rate.
  if (canEnhanceET) {
    // Check if emission weight was enhanced. Get enhance weight factor.
    bool foundEnhance = false;
    // Move backwards as last elements have highest pT, thus are chosen
    // splittings.
    for ( map<double,pair<string,double> >::reverse_iterator
          it = enhanceFactors.rbegin();
          it != enhanceFactors.rend(); ++it ){
      if (splittingNameSel.find(it->second.first) != string::npos
        && abs(it->second.second-1.0) > 1e-9) {
        foundEnhance = true;
        weight       = it->second.second;
        vp           = userHooksPtr->vetoProbability(splittingNameSel);
        break;
      }
    }

    // Check emission veto.
    if (foundEnhance && rndmPtr->flat() < vp ) vetoedEnhancedEmission = true;
    // Calculate new event weight.
    double rwgt = 1.;
    if (foundEnhance && vetoedEnhancedEmission) rwgt *= (1.-1./weight)/vp;
    else if (foundEnhance) rwgt *= 1./((1.-vp)*weight);

    // Reset enhance factors after usage.
    enhanceFactors.clear();

    // Set events weights, so that these could be used externally.
    double wtOld = userHooksPtr->getEnhancedEventWeight();
    if (!doTrialNow && canEnhanceEmission && !doUncertaintiesNow)
      userHooksPtr->setEnhancedEventWeight(wtOld*rwgt);
    if ( doTrialNow && canEnhanceTrial)
      userHooksPtr->setEnhancedTrial(sqrt(dipSel->pT2), weight);
    // Increment counter to handle counting of rejected emissions.
    if (vetoedEnhancedEmission && canEnhanceEmission) infoPtr->addCounter(40);
  }

  // Emission veto is a phase space restriction, and should not be included
  // in the uncertainty calculation.
  if (vetoedEnhancedEmission) acceptEvent = false;
  if (doUncertaintiesNow) calcUncertainties( acceptEvent, pAccept, weight, vp,
    dipSel, &rad, &emt, &rec);

  // Return false if we decided to reject this branching.
  // Veto if necessary.
  if ( (vetoedEnhancedEmission && canEnhanceEmission) || !acceptEvent) {
    event.popBack( event.size() - eventSizeOld);
    event[iRadBef].status( iRadStatusV);
    event[iRadBef].daughters( iRadDau1V, iRadDau2V);
    if (useLocalRecoilNow && isrTypeNow == 0) {
      event[iRecBef].status( iRecStatusV);
      event[iRecBef].daughters( iRecDau1V, iRecDau2V);
    } else if (useLocalRecoilNow) {
      event[iRecBef].mothers( iRecMot1V, iRecMot2V);
      if (iRecMot1V == beamOff1) event[beamOff1].daughter1( ev1Dau1V);
      if (iRecMot1V == beamOff2) event[beamOff2].daughter1( ev2Dau1V);
    } else {
      for (int iG = 0; iG < int(iGRecBef.size()); ++iG) {
        event[iGRecBef[iG]].statusPos();
        event[iGRecBef[iG]].daughters( 0, 0);
      }
    }
    return false;
  }

  // For global recoil restore the one nominal recoiler, for bookkeeping.
  if (!useLocalRecoilNow) {
    iRec = iRecBef;
    for (int iG = 0; iG < int(iGRecBef.size()); ++iG)
    if (iGRecBef[iG] == iRecBef) iRec = iGRec[iG];
  }

  // For initial-state recoiler also update beam and sHat info.
  if (isrTypeNow != 0) {
    BeamParticle& beamRec = (isrTypeNow == 1) ? *beamAPtr : *beamBPtr;
    double xOld = beamRec[iSysSelRec].x();
    double xRec = 2. * pRec.e() / (beamAPtr->e() + beamBPtr->e());
    beamRec[iSysSelRec].iPos( iRec);
    beamRec[iSysSelRec].x( xRec);
    partonSystemsPtr->setSHat( iSysSelRec,
    partonSystemsPtr->getSHat(iSysSelRec) * xRec / xOld);
  }

  // For global recoil: if everything went as expected, remove the line
  // from the list of "hard lines" that are allowed to use global recoil.
  if ( !useLocalRecoilNow || nGlobal >= nMaxGlobalBranch) {
    bool doRemove=true;
    while ( doRemove ) {
      bool hasRemoved = false;
      for (int iHard = 0; iHard < int(hardPartons.size()); ++iHard)
        if ( event[dipSel->iRadiator].isAncestor(hardPartons[iHard]) ) {
          hardPartons.erase( hardPartons.begin() + iHard );
          hasRemoved = true;
          break;
        }
      doRemove = hasRemoved;
    }
  }

  // Update number of splittings that have been produced with global recoil.
  if ( !useLocalRecoilNow ) ++nGlobal;

  // Photon emission: update to new dipole ends; add new photon "dipole".
  if (dipSel->flavour == 22) {
    dipSel->iRadiator = iRad;
    dipSel->iRecoiler = iRec;
    // When recoiler was uncharged particle, in resonance decays,
    // assign recoil to emitted photons.
    if (recoilToColoured && inResonance && event[iRec].chargeType() == 0)
      dipSel->iRecoiler = iEmt;
    dipSel->pTmax = pTsel;
    if (doQEDshowerByGamma) dipEnd.push_back( TimeDipoleEnd(iEmt, iRad,
      pTsel, 0, 0, 1, 0, 0, iSysSel, 0) );

  // Gluon emission: update both dipole ends and add two new ones.
  } else if (dipSel->flavour == 21) {
    dipSel->iRadiator  = iRad;
    dipSel->iRecoiler  = iEmt;
    dipSel->systemRec  = iSysSel;
    dipSel->isrType    = 0;
    dipSel->pTmax      = pTsel;
    // Optionally also kill ME corrections after first emission.
    if (!doMEafterFirst) dipSel->MEtype = 0;
    // PS dec 2010: check normalization of radiating dipole
    // Dipole corresponding to the newly created colour tag has normal strength
    double flexFactor  = (isFlexible) ? dipSel->flexFactor : 1.0;
    dipSel->isFlexible = false;
    for (int i = 0; i < int(dipEnd.size()); ++i) {
      if (dipEnd[i].iRadiator == iRecBef && dipEnd[i].iRecoiler == iRadBef
        && dipEnd[i].colType != 0) {
        dipEnd[i].iRadiator = iRec;
        dipEnd[i].iRecoiler = iEmt;
        // Optionally also kill ME corrections after first emission.
        if (!doMEafterFirst) dipEnd[i].MEtype = 0;
        // Strive to match colour to anticolour inside closed system.
        if ( !isFlexible && dipEnd[i].colType * dipSel->colType > 0)
          dipEnd[i].iRecoiler = iRad;
        dipEnd[i].pTmax = pTsel;
        // PS dec 2010: if the (iRadBef,iRecBef) dipole was flexible, the
        // same should be true for this (opposite) end. If so, this end keeps
        // the modified normalization, so we shouldn't need to do anything.
      }
      // Weak shower can have gluons as recoiler. Always choose
      // the outgoing gluon that produces the highest invariant mass.
      if (event[iRadBef].id() == 21 && dipEnd[i].iRecoiler == iRadBef
         && dipEnd[i].weakType != 0) {
        double m1 = (event[iRad].p()+event[dipEnd[i].iRadiator].p()).m2Calc();
        double m2 = (event[iEmt].p()+event[dipEnd[i].iRadiator].p()).m2Calc();
        dipEnd[i].iRecoiler = (m1 > m2) ? iRad : iEmt;
        dipEnd[i].iMEpartner = dipEnd[i].iRecoiler;
      }
    }
    int colType = (dipSel->colType > 0) ? 2 : -2 ;
    // When recoiler was uncoloured particle, in resonance decays,
    // assign recoil to coloured particle.
    int iRecMod = iRec;
    if (recoilToColoured && inResonance && event[iRec].col() == 0
      && event[iRec].acol() == 0) iRecMod = iRad;
    dipEnd.push_back( TimeDipoleEnd(iEmt, iRecMod, pTsel,
       colType, 0, 0, 0, isrTypeSave, iSysSel, 0));
    dipEnd.back().systemRec = iSysSelRec;
    // PS dec 2010: the (iEmt,iRec) dipole "inherits" flexible normalization
    if (isFlexible) {
      dipEnd.back().isFlexible = true;
      dipEnd.back().flexFactor = flexFactor;
    }
    dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel,
      -colType, 0, 0, 0, 0, iSysSel, 0));

  // Gluon branching to q qbar: update current dipole and other of gluon.
  } else if (dipSel->colType != 0) {
    for (int i = 0; i < int(dipEnd.size()); ++i) {
      // Strive to match colour to anticolour inside closed system.
      if ( !isFlexible && dipEnd[i].iRecoiler == iRadBef
        && dipEnd[i].colType * dipSel->colType < 0 )
        dipEnd[i].iRecoiler = iEmt;
      if (dipEnd[i].iRadiator == iRadBef && abs(dipEnd[i].colType) == 2) {
        dipEnd[i].colType /= 2;
        if (dipEnd[i].system != dipEnd[i].systemRec) continue;

        // Note: gluino -> quark + squark gives a deeper radiation dip than
        // the more obvious alternative photon decay, so is more realistic.
        dipEnd[i].MEtype = (doMEcorrections && doMEafterFirst) ? 66 : 0;
        if (&dipEnd[i] == dipSel) dipEnd[i].iMEpartner = iRad;
        else                      dipEnd[i].iMEpartner = iEmt;
      }
      // Choose recoiler to Z/W to get largest mass.
      if (event[iRadBef].id() == 21 && dipEnd[i].iRecoiler == iRadBef
         && dipEnd[i].weakType != 0) {
        double m1 = (event[iRad].p()+event[dipEnd[i].iRadiator].p()).m2Calc();
        double m2 = (event[iEmt].p()+event[dipEnd[i].iRadiator].p()).m2Calc();
        dipEnd[i].iRecoiler = (m1 > m2) ? iRad : iEmt;
        dipEnd[i].iMEpartner = dipEnd[i].iRecoiler;
      }
    }
    dipSel->iRadiator = iEmt;
    dipSel->iRecoiler = iRec;
    dipSel->pTmax     = pTsel;

    // Gluon branching to q qbar: also add two charge dipole ends.
    // Note: gluino -> quark + squark gives a deeper radiation dip than
    // the more obvious alternative photon decay, so is more realistic.
    if (doQEDshowerByQ) {
      int chgType = event[iRad].chargeType();
      int meType = (doMEcorrections && doMEafterFirst) ? 66 : 0;
      dipEnd.push_back( TimeDipoleEnd(iRad, iEmt, pTsel,
        0,  chgType, 0, 0, 0, iSysSel, meType, iEmt));
      dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel,
        0, -chgType, 0, 0, 0, iSysSel, meType, iRad));
    }

    // Gluon branching to q qbar: also add weak dipoles.
    // Randomly decided whether to use left or right quarks.
    if (doWeakShower && iSysSel == 0 &&
      !(hasWeaklyRadiated && singleWeakEmission)) {
      int weakPol = (rndmPtr->flat() > 0.5) ? -1 : 1;
      event[iRad].pol(weakPol);
      event[iEmt].pol(weakPol);
      if ((weakMode == 0 || weakMode == 1) && weakPol == -1) {
        dipEnd.push_back( TimeDipoleEnd(iRad, iEmt, pTsel,
          0, 0, 0, 1, 0, iSysSel, 200, iEmt, weakPol) );
        dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel,
          0, 0, 0, 1, 0, iSysSel, 200, iRad, weakPol) );
      }
      if (weakMode == 0 || weakMode == 2) {
        dipEnd.push_back( TimeDipoleEnd(iRad, iEmt, pTsel,
          0, 0, 0, 2, 0, iSysSel, 205, iEmt, weakPol) );
        dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel,
          0, 0, 0, 2, 0, iSysSel, 205, iRad, weakPol) );
      }
    }

  // Photon branching to f fbar: inactivate photon "dipole";
  // optionally add new charge and colour dipole ends.
  // (But not W or Z ends, since W/Z are heavier than gamma*.)
  } else if (dipSel->gamType != 0) {
    dipSel->gamType = 0;
    int chgType = event[iRad].chargeType();
    int colType = event[iRad].colType();
    // MEtype = 102 for charge in vector decay.
    if ( chgType != 0 && ( ( doQEDshowerByQ && colType != 0 )
      || ( doQEDshowerByL && colType == 0 ) ) ) {
      int MEtype = (doMEcorrections && doMEafterFirst) ? 102 : 0;
      dipEnd.push_back( TimeDipoleEnd(iRad, iEmt, pTsel,
        0,  chgType, 0, 0, 0, iSysSel, MEtype, iEmt) );
      dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel,
        0, -chgType, 0, 0, 0, iSysSel, MEtype, iRad) );
    }
    // MEtype = 11 for colour in vector decay.
    if (colType != 0 && doQCDshower) {
      int MEtype = (doMEcorrections && doMEafterFirst) ? 11 : 0;
      dipEnd.push_back( TimeDipoleEnd(iRad, iEmt, pTsel,
         colType, 0, 0, 0, 0, iSysSel, MEtype, iEmt) );
      dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel,
        -colType, 0, 0, 0, 0, iSysSel, MEtype, iRad) );
    }

  // Photon_HV emission: update to new dipole ends.
  } else if (dipSel->flavour == 4900022) {
    dipSel->iRadiator = iRad;
    dipSel->iRecoiler = iRec;
    dipSel->pTmax = pTsel;

  // Gluon_HV emission: update to new dipole ends.
  } else if (dipSel->flavour == 4900021) {
    dipSel->iRadiator = iRad;
    dipSel->iRecoiler = iEmt;
    dipSel->pTmax     = pTsel;
    for (int i = 0; i < int(dipEnd.size()); ++i)
    if (dipEnd[i].iRadiator == iRecBef && dipEnd[i].iRecoiler == iRadBef
      && dipEnd[i].isHiddenValley) {
      dipEnd[i].iRadiator = iRec;
      dipEnd[i].iRecoiler = iEmt;
      dipEnd[i].pTmax     = pTsel;
    }
    int colvType = (dipSel->colvType > 0) ? 2 : -2 ;
    dipEnd.push_back( TimeDipoleEnd(iEmt, iRec, pTsel,
      0, 0, 0, 0, isrTypeSave, iSysSel, 0, -1, 0, false, true, colvType) );
    dipEnd.back().systemRec = iSysSelRec;
    dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel,
      0, 0, 0, 0, 0, iSysSel, 0, -1, 0, false, true, -colvType) );

  // W/Z emission, if only a single weak emission is allowed.
  } else if (dipSel->weakType != 0) {
    hasWeaklyRadiated = true;
    if (singleWeakEmission)
      for (int i = 0; i < int(dipEnd.size()); ++i) dipEnd[i].weakType = 0;
  }

  // Copy or set lifetime for new final state.
  if (event[iRad].id() == event[iRadBef].id()) {
    event[iRad].tau( event[iRadBef].tau() );
  } else {
    event[iRad].tau( event[iRad].tau0() * rndmPtr->exp() );
  }
  event[iRec].tau( event[iRecBef].tau() );
  event[iEmt].tau( event[iEmt].tau0() * rndmPtr->exp() );

  // Now update other dipoles that also involved the radiator or recoiler.
  for (int i = 0; i < int(dipEnd.size()); ++i) {
    // PS dec 2010: if radiator was flexible and now is normal, there may
    // be other flexible dipoles that need updating.
    if (isFlexible && !dipSel->isFlexible && dipEnd[i].isFlexible) {
      if (dipEnd[i].iRecoiler  == iRadBef) dipEnd[i].iRecoiler = iEmt;
      if (dipEnd[i].iRadiator  == iRadBef) {
        dipEnd[i].iRadiator = iEmt;
        if (dipEnd[i].colType == 1 && dipSel->flavour == 21)
          dipEnd[i].colType = 2;
        if (dipEnd[i].colType ==-1 && dipSel->flavour == 21)
          dipEnd[i].colType =-2;
      }
    }
    if (dipEnd[i].iRadiator  == iRadBef) dipEnd[i].iRadiator  = iRad;
    if (dipEnd[i].iRecoiler  == iRadBef) dipEnd[i].iRecoiler  = iRad;
    if (dipEnd[i].iMEpartner == iRadBef) dipEnd[i].iMEpartner = iRad;
    if (useLocalRecoilNow) {
      if (dipEnd[i].iRadiator  == iRecBef) dipEnd[i].iRadiator  = iRec;
      if (dipEnd[i].iRecoiler  == iRecBef) dipEnd[i].iRecoiler  = iRec;
      if (dipEnd[i].iMEpartner == iRecBef) dipEnd[i].iMEpartner = iRec;
    } else {
      for (int iG = 0; iG < int(iGRecBef.size()); ++iG) {
        if (dipEnd[i].iRadiator  == iGRecBef[iG])
            dipEnd[i].iRadiator  =  iGRec[iG];
        if (dipEnd[i].iRecoiler  == iGRecBef[iG])
            dipEnd[i].iRecoiler  =  iGRec[iG];
        if (dipEnd[i].iMEpartner == iGRecBef[iG])
            dipEnd[i].iMEpartner =  iGRec[iG];
      }
    }
  }

  // PS Apr 2011
  // Update any junctions downstream of this branching (if necessary)
  // (This happens, e.g., via LHEF, when adding showers to intermediate
  //  coloured resonances whose decays involved junctions)
  for (int iJun = 0; iJun < event.sizeJunction(); iJun++) {
    // Number of incoming colour lines for this junction.
    int nIncoming = (event.kindJunction(iJun)-1)/2;
    // Check radiator colour or anticolour, depending on junction kind
    // (if junction, incoming = anticolours, and vice versa)
    int colChk = 0;
    colChk = ( event.kindJunction(iJun) % 2 == 0 )
           ? event[iRadBef].col() : event[iRadBef].acol();
    // Loop over incoming junction ends
    for (int iCol = 0; iCol < nIncoming; iCol++) {
      int colJun = event.colJunction( iJun, iCol);
      // If match, update junction end with new upstream (anti)colour
      if (colJun == colChk) {
        int colNew = 0;
        if ( event.kindJunction(iJun) % 2 == 0 ) colNew = colRad;
        else colNew = acolRad;
        event.colJunction( iJun, iCol, colNew );
      }
    }
  }

  // Finally update the list of all partons in all systems.
  partonSystemsPtr->replace(iSysSel, iRadBef, iRad);
  partonSystemsPtr->addOut(iSysSel, iEmt);
  if (useLocalRecoilNow)
    partonSystemsPtr->replace(iSysSelRec, iRecBef, iRec);
  else {
    for (int iG = 0; iG < int(iGRecBef.size()); ++iG)
    partonSystemsPtr->replace(iSysSel, iGRecBef[iG], iGRec[iG]);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Initialize the choices of uncertainty variations of the shower.

bool SimpleTimeShower::initUncertainties() {

  if( infoPtr->nWeights() > 1 ) return(nUncertaintyVariations);

  // Populate lists of uncertainty variations for SimpleTimeShower, by keyword.
  uVarMuSoftCorr = settingsPtr->flag("UncertaintyBands:muSoftCorr");
  dASmax         = settingsPtr->parm("UncertaintyBands:deltaAlphaSmax");
  // Variations handled by SpaceShower.
  varPDFplus    = &infoPtr->varPDFplus;
  varPDFminus   = &infoPtr->varPDFminus;
  varPDFmember  = &infoPtr->varPDFmember;

  // Reset uncertainty variation maps.
  varG2GGmuRfac.clear();    varG2GGcNS.clear();
  varQ2QGmuRfac.clear();    varQ2QGcNS.clear();
  varX2XGmuRfac.clear();    varX2XGcNS.clear();
  varG2QQmuRfac.clear();    varG2QQcNS.clear();

  vector<string> keys;
  keys.push_back("fsr:murfac");
  keys.push_back("fsr:g2gg:murfac");
  keys.push_back("fsr:q2qg:murfac");
  keys.push_back("fsr:x2xg:murfac");
  keys.push_back("fsr:g2qq:murfac");
  keys.push_back("fsr:cns");
  keys.push_back("fsr:g2gg:cns");
  keys.push_back("fsr:q2qg:cns");
  keys.push_back("fsr:x2xg:cns");
  keys.push_back("fsr:g2qq:cns");
  // Store number of QCD variations (as separator to QED ones).
  int nKeysQCD=keys.size();

  // Get uncertainty variations from Settings (as list of strings to parse).
  vector<string> uVars = settingsPtr->wvec("UncertaintyBands:List");
  size_t varSize = uVars.size();
  nUncertaintyVariations = int(uVars.size());
  if (nUncertaintyVariations == 0) return false;
  vector<string> uniqueVars;

  // Parse each string in uVars to look for recognized keywords.
  for (size_t iWeight = 0; iWeight < varSize; ++iWeight) {
    // Convert to lowercase (to be case-insensitive). Also remove "=" signs
    // and extra spaces, so "key=value", "key = value" mapped to "key value"
    string uVarString = toLower(uVars[iWeight]);
    while (uVarString.find(" ") == 0) uVarString.erase( 0, 1);
    int iEnd = uVarString.find(" ", 0);
    uVarString.erase(0,iEnd+1);
    while (uVarString.find("=") != string::npos) {
      int firstEqual = uVarString.find_first_of("=");
      string testString = uVarString.substr(0, firstEqual);
      iEnd = uVarString.find_first_of(" ", 0);
      if( iEnd<0 ) iEnd = uVarString.length();
      string insertString = uVarString.substr(0,iEnd);
      // does the key match an fsr one?
      if( find(keys.begin(), keys.end(), testString) != keys.end() ) {
        if( uniqueVars.size() == 0 ) {
          uniqueVars.push_back(insertString);
        } else if ( find(uniqueVars.begin(), uniqueVars.end(), insertString)
        == uniqueVars.end() ) {
          uniqueVars.push_back(insertString);
        }
      }
      uVarString.erase(0,iEnd+1);
    }
  }

  nUncertaintyVariations = int(uniqueVars.size());

  // Only perform for the first call to Timeshower
  if (infoPtr->nWeights() <= 1.) {
    infoPtr->setNWeights( nUncertaintyVariations + 1 );
    infoPtr->setWeightLabel( 0, "Baseline");
    for(int iWeight = 1; iWeight <= nUncertaintyVariations; ++iWeight) {
      string uVarString = uniqueVars[iWeight-1];
      infoPtr->setWeightLabel(iWeight, uVarString);

      while (uVarString.find("=") != string::npos) {
        int firstEqual = uVarString.find_first_of("=");
        uVarString.replace(firstEqual, 1, " ");
      }
      while (uVarString.find("  ") != string::npos)
        uVarString.erase( uVarString.find("  "), 1);
      if (uVarString == "" || uVarString == " ") continue;

      // Loop over all keywords.
      int nRecognizedQCD = 0;
      for (int iWord = 0; iWord < int(keys.size()); ++iWord) {
        // Transform string to lowercase to avoid case-dependence.
        string key = toLower(keys[iWord]);
        // Skip if empty or keyword not found.
        if (uVarString.find(key) == string::npos) continue;
        // Extract variation value/factor.
        int iKey = uVarString.find(key);
        int iBeg = uVarString.find(" ", iKey) + 1;
        int iEnd = uVarString.find(" ", iBeg);
        string valueString = uVarString.substr(iBeg, iEnd - iBeg);
        stringstream ss(valueString);
        double value;
        ss >> value;
        if (!ss) continue;

        // Store (iWeight,value) pairs
        // RECALL: use lowercase for all keys here (since converted above).
        if (key == "fsr:murfac" || key == "fsr:g2gg:murfac")
          varG2GGmuRfac[iWeight] = value;
        if (key == "fsr:murfac" || key == "fsr:q2qg:murfac")
          varQ2QGmuRfac[iWeight] = value;
        if (key == "fsr:murfac" || key == "fsr:x2xg:murfac")
          varX2XGmuRfac[iWeight] = value;
        if (key == "fsr:murfac" || key == "fsr:g2qq:murfac")
          varG2QQmuRfac[iWeight] = value;
        if (key == "fsr:cns" || key == "fsr:g2gg:cns")
          varG2GGcNS[iWeight] = value;
        if (key == "fsr:cns" || key == "fsr:q2qg:cns")
          varQ2QGcNS[iWeight] = value;
        if (key == "fsr:cns" || key == "fsr:x2xg:cns")
          varX2XGcNS[iWeight] = value;
        if (key == "fsr:cns" || key == "fsr:g2qq:cns")
          varG2QQcNS[iWeight] = value;
        // Tell that we found at least one recognized and parseable keyword.
        if (iWord < nKeysQCD) nRecognizedQCD++;
      } // End loop over QCD keywords

      // Tell whether this uncertainty variation contained >= 1 QCD variation.
      if (nRecognizedQCD > 0) ++nVarQCD;
    } // End loop over UVars.
  }
  infoPtr->initUncertainties(&uVars);
  // Let the calling function know if we found anything.
  return (nUncertaintyVariations > 0);
}


//==========================================================================

// Calculate uncertainties for the current event.

void SimpleTimeShower::calcUncertainties(bool accept, double pAccept,
  double enhance, double vp, TimeDipoleEnd* dip, Particle* radPtr,
  Particle* emtPtr, Particle* recPtr) {

  // Sanity check.
  if (!doUncertainties || !doUncertaintiesNow || nUncertaintyVariations <= 0)
    return;

  // Define pointer and iterator to loop over the contents of each
  // (iWeight,value) map.
  map<int,double>* varPtr=0;
  map<int,double>::iterator itVar;
  // Make sure we have a dummy to point to if no map to be used.
  map<int,double> dummy;     dummy.clear();

  int numWeights = infoPtr->nWeights();
  // Store uncertainty variation factors, initialised to unity.
  // Make vector sizes + 1 since 0 = default and variations start at 1.
  vector<double> uVarFac(numWeights, 1.0);
  vector<bool> doVar(numWeights, false);

  // For the case of biasing, the nominal weight might not be unity.
  doVar[0] = true;
  uVarFac[0] = 1.0;

  // Extract relevant quantities.
  int idEmt = emtPtr->id();
  int idRad = radPtr->id();

  // QCD variations.
  if (dip->colType != 0) {

    // QCD renormalization-scale variations.
    if (alphaSorder == 0) varPtr = &dummy;
    else if (idEmt == 21 && idRad == 21) varPtr = &varG2GGmuRfac;
    else if (idEmt == 21 && abs(idRad) <= uVarNflavQ)
      varPtr = &varQ2QGmuRfac;
    else if (idEmt == 21) varPtr = &varX2XGmuRfac;
    else if (abs(idRad) <= nGluonToQuark && abs(idEmt) <= nGluonToQuark)
      varPtr = &varG2QQmuRfac;
    else varPtr = &dummy;
    for (itVar = varPtr->begin(); itVar != varPtr->end(); ++itVar) {
      int iWeight   = itVar->first;
      double valFac = itVar->second;
      double muR2 = renormMultFac * dip->pT2;
      double alphaSbaseline = alphaS.alphaS(muR2);
      // Correction-factor alphaS.
      double muR2var = max(1.1 * Lambda3flav2, pow2(valFac) * muR2);
      double alphaSratio = alphaS.alphaS(muR2var) / alphaSbaseline;
      // Apply soft correction factor to X2XG.
      double facCorr = 1.;
      if (idEmt == 21 && uVarMuSoftCorr) {
        // Use smallest alphaS and b0, to make the compensation conservative.
        int nf = 5;
        if (dip->pT2 < pow2(mc)) nf = 3;
        else if (dip->pT2 < pow2(mb)) nf = 4;
        double alphaScorr = alphaS.alphaS(dip->m2Dip);
        double facSoft    = alphaScorr * (33. - 2. * nf) / (6. * M_PI);
        double zeta = 1. - dip->z;
        if (idRad == 21) zeta = min(dip->z, 1. - dip->z);
        facCorr = 1. + (1. - zeta) * facSoft * log(valFac);
      }
      // Apply correction factor here for emission processes.
      double alphaSfac   = alphaSratio * facCorr;
      // Limit absolute variation to +/- deltaAlphaSmax.
      if (alphaSfac > 1.)
        alphaSfac = min(alphaSfac, (alphaSbaseline + dASmax) / alphaSbaseline);
      else if (alphaSbaseline > dASmax)
        alphaSfac = max(alphaSfac, (alphaSbaseline - dASmax) / alphaSbaseline);
      uVarFac[iWeight] *= alphaSfac;
      doVar[iWeight] = true;
    }

    // QCD finite-term variations (only when no MECs and above pT threshold).
    if (dip->MEtype != 0 || dip->pT2 < pow2(cNSpTmin) ) varPtr = &dummy;
    else if (idEmt == 21 && idRad == 21) varPtr = &varG2GGcNS;
    else if (idEmt == 21 && abs(idRad) <= uVarNflavQ) varPtr = &varQ2QGcNS;
    else if (idEmt == 21) varPtr = &varX2XGcNS;
    else if (abs(idRad) <= nGluonToQuark && abs(idEmt) <= nGluonToQuark)
      varPtr = &varG2QQcNS;
    else varPtr = &dummy;
    for (itVar = varPtr->begin(); itVar != varPtr->end(); ++itVar) {
      int iWeight   = itVar->first;
      double valFac = itVar->second;
      // Correction-factor alphaS.
      double z   = dip->z;
      double Q2  = dip->m2;
      // Virtuality for massive radiators.
      if (abs(idRad) >= 4 && idRad != 21) Q2 = max(1., Q2-radPtr->m2());
      double yQ  = Q2 / dip->m2Dip;
      double num = yQ * valFac;
      double denom = 1.;
      // G->GG.
      if (idEmt == 21 && idRad == 21)
        denom = pow2(1. - z * (1.-z)) / (z*(1.-z));
      // Q->QG.
      else if (idEmt == 21)
          denom = (1. + pow2(z)) / (1. - z);
      // G->QQ.
      else
          denom = pow2(z) + pow2(1. - z);
      // Compute reweight ratio.
      uVarFac[iWeight] *= 1. + num / denom;
      doVar[iWeight] = true;
    }

    // PDF variations for dipoles that connect to the initial state.
    if ( dip->isrType != 0 ){
      if ( !varPDFplus->empty() || !varPDFminus->empty()
        || !varPDFmember->empty() ) {
        // Evaluation of new daughter and mother PDF's.
        double scale2 = (useFixedFacScale) ? fixedFacScale2
          : factorMultFac * dip->pT2;
        BeamParticle& beam  = (dip->isrType == 1) ? *beamAPtr : *beamBPtr;
        int iSysRec   = dip->systemRec;
        double xOld   = beam[iSysRec].x();
        double xNew   = xOld * (1. + (dip->m2 - dip->m2Rad)
                             / (dip->m2Dip - dip->m2Rad));
        int idRec     = recPtr->id();
        int valSea = (beam[iSysSel].isValence()) ? 1 : 0;
        if( beam[iSysSel].isUnmatched() ) valSea = 2;
        beam.calcPDFEnvelope( make_pair(idRec,idRec),
                              make_pair(xNew,xOld), scale2, valSea);
        PDF::PDFEnvelope ratioPDFEnv = beam.getPDFEnvelope( );
        //
        varPtr = varPDFplus;
        for (itVar = varPtr->begin(); itVar != varPtr->end(); ++itVar) {
          int iWeight   = itVar->first;
          uVarFac[iWeight] *= 1.0 + min(ratioPDFEnv.errplusPDF
            / ratioPDFEnv.centralPDF,0.5);
          doVar[iWeight] = true;
        }
        //
        varPtr = varPDFminus;
        for (itVar = varPtr->begin(); itVar != varPtr->end(); ++itVar) {
          int iWeight   = itVar->first;
          uVarFac[iWeight] *= max(.01,1.0 - min(ratioPDFEnv.errminusPDF
            / ratioPDFEnv.centralPDF,0.5));
          doVar[iWeight] = true;
        }
        varPtr = varPDFmember;
        for (itVar = varPtr->begin(); itVar != varPtr->end(); ++itVar) {
          int iWeight   = itVar->first;
          int member    = int( itVar->second );
          uVarFac[iWeight] *= max(.01,ratioPDFEnv.pdfMemberVars[member]
            / ratioPDFEnv.centralPDF);
          doVar[iWeight] = true;
        }
      }
    }

  }

  // Ensure 0 < PacceptPrime < 1 (with small margins).
  // Skip the central weight, so as to avoid confusion
  for (int iWeight = 1; iWeight<=nUncertaintyVariations; ++iWeight) {
    if (!doVar[iWeight]) continue;
    double pAcceptPrime = pAccept * uVarFac[iWeight];
    if (pAcceptPrime > PROBLIMIT && dip->colType != 0) {
      uVarFac[iWeight] *= PROBLIMIT / pAcceptPrime;
    }
  }

  // Apply reject or accept reweighting factors according to input decision.
  for (int iWeight = 0; iWeight <= nUncertaintyVariations; ++iWeight) {
    if (!doVar[iWeight]) continue;
    // If trial accepted: apply ratio of accept probabilities.
    if (accept) infoPtr->reWeight(iWeight,
      uVarFac[iWeight] / ((1.0 - vp) * enhance) );
    // If trial rejected : apply Sudakov reweightings.
    else {
      // Check for near-singular denominators (indicates too few failures,
      // and hence would need to increase headroom).
      double denom = 1. - pAccept*(1.0 - vp);
      if (denom < REJECTFACTOR) {
        stringstream message;
        message << iWeight;
        infoPtr->errorMsg("Warning in SimpleTimeShower: reject denom for "
          "iWeight = ", message.str());
      }
      // Force reweighting factor > 0.
      double reWtFail = max(0.01, (1. - uVarFac[iWeight] * pAccept / enhance)
        / denom);
      infoPtr->reWeight(iWeight, reWtFail);
    }
  }
}

//==========================================================================

// Rescatter: If a dipole stretches between two different systems, those
//            systems will no longer locally conserve momentum. These
//            imbalances become problematic when ISR or primordial kT
//            is switched on as these steps involve Lorentz boosts.
//
//            'rescatterPropagateRecoil' tries to fix momentum in all
//            systems by propogating recoil momentum through all
//            intermediate systems. As the momentum transfer is already
//            defined, this can lead to internal lines gaining a
//            virtuality.

// Useful definitions for a pair of integers and a vector of pairs
typedef pair < int, int >  pairInt;
typedef vector < pairInt > vectorPairInt;

//--------------------------------------------------------------------------

// findParentSystems
//  Utility routine to find all parent systems of a given system
//  Returns a vector of pairs of integers with:
//   a) The system index, including the starting system (negative
//      if (b) points to a parent system, positive if (b) points
//      to a daughter system
//   b) The event record index that is the path out of the system
//      (if forwards == false, this is an incoming parton to the
//      system, and is +ve if side A or -ve if side B,
//      if forwards == true, this is an outgoing parton from the
//      system).
//  Returns as empty vector on failure
//  Note: this assumes single rescattering only and therefore only
//        one possible parent system

inline vectorPairInt findParentSystems(const int sys,
  Event& event, PartonSystems* partonSystemsPtr, bool forwards) {

  vectorPairInt parentSystems;
  parentSystems.reserve(10);

  int iSysCur = sys;
  while (true) {
    // Get two incoming partons
    int iInA = partonSystemsPtr->getInA(iSysCur);
    int iInB = partonSystemsPtr->getInB(iSysCur);

    // Check if either of these links to another system
    int iIn = 0;
    if (event[iInA].isRescatteredIncoming()) iIn =  iInA;
    if (event[iInB].isRescatteredIncoming()) iIn = -iInB;

    // Save the current system to the vector
    parentSystems.push_back( pairInt(-iSysCur, iIn) );
    if (iIn == 0) break;

    int iInAbs  = abs(iIn);
    int iMother = event[iInAbs].mother1();
    iSysCur     = partonSystemsPtr->getSystemOf(iMother);
    if (iSysCur == -1) {
      parentSystems.clear();
      break;
    }
  } // while (true)

  // If forwards is set, change all event record indices to go to daughter
  // systems rather than parent systems
  if (forwards) {
    vectorPairInt::reverse_iterator rit;
    for (rit = parentSystems.rbegin(); rit < (parentSystems.rend() - 1);
         ++rit) {
      pairInt &cur  = *rit;
      pairInt &next = *(rit + 1);
      cur.first     = -cur.first;
      cur.second    = (next.second < 0) ? -event[abs(next.second)].mother1() :
                                           event[abs(next.second)].mother1();
    }
  }

  return parentSystems;
}

//--------------------------------------------------------------------------

// rescatterPropagateRecoil
//  Fix up momentum in all intermediate systems when radiator and recoiler
//  systems are different. The strategy is to look at all parent systems
//  from the radiator system and the recoiler system and find where they
//  intersect.

bool SimpleTimeShower::rescatterPropagateRecoil( Event& event, Vec4& pNew) {

  // Some useful variables for later
  int  iRadBef    = dipSel->iRadiator;
       iSysSel    = dipSel->system;
  int  iSysSelRec = dipSel->systemRec;
  Vec4 pImbal     = pNew - event[iRadBef].p();

  // Store changes locally at first in case we veto the branching
  // eventMod stores index into the event record and new 4-vector
  vector < pair < int, Vec4 > > eventMod;
  eventMod.reserve(10);
  // systemMod stores system index (iSys) and system-parton index (iMem)
  //   iMem >=  0 - index into outgoing partons (iOut)
  //   iMem == -1 - incoming A
  //   iMem == -2 - incoming B
  vectorPairInt systemMod;
  systemMod.reserve(10);

  // Find all parent systems from radiating and recoiling systems
  vectorPairInt radParent = findParentSystems(iSysSel, event,
                                              partonSystemsPtr, false);
  vectorPairInt recParent = findParentSystems(iSysSelRec, event,
                                              partonSystemsPtr, true);
  if (radParent.size() == 0 || recParent.size() == 0) {
    // This should never happen
    infoPtr->errorMsg("Error in SimpleTimeShower::rescatterPropagateRecoil: "
      "couldn't find parent system; branching vetoed");
    return false;
  }
  // Find the system that connects radiating and recoiling system
  bool foundPath = false;
  unsigned int iRadP = 0;
  unsigned int iRecP = 0;
  for (iRadP = 0; iRadP < radParent.size(); iRadP++) {
    for (iRecP = 0; iRecP < recParent.size(); iRecP++)
      if (abs(radParent[iRadP].first) == abs(recParent[iRecP].first)) {
        foundPath = true;
        break;
      }
    if (foundPath) break;
  }
  if (!foundPath) {
    // Can fail e.g. for QED dipoles where there is no connection
    // between radiator and recoiler systems
    infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterPropagateRecoil: "
      "couldn't find recoil path; branching vetoed");
    return false;
  }

  // Join together to form complete path from radiating system
  // to recoiling system
  vectorPairInt path;
  if (radParent.size() > 1)
    path.assign(radParent.begin(), radParent.begin() + iRadP);
  if (recParent.size() > 1)
    path.insert(path.end(), recParent.rend() - iRecP - 1,
                recParent.rend() - 1);

  // Follow the path fixing up momenta as we go
  for (unsigned int i = 0; i < path.size(); i++) {
    // Line out of the current system
    bool isIncoming  = (path[i].first < 0) ? true : false;
    int  iSysCur     = abs(path[i].first);
    bool isIncomingA = (path[i].second > 0) ? true : false;
    int  iLink       = abs(path[i].second);

    int iMemCur;
    if (isIncoming) iMemCur = (isIncomingA) ? -1 : -2;
    else {
      iMemCur = -1;
      for (int j = 0; j < partonSystemsPtr->sizeOut(iSysCur); j++)
        if (partonSystemsPtr->getOut(iSysCur, j) == iLink) {
          iMemCur = j;
          break;
        }
      if (iMemCur == -1) {
        // This should never happen
        infoPtr->errorMsg("Error in SimpleTimeShower::rescatterPropagateRecoil"
          ": couldn't find parton system; branching vetoed");
        return false;
      }
    }

    Vec4 pMod = (isIncoming) ? event[iLink].p() + pImbal :
                               event[iLink].p() - pImbal;
    eventMod.push_back(pair <int, Vec4> (iLink, pMod));
    systemMod.push_back(pairInt(iSysCur, iMemCur));

    // Calculate sHat of iSysCur
    int  iInCurA = partonSystemsPtr->getInA(iSysCur);
    int  iInCurB = partonSystemsPtr->getInB(iSysCur);
    Vec4 pTotCur = event[iInCurA].p() + event[iInCurB].p();

    // If iMemCur is -1 or -2, then we must have changed the sHat of iSysCur
    if (iMemCur < 0) pTotCur += (isIncoming) ? pImbal : -pImbal;
    double sHatCur = pTotCur.m2Calc();

    // The fixed-up incoming and outgoing partons should not have
    // too large a virtuality in relation to the system mass-square.
    if (abs(pMod.m2Calc()) > MAXVIRTUALITYFRACTION * sHatCur) {
      infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterPropagateRecoil"
        ": virtuality much larger than sHat; branching vetoed");
      return false;
    }

    // Outgoing ones should also not have too large negative energy
    // in the rest frame of the system.
    if (!isIncoming && pMod * pTotCur < -MAXNEGENERGYFRACTION * sHatCur) {
      infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterPropagateRecoil"
        ": rest frame energy too negative; branching vetoed");
      return false;
    }

    // Veto negative sHat.
    if (sHatCur < 0.0) {
      infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterPropagateRecoil"
        ": sHat became negative; branching vetoed");
      return false;
    }

    // Line into the new current system
    iLink   = (isIncoming) ? event[iLink].mother1()  :
                             event[iLink].daughter1();
    iSysCur = partonSystemsPtr->getSystemOf(iLink, true);

    if (!isIncoming) iMemCur = (isIncomingA) ? -1 : -2;
    else {
      iMemCur = -1;
      for (int j = 0; j < partonSystemsPtr->sizeOut(iSysCur); j++)
        if (partonSystemsPtr->getOut(iSysCur, j) == iLink) {
          iMemCur = j;
          break;
        }
      if (iMemCur == -1) {
        // This should never happen
        infoPtr->errorMsg("Error in SimpleTimeShower::rescatterPropagateRecoil"
          ": couldn't find parton system; branching vetoed");
        return false;
      }
    }

    pMod = (isIncoming) ? event[iLink].p() + pImbal :
                          event[iLink].p() - pImbal;
    eventMod.push_back(pair <int, Vec4> (iLink, pMod));
    systemMod.push_back(pairInt(iSysCur, iMemCur));

    // Calculate sHat of iSysCur
    iInCurA = partonSystemsPtr->getInA(iSysCur);
    iInCurB = partonSystemsPtr->getInB(iSysCur);
    pTotCur = event[iInCurA].p() + event[iInCurB].p();

    // If iMemCur is -1 or -2, then we must have changed the sHat of iSysCur
    if (iMemCur < 0) pTotCur += (isIncoming) ? pImbal : -pImbal;
    sHatCur = pTotCur.m2Calc();

    // The fixed-up incoming and outgoing partons should not have
    // too large a virtuality in relation to the system mass-square.
    if (abs(pMod.m2Calc()) > MAXVIRTUALITYFRACTION * sHatCur) {
      infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterPropagateRecoil"
        ": virtuality much larger than sHat; branching vetoed");
      return false;
    }

    // Outgoing ones should also not have too large negative energy
    // in the rest frame of the system.
    if (!isIncoming && pMod * pTotCur < -MAXNEGENERGYFRACTION * sHatCur) {
      infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterPropagateRecoil"
        ": rest frame energy too negative; branching vetoed");
      return false;
    }

    // Veto negative sHat
    if (sHatCur < 0.0) {
      infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterPropagateRecoil"
        ": sHat became negative; branching vetoed");
      return false;
    }

    // Do negative energy veto
    if (VETONEGENERGY && pMod.e() < 0.0) {
      infoPtr->errorMsg("Warning in SimpleTimeShower::rescatterPropagateRecoil"
        ": energy became negative; branching vetoed");
      return false;
    }

  } // for (unsigned int i = 0; i < path.size(); i++)

  // If no vetos by this point, apply the changes to the event record
  // An incoming particle with changed momentum is given status code -54,
  // an outgoing particle with changed momentum is given status code -55
  for (unsigned int i = 0; i < eventMod.size(); i++) {
    int idx    = eventMod[i].first;
    Vec4 &pMod = eventMod[i].second;
    int iSys   = systemMod[i].first;
    int iMem   = systemMod[i].second;

    // If incoming to a process then set the copy to be the mother
    if (event[idx].isRescatteredIncoming()) {
      int mother1 = event[idx].mother1();
      idx = event.copy(idx, -54);
      event[mother1].daughters(idx, idx);

      // Update beam information if necessary
      double eCM = sqrt(m2( beamAPtr->p(), beamBPtr->p()));
      if        (iMem == -1) {
        partonSystemsPtr->setInA(iSys, idx);
        (*beamAPtr)[iSys].x((pMod.e() + pMod.pz()) / eCM);
        (*beamAPtr)[iSys].m(pMod.mCalc());
        (*beamAPtr)[iSys].p(pMod);
        (*beamAPtr)[iSys].iPos(idx);
      } else if (iMem == -2) {
        partonSystemsPtr->setInB(iSys, idx);
        (*beamBPtr)[iSys].x((pMod.e() - pMod.pz()) / eCM);
        (*beamBPtr)[iSys].m(pMod.mCalc());
        (*beamBPtr)[iSys].p(pMod);
        (*beamBPtr)[iSys].iPos(idx);
      } else {
        // This should never happen
        infoPtr->errorMsg("Error in SimpleTimeShower::rescatterPropagateRecoil"
        ": internal bookeeping error");
      }

    // Otherwise set the new event record entry to be the daughter
    } else {
      int daughter1 = event[idx].daughter1();
      idx = event.copy(idx, 55);
      event[idx].statusNeg();
      event[daughter1].mothers(idx, idx);

      partonSystemsPtr->setOut(iSys, iMem, idx);
    }

    event[idx].p( eventMod[i].second );
    event[idx].m( event[idx].mCalc() );
  }

  return true;
}


//--------------------------------------------------------------------------

// Find class of QCD ME correction.
// MEtype classification follow codes in Norrbin article,
// additionally -1 = try to find type, 0 = no ME corrections.

void SimpleTimeShower::findMEtype( Event& event, TimeDipoleEnd& dip) {

  // Initial value. Mark if no ME corrections to be applied.
  bool setME = doMEcorrections;
  int iMother  = event[dip.iRadiator].mother1();
  int iMother2 = event[dip.iRadiator].mother2();

  // Allow ME corrections for Hidden Valley pair in 2 -> 2.
  if (dip.isHiddenValley && event[dip.iRecoiler].id()
    == -event[dip.iRadiator].id());

  // Allow ME corrections for all weak branchings.
  else if (dip.weakType != 0);

  // Else optionally no ME corrections in 2 -> n processes.
  else if (!doMEextended) {
    if (iMother2 != iMother && iMother2 != 0) setME = false;
    if (event[dip.iRecoiler].mother1() != iMother)  setME = false;
    if (event[dip.iRecoiler].mother2() != iMother2) setME = false;
  }

  // Optionally no ME corrections for recoiler in initial state.
  if (event[dip.iRecoiler].status() < 0) setME = doMEextended;

  // No ME corrections for recoiler not in same system
  if (dip.system != dip.systemRec) setME = false;

  // Done if no ME to be set.
  if (!setME) {
    dip.MEtype = 0;
    return;
  }

  // Pair "rare" particles, if possible.
  if (dip.iMEpartner < 0) {
    int idAbs1   = event[dip.iRadiator].idAbs();
    int idAbs2   = event[dip.iRecoiler].idAbs();
    bool isRare1 = (idAbs1 > 5 && idAbs1 < 11) || (idAbs1 > 16 && idAbs1 < 21)
                 || idAbs1 > 22;
    bool isRare2 = (idAbs2 > 5 && idAbs2 < 11) || (idAbs2 > 16 && idAbs2 < 21)
                 || idAbs2 > 22;
    if (isRare1 && !isRare2) {
      vector<int> iSis = event[dip.iRadiator].sisterList();
      // Prio on particle-(anti)particle pairs, else other rare.
      for (int iS = 0; iS < int(iSis.size()); ++iS) {
        idAbs2   = event[iSis[iS]].idAbs();
        isRare2 = (idAbs2 > 5 && idAbs2 < 11) || (idAbs2 > 16 && idAbs2 < 21)
                || idAbs2 > 22;
        if (idAbs2 == idAbs1) dip.iMEpartner = iSis[iS];
        if (isRare2 && dip.iMEpartner < 0) dip.iMEpartner = iSis[iS];
      }
    }
  }

  // If no ME partner set, assume it is the recoiler.
  if (dip.iMEpartner < 0) dip.iMEpartner = dip.iRecoiler;

  // If ME already set, assume everything is in order.
  if (dip.MEtype != -1) return;

  // Now begin processing of colour dipole, including Hidden Valley.
  if (dip.colType != 0 || dip.colvType != 0) {
    bool isHiddenColour = (dip.colvType != 0);

    // Find daughter types (may or may not be used later on).
    int idDau1      = event[dip.iRadiator].id();
    int idDau2      = event[dip.iMEpartner].id();
    int dau1Type    = findMEparticle(idDau1, isHiddenColour);
    int dau2Type    = findMEparticle(idDau2, isHiddenColour);
    int minDauType  = min(dau1Type, dau2Type);
    int maxDauType  = max(dau1Type, dau2Type);

    // Reorder dipole ends in kinematics. Split ME expression in two sides.
    dip.MEorder     = (dau2Type >= dau1Type);
    dip.MEsplit     = (maxDauType <= 6);
    dip.MEgluinoRec = false;

    // If type already set (or set not to have) then done.
    if (minDauType == 0 && dip.MEtype < 0) dip.MEtype = 0;
    if (dip.MEtype >= 0) return;
    dip.MEtype = 0;

    // For H -> gg -> ggg we found that DGLAP kernels do better than eikonal.
    if (dau1Type == 4 && dau2Type == 4) return;

    // Find mother type.
    int idMother = 0;
    if ( event[dip.iRecoiler].mother1() == iMother && iMother >= 0
      && (iMother2 == 0 || iMother2 == iMother) )
      idMother = event[iMother].id();
    int motherType = (idMother != 0)
      ? findMEparticle(idMother, isHiddenColour) : 0;

    // When a mother is not known then use colour and spin content to guess.
    if (motherType == 0) {
      int col1  = event[dip.iRadiator].col();
      int acol1 = event[dip.iRadiator].acol();
      int col2  = event[dip.iMEpartner].col();
      int acol2 = event[dip.iMEpartner].acol();
      // spinT = 0/1 = integer or half-integer.
      int spinT = ( event[dip.iRadiator].spinType()
                + event[dip.iMEpartner].spinType() )%2;
      // Colour singlet mother.
      if ( col1 == acol2 && acol1 == col2 )
        motherType = (spinT == 0) ? 7 : 9;
      // Colour octet mother.
      else if ( (col1 == acol2 && acol1 != 0 && col2 != 0)
        || (acol1 == col2 && col1 != 0 && acol2 != 0) )
        motherType = (spinT == 0) ? 4 : 5;
      // Colour triplet mother.
      else if ( (col1 == acol2 && acol1 != col2)
        || (acol1 == col2 && col1 != acol2) )
        motherType = (spinT == 0) ? 2 : 1;
      // If no colours are matched then cannot have common mother, so done.
      else return;
    }

    // Now start from default, which is eikonal ME corrections,
    // and try to find matching ME cases below.
    int MEkind = 0;
    int MEcombi = 4;
    dip.MEmix = 0.5;

    // Hidden Valley with massive gamma_v covered by two special cases.
    if (isHiddenColour && brokenHVsym) {
      MEkind = (dau2Type == 0 || dau2Type > 6) ? 30 : 31;
      dip.MEtype = 5 * MEkind + 1;
      return;
    }

    // Triplet recoiling against gluino needs enhanced radiation
    // to match to matrix elements.
    dip.MEgluinoRec = (dau1Type >= 1 && dau1Type <= 3 && dau2Type == 5);

    // Vector/axial vector -> q + qbar.
    if (minDauType == 1 && maxDauType == 1 &&
      (motherType == 4 || motherType == 7) ) {
      MEkind = 2;
      if (idMother == 21 || idMother == 22 || motherType == 4) MEcombi = 1;
      else if (idMother == 23 || idDau1 + idDau2 == 0) {
        MEcombi = 3;
        dip.MEmix = gammaZmix( event, iMother, dip.iRadiator, dip.iRecoiler );
      }
      else if (idMother == 24) MEcombi = 4;
    }
    // For chi -> chi q qbar, use V/A -> q qbar as first approximation.
    else if (minDauType == 1 && maxDauType == 1 && motherType == 9)
      MEkind = 2;

    // q -> q + V.
    else if (minDauType == 1 && maxDauType == 7 && motherType == 1) {
      MEkind = 3;
      if (idDau1 == 22 || idDau2 == 22) MEcombi = 1;

    // Scalar/pseudoscalar -> q + qbar; q -> q + S.
    } else if (minDauType == 1 && maxDauType == 1 && motherType == 8) {
      MEkind = 4;
      if (idMother == 25 || idMother == 35 || idMother == 37) MEcombi = 1;
      else if (idMother == 36) MEcombi = 2;
    }
    else if (minDauType == 1 && maxDauType == 8 && motherType == 1)
      MEkind = 5;

    // V -> ~q + ~qbar; ~q -> ~q + V; S -> ~q + ~qbar; ~q -> ~q + S.
    else if (minDauType == 2 && maxDauType == 2 && (motherType == 4
      || motherType == 7) ) MEkind = 6;
    else if (minDauType == 2 && (maxDauType == 4 || maxDauType == 7)
      && motherType == 2) MEkind = 7;
    else if (minDauType == 2 && maxDauType == 2 && motherType == 8)
      MEkind = 8;
    else if (minDauType == 2 && maxDauType == 8 && motherType == 2)
      MEkind = 9;

    // chi -> q + ~qbar; ~q -> q + chi; q -> ~q + chi.
    else if (minDauType == 1 && maxDauType == 2 && motherType == 9)
      MEkind = 10;
    else if (minDauType == 1 && maxDauType == 9 && motherType == 2)
      MEkind = 11;
    else if (minDauType == 2 && maxDauType == 9 && motherType == 1)
      MEkind = 12;

    // ~g -> q + ~qbar; ~q -> q + ~g; q -> ~q + ~g.
    else if (minDauType == 1 && maxDauType == 2 && motherType == 5)
      MEkind = 13;
    else if (minDauType == 1 && maxDauType == 5 && motherType == 2)
      MEkind = 14;
    else if (minDauType == 2 && maxDauType == 5 && motherType == 1)
      MEkind = 15;

    // In cases where coloured spin 1 particle involved use spin 0.
    // V_coloured -> q + l.
    else if (minDauType == 1 && maxDauType == 9 && motherType == 3)
      MEkind = 11;
    // q -> V_coloured + l;
    else if (minDauType == 3 && maxDauType == 9 && motherType == 1)
      MEkind = 12;

    // g (+V, S) -> ~g + ~g (eikonal approximation).
    else if (minDauType == 5 && maxDauType == 5) MEkind = 16;

    // Save ME type and gamma_5 admixture.
    dip.MEtype = 5 * MEkind + MEcombi;

  // Now begin processing of charge dipole - still primitive.
  } else if (dip.chgType != 0) {

    // Set defaults for QED case; then possibly done.
    dip.MEorder = true;
    dip.MEsplit = true;
    if (dip.MEtype >= 0) return;

    // So far only ME corrections for q qbar or l lbar.
    int idDau1 = event[dip.iRadiator].id();
    int idDau2 = event[dip.iMEpartner].id();
    if (abs(idDau1) < 9 && abs(idDau2) < 9 && idDau1 * idDau2 < 0) ;
    else if (abs(idDau1) > 10 && abs(idDau1) < 19 && abs(idDau2) > 10
      && abs(idDau2) < 19 && idDau1 * idDau2 < 0) ;
    else { dip.MEtype = 0; return; }

    // Distinguish charge sum != 0 or = 0; in latter assume vector source.
    dip.MEtype = 101;
    if (idDau1 + idDau2 == 0) dip.MEtype = 102;
    dip.MEmix = 1.;
  }

  // Identify 2 -> 2 processes for weak corrections.
  else if (dip.weakType == 1) {
    if (event[dip.iRadiator].id() == -event[dip.iRecoiler].id()
      || event[event[dip.iRadiator].mother1()].idAbs() == 24
        || infoPtr->nFinal() != 2) dip.MEtype = 200;
    else if (event[dip.iRadiator].idAbs() == 21
      || event[dip.iRecoiler].idAbs() == 21) dip.MEtype = 201;
    else if (event[dip.iRadiator].id() == event[dip.iRecoiler].id())
      dip.MEtype = 202;
    else dip.MEtype = 203;
  } else if (dip.weakType == 2) {
    if (event[dip.iRadiator].id() == -event[dip.iRecoiler].id()
      || event[event[dip.iRadiator].mother1()].idAbs() == 24) dip.MEtype = 205;
    else if (event[dip.iRadiator].idAbs() == 21
      || event[dip.iRecoiler].idAbs() == 21) dip.MEtype = 206;
    else if (event[dip.iRadiator].id() == event[dip.iRecoiler].id())
      dip.MEtype = 207;
    else dip.MEtype = 208;
  }

}

//--------------------------------------------------------------------------

// Find type of particle for ME type: 0 = unknown, 1 = quark, 2 = squark,
// 3 = spare triplet, 4 = gluon, 5 = gluino, 6 = spare octet,
// 7 = vector boson, 8 = colourless scalar, 9 = colourless spin 1/2.

int SimpleTimeShower::findMEparticle( int id, bool isHiddenColour) {

  // find colour and spin of particle.
  int type = 0;
  int colType = abs(particleDataPtr->colType(id));
  int spinType = particleDataPtr->spinType(id);

  // For hidden valley particle treat HV colour as normal one.
  // Note: no need to assign gv/gammav since not in ME.
  if (isHiddenColour) {
    colType = 0;
    int idAbs = abs(id);
    if (  (idAbs > 4900000 && idAbs < 4900007)
       || (idAbs > 4900010 && idAbs < 4900017)
       || (idAbs > 4900100 && idAbs < 4900109) ) colType = 1;
  }

  // Find particle type from colour and spin.
  if      (colType == 1 && spinType == 2) type = 1;
  else if (colType == 1 && spinType == 1) type = 2;
  else if (colType == 1)                  type = 3;
  else if (colType == 2 && spinType == 3) type = 4;
  else if (colType == 2 && spinType == 2) type = 5;
  else if (colType == 2)                  type = 6;
  else if (colType == 0 && spinType == 3) type = 7;
  else if (colType == 0 && spinType == 1) type = 8;
  else if (colType == 0 && spinType == 2) type = 9;

  // Done.
  return type;

}

//--------------------------------------------------------------------------

// Find mixture of V and A in gamma/Z: energy- and flavour-dependent.

double SimpleTimeShower::gammaZmix( Event& event, int iRes, int iDau1,
  int iDau2) {

  // Try to identify initial flavours; use e+e- as default.
  int idIn1 = -11;
  int idIn2 = 11;
  int iIn1  = (iRes >= 0) ? event[iRes].mother1() : -1;
  int iIn2  = (iRes >= 0) ? event[iRes].mother2() : -1;
  if (iIn1 > 0 && iIn2 <= 0 && event[iDau1].mother2() > 0)
    iIn2 = event[event[iDau1].mother2()].mother1();
  if (iIn1 >=0) idIn1 = event[iIn1].id();
  if (iIn2 >=0) idIn2 = event[iIn2].id();

  // In processes f + g/gamma -> f + Z only need find one fermion.
  if (idIn1 == 21 || idIn1 == 22) idIn1 = -idIn2;
  if (idIn2 == 21 || idIn2 == 22) idIn2 = -idIn1;

  // Initial flavours and couplings; return if don't make sense.
  if (idIn1 + idIn2 != 0 ) return 0.5;
  int idInAbs = abs(idIn1);
  if (idInAbs == 0 || idInAbs > 18 ) return 0.5;
  double ei = coupSMPtr->ef(idInAbs);
  double vi = coupSMPtr->vf(idInAbs);
  double ai = coupSMPtr->af(idInAbs);

  // Final flavours and couplings; return if don't make sense.
  if (event[iDau1].id() + event[iDau2].id() != 0) return 0.5;
  int idOutAbs = abs(event[iDau1].id());
  if (idOutAbs == 0 || idOutAbs >18 ) return 0.5;
  double ef = coupSMPtr->ef(idOutAbs);
  double vf = coupSMPtr->vf(idOutAbs);
  double af = coupSMPtr->af(idOutAbs);

  // Calculate prefactors for interference and resonance part.
  Vec4 psum = event[iDau1].p() + event[iDau2].p();
  double sH = psum.m2Calc();
  double intNorm = 2. * thetaWRat * sH * (sH - mZ*mZ)
    / ( pow2(sH - mZ*mZ) + pow2(sH * gammaZ / mZ) );
  double resNorm = pow2(thetaWRat * sH)
    / ( pow2(sH - mZ*mZ) + pow2(sH * gammaZ / mZ) );

  // Calculate vector and axial expressions and find mix.
  double vect = ei*ei * ef*ef + ei*vi * intNorm * ef*vf
    + (vi*vi + ai*ai) * resNorm * vf*vf;
  double axiv = (vi*vi + ai*ai) * resNorm * af*af;
  return vect / (vect + axiv);
}

//--------------------------------------------------------------------------

// Set up to calculate QCD ME correction with calcMEcorr.
// Normally for primary particles, but also from g/gamma -> f fbar.

double SimpleTimeShower::findMEcorr(TimeDipoleEnd* dip, Particle& rad,
  Particle& partner, Particle& emt, bool cutEdge) {

  // Initial values and matrix element kind.
  double wtME    = 1.;
  double wtPS    = 1.;
  int    MEkind  = dip->MEtype / 5;
  int    MEcombi = dip->MEtype % 5;

  // Construct ME variables.
  Vec4   sum     = rad.p() + partner.p() + emt.p();
  double eCMME   = sum.mCalc();
  double x1      = 2. * (sum * rad.p()) / pow2(eCMME);
  double x2      = 2. * (sum * partner.p()) / pow2(eCMME);
  double r1      = rad.m() / eCMME;
  double r2      = partner.m() / eCMME;
  double r3      = 0.;

  // Evaluate kinematics for Hidden Valley with massive gamma_v.
  double gammavCorr = 1.;
  if (dip->colvType != 0 && brokenHVsym) {
    r3              = emt.m() / eCMME;
    double x3Tmp    = 2. - x1 - x2;
    gammavCorr      = x3Tmp / (x3Tmp - kRad * (x1 + x3Tmp));
    // For Q_v Qbar_v pair correct kinematics to common average mass.
    if (MEkind == 31) {
      double m2Pair = (rad.p() + partner.p()).m2Calc();
      double m2Avg  = 0.5 * (rad.m2() + partner.m2())
                    - 0.25 * pow2(rad.m2() - partner.m2()) / m2Pair;
      r1            = sqrt(m2Avg) / eCMME;
      r2            = r1;
      double xShift = 0.5 * (x1 + x2) * (partner.m2() - rad.m2()) / m2Pair;
      x1           += xShift;
      x2           -= xShift;
    }
  }

  // Derived ME variables, suitably protected.
  double x1minus, x2minus, x3;
  if (cutEdge) {
    x1minus = max(XMARGIN, 1. + r1*r1 - r2*r2 - x1);
    x2minus = max(XMARGIN, 1. + r2*r2 - r1*r1 - x2) ;
    x3      = max(XMARGIN, 2. - x1 - x2);
  } else {
    x1minus = max(XMARGIN*XMARGIN, 1. + r1*r1 - r2*r2 - x1);
    x2minus = max(XMARGIN*XMARGIN, 1. + r2*r2 - r1*r1 - x2) ;
    x3      = max(XMARGIN*XMARGIN, 2. - x1 - x2);
  }

  // Begin processing of QCD dipoles.
  if (dip->colType !=0 || dip->colvType != 0) {

    // Evaluate normal ME, for proper order of particles.
    if (dip->MEorder) wtME = calcMEcorr(MEkind, MEcombi, dip->MEmix,
      x1, x2, r1, r2, r3, cutEdge);
    else wtME = calcMEcorr(MEkind, MEcombi, dip->MEmix,
      x2, x1, r2, r1, r3, cutEdge);

    // Split up total ME when two radiating particles.
    if (dip->MEsplit) wtME = wtME * x1minus / x3;

    // Evaluate shower rate to be compared with.
    wtPS = 2. / ( x3 * x2minus );
    if (dip->MEgluinoRec) wtPS *= 9./4.;
    if (dip->colvType != 0 && brokenHVsym) wtPS *= gammavCorr;

  // For generic charge combination currently only massless expression.
  // (Masses included only to respect phase space boundaries.)
  } else if (dip->chgType !=0 && dip->MEtype == 101) {
    double chg1 = particleDataPtr->charge(rad.id());
    double chg2 = particleDataPtr->charge(partner.id());
    wtME = (x1*x1 + x2*x2) * pow2( chg1 * x1minus / x3
      - chg2 * x2minus / x3 );
    wtPS = 2. * ( chg1*chg1 * x1minus / x3 + chg2*chg2 * x2minus / x3 );

  // For flavour neutral system assume vector source and include masses.
  } else if (dip->chgType !=0 && dip->MEtype == 102) {
    wtME = calcMEcorr(2, 1, dip->MEmix, x1, x2, r1, r2, 0., cutEdge)
      * x1minus / x3;
    wtPS = 2. / ( x3 * x2minus );
  }

  // Weak W and Z emissions, currently using same matrix element.
  // The s-channel corrections are handled with simple MEs.
  else if (dip->MEtype == 200 || dip->MEtype == 205) {
    r3   = emt.m() / eCMME;
    wtME = calcMEcorr(32, 1, dip->MEmix, x1, x2, r1, r2, r3, cutEdge)
      * x1minus / x3;
    wtPS = 8. / (x3 * x2minus);
    wtPS *= x3 / (x3 - kRad * (x1 + x3));
  }
  // The t-channel corrections are handled separately in findMEweak.
  else if (dip->MEtype == 201 || dip->MEtype == 202
       ||  dip->MEtype == 203 || dip->MEtype == 205
       ||  dip->MEtype == 206 || dip->MEtype == 207) return 1.;

  // Return ratio of actual ME to assumed PS rate of emission.
  if (wtME > wtPS) infoPtr->errorMsg("Warning in SimpleTimeShower::findMEcorr:"
    " ME weight above PS one");
  return wtME / wtPS;

}

//--------------------------------------------------------------------------

// Matrix elements for gluon (or photon) emission from
// a two-body state; to be used by the parton shower routine.
// Here x_i = 2 E_i/E_cm, r_i = m_i/E_cm and
// 1/sigma_0 d(sigma)/d(x_1)d(x_2) = (alpha-strong/2 pi) * C_F * (this),
// i.e. normalization is such that one recovers the familiar
// (x_1^2 + x_2^2)/((1-x_1)*(1-x_2)) for the massless case.
// Coupling structure:
// kind =  1 : eikonal soft-gluon expression (spin-independent)
//      =  2 : V -> q qbar (V = vector/axial vector colour singlet)
//      =  3 : q -> q V
//      =  4 : S -> q qbar (S = scalar/pseudoscalar colour singlet)
//      =  5 : q -> q S
//      =  6 : V -> ~q ~qbar (~q = squark)
//      =  7 : ~q -> ~q V
//      =  8 : S -> ~q ~qbar
//      =  9 : ~q -> ~q S
//      = 10 : chi -> q ~qbar (chi = neutralino/chargino)
//      = 11 : ~q -> q chi
//      = 12 : q -> ~q chi
//      = 13 : ~g -> q ~qbar
//      = 14 : ~q -> q ~g
//      = 15 : q -> ~q ~g
//      = 16 : (9/4)*(eikonal) for gg -> ~g ~g
//      = 30 : Dv -> d qv     (Dv= hidden valley fermion, qv= valley scalar)
//      = 31 : S  -> Dv Dvbar (S=scalar colour singlet)
// Note that the order of the decay products is important.
// combi = 1 : pure non-gamma5, i.e. vector/scalar/...
//       = 2 : pure gamma5, i.e. axial vector/pseudoscalar/....
//       = 3 : mixture mix*(combi=1) + (1-mix)*(combi=2)
//       = 4 : mixture (combi=1) +- (combi=2)

double SimpleTimeShower::calcMEcorr( int kind, int combiIn, double mixIn,
  double x1, double x2, double r1, double r2, double r3, bool cutEdge) {

  // Frequent variable combinations.
  double x3     = 2. - x1 - x2;
  double x1s    = x1 * x1;
  double x2s    = x2 * x2;
  double x3s    = x3 * x3;
  double x1c    = x1 * x1s;
  double x2c    = x2 * x2s;
  double x3c    = x3 * x3s;
  double r1s    = r1 * r1;
  double r2s    = r2 * r2;
  double r1c    = r1 * r1s;
  double r2c    = r2 * r2s;
  double r1q    = r1s * r1s;
  double r2q    = r2s * r2s;
  double prop1  = 1. + r1s - r2s - x1;
  double prop2  = 1. + r2s - r1s - x2;
  double prop1s = prop1 * prop1;
  double prop2s = prop2 * prop2;
  double prop12 = prop1 * prop2;
  double prop13 = prop1 * x3;
  double prop23 = prop2 * x3;

  // Special case: Hidden-Valley massive photon.
  double r3s    = r3 * r3;
  double prop3  = r3s - x3;
  double prop3s = prop3 * prop3;
  if (kind == 30) prop13 = prop1 * prop3;

  // Check input values. Return zero outside allowed phase space.
  if (cutEdge) {
    if (x1 - 2.*r1 < XMARGIN || prop1 < XMARGIN) return 0.;
    if (x2 - 2.*r2 < XMARGIN || prop2 < XMARGIN) return 0.;
    // Limits not worked out for r3 > 0.
    if (kind != 30 && kind != 31) {
      if (x1 + x2 - 1. - pow2(r1+r2) < XMARGIN) return 0.;
      // Note: equivalent rewritten form 4. * ( (1. - x1) * (1. - x2)
      // * (1. - r1s - r2s - x3) + r1s * (1. - x2s - x3) + r2s
      // * (1. - x1s - x3) - pow2(r1s - r2s) ) gives about same result.
      if ( (x1s - 4.*r1s) * (x2s - 4.*r2s)
        - pow2( 2. * (1. - x1 - x2 + r1s + r2s) + x1*x2 )
        < XMARGIN * (XMARGINCOMB + r1 + r2) ) return 0.;
    }
  }

  // Initial values; phase space.
  int combi   = max(1, min(4, combiIn) );
  double mix  = max(0., min(1., mixIn) );
  bool isSet1 = false;
  bool isSet2 = false;
  bool isSet4 = false;
  double ps = sqrtpos( pow2(1. - r1*r1 - r2*r2) - pow2(2. * r1 * r2) );
  double rLO = 0., rFO = 0., rLO1 = 0., rFO1 = 0., rLO2 = 0.,
    rFO2 = 0., rLO4 = 0., rFO4 = 0.;
  double offset = 0;

  // Select which kind of ME to use.
  switch (kind) {

    // case 1 is equal to default, i.e. eikonal expression.

    // V -> q qbar (V = gamma*/Z0/W+-/...).
    case 2:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(2.-r1s-r1q+6.*r1*r2-r2s+2.*r1s*r2s-r2q)/2.;
        rFO1 = -(3.+6.*r1s+r1q-6.*r1*r2+6.*r1c*r2-2.*r2s-6.*r1s*r2s
        +6.*r1*r2c+r2q-3.*x1+6.*r1*r2*x1+2.*r2s*x1+x1s-2.*r1s*x1s
        +3.*r1s*x3+6.*r1*r2*x3-r2s*x3-2.*x1*x3-5.*r1s*x1*x3
        +r2s*x1*x3+x1s*x3-3.*x3s-3.*r1s*x3s+r2s*x3s
        +2.*x1*x3s+x3c-x2)
        /prop2s
        -2.*(-3.+r1s-6.*r1*r2+6.*r1c*r2+3.*r2s-4.*r1s*r2s
        +6.*r1*r2c+2.*x1+3.*r1s*x1+r2s*x1-x1s-r1s*x1s
        -r2s*x1s+4.*x3+2.*r1s*x3+3.*r1*r2*x3-r2s*x3-3.*x1*x3
        -2.*r1s*x1*x3+x1s*x3-x3s-r1s*x3s+r1*r2*x3s+x1*x3s)
        /prop12
        -(-1.+2.*r1s+r1q+6.*r1*r2+6.*r1c*r2-2.*r2s-6.*r1s*r2s
        +6.*r1*r2c+r2q-x1-2.*r1s*x1-6.*r1*r2*x1+8.*r2s*x1+x1s
        -2.*r2s*x1s-r1s*x3+r2s*x3-r1s*x1*x3+r2s*x1*x3+x1s*x3+x2)
        /prop1s;
        rFO1 = rFO1/2.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(2.-r1s-r1q-6.*r1*r2-r2s+2.*r1s*r2s-r2q)/2.;
        rFO2 = -(3.+6.*r1s+r1q+6.*r1*r2-6.*r1c*r2-2.*r2s-6.*r1s*r2s
        -6.*r1*r2c+r2q-3.*x1-6.*r1*r2*x1+2.*r2s*x1+x1s-2.*r1s*x1s
        +3.*r1s*x3-6.*r1*r2*x3-r2s*x3-2.*x1*x3-5.*r1s*x1*x3
        +r2s*x1*x3+x1s*x3-3.*x3s-3.*r1s*x3s+r2s*x3s+2.*x1*x3s+x3c-x2)
        /prop2s
        -2.*(-3+r1s+6.*r1*r2-6.*r1c*r2+3.*r2s-4.*r1s*r2s-6.*r1*r2c
        +2.*x1+3.*r1s*x1+r2s*x1-x1s-r1s*x1s-r2s*x1s+4.*x3+2.*r1s*x3
        -3.*r1*r2*x3-r2s*x3-3.*x1*x3-2.*r1s*x1*x3+x1s*x3-x3s-r1s*x3s
        -r1*r2*x3s+x1*x3s)
        /prop12
        -(-1.+2.*r1s+r1q-6.*r1*r2-6.*r1c*r2-2.*r2s-6.*r1s*r2s
        -6.*r1*r2c+r2q-x1-2.*r1s*x1+6.*r1*r2*x1+8.*r2s*x1+x1s
        -2.*r2s*x1s-r1s*x3+r2s*x3-r1s*x1*x3+r2s*x1*x3+x1s*x3+x2)
        /prop1s;
        rFO2 = rFO2/2.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(2.-r1s-r1q-r2s+2.*r1s*r2s-r2q)/2.;
        rFO4 = (1.-r1q+6.*r1s*r2s-r2q+x1+3.*r1s*x1-9.*r2s*x1-3.*x1s
        -r1s*x1s+3.*r2s*x1s+x1c-x2-r1s*x2+r2s*x2-r1s*x1*x2+r2s*x1*x2
        +x1s*x2)
        /prop1s
        -2.*(1.+r1s+r2s-4.*r1s*r2s+r1s*x1+2.*r2s*x1-x1s-r2s*x1s
        +2.*r1s*x2+r2s*x2-3.*x1*x2+x1s*x2-x2s-r1s*x2s+x1*x2s)
        /prop12
        +(1.-r1q+6.*r1s*r2s-r2q-x1+r1s*x1-r2s*x1+x2-9.*r1s*x2
        +3.*r2s*x2+r1s*x1*x2-r2s*x1*x2-3.*x2s+3.*r1s*x2s-r2s*x2s
        +x1*x2s+x2c)
        /prop2s;
        rFO4 = rFO4/2.;
        isSet4 = true;
      }
      break;

    // q -> q V.
    case 3:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-2.*r1s+r1q+r2s-6.*r1*r2s+r1s*r2s-2.*r2q);
        rFO1 = -2.*(-1.+r1-2.*r1s+2.*r1c-r1q+pow5(r1)-r2s+r1*r2s
        -5.*r1s*r2s+r1c*r2s-2.*r1*r2q+2.*x1-2.*r1*x1+2.*r1s*x1
        -2.*r1c*x1+2.*r2s*x1+5.*r1*r2s*x1+r1s*r2s*x1+2.*r2q*x1
        -x1s+r1*x1s-r2s*x1s+3.*x2+4.*r1s*x2+r1q*x2+2.*r2s*x2
        +2.*r1s*r2s*x2-4.*x1*x2-2.*r1s*x1*x2-r2s*x1*x2+x1s*x2
        -2.*x2s-2.*r1s*x2s+x1*x2s)
        /prop23
        +(2.*r2s+6.*r1*r2s-6.*r1s*r2s+6.*r1c*r2s+2.*r2q+6.*r1*r2q
        -r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2-3.*r2s*x2-6.*r1*r2s*x2
        +9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1+6.*r1*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s-
        2.*r1s*x1s+x1c+7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2+6.*r1*r2s*x2
        +r1s*r2s*x2-2.*r2q*x2-9.*x1*x2-3.*r1s*x1*x2+2.*r2s*x1*x2
        +2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s+x1*x2s)
        /x3s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-2.*r1s+r1q+r2s+6.*r1*r2s+r1s*r2s-2.*r2q);
        rFO2 = 2*(1.+r1+2.*r1s+2.*r1c+r1q+pow5(r1)+r2s+r1*r2s
        +5.*r1s*r2s+r1c*r2s-2.*r1*r2q-2.*x1-2.*r1*x1-2.*r1s*x1
        -2.*r1c*x1-2.*r2s*x1+5.*r1*r2s*x1-r1s*r2s*x1-2.*r2q*x1+x1s
        +r1*x1s+r2s*x1s-3.*x2-4.*r1s*x2-r1q*x2-2.*r2s*x2
        -2.*r1s*r2s*x2+4.*x1*x2+2.*r1s*x1*x2+r2s*x1*x2-x1s*x2
        +2.*x2s+2.*r1s*x2s-x1*x2s)
        /prop23
        +(2.*r2s-6.*r1*r2s-6.*r1s*r2s-6.*r1c*r2s+2.*r2q-6.*r1*r2q
        -r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2-3.*r2s*x2+6.*r1*r2s*x2
        +9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1-6.*r1*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s
        -2.*r1s*x1s+x1c+7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2-6.*r1*r2s*x2
        +r1s*r2s*x2-2.*r2q*x2-9.*x1*x2-3.*r1s*x1*x2+2.*r2s*x1*x2
        +2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s+x1*x2s)
        /x3s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-2.*r1s+r1q+r2s+r1s*r2s-2.*r2q);
        rFO4 = 2*(1.+2.*r1s+r1q+r2s+5.*r1s*r2s-2.*x1-2.*r1s*x1
        -2.*r2s*x1-r1s*r2s*x1-2.*r2q*x1+x1s+r2s*x1s-3.*x2-4.*r1s*x2
        -r1q*x2-2.*r2s*x2-2.*r1s*r2s*x2+4.*x1*x2+2.*r1s*x1*x2+r2s*x1*x2
        -x1s*x2+2.*x2s+2.*r1s*x2s-x1*x2s)
        /prop23
        +(2.*r2s-6.*r1s*r2s+2.*r2q-r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2
        -3.*r2s*x2+9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s-2.*r1s*x1s+x1c
        +7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2+r1s*r2s*x2-2.*r2q*x2-9.*x1*x2
        -3.*r1s*x1*x2+2.*r2s*x1*x2+2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s
        +x1*x2s)
        /x3s;
        isSet4 = true;
      }
      break;

    // S -> q qbar    (S = h0/H0/A0/H+-/...).
    case 4:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s-r2s-2.*r1*r2);
        rFO1 = -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -2.*(r1s+r1q-2.*r1c*r2+r2s-6.*r1s*r2s-2.*r1*r2c+r2q-r1s*x1
        +r1*r2*x1+2.*r2s*x1+2.*r1s*x2+r1*r2*x2-r2s*x2-x1*x2)
        /prop12
        -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1-r1s*x1
        +r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s-r2s+2.*r1*r2);
        rFO2 = -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1-2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +2.*(-r1s-r1q-2.*r1c*r2-r2s+6.*r1s*r2s-2.*r1*r2c-r2q+r1s*x1
        +r1*r2*x1-2.*r2s*x1-2.*r1s*x2+r1*r2*x2+r2s*x2+x1*x2)
        /prop12;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+3.*r2s*x1+x2
        +r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -2.*(r1s+r1q+r2s-6.*r1s*r2s+r2q-r1s*x1
        +2.*r2s*x1+2.*r1s*x2-r2s*x2-x1*x2)
        /prop12
        -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1
        +x2+3.*r1s*x2-r2s*x2-x1*x2)
        /prop2s;
        isSet4 = true;
      }
      break;

    // q -> q S.
    case 5:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = (4.-4.*r1s+4.*r2s-3.*x1-2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        -2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.-r1-5.*r1s-r1c+3.*r2s+r1*r2s-2.*x1-r1*x1
        +r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.+r1s-r2s-2.*r1);
        rFO2 = (4.-4.*r1s+4.*r2s-3.*x1+2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        +2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.+r1-5.*r1s+r1c+3.*r2s-r1*r2s-2.*x1+r1*x1
        +r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = (4.-4.*r1s+4.*r2s-3.*x1+r1s*x1-r2s*x1-5.*x2+r1s*x2
        -r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.-5.*r1s+3.*r2s-2.*x1+r1s*x1-4.*x2+2.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop23
        +(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet4 = true;
      }
      break;

    // V -> ~q ~qbar  (~q = squark).
    case 6:
      rLO1 = ps*(1.-2.*r1s+r1q-2.*r2s-2.*r1s*r2s+r2q);
      rFO1 = 2.*3.+(1.+r1s+r2s-x1)*(4.*r1s-x1s)
      /prop1s
      +2.*(-1.-3.*r1s-r2s+x1+x1s*0.5+x2-x1*x2*0.5)
      /prop1
      +(1.+r1s+r2s-x2)*(4.*r2s-x2s)
      /prop2s
      +2.*(-1.-r1s-3.*r2s+x1+x2-x1*x2*0.5+x2s*0.5)
      /prop2
      -(-4.*r1s-4.*r1q-4.*r2s-8.*r1s*r2s-4.*r2q+2.*x1+6.*r1s*x1
      +6.*r2s*x1-2.*x1s+2.*x2+6.*r1s*x2+6.*r2s*x2-4.*x1*x2
      -2.*r1s*x1*x2-2.*r2s*x1*x2+x1s*x2-2.*x2s+x1*x2s)
      /prop12;
      isSet1 = true;
      break;

    // ~q -> ~q V.
    case 7:
      rLO1 = ps*(1.-2.*r1s+r1q-2.*r2s-2.*r1s*r2s+r2q);
      rFO1 = 16.*r2s-8.*(4.*r2s+2.*r2s*x1+x2+r1s*x2+r2s*x2-x1*x2
      -2.*x2s)
      /(3.*prop2)
      +8.*(1.+r1s+r2s-x2)*(4.*r2s-x2s)
      /(3.*prop2s)
      +8.*(x1+x2)*(-1.-2.*r1s-r1q-2.*r2s+2.*r1s*r2s-r2q+2.*x1
      +2.*r1s*x1+2.*r2s*x1-x1s+2.*x2+2.*r1s*x2+2.*r2s*x2-2.*x1*x2-x2s)
      /(3.*x3s)
      +8.*(-1.-r1s+r2s-x1)*(2.*r2s*x1+x2+r1s*x2+r2s*x2-x1*x2-x2s)
      /(3.*prop2*x3)
      -8.*(1.+2.*r1s+r1q+2.*r2s-2.*r1s*r2s+r2q-2.*x1-2.*r1s*x1
      -4.*r2s*x1+x1s-3.*x2-3.*r1s*x2-3.*r2s*x2+3.*x1*x2+2.*x2s)
      /(3.*x3);
      rFO1 = 3.*rFO1/8.;
      isSet1 = true;
      break;

    // S -> ~q ~qbar.
    case 8:
      rLO1 = ps;
      rFO1 = (-1.-2.*r1s-r1q-2.*r2s+2.*r1s*r2s-r2q+2.*x1+2.*r1s*x1
      +2.*r2s*x1-x1s-r2s*x1s+2.*x2+2.*r1s*x2+2.*r2s*x2-3.*x1*x2
      -r1s*x1*x2-r2s*x1*x2+x1s*x2-x2s-r1s*x2s+x1*x2s)
      /(prop1s*prop2s);
      rFO1 = 2.*rFO1;
      isSet1 = true;
      break;

    // ~q -> ~q S.
    case 9:
      rLO1 = ps;
      rFO1 = (-1.-r1s-r2s+x2)
      /prop2s
      +(1.+r1s-r2s+x1)
      /prop23
      -(x1+x2)
      /x3s;
      isSet1 = true;
      break;

    // chi -> q ~qbar   (chi = neutralino/chargino).
    case 10:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = (2.*r1+x1)*(-1.-r1s-r2s+x1)
        /prop1s
        +2.*(-1.-r1s-2.*r1c-r2s-2.*r1*r2s+3.*x1*0.5+r1*x1
        -r1s*x1*0.5-r2s*x1*0.5+x2+r1*x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-2.*r1+r1s-r2s);
        rFO2 = (2.*r1-x1)*(1.+r1s+r2s-x1)
        /prop1s
        +2.*(-1.-r1s+2.*r1c-r2s+2.*r1*r2s+3.*x1*0.5-r1*x1
        -r1s*x1*0.5-r2s*x1*0.5+x2-r1*x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)/
        prop2s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = x1*(-1.-r1s-r2s+x1)
        /prop1s
        +2.*(-1.-r1s-r2s+3.*x1*0.5-r1s*x1*0.5-r2s*x1*0.5
        +x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet4 = true;
      }
      break;

    // ~q -> q chi.
    case 11:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-pow2(r1+r2));
        rFO1 = (1.+r1s+2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q-2.*r1*r2-2.*r1c*r2+2.*r1*r2c+r2q+x1+r1s*x1
        -2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-pow2(r1-r2));
        rFO2 = (1.+r1s-2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q+2.*r1*r2+2.*r1c*r2-2.*r1*r2c+r2q+x1+r1s*x1
        +2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = (1.+r1s+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1+x2
        +3.*r1s*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q+r2q+x1+r1s*x1-3.*r2s*x1
        +2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet4 = true;
      }
      break;

    // q -> ~q chi.
    case 12:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s+r2s+2.*r2);
        rFO1 = (2.*r2+x2)*(-1.-r1s-r2s+x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1-2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2-2.*r2*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s+r2+r1s*r2-r2s-r2c+x1+r2*x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s+r2s-2.*r2);
        rFO2 = (2.*r2-x2)*(1.+r1s+r2s-x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2+2.*r2*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s-r2-r1s*r2-r2s+r2c+x1-r2*x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s+r2s);
        rFO4 = x2*(-1.-r1s-r2s+x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+r2s*x1+x1s
        -3.*x2-r1s*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s-r2s+x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet4 = true;
      }
      break;

    // ~g -> q ~qbar.
    case 13:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = 4.*(2.*r1+x1)*(-1.-r1s-r2s+x1)
        /(3.*prop1s)
        -(-1.-r1s-2.*r1c-r2s-2.*r1*r2s+3.*x1*0.5+r1*x1-r1s*x1*0.5
        -r2s*x1*0.5+x2+r1*x2+r1s*x2-x1*x2*0.5)
        /(3.*prop12)
        +3.*(-1.+r1-r1s-r1c-r2s+r1*r2s+2.*x1+r2s*x1-x1s*0.5+x2+r1*x2
        +r1s*x2-x1*x2*0.5)
        /prop13
        +3.*(4.-4.*r1s+4.*r2s-3.*x1-2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        -2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -3.*(3.-r1-5.*r1s-r1c+3.*r2s+r1*r2s-2.*x1-r1*x1+r1s*x1
        -4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +4.*(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1-r2s*x1
        -3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /(3.*prop2s);
        rFO1 = 3.*rFO1/4.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.+r1s-r2s-2.*r1);
        rFO2 = 4.*(2.*r1-x1)*(1.+r1s+r2s-x1)
        /(3.*prop1s)
        +3.*(-1.-r1-r1s+r1c-r2s-r1*r2s+2.*x1+r2s*x1-x1s*0.5
        +x2-r1*x2+r1s*x2-x1*x2*0.5)
        /prop13
        +(2.+2.*r1s-4.*r1c+2.*r2s-4.*r1*r2s-3.*x1+2.*r1*x1
        +r1s*x1+r2s*x1-2.*x2+2.*r1*x2-2.*r1s*x2+x1*x2)
        /(6.*prop12)
        +3.*(4.-4.*r1s+4.*r2s-3.*x1+2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        +2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -3.*(3.+r1-5.*r1s+r1c+3.*r2s-r1*r2s-2.*x1+r1*x1+r1s*x1-4.*x2
        +2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +4.*(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1-r2s*x1
        -3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /(3.*prop2s);
        rFO2 = 3.*rFO2/4.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = 8.*x1*(-1.-r1s-r2s+x1)
        /(3.*prop1s)
        +6.*(-1-r1s-r2s+2.*x1+r2s*x1-x1s*0.5+x2+r1s*x2-x1*x2*0.5)
        /prop13
        +(2.+2.*r1s+2.*r2s-3.*x1+r1s*x1+r2s*x1-2.*x2-2.*r1s*x2+x1*x2)
        /(3.*prop12)
        +6.*(4.-4.*r1s+4.*r2s-3.*x1+r1s*x1-r2s*x1-5.*x2+r1s*x2-r2s*x2
        +x1*x2+x2s)
        /x3s
        -6.*(3.-5.*r1s+3.*r2s-2.*x1+r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +8.*(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2-r2s*x2
        +x1*x2+x2s)
        /(3.*prop2s);
        rFO4 = 3.*rFO4/8.;
        isSet4 = true;
      }
      break;

    // ~q -> q ~g.
    case 14:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s-r2s-2.*r1*r2);
        rFO1 = 64.*(1.+r1s+2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -16.*(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q
        +x1-r1s*x1+2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -16.*(r1s+r1q-2.*r1c*r2+r2s-6.*r1s*r2s-2.*r1*r2c+r2q-r1s*x1
        +r1*r2*x1+2.*r2s*x1+2.*r1s*x2+r1*r2*x2-r2s*x2-x1*x2)
        /prop12
        -64.*(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /(9.*prop2s)
        +8.*(-1.+r1q-2.*r1*r2+2.*r1c*r2-2.*r2s-2.*r1*r2c-r2q-2.*r1s*x1
        +2.*r2s*x1+x1s+x2-3.*r1s*x2-2.*r1*r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-2.*r1s-r1q-2.*r1*r2-2.*r1c*r2+2.*r1*r2c+r2q+x1+r1s*x1
        -2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO1 = 9.*rFO1/64.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s-r2s+2.*r1*r2);
        rFO2 = 64.*(1.+r1s-2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -16.*(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1-2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -64.*(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /(9.*prop2s)
        +16.*(-r1s-r1q-2.*r1c*r2-r2s+6.*r1s*r2s-2.*r1*r2c-r2q+r1s*x1
        +r1*r2*x1-2.*r2s*x1-2.*r1s*x2+r1*r2*x2+r2s*x2+x1*x2)
        /prop12
        +8.*(-1.+r1q+2.*r1*r2-2.*r1c*r2-2.*r2s+2.*r1*r2c-r2q-2.*r1s*x1
        +2.*r2s*x1+x1s+x2-3.*r1s*x2+2.*r1*r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-2.*r1s-r1q+2.*r1*r2+2.*r1c*r2-2.*r1*r2c+r2q+x1+r1s*x1+
        2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO2 = 9.*rFO2/64.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = 128.*(1.+r1s+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -32*(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+3.*r2s*x1+x2
        +r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -32.*(r1s+r1q+r2s-6.*r1s*r2s+r2q-r1s*x1+2.*r2s*x1+2.*r1s*x2
        -r2s*x2-x1*x2)
        /prop12
        -128.*(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1+x2+3.*r1s*x2
        -r2s*x2-x1*x2)
        /(9.*prop2s)
        +16.*(-1.+r1q-2.*r2s-r2q-2.*r1s*x1+2.*r2s*x1+x1s
        +x2-3.*r1s*x2+r2s*x2+x1*x2)
        /prop13
        -16.*(-1.-2.*r1s-r1q+r2q+x1+r1s*x1-3.*r2s*x1
        +2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO4 = 9.*rFO4/128.;
        isSet4 = true;
      }
      break;

    // q -> ~q ~g.
    case 15:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s+r2s+2.*r2);
        rFO1 = 32*(2.*r2+x2)*(-1.-r1s-r2s+x2)
        /(9.*prop2s)
        +8.*(-1.-r1s-2.*r1s*r2-r2s-2.*r2c+x1+r2*x1+r2s*x1
        +3.*x2*0.5-r1s*x2*0.5+r2*x2-r2s*x2*0.5-x1*x2*0.5)
        /prop12
        +8.*(2.+2.*r1s-2.*r2-2.*r1s*r2-6.*r2s-2.*r2c-3.*x1-r1s*x1
        +2.*r2*x1+3.*r2s*x1+x1s-x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        +32.*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1-2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2-2.*r2*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        -8.*(3.+3.*r1s-r2+r1s*r2-5.*r2s-r2c-4.*x1-r1s*x1
        +2.*r2s*x1+x1s-2.*x2-r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-r1s+r2+r1s*r2-r2s-r2c+x1+r2*x1+r2s*x1+2.*x2+r1s*x2
        -x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO1 = 9.*rFO1/32.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s+r2s-2.*r2);
        rFO2 = 32*(2.*r2-x2)*(1.+r1s+r2s-x2)
        /(9.*prop2s)
        +8.*(-1.-r1s+2.*r1s*r2-r2s+2.*r2c+x1-r2*x1+r2s*x1
        +3.*x2*0.5-r1s*x2*0.5-r2*x2-r2s*x2*0.5-x1*x2*0.5)
        /prop12
        +8.*(2.+2.*r1s+2.*r2+2.*r1s*r2-6.*r2s+2.*r2c-3.*x1-r1s*x1
        -2.*r2*x1+3.*r2s*x1+x1s-x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        -8.*(3.+3.*r1s+r2-r1s*r2-5.*r2s+r2c-4.*x1-r1s*x1+2.*r2s*x1+x1s
        -2.*x2+r2*x2+r2s*x2+x1*x2)
        /prop13
        +32*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+2.*r2*x1+r2s*x1
        +x1s-3.*x2-r1s*x2+2.*r2*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        -8.*(-1.-r1s-r2-r1s*r2-r2s+r2c+x1-r2*x1+r2s*x1+2.*x2+r1s*x2
        -x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO2 = 9.*rFO2/32.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s+r2s);
        rFO4 = 64.*x2*(-1.-r1s-r2s+x2)
        /(9.*prop2s)
        +16.*(-1.-r1s-r2s+x1+r2s*x1+3.*x2*0.5-r1s*x2*0.5
        -r2s*x2*0.5-x1*x2*0.5)
        /prop12
        -16.*(3.+3.*r1s-5.*r2s-4.*x1-r1s*x1+2.*r2s*x1+x1s-2.*x2+r2s*x2
        +x1*x2)
        /prop13
        +64.*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+r2s*x1+x1s-3.*x2
        -r1s*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        +16.*(2.+2.*r1s-6.*r2s-3.*x1-r1s*x1+3.*r2s*x1+x1s
        -x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        -16.*(-1.-r1s-r2s+x1+r2s*x1+2.*x2+r1s*x2-x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO4 = 9.*rFO4/64.;
        isSet4 = true;
      }
      break;

    // g -> ~g ~g. Use (9/4)*eikonal. May be changed in the future.
    case 16:
      rLO = ps;
      if      (combi == 2) offset = x3s;
      else if (combi == 3) offset = mix * x3s;
      else if (combi == 4) offset = 0.5 * x3s;
      rFO = ps * 4.5 * ( (x1+x2-1.+offset-r1s-r2s)/prop12
      - r1s/prop2s - r2s/prop1s );
      break;

    // Dv -> qv d.
    case 30:
      rLO = ps*(1.-r1s+r2s+2.*r2);
      rFO = ( 0.5*r3s + 2.*r1q + 0.5*r2s*r3s + r2*r3s - 2.*r1s
             - 0.5*r1s*r3s - 2.*r1s*r2s - 4.*r1s*r2 ) / prop2s
          + ( -2. + 2.*r2q + 2.*r1q + 2.*r2s*r3s - 4.*r2 + 2.*r2*r3s
             + 4.*r2*r2s - 4.*r1s*r2s - 4.*r1s*r2 ) /prop23
          + ( -2. - 0.5*r3s - 2.*r2s - 4.*r2 + 2.*r1s ) / prop2
          + ( -2. - r3s - 2.*r2s - r2s*r3s - 4.*r2 - 2.*r2*r3s
             + 2.*r1s + r1s*r3s ) / prop3s
          + ( -1. - r3s - r2s - 4.*r2 + r1s - x2 ) / prop3
          + 1.;
      break;

    // S -> Dv Dvbar
    case 31:
      rLO = ps*(1.-4.*r1s);
      rFO = (r3s + 2.*r1s) * (-1. + 4.*r1s) * (1./prop1s + 1./prop2s)
          + (-1. + 8.*r1s - x2) / prop1
          + (-1. + 8.*r1s - x1) / prop2
          + 2. * (1. - 6.*r1s + 8.*r1q + 4.*r3s*r1s) / prop12
          + 2.;
      break;

     // q -> q~ W
    case 32:
      rLO = 1.;
      rFO = (2. * r3s * r3s + 2. * r3s * (x1 + x2) + x1s + x2s) / prop12
          - r3s / prop1s - r3s / prop2s;
      break;

    // q -> q Z
    case 33:
      rLO = 1.;
      rFO = (2. * r3s * r3s + 2. * r3s * (x1 + x2) + x1s + x2s) / prop12
          - r3s / prop1s - r3s / prop2s;
      break;

    // Eikonal expression for kind == 1; also acts as default.
    default:
      rLO = ps;
      if      (combi == 2) offset = x3s;
      else if (combi == 3) offset = mix * x3s;
      else if (combi == 4) offset = 0.5 * x3s;
      rFO = ps * 2. * ( (x1+x2-1.+offset-r1s-r2s)/prop12
      - r1s/prop2s - r2s/prop1s );
      break;

  // End of ME cases.
  }

  // Find relevant leading and first order expressions.
  if      (combi == 1 && isSet1) {
    rLO = rLO1;
    rFO = rFO1; }
  else if (combi == 2 && isSet2) {
    rLO = rLO2;
    rFO = rFO2; }
  else if (combi == 3 && isSet1 && isSet2) {
    rLO = mix * rLO1 + (1.-mix) * rLO2;
    rFO = mix * rFO1 + (1.-mix) * rFO2; }
  else if (isSet4) {
    rLO = rLO4;
    rFO = rFO4; }
  else if (combi == 4 && isSet1 && isSet2) {
    rLO = 0.5 * (rLO1 + rLO2);
    rFO = 0.5 * (rFO1 + rFO2); }
  else if (isSet1) {
    rLO = rLO1;
    rFO = rFO1; }

  // Return ratio of first to leading order cross section.
  return rFO / rLO;
}

//--------------------------------------------------------------------------

// Return the ME corrections for weak t-channel processes.

double SimpleTimeShower::findMEcorrWeak(TimeDipoleEnd* dip,Vec4 rad,
  Vec4 rec, Vec4 emt,Vec4 p3,Vec4 p4,Vec4 radBef, Vec4 recBef) {

  // Check that it is weak emission.
  if (dip->MEtype > 210 || dip->MEtype < 200) return 1.;

  // Remove double counting. Only implemented for QCD hard processes
  // and for the first emission.
  bool cut = false;
  if (infoPtr->nISR() + infoPtr->nFSRinProc() == 0
   && infoPtr->code() > 110 && infoPtr->code() < 130
   && vetoWeakJets) {
    double d = emt.pT2();
    if (rad.pT2() < d) {d = rad.pT2(); cut = true;}
    if (rec.pT2() < d) {d = rec.pT2(); cut = true;}

    // Always check for combination of radiator and emitted.
    double dij = min(rad.pT2(),emt.pT2())
      * pow2(RRapPhi(rad,emt)) / vetoWeakDeltaR2;
    if (dij < d) {
      d = dij;
      cut = false;
    }

    // Check for angle between recoiler and radiator, if quark anti-quark pair,
    // or if the recoiler is a gluon.
    if (dip->MEtype == 200 || dip->MEtype == 205 ||
        dip->MEtype == 201 || dip->MEtype == 206) {
      dij = min(rad.pT2(),rec.pT2()) * pow2(RRapPhi(rad,rec))
          / vetoWeakDeltaR2;
      if (dij < d) {
        d = dij;
        cut = true;
      }
    }
    // Check for angle between recoiler and emitted, if recoiler is a quark.
    if (dip->MEtype == 200 || dip->MEtype == 205 ||
        dip->MEtype == 202 || dip->MEtype == 207 ||
        dip->MEtype == 203 || dip->MEtype == 208) {
      dij = min(emt.pT2(),rec.pT2()) * pow2(RRapPhi(emt,rec))
          / vetoWeakDeltaR2;
      if (dij < d) {
        d = dij;
        cut = false;
      }
    }
    if (cut) return 0.;
  }

  // Check that MEtype is t-channel weak emission.
  if ( dip->MEtype != 201 && dip->MEtype != 202 && dip->MEtype != 203
    && dip->MEtype != 206 && dip->MEtype != 207 && dip->MEtype != 208)
    return 1;

  // Rescaling of incoming partons p3 and p4.
  double scaleFactor2 = (rad + rec + emt).m2Calc() / (p3 + p4).m2Calc();
  double scaleFactor = sqrt(scaleFactor2);
  p3 *= scaleFactor;
  p4 *= scaleFactor;

  // Longitudinal boost to rest frame of incoming partons of hard interaction.
  RotBstMatrix rot2to2frame;
  rot2to2frame.bstback(p3 + p4);
  p3.rotbst(rot2to2frame);
  p4.rotbst(rot2to2frame);
  rad.rotbst(rot2to2frame);
  emt.rotbst(rot2to2frame);
  rec.rotbst(rot2to2frame);
  recBef.rotbst(rot2to2frame);
  radBef.rotbst(rot2to2frame);

  // Further boost to rest frame of outgoing state.
  RotBstMatrix rot2to3frame;
  rot2to3frame.bstback(rad + emt + rec);
  rad.rotbst(rot2to3frame);
  emt.rotbst(rot2to3frame);
  rec.rotbst(rot2to3frame);
  recBef.rotbst(rot2to3frame);
  radBef.rotbst(rot2to3frame);

  // Kinematical quantities.
  double sHat = (p3 + p4).m2Calc();
  double tHat = (radBef - p3).m2Calc();
  double uHat = (recBef - p3).m2Calc();
  double z    = dip->z;
  double pT2  = dip->pT2;
  double Q2   = pT2 / (z*(1.-z));

  // ME weight. Prefactor mainly from Jacobian.
  double wt = 2. * pT2 / z * (Q2+sHat)/sHat * (1. - kRad - kEmt) / 4.;
  if (dip->MEtype == 201 || dip->MEtype == 206)
    wt *= simpleWeakShowerMEs.getMEqg2qgZ( p3, p4, rec, emt, rad)
        / simpleWeakShowerMEs.getMEqg2qg( sHat, tHat, uHat);
  else if (dip->MEtype == 202 || dip->MEtype == 207)
    wt *= simpleWeakShowerMEs.getMEqq2qqZ( p3, p4, emt, rec, rad)
        / simpleWeakShowerMEs.getMEqq2qq( sHat, tHat, uHat, true);
  else if (dip->MEtype == 203 || dip->MEtype == 208)
    wt *= simpleWeakShowerMEs.getMEqq2qqZ( p3, p4, emt, rec, rad)
        / simpleWeakShowerMEs.getMEqq2qq( sHat, tHat, uHat, false);

  // Split of ME into an ISR part and FSR part.
  wt *= abs((-emt + p3).m2Calc()) / ((emt + rad).m2Calc()
      + abs((-p3 + emt).m2Calc()));

  // Correction for previous fudge-factor enhancement of weak emission rate.
  wt /= WEAKPSWEIGHT;
  if (wt > 1.) infoPtr->errorMsg("Warning in SimpleTimeShower::findMEcorrWeak:"
    " weight is above unity");
  return wt;

}

//--------------------------------------------------------------------------

// Find coefficient of azimuthal asymmetry from gluon polarization.

void SimpleTimeShower::findAsymPol( Event& event, TimeDipoleEnd* dip) {

  // Default is no asymmetry. Only gluons are studied.
  dip->asymPol = 0.;
  dip->iAunt = 0;
  int iRad = dip->iRadiator;
  if (!doPhiPolAsym || event[iRad].id() != 21) return;

  // Trace grandmother via possibly intermediate recoil copies.
  int iMother = event[iRad].iTopCopy();
  int iGrandM = event[iMother].mother1();

  // If grandmother in initial state of hard scattering,
  // then at most keep only gg and qq initial states.
  int statusGrandM = event[iGrandM].status();
  bool isHardProc  = (statusGrandM == -21 || statusGrandM == -31);
  if (isHardProc) {
    if (!doPhiPolAsymHard) return;
    if (event[iGrandM + 1].status() != statusGrandM) return;
    if (event[iGrandM].isGluon() && event[iGrandM + 1].isGluon());
    else if (event[iGrandM].isQuark() && event[iGrandM + 1].isQuark());
    else return;
  }

  // Set aunt by history or, for hard scattering, by colour flow.
  if (isHardProc) dip->iAunt = dip->iRecoiler;
  else dip->iAunt = (event[iGrandM].daughter1() == iMother)
    ? event[iGrandM].daughter2() : event[iGrandM].daughter1();

  // Coefficient from gluon production (approximate z by energy).
  // For hard process arbitrarily put z = 1/2.
  double zProd = (isHardProc) ? 0.5 : event[iRad].e()
    / (event[iRad].e() + event[dip->iAunt].e());
  if (event[iGrandM].isGluon()) dip->asymPol = pow2( (1. - zProd)
    / (1. - zProd * (1. - zProd) ) );
  else dip->asymPol = 2. * (1. - zProd) / (1. + pow2(1. - zProd) );

  // Coefficients from gluon decay.
  if (dip->flavour == 21) dip->asymPol *= pow2( dip->z * (1. - dip->z)
    / (1. - dip->z * (1. - dip->z) ) );
  else  dip->asymPol *= -2. * dip->z * ( 1. - dip->z )
    / (1. - 2. * dip->z * (1. - dip->z) );

}

//--------------------------------------------------------------------------

// Print the list of dipoles.

void SimpleTimeShower::list() const {

  // Header.
  cout << "\n --------  PYTHIA SimpleTimeShower Dipole Listing  -----------"
       << "------------------------------------------------------- \n \n  "
       << "  i    rad    rec       pTmax  col  chg  gam weak  oni   hv  is"
       << "r  sys sysR type  MErec     mix  ord  spl  ~gR  pol \n"
       << fixed << setprecision(3);

  // Loop over dipole list and print it.
  for (int i = 0; i < int(dipEnd.size()); ++i)
  cout << setw(5) << i                     << setw(7) << dipEnd[i].iRadiator
       << setw(7) << dipEnd[i].iRecoiler   << setw(12) << dipEnd[i].pTmax
       << setw(5) << dipEnd[i].colType     << setw(5) << dipEnd[i].chgType
       << setw(5) << dipEnd[i].gamType     << setw(5) << dipEnd[i].weakType
       << setw(5) << dipEnd[i].isOctetOnium
       << setw(5) << dipEnd[i].isHiddenValley << setw(5) << dipEnd[i].isrType
       << setw(5) << dipEnd[i].system      << setw(5) << dipEnd[i].systemRec
       << setw(5) << dipEnd[i].MEtype      << setw(7) << dipEnd[i].iMEpartner
       << setw(8) << dipEnd[i].MEmix       << setw(5) << dipEnd[i].MEorder
       << setw(5) << dipEnd[i].MEsplit     << setw(5) << dipEnd[i].MEgluinoRec
       << setw(5) << dipEnd[i].weakPol << "\n";

  // Done.
  cout << "\n --------  End PYTHIA SimpleTimeShower Dipole Listing  -------"
       << "-------------------------------------------------------" << endl;

}

//==========================================================================

} // end namespace Pythia8
