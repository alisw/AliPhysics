// PartonLevel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Hard diffraction added by Christine Rasmussen.

// Function definitions (not found in the header) for the PartonLevel class.

#include "Pythia8/PartonLevel.h"

namespace Pythia8 {

//==========================================================================

// The PartonLevel class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to produce parton level from given input.
const int PartonLevel::NTRY = 10;

//--------------------------------------------------------------------------

// Main routine to initialize the parton-level generation process.

bool PartonLevel::init( Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
  BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  BeamParticle* beamPomAPtrIn, BeamParticle* beamPomBPtrIn,
  BeamParticle* beamGamAPtrIn, BeamParticle* beamGamBPtrIn,
  BeamParticle* beamVMDAPtrIn, BeamParticle* beamVMDBPtrIn,
  Couplings* couplingsPtrIn, PartonSystems* partonSystemsPtrIn,
  SigmaTotal* sigmaTotPtr, TimeShower* timesDecPtrIn, TimeShower* timesPtrIn,
  SpaceShower* spacePtrIn, RHadrons* rHadronsPtrIn, UserHooks* userHooksPtrIn,
  MergingHooks* mergingHooksPtrIn, PartonVertex* partonVertexPtrIn,
  bool useAsTrial ) {

  // Store input pointers and modes for future use.
  infoPtr            = infoPtrIn;
  particleDataPtr    = particleDataPtrIn;
  rndmPtr            = rndmPtrIn;
  beamAPtr           = beamAPtrIn;
  beamBPtr           = beamBPtrIn;
  beamHadAPtr        = beamAPtr;
  beamHadBPtr        = beamBPtr;
  beamPomAPtr        = beamPomAPtrIn;
  beamPomBPtr        = beamPomBPtrIn;
  beamGamAPtr        = beamGamAPtrIn;
  beamGamBPtr        = beamGamBPtrIn;
  beamVMDAPtr        = beamVMDAPtrIn;
  beamVMDBPtr        = beamVMDBPtrIn;
  couplingsPtr       = couplingsPtrIn;
  partonSystemsPtr   = partonSystemsPtrIn;
  timesDecPtr        = timesDecPtrIn;
  timesPtr           = timesPtrIn;
  spacePtr           = spacePtrIn;
  rHadronsPtr        = rHadronsPtrIn;
  userHooksPtr       = userHooksPtrIn;
  mergingHooksPtr    = mergingHooksPtrIn;
  partonVertexPtr    = partonVertexPtrIn;

  // Min bias and diffraction processes need special treatment.
  bool doSQ          = settings.flag("SoftQCD:all")
                    || settings.flag("SoftQCD:inelastic");
  bool doND          = settings.flag("SoftQCD:nonDiffractive");
  bool doSD          = settings.flag("SoftQCD:singleDiffractive");
  bool doDD          = settings.flag("SoftQCD:doubleDiffractive");
  bool doCD          = settings.flag("SoftQCD:centralDiffractive");
  doNonDiff          = doSQ || doND;
  doDiffraction      = doSQ || doSD || doDD || doCD;
  doHardDiff         = settings.flag("Diffraction:doHard");
  hardDiffSide       = (doHardDiff) ? settings.mode("Diffraction:hardDiffSide")
                     : 0;
  sampleTypeDiff     = (doHardDiff) ? settings.mode("Diffraction:sampleType")
                     : 0;

  // Separate low-mass (unresolved) and high-mass (perturbative) diffraction.
  mMinDiff           = settings.parm("Diffraction:mMinPert");
  mWidthDiff         = settings.parm("Diffraction:mWidthPert");
  pMaxDiff           = settings.parm("Diffraction:probMaxPert");
  if (mMinDiff > infoPtr->eCM()) doDiffraction = false;

  // Set whether photon inside lepton. Mode updated event-by-event.
  gammaMode          = settings.mode("Photon:ProcessType");
  gammaModeEvent     = 0;
  beamHasGamma       = settings.flag("PDF:lepton2gamma")
                    && (beamAPtr != 0) && (beamBPtr != 0);
  hasGammaA          = false;
  hasGammaB          = false;
  beamAisGamma       = (beamAPtr != 0) ? beamAPtr->isGamma() : false;
  beamBisGamma       = (beamBPtr != 0) ? beamBPtr->isGamma() : false;
  beamAhasGamma      = (beamHasGamma && beamAPtr->isLepton());
  beamBhasGamma      = (beamHasGamma && beamBPtr->isLepton());
  beamAhasResGamma   = (beamAPtr != 0) ? beamAPtr->hasResGamma() : false;
  beamBhasResGamma   = (beamBPtr != 0) ? beamBPtr->hasResGamma() : false;
  beamHasResGamma    = (gammaMode < 4) && beamHasGamma;
  isGammaHadronDir   = false;
  bool isAgamBhad    = (beamAisGamma || beamAhasGamma) && beamBPtr->isHadron();
  bool isBgamAhad    = (beamBisGamma || beamBhasGamma) && beamAPtr->isHadron();
  bool isAgamBgam    = (beamAisGamma || beamAhasGamma)
                    && (beamBisGamma || beamBhasGamma);
  bool onlyDirGamma  = false;
  if ( gammaMode == 4 || ( gammaMode == 3 && isAgamBhad )
    || ( gammaMode == 2 && isBgamAhad ) || ( gammaMode > 1 && isAgamBgam) )
    onlyDirGamma = true;

  // Show the copies of beam photon if found in ISR.
  showUnresGamma     = settings.flag("Photon:showUnres");

  // Need MPI initialization for soft QCD processes, even if only first MPI.
  // But no need to initialize MPI if never going to use it.
  doMPI              = settings.flag("PartonLevel:MPI");
  doMPIMB            = doMPI;
  doMPISDA           = doMPI;
  doMPISDB           = doMPI;
  doMPICD            = doMPI;
  doMPIinit          = doMPI;
  doMPIgmgm          = doMPI;
  if (doNonDiff || doDiffraction)        doMPIinit = true;
  if (!settings.flag("PartonLevel:all")) doMPIinit = false;

  // Nature of MPI matching also used here for one case.
  pTmaxMatchMPI      = settings.mode("MultipartonInteractions:pTmaxMatch");

  // Initialise trial shower switch.
  doTrial            = useAsTrial;
  // Merging initialization.
  bool hasMergingHooks = (mergingHooksPtr != 0);
  canRemoveEvent       = !doTrial && hasMergingHooks
    && ( mergingHooksPtr->doCKKWLMerging() || mergingHooksPtr->doNL3Merging());
  canRemoveEmission    = !doTrial && hasMergingHooks
    && ( mergingHooksPtr->doUMEPSMerging() || mergingHooksPtr->doNL3Merging()
      || mergingHooksPtr->doUNLOPSMerging() );
  nTrialEmissions    = 1;
  pTLastBranch       = 0.0;
  typeLastBranch     = 0;

  // Flags for showers: ISR and FSR.
  doISR              = settings.flag("PartonLevel:ISR");
  bool FSR           = settings.flag("PartonLevel:FSR");
  bool FSRinProcess  = settings.flag("PartonLevel:FSRinProcess");
  bool interleaveFSR = settings.flag("TimeShower:interleave");
  doFSRduringProcess = FSR && FSRinProcess &&  interleaveFSR;
  doFSRafterProcess  = FSR && FSRinProcess && !interleaveFSR;
  doFSRinResonances  = FSR && settings.flag("PartonLevel:FSRinResonances");

  // Flags for colour reconnection.
  doReconnect        = settings.flag("ColourReconnection:reconnect");
  reconnectMode      = settings.mode("ColourReconnection:mode");
  forceResonanceCR   = settings.flag("ColourReconnection:forceResonance");

  // Some other flags.
  doRemnants         = settings.flag("PartonLevel:Remnants");
  doSecondHard       = settings.flag("SecondHard:generate");
  twoHard            = doSecondHard;
  earlyResDec        = settings.flag("PartonLevel:earlyResDec");

  // Allow R-hadron formation.
  allowRH            = settings.flag("RHadrons:allow");

  // Possibility to allow user veto during evolution.
  canVetoPT          = (userHooksPtr != 0)
                     ? userHooksPtr->canVetoPT()   : false;
  pTvetoPT           = (canVetoPT)
                     ? userHooksPtr->scaleVetoPT() : -1.;
  canVetoStep        = (userHooksPtr != 0)
                     ? userHooksPtr->canVetoStep() : false;
  nVetoStep          = (canVetoStep)
                     ? userHooksPtr->numberVetoStep() : -1;
  canVetoMPIStep     = (userHooksPtr != 0)
                     ? userHooksPtr->canVetoMPIStep() : false;
  nVetoMPIStep       = (canVetoMPIStep)
                     ? userHooksPtr->numberVetoMPIStep() : -1;
  canVetoEarly       = (userHooksPtr != 0)
                     ? userHooksPtr->canVetoPartonLevelEarly() : false;

  // Settings for vetoing of QCD emission for Drell-Yan weak boson production.
  vetoWeakJets       = settings.flag("WeakShower:vetoQCDjets");
  vetoWeakDeltaR2    = pow2(settings.parm("WeakShower:vetoWeakDeltaR"));

  // Possibility to set maximal shower scale in resonance decays.
  canSetScale        = (userHooksPtr != 0)
                     ? userHooksPtr->canSetResonanceScale() : false;

  // Possibility to reconnect specifically for resonance decays.
  canReconResSys     = (userHooksPtr != 0)
                     ? userHooksPtr->canReconnectResonanceSystems() : false;

  // Done with initialization only for FSR in resonance decays.
  if (beamAPtr == 0 || beamBPtr == 0) return true;

  // Make sure that photons are in resolved mode when mixing with unresolved
  // before initializing MPIs.
  if ( (beamAisGamma || beamAhasGamma) && gammaMode == 0) {
    beamAPtr->setGammaMode(1);
    if (beamAhasGamma) beamGamAPtr->setGammaMode(1);
  }
  if ( (beamBisGamma || beamBhasGamma) && gammaMode == 0) {
    beamBPtr->setGammaMode(1);
    if (beamBhasGamma) beamGamBPtr->setGammaMode(1);
  }

  // If only direct photons do not initialize diffractive MPI systems.
  if ( ( (gammaMode == 3) &&
    ( (beamAPtr->isGamma() || beamAhasGamma) && beamBPtr->isHadron() ) )
    || ( (gammaMode == 2) &&
    ( (beamBPtr->isGamma() || beamBhasGamma) && beamAPtr->isHadron() ) )
    || ( (gammaMode > 1) &&
    (  (beamAPtr->isGamma() || beamAhasGamma)
    && (beamBPtr->isGamma() || beamBhasGamma) ) ) )
    onlyDirGamma = true;

  // Flag if lepton beams, and if non-resolved ones. May change main flags.
  hasTwoLeptonBeams  =  beamAPtr->isLepton() && beamBPtr->isLepton();
  hasOneLeptonBeam   = (beamAPtr->isLepton() || beamBPtr->isLepton())
                    && !hasTwoLeptonBeams;
  hasPointLeptons    = (hasOneLeptonBeam || hasTwoLeptonBeams)
    && (beamAPtr->isUnresolved() || beamBPtr->isUnresolved());
  if ( (hasOneLeptonBeam || hasTwoLeptonBeams) && !beamHasResGamma ) {
    doMPIMB          = false;
    doMPISDA         = false;
    doMPISDB         = false;
    doMPICD          = false;
    doMPIinit        = false;
    doMPIgmgm        = false;
  }
  if (hasTwoLeptonBeams && hasPointLeptons) {
    doISR            = false;
    doRemnants       = false;
  }

  // For ND events in lepton->gamma events no need to initialize MPIs for l+l-.
  doNDgamma = false;
  if (beamHasResGamma)         doMPIinit = false;
  if (beamHasResGamma && doND) doNDgamma = true;

  // Set info and initialize the respective program elements.
  timesPtr->init( beamAPtr, beamBPtr);
  if (doISR) spacePtr->init( beamAPtr, beamBPtr);

  doMPIMB  =  multiMB.init( doMPIinit, 0, infoPtr, settings, particleDataPtr,
    rndmPtr, beamAPtr, beamBPtr, couplingsPtr,
    partonSystemsPtr, sigmaTotPtr, userHooksPtr, partonVertexPtr);

  // Initialize MPIs for diffractive system, possibly photon beam from
  // lepton, possibly VMD from photon.
  if (doSD || doDD || doSQ || ( doHardDiff && (hardDiffSide == 0
    || hardDiffSide == 1) && beamBPtr->getGammaMode() < 2 ) ) {
    BeamParticle* tmpBeamA = (beamAhasGamma) ? beamGamAPtr : beamAPtr;
    if (infoPtr->isVMDstateA()) tmpBeamA = beamVMDAPtr;
    doMPISDA = multiSDA.init( !onlyDirGamma, 1, infoPtr, settings,
      particleDataPtr, rndmPtr, tmpBeamA, beamPomBPtr, couplingsPtr,
      partonSystemsPtr, sigmaTotPtr, userHooksPtr, partonVertexPtr,
      (beamAisGamma || beamAhasGamma) );
  }
  if (doSD || doDD || doSQ || ( doHardDiff && (hardDiffSide == 0
    || hardDiffSide == 2) && beamAPtr->getGammaMode() < 2 ) ) {
    BeamParticle* tmpBeamB = (beamBhasGamma) ? beamGamBPtr : beamBPtr;
    if (infoPtr->isVMDstateB()) tmpBeamB = beamVMDBPtr;
    doMPISDB = multiSDB.init( !onlyDirGamma, 2, infoPtr, settings,
      particleDataPtr, rndmPtr, beamPomAPtr, tmpBeamB, couplingsPtr,
      partonSystemsPtr, sigmaTotPtr, userHooksPtr, partonVertexPtr,
      (beamBisGamma || beamBhasGamma) );
  }
  if (doCD || doSQ) doMPICD = multiCD.init( doMPIinit, 3, infoPtr, settings,
    particleDataPtr, rndmPtr, beamPomAPtr, beamPomBPtr, couplingsPtr,
    partonSystemsPtr, sigmaTotPtr, userHooksPtr, partonVertexPtr);
  if (!remnants.init( infoPtr, settings, rndmPtr, beamAPtr, beamBPtr,
    partonSystemsPtr, partonVertexPtr, particleDataPtr, &colourReconnection))
    return false;
  resonanceDecays.init( infoPtr, particleDataPtr, rndmPtr);
  colourReconnection.init( infoPtr, settings, rndmPtr, particleDataPtr,
    beamAPtr, beamBPtr, partonSystemsPtr);
  junctionSplitting.init(infoPtr, settings, rndmPtr, particleDataPtr);

  // Initialize hard diffraction with possible photon beams from leptons.
  if (doHardDiff && (gammaMode !=4) )
    hardDiffraction.init(infoPtr, settings, rndmPtr, beamAhasGamma ?
      beamGamAPtr : beamAPtr, beamBhasGamma ? beamGamBPtr : beamBPtr,
      beamPomAPtr, beamPomBPtr, sigmaTotPtr);

  // Initialize an MPI instance for photons from leptons.
  if ( beamHasResGamma && (doMPI || doNDgamma) ) {
    doMPIinit = true;
    // Lepton-hadron.
    if (beamAPtr->isLepton() && beamBPtr->isHadron() ) {
      doMPIgmgm = multiGmGm.init( doMPIinit, 0, infoPtr, settings,
        particleDataPtr, rndmPtr, beamGamAPtr, beamBPtr, couplingsPtr,
        partonSystemsPtr, sigmaTotPtr, userHooksPtr, partonVertexPtr, true);
    // Hadron-lepton.
    } else if (beamBPtr->isLepton() && beamAPtr->isHadron() ) {
      doMPIgmgm = multiGmGm.init( doMPIinit, 0, infoPtr, settings,
        particleDataPtr, rndmPtr, beamAPtr, beamGamBPtr, couplingsPtr,
        partonSystemsPtr, sigmaTotPtr, userHooksPtr, partonVertexPtr, true);
    // Lepton-lepton.
    } else {
      doMPIgmgm = multiGmGm.init( doMPIinit, 0, infoPtr, settings,
        particleDataPtr, rndmPtr, beamGamAPtr, beamGamBPtr, couplingsPtr,
        partonSystemsPtr, sigmaTotPtr, userHooksPtr, partonVertexPtr, true);
    }
    doMPIMB = doMPIgmgm;
  }

  // Succeeded, or not.
  multiPtr       = &multiMB;
  if (doMPIinit && !doMPIMB) return false;
  if (doMPIinit && (doSD || doDD || doSQ) && (!doMPISDA || !doMPISDB))
     return false;
  if (doMPIinit && (doCD || doSQ) && !doMPICD) return false;
  if (!doMPIMB || !doMPISDA || !doMPISDB || !doMPICD) doMPI = false;
  return true;

}

//--------------------------------------------------------------------------

// Function to reset PartonLevel object for trial shower usage.

void PartonLevel::resetTrial() {

  // Clear input pointers.
  partonSystemsPtr->clear();
  beamAPtr->clear();
  beamBPtr->clear();
  beamHadAPtr->clear();
  beamHadBPtr->clear();
  beamPomAPtr->clear();
  beamPomBPtr->clear();
  beamGamAPtr->clear();
  beamGamBPtr->clear();
  beamVMDAPtr->clear();
  beamVMDBPtr->clear();

  // Clear last branching return values.
  pTLastBranch   = 0.0;
  typeLastBranch = 0;

}

//--------------------------------------------------------------------------

// Main routine to do the parton-level evolution.

bool PartonLevel::next( Event& process, Event& event) {

  // Current event classification.
  isResolved        = infoPtr->isResolved();
  isResolvedA       = isResolved;
  isResolvedB       = isResolved;
  isResolvedC       = isResolved;
  isDiffA           = infoPtr->isDiffractiveA();
  isDiffB           = infoPtr->isDiffractiveB();
  isDiffC           = infoPtr->isDiffractiveC();
  isDiff            = isDiffA || isDiffB || isDiffC;
  isCentralDiff     = isDiffC;
  isDoubleDiff      = isDiffA && isDiffB;
  isSingleDiff      = isDiff && !isDoubleDiff  && !isCentralDiff;
  isNonDiff         = infoPtr->isNonDiffractive();
  isElastic         = infoPtr->isElastic();

  // Default values for what is to come with diffraction.
  isHardDiffA       = false;
  isHardDiffB       = false;
  isHardDiff        = false;
  doDiffVeto        = false;
  // Mark hard diffractive events to handle CR correctly.
  bool doDiffCR     = false;
  // The setup of the diffractive events can come after the first evolution.
  int nHardDiffLoop = 1;
  // Flag to check whether hard diffraction system is set up for the process.
  hardDiffSet       = false;
  // Offset for the initiator position when photons from leptons.
  gammaOffset = 0;

  // Parton-level vetoes for matching and merging.
  doVeto            = false;
  infoPtr->setAbortPartonLevel(false);

  // Update photon state according to beams set in processContainer.
  beamAhasResGamma  = beamAPtr->hasResGamma();
  beamBhasResGamma  = beamBPtr->hasResGamma();
  beamHasResGamma   = beamAhasResGamma || beamBhasResGamma;

  // Save if photoproduction from either side.
  hasGammaA         = (beamAPtr != 0) ? beamAPtr->getGammaMode() > 0 : false;
  hasGammaB         = (beamBPtr != 0) ? beamBPtr->getGammaMode() > 0 : false;

  // Save current photon mode when mixing processes.
  gammaModeEvent = gammaMode;
  if ( hasGammaA || hasGammaB ) {
    if (beamAPtr->getGammaMode() < 2 && beamBPtr->getGammaMode() < 2)
      gammaModeEvent = 1;
    if (beamAPtr->getGammaMode() < 2 && beamBPtr->getGammaMode() == 2)
      gammaModeEvent = 2;
    if (beamAPtr->getGammaMode() == 2 && beamBPtr->getGammaMode() < 2)
      gammaModeEvent = 3;
    if (beamAPtr->getGammaMode() == 2 && beamBPtr->getGammaMode() == 2)
      gammaModeEvent = 4;
  }

  // Check if direct-photon + hadron to set hard system correctly.
  isGammaHadronDir = !(hasGammaA || hasGammaB) ? false :
    ( (beamAPtr->getGammaMode() == 2) && (beamBPtr->getGammaMode() == 0) )
    || ( (beamAPtr->getGammaMode() == 0) && (beamBPtr->getGammaMode() == 2) );

  // Set up gamma+gamma/p subcollision. May fail due to extreme kinematics.
  if (beamHasGamma) {
    if ( !setupResolvedLeptonGamma( process) ) return false;
  }

  // Prepare for a potential hard diffractive event.
  if (doHardDiff) {

    // Switch for the diffractive side to be considered.
    bool checkSideA = (hardDiffSide < 2) && (beamBPtr->isHadron() ||
      ((beamBisGamma || beamBhasGamma) && beamBPtr->getGammaMode() == 1));
    bool checkSideB = (hardDiffSide%2 == 0) && (beamAPtr->isHadron() ||
      ((beamAisGamma || beamAhasGamma) && beamAPtr->getGammaMode() == 1));

    // Preliminary decision based on diffractive-to-inclusive PDF ratio.
    // If Pomeron taken from side A(=1), then B is the diffractive system.
    // If Pomeron taken from side B(=2), then A is the diffractive system.
    // Do not check for diffraction if direct photons.
    // Rescale the x values when photon from lepton and calculate PDFs
    // wrt. photon beam.
    if ( checkSideA ) {
      double xGammaB    = beamBhasGamma ? beamHadBPtr->xGamma() : 1.;
      double xBrescaled = infoPtr->x2pdf() / xGammaB;
      double pdfB       = beamBhasGamma ? beamGamBPtr->xf(infoPtr->id2pdf(),
        xBrescaled, infoPtr->Q2Fac()) : infoPtr->pdf2();
      isHardDiffA = hardDiffraction.isDiffractive(2, infoPtr->id2pdf(),
        xBrescaled, infoPtr->Q2Fac(), pdfB);
    }
    if ( checkSideB ) {
      double xGammaA    = beamAhasGamma ? beamHadAPtr->xGamma() : 1.;
      double xArescaled = infoPtr->x1pdf() / xGammaA;
      double pdfA       = beamAhasGamma ? beamGamAPtr->xf(infoPtr->id1pdf(),
        xArescaled, infoPtr->Q2Fac()) : infoPtr->pdf1();
      isHardDiffB = hardDiffraction.isDiffractive(1, infoPtr->id1pdf(),
        xArescaled, infoPtr->Q2Fac(), pdfA);
    }

    // No hard double diffraction yet, so randomly choose one of the sides.
    if (isHardDiffA && isHardDiffB) {
      if (rndmPtr->flat() < 0.5) isHardDiffA = false;
      else isHardDiffB = false;
    }
    isHardDiff = isHardDiffA || isHardDiffB;

    // Save diffractive values.
    double xPomA = (isHardDiffB) ? hardDiffraction.getXPomeronA() : 0.;
    double xPomB = (isHardDiffA) ? hardDiffraction.getXPomeronB() : 0.;
    double tPomA = (isHardDiffB) ? hardDiffraction.getTPomeronA() : 0.;
    double tPomB = (isHardDiffA) ? hardDiffraction.getTPomeronB() : 0.;
    infoPtr->setHardDiff( false, false, isHardDiffA, isHardDiffB,
      xPomA, xPomB, tPomA, tPomB);

    // Discard all nondiffractive events if only diffractive sample is wanted.
    if (!isHardDiff && sampleTypeDiff > 2) {
      doDiffVeto = true;
      // Reset beam pointers and set event back to original frame.
      if (beamHasGamma) leaveResolvedLeptonGamma( process, event, false);
      return false;
    }

    if (isHardDiff) {
      // Set up the diffractive system if run without MPI veto.
      if (sampleTypeDiff%2 == 1) setupHardDiff( process);
      // Allow for second loop if run with MPI veto.
      else nHardDiffLoop = 2;
    }
  }

  // Check if two subcollisions, from secondHard or Les Houches input.
  int n21 = 0;
  for (int i = 1; i < process.size(); ++i)
    if (process[i].status() == -21) ++n21;
  twoHard = (n21 == 4);

  // nHardLoop counts how many hard-scattering subsystems are to be processed.
  // Almost always 1, but elastic and low-mass diffraction gives 0, while
  // double diffraction can give up to 2. Not to be confused with SecondHard.
  int nHardLoop  = 1;
  if (!isResolved) nHardLoop = (isDiff) ? decideResolvedDiff( process) : 0;

  // Handle unresolved subsystems. Done if no resolved ones.
  sizeProcess    = 0;
  sizeEvent      = 0;
  if (!isResolvedA || !isResolvedB || !isResolvedC) {
    bool physical = setupUnresolvedSys( process, event);
    // Set up the scattered photon if present before exiting.
    if (!physical || nHardLoop == 0) {
      if (beamHasGamma) leaveResolvedLeptonGamma( process, event, physical);
      return physical;
    }
    sizeProcess  = process.size();
    sizeEvent    = event.size();
  }

  // Number of actual branchings.
  int nBranch        = 0;
  // Number of desired branchings, negative value means no restriction.
  int nBranchMax     = (doTrial) ? nTrialEmissions : -1;

  // Store merging weight.
  bool hasMergingHooks = (mergingHooksPtr != 0);
  if ( hasMergingHooks && canRemoveEvent )
    mergingHooksPtr->storeWeights(infoPtr->getWeightCKKWL());

  // Reset event weight coming from enhanced branchings.
  if (userHooksPtr != 0) userHooksPtr->setEnhancedEventWeight(1.);

  // Loop to set up diffractive system if run with MPI veto.
  for (int iHardDiffLoop = 1; iHardDiffLoop <= nHardDiffLoop;
    ++iHardDiffLoop) {

  // Big outer loop to handle up to two systems (in double diffraction),
  // but normally one. (Not indented in following, but end clearly marked.)
  for (int iHardLoop = 1; iHardLoop <= nHardLoop; ++iHardLoop) {
    infoPtr->setCounter(20, iHardLoop);
    infoPtr->setCounter(21);

  // Classification of diffractive system: 1 = A, 2 = B, 3 = central.
  iDS = 0;
  if (isDiffA || isDiffB) iDS = (iHardLoop == 2 || !isResolvedA) ? 2 : 1;
  if (isDiffC) iDS = 3;

  // Process and event records can be out of step for diffraction.
  if (iHardLoop == 2) {
    sizeProcess = process.size();
    sizeEvent   = event.size();
    partonSystemsPtr->clear();
    if (event.lastColTag() > process.lastColTag())
      process.initColTag(event.lastColTag());
  }

  // If you need to restore then do not throw existing diffractive system.
  if (isDiff) {
    event.saveSize();
    event.saveJunctionSize();

    // Allow special treatment of diffractive systems.
    setupResolvedDiff( process);
  }

  // Prepare to do multiparton interactions; at new mass for diffraction.
  if (doMPIinit || doDiffraction) multiPtr->reset();

  // Special case if nondiffractive: do hardest interaction.
  if (isNonDiff || isDiff) {
    multiPtr->pTfirst();
    multiPtr->setupFirstSys( process);
  }

  // Allow up to ten tries; failure possible for beam remnants.
  // Main cause: inconsistent colour flow at the end of the day.
  bool physical = true;
  int  nRad     = 0;
  for (int iTry = 0; iTry < NTRY; ++ iTry) {
    infoPtr->addCounter(21);
    for (int i = 22; i < 32; ++i) infoPtr->setCounter(i);

    // Reset flag, counters and max scales.
    physical   = true;
    nMPI       = (twoHard) ? 2 : 1;
    nISR       = 0;
    nFSRinProc = 0;
    nFSRinRes  = 0;
    nISRhard   = 0;
    nFSRhard   = 0;
    pTsaveMPI  = 0.;
    pTsaveISR  = 0.;
    pTsaveFSR  = 0.;

    // Reset nMPI and nISR for showers.
    infoPtr->setPartEvolved(nMPI, nISR);

    // Reset parameters related to valence content and remnants of photon
    // beam if ISR or MPI generated.
    if (beamAPtr->isGamma() && ( doMPI || doISR ) ) beamAPtr->resetGamma();
    if (beamBPtr->isGamma() && ( doMPI || doISR ) ) beamBPtr->resetGamma();

    // Identify hard interaction system for showers.
    setupHardSys( process, event);

    // Optionally check for a veto after the hardest interaction.
    if (canVetoMPIStep) {
      doVeto = userHooksPtr->doVetoMPIStep( 1, event);
      // Abort event if vetoed.
      if (doVeto) {
        if (isDiff) leaveResolvedDiff( iHardLoop, process, event);
        if (beamHasResGamma) leaveResolvedLeptonGamma( process, event, false);
        return false;
      }
    }

    // Check matching of process scale to maximum ISR/FSR/MPI scales.
    double Q2Fac       = infoPtr->Q2Fac();
    double Q2Ren       = infoPtr->Q2Ren();
    bool limitPTmaxISR = (doISR)
      ? spacePtr->limitPTmax( event, Q2Fac, Q2Ren) : false;
    bool limitPTmaxFSR = (doFSRduringProcess)
      ? timesPtr->limitPTmax( event, Q2Fac, Q2Ren) : false;
    bool limitPTmaxMPI = (doMPI)  ? multiPtr->limitPTmax( event) : false;

    // Global recoil: reset counters and store locations of outgoing partons.
    timesPtr->prepareGlobal( event);
    bool isFirstTrial = true;

    // Set hard scale, maximum for showers and multiparton interactions.
    double pTscaleRad  = process.scale();
    double pTscaleMPI  = (doMPI && pTmaxMatchMPI == 3)
                       ? multiPtr->scaleLimitPT() : pTscaleRad;
    if (twoHard) {
      pTscaleRad       = max( pTscaleRad, process.scaleSecond() );
      pTscaleMPI       = min( pTscaleMPI, process.scaleSecond() );
    }
    double pTmaxMPI = (limitPTmaxMPI) ? pTscaleMPI : infoPtr->eCM();
    double pTmaxISR = (limitPTmaxISR) ? spacePtr->enhancePTmax() * pTscaleRad
                                      : infoPtr->eCM();
    double pTmaxFSR = (limitPTmaxFSR) ? timesPtr->enhancePTmax() * pTscaleRad
                                      : infoPtr->eCM();

    // Store the starting scale to use it for valence selection for gamma
    // beam. In case of MPIs use the pT scale of the last one.
    beamAPtr->pTMPI( process.scale() );
    beamBPtr->pTMPI( process.scale() );

    // Potentially reset starting scales for matrix element merging.
    if ( hasMergingHooks && (doTrial || canRemoveEvent || canRemoveEmission) )
      mergingHooksPtr->setShowerStartingScales( doTrial,
        (canRemoveEvent || canRemoveEmission), pTscaleRad, process, pTmaxFSR,
        limitPTmaxFSR, pTmaxISR, limitPTmaxISR, pTmaxMPI, limitPTmaxMPI );
    double pTmax    = max( pTmaxMPI, max( pTmaxISR, pTmaxFSR) );
    pTsaveMPI       = pTmaxMPI;
    pTsaveISR       = pTmaxISR;
    pTsaveFSR       = pTmaxFSR;

    // Prepare the classes to begin the generation.
    if (doMPI) multiPtr->prepare( event, pTmaxMPI, (iHardDiffLoop == 2) );
    if (doISR) spacePtr->prepare( 0, event, limitPTmaxISR);
    if (doFSRduringProcess) timesPtr->prepare( 0, event, limitPTmaxFSR);
    if (twoHard && doISR) spacePtr->prepare( 1, event, limitPTmaxISR);
    if (twoHard && doFSRduringProcess) timesPtr->prepare( 1, event,
       limitPTmaxFSR);

    // Impact parameter has now been chosen, but not usefully for diffraction.
    // (Never chosen for low-mass diffraction, twice for double diffraction.)
    if (!isDiff) infoPtr->setImpact( multiPtr->bMPI(), multiPtr->enhanceMPI(),
      multiPtr->enhanceMPIavg(), true, (iHardDiffLoop == 2) );

    // Set up initial veto scale.
    doVeto        = false;
    double pTveto = pTvetoPT;
    typeLatest    = 0;

    // Begin evolution down in pT from hard pT scale.
    do {
      infoPtr->addCounter(22);
      typeVetoStep = 0;
      nRad         =  nISR + nFSRinProc;

      // Check whether the beam photon has unresolved during the evolution.
      bool unresolvedGammaA = (beamAPtr->isGamma()
        && !(beamAPtr->resolvedGamma()) );
      bool unresolvedGammaB = (beamBPtr->isGamma()
        && !(beamBPtr->resolvedGamma()) );
      bool unresolvedGamma  = (unresolvedGammaA || unresolvedGammaB)
        || gammaModeEvent == 4;

      // Find next pT value for FSR, MPI and ISR.
      // Order calls to minimize time expenditure.
      double pTgen = 0.;

      // Potentially increase shower stopping scale for trial showers, to
      // avoid accumulating low-pT emissions (and weights thereof)
      if ( hasMergingHooks && doTrial)
        pTgen = max( pTgen, mergingHooksPtr->getShowerStoppingScale() );

      double pTtimes = (doFSRduringProcess)
        ? timesPtr->pTnext( event, pTmaxFSR, pTgen, isFirstTrial, doTrial)
        : -1.;
      pTgen = max( pTgen, pTtimes);
      // No MPIs for unresolved photons.
      double pTmulti = (doMPI && !unresolvedGamma)
        ? multiPtr->pTnext( pTmaxMPI, pTgen, event) : -1.;
      pTgen = max( pTgen, pTmulti);
      double pTspace = (doISR)
        ? spacePtr->pTnext( event, pTmaxISR, pTgen, nRad, doTrial) : -1.;
      double pTnow = max( pTtimes, max( pTmulti, pTspace));

      // Update information.
      infoPtr->setPTnow( pTnow);
      isFirstTrial = false;

      // Allow a user veto. Only do it once, so remember to change pTveto.
      if (pTveto > 0. && pTveto > pTnow) {
        pTveto = -1.;
        doVeto = userHooksPtr->doVetoPT( typeLatest, event);
        // Abort event if vetoed.
        if (doVeto) {
          if (isDiff) leaveResolvedDiff( iHardLoop, process, event);
          if (beamHasResGamma)
            leaveResolvedLeptonGamma( process, event, false);
          return false;
        }
      }

      // Do a multiparton interaction (if allowed).
      if (pTmulti > 0. && pTmulti > pTspace && pTmulti > pTtimes) {
        infoPtr->addCounter(23);
        if (multiPtr->scatter( event)) {
          typeLatest = 1;
          ++nMPI;
          if (canVetoMPIStep && nMPI <= nVetoMPIStep) typeVetoStep = 1;

          // Break for hard diffraction with MPI veto.
          if (isHardDiff && sampleTypeDiff == 4 && iHardDiffLoop == 1) {
            infoPtr->setHardDiff( false, false, false, false, 0., 0., 0., 0.);
            doDiffVeto = true;
            // Reset beam pointers and set event back to original frame.
            if (beamHasResGamma)
              leaveResolvedLeptonGamma( process, event, false);
            return false;
          }

          // Update ISR and FSR dipoles.
          if (doISR)              spacePtr->prepare( nMPI - 1, event);
          if (doFSRduringProcess) timesPtr->prepare( nMPI - 1, event);
          nBranch++;
          pTLastBranch = pTmulti;
          typeLastBranch = 1;
        }

        // Set maximal scales for next pT to pick.
        pTmaxMPI = pTmulti;
        pTmaxISR = min(pTmulti, pTmaxISR);
        pTmaxFSR = min(pTmulti, pTmaxFSR);
        pTmax    = pTmulti;
      }

      // Do an initial-state emission (if allowed).
      else if (pTspace > 0. && pTspace > pTtimes) {
        infoPtr->addCounter(24);

        // If MPIs, construct the gamma->qqbar branching in beamRemnants.
        if (spacePtr->branch( event)
            && ( !(nMPI > 1 && spacePtr->wasGamma2qqbar()) ) ) {
          typeLatest = 2;
          iSysNow = spacePtr->system();
          ++nISR;
          if (iSysNow == 0) ++nISRhard;
          if (canVetoStep && iSysNow == 0 && nISRhard <= nVetoStep)
            typeVetoStep = 2;

          // Update FSR dipoles.
          if (doFSRduringProcess) timesPtr->update( iSysNow, event,
            spacePtr->getHasWeaklyRadiated());
          nBranch++;
          pTLastBranch = pTspace;
          typeLastBranch = 2;

        // Rescatter: it is possible for kinematics to fail, in which
        //            case we need to restart the parton level processing.
        } else if (spacePtr->doRestart()) {
          physical = false;
          break;
        }

        // Set maximal scales for next pT to pick.
        pTmaxMPI = min( min(pTspace,pTmaxISR), pTmaxMPI);
        pTmaxISR = min(pTspace,pTmaxISR);
        pTmaxFSR = min( min(pTspace,pTmaxISR), pTmaxFSR);
        pTmax    = pTspace;
      }

      // Do a final-state emission (if allowed).
      else if (pTtimes > 0.) {
        infoPtr->addCounter(25);
        if (timesPtr->branch( event, true)) {
          typeLatest = 3;
          iSysNow = timesPtr->system();
          ++nFSRinProc;
          if (iSysNow == 0) ++nFSRhard;
          if (canVetoStep && iSysNow == 0 && nFSRhard <= nVetoStep)
            typeVetoStep = 3;

          // Update ISR dipoles.
          if (doISR) spacePtr->update( iSysNow, event,
            timesPtr->getHasWeaklyRadiated());
          nBranch++;
          pTLastBranch = pTtimes;
          typeLastBranch = 3;

        }

        // Set maximal scales for next pT to pick.
        pTmaxMPI = min( min(pTtimes,pTmaxFSR), pTmaxMPI);
        pTmaxISR = min( min(pTtimes,pTmaxFSR), pTmaxISR);
        pTmaxFSR = min(pTtimes, pTmaxFSR);
        pTmax    = pTtimes;
      }

      // If no pT scales above zero then nothing to be done.
      else pTmax = 0.;

      // Check for double counting for Drell-Yan weak production.
      // Only look at the second emission.
      if ( (infoPtr->code() == 221 || infoPtr->code() == 222) &&
            nISRhard + nFSRhard == 2 && vetoWeakJets) {
        int id1 = event[partonSystemsPtr->getOut(0,0)].id();
        int id2 = event[partonSystemsPtr->getOut(0,1)].id();
        int id3 = event[partonSystemsPtr->getOut(0,2)].id();
        Vec4 p1 = event[partonSystemsPtr->getOut(0,0)].p();
        Vec4 p2 = event[partonSystemsPtr->getOut(0,1)].p();
        Vec4 p3 = event[partonSystemsPtr->getOut(0,2)].p();

        // Make sure id1 is weak boson, and check that there
        // only is a single weak boson and no photons.
        bool doubleCountEvent = true;
        if (abs(id1) == 24 || abs(id1) == 23) {
          if (abs(id2) > 21 || abs(id3) > 21)
            doubleCountEvent = false;
        } else if (abs(id2) == 24 || abs(id2) == 23) {
          swap(id1,id2);
          swap(p1,p2);
          if (abs(id3) > 21)
            doubleCountEvent = false;
        } else if ( abs(id3) == 24 || abs(id3) == 23) {
          swap(id1,id3);
          swap(p1,p3);
        }

        if (doubleCountEvent) {
          double d = p1.pT2();
          bool cut = true;
          if (p2.pT2() < d) {d = p2.pT2(); cut = false;}
          if (p3.pT2() < d) {d = p3.pT2(); cut = false;}

          // Check for angle between weak boson and quarks.
          // (require final state particle to be a fermion)
          if (abs(id2) < 20) {
            double dij = min(p1.pT2(),p2.pT2())
              * pow2(RRapPhi(p1,p2)) / vetoWeakDeltaR2;
            if (dij < d) {
              d = dij;
              cut = true;
            }
          }

          if (abs(id3) < 20) {
            double dij = min(p1.pT2(),p3.pT2())
              * pow2(RRapPhi(p1,p3)) / vetoWeakDeltaR2;
            if (dij < d) {
              d = dij;
              cut = true;
            }
          }

          // Check for angle between recoiler and radiator,
          // if it is a quark anti-quark pair
          // or if the recoiler is a gluon.
          if (abs(id2) == 21 || abs(id3) == 21 || id2 == - id3) {
            double dij = min(p2.pT2(),p3.pT2())
              * pow2(RRapPhi(p2,p3)) / vetoWeakDeltaR2;
            if (dij < d) {
              d = dij;
              cut = false;
            }
          }

          // Veto event if it does not belong to Drell-Yan production.
          if (cut) return false;
        }
      }

      // Optionally check for a veto after the first few interactions,
      // or after the first few emissions, ISR or FSR, in the hardest system.
      if (typeVetoStep == 1) {
        doVeto = userHooksPtr->doVetoMPIStep( nMPI, event);
      } else if (typeVetoStep > 1) {
        doVeto = userHooksPtr->doVetoStep( typeVetoStep, nISRhard,
          nFSRhard, event);
      }

      // Abort event if vetoed.
      if (doVeto) {
        if (isDiff) leaveResolvedDiff( iHardLoop, process, event);
        if (beamHasResGamma) leaveResolvedLeptonGamma( process, event);
        return false;
      }

      // Keep on evolving until nothing is left to be done.
      if (typeLatest > 0 && typeLatest < 4)
        infoPtr->addCounter(25 + typeLatest);
      if (!isDiff) infoPtr->setPartEvolved( nMPI, nISR);

      // Handle potential merging veto.
      if ( canRemoveEvent && nISRhard + nFSRhard == 1 ) {
        // Simply check, and possibly reset weights.
        mergingHooksPtr->doVetoStep( process, event );
      }

    // End loop evolution down in pT from hard pT scale.
    } while (pTmax > 0.  && (nBranchMax <= 0 || nBranch < nBranchMax) );

    // Do all final-state emissions if not already considered above.
    if (doFSRafterProcess && (nBranchMax <= 0 || nBranch < nBranchMax) ) {

      // Find largest scale for final partons.
      pTmax = 0.;
      for (int i = 0; i < event.size(); ++i)
        if (event[i].isFinal() && event[i].scale() > pTmax)
          pTmax = event[i].scale();
      pTsaveFSR = pTmax;

      // Prepare all subsystems for evolution.
      for (int iSys = 0; iSys < partonSystemsPtr->sizeSys(); ++iSys)
        timesPtr->prepare( iSys, event);

      // Set up initial veto scale.
      doVeto = false;
      pTveto = pTvetoPT;

      // Begin evolution down in pT from hard pT scale.
      do {
        infoPtr->addCounter(29);
        typeVetoStep = 0;
        double pTtimes = timesPtr->pTnext( event, pTmax, 0.);
        infoPtr->setPTnow( pTtimes);

        // Allow a user veto. Only do it once, so remember to change pTveto.
        if (pTveto > 0. && pTveto > pTtimes) {
          pTveto = -1.;
          doVeto = userHooksPtr->doVetoPT( 4, event);
          // Abort event if vetoed.
          if (doVeto) {
            if (isDiff) leaveResolvedDiff( iHardLoop, process, event);
            if (beamHasResGamma) leaveResolvedLeptonGamma( process, event);
            return false;
          }
        }

        // Do a final-state emission (if allowed).
        if (pTtimes > 0.) {
          infoPtr->addCounter(30);
          if (timesPtr->branch( event, true)) {
            iSysNow = timesPtr->system();
            ++nFSRinProc;
            if (iSysNow == 0) ++nFSRhard;
            if (canVetoStep && iSysNow == 0 && nFSRhard <= nVetoStep)
            typeVetoStep = 4;

            nBranch++;
            pTLastBranch = pTtimes;
            typeLastBranch = 4;

          }
          pTmax = pTtimes;
        }

        // If no pT scales above zero then nothing to be done.
        else pTmax = 0.;

        // Optionally check for a veto after the first few emissions.
        if (typeVetoStep > 0) {
          doVeto = userHooksPtr->doVetoStep( typeVetoStep, nISRhard,
            nFSRhard, event);
          // Abort event if vetoed.
          if (doVeto) {
            if (isDiff) leaveResolvedDiff( iHardLoop, process, event);
            if (beamHasResGamma) leaveResolvedLeptonGamma( process, event);
            return false;
          }
        }

        // Handle potential merging veto.
        if ( canRemoveEvent && nISRhard + nFSRhard == 1 ) {
          // Simply check, and possibly reset weights.
          mergingHooksPtr->doVetoStep( process, event );
        }

        // Keep on evolving until nothing is left to be done.
        infoPtr->addCounter(31);

      } while (pTmax > 0.  && (nBranchMax <= 0 || nBranch < nBranchMax) );
    }

    // Handle veto after ISR + FSR + MPI, but before beam remnants
    // and resonance decays, e.g. for MLM matching.
    if (canVetoEarly && userHooksPtr->doVetoPartonLevelEarly( event)) {
      doVeto = true;
      if (isDiff) leaveResolvedDiff( iHardLoop, process, event);
      if (beamHasResGamma) leaveResolvedLeptonGamma( process, event);
      return false;
    }

    // Perform showers in resonance decay chains before beams & reconnection.
    if (earlyResDec) {
      int oldSizeEvt = event.size();
      int oldSizeSys = partonSystemsPtr->sizeSys();
      if (nBranchMax <= 0 || nBranch < nBranchMax)
        doVeto = !resonanceShowers( process, event, true);
      // Abort event if vetoed.
      if (doVeto) return false;

      // Reassign new decay products to original system, or else beam remnant
      // handling would be confused by systems without incoming partons.
      for (int iSys = oldSizeSys; iSys < partonSystemsPtr->sizeSys(); ++iSys)
        for (int iOut = 0; iOut < partonSystemsPtr->sizeOut(iSys); ++iOut)
          partonSystemsPtr->addOut(0, partonSystemsPtr->getOut( iSys, iOut) );
      partonSystemsPtr->setSizeSys( oldSizeSys);

      // Perform decays and showers of W and Z emitted in shower.
      // To do: check if W/Z emission is on in ISR or FSR??
      if (!wzDecayShowers( event)) return false;

      // User hook to reconnect colours specifically in resonance decays.
      if (canReconResSys && !userHooksPtr->doReconnectResonanceSystems(
        oldSizeEvt, event)) return false;
    }

    // Find the first particle in the current diffractive system.
    int  iFirst = 0;
    if (isDiff) {
      doDiffCR = isDiff;
      iFirst   = (iHardLoop == 1) ? 5 + sizeEvent - sizeProcess : sizeEvent;
      if (isDiffC) iFirst = 6 + sizeEvent - sizeProcess;
    }

    // Change the first particle for hard diffraction.
    if (infoPtr->hasPomPsystem()) {
      doDiffCR = true;
      iFirst   = 5;
    }

    // Add beam remnants, including primordial kT kick and colour tracing.
    if (!doTrial && physical && doRemnants
      && (!beamHasGamma || gammaModeEvent != 4)
      && !remnants.add( event, iFirst, doDiffCR)) physical = false;

    // If no problems then done.
    if (physical) break;

    // Else restore and loop, but do not throw existing diffractive system.
    if (!isDiff) event.clear();
    else {
      event.restoreSize();
      event.restoreJunctionSize();
    }
    beamAPtr->clear();
    beamBPtr->clear();
    partonSystemsPtr->clear();

    // Restore also the lepton beams if include photons.
    if (beamAhasResGamma) beamHadAPtr->clear();
    if (beamBhasResGamma) beamHadBPtr->clear();
    if (infoPtr->isVMDstateA()) beamVMDAPtr->clear();
    if (infoPtr->isVMDstateB()) beamVMDBPtr->clear();

  // End loop over ten tries. Restore from diffraction. Hopefully it worked.
  }
  if (isDiff) leaveResolvedDiff( iHardLoop, process, event);

  if (!physical) {
    // Leave hard diffractive system properly if beam remnant failed.
    if (infoPtr->hasPomPsystem()) leaveHardDiff( process, event, false);
    // Leave also photon from lepton framework.
    if (beamHasGamma) leaveResolvedLeptonGamma( process, event, false);
    return false;
  }

  // End big outer loop to handle two systems in double diffraction.
  }

  // If no additional MPI has been found then set up the diffractive
  // system the first time around.
  if (isHardDiff && sampleTypeDiff%2 == 0 && iHardDiffLoop == 1 && nMPI == 1) {
    event.clear();
    beamAPtr->clear();
    beamBPtr->clear();
    partonSystemsPtr->clear();
    setupHardDiff( process);
    continue;
  }

  // Do colour reconnection for non-diffractive events before resonance decays.
  if (doReconnect && !doDiffCR && reconnectMode > 0) {
    Event eventSave = event;
    bool colCorrect = false;
    for (int i = 0; i < 10; ++i) {
      colourReconnection.next(event, 0);
      if (junctionSplitting.checkColours(event)) {
        colCorrect = true;
        break;
      }
      else event = eventSave;
    }
    if (!colCorrect) {
      infoPtr->errorMsg("Error in PartonLevel::next: "
        "Colour reconnection failed.");
      return false;
    }
  }

  // Perform showers in resonance decay chains after beams & reconnection.
  int oldSizeEvt = event.size();
  if (!earlyResDec) {
    if (nBranchMax <= 0 || nBranch < nBranchMax)
      doVeto = !resonanceShowers( process, event, true);
    // Abort event if vetoed.
    if (doVeto) return false;

    // Perform decays and showers of W and Z emitted in shower.
    // To do:check if W/Z emission is on in ISR or FSR??
    if (!wzDecayShowers( event)) return false;

    // User hook to reconnect colours specifically in resonance decays.
    if (canReconResSys && !userHooksPtr->doReconnectResonanceSystems(
      oldSizeEvt, event)) return false;
  }

  // Store event properties. Not available for diffraction.
  if (!isDiff) infoPtr->setEvolution( pTsaveMPI, pTsaveISR, pTsaveFSR,
    nMPI, nISR, nFSRinProc, nFSRinRes);
  if (isDiff) {
    multiPtr->setEmpty();
    infoPtr->setImpact( multiPtr->bMPI(), multiPtr->enhanceMPI(),
      multiPtr->enhanceMPIavg(), false);
  }

  // Do colour reconnection for resonance decays.
  if (!earlyResDec && forceResonanceCR && doReconnect &&
      !doDiffCR && reconnectMode != 0) {
    Event eventSave = event;
    bool colCorrect = false;
    for (int i = 0; i < 10; ++i) {
      colourReconnection.next(event, oldSizeEvt);
      if (junctionSplitting.checkColours(event)) {
        colCorrect = true;
        break;
      }
      else event = eventSave;
    }
    if (!colCorrect) {
      infoPtr->errorMsg("Error in PartonLevel::next: "
        "Colour reconnection failed.");
      return false;
    }
  }

  // Leave diffractive events.
  if (isHardDiff) {

    // If inclusive sample wanted for MPI veto and nMPI > 1
    // then event is non-diffractive and we can break the loop.
    if (sampleTypeDiff == 2 && iHardDiffLoop == 1 && nMPI > 1) {
      infoPtr->setHardDiff( false, false, false, false, 0., 0., 0., 0.);
      break;
    }

    // Leave diffractive system properly.
    if (infoPtr->hasPomPsystem()) leaveHardDiff( process, event);
  }

  // End big outer loop to handle the setup of the diffractive system.
  }

  // If beam photon unresolved during evolution remove the copies of the
  // beam particle from the event record.
  if ( ( beamAPtr->isGamma() || beamBPtr->isGamma() )
       && ( !beamAPtr->resolvedGamma() || !beamBPtr->resolvedGamma() ) ) {
    if (!showUnresGamma && (gammaModeEvent != 4) ) cleanEventFromGamma( event);
  }

  // After parton level generation, add scattered leptons, restore the event.
  if (beamHasGamma) leaveResolvedLeptonGamma( process, event, true);

  // Allow option for QED shower phase after remnants but before hadronisation
  // Note: incoming leptons or photons may evolve as part of this phase.
  timesPtr->showerQEDafterRemnants(event);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Decide which diffractive subsystems are resolved (= perturbative).

int PartonLevel::decideResolvedDiff( Event& process) {

  // Loop over two systems.
  int nHighMass = 0;
  int iDSmin    = (isDiffC) ? 3 : 1;
  int iDSmax    = (isDiffC) ? 3 : 2;
  for (int iDSnow = iDSmin; iDSnow <= iDSmax; ++iDSnow) {

    // Offset for diffractive system when photons from lepton beams.
    int iDiffMot = iDSnow + 2 + gammaOffset;

    // Only high-mass diffractive systems should be resolved.
    double mDiff = process[iDiffMot].m();
    bool isHighMass = ( mDiff > mMinDiff && rndmPtr->flat()
      < pMaxDiff * ( 1. - exp( -(mDiff - mMinDiff) / mWidthDiff ) ) );

    // Set outcome and done.
    if (isHighMass) ++nHighMass;
    if (iDSnow == 1) isResolvedA = isHighMass;
    if (iDSnow == 2) isResolvedB = isHighMass;
    if (iDSnow == 3) isResolvedC = isHighMass;
  }
  return nHighMass;

}

//--------------------------------------------------------------------------

// Set up an unresolved process, i.e. elastic or diffractive.

bool PartonLevel::setupUnresolvedSys( Event& process, Event& event) {

  // No hard scale in event.
  process.scale( 0.);

  // Copy particles from process to event.
  for (int i = 0; i < process.size(); ++ i) event.append( process[i]);

  // Loop to find diffractively excited beams.
  for (iDS = 1; iDS < 4; ++iDS)
  if ( (iDS == 1 && isDiffA && !isResolvedA)
    || (iDS == 2 && isDiffB && !isResolvedB)
    || (iDS == 3 && isDiffC && !isResolvedC) ) {
    int iBeam = iDS + 2 + gammaOffset;

    // Diffractive mass. Boost and rotation from diffractive system
    // rest frame, aligned along z axis, to event cm frame.
    // Or to photon-photon(hadron) cm frame with photoproduction.
    double mDiff  = process[iBeam].m();
    double m2Diff = mDiff * mDiff;
    Vec4 pDiffA   = (iDS == 1) ? process[1 + gammaOffset].p()
      : process[1 + gammaOffset].p() - process[3 + gammaOffset].p();
    Vec4 pDiffB   = (iDS == 2) ? process[2 + gammaOffset].p()
      : process[2 + gammaOffset].p() - process[4 + gammaOffset].p();
    RotBstMatrix MtoCM;
    MtoCM.fromCMframe( pDiffA, pDiffB);

    // Beam Particle used for flavour content kicked out by Pomeron.
    // Randomize for central diffraction; misses closed gluon loop case.
    bool beamSideA = (iDS == 1 || (iDS == 3 && rndmPtr->flat() < 0.5));
    BeamParticle* beamPtr = (beamSideA) ? beamAPtr    : beamBPtr;
    if (iDS == 3) beamPtr = (beamSideA) ? beamPomAPtr : beamPomBPtr;

    // Pick quark or gluon kicked out and flavour subdivision.
    beamPtr->newValenceContent();
    bool gluonIsKicked = beamPtr->pickGluon(mDiff);
    int id1 = beamPtr->pickValence();
    int id2 = beamPtr->pickRemnant();

    // Find flavour masses. Scale them down if too big.
    double m1 = particleDataPtr->constituentMass(id1);
    double m2 = particleDataPtr->constituentMass(id2);
    if (m1 + m2 > 0.5 * mDiff) {
      double reduce = 0.5 * mDiff / (m1 + m2);
      m1 *= reduce;
      m2 *= reduce;
    }

    // If quark is kicked out, then trivial kinematics in rest frame.
    if (!gluonIsKicked) {
      double pAbs = sqrt( pow2(m2Diff - m1*m1 - m2*m2)
        - pow2(2. * m1 * m2) ) / (2. * mDiff);
      if (!beamSideA) pAbs = -pAbs;
      double e1 = (m2Diff + m1*m1 - m2*m2) / (2. * mDiff);
      double e2 = (m2Diff + m2*m2 - m1*m1) / (2. * mDiff);
      Vec4 p1( 0., 0., -pAbs, e1);
      Vec4 p2( 0., 0.,  pAbs, e2);

      // Boost and rotate to event cm frame.
      p1.rotbst( MtoCM);
      p2.rotbst( MtoCM);

      // Set colours.
      int col1, acol1, col2, acol2;
      if (particleDataPtr->colType(id1) == 1) {
        col1  = event.nextColTag();
        acol1 = 0;
        col2  = 0;
        acol2 = col1;
      } else {
        col1  = 0;
        acol1 = event.nextColTag();
        col2  = acol1;
        acol2 = 0;
      }
     // Update process colours to stay in step.
      process.nextColTag();

      // Store partons of diffractive system and mark system decayed.
      int iDauBeg = event.append( id1, 24, iBeam, 0, 0, 0, col1, acol1,
        p1, m1);
      int iDauEnd = event.append( id2, 63, iBeam, 0, 0, 0, col2, acol2,
        p2, m2);
      event[iBeam].statusNeg();
      event[iBeam].daughters(iDauBeg, iDauEnd);

    // If gluon is kicked out: share momentum between two remnants.
    } else {
      double m2Sys, zSys, pxSys, pySys, mTS1, mTS2;
      zSys = beamPtr->zShare(mDiff, m1, m2);

      // Provide relative pT kick in remnant. Construct (transverse) masses.
      pxSys = beamPtr->pxShare();
      pySys = beamPtr->pyShare();
      mTS1  = m1*m1 + pxSys*pxSys + pySys*pySys;
      mTS2  = m2*m2 + pxSys*pxSys + pySys*pySys;
      m2Sys = mTS1 / zSys + mTS2 / (1. - zSys);

      // Momentum of kicked-out massless gluon in diffractive rest frame.
      double pAbs  = (m2Diff - m2Sys) / (2. * mDiff);
      double pLRem = (beamSideA) ? pAbs : -pAbs;
      Vec4 pG(  0., 0., -pLRem, pAbs);
      Vec4 pRem(0., 0.,  pLRem, mDiff - pAbs);

      // Momenta of the two beam remnant flavours. (Lightcone p+ = m_diff!)
      double e1 = 0.5 * (zSys * mDiff + mTS1 / (zSys * mDiff));
      double pL1 = 0.5 * (zSys * mDiff - mTS1 / (zSys * mDiff));
      if (!beamSideA) pL1 = -pL1;
      Vec4 p1(pxSys, pySys, pL1, e1);
      Vec4 p2 = pRem - p1;

      // Boost and rotate to event cm frame. Improve precision.
      pG.rotbst( MtoCM);
      p1.rotbst( MtoCM);
      p2.rotbst( MtoCM);
      pG.e( pG.pAbs());

      // Set colours.
      int colG, acolG, col1, acol1, col2, acol2;
      if (particleDataPtr->colType(id1) == 1) {
        col1  = event.nextColTag();
        acol1 = 0;
        colG  = event.nextColTag();
        acolG = col1;
        col2  = 0;
        acol2 = colG;
      } else {
        col1  = 0;
        acol1 = event.nextColTag();
        colG  = acol1;
        acolG = event.nextColTag();
        col2  = acolG;
        acol2 = 0;
      }
      // Update process colours to stay in step.
      process.nextColTag();
      process.nextColTag();

      // Store partons of diffractive system and mark system decayed.
      int iDauBeg = event.append( 21, 24, iBeam, 0, 0, 0, colG, acolG, pG, 0.);
      event.append( id1, 63, iBeam, 0, 0, 0, col1, acol1, p1, m1);
      int iDauEnd = event.append( id2, 63, iBeam, 0, 0, 0, col2, acol2,
        p2, m2);
      event[iBeam].statusNeg();
      event[iBeam].daughters(iDauBeg, iDauEnd);
    }

  // End loop over beams. Done.
  }
  return true;

}

//--------------------------------------------------------------------------

// Set up the hard process(es), excluding subsequent resonance decays.

void PartonLevel::setupHardSys( Event& process, Event& event) {

  // Incoming partons to hard process are stored in slots 3 and 4.
  int inS = 0;
  int inP = 3;
  int inM = 4;

  // Mother and last entry of diffractive system. Offset.
  int iDiffMot = iDS + 2;
  int iDiffDau = process.size() - 1;
  int nOffset  = sizeEvent - sizeProcess;

  // Corrected information for hard diffraction.
  if (infoPtr->hasPomPsystem()) {
    iDiffMot = (isHardDiffB) ? 4 : 3;
    inS      = iDiffMot;
    inP      = 7;
    inM      = 8;
  }

  // If photons inside leptons more entries in event.
  // Further offset if hard diffraction set up.
  int nGammaOffset     = 0;
  int nGammaDiffOffset = 0;
  if ( beamHasResGamma || (beamHasGamma && isGammaHadronDir) ) {
    nGammaOffset = 2;
    if (!isDiff) {
      inP += nGammaOffset;
      inM += nGammaOffset;
    } else {
      iDiffMot += nGammaOffset;
      inS      += nGammaOffset;
    }
    if ( (isHardDiffA || isHardDiffB) && hardDiffSet ) nGammaDiffOffset = 4;
  }

  // Resolved diffraction means more entries.
  if (isDiff) {
    inS   = iDiffMot;
    inP   = iDiffDau - 3;
    inM   = iDiffDau - 2;

    // Diffractively excited particle described as Pomeron-hadron beams.
    event[inS].statusNeg();
    event[inS].daughters( inP - 2 + nOffset, inM - 2 + nOffset);
  }

  // If two hard interactions then find where second begins.
  int iBeginSecond = process.size();
  if (twoHard) {
    iBeginSecond = 5;
    while (process[iBeginSecond].status() != -21) ++iBeginSecond;
  }

  // If incoming partons are massive then recalculate to put them massless.
  if (process[inP].m() != 0. || process[inM].m() != 0.) {
    double pPos = process[inP].pPos() + process[inM].pPos();
    double pNeg = process[inP].pNeg() + process[inM].pNeg();
    process[inP].pz( 0.5 * pPos);
    process[inP].e(  0.5 * pPos);
    process[inP].m(  0.);
    process[inM].pz(-0.5 * pNeg);
    process[inM].e(  0.5 * pNeg);
    process[inM].m(  0.);
  }

  // Add incoming hard-scattering partons to list in beam remnants.
  double x1 = process[inP].pPos() / process[inS].m();
  double x2 = process[inM].pNeg() / process[inS].m();

  // If photon inside gamma calculate x with respect to photon beams.
  // With diffraction the beam offsets already taken care of.
  if ( (beamHasResGamma || (beamHasGamma && isGammaHadronDir) )
      && !isDiff ) {
    if ( hasGammaA)
      beamHadAPtr->append( 3, process[3].id(), beamHadAPtr->xGamma() );
    if ( hasGammaB)
      beamHadBPtr->append( 4, process[4].id(), beamHadBPtr->xGamma() );
    x1 = process[inP].pPos() / ( process[3 + nGammaDiffOffset].p()
       + process[4 + nGammaDiffOffset].p() ).mCalc();
    x2 = process[inM].pNeg() / ( process[3 + nGammaDiffOffset].p()
       + process[4 + nGammaDiffOffset].p() ).mCalc();
  }
  beamAPtr->append( inP + nOffset, process[inP].id(), x1);
  beamBPtr->append( inM + nOffset, process[inM].id(), x2);

  // Scale. Find whether incoming partons are valence or sea. Store.
  // When an x-dependent matter profile is used with nonDiffractive,
  // trial interactions mean that the valence/sea choice has already
  // been made and should be restored here.
  double scale = process.scale();
  int vsc1, vsc2;
  beamAPtr->xfISR( 0, process[inP].id(), x1, scale*scale);
  if (isNonDiff && (vsc1 = multiPtr->getVSC1()) != 0)
    (*beamAPtr)[0].companion(vsc1);
  else if (beamAPtr->isGamma()) vsc1 = beamAPtr->gammaValSeaComp(0);
  else vsc1 = beamAPtr->pickValSeaComp();
  beamBPtr->xfISR( 0, process[inM].id(), x2, scale*scale);
  if (isNonDiff && (vsc2 = multiPtr->getVSC2()) != 0)
    (*beamBPtr)[0].companion(vsc2);
  else if (beamBPtr->isGamma()) vsc2 = beamBPtr->gammaValSeaComp(0);
  else vsc2 = beamBPtr->pickValSeaComp();
  bool isVal1 = (vsc1 == -3);
  bool isVal2 = (vsc2 == -3);
  infoPtr->setValence( isVal1, isVal2);

  // Initialize info needed for subsequent sequential decays + showers.
  nHardDone = sizeProcess;
  iPosBefShow.resize( process.size() );
  fill( iPosBefShow.begin(), iPosBefShow.end(), 0);

  // Add the beam and hard subprocess partons to the event record.
  for (int i = sizeProcess; i < iBeginSecond; ++i) {
    if (process[i].mother1() > inM) break;
    int j = event.append(process[i]);
    iPosBefShow[i] = i;

    // Offset history if process and event not in step.
    if (nOffset != 0) {
      int iOrd = i - iBeginSecond + 7;
      if (iOrd == 1 || iOrd == 2)
        event[j].offsetHistory( 0, 0, 0, nOffset);
      else if (iOrd == 3 || iOrd == 4)
        event[j].offsetHistory( 0, nOffset, 0, nOffset);
      else if (iOrd == 5 || iOrd == 6)
        event[j].offsetHistory( 0, nOffset, 0, 0);
    }

    // Currently outgoing ones should not count as decayed.
    if (event[j].status() == -22) {
      event[j].statusPos();
      event[j].daughters(0, 0);
    }

    // Complete task of copying hard subsystem into event record.
    ++nHardDone;
  }

  // Store participating partons as first set in list of all systems.
  partonSystemsPtr->addSys();
  partonSystemsPtr->setInA(0, inP + nOffset);
  partonSystemsPtr->setInB(0, inM + nOffset);
  for (int i = inM + 1; i < nHardDone; ++i)
    partonSystemsPtr->addOut(0, i + nOffset);
  partonSystemsPtr->setSHat( 0,
    (event[inP + nOffset].p() + event[inM + nOffset].p()).m2Calc() );
  partonSystemsPtr->setPTHat( 0, scale);

  // Identify second hard process where applicable.
  // Since internally generated incoming partons are guaranteed massless.
  if (twoHard) {
    int inP2 = iBeginSecond;
    int inM2 = iBeginSecond + 1;

    // Add incoming hard-scattering partons to list in beam remnants.
    // Not valid if not in rest frame??
    x1 = process[inP2].pPos() / process[0].e();
    beamAPtr->append( inP2, process[inP2].id(), x1);
    x2 = process[inM2].pNeg() / process[0].e();
    beamBPtr->append( inM2, process[inM2].id(), x2);

    // Find whether incoming partons are valence or sea.
    scale = process.scaleSecond();
    beamAPtr->xfISR( 1, process[inP2].id(), x1, scale*scale);
    beamAPtr->pickValSeaComp();
    beamBPtr->xfISR( 1, process[inM2].id(), x2, scale*scale);
    beamBPtr->pickValSeaComp();

    // Add the beam and hard subprocess partons to the event record.
    for (int i = inP2; i < process.size(); ++ i) {
      int mother = process[i].mother1();
      if ( (mother > 2 && mother < inP2) || mother > inM2 ) break;
      event.append(process[i]);
      iPosBefShow[i] = i;

      // Currently outgoing ones should not count as decayed.
      if (event[i].status() == -22) {
        event[i].statusPos();
        event[i].daughters(0, 0);
      }

      // Complete task of copying hard subsystem into event record.
      ++nHardDone;
    }

    // Store participating partons as second set in list of all systems.
    partonSystemsPtr->addSys();
    partonSystemsPtr->setInA(1, inP2);
    partonSystemsPtr->setInB(1, inM2);
    for (int i = inM2 + 1; i < nHardDone; ++i)
      partonSystemsPtr->addOut(1, i);
    partonSystemsPtr->setSHat( 1,
      (event[inP2].p() + event[inM2].p()).m2Calc() );
    partonSystemsPtr->setPTHat( 1, scale);

  // End code for second hard process.
  }

  // Update event colour tag to maximum in whole process.
  int maxColTag = 0;
  for (int i = 0; i < process.size(); ++ i) {
    if (process[i].col() > maxColTag) maxColTag = process[i].col();
    if (process[i].acol() > maxColTag) maxColTag = process[i].acol();
  }
  event.initColTag(maxColTag);

  // Copy junctions from process to event.
  for (int iJun = 0; iJun < process.sizeJunction(); ++iJun) {
    // Resonance decay products may not have been copied from process to
    // event yet. If so, do not add junctions associated with decays yet.
    int kindJunction = process.kindJunction(iJun);
    bool doCopy = true;
    // For junction types <= 4, check if final-state legs were copied.
    if (kindJunction <= 4) {
      int iLegF1 = (kindJunction - 1) / 2;
      for (int iLeg = iLegF1; iLeg <= 2; ++iLeg) {
        bool colFound = false;
        for (int i = inM + 1; i < event.size(); ++i) {
          int col = (kindJunction % 2 == 1) ? event[i].col() : event[i].acol();
          if (col == process.colJunction(iJun,iLeg)) colFound = true;
        }
        if (!colFound) doCopy = false;
      }
    }
    if (doCopy) {
      event.appendJunction( process.getJunction(iJun));
    }
  }

  // Done.
}

//--------------------------------------------------------------------------

// Set up the event for subsequent resonance decays and showers.

void PartonLevel::setupShowerSys( Event& process, Event& event) {

  // Reset event record to only contain line 0.
  event.clear();
  event.append( process[0]);

  // Initialize info needed for subsequent sequential decays + showers.
  nHardDone = 1;
  iPosBefShow.resize( process.size());
  fill( iPosBefShow.begin(), iPosBefShow.end(), 0);

  // Add the hard subprocess partons to the event record.
  for (int i = 1; i < process.size(); ++i) {
    if (process[i].mother1() > 0) break;
    int j = event.append(process[i]);
    iPosBefShow[i] = i;

    // Currently outgoing ones should not count as decayed.
    if (event[j].status() == -22) {
      event[j].statusPos();
      event[j].daughters(0, 0);
    }

    // Complete task of copying hard subsystem into event record.
    ++nHardDone;
  }

  // Store participating partons as first set in list of all systems.
  partonSystemsPtr->clear();
  partonSystemsPtr->addSys();
  for (int i = 1; i < nHardDone; ++i) partonSystemsPtr->addOut(0, i);
  partonSystemsPtr->setSHat( 0, pow2(process[0].m()) );
  partonSystemsPtr->setPTHat( 0, 0.5 * process[0].m() );

  // Copy junctions from process to event.
  for (int iJun = 0; iJun < process.sizeJunction(); ++iJun) {
    // Resonance decay products may not have been copied from process to
    // event yet. If so, do not add junctions associated with decays yet.
    int kindJunction = process.kindJunction(iJun);
    bool doCopy = true;
    // For junction types <= 4, check if final-state legs were copied.
    if (kindJunction <= 4) {
      int iLegF1 = (kindJunction - 1) / 2;
      for (int iLeg = iLegF1; iLeg <= 2; ++iLeg) {
        bool colFound = false;
        for (int i = 1; i < event.size(); ++i) {
          int col = (kindJunction % 2 == 1) ? event[i].col() : event[i].acol();
          if (col == process.colJunction(iJun,iLeg)) colFound = true;
        }
        if (!colFound) doCopy = false;
      }
    }
    if (doCopy) {
      event.appendJunction( process.getJunction(iJun));
    }
  }

  // Done.
}

//--------------------------------------------------------------------------

// Resolved diffraction: replace full event with diffractive subsystem.

void PartonLevel::setupResolvedDiff( Event& process) {

  // Mother and last entry of diffractive system.
  int iDiffMot     = iDS + 2 + gammaOffset;
  int iDiffDau     = process.size() - 1;

  // Diffractively excited particle to be replaced by Pomeron-hadron beams
  // (or Pomeron-Pomeron beams for central diffraction).
  process[iDiffMot].statusNeg();
  process[iDiffMot].daughters( iDiffDau + 1, iDiffDau + 2);

  // Diffractive system mass.
  double mDiff   = process[iDiffMot].m();
  double m2Diff  = mDiff * mDiff;

  // Set up Pomeron-proton or Pomeron-Pomeron system as if it were
  // the complete collision. Set Pomeron "beam particle" massless.
  int idDiffA    = (iDS == 1) ? process[1 + gammaOffset].id() : 990;
  int idDiffB    = (iDS == 2) ? process[2 + gammaOffset].id() : 990;
  double mDiffA  = (iDS == 1) ? process[1 + gammaOffset].m() : 0.;
  double mDiffB  = (iDS == 2) ? process[2 + gammaOffset].m() : 0.;
  if (idDiffA == 22 && infoPtr->isVMDstateA()) {
    idDiffA = (iDS == 1) ? infoPtr->idVMDA() : 990;
    mDiffA  = (iDS == 1) ? infoPtr->mVMDA()  : 0.;
  }
  if (idDiffB == 22 && infoPtr->isVMDstateB()) {
    idDiffB = (iDS == 2) ? infoPtr->idVMDB() : 990;
    mDiffB  = (iDS == 2) ? infoPtr->mVMDB()  : 0.;
  }
  double m2DiffA = mDiffA * mDiffA;
  double m2DiffB = mDiffB * mDiffB;
  double eDiffA  = 0.5 * (m2Diff + m2DiffA - m2DiffB) / mDiff;
  double eDiffB  = 0.5 * (m2Diff + m2DiffB - m2DiffA) / mDiff;
  double pzDiff  = 0.5 * sqrtpos( pow2(m2Diff - m2DiffA - m2DiffB)
                 - 4. * m2DiffA * m2DiffB ) / mDiff;
  process.append( idDiffA, 13, iDiffMot, 0, 0, 0, 0, 0,
                   0., 0.,  pzDiff, eDiffA, mDiffA);
  process.append( idDiffB, 13, iDiffMot, 0, 0, 0, 0, 0,
                   0., 0., -pzDiff, eDiffB, mDiffB);

  // Reassign beam pointers to refer to subsystem effective beams.
  // If VMD state in gamma, then use VMD instead of gamma as beam pointer.
  beamAPtr       = (iDS == 1) ? beamHadAPtr : beamPomAPtr;
  beamBPtr       = (iDS == 2) ? beamHadBPtr : beamPomBPtr;
  if (infoPtr->isVMDstateA())
    beamAPtr = (iDS == 1) ? beamVMDAPtr : beamPomAPtr;
  if (infoPtr->isVMDstateB())
    beamBPtr = (iDS == 2) ? beamVMDBPtr : beamPomBPtr;

  // Pretend that the diffractive system is the whole collision.
  eCMsave = infoPtr->eCM();
  infoPtr->setECM( mDiff);
  beamAPtr->newPzE(  pzDiff, eDiffA);
  beamBPtr->newPzE( -pzDiff, eDiffB);

  // Keep track of pomeron momentum fraction.
  if ( beamAPtr->id() == 990 )
    beamAPtr->xPom(pow2(mDiff/eCMsave));
  if ( beamBPtr->id() == 990 )
    beamBPtr->xPom(pow2(mDiff/eCMsave));

  // Beams not found in normal slots 1 and 2.
  int beamOffset = (sizeEvent > 0) ? sizeEvent - 1 : 4;
  int beamRemnOffset = iDS;
  if (beamAPtr->isGamma() || beamBPtr->isGamma()) beamRemnOffset = 4;

  // Reassign beam pointers in other classes.
  timesPtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  timesDecPtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  spacePtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  remnants.reassignBeamPtrs(  beamAPtr, beamBPtr, beamRemnOffset);
  colourReconnection.reassignBeamPtrs(  beamAPtr, beamBPtr);


  // Reassign multiparton interactions pointer to right object.
  if      (iDS == 1) multiPtr = &multiSDA;
  else if (iDS == 2) multiPtr = &multiSDB;
  else               multiPtr = &multiCD;

}

//--------------------------------------------------------------------------

// Resolved diffraction: restore to original behaviour.

void PartonLevel::leaveResolvedDiff( int iHardLoop, Event& process,
  Event& event) {

  // Reconstruct boost and rotation to event cm frame.
  // In case of photoproduction boost to photon-photon(hadrom) cm frame.
  Vec4 pDiffA = (iDS == 1) ? process[1 + gammaOffset].p()
              : process[1 + gammaOffset].p() - process[3 + gammaOffset].p();
  Vec4 pDiffB = (iDS == 2) ? process[2 + gammaOffset].p()
              : process[2 + gammaOffset].p() - process[4 + gammaOffset].p();
  RotBstMatrix MtoCM;
  MtoCM.fromCMframe( pDiffA, pDiffB);

  // Perform rotation and boost on diffractive system.
  for (int i = sizeProcess; i < process.size(); ++i)
    process[i].rotbst( MtoCM);
  int iFirst = (iHardLoop == 1) ? 5 + sizeEvent - sizeProcess + gammaOffset
    : sizeEvent;
  if (isDiffC) iFirst = 6 + sizeEvent - sizeProcess;
  for (int i = iFirst; i < event.size(); ++i)
    event[i].rotbst( MtoCM);

  // Restore cm energy.
  infoPtr->setECM( eCMsave);
  beamAPtr->newPzE( event[1].pz(), event[1].e());
  beamBPtr->newPzE( event[2].pz(), event[2].e());
  // Keeping track of pomeron momentum fraction.
  beamAPtr->xPom();
  beamBPtr->xPom();

  // Restore beam pointers to incoming hadrons.
  beamAPtr = beamHadAPtr;
  beamBPtr = beamHadBPtr;

  // Reassign beam pointers in other classes.
  timesPtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  timesDecPtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  spacePtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  remnants.reassignBeamPtrs(  beamAPtr, beamBPtr, 0);
  colourReconnection.reassignBeamPtrs(  beamAPtr, beamBPtr);

  // Restore multiparton interactions pointer to default object.
  multiPtr = &multiMB;

}

//--------------------------------------------------------------------------

// Set up special handling of hard diffraction.

void PartonLevel::setupHardDiff( Event& process) {

  // Create a temporary event record holding the info of the hard process.
  Event tmpProcess = process;
  process.clear();
  process.scale(tmpProcess.scale());

  // Add the first three entries: system + incoming beams.
  for (int iEntry = 0; iEntry < 3; ++iEntry)
    process.append( tmpProcess[iEntry]);

  // Get system info and calculate diffractive system mass.
  double eCM      = infoPtr->eCM();
  double sNow     = eCM * eCM;
  double xPom     = (isHardDiffB) ? infoPtr->xPomeronA()
                  : infoPtr->xPomeronB();
  double phiPom   = 2. * M_PI * rndmPtr->flat();
  double thetaPom = (isHardDiffB) ? hardDiffraction.getThetaPomeronA()
                  :  hardDiffraction.getThetaPomeronB();
  double m2Diff   = xPom * sNow;
  double mDiff    = sqrt(m2Diff);

  // Add possible photon beams from leptons.
  if (beamAhasGamma || beamBhasGamma) {
    process.append( tmpProcess[3]);
    process.append( tmpProcess[4]);
  }

  // Particle masses and ids.
  // If isHardDiffB (or A) then m1 + m2 -> m1 + mDiff (or mDiff + m2).
  // Change an elastically scattered gamma to a rho, since a Pomeron only
  // can be taken from a VMD particle, not from the gamma proper.
  int id1   = process[1 + gammaOffset].id();
  int id2   = process[2 + gammaOffset].id();
  int idEl  = (isHardDiffB) ? id1 : id2;
  if (idEl == 22) idEl = 113;
  int idX   = (isHardDiffB) ? id2 : id1;
  double m1 = process[1 + gammaOffset].m();
  double m2 = process[2 + gammaOffset].m();
  double mEl = (isHardDiffB) ? m1 : m2;

  // Evaluate momenta of outgoing particles, initially along beam axis.
  // Special handling for rho with too high a Breit-Wigner-selected mass,
  // where kinematics failures leads to retries.
  double m3, m4, s3, s4, lambda34;
  if (idEl == 113) {
    int nTry = 0;
    do {
      mEl      = particleDataPtr->mSel(113);
      m3       = (isHardDiffB) ? mEl : mDiff;
      m4       = (isHardDiffA) ? mEl : mDiff;
      s3       = pow2(m3);
      s4       = pow2(m4);
      lambda34 = sqrtpos( pow2( sNow - s3 - s4) - 4. * s3 * s4 );
      ++nTry;
    } while (lambda34 <= 0. && nTry < 10);
  } else {
    m3       = (isHardDiffB) ? mEl : mDiff;
    m4       = (isHardDiffA) ? mEl : mDiff;
    s3       = pow2(m3);
    s4       = pow2(m4);
    lambda34 = sqrtpos( pow2( sNow - s3 - s4) - 4. * s3 * s4 );
  }
  double pAbs = 0.5 * lambda34 / eCM;
  Vec4 p3 = Vec4( 0., 0.,  pAbs, 0.5 * (sNow + s3 - s4) / eCM);
  Vec4 p4 = Vec4( 0., 0., -pAbs, 0.5 * (sNow + s4 - s3) / eCM);

  // Take copies for later longitudinal boost; then rotate outgoing beams.
  Vec4 pD3 = p3;
  Vec4 pD4 = p4;
  p3.rot( thetaPom, phiPom);
  p4.rot( thetaPom, phiPom);

  //  Append intermediate states to the event record.
  int status3 = (isHardDiffB) ? 14 : 15;
  int status4 = (isHardDiffA) ? 14 : 15;
  int sign    = (idX > 0) ? 1 : -1;
  int idDiff  = (idX == 22) ? 9900020 : sign * 9902210;
  int id3     = (isHardDiffB) ? idEl : idDiff;
  int id4     = (isHardDiffA) ? idEl : idDiff;
  process.append( id3, status3, 1 + gammaOffset, 0, 0, 0, 0, 0, p3, m3);
  process.append( id4, status4, 2 + gammaOffset, 0, 0, 0, 0, 0, p4, m4);

  // Correct event record history accordingly.
  process[1 + gammaOffset].daughters(3 + gammaOffset, 0);
  process[2 + gammaOffset].daughters(4 + gammaOffset, 0);
  int iDiffMot     = (isHardDiffB) ? 4 + gammaOffset: 3 + gammaOffset;
  int iDiffRad     = process.size() - 1;
  process[iDiffMot].statusNeg();
  process[iDiffMot].daughters( iDiffRad + 1, iDiffRad + 2);

  // Set up Pomeron-particle system as if it were the complete collision.
  // Set Pomeron "beam particle" massless.
  int idDiffA    = (isHardDiffB) ? 990 : process[1 + gammaOffset].id();
  int idDiffB    = (isHardDiffA) ? 990 : process[2 + gammaOffset].id();
  double mDiffA  = (isHardDiffB) ? 0. : process[1 + gammaOffset].m();
  double mDiffB  = (isHardDiffA) ? 0. : process[2 + gammaOffset].m();
  double m2DiffA = mDiffA * mDiffA;
  double m2DiffB = mDiffB * mDiffB;
  double eDiffA  = 0.5 * (m2Diff + m2DiffA - m2DiffB) / mDiff;
  double eDiffB  = 0.5 * (m2Diff + m2DiffB - m2DiffA) / mDiff;
  double pzDiff  = 0.5 * sqrtpos( pow2(m2Diff - m2DiffA - m2DiffB)
                 - 4. * m2DiffA * m2DiffB ) / mDiff;
  process.append( idDiffA, 13, iDiffMot, 0, 0, 0, 0, 0,
                   0., 0.,  pzDiff, eDiffA, mDiffA);
  process.append( idDiffB, 13, iDiffMot, 0, 0, 0, 0, 0,
                   0., 0., -pzDiff, eDiffB, mDiffB);

  // Append hard process.
  vector<int> hardParton;
  for (int iHard = 3 + gammaOffset; iHard < tmpProcess.size(); ++iHard)
    hardParton.push_back( process.append(tmpProcess[iHard]) );

  // Boost the hard partons in z-direction (from pp to Pp system).
  Vec4 pDiffA = (isHardDiffA) ? process[1 + gammaOffset].p()
                              : process[1 + gammaOffset].p() - pD3;
  Vec4 pDiffB = (isHardDiffB) ? process[2 + gammaOffset].p()
                              : process[2 + gammaOffset].p() - pD4;
  RotBstMatrix MtoCM;
  MtoCM.toCMframe( pDiffA, pDiffB);
  for (unsigned int i = 0; i < hardParton.size(); ++i)
    process[hardParton[i]].rotbst(MtoCM);

  // Change mothers and daughters after appending hard process.
  for (unsigned int j = 0; j < hardParton.size(); ++j) {
    int mother1 = (tmpProcess[j + 3 + gammaOffset].mother1() == 0)
      ? 0 : tmpProcess[j + 3 + gammaOffset].mother1() + 4;
    int mother2 = (tmpProcess[j + 3 + gammaOffset].mother2() == 0)
      ? 0 : tmpProcess[j + 3 + gammaOffset].mother2() + 4;
    int daughter1 = (tmpProcess[j + 3 + gammaOffset].daughter1() == 0)
      ? 0 : tmpProcess[j + 3 + gammaOffset].daughter1() + 4;
    int daughter2 = (tmpProcess[j + 3 + gammaOffset].daughter2() == 0)
      ? 0 : tmpProcess[j + 3 + gammaOffset].daughter2() + 4;
    process[hardParton[j]].mothers( mother1,mother2);
    process[hardParton[j]].daughters( daughter1, daughter2);
  }

  // Search for pomeron and particle with status codes 13 (beam-inside-beam).
  int iPomeron = 0;
  int iPrtcl   = 0;
  for (int i = 0; i < process.size(); ++i) {
    if (process[i].id() == 990 && process[i].status() == 13) iPomeron = i;
    if (process[i].idAbs() == idX && process[i].status() == 13) iPrtcl = i;
  }

  if (isHardDiffB) {
    process[iPomeron].daughters(hardParton[0], 0);
    process[iPrtcl].daughters(hardParton[1],0);
    process[hardParton[0]].mothers(iPomeron,0);
    process[hardParton[1]].mothers(iPrtcl, 0);
  } else {
    process[iPomeron].daughters(hardParton[1], 0);
    process[iPrtcl].daughters(hardParton[0],0);
    process[hardParton[1]].mothers(iPomeron,0);
    process[hardParton[0]].mothers(iPrtcl, 0);
  }

  // Negate status of Pomeron and proton
  process[iPomeron].statusNeg();
  process[iPrtcl].statusNeg();

  // Change state of system to unresolved to avoid aborting from Pythia.
  infoPtr->setHasUnresolvedBeams( true);

  // Reassign beam pointers to refer to subsystem effective beams.
  // If gamma-in-lepton change to beamGamPtr.
  beamAPtr = (isHardDiffB) ? beamPomAPtr
    : (beamAhasGamma ? beamGamAPtr : beamHadAPtr);
  beamBPtr = (isHardDiffA) ? beamPomBPtr
    : (beamBhasGamma ? beamGamBPtr : beamHadBPtr);

  // Pretend that the diffractive system is the whole collision.
  eCMsave = infoPtr->eCM();
  infoPtr->setECM( mDiff);
  beamAPtr->newPzE(  pzDiff, eDiffA);
  beamBPtr->newPzE( -pzDiff, eDiffB);

  // Beams not found in normal slots 1 and 2.
  int beamOffset = 4 + gammaOffset;
  int beamRemnOffset = (isHardDiffB) ? 2 : 1;
  if (beamAPtr->isGamma() || beamBPtr->isGamma())
    beamRemnOffset = 4 + gammaOffset;

  // Reassign beam pointers in other classes.
  timesPtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  timesDecPtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  spacePtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  remnants.reassignBeamPtrs(  beamAPtr, beamBPtr, beamRemnOffset);
  colourReconnection.reassignBeamPtrs(  beamAPtr, beamBPtr);

  // Reassign multiparton interactions pointer to right object.
  if      (isHardDiffA) multiPtr = &multiSDA;
  else if (isHardDiffB) multiPtr = &multiSDB;

  // Set the beam offset for MPIs.
  multiPtr->setBeamOffset(beamOffset);

  // Hard diffractive system set.
  hardDiffSet = true;

  // Done.
  infoPtr->setHasPomPsystem( true);

}

//--------------------------------------------------------------------------

// Leave special handling of hard diffraction.

void PartonLevel::leaveHardDiff( Event& process, Event& event,
  bool physical) {

  // Calculate kinematics only for physical events.
  if (physical) {

    // Reconstruct boost and rotation to event cm frame.
    Vec4 pDiffA = (isHardDiffA) ? process[1 + gammaOffset].p()
                : process[1 + gammaOffset].p() - process[3 + gammaOffset].p();
    Vec4 pDiffB = (isHardDiffB) ? process[2 + gammaOffset].p()
                : process[2 + gammaOffset].p() - process[4 + gammaOffset].p();
    RotBstMatrix MtoCM;
    MtoCM.fromCMframe( pDiffA, pDiffB);

    // Perform rotation and boost on diffractive system.
    for (int i = 5 + gammaOffset; i < process.size(); ++i)
      process[i].rotbst( MtoCM);
    for (int i = 5 + gammaOffset; i < event.size(); ++i)
      event[i].rotbst( MtoCM);

    // Reset beam energies.
    beamAPtr->newPzE( event[1 + gammaOffset].pz(), event[1 + gammaOffset].e());
    beamBPtr->newPzE( event[2 + gammaOffset].pz(), event[2 + gammaOffset].e());

  }

  // Clear diffractive info.
  isHardDiffA = isHardDiffB = isHardDiff = false;

  // Restore cm energy.
  infoPtr->setECM( eCMsave);

  // Restore beam pointers to incoming hadrons.
  // If gamma-in-lepton need to change to beamGamPtr.
  beamAPtr = beamAhasGamma ? beamGamAPtr : beamHadAPtr;
  beamBPtr = beamBhasGamma ? beamGamBPtr : beamHadBPtr;

  // Reassign beam pointers in other classes.
  timesPtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  timesDecPtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  spacePtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  remnants.reassignBeamPtrs(  beamAPtr, beamBPtr, 0);
  colourReconnection.reassignBeamPtrs(  beamAPtr, beamBPtr);

  // Reset the beam offset to normal.
  multiPtr->setBeamOffset(0);

  // Restore multiparton interactions pointer to default object.
  multiPtr = &multiMB;

}

//--------------------------------------------------------------------------

// Replace full event with photon-photon/hadron subsystem.
// Use the sampled kinematics to construct the momenta.

bool PartonLevel::setupResolvedLeptonGamma( Event& process) {

  // Save the collision energy of the lepton system.
  eCMsaveGamma = infoPtr->eCM();

  // Beams not found in normal slots 1 and 2 but 2 step forward.
  int beamOffset = 2;
  int iBeamA     = 3;
  int iBeamB     = 4;
  gammaOffset    = beamOffset;

  // Retrieve the information set on GammaKinematics. Use sHat for 2->1.
  double mGmGm  = (infoPtr->nFinal() > 1) || (gammaModeEvent != 4)
                ? infoPtr->eCMsub() : sqrt( infoPtr->sHat());
  double m2GmGm = pow2(mGmGm);

  // Massless photons here, virtualities added after parton level evolution.
  double m2Gamma1 = hasGammaA ? 0. : pow2(beamAPtr->m());
  double m2Gamma2 = hasGammaB ? 0. : pow2(beamBPtr->m());

  // Derive the new momenta in the CM frame of the gamma-gamma system.
  double eGamA    = 0.5 * (m2GmGm + m2Gamma1 - m2Gamma2) / mGmGm;
  double eGamB    = 0.5 * (m2GmGm + m2Gamma2 - m2Gamma1) / mGmGm;
  double pzGam    = 0.5 * sqrtpos( pow2(m2GmGm - m2Gamma1 - m2Gamma2)
                  - 4. * m2Gamma1 * m2Gamma2 ) / mGmGm;
  Vec4 pGammaANew( 0, 0,  pzGam, eGamA);
  Vec4 pGammaBNew( 0, 0, -pzGam, eGamB);

  // Set the beam momenta to new rest frame of gamma-gamma.
  beamGamAPtr->newPzE(  pzGam, eGamA);
  beamGamBPtr->newPzE( -pzGam, eGamB);

  // Save the original photon momenta.
  Vec4 pGammaA = process[iBeamA].p();
  Vec4 pGammaB = process[iBeamB].p();

  // Boost the process to gamma-gamma rest frame (no kT added yet).
  // For diffractive and elastic boost only final state particles.
  RotBstMatrix MtoGammaGamma;
  MtoGammaGamma.toCMframe( pGammaA, pGammaB);
  if (!isDiff && !isElastic) process.rotbst(MtoGammaGamma);
  else {
    for (int i = 0; i < 5; ++i) process[i].rotbst(MtoGammaGamma);
  }

  // Set momenta of photons according to new m2GmGm.
  process[iBeamA].p(pGammaANew);
  process[iBeamB].p(pGammaBNew);

  // If another beam not a photon set the mass as well.
  if ( !hasGammaA && !(beamBPtr->getGammaMode() == 2) )
    process[iBeamA].m( sqrt(m2Gamma1) );
  if ( !hasGammaB && !(beamAPtr->getGammaMode() == 2) )
    process[iBeamB].m( sqrt(m2Gamma2) );

  // Done for direct-direct and elastic processes since no need to reassign
  // beams.
  if ( (gammaModeEvent == 4) || isElastic ) return true;

  // Copy sampled VMD states to photon beams when necessary.
  if (infoPtr->isVMDstateA() ) beamGamAPtr->setVMDstate(true,
    infoPtr->idVMDA(), infoPtr->mVMDA(), infoPtr->scaleVMDA(), false);
  if (infoPtr->isVMDstateB() ) beamGamBPtr->setVMDstate(true,
    infoPtr->idVMDB(), infoPtr->mVMDB(), infoPtr->scaleVMDB(), false);

  // Reassign beam pointers to refer to subsystem effective beams.
  if ( hasGammaA) beamAPtr = beamGamAPtr;
  else beamAPtr->newPzE(  pzGam, eGamA);
  if ( hasGammaB) beamBPtr = beamGamBPtr;
  else beamBPtr->newPzE( -pzGam, eGamB);

  // Change state of system to unresolved to avoid aborting from Pythia.
  if ( (beamAhasResGamma && (!beamBhasResGamma && hasGammaB))
       || ( (!beamAhasResGamma && hasGammaA) && beamBhasResGamma) )
    infoPtr->setHasUnresolvedBeams( true);

  // Pretend that the gamma-gamma system is the whole collision.
  infoPtr->setECM( mGmGm);

  // Reassign beam pointers in other classes.
  timesPtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  timesDecPtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  spacePtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  remnants.reassignBeamPtrs(  beamAPtr, beamBPtr, beamOffset);
  colourReconnection.reassignBeamPtrs(  beamAPtr, beamBPtr);

  // Set the MPI to point the gamma-gamma system.
  multiPtr = &multiGmGm;
  multiPtr->setBeamOffset(beamOffset);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Move back to CM-frame of colliding leptons. Add also the scatterd lepton
// if remnants are constructed.

void PartonLevel::leaveResolvedLeptonGamma( Event& process, Event& event,
  bool physical) {

  // Restore beam pointers to incoming leptons.
  if ( hasGammaA) beamAPtr = beamHadAPtr;
  if ( hasGammaB) beamBPtr = beamHadBPtr;

  // Restore cm energy.
  infoPtr->setECM( eCMsaveGamma);

  // Set the new kinematics only for physical events.
  if (physical) {

    // Find the momenta of incoming leptons and photons.
    Vec4 pLeptonA =  process[1].p();
    Vec4 pLeptonB =  process[2].p();

    // Momenta for scattered leptons to be set later according to final state.
    Vec4 pLepton1scat;
    Vec4 pLepton2scat;

    // Skip the boost of the extra intermediate photon in hard diffraction.
    int iSkipHardDiff = -1;

    // Reconstruct boost and rotation to event cm frame.
    RotBstMatrix MtoLeptonLepton;
    MtoLeptonLepton.toCMframe( pLeptonA, pLeptonB);

    // Boost event to cm frame.
    process.rotbst(MtoLeptonLepton);
    event.rotbst(MtoLeptonLepton);

    // Get the masses of beam particles.
    double m2BeamA  = pow2( beamAPtr->m());
    double m2BeamB  = pow2( beamBPtr->m());

    // Get the original collision energy and derive the lepton energies in CMS.
    double sCM      = infoPtr->s();
    double eCM2A    = 0.25 * pow2(sCM + m2BeamA - m2BeamB) / sCM;
    double eCM2B    = 0.25 * pow2(sCM - m2BeamA + m2BeamB) / sCM;

    // Find the current momenta of photons.
    Vec4 pGamma1Orig = process[3].p();
    Vec4 pGamma2Orig = process[4].p();
    Vec4 pGamma1     = pGamma1Orig;
    Vec4 pGamma2     = pGamma2Orig;
    double mGamma1   = sqrt(m2BeamA);
    double mGamma2   = sqrt(m2BeamB);

    // Separate treatment for 2 -> 1 processes with two direct gammas
    // to preserve sampled invariant mass. Corresponds to the usual
    // beam remnant handling.
    if ( (infoPtr->nFinal() == 1) && (gammaModeEvent == 4)
      && hasGammaA && hasGammaB ) {

      // Save the virtualities and transverse mometum of the photons.
      mGamma1     = -sqrt( beamAPtr->Q2Gamma());
      mGamma2     = -sqrt( beamBPtr->Q2Gamma());
      double kTxA = beamAPtr->gammaKTx();
      double kTyA = beamAPtr->gammaKTy();
      double kTxB = beamBPtr->gammaKTx();
      double kTyB = beamBPtr->gammaKTy();

      // Derive transverse mass and new sHat.
      double mT2A    = pow2( beamAPtr->gammaKT() ) - pow2(mGamma1);
      double mT2B    = pow2( beamBPtr->gammaKT() ) - pow2(mGamma2);
      double sHatBef = infoPtr->sHat();
      double sHatAft = infoPtr->sHat() + pow2(kTxA + kTxB) + pow2(kTyA + kTyB);
      double lambda  = pow2(sHatAft) + pow2(mT2A) + pow2(mT2B)
                     - 2. * ( sHatAft * (mT2A + mT2B) + mT2A * mT2B );
      double kz2New  = 0.25 * lambda / sHatAft;

      // Photon momenta in CM frame.
      Vec4 pGamma1New(kTxA, kTyA,  sqrt(kz2New), sqrt(kz2New + mT2A) );
      Vec4 pGamma2New(kTxB, kTyB, -sqrt(kz2New), sqrt(kz2New + mT2B) );

      // Boost and rotate the new photon momenta.
      RotBstMatrix MfromGmGmOrig;
      MfromGmGmOrig.toCMframe( pGamma1Orig, pGamma2Orig);
      MfromGmGmOrig.invert();
      pGamma1New.rotbst(MfromGmGmOrig);
      pGamma2New.rotbst(MfromGmGmOrig);

      // Set the new photon momenta.
      pGamma1 = pGamma1New;
      pGamma2 = pGamma2New;

      // Modify the mass and momenta in event record.
      event[3].p( pGamma1);
      event[3].m( mGamma1);
      event[4].p( pGamma2);
      event[4].m( mGamma2);

      // Light-cone momentum removed by the photons.
      double rescale = sqrt(sHatAft / sHatBef);
      double wPosBef = beamAPtr->xGamma() * infoPtr->eCM();
      double wNegBef = beamBPtr->xGamma() * infoPtr->eCM();

      // Transverse mass of the leptons and invariant mass left for remnants.
      double w2BeamA = pow2(beamAPtr->gammaKT()) + m2BeamA;
      double w2BeamB = pow2(beamBPtr->gammaKT()) + m2BeamB;
      double wPosRem = infoPtr->eCM() - rescale * wPosBef;
      double wNegRem = infoPtr->eCM() - rescale * wNegBef;
      double w2Rem   = wPosRem * wNegRem;
      double xLepA   = 1. - beamAPtr->xGamma();
      double xLepB   = 1. - beamBPtr->xGamma();

      // Rescaling factors for leptons.
      double lambdaRoot = sqrtpos( pow2(w2Rem - w2BeamA - w2BeamB)
                        - 4. * w2BeamA * w2BeamB );
      double rescaleA   = (w2Rem + w2BeamA - w2BeamB + lambdaRoot)
                        / (2. * w2Rem * xLepA);
      double rescaleB   = (w2Rem + w2BeamB - w2BeamA + lambdaRoot)
                        / (2. * w2Rem * xLepB);

      // Momenta of the scattered leptons.
      double pPosA  = rescaleA * xLepA * wPosRem;
      double pNegA  = w2BeamA / pPosA;
      double eLepA  = 0.5 * (pPosA + pNegA);
      double pzLepA = 0.5 * (pPosA - pNegA);
      double pNegB  = rescaleB * xLepB * wNegRem;
      double pPosB  = w2BeamB / pNegB;
      double eLepB  = 0.5 * (pPosB + pNegB);
      double pzLepB = 0.5 * (pPosB - pNegB);
      Vec4 pLeptonANew( -kTxA, -kTyA, pzLepA, eLepA);
      Vec4 pLeptonBNew( -kTxB, -kTyB, pzLepB, eLepB);

      // Save the momenta of scattered leptons.
      pLepton1scat = pLeptonANew;
      pLepton2scat = pLeptonBNew;

      // Save the scattering angles of the leptons.
      double theta1 = pLepton1scat.theta();
      double theta2 = M_PI - pLepton2scat.theta();
      infoPtr->setTheta1(theta1);
      infoPtr->setTheta2(theta2);
      infoPtr->setECMsub(sqrt(sHatBef));
      infoPtr->setsHatNew(sHatBef);

    // Otherwise derive the kinematics according to sampled virtualities.
    } else {

      // Calculate the new momentum for virtual photon if present for side A.
      if ( hasGammaA ) {

        // Get the x_gamma and virtuality and derive mass.
        double xGamma1  = beamAPtr->xGamma();
        double Q2gamma1 = beamAPtr->Q2Gamma();
        mGamma1         = -sqrt( Q2gamma1);
        beamGamAPtr->newM( mGamma1);

        // Derive the kinematics with virtuality and kT.
        double eGamma1  = xGamma1 * sqrt( eCM2A);
        double kz1      = (eCM2A * xGamma1 + 0.5 * Q2gamma1)
                        / sqrt(eCM2A - m2BeamA);

        // Set the new momemtum and mass.
        pGamma1 = Vec4( beamAPtr->gammaKTx(), beamAPtr->gammaKTy(),
          kz1, eGamma1);
        event[3].p( pGamma1);
        event[3].m( mGamma1);

        // For hard diffraction set the new kinematics also for the copy.
        if (infoPtr->isHardDiffractiveA()) {
          event[7].p( pGamma1);
          event[7].m( mGamma1);
          iSkipHardDiff = 7;
        }
      }

      // Calculate the new momentum for virtual photon if present for side B.
      if ( hasGammaB ) {

        // Get the x_gamma and virtuality and derive mass.
        double xGamma2  = beamBPtr->xGamma();
        double Q2gamma2 = beamBPtr->Q2Gamma();
        mGamma2         = -sqrt( Q2gamma2);
        beamGamBPtr->newM( mGamma2);

        // Derive the kinematics with virtuality and kT.
        double eGamma2  = xGamma2 * sqrt( eCM2B);
        double kz2      = (eCM2B * xGamma2 + 0.5 * Q2gamma2)
                        / sqrt(eCM2B - m2BeamB);

        // Save the 4-momentum of photons with sampled kT.
        pGamma2 = Vec4( beamBPtr->gammaKTx(), beamBPtr->gammaKTy(),
          -kz2, eGamma2);
        event[4].p( pGamma2);
        event[4].m( mGamma2);

        // For hard diffraction set the new kinematics also for the copy.
        if (infoPtr->isHardDiffractiveB()) {
          event[8].p( pGamma2);
          event[8].m( mGamma2);
          iSkipHardDiff = 8;
        }
      }

      // Find momenta for scattered lepton.
      Vec4 pLepton1 = process[1].p();
      Vec4 pLepton2 = process[2].p();
      pLepton1scat  = pLepton1 - pGamma1;
      pLepton2scat  = pLepton2 - pGamma2;
    }

    // For photon-hadron cases use the 4-momentum of the original beam particle
    // since virtual photon kinematics are derived in this frame.
    // Effect negligible when mass of the beam providing photon is small.
    if      ( hasGammaA && !hasGammaB) pGamma2 = event[2].p();
    else if (!hasGammaA &&  hasGammaB) pGamma1 = event[1].p();

    // Find the boost from rest frame of collinear photons to rest frame of
    // photons with kT.
    RotBstMatrix MfromGmGm;
    MfromGmGm.toCMframe( pGamma1, pGamma2);
    MfromGmGm.fromCMframe( pGamma1Orig, pGamma2Orig);
    MfromGmGm.invert();

    // Copy the momentum and mass of the unresolved photon for direct-resolved
    // processes to have correct virtualities in the event record.
    int iSkipGamma = -1;
    if ( gammaModeEvent == 3 ) {
      iSkipGamma = 5;
      event[iSkipGamma].m( mGamma1);
      event[iSkipGamma].p( pGamma1);
    } else if ( gammaModeEvent == 2 ) {
      iSkipGamma = 6;
      event[iSkipGamma].m( mGamma2);
      event[iSkipGamma].p( pGamma2);
    }

    // Boost scattered system to frame where photon beam has non-zero kT.
    // Do not boost photons which has already correct four momenta.
    // For elastic events skip both photons.
    int iProcessBegin = 5;
    int iSkipGammaEl = -1;
    if (isElastic) {
      iProcessBegin = 3;
      if (hasGammaA) iSkipGamma   = 3;
      if (hasGammaB) iSkipGammaEl = 4;
    }
    for (int i = iProcessBegin; i < event.size(); ++i) {
      if ( (i != iSkipGamma) && (i != iSkipHardDiff) && (i != iSkipGammaEl) )
        event[i].rotbst( MfromGmGm);
    }

    // For photon-hadron cases set the beam particle copy originally derived
    // using zero virtuality for photon to match the derived kinematics.
    if      ( hasGammaA && !hasGammaB) event[4].p(pGamma2);
    else if (!hasGammaA &&  hasGammaB) event[3].p(pGamma1);

    // Add the scattered leptons if remnants are constructed.
    if (doRemnants) {

      // Add scattered leptons and fix the daughter codes.
      if ( hasGammaA) {
        int iPosLepton1 = event.append( beamAPtr->id(), 63, 1, 0, 0, 0, 0, 0,
          pLepton1scat, beamAPtr->m());
        event[1].daughter2( event[1].daughter1());
        event[1].daughter1( iPosLepton1);
      }
      if ( hasGammaB) {
        int iPosLepton2 = event.append( beamBPtr->id(), 63, 2, 0, 0, 0, 0, 0,
          pLepton2scat, beamBPtr->m());
        event[2].daughter2( event[2].daughter1());
        event[2].daughter1( iPosLepton2);
      }
    }

  // Reset all pointers also for non-physical events.
  }

  // Done for direct-direct and elastic processes.
  if ( (gammaModeEvent == 4) || isElastic ) return;

  // Reassign beam pointers in other classes.
  timesPtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  timesDecPtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  spacePtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  remnants.reassignBeamPtrs(  beamAPtr, beamBPtr, 0);
  colourReconnection.reassignBeamPtrs(  beamAPtr, beamBPtr);

  // Set the MPI pointer back to the original collisions.
  multiPtr = &multiMB;
  multiPtr->setBeamOffset(0);

}

//--------------------------------------------------------------------------

// Remove the copies of the beam photon from the event record.

void PartonLevel::cleanEventFromGamma( Event& event) {

  // Offset to normal beam position when photons emitted from a lepton beam.
  int beamOffset = 2;
  // If hard diffraction more offset for sub-system beams.
  if (infoPtr->isHardDiffractiveA() || infoPtr->isHardDiffractiveB())
    beamOffset += 4;
  int iPosBeam1  = 1 + beamOffset;
  int iPosBeam2  = 2 + beamOffset;

  // Go through the event record from the end and find the copies.
  int iPosGamma1 = 0;
  int iPosGamma2 = 0;
  for (int i = event.size() - 1; i > 0; --i) {
    if ( (event[i].id() == 22 && event[i].mother1() == iPosBeam1)
      && beamAhasResGamma ) iPosGamma1 = i;
    if ( (event[i].id() == 22 && event[i].mother1() == iPosBeam2)
      && beamBhasResGamma ) iPosGamma2 = i;
  }

  // Check how many unresolved photons are present in the event.
  int nGamma = 0;
  if (iPosGamma1 > 0) ++nGamma;
  if (iPosGamma2 > 0) ++nGamma;

  // Exit if no copies found.
  if ( nGamma == 0 ) return;

  // Loop over two beams if required.
  for (int i = 0; i < nGamma; ++i) {

    // Set the positions to match the beam.
    int iPosGamma = (iPosGamma1 > 0 && i == 0) ? iPosGamma1 : iPosGamma2;
    int iPosBeam  = (iPosGamma1 > 0 && i == 0) ? iPosBeam1  : iPosBeam2;

    // Go through the history of the beam photon.
    while ( iPosGamma > iPosBeam ) {
      int iDaughter1 = event[iPosGamma].daughter1();
      int iDaughter2 = event[iPosGamma].daughter2();
      int iMother1   = event[iPosGamma].mother1();
      int iMother2   = event[iPosGamma].mother2();

      // If equal daughters photon just a carbon copy.
      if ( iDaughter1 == iDaughter2 ) {
        event[iDaughter1].mothers( iMother1, iMother2 );
        event.remove( iPosGamma, iPosGamma);
        iPosGamma = iDaughter1;

      // If non-equal daughters the photon from ISR branching.
      } else {
        event[iMother1].daughters( iDaughter1, iDaughter2 );
        event[iDaughter1].mother1( iMother1 );
        event[iDaughter2].mother1( iMother1 );
        event.remove( iPosGamma, iPosGamma);
        iPosGamma = iMother1;
      }

      // If both beams unresolved fix the position of the latter one.
      if ( (i == 0 && nGamma > 1) && iPosGamma2 > iPosGamma ) --iPosGamma2;

    }
  }
}

//--------------------------------------------------------------------------

// Handle showers in successive resonance decays.

bool PartonLevel::resonanceShowers( Event& process, Event& event,
  bool skipForR) {

  // Prepare to start over from beginning for R-hadron decays.
  if (allowRH) {
    if (skipForR) {
      nHardDoneRHad = nHardDone;
      inRHadDecay.resize(0);
      for (int i = 0; i < process.size(); ++i)
        inRHadDecay.push_back( false);
    } else nHardDone = nHardDoneRHad;
  }

  // Isolate next system to be processed, if anything remains.
  int nRes    = 0;
  int nFSRres = 0;
  // Number of desired branchings, negative value means no restriction.
  int nBranchMax = (doTrial) ? nTrialEmissions : -1;
  // Vector to tell which junctions have already been copied
  vector<int> iJunCopied;

  while (nHardDone < process.size()) {
    ++nRes;
    int iBegin = nHardDone;

    // In first call (skipForR = true) skip over daughters
    // of resonances that should form R-hadrons
    if (allowRH) {
      if (skipForR) {
        bool comesFromR = false;
        int iTraceUp = process[iBegin].mother1();
        do {
          if ( rHadronsPtr->givesRHadron(process[iTraceUp].id()) )
            comesFromR = true;
          iTraceUp = process[iTraceUp].mother1();
        } while (iTraceUp > 0 && !comesFromR);
        if (comesFromR) {
          inRHadDecay[iBegin] = true;
          ++nHardDone;
          continue;
        }

      // In optional second call (skipForR = false) process decay chains
      // inside R-hadrons.
      } else if (!inRHadDecay[iBegin]) {
        ++nHardDone;
        continue;
      }
    }

    // Mother in hard process and in complete event (after shower).
    int iHardMother      = process[iBegin].mother1();
    Particle& hardMother = process[iHardMother];
    int iBefMother       = iPosBefShow[iHardMother];
    int iAftMother       = event[iBefMother].iBotCopyId();
    // Possibly trace across intermediate R-hadron state.
    if (allowRH) {
      int iTraceRHadron    = rHadronsPtr->trace( iAftMother);
      if (iTraceRHadron > 0) iAftMother = iTraceRHadron;
    }
    Particle& aftMother  = event[iAftMother];

    // From now on mother counts as decayed.
    aftMother.statusNeg();

    // Mother can have been moved by showering (in any of previous steps),
    // so prepare to update colour and momentum information for system.
    int colBef  = hardMother.col();
    int acolBef = hardMother.acol();
    int colAft  = aftMother.col();
    int acolAft = aftMother.acol();
    RotBstMatrix M;
    M.bst( hardMother.p(), aftMother.p());

    // New colour reconnection can not handle late resonance decay
    // of coloured particles so abort event.
    if ( (colBef != 0 || acolBef != 0) && doReconnect && reconnectMode == 1
      && forceResonanceCR && !earlyResDec) {
      infoPtr->errorMsg("Abort from PartonLevel::resonanceShower: "
        "new CR can't handle separate CR for coloured resonance decays");
      infoPtr->setAbortPartonLevel(true);
      return false;
    }

    // Extract next partons from hard event into normal event record.
    vector<bool> doCopyJun;
    for (int i = iBegin; i < process.size(); ++i) {
      if (process[i].mother1() != iHardMother) break;
      int iNow = event.append( process[i] );
      iPosBefShow[i] = iNow;
      Particle& now = event.back();

      // Currently outgoing ones should not count as decayed.
      if (now.status() == -22) {
        now.statusPos();
        now.daughters(0, 0);
      }

      // Update daughter and mother information.
      if (i == iBegin) event[iAftMother].daughter1( iNow);
      else             event[iAftMother].daughter2( iNow);
      now.mother1(iAftMother);

      // Check if this parton came from a BNV (junction) decay in hard event.
      for (int iJun = 0; iJun < process.sizeJunction(); ++iJun) {
        if (iJun >= int(doCopyJun.size())) doCopyJun.push_back(false);
        // Skip if we already decided we're going to copy this junction
        if (doCopyJun[iJun]) continue;
        // Only consider junctions that can appear in decays.
        int kindJunction = process.kindJunction(iJun);
        if (kindJunction <= 2) {
          // Junction Kinds 1 and 2: all legs in final state
          int nMatch = 0;
          // Loop over resonance-decay daughters
          int iDau1 = hardMother.daughter1();
          int iDau2 = hardMother.daughter2();
          // Must have at least 3 decay products
          if (iDau1 == 0 || iDau2 - iDau1 < 2) continue;
          for (int iDau=iDau1; iDau<=iDau2; ++iDau) {
            int colDau = (kindJunction == 1 ? process[iDau].col()
              : process[iDau].acol());
            for (int iLeg = 0; iLeg <= 2; ++iLeg)
              if ( process.colJunction(iJun,iLeg) == colDau ) nMatch += 1;
          }
          // If three legs match
          if (nMatch == 3) doCopyJun[iJun] = true;
        } else if (kindJunction <= 4) {
          // Junction Kinds 3 and 4: copy if initial-state leg matches
          // this resonance.
          int col = (kindJunction == 3 ? hardMother.acol() : hardMother.col());
          if ( process.colJunction(iJun,0) == col ) doCopyJun[iJun] = true;
        }
        // Extra safety: Check if this junction has already been copied
        // (e.g., in setupHardSys). If so, do not copy again.
        for (int kJun = 0; kJun < event.sizeJunction(); ++kJun) {
          int nMatch = 0;
          for (int iLeg = 0; iLeg <= 2; ++iLeg)
            if (event.colJunction(kJun,iLeg) == process.colJunction(iJun,iLeg))
              ++nMatch;
          if (nMatch == 3) doCopyJun[iJun] = false;
        }
      }

      // Update colour and momentum information.
      if (now.col() == colBef) now.col( colAft);
      if (now.acol() == acolBef) now.acol( acolAft);
      // Sextet mothers have additional (negative) tag
      if (now.col() == -acolBef) now.col( -acolAft);
      if (now.acol() == -colBef) now.acol( -colAft);
      now.rotbst( M);

      // Update vertex information.
      if (now.hasVertex()) now.vProd( event[iAftMother].vDec() );

      // Complete task of copying next subsystem into event record.
      ++nHardDone;
    }
    int iEnd = nHardDone - 1;

    // Copy down junctions from hard event into normal event record.
    for (int iJun = 0; iJun < int(doCopyJun.size()); ++iJun) {
      // Check if this junction was already copied
      for (int jJun = 0; jJun < int(iJunCopied.size()); ++jJun)
        if (iJunCopied[jJun] == iJun) doCopyJun[iJun] = false;
      // Skip if not doing anything
      if (!doCopyJun[iJun]) continue;
      // Check for changed colors and update as necessary.
      Junction junCopy = process.getJunction(iJun);
      for (int iLeg = 0; iLeg <= 2; ++iLeg) {
        int colLeg = junCopy.col(iLeg);
        if (colLeg == colBef) junCopy.col(iLeg, colAft);
        if (colLeg == acolBef) junCopy.col(iLeg, acolAft);
      }
      event.appendJunction(junCopy);
      // Mark junction as copied (to avoid later recopying)
      iJunCopied.push_back(iJun);
    }

    // Reset pT of last branching
    pTLastBranch = 0.0;

    // Add new resonance-decay system (with empty beam slots).
    int iSys = partonSystemsPtr->addSys();
    partonSystemsPtr->setInRes( iSys,iAftMother);
    partonSystemsPtr->setSHat( iSys, pow2(hardMother.m()) );
    partonSystemsPtr->setPTHat( iSys, 0.5 * hardMother.m() );

    // Loop over allowed range to find all final-state particles.
    for (int i = iPosBefShow[iBegin]; i <= iPosBefShow[iEnd]; ++i)
    if (event[i].isFinal()) partonSystemsPtr->addOut( iSys, i);

    // Do parton showers inside subsystem: maximum scale by mother mass.
    if (doFSRinResonances) {
      double pTmax = 0.5 * hardMother.m();
      if (canSetScale) pTmax
        = userHooksPtr->scaleResonance( iAftMother, event);
      nFSRhard     = 0;

      // Set correct scale for trial showers.
      if (doTrial) pTmax = process.scale();

      // Set correct scale for showers off multi-parton events when
      // merging e+e- -> V -> jets.
      int iMother1 = hardMother.mother1();
      int iMother2 = hardMother.mother2();
      if ( canRemoveEvent
        && event[iMother1].colType() == 0 && event[iMother2].colType() == 0)
        pTmax = process.scale();

      // Let prepare routine do the setup.
      timesDecPtr->prepare( iSys, event);

       // Number of actual branchings
      int nBranch = 0;

      // Set up initial veto scale.
      doVeto        = false;
      double pTveto = pTvetoPT;

      // Begin evolution down in pT from hard pT scale.
      do {
        typeVetoStep = 0;
        double pTtimes = timesDecPtr->pTnext( event, pTmax, 0.);

        // Allow a user veto. Only do it once, so remember to change pTveto.
        if (pTveto > 0. && pTveto > pTtimes) {
          pTveto = -1.;
          doVeto = userHooksPtr->doVetoPT( 5, event);
          // Abort event if vetoed.
          if (doVeto) return false;
        }

        // Do a final-state emission (if allowed).
        if (pTtimes > 0.) {
          if (timesDecPtr->branch( event)) {
            ++nFSRres;
            ++nFSRhard;
            if (canVetoStep && nFSRhard <= nVetoStep) typeVetoStep = 5;

            nBranch++;
            pTLastBranch = pTtimes;
            typeLastBranch = 5;

          }
          pTmax = pTtimes;
        }

        // If no pT scales above zero then nothing to be done.
        else pTmax = 0.;

        // Optionally check for a veto after the first few emissions.
        if (typeVetoStep > 0) {
          doVeto = userHooksPtr->doVetoStep( typeVetoStep, 0, nFSRhard,
            event);
          // Abort event if vetoed.
          if (doVeto) return false;
        }

        // Handle potential merging veto.
        if ( canRemoveEvent && nFSRhard == 1 ) {
          // Simply check, and possibly reset weights.
          mergingHooksPtr->doVetoStep( process, event, true );
        }

      // Keep on evolving until nothing is left to be done.
      } while (pTmax > 0.  && (nBranchMax <= 0 || nBranch < nBranchMax) );

    }

  // No more systems to be processed. Set total number of emissions.
  }
  if (skipForR) nFSRinRes = nFSRres;
  return true;

}

//--------------------------------------------------------------------------

// Perform decays and showers of W and Z emitted in shower.

bool PartonLevel::wzDecayShowers( Event& event) {

  // Identify W/Z produced by a parton shower.
  for (int iWZ = 0; iWZ < event.size(); ++iWZ)
  if (event[iWZ].isFinal()
    && (event[iWZ].id() == 23 || event[iWZ].idAbs() == 24) ) {

    // Do nothing if particle should not be decayed.
    if ( event[iWZ].canDecay() && event[iWZ].mayDecay()
      && event[iWZ].isResonance() ) ;
    else continue;

    int iWZtop = event[iWZ].iTopCopy();
    int typeWZ = 0;
    if (event[iWZtop].statusAbs() == 56) typeWZ = 1;
    if (event[iWZtop].statusAbs() == 47) typeWZ = 2;
    int sizeSave = event.size();

    // Map id_Z = 23 -> 93 and id_W = 24 -> 94, for separate decay settings.
    // Let W/Z resonance decay. Restore correct identity and status codes.
    if (typeWZ > 0) {
      int idSave     = event[iWZ].id();
      event[iWZ].id( (idSave > 0) ? idSave + 70 : idSave - 70);
      int statusSave = event[iWZ].status();
      resonanceDecays.next( event, iWZ);
      event[iWZ].id( idSave);
      if (event.size() - sizeSave != 2) return true;
      event[iWZ].status( -statusSave);
    }

    // FSR decays.
    if (typeWZ == 1) {

      // Identify fermion after W/Z emission.
      vector<int>  iSisters = event[iWZtop].sisterList();
      if (iSisters.size() != 1) {
        infoPtr->errorMsg("Error in PartonLevel::wzDecayShowers: "
          "Not able to find a single sister particle");
        return false;
      }
      int iEmitter = iSisters[0];

      // Boosts to study decay in W/Z rest frame.
      RotBstMatrix MtoNew, MtoRest, MtoCM;
      MtoNew.bst( event[iWZtop].p(), event[iWZ].p());
      MtoRest.bstback( event[iWZ].p());
      MtoCM.bst( event[iWZ].p());

      // Emitter and recoiler in W/Z rest frame.
      Vec4 pEmitter = event[iEmitter].p();
      pEmitter.rotbst( MtoNew);
      pEmitter.rotbst( MtoRest);
      if (event[iWZtop + 1].statusAbs() != 52) {
        infoPtr->errorMsg("Error in PartonLevel::wzDecayShowers: "
          "Found wrong recoiler");
        return false;
      }
      Vec4 pRecoiler = event[iWZtop + 1].p();
      pRecoiler.rotbst( MtoNew);
      pRecoiler.rotbst( MtoRest);
      Vec4 pWZRest = event[iWZ].p();
      pWZRest.rotbst( MtoRest);

      // Always choose p4 as the particle and p5 as the anti-particle.
      Vec4 p4 = pEmitter;
      Vec4 p5 = pRecoiler;
      if (event[iEmitter].id() < 0) swap( p4, p5);

      // Decay daughters in W/Z rest frame.
      // Always choose pDec1 as the particle and p2Dec as the anti-particle.
      Vec4 pDec1 = event[sizeSave].p();
      Vec4 pDec2 = event[sizeSave + 1].p();
      if (event[sizeSave].id() < 0) swap( pDec1, pDec2);
      pDec1.rotbst( MtoRest);
      pDec2.rotbst( MtoRest);

      // Couplings.
      double li2, ri2, lf2, rf2;
      // Z couplings: make use of initial fermion polarization if set.
      if (event[iWZ].id() == 23) {
        li2 = pow2(couplingsPtr->lf( event[iEmitter].idAbs() ));
        ri2 = pow2(couplingsPtr->rf( event[iEmitter].idAbs() ));
        lf2 = pow2(couplingsPtr->lf( event[sizeSave].idAbs() ));
        rf2 = pow2(couplingsPtr->rf( event[sizeSave].idAbs() ));
        if ( abs( event[iEmitter].pol() + 1.) < 0.1) ri2 = 0.;
        if ( abs( event[iEmitter].pol() - 1.) < 0.1) li2 = 0.;
      // W couplings.
      } else {
        li2 = 1.;
        ri2 = 0.;
        lf2 = 1.;
        rf2 = 0.;
      }

      // Different needed kinematic variables.
      double sWZER = (p4 + pWZRest + p5).m2Calc();
      double x1    = 2. * p4 * (p4 + pWZRest + p5) / sWZER;
      double x2    = 2. * p5 * (p4 + pWZRest + p5) / sWZER;
      double x1s   = x1 * x1;
      double x2s   = x2 * x2;
      double m2Sel = pWZRest.m2Calc();
      double rWZER = m2Sel / sWZER;

      // Calculate constants needed in correction.
      double con[9];
      con[0] = 2. * m2Sel * (1.-x1) * ((x2s+1.-x1-x2) - rWZER * (1.-x2))
             * (li2 * lf2 + ri2 * rf2);
      con[1] = 2. * m2Sel * (1.-x2) * ((x1s+1.-x1-x2) - rWZER * (1.-x1))
             * (li2 * rf2 + ri2 * lf2);
      con[2] = 2. * m2Sel * (1.-x1) * ((x2s+1.-x1-x2) - rWZER * (1.-x2))
             * (li2 * rf2 + ri2 * lf2);
      con[3] = 2. * m2Sel * (1.-x2) * ((x1s+1.-x1-x2) - rWZER * (1.-x1))
             * (li2 * lf2 + ri2 * rf2);
      con[4] = m2Sel * sWZER * (1.-x1) * (1.-x2) * ((x1+x2-1.) + rWZER)
             * (li2 + ri2) * (lf2 + rf2);
      con[5] = -4. * (1.-x1) * (1.-x2) * (li2 + ri2) * (lf2 + rf2);
      con[6] = -4. * (1.-x1) * (1.-x2) * (li2 + ri2) * (lf2 + rf2);
      con[7] = 4. * (1.-x1) * ((1.-x1) - rWZER * (1.-x2))
             * (li2 + ri2) * (lf2 + rf2);
      con[8] = 4. * (1.-x2) * ((1.-x2) - rWZER * (1.-x1))
             * (li2 + ri2) * (lf2 + rf2);

      // Find maximum value: try pDec1 and pDec2 = -pDec1 along +-x, +-y, +-z.
      double wtMax  = 0.;
      double pAbs12 = pDec1.pAbs();
      for (int j = 0; j < 6; ++j) {
        Vec4 pDec1Test( 0., 0., 0., pDec1.e());
        Vec4 pDec2Test( 0., 0., 0., pDec2.e());
        if      (j == 0) { pDec1Test.px(  pAbs12);  pDec1Test.px( -pAbs12);}
        else if (j == 1) { pDec1Test.px( -pAbs12);  pDec1Test.px(  pAbs12);}
        else if (j == 2) { pDec1Test.py(  pAbs12);  pDec1Test.py( -pAbs12);}
        else if (j == 3) { pDec1Test.py( -pAbs12);  pDec1Test.py(  pAbs12);}
        else if (j == 4) { pDec1Test.pz(  pAbs12);  pDec1Test.pz( -pAbs12);}
        else if (j == 5) { pDec1Test.pz( -pAbs12);  pDec1Test.pz(  pAbs12);}

        // Evaluate matrix element and compare with current maximum.
        double p2p4Test = p4 * pDec1Test;
        double p3p4Test = p4 * pDec2Test;
        double p2p5Test = p5 * pDec1Test;
        double p3p5Test = p5 * pDec2Test;
        double testValues[9] = { p2p4Test, p2p5Test, p3p4Test, p3p5Test, 1.,
          p2p5Test * p3p4Test, p2p4Test * p3p5Test, p2p4Test * p3p4Test,
          p2p5Test * p3p5Test};
        double wtTest = 0.;
        for (int i = 0; i < 9; ++i) wtTest += con[i] * testValues[i];
        if (wtTest > wtMax) wtMax = wtTest;
      }

      // Multiply by four to ensure maximum is an overestimate.
      wtMax *= 4.;

      // Iterate with new angles until weighting succeeds.
      int nRot = -1;
      double wt = 0.;
      do {
        ++nRot;
        if (nRot > 0) {
          RotBstMatrix MrndmRot;
          MrndmRot.rot( acos(2. * rndmPtr->flat() - 1.),
            2. * M_PI * rndmPtr->flat());
          pDec1.rotbst(MrndmRot);
          pDec2.rotbst(MrndmRot);
        }

        // p2 is decay product, p3 is anti decay product,
        // p4 is dipole particle, p5 is dipole anti particle.
        // So far assumed that we always have qQ-dipole.
        double p2p4 = p4 * pDec1;
        double p3p4 = p4 * pDec2;
        double p2p5 = p5 * pDec1;
        double p3p5 = p5 * pDec2;

        // Calculate weight and compare with maximum weight.
        double wtValues[9] = { p2p4, p2p5, p3p4, p3p5, 1., p2p5 * p3p4,
          p2p4 * p3p5, p2p4 * p3p4, p2p5 * p3p5};
        wt =  0.;
        for (int i = 0; i < 9; ++i) wt += con[i] * wtValues[i];
        if (wt > wtMax || wt < 0.) {
          infoPtr->errorMsg("Error in PartonLevel::wzDecayShowers: "
          "wt bigger than wtMax or less than zero");
          return false;
        }
      } while (wt < wtMax * rndmPtr->flat());

      // If momenta rotated then store new ones.
      if (nRot > 0) {
        pDec1.rotbst( MtoCM);
        pDec2.rotbst( MtoCM);
        if (event[sizeSave].id() > 0) {
          event[sizeSave].p( pDec1);
          event[sizeSave + 1].p( pDec2);
        }
        else {
          event[sizeSave].p( pDec2);
          event[sizeSave + 1].p( pDec1);
        }
      }
    }

    // ISR decays.
    if (typeWZ == 2) {

      // Identify mother of W/Z emission.
      double iMother = event[iWZtop].mother1();

      // Boosts to study decay in W/Z rest frame.
      RotBstMatrix MtoNew, MtoRest, MtoCM;
      MtoNew.bst( event[iWZtop].p(), event[iWZ].p());
      MtoRest.bstback( event[iWZ].p());
      MtoCM.bst( event[iWZ].p());

      // Find recoiler.
      int iRecoiler;
      if (event[iMother + 1].statusAbs() == 42) iRecoiler = iMother + 1;
      else if (event[iMother - 1].statusAbs() == 42) iRecoiler = iMother - 1;
      else {
        infoPtr->errorMsg("Error in PartonLevel::wzDecayShowers: "
          "Found wrong recoiler");
        return false;
      }

      // Emitter and recoiler in W/Z rest frame.
      Vec4 pMother = event[iMother].p();
      pMother.rotbst( MtoNew);
      pMother.rotbst( MtoRest);
      Vec4 pRecoiler = event[iRecoiler].p();
      pRecoiler.rotbst( MtoNew);
      pRecoiler.rotbst( MtoRest);
      Vec4 pWZRest = event[iWZ].p();
      pWZRest.rotbst( MtoRest);

      // Always choose p1 as the particle and p2 as the anti-particle.
      // If no anti-particles just continue.
      Vec4 p1 = pMother;
      Vec4 p2 = pRecoiler;
      if (event[iMother].id() < 0) swap( p1, p2);

      // Decay daughters in W/Z rest frame.
      // Always choose pDec1 as the particle and p2Dec as the anti-particle.
      Vec4 pDec1 = event[sizeSave].p();
      Vec4 pDec2 = event[sizeSave + 1].p();
      if (event[sizeSave].id() < 0) swap( pDec1, pDec2);
      pDec1.rotbst( MtoRest);
      pDec2.rotbst( MtoRest);

      // Couplings.
      double li2, ri2, lf2, rf2;
      // Z couplings: make use of initial fermion polarization if set.
      if (event[iWZ].id() == 23) {
        li2 = pow2(couplingsPtr->lf( event[iMother].idAbs() ));
        ri2 = pow2(couplingsPtr->rf( event[iMother].idAbs() ));
        lf2 = pow2(couplingsPtr->lf( event[sizeSave].idAbs() ));
        rf2 = pow2(couplingsPtr->rf( event[sizeSave].idAbs() ));
        if ( abs( event[iMother].pol() + 1.) < 0.1) ri2 = 0.;
        if ( abs( event[iMother].pol() - 1.) < 0.1) li2 = 0.;
      // W couplings.
      } else {
        li2 = 1.;
        ri2 = 0.;
        lf2 = 1.;
        rf2 = 0.;
      }

      // Different needed kinematic variables.
      double s = (p1 + p2).m2Calc();
      double t = (p1-pWZRest).m2Calc();
      double u = (p2-pWZRest).m2Calc();
      double m2Sel = pWZRest.m2Calc();
      double m2Final = t + u + s - m2Sel;

      // Calculate constants needed in correction.
      double con[9];
      con[0] = 2. * m2Sel * (u * (1. - m2Final / t) + s)
        * (li2 * rf2 + ri2 * lf2);
      con[1] = 2. * m2Sel * (u * (1. - m2Final / t) + s)
        * (li2 * lf2 + ri2 * rf2);
      con[2] = 2. * m2Sel * (t * (1. - m2Final / u) + s)
        * (li2 * lf2 + ri2 * rf2);
      con[3] = 2. * m2Sel * (t * (1. - m2Final / u) + s)
        * (li2 * rf2 + ri2 * lf2);
      con[4] = m2Sel * m2Final * s * (li2 + ri2) * (lf2 + rf2);
      con[5] = -4. * m2Final * (li2 + ri2) * (lf2 + rf2);
      con[6] = -4. * m2Final * (li2 + ri2) * (lf2 + rf2);
      con[7] = 4. * (m2Final * u / t - m2Sel) * (li2 + ri2) * (lf2 + rf2);
      con[8] = 4. * (m2Final * t / u - m2Sel) * (li2 + ri2) * (lf2 + rf2);

      // Find maximum value: try pDec1 and pDec2 = -pDec1 along +-x, +-y, +-z.
      double wtMax  = 0.;
      double pAbs12 = pDec1.pAbs();
      for (int j = 0; j < 6; ++j) {
        Vec4 pDec1Test( 0., 0., 0., pDec1.e());
        Vec4 pDec2Test( 0., 0., 0., pDec2.e());
        if      (j == 0) { pDec1Test.px(  pAbs12);  pDec1Test.px( -pAbs12);}
        else if (j == 1) { pDec1Test.px( -pAbs12);  pDec1Test.px(  pAbs12);}
        else if (j == 2) { pDec1Test.py(  pAbs12);  pDec1Test.py( -pAbs12);}
        else if (j == 3) { pDec1Test.py( -pAbs12);  pDec1Test.py(  pAbs12);}
        else if (j == 4) { pDec1Test.pz(  pAbs12);  pDec1Test.pz( -pAbs12);}
        else if (j == 5) { pDec1Test.pz( -pAbs12);  pDec1Test.pz(  pAbs12);}

        // Evaluate matrix element and compare with current maximum.
        double p1p4Test = p1 * pDec1Test;
        double p1p5Test = p1 * pDec2Test;
        double p2p4Test = p2 * pDec1Test;
        double p2p5Test = p2 * pDec2Test;
        double testValues[9] = { p1p4Test, p1p5Test, p2p4Test, p2p5Test, 1.,
          p1p4Test * p2p5Test, p1p5Test * p2p4Test, p1p4Test * p1p5Test,
          p2p4Test * p2p5Test};
        double wtTest = 0.;
        for (int i = 0; i < 9; ++i) wtTest += con[i] * testValues[i];
        if (wtTest > wtMax) wtMax = wtTest;
      }

      // Multiply by four to ensure maximum is an overestimate.
      wtMax *= 4.;

      // Iterate with new angles until weighting succeeds.
      int nRot = -1;
      double wt = 0.;
      do {
        ++nRot;
        if (nRot > 0) {
          RotBstMatrix MrndmRot;
          MrndmRot.rot( acos(2. * rndmPtr->flat() - 1.),
            2. * M_PI * rndmPtr->flat());
          pDec1.rotbst(MrndmRot);
          pDec2.rotbst(MrndmRot);
        }

        // p2 is decay product, p3 is anti decay product,
        // p4 is dipole particle, p5 is dipole anti particle.
        // So far assumed that we always have qQ-dipole.
        double p1p4 = p1 * pDec1;
        double p1p5 = p1 * pDec2;
        double p2p4 = p2 * pDec1;
        double p2p5 = p2 * pDec2;

        // Calculate weight and compare with maximum weight.
        double wtValues[9] = { p1p4, p1p5, p2p4, p2p5, 1., p1p4 * p2p5,
          p1p5 * p2p4, p1p4 * p1p5, p2p4 * p2p5};
        wt =  0.;
        for (int i = 0; i < 9; ++i) wt += con[i] * wtValues[i];
        if (wt > wtMax || wt < 0.) {
          infoPtr->errorMsg("Error in PartonLevel::wzDecayShowers: "
            "wt bigger than wtMax or less than zero");
          return false;
        }
      } while (wt < wtMax * rndmPtr->flat());

      // If momenta rotated then store new ones.
      if (nRot > 0) {
        pDec1.rotbst( MtoCM);
        pDec2.rotbst( MtoCM);
        if (event[sizeSave].id() > 0) {
          event[sizeSave].p( pDec1);
          event[sizeSave + 1].p( pDec2);
        }
        else {
          event[sizeSave].p( pDec2);
          event[sizeSave + 1].p( pDec1);
        }
      }
    }

    // Add new resonance-decay system (with empty beam slots).
    if (typeWZ > 0) {
      // Maximum shower scale set by mother mass.
      double pTmax = 0.5 * event[iWZ].m();
      int iSys = partonSystemsPtr->addSys();
      partonSystemsPtr->setInRes( iSys, iWZ);
      partonSystemsPtr->setSHat( iSys, pow2(event[iWZ].m()) );
      partonSystemsPtr->setPTHat( iSys, pTmax);
      for (int i = sizeSave; i < event.size(); ++i)
        partonSystemsPtr->addOut( iSys, i);

      // Do parton showers inside subsystem.
      if (doFSRinResonances) {

        // Let prepare routine do the setup.
        timesDecPtr->prepare( iSys, event);

        // Begin evolution down in pT from hard pT scale.
        do {
          double pTtimes = timesDecPtr->pTnext( event, pTmax, 0.);

          // Do a final-state emission (if allowed).
          if (pTtimes > 0.) {
            timesDecPtr->branch( event);
            pTmax = pTtimes;
          }

          // If no pT scales above zero then nothing to be done.
          else pTmax = 0.;

          // Keep on evolving until nothing is left to be done.
        } while (pTmax > 0.);
      }
    }

  // End loop over event to find W/Z gauge bosons.
  }

  // Done.
  return true;

}

//==========================================================================

} // end namespace Pythia8
