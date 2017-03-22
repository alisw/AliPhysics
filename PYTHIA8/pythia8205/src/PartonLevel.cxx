// PartonLevel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

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
  Couplings* couplingsPtrIn, PartonSystems* partonSystemsPtrIn,
  SigmaTotal* sigmaTotPtr, TimeShower* timesDecPtrIn, TimeShower* timesPtrIn,
  SpaceShower* spacePtrIn, RHadrons* rHadronsPtrIn, UserHooks* userHooksPtrIn,
  MergingHooks* mergingHooksPtrIn, bool useAsTrial ) {

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
  couplingsPtr       = couplingsPtrIn;
  partonSystemsPtr   = partonSystemsPtrIn;
  timesDecPtr        = timesDecPtrIn;
  timesPtr           = timesPtrIn;
  spacePtr           = spacePtrIn;
  rHadronsPtr        = rHadronsPtrIn;
  userHooksPtr       = userHooksPtrIn;
  mergingHooksPtr    = mergingHooksPtrIn;

  // Min bias and diffraction processes need special treatment.
  bool doSQ          = settings.flag("SoftQCD:all")
                    || settings.flag("SoftQCD:inelastic");
  bool doND          = settings.flag("SoftQCD:nonDiffractive");
  bool doSD          = settings.flag("SoftQCD:singleDiffractive");
  bool doDD          = settings.flag("SoftQCD:doubleDiffractive");
  bool doCD          = settings.flag("SoftQCD:centralDiffractive");
  doNonDiff          =  doSQ || doND;
  doDiffraction      =  doSQ || doSD || doDD || doCD;

  // Separate low-mass (unresolved) and high-mass (perturbative) diffraction.
  mMinDiff           = settings.parm("Diffraction:mMinPert");
  mWidthDiff         = settings.parm("Diffraction:mWidthPert");
  pMaxDiff           = settings.parm("Diffraction:probMaxPert");
  if (mMinDiff > infoPtr->eCM()) doDiffraction = false;

  // Need MPI initialization for soft QCD processes, even if only first MPI.
  // But no need to initialize MPI if never going to use it.
  doMPI              = settings.flag("PartonLevel:MPI");
  doMPIMB            = doMPI;
  doMPISDA           = doMPI;
  doMPISDB           = doMPI;
  doMPICD            = doMPI;
  doMPIinit          = doMPI;
  if (doNonDiff || doDiffraction)        doMPIinit = true;
  if (!settings.flag("PartonLevel:all")) doMPIinit = false;

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

  // Some other flags.
  doRemnants         = settings.flag("PartonLevel:Remnants");
  doSecondHard       = settings.flag("SecondHard:generate");
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

  // Flag if lepton beams, and if non-resolved ones. May change main flags.
  hasLeptonBeams     = ( beamAPtr->isLepton() || beamBPtr->isLepton() );
  hasPointLeptons    = ( hasLeptonBeams
    && (beamAPtr->isUnresolved() || beamBPtr->isUnresolved() ) );
  if (hasLeptonBeams) {
    doMPIMB          = false;
    doMPISDA         = false;
    doMPISDB         = false;
    doMPICD          = false;
    doMPIinit        = false;
  }
  if (hasPointLeptons) {
    doISR            = false;
    doRemnants       = false;
  }

  // Set info and initialize the respective program elements.
  timesPtr->init( beamAPtr, beamBPtr);
  if (doISR) spacePtr->init( beamAPtr, beamBPtr);
  doMPIMB  =  multiMB.init( doMPIinit, 0, infoPtr, settings, particleDataPtr,
    rndmPtr, beamAPtr, beamBPtr, couplingsPtr, partonSystemsPtr, sigmaTotPtr,
    userHooksPtr);
  if (doSD || doDD || doSQ) doMPISDA = multiSDA.init( doMPIinit, 1, infoPtr,
    settings, particleDataPtr, rndmPtr, beamAPtr, beamPomBPtr, couplingsPtr,
    partonSystemsPtr, sigmaTotPtr, userHooksPtr);
  if (doSD || doDD || doSQ) doMPISDB = multiSDB.init( doMPIinit, 2, infoPtr,
    settings, particleDataPtr, rndmPtr, beamPomAPtr, beamBPtr, couplingsPtr,
    partonSystemsPtr, sigmaTotPtr, userHooksPtr);
  if (doCD || doSQ) doMPICD = multiCD.init( doMPIinit, 3, infoPtr, settings,
    particleDataPtr, rndmPtr, beamPomAPtr, beamPomBPtr, couplingsPtr,
    partonSystemsPtr, sigmaTotPtr, userHooksPtr);
  if (!remnants.init( infoPtr, settings, rndmPtr, beamAPtr, beamBPtr,
    partonSystemsPtr, particleDataPtr, &colourReconnection)) return false;
  resonanceDecays.init( infoPtr, particleDataPtr, rndmPtr);
  colourReconnection.init( infoPtr, settings, rndmPtr, beamAPtr, beamBPtr,
    partonSystemsPtr);
  junctionSplitting.init(infoPtr, settings, rndmPtr, particleDataPtr);

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

  // Clear last branching return values.
  pTLastBranch   = 0.0;
  typeLastBranch = 0;

}

//--------------------------------------------------------------------------

// Main routine to do the parton-level evolution.

bool PartonLevel::next( Event& process, Event& event) {

  // Current event classification.
  isResolved     = infoPtr->isResolved();
  isResolvedA    = isResolved;
  isResolvedB    = isResolved;
  isResolvedC    = isResolved;
  isDiffA        = infoPtr->isDiffractiveA();
  isDiffB        = infoPtr->isDiffractiveB();
  isDiffC        = infoPtr->isDiffractiveC();
  isDiff         = isDiffA || isDiffB || isDiffC;
  isCentralDiff  = isDiffC;
  isDoubleDiff   = isDiffA && isDiffB;
  isSingleDiff   = isDiff && !isDoubleDiff  && !isCentralDiff;
  isNonDiff      = infoPtr->isNonDiffractive();
  doVeto         = false;
  infoPtr->setAbortPartonLevel(false);

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
    if (!physical || nHardLoop == 0) return physical;
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
  if (doMPIinit) multiPtr->reset();

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
    nMPI       = (doSecondHard) ? 2 : 1;
    nISR       = 0;
    nFSRinProc = 0;
    nFSRinRes  = 0;
    nISRhard   = 0;
    nFSRhard   = 0;
    pTsaveMPI  = 0.;
    pTsaveISR  = 0.;
    pTsaveFSR  = 0.;

    // Identify hard interaction system for showers.
    setupHardSys( process, event);

    // Optionally check for a veto after the hardest interaction.
    if (canVetoMPIStep) {
      doVeto = userHooksPtr->doVetoMPIStep( 1, event);
      // Abort event if vetoed.
      if (doVeto) {
        if (isDiff) leaveResolvedDiff( iHardLoop, process, event);
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
    double pTscaleMPI  = pTscaleRad;
    if (doSecondHard) {
      pTscaleRad       = max( pTscaleRad, process.scaleSecond() );
      pTscaleMPI       = min( pTscaleMPI, process.scaleSecond() );
    }
    double pTmaxMPI = (limitPTmaxMPI)  ? pTscaleMPI : infoPtr->eCM();
    double pTmaxISR = (limitPTmaxISR) ? spacePtr->enhancePTmax() * pTscaleRad
                                      : infoPtr->eCM();
    double pTmaxFSR = (limitPTmaxFSR) ? timesPtr->enhancePTmax() * pTscaleRad
                                      : infoPtr->eCM();

    // Potentially reset up starting scales for matrix element merging.
    if ( hasMergingHooks && (doTrial || canRemoveEvent || canRemoveEmission) )
      mergingHooksPtr->setShowerStartingScales( doTrial,
        (canRemoveEvent || canRemoveEmission), pTscaleRad, process, pTmaxFSR,
        limitPTmaxFSR, pTmaxISR, limitPTmaxISR, pTmaxMPI, limitPTmaxMPI );

    double pTmax    = max( pTmaxMPI, max( pTmaxISR, pTmaxFSR) );
    pTsaveMPI       = pTmaxMPI;
    pTsaveISR       = pTmaxISR;
    pTsaveFSR       = pTmaxFSR;
    // Prepare the classes to begin the generation.
    if (doMPI) multiPtr->prepare( event, pTmaxMPI);
    if (doISR) spacePtr->prepare( 0, event, limitPTmaxISR);
    if (doFSRduringProcess) timesPtr->prepare( 0, event, limitPTmaxFSR);
    if (doSecondHard && doISR) spacePtr->prepare( 1, event, limitPTmaxISR);
    if (doSecondHard && doFSRduringProcess) timesPtr->prepare( 1, event,
       limitPTmaxFSR);

    // Impact parameter has now been chosen, except for diffraction.
    if (!isDiff) infoPtr->setImpact( multiPtr->bMPI(),
      multiPtr->enhanceMPI(), true);
    // Set up initial veto scale.
    doVeto        = false;
    double pTveto = pTvetoPT;
    typeLatest    = 0;

    // Begin evolution down in pT from hard pT scale.
    do {
      infoPtr->addCounter(22);
      typeVetoStep = 0;
      nRad         =  nISR + nFSRinProc;

      // Find next pT value for FSR, MPI and ISR.
      // Order calls to minimize time expenditure.
      double pTgen = 0.;
      double pTtimes = (doFSRduringProcess)
        ? timesPtr->pTnext( event, pTmaxFSR, pTgen, isFirstTrial) : -1.;
      pTgen = max( pTgen, pTtimes);
      double pTmulti = (doMPI)
        ? multiPtr->pTnext( pTmaxMPI, pTgen, event) : -1.;
      pTgen = max( pTgen, pTmulti);
      double pTspace = (doISR)
        ? spacePtr->pTnext( event, pTmaxISR, pTgen, nRad) : -1.;
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

          // Update ISR and FSR dipoles.
          if (doISR)              spacePtr->prepare( nMPI - 1, event);
          if (doFSRduringProcess) timesPtr->prepare( nMPI - 1, event);
        }

        // Set maximal scales for next pT to pick.
        pTmaxMPI = pTmulti;
        pTmaxISR = min(pTmulti, pTmaxISR);
        pTmaxFSR = min(pTmulti, pTmaxFSR);
        pTmax    = pTmulti;
        nBranch++;
        pTLastBranch = pTmulti;
        typeLastBranch = 1;
      }

      // Do an initial-state emission (if allowed).
      else if (pTspace > 0. && pTspace > pTtimes) {
        infoPtr->addCounter(24);
        if (spacePtr->branch( event)) {
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
        pTmaxMPI = min(pTspace, pTmaxMPI);
        pTmaxISR = pTspace;
        pTmaxFSR = min(pTspace, pTmaxFSR);
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
        pTmaxMPI = min(pTtimes, pTmaxMPI);
        pTmaxISR = min(pTtimes, pTmaxISR);
        pTmaxFSR = pTtimes;
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
        return false;
      }

      // Keep on evolving until nothing is left to be done.
      if (typeLatest > 0 && typeLatest < 4)
        infoPtr->addCounter(25 + typeLatest);
      if (!isDiff) infoPtr->setPartEvolved( nMPI, nISR);

      // Handle potential merging veto.
      if ( canRemoveEvent && nISRhard + nFSRhard == 1 ){
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
            return false;
          }
        }

        // Handle potential merging veto.
        if ( canRemoveEvent && nISRhard + nFSRhard == 1 ){
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

      // Reassign new decay products to original system.
      for (int iSys = oldSizeSys; iSys < partonSystemsPtr->sizeSys(); ++iSys)
        for (int iOut = 0; iOut < partonSystemsPtr->sizeOut(iSys); ++iOut)
          partonSystemsPtr->addOut(0, partonSystemsPtr->getOut( iSys, iOut) );
      partonSystemsPtr->setSizeSys( oldSizeSys);

      // Perform decays and showers of W and Z emitted in shower.
      // To do:check if W/Z emission is on in ISR or FSR??
      if (!wzDecayShowers( event)) return false;

      // User hook to reconnect colours specifically in resonance decays.
      if (canReconResSys && !userHooksPtr->doReconnectResonanceSystems(
        oldSizeEvt, event)) return false;
    }

    // Find the first particle in the current diffractive system.
    int iFirst = 0;
    if (isDiff) {
      iFirst = (iHardLoop == 1) ? 5 + sizeEvent - sizeProcess : sizeEvent;
      if (isDiffC)  iFirst = 6 + sizeEvent - sizeProcess;
    }

    // Add beam remnants, including primordial kT kick and colour tracing.
    if (!doTrial && physical && doRemnants
      && !remnants.add( event, iFirst, isDiff)) physical = false;

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

  // End loop over ten tries. Restore from diffraction. Hopefully it worked.
  }
  if (isDiff) leaveResolvedDiff( iHardLoop, process, event);
  if (!physical) return false;

  // End big outer loop to handle two systems in double diffraction.
  }

  // Do colour reconnection for non-diffractive events before resonance decays.
  if (doReconnect && !isDiff && reconnectMode != 0) {
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
    infoPtr->setImpact( multiPtr->bMPI(), multiPtr->enhanceMPI(), false);
  }

  // Do colour reconnection for resonance decays.
  if (doReconnect && !isDiff && reconnectMode != 0) {
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
    int iDiffMot = iDSnow + 2;

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
    int iBeam = iDS + 2;

    // Diffractive mass. Boost and rotation from diffractive system
    // rest frame, aligned along z axis, to event cm frame.
    double mDiff  = process[iBeam].m();
    double m2Diff = mDiff * mDiff;
    Vec4 pDiffA   = (iDS == 1) ? process[1].p()
                               : process[1].p() - process[3].p();
    Vec4 pDiffB   = (iDS == 2) ? process[2].p()
                               : process[2].p() - process[4].p();
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
  if (doSecondHard) {
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
  beamAPtr->append( inP + nOffset, process[inP].id(), x1);
  double x2 = process[inM].pNeg() / process[inS].m();
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
  else vsc1 = beamAPtr->pickValSeaComp();
  beamBPtr->xfISR( 0, process[inM].id(), x2, scale*scale);
  if (isNonDiff && (vsc2 = multiPtr->getVSC2()) != 0)
    (*beamBPtr)[0].companion(vsc2);
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
  if (doSecondHard) {
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
    if (doCopy) event.appendJunction( process.getJunction(iJun));
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
    if (doCopy) event.appendJunction( process.getJunction(iJun));
  }

  // Done.
}

//--------------------------------------------------------------------------

// Resolved diffraction: replace full event with diffractive subsystem.

void PartonLevel::setupResolvedDiff( Event& process) {

  // Mother and last entry of diffractive system.
  int iDiffMot     = iDS + 2;
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
  int idDiffA    = (iDS == 1) ? process[1].id() : 990;
  int idDiffB    = (iDS == 2) ? process[2].id() : 990;
  double mDiffA  = (iDS == 1) ? process[1].m() : 0.;
  double mDiffB  = (iDS == 2) ? process[2].m() : 0.;
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
  beamAPtr       = (iDS == 1) ? beamHadAPtr : beamPomAPtr;
  beamBPtr       = (iDS == 2) ? beamHadBPtr : beamPomBPtr;

  // Pretend that the diffractive system is the whole collision.
  eCMsave = infoPtr->eCM();
  infoPtr->setECM( mDiff);
  beamAPtr->newPzE(  pzDiff, eDiffA);
  beamBPtr->newPzE( -pzDiff, eDiffB);

  // Beams not found in normal slots 1 and 2.
  int beamOffset = (sizeEvent > 0) ? sizeEvent - 1 : 4;

  // Reassign beam pointers in other classes.
  timesPtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  spacePtr->reassignBeamPtrs( beamAPtr, beamBPtr, beamOffset);
  remnants.reassignBeamPtrs(  beamAPtr, beamBPtr, iDS);
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
  Vec4 pDiffA = (iDS == 1) ? process[1].p()
                           : process[1].p() - process[3].p();
  Vec4 pDiffB = (iDS == 2) ? process[2].p()
                           : process[2].p() - process[4].p();
  RotBstMatrix MtoCM;
  MtoCM.fromCMframe( pDiffA, pDiffB);

  // Perform rotation and boost on diffractive system.
  for (int i = sizeProcess; i < process.size(); ++i)
    process[i].rotbst( MtoCM);
  int iFirst = (iHardLoop == 1) ? 5 + sizeEvent - sizeProcess : sizeEvent;
  if(isDiffC)  iFirst = 6 + sizeEvent - sizeProcess;
  for (int i = iFirst; i < event.size(); ++i)
    event[i].rotbst( MtoCM);

  // Restore cm energy.
  infoPtr->setECM( eCMsave);
  beamAPtr->newPzE( event[1].pz(), event[1].e());
  beamBPtr->newPzE( event[2].pz(), event[2].e());

  // Restore beam pointers to incoming hadrons.
  beamAPtr = beamHadAPtr;
  beamBPtr = beamHadBPtr;

  // Reassign beam pointers in other classes.
  timesPtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  spacePtr->reassignBeamPtrs( beamAPtr, beamBPtr, 0);
  remnants.reassignBeamPtrs(  beamAPtr, beamBPtr, 0);
  colourReconnection.reassignBeamPtrs(  beamAPtr, beamBPtr);

  // Restore multiparton interactions pointer to default object.
  multiPtr = &multiMB;

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

    // In first call (skipForR = true) skip over resonances
    // that should form R-hadrons, and their daughters.
    if (allowRH) {
      if (skipForR) {
        bool comesFromR = false;
        int iTraceUp = iBegin;
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
      && !earlyResDec) {
      infoPtr->errorMsg("Abort in PartonLevel::resonanceShower: "
        "new CR can't handle late resonance decay of coloured particles");
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

      // Check if this parton carries a junction color in hard event.
      for (int iJun = 0; iJun < process.sizeJunction(); ++iJun) {
        if (iJun >= int(doCopyJun.size())) doCopyJun.push_back(false);
        // Only consider junctions that can appear in decays.
        int kindJunction = process.kindJunction(iJun);
        if (kindJunction >= 5) continue;
        int col = (kindJunction % 2 == 1) ? now.col() : now.acol();
        int iLegF1 = (kindJunction - 1) / 2;
        for (int iLeg = iLegF1; iLeg <= 2; ++iLeg)
        if (col == process.colJunction(iJun,iLeg)) doCopyJun[iJun] = true;
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

    // Add new system, automatically with two empty beam slots.
    int iSys = partonSystemsPtr->addSys();
    partonSystemsPtr->setSHat(iSys, pow2(hardMother.m()) );
    partonSystemsPtr->setPTHat(iSys, 0.5 * hardMother.m() );

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

    // Add new system, automatically with two empty beam slots.
    if (typeWZ > 0) {
      // Maximum shower scale set by mother mass.
      double pTmax = 0.5 * event[iWZ].m();
      int iSys = partonSystemsPtr->addSys();
      partonSystemsPtr->setSHat(iSys, pow2(event[iWZ].m()) );
      partonSystemsPtr->setPTHat(iSys, pTmax );
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
