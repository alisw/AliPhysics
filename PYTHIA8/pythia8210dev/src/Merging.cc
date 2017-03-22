// MergingHooks.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel.
// Function definitions (not found in the header) for the Merging class.

#include "Pythia8/Merging.h"

namespace Pythia8 {

//==========================================================================

// The Merging class.

//--------------------------------------------------------------------------

// Factor by which the maximal value of the merging scale can deviate before
// a warning is printed.
const double Merging::TMSMISMATCH = 1.5;

//--------------------------------------------------------------------------

// Initialise Merging class

void Merging::init( Settings* settingsPtrIn, Info* infoPtrIn,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
  BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  MergingHooks* mergingHooksPtrIn, PartonLevel* trialPartonLevelPtrIn ){

  // Save pointers.
  settingsPtr           = settingsPtrIn;
  infoPtr               = infoPtrIn;
  particleDataPtr       = particleDataPtrIn;
  rndmPtr               = rndmPtrIn;
  mergingHooksPtr       = mergingHooksPtrIn;
  trialPartonLevelPtr   = trialPartonLevelPtrIn;
  beamAPtr              = beamAPtrIn;
  beamBPtr              = beamBPtrIn;
  // Reset minimal tms value.
  tmsNowMin             = infoPtr->eCM();

}

//--------------------------------------------------------------------------

// Function to print information.
void Merging::statistics( ostream& os ) {

  // Recall switch to enfore merging scale cut.
  bool enforceCutOnLHE  = settingsPtr->flag("Merging:enforceCutOnLHE");
  // Recall merging scale value.
  double tmsval         = mergingHooksPtr->tms();
  bool printBanner      = enforceCutOnLHE && tmsNowMin > TMSMISMATCH*tmsval;
  // Reset minimal tms value.
  tmsNowMin             = infoPtr->eCM();

  if (!printBanner) return;

  // Header.
  os << "\n *-------  PYTHIA Matrix Element Merging Information  ------"
     << "-------------------------------------------------------*\n"
     << " |                                                            "
     << "                                                     |\n";
  // Print warning if the minimal tms value of any event was significantly
  // above the desired merging scale value.
  os << " | Warning in Merging::statistics: All Les Houches events"
     << " significantly above Merging:TMS cut. Please check.       |\n";

  // Listing finished.
  os << " |                                                            "
     << "                                                     |\n"
     << " *-------  End PYTHIA Matrix Element Merging Information -----"
     << "-----------------------------------------------------*" << endl;
}

//--------------------------------------------------------------------------

// Function to steer different merging prescriptions.

int Merging::mergeProcess(Event& process){

  int vetoCode = 1;

  // Reinitialise hard process.
  mergingHooksPtr->hardProcess.clear();
  mergingHooksPtr->processSave = settingsPtr->word("Merging:Process");
  mergingHooksPtr->hardProcess.initOnProcess(
    settingsPtr->word("Merging:Process"), particleDataPtr);

  mergingHooksPtr->doUserMergingSave
    = settingsPtr->flag("Merging:doUserMerging");
  mergingHooksPtr->doMGMergingSave
    = settingsPtr->flag("Merging:doMGMerging");
  mergingHooksPtr->doKTMergingSave
    = settingsPtr->flag("Merging:doKTMerging");
  mergingHooksPtr->doPTLundMergingSave
    = settingsPtr->flag("Merging:doPTLundMerging");
  mergingHooksPtr->doCutBasedMergingSave
    = settingsPtr->flag("Merging:doCutBasedMerging");
  mergingHooksPtr->doNL3TreeSave
    = settingsPtr->flag("Merging:doNL3Tree");
  mergingHooksPtr->doNL3LoopSave
    = settingsPtr->flag("Merging:doNL3Loop");
  mergingHooksPtr->doNL3SubtSave
    = settingsPtr->flag("Merging:doNL3Subt");
  mergingHooksPtr->doUNLOPSTreeSave
    = settingsPtr->flag("Merging:doUNLOPSTree");
  mergingHooksPtr->doUNLOPSLoopSave
    = settingsPtr->flag("Merging:doUNLOPSLoop");
  mergingHooksPtr->doUNLOPSSubtSave
    = settingsPtr->flag("Merging:doUNLOPSSubt");
  mergingHooksPtr->doUNLOPSSubtNLOSave
    = settingsPtr->flag("Merging:doUNLOPSSubtNLO");
  mergingHooksPtr->doUMEPSTreeSave
    = settingsPtr->flag("Merging:doUMEPSTree");
  mergingHooksPtr->doUMEPSSubtSave
    = settingsPtr->flag("Merging:doUMEPSSubt");
  mergingHooksPtr->nReclusterSave
    = settingsPtr->mode("Merging:nRecluster");

  // Possibility to apply merging scale to an input event.
  bool applyTMSCut = settingsPtr->flag("Merging:doXSectionEstimate");
  if ( applyTMSCut && cutOnProcess(process) ) return -1;
  // Done if only a cut should be applied.
  if ( applyTMSCut ) return 1;

  // Possibility to perform CKKW-L merging on this event.
  if ( mergingHooksPtr->doCKKWLMerging() )
    vetoCode = mergeProcessCKKWL(process);

  // Possibility to perform UMEPS merging on this event.
  if ( mergingHooksPtr->doUMEPSMerging() )
     vetoCode = mergeProcessUMEPS(process);

  // Possibility to perform NL3 NLO merging on this event.
  if ( mergingHooksPtr->doNL3Merging() )
    vetoCode = mergeProcessNL3(process);

  // Possibility to perform UNLOPS merging on this event.
  if ( mergingHooksPtr->doUNLOPSMerging() )
    vetoCode = mergeProcessUNLOPS(process);

  return vetoCode;

}

//--------------------------------------------------------------------------

// Function to perform CKKW-L merging on this event.

int Merging::mergeProcessCKKWL( Event& process) {

  // Ensure that merging hooks to not veto events in the trial showers.
  mergingHooksPtr->doIgnoreStep(true);
  // For pp > h, allow cut on state, so that underlying processes
  // can be clustered to gg > h
  if ( mergingHooksPtr->getProcessString().compare("pp>h") == 0 )
    mergingHooksPtr->allowCutOnRecState(true);
  // For now, prefer construction of ordered histories.
  mergingHooksPtr->orderHistories(true);

  // Reset weight of the event.
  double wgt = 1.0;
  mergingHooksPtr->setWeightCKKWL(1.);
  mergingHooksPtr->muMI(-1.);

  // Prepare process record for merging. If Pythia has already decayed
  // resonances used to define the hard process, remove resonance decay
  // products.
  Event newProcess( mergingHooksPtr->bareEvent( process, true) );
  // Store candidates for the splitting V -> qqbar'.
  mergingHooksPtr->storeHardProcessCandidates( newProcess);

  // Check if event passes the merging scale cut.
  double tmsval = mergingHooksPtr->tms();
  // Get merging scale in current event.
  double tmsnow = mergingHooksPtr->tmsNow( newProcess );
  // Calculate number of clustering steps.
  int nSteps    = mergingHooksPtr->getNumberOfClusteringSteps( newProcess);

  // Too few steps can be possible if a chain of resonance decays has been
  // removed. In this case, reject this event, since it will be handled in
  // lower-multiplicity samples.
  int nRequested = settingsPtr->mode("Merging:nRequested");
  if (nSteps < nRequested) {
    mergingHooksPtr->setWeightCKKWL(0.);
    return -1;
  }

  // Reset the minimal tms value, if necessary.
  tmsNowMin     = (nSteps == 0) ? 0. : min(tmsNowMin, tmsnow);

  // Get random number to choose a path.
  double RN = rndmPtr->flat();
  // Set dummy process scale.
  newProcess.scale(0.0);
  // Generate all histories.
  History FullHistory( nSteps, 0.0, newProcess, Clustering(), mergingHooksPtr,
            (*beamAPtr), (*beamBPtr), particleDataPtr, infoPtr, true, true,
            true, true, 1.0, 0);
  // Project histories onto desired branches, e.g. only ordered paths.
  FullHistory.projectOntoDesiredHistories();

  // Do not apply cut if the configuration could not be projected onto an
  // underlying born configuration.
  bool applyCut = nSteps > 0 && FullHistory.select(RN)->nClusterings() > 0;

  // Enfore merging scale cut if the event did not pass the merging scale
  // criterion.
  bool enforceCutOnLHE  = settingsPtr->flag("Merging:enforceCutOnLHE");
  if ( enforceCutOnLHE && applyCut && tmsnow < tmsval ) {
    string message="Warning in Merging::mergeProcessCKKWL: Les Houches Event";
    message+=" fails merging scale cut. Reject event.";
    infoPtr->errorMsg(message);
    mergingHooksPtr->setWeightCKKWL(0.);
    return -1;
  }

  if ( FullHistory.select(RN)->nClusterings() < nSteps) {
    string message="Warning in Merging::mergeProcessCKKWL: No clusterings";
    message+=" found. History incomplete.";
    infoPtr->errorMsg(message);
  }

  // Calculate CKKWL weight:
  // Perform reweighting with Sudakov factors, save alpha_s ratios and
  // PDF ratio weights.
  wgt = FullHistory.weightTREE( trialPartonLevelPtr,
          mergingHooksPtr->AlphaS_FSR(), mergingHooksPtr->AlphaS_ISR(), RN);

  // Event with production scales set for further (trial) showering
  // and starting conditions for the shower.
  FullHistory.getStartingConditions( RN, process );
  // If necessary, reattach resonance decay products.
  mergingHooksPtr->reattachResonanceDecays(process);

  // Allow to dampen histories in which the lowest multiplicity reclustered
  // state does not pass the lowest multiplicity cut of the matrix element.
  double dampWeight = mergingHooksPtr->dampenIfFailCuts(
           FullHistory.lowestMultProc(RN) );
  // Save the weight of the event for histogramming. Only change the
  // event weight after trial shower on the matrix element
  // multiplicity event (= in doVetoStep).
  wgt *= dampWeight;

  // Set QCD 2->2 starting scale different from arbitrary scale in LHEF!
  // --> Set to minimal mT of partons.
  int nFinal = 0;
  double muf = process[0].e();
  for ( int i=0; i < process.size(); ++i )
  if ( process[i].isFinal()
    && (process[i].colType() != 0 || process[i].id() == 22 ) ) {
    nFinal++;
    muf = min( muf, abs(process[i].mT()) );
  }
  // For pure QCD dijet events (only!), set the process scale to the
  // transverse momentum of the outgoing partons.
  if ( nSteps == 0 && nFinal == 2
    && ( mergingHooksPtr->getProcessString().compare("pp>jj") == 0
      || mergingHooksPtr->getProcessString().compare("pp>aj") == 0) )
    process.scale(muf);

  // Save the weight of the event for histogramming.
  mergingHooksPtr->setWeightCKKWL(wgt);

  // Allow merging hooks to veto events from now on.
  mergingHooksPtr->doIgnoreStep(false);

  // If no-emission probability is zero.
  if ( wgt == 0. ) return 0;

  // Done
  return 1;

}

//--------------------------------------------------------------------------

// Function to perform UMEPS merging on this event.

int Merging::mergeProcessUMEPS( Event& process) {

  // Initialise which part of UMEPS merging is applied.
  bool doUMEPSTree                = settingsPtr->flag("Merging:doUMEPSTree");
  bool doUMEPSSubt                = settingsPtr->flag("Merging:doUMEPSSubt");
  // Save number of looping steps
  mergingHooksPtr->nReclusterSave = settingsPtr->mode("Merging:nRecluster");
  int nRecluster                  = settingsPtr->mode("Merging:nRecluster");

  // Ensure that merging hooks does not remove emissions.
  mergingHooksPtr->doIgnoreEmissions(true);
  // For pp > h, allow cut on state, so that underlying processes
  // can be clustered to gg > h
  if ( mergingHooksPtr->getProcessString().compare("pp>h") == 0 )
    mergingHooksPtr->allowCutOnRecState(true);
  // For now, prefer construction of ordered histories.
  mergingHooksPtr->orderHistories(true);

  // Reset weights of the event.
  double wgt   = 1.;
  mergingHooksPtr->setWeightCKKWL(1.);
  mergingHooksPtr->muMI(-1.);

  // Prepare process record for merging. If Pythia has already decayed
  // resonances used to define the hard process, remove resonance decay
  // products.
  Event newProcess( mergingHooksPtr->bareEvent( process, true) );
  // Store candidates for the splitting V -> qqbar'.
  mergingHooksPtr->storeHardProcessCandidates( newProcess );

  // Check if event passes the merging scale cut.
  double tmsval   = mergingHooksPtr->tms();
  // Get merging scale in current event.
  double tmsnow  = mergingHooksPtr->tmsNow( newProcess );
  // Calculate number of clustering steps.
  int nSteps     = mergingHooksPtr->getNumberOfClusteringSteps( newProcess );

  // Too few steps can be possible if a chain of resonance decays has been
  // removed. In this case, reject this event, since it will be handled in
  // lower-multiplicity samples.
  int nRequested = settingsPtr->mode("Merging:nRequested");
  if (nSteps < nRequested) {
    mergingHooksPtr->setWeightCKKWL(0.);
    return -1;
  }

  // Reset the minimal tms value, if necessary.
  tmsNowMin      = (nSteps == 0) ? 0. : min(tmsNowMin, tmsnow);

  // Get random number to choose a path.
  double RN = rndmPtr->flat();
  // Set dummy process scale.
  newProcess.scale(0.0);
  // Generate all histories.
  History FullHistory( nSteps, 0.0, newProcess, Clustering(), mergingHooksPtr,
            (*beamAPtr), (*beamBPtr), particleDataPtr, infoPtr, true, true,
            true, true, 1.0, 0);
  // Project histories onto desired branches, e.g. only ordered paths.
  FullHistory.projectOntoDesiredHistories();

  // Do not apply cut if the configuration could not be projected onto an
  // underlying born configuration.
  bool applyCut = nSteps > 0 && FullHistory.select(RN)->nClusterings() > 0;

  // Enfore merging scale cut if the event did not pass the merging scale
  // criterion.
  bool enforceCutOnLHE  = settingsPtr->flag("Merging:enforceCutOnLHE");
  if ( enforceCutOnLHE && applyCut && tmsnow < tmsval ) {
    string message="Warning in Merging::mergeProcessUMEPS: Les Houches Event";
    message+=" fails merging scale cut. Reject event.";
    infoPtr->errorMsg(message);
    mergingHooksPtr->setWeightCKKWL(0.);
    return -1;
  }

  // Check reclustering steps to correctly apply MPI.
  int nPerformed = 0;
  if ( nSteps > 0 && doUMEPSSubt
    && !FullHistory.getFirstClusteredEventAboveTMS( RN, nRecluster,
          newProcess, nPerformed, false ) ) {
    // Discard if the state could not be reclustered to a state above TMS.
    mergingHooksPtr->setWeightCKKWL(0.);
    return -1;
  }

  mergingHooksPtr->nMinMPI(nSteps - nPerformed);

  // Calculate CKKWL weight:
  // Perform reweighting with Sudakov factors, save alpha_s ratios and
  // PDF ratio weights.
  if ( doUMEPSTree ) {
    wgt = FullHistory.weight_UMEPS_TREE( trialPartonLevelPtr,
      mergingHooksPtr->AlphaS_FSR(), mergingHooksPtr->AlphaS_ISR(), RN );
  } else {
    wgt = FullHistory.weight_UMEPS_SUBT( trialPartonLevelPtr,
      mergingHooksPtr->AlphaS_FSR(), mergingHooksPtr->AlphaS_ISR(), RN );
  }

  // Event with production scales set for further (trial) showering
  // and starting conditions for the shower.
  if ( doUMEPSTree ) FullHistory.getStartingConditions( RN, process );
  // Do reclustering (looping) steps.
  else FullHistory.getFirstClusteredEventAboveTMS( RN, nRecluster, process,
    nPerformed, true );

  // Allow to dampen histories in which the lowest multiplicity reclustered
  // state does not pass the lowest multiplicity cut of the matrix element
  double dampWeight = mergingHooksPtr->dampenIfFailCuts(
           FullHistory.lowestMultProc(RN) );
  // Save the weight of the event for histogramming. Only change the
  // event weight after trial shower on the matrix element
  // multiplicity event (= in doVetoStep)
  wgt *= dampWeight;

  // Save the weight of the event for histogramming.
  mergingHooksPtr->setWeightCKKWL(wgt);

  // Set QCD 2->2 starting scale different from arbitrary scale in LHEF!
  // --> Set to minimal mT of partons.
  int nFinal = 0;
  double muf = process[0].e();
  for ( int i=0; i < process.size(); ++i )
  if ( process[i].isFinal()
    && (process[i].colType() != 0 || process[i].id() == 22 ) ) {
    nFinal++;
    muf = min( muf, abs(process[i].mT()) );
  }

  // For pure QCD dijet events (only!), set the process scale to the
  // transverse momentum of the outgoing partons.
  // Calculate number of clustering steps.
  int nStepsNew = mergingHooksPtr->getNumberOfClusteringSteps( process );
  if ( nStepsNew == 0
    && ( mergingHooksPtr->getProcessString().compare("pp>jj") == 0
      || mergingHooksPtr->getProcessString().compare("pp>aj") == 0) )
    process.scale(muf);

  // Reset hard process candidates (changed after clustering a parton).
  mergingHooksPtr->storeHardProcessCandidates( process );
  // If necessary, reattach resonance decay products.
  mergingHooksPtr->reattachResonanceDecays(process);

  // Allow merging hooks to remove emissions from now on.
  mergingHooksPtr->doIgnoreEmissions(false);

  // If no-emission probability is zero.
  if ( wgt == 0. ) return 0;

  // Done
  return 1;

}

//--------------------------------------------------------------------------

// Function to perform NL3 NLO merging on this event.

int Merging::mergeProcessNL3( Event& process) {

  // Initialise which part of NL3 merging is applied.
  bool doNL3Tree = settingsPtr->flag("Merging:doNL3Tree");
  bool doNL3Loop = settingsPtr->flag("Merging:doNL3Loop");
  bool doNL3Subt = settingsPtr->flag("Merging:doNL3Subt");
  // Save number of looping steps.
  int nRequested = settingsPtr->mode("Merging:nRequested");

  // Ensure that hooks (NL3 part) to not remove emissions.
  mergingHooksPtr->doIgnoreEmissions(true);
  // Ensure that hooks (CKKWL part) to not veto events in trial showers.
  mergingHooksPtr->doIgnoreStep(true);
  // For pp > h, allow cut on state, so that underlying processes
  // can be clustered to gg > h
  if ( mergingHooksPtr->getProcessString().compare("pp>h") == 0)
    mergingHooksPtr->allowCutOnRecState(true);
  // For now, prefer construction of ordered histories.
  mergingHooksPtr->orderHistories(true);

  // Reset weight of the event
  double wgt      = 1.;
  mergingHooksPtr->setWeightCKKWL(1.);
  // Reset the O(alphaS)-term of the CKKW-L weight.
  double wgtFIRST = 0.;
  mergingHooksPtr->setWeightFIRST(0.);
  mergingHooksPtr->muMI(-1.);

  // Prepare process record for merging. If Pythia has already decayed
  // resonances used to define the hard process, remove resonance decay
  // products.
  Event newProcess( mergingHooksPtr->bareEvent( process, true) );
  // Store candidates for the splitting V -> qqbar'
  mergingHooksPtr->storeHardProcessCandidates( newProcess);

  // Check if event passes the merging scale cut.
  double tmsval  = mergingHooksPtr->tms();
  // Get merging scale in current event.
  double tmsnow  = mergingHooksPtr->tmsNow( newProcess );
  // Calculate number of clustering steps
  int nSteps   = mergingHooksPtr->getNumberOfClusteringSteps( newProcess);

  // Too few steps can be possible if a chain of resonance decays has been
  // removed. In this case, reject this event, since it will be handled in
  // lower-multiplicity samples.
  if (nSteps < nRequested) {
    mergingHooksPtr->setWeightCKKWL(0.);
    mergingHooksPtr->setWeightFIRST(0.);
    return -1;
  }

  // Reset the minimal tms value, if necessary.
  tmsNowMin = (nSteps == 0) ? 0. : min(tmsNowMin, tmsnow);

  // Enfore merging scale cut if the event did not pass the merging scale
  // criterion.
  bool enforceCutOnLHE  = settingsPtr->flag("Merging:enforceCutOnLHE");
  if ( enforceCutOnLHE && nSteps > 0 && nSteps == nRequested
    && tmsnow < tmsval ) {
    string message="Warning in Merging::mergeProcessNL3: Les Houches Event";
    message+=" fails merging scale cut. Reject event.";
    infoPtr->errorMsg(message);
    mergingHooksPtr->setWeightCKKWL(0.);
    mergingHooksPtr->setWeightFIRST(0.);
    return -1;
  }

  // Get random number to choose a path.
  double RN = rndmPtr->flat();
  // Set dummy process scale.
  newProcess.scale(0.0);
  // Generate all histories
  History FullHistory( nSteps, 0.0, newProcess, Clustering(), mergingHooksPtr,
            (*beamAPtr), (*beamBPtr), particleDataPtr, infoPtr, true, true,
            true, true, 1.0, 0);
  // Project histories onto desired branches, e.g. only ordered paths.
  FullHistory.projectOntoDesiredHistories();

  // Discard states that cannot be projected unto a state with one less jet.
  if ( nSteps > 0 && doNL3Subt
    && FullHistory.select(RN)->nClusterings() == 0 ){
    mergingHooksPtr->setWeightCKKWL(0.);
    mergingHooksPtr->setWeightFIRST(0.);
    return -1;
  }

  // Potentially recluster real emission jets for powheg input containing
  // "too many" jets, i.e. real-emission kinematics.
  bool containsRealKin = nSteps > nRequested && nSteps > 0;

  // Perform one reclustering for real emission kinematics, then apply merging
  // scale cut on underlying Born kinematics.
  if ( containsRealKin ) {
    Event dummy = Event();
    // Initialise temporary output of reclustering.
    dummy.clear();
    dummy.init( "(hard process-modified)", particleDataPtr );
    dummy.clear();
    // Recluster once.
    if ( !FullHistory.getClusteredEvent( RN, nSteps, dummy )) {
      mergingHooksPtr->setWeightCKKWL(0.);
      mergingHooksPtr->setWeightFIRST(0.);
      return -1;
    }
    double tnowNew  = mergingHooksPtr->tmsNow( dummy );
    // Veto if underlying Born kinematics do not pass merging scale cut.
    if ( enforceCutOnLHE && nSteps > 0 && nRequested > 0
      && tnowNew < tmsval ) {
      mergingHooksPtr->setWeightCKKWL(0.);
      mergingHooksPtr->setWeightFIRST(0.);
      return -1;
    }
  }

  // Remember number of jets, to include correct MPI no-emission probabilities.
  if ( doNL3Subt || containsRealKin ) mergingHooksPtr->nMinMPI(nSteps - 1);
  else mergingHooksPtr->nMinMPI(nSteps);

  // Calculate weight
  // Do LO or first part of NLO tree-level reweighting
  if( doNL3Tree ) {
    // Perform reweighting with Sudakov factors, save as ratios and
    // PDF ratio weights
    wgt = FullHistory.weightTREE( trialPartonLevelPtr,
      mergingHooksPtr->AlphaS_FSR(), mergingHooksPtr->AlphaS_ISR(), RN);
  } else if( doNL3Loop || doNL3Subt ) {
    // No reweighting, just set event scales properly and incorporate MPI
    // no-emission probabilities.
    wgt = FullHistory.weightLOOP( trialPartonLevelPtr, RN);
  }

  // Event with production scales set for further (trial) showering
  // and starting conditions for the shower
  if ( !doNL3Subt && !containsRealKin )
    FullHistory.getStartingConditions(RN, process);
  // For sutraction of nSteps-additional resolved partons from
  // the nSteps-1 parton phase space, recluster the last parton
  // in nSteps-parton events, and sutract later
  else {
    // Function to return the reclustered event
    if ( !FullHistory.getClusteredEvent( RN, nSteps, process )) {
      mergingHooksPtr->setWeightCKKWL(0.);
      mergingHooksPtr->setWeightFIRST(0.);
      return -1;
    }
  }

  // Allow to dampen histories in which the lowest multiplicity reclustered
  // state does not pass the lowest multiplicity cut of the matrix element
  double dampWeight = mergingHooksPtr->dampenIfFailCuts(
           FullHistory.lowestMultProc(RN) );
  // Save the weight of the event for histogramming. Only change the
  // event weight after trial shower on the matrix element
  // multiplicity event (= in doVetoStep)
  wgt *= dampWeight;

  // For tree level samples in NL3, rescale with k-Factor
  if (doNL3Tree ){
    // Find k-factor
    double kFactor = 1.;
    if( nSteps > mergingHooksPtr->nMaxJetsNLO() )
      kFactor = mergingHooksPtr->kFactor( mergingHooksPtr->nMaxJetsNLO() );
    else kFactor = mergingHooksPtr->kFactor(nSteps);
    // For NLO merging, rescale CKKW-L weight with k-factor
    wgt *= kFactor;
  }

  // Save the weight of the event for histogramming
  mergingHooksPtr->setWeightCKKWL(wgt);

  // Check if we need to subtract the O(\alpha_s)-term. If the number
  // of additional partons is larger than the number of jets for
  // which loop matrix elements are available, do standard CKKW-L
  bool doOASTree = doNL3Tree && nSteps <= mergingHooksPtr->nMaxJetsNLO();

  // Now begin NLO part for tree-level events
  if ( doOASTree ) {
    // Calculate the O(\alpha_s)-term of the CKKWL weight
    wgtFIRST = FullHistory.weightFIRST( trialPartonLevelPtr,
      mergingHooksPtr->AlphaS_FSR(), mergingHooksPtr->AlphaS_ISR(), RN,
      rndmPtr );
    // If necessary, also dampen the O(\alpha_s)-term
    wgtFIRST *= dampWeight;
    // Set the subtractive weight to the value calculated so far
    mergingHooksPtr->setWeightFIRST(wgtFIRST);
    // Subtract the O(\alpha_s)-term from the CKKW-L weight
    // If PDF contributions have not been included, subtract these later
    wgt = wgt - wgtFIRST;
  }

  // Set qcd 2->2 starting scale different from arbirtrary scale in LHEF!
  // --> Set to pT of partons
  double pT = 0.;
  for( int i=0; i < process.size(); ++i)
    if(process[i].isFinal() && process[i].colType() != 0) {
      pT = sqrt(pow(process[i].px(),2) + pow(process[i].py(),2));
      break;
    }
  // For pure QCD dijet events (only!), set the process scale to the
  // transverse momentum of the outgoing partons.
  if ( nSteps == 0
    && mergingHooksPtr->getProcessString().compare("pp>jj") == 0)
    process.scale(pT);

  // Reset hard process candidates (changed after clustering a parton).
  mergingHooksPtr->storeHardProcessCandidates( process );
  // If necessary, reattach resonance decay products.
  mergingHooksPtr->reattachResonanceDecays(process);

  // Allow merging hooks (NL3 part) to remove emissions from now on.
  mergingHooksPtr->doIgnoreEmissions(false);
  // Allow merging hooks (CKKWL part) to veto events from now on.
  mergingHooksPtr->doIgnoreStep(false);

  // Done
  return 1;

}

//--------------------------------------------------------------------------

// Function to perform UNLOPS merging on this event.

int Merging::mergeProcessUNLOPS( Event& process) {

  // Initialise which part of UNLOPS merging is applied.
  bool nloTilde         = settingsPtr->flag("Merging:doUNLOPSTilde");
  bool doUNLOPSTree     = settingsPtr->flag("Merging:doUNLOPSTree");
  bool doUNLOPSLoop     = settingsPtr->flag("Merging:doUNLOPSLoop");
  bool doUNLOPSSubt     = settingsPtr->flag("Merging:doUNLOPSSubt");
  bool doUNLOPSSubtNLO  = settingsPtr->flag("Merging:doUNLOPSSubtNLO");
  // Save number of looping steps
  mergingHooksPtr->nReclusterSave = settingsPtr->mode("Merging:nRecluster");
  int nRecluster        = settingsPtr->mode("Merging:nRecluster");
  int nRequested        = settingsPtr->mode("Merging:nRequested");

  // Ensure that merging hooks to not remove emissions
  mergingHooksPtr->doIgnoreEmissions(true);
  // For now, prefer construction of ordered histories.
  mergingHooksPtr->orderHistories(true);
  // For pp > h, allow cut on state, so that underlying processes
  // can be clustered to gg > h
  if ( mergingHooksPtr->getProcessString().compare("pp>h") == 0)
    mergingHooksPtr->allowCutOnRecState(true);

  // Reset weight of the event.
  double wgt      = 1.;
  mergingHooksPtr->setWeightCKKWL(1.);
  // Reset the O(alphaS)-term of the UMEPS weight.
  double wgtFIRST = 0.;
  mergingHooksPtr->setWeightFIRST(0.);
  mergingHooksPtr->muMI(-1.);

  // Prepare process record for merging. If Pythia has already decayed
  // resonances used to define the hard process, remove resonance decay
  // products.
  Event newProcess( mergingHooksPtr->bareEvent( process, true) );
  // Store candidates for the splitting V -> qqbar'
  mergingHooksPtr->storeHardProcessCandidates( newProcess );

  // Check if event passes the merging scale cut.
  double tmsval  = mergingHooksPtr->tms();
  // Get merging scale in current event.
  double tmsnow  = mergingHooksPtr->tmsNow( newProcess );
  // Calculate number of clustering steps
  int nSteps     = mergingHooksPtr->getNumberOfClusteringSteps( newProcess);

  // Too few steps can be possible if a chain of resonance decays has been
  // removed. In this case, reject this event, since it will be handled in
  // lower-multiplicity samples.
  if (nSteps < nRequested) {
    string message="Warning in Merging::mergeProcessUNLOPS: Les Houches Event";
    message+=" after removing decay products does not contain enough partons.";
    infoPtr->errorMsg(message);
    mergingHooksPtr->setWeightCKKWL(0.);
    mergingHooksPtr->setWeightFIRST(0.);
    return -1;
  }

  // Reset the minimal tms value, if necessary.
  tmsNowMin = (nSteps == 0) ? 0. : min(tmsNowMin, tmsnow);

  // Get random number to choose a path.
  double RN = rndmPtr->flat();
  // Set dummy process scale.
  newProcess.scale(0.0);
  // Generate all histories
  History FullHistory( nSteps, 0.0, newProcess, Clustering(), mergingHooksPtr,
            (*beamAPtr), (*beamBPtr), particleDataPtr, infoPtr, true, true,
            true, true, 1.0, 0);
  // Project histories onto desired branches, e.g. only ordered paths.
  FullHistory.projectOntoDesiredHistories();

  // Do not apply cut if the configuration could not be projected onto an
  // underlying born configuration.
  bool applyCut = nSteps > 0 && FullHistory.select(RN)->nClusterings() > 0;

  // Enfore merging scale cut if the event did not pass the merging scale
  // criterion.
  bool enforceCutOnLHE  = settingsPtr->flag("Merging:enforceCutOnLHE");
  if ( enforceCutOnLHE && applyCut && nSteps == nRequested
    && tmsnow < tmsval ) {
    string message="Warning in Merging::mergeProcessUNLOPS: Les Houches Event";
    message+=" fails merging scale cut. Reject event.";
    infoPtr->errorMsg(message);
    mergingHooksPtr->setWeightCKKWL(0.);
    mergingHooksPtr->setWeightFIRST(0.);
    return -1;
  }

  // Potentially recluster real emission jets for powheg input containing
  // "too many" jets, i.e. real-emission kinematics.
  bool containsRealKin = nSteps > nRequested && nSteps > 0;
  if ( containsRealKin ) nRecluster += nSteps - nRequested;

  // Remove real emission events without underlying Born configuration from
  // the loop sample, since such states will be taken care of by tree-level
  // samples.
  bool allowIncompleteReal =
    settingsPtr->flag("Merging:allowIncompleteHistoriesInReal");
  if ( doUNLOPSLoop && containsRealKin && !allowIncompleteReal
    && FullHistory.select(RN)->nClusterings() == 0 ) {
    mergingHooksPtr->setWeightCKKWL(0.);
    mergingHooksPtr->setWeightFIRST(0.);
    return -1;
  }

  // Discard if the state could not be reclustered to any state above TMS.
  int nPerformed = 0;
  if ( nSteps > 0 && !allowIncompleteReal
    && ( doUNLOPSSubt || doUNLOPSSubtNLO || containsRealKin )
    && !FullHistory.getFirstClusteredEventAboveTMS( RN, nRecluster,
          newProcess, nPerformed, false ) ) {
    mergingHooksPtr->setWeightCKKWL(0.);
    mergingHooksPtr->setWeightFIRST(0.);
    return -1;
  }
  // Check reclustering steps to correctly apply MPI.
  mergingHooksPtr->nMinMPI(nSteps - nPerformed);

  // Perform one reclustering for real emission kinematics, then apply
  // merging scale cut on underlying Born kinematics.
  if ( containsRealKin ) {
    Event dummy = Event();
    // Initialise temporary output of reclustering.
    dummy.clear();
    dummy.init( "(hard process-modified)", particleDataPtr );
    dummy.clear();
    // Recluster once.
    FullHistory.getClusteredEvent( RN, nSteps, dummy );
    double tnowNew  = mergingHooksPtr->tmsNow( dummy );
    // Veto if underlying Born kinematics do not pass merging scale cut.
    if ( enforceCutOnLHE && nSteps > 0 && nRequested > 0
      && tnowNew < tmsval ) {
      mergingHooksPtr->setWeightCKKWL(0.);
      mergingHooksPtr->setWeightFIRST(0.);
      return -1;
    }
  }

  // Calculate weights.
  // Do LO or first part of NLO tree-level reweighting
  if( doUNLOPSTree ) {
    // Perform reweighting with Sudakov factors, save as ratios and
    // PDF ratio weights
    wgt = FullHistory.weight_UNLOPS_TREE( trialPartonLevelPtr,
            mergingHooksPtr->AlphaS_FSR(), mergingHooksPtr->AlphaS_ISR(), RN);
  } else if( doUNLOPSLoop ) {
    // No reweighting, just set event scales properly
    wgt = FullHistory.weight_UNLOPS_LOOP( trialPartonLevelPtr, RN);
  } else if( doUNLOPSSubtNLO ) {
    // The standard prescripition contains no real-emission parts
    // No reweighting, just set event scales properly
    wgt = FullHistory.weight_UNLOPS_SUBTNLO( trialPartonLevelPtr, RN);
  } else if( doUNLOPSSubt ) {
    // The standard prescripition contains no subtraction parts
    // No reweighting, just set event scales properly
    wgt = FullHistory.weight_UNLOPS_SUBT( trialPartonLevelPtr,
            mergingHooksPtr->AlphaS_FSR(), mergingHooksPtr->AlphaS_ISR(), RN);
  }

  // Event with production scales set for further (trial) showering
  // and starting conditions for the shower.
  if (!doUNLOPSSubt && !doUNLOPSSubtNLO && !containsRealKin )
    FullHistory.getStartingConditions(RN, process);
  // Do reclustering (looping) steps.
  else FullHistory.getFirstClusteredEventAboveTMS( RN, nRecluster, process,
    nPerformed, true );

  // Allow to dampen histories in which the lowest multiplicity reclustered
  // state does not pass the lowest multiplicity cut of the matrix element
  double dampWeight = mergingHooksPtr->dampenIfFailCuts(
           FullHistory.lowestMultProc(RN) );
  // Save the weight of the event for histogramming. Only change the
  // event weight after trial shower on the matrix element
  // multiplicity event (= in doVetoStep)
  wgt *= dampWeight;

  // For tree-level or subtractive sammples, rescale with k-Factor
  if ( doUNLOPSTree || doUNLOPSSubt ){
    // Find k-factor
    double kFactor = 1.;
    if ( nSteps > mergingHooksPtr->nMaxJetsNLO() )
      kFactor = mergingHooksPtr->kFactor( mergingHooksPtr->nMaxJetsNLO() );
    else kFactor = mergingHooksPtr->kFactor(nSteps);
    // For NLO merging, rescale CKKW-L weight with k-factor
    wgt *= (nRecluster == 2 && nloTilde) ? 1. : kFactor;
  }

  // Save the weight of the event for histogramming
  mergingHooksPtr->setWeightCKKWL(wgt);

  // Check if we need to subtract the O(\alpha_s)-term. If the number
  // of additional partons is larger than the number of jets for
  // which loop matrix elements are available, do standard UMEPS.
  int nMaxNLO     = mergingHooksPtr->nMaxJetsNLO();
  bool doOASTree  = doUNLOPSTree && nSteps <= nMaxNLO;
  bool doOASSubt  = doUNLOPSSubt && nSteps <= nMaxNLO+1 && nSteps > 0;

  // Now begin NLO part for tree-level events
  if ( doOASTree || doOASSubt ) {

    // Decide on which order to expand to.
    int order = ( nSteps > 0 && nSteps <= nMaxNLO) ? 1 : -1;

    // Exclusive inputs:
    // Subtract only the O(\alpha_s^{n+0})-term from the tree-level
    // subtraction, if we're at the highest NLO multiplicity (nMaxNLO).
    if ( nloTilde && doUNLOPSSubt && nRecluster == 1
      && nSteps == nMaxNLO+1 ) order = 0;

    // Exclusive inputs:
    // Do not remove the O(as)-term if the number of reclusterings
    // exceeds the number of NLO jets, or if more clusterings have
    // been performed.
    if (nloTilde && doUNLOPSSubt && ( nSteps > nMaxNLO+1
      || (nSteps == nMaxNLO+1 && nPerformed != nRecluster) ))
        order = -1;

    // Calculate terms in expansion of the CKKW-L weight.
    wgtFIRST = FullHistory.weight_UNLOPS_CORRECTION( order,
      trialPartonLevelPtr, mergingHooksPtr->AlphaS_FSR(),
      mergingHooksPtr->AlphaS_ISR(), RN, rndmPtr );

    // Exclusive inputs:
    // Subtract the O(\alpha_s^{n+1})-term from the tree-level
    // subtraction, not the O(\alpha_s^{n+0})-terms.
    if ( nloTilde && doUNLOPSSubt && nRecluster == 1
      && nPerformed == nRecluster && nSteps <= nMaxNLO )
      wgtFIRST += 1.;

    // If necessary, also dampen the O(\alpha_s)-term
    wgtFIRST *= dampWeight;
    // Set the subtractive weight to the value calculated so far
    mergingHooksPtr->setWeightFIRST(wgtFIRST);
    // Subtract the O(\alpha_s)-term from the CKKW-L weight
    // If PDF contributions have not been included, subtract these later
    wgt = wgt - wgtFIRST;

  }

  // Set QCD 2->2 starting scale different from arbitrary scale in LHEF!
  // --> Set to minimal mT of partons.
  int nFinal = 0;
  double muf = process[0].e();
  for ( int i=0; i < process.size(); ++i )
  if ( process[i].isFinal()
    && (process[i].colType() != 0 || process[i].id() == 22 ) ) {
    nFinal++;
    muf = min( muf, abs(process[i].mT()) );
  }
  // For pure QCD dijet events (only!), set the process scale to the
  // transverse momentum of the outgoing partons.
  if ( nSteps == 0 && nFinal == 2
    && ( mergingHooksPtr->getProcessString().compare("pp>jj") == 0
      || mergingHooksPtr->getProcessString().compare("pp>aj") == 0) )
    process.scale(muf);

  // Reset hard process candidates (changed after clustering a parton).
  mergingHooksPtr->storeHardProcessCandidates( process );

  // Check if resonance structure has been changed
  //  (e.g. because of clustering W/Z/gluino)
  vector <int> oldResonance;
  for ( int i=0; i < newProcess.size(); ++i )
    if ( newProcess[i].status() == 22 )
      oldResonance.push_back(newProcess[i].id());
  vector <int> newResonance;
  for ( int i=0; i < process.size(); ++i )
    if ( process[i].status() == 22 )
      newResonance.push_back(process[i].id());
  // Compare old and new resonances
  for ( int i=0; i < int(oldResonance.size()); ++i )
    for ( int j=0; j < int(newResonance.size()); ++j )
      if ( newResonance[j] == oldResonance[i] ) {
        oldResonance[i] = 99;
        break;
      }
  bool hasNewResonances = (newResonance.size() != oldResonance.size());
  for ( int i=0; i < int(oldResonance.size()); ++i )
    hasNewResonances = (oldResonance[i] != 99);

  // If necessary, reattach resonance decay products.
  if (!hasNewResonances) mergingHooksPtr->reattachResonanceDecays(process);

  // Allow merging hooks to remove emissions from now on.
  mergingHooksPtr->doIgnoreEmissions(false);

  // If no-emission probability is zero.
  if ( wgt == 0. ) return 0;

  // If the resonance structure of the process has changed due to reclustering,
  // redo the resonance decays in Pythia::next()
  if (hasNewResonances) return 2;

  // Done
  return 1;

}

//--------------------------------------------------------------------------

// Function to apply the merging scale cut on an input event.

bool Merging::cutOnProcess( Event& process) {

  // Save number of looping steps
  mergingHooksPtr->nReclusterSave = settingsPtr->mode("Merging:nRecluster");
  int nRequested        = settingsPtr->mode("Merging:nRequested");

  // For now, prefer construction of ordered histories.
  mergingHooksPtr->orderHistories(true);
  // For pp > h, allow cut on state, so that underlying processes
  // can be clustered to gg > h
  if ( mergingHooksPtr->getProcessString().compare("pp>h") == 0)
    mergingHooksPtr->allowCutOnRecState(true);

  // Prepare process record for merging. If Pythia has already decayed
  // resonances used to define the hard process, remove resonance decay
  // products.
  Event newProcess( mergingHooksPtr->bareEvent( process, true) );
  // Store candidates for the splitting V -> qqbar'
  mergingHooksPtr->storeHardProcessCandidates( newProcess );

  // Check if event passes the merging scale cut.
  double tmsval  = mergingHooksPtr->tms();
  // Get merging scale in current event.
  double tmsnow  = mergingHooksPtr->tmsNow( newProcess );
  // Calculate number of clustering steps
  int nSteps     = mergingHooksPtr->getNumberOfClusteringSteps( newProcess);

  // Too few steps can be possible if a chain of resonance decays has been
  // removed. In this case, reject this event, since it will be handled in
  // lower-multiplicity samples.
  if (nSteps < nRequested) return true;

  // Reset the minimal tms value, if necessary.
  tmsNowMin = (nSteps == 0) ? 0. : min(tmsNowMin, tmsnow);

  // Potentially recluster real emission jets for powheg input containing
  // "too many" jets, i.e. real-emission kinematics.
  bool containsRealKin = nSteps > nRequested && nSteps > 0;

  // Get random number to choose a path.
  double RN = rndmPtr->flat();
  // Set dummy process scale.
  newProcess.scale(0.0);
  // Generate all histories
  History FullHistory( nSteps, 0.0, newProcess, Clustering(), mergingHooksPtr,
            (*beamAPtr), (*beamBPtr), particleDataPtr, infoPtr, true, true,
            true, true, 1.0, 0);
  // Project histories onto desired branches, e.g. only ordered paths.
  FullHistory.projectOntoDesiredHistories();

  // Remove real emission events without underlying Born configuration from
  // the loop sample, since such states will be taken care of by tree-level
  // samples.
  bool allowIncompleteReal =
    settingsPtr->flag("Merging:allowIncompleteHistoriesInReal");
  if ( containsRealKin && !allowIncompleteReal
    && FullHistory.select(RN)->nClusterings() == 0 )
    return true;

  // Cut if no history passes the cut on the lowest-multiplicity state.
  double dampWeight = mergingHooksPtr->dampenIfFailCuts(
           FullHistory.lowestMultProc(RN) );
  if ( dampWeight == 0. ) return true;

  // Do not apply cut if the configuration could not be projected onto an
  // underlying born configuration.
  if ( nSteps > 0 && FullHistory.select(RN)->nClusterings() == 0 )
    return false;

  // Now enfore merging scale cut if the event did not pass the merging scale
  // criterion.
  if ( nSteps > 0 && nSteps == nRequested && tmsnow < tmsval ) {
    string message="Warning in Merging::cutOnProcess: Les Houches Event";
    message+=" fails merging scale cut. Reject event.";
    infoPtr->errorMsg(message);
    return true;
  }

  if ( FullHistory.select(RN)->nClusterings() < nSteps) {
    string message="Warning in Merging::cutOnProcess: No clusterings";
    message+=" found. History incomplete.";
    infoPtr->errorMsg(message);
  }

  // Done if no real-emission jets are present.
  if ( !containsRealKin ) return false;

  // Now cut on events that contain an additional real-emission jet.
  // Perform one reclustering for real emission kinematics, then apply merging
  // scale cut on underlying Born kinematics.
  if ( containsRealKin ) {
    Event dummy = Event();
    // Initialise temporary output of reclustering.
    dummy.clear();
    dummy.init( "(hard process-modified)", particleDataPtr );
    dummy.clear();
    // Recluster once.
    FullHistory.getClusteredEvent( RN, nSteps, dummy );
    double tnowNew  = mergingHooksPtr->tmsNow( dummy );
    // Veto if underlying Born kinematics do not pass merging scale cut.
    if ( nSteps > 0 && nRequested > 0 && tnowNew < tmsval ) return true;
  }

  // Done if only interested in cross section estimate after cuts.
  return false;

}

//==========================================================================

} // end namespace Pythia8
