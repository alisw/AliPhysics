// MergingHooks.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel.
// Function definitions (not found in the header) for the Merging class.

#include "Merging.h"

namespace Pythia8 {
 
//==========================================================================

// The Merging class.

//--------------------------------------------------------------------------

// Number of trial emission to use for calculating the average number of 
// emissions
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

  // Header.
  os << "\n *-------  PYTHIA Matrix Element Merging Information  ------"
     << "-------------------------------------------------------*\n"
     << " |                                                            "
     << "                                                     |\n";

  // Recall switch to enfore merging scale cut.
  bool enforceCutOnLHE  = settingsPtr->flag("Merging:enforceCutOnLHE");
  // Recall merging scale value.
  double tmsval = mergingHooksPtr->tms();

  // Print warning if the minimal tms value of any event was significantly
  // above the desired merging scale value.
  if ( enforceCutOnLHE && tmsNowMin > TMSMISMATCH*tmsval )
    os << " | Warning in Merging::statistics: All Les Houches events"
       << " significantly above Merging:TMS cut. Please check.\n";
  else
    os << " |      0   no warnings to report              \n";

  // Listing finished.
  os << " |                                                            "
     << "                                                     |\n"
     << " *-------  End PYTHIA Matrix Element Merging Information -----"
     << "-----------------------------------------------------*" << endl;

  // Reset minimal tms value. 
  tmsNowMin             = infoPtr->eCM();

} 

//--------------------------------------------------------------------------

// Function to steer different merging prescriptions.

int Merging::mergeProcess(Event& process){

  int vetoCode = 1;

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
  int nSteps   = mergingHooksPtr->getNumberOfClusteringSteps( newProcess);

  // Reset the minimal tms value, if necessary.
  tmsNowMin = (nSteps == 0) ? 0. : min(tmsNowMin, tmsnow);

  // Enfore merging scale cut if the event did not pass the merging scale
  // criterion.
  bool enforceCutOnLHE  = settingsPtr->flag("Merging:enforceCutOnLHE");
  if ( enforceCutOnLHE && nSteps > 0 && tmsnow < tmsval ) {
    string message="Warning in Merging::mergeProcessCKKWL: Les Houches Event";
    message+=" fails merging scale cut. Cut by rejecting event.";
    infoPtr->errorMsg(message);
    mergingHooksPtr->setWeightCKKWL(0.);
    return -1;
  }

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
  int nSteps   = mergingHooksPtr->getNumberOfClusteringSteps( newProcess );

  // Reset the minimal tms value, if necessary.
  tmsNowMin = (nSteps == 0) ? 0. : min(tmsNowMin, tmsnow);

  // Enfore merging scale cut if the event did not pass the merging scale
  // criterion.
  bool enforceCutOnLHE  = settingsPtr->flag("Merging:enforceCutOnLHE");
  if ( enforceCutOnLHE && nSteps > 0 && tmsnow < tmsval ) {
    string message="Warning in Pythia::mergeProcessUMEPS: Les Houches Event";
    message+=" fails merging scale cut. Cut by rejecting event.";
    infoPtr->errorMsg(message);
    mergingHooksPtr->setWeightCKKWL(0.);
    return -1;
  }

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

  // Discard incomplete histories.
  if ( nSteps > 0 && doUMEPSSubt
     && !FullHistory.foundCompleteHistories() ){
    mergingHooksPtr->setWeightCKKWL(0.);
    return -1;
  }

  // Check reclustering steps to correctly apply MPI.
  int nPerformed = 0;
  if ( doUMEPSSubt ) FullHistory.getFirstClusteredEventAboveTMS( RN, 
    nRecluster, newProcess, nPerformed, false );
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


  // If necessary, reattach resonance decay products.
  mergingHooksPtr->reattachResonanceDecays(process); 

  if ( doUMEPSSubt && nStepsNew > 0 )
    mergingHooksPtr->muMI( process.scale() );
  else if ( doUMEPSSubt ) mergingHooksPtr->muMI( infoPtr->eCM() );

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

  // Reset the minimal tms value, if necessary.
  tmsNowMin = (nSteps == 0) ? 0. : min(tmsNowMin, tmsnow);

  // Enfore merging scale cut if the event did not pass the merging scale
  // criterion.
  bool enforceCutOnLHE  = settingsPtr->flag("Merging:enforceCutOnLHE");
  if ( enforceCutOnLHE && nSteps > 0 && nSteps == nRequested
    && tmsnow < tmsval ) {
    string message="Warning in Pythia::mergeProcessNL3: Les Houches Event";
    message+=" fails merging scale cut. Cut by rejecting event.";
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

  // Discard incomplete histories when reclustering.
  if ( nSteps > 0 && doNL3Subt
    && !FullHistory.foundCompleteHistories() ){
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
      return true;
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

  // Set shower starting scale for non-highest multiplicities to merging
  // scale.
  int nStepsFin = mergingHooksPtr->getNumberOfClusteringSteps( process);

  // If necessary, reattach resonance decay products.
  mergingHooksPtr->reattachResonanceDecays(process); 

  if ( ( doNL3Subt || containsRealKin) && nStepsFin > 0 )
    mergingHooksPtr->muMI( process.scale() );
  else if ( doNL3Subt || containsRealKin )
    mergingHooksPtr->muMI( infoPtr->eCM() );

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

  // Reset the minimal tms value, if necessary.
  tmsNowMin = (nSteps == 0) ? 0. : min(tmsNowMin, tmsnow);

  // Enfore merging scale cut if the event did not pass the merging scale
  // criterion.
  bool enforceCutOnLHE  = settingsPtr->flag("Merging:enforceCutOnLHE");
  if ( enforceCutOnLHE && nSteps > 0 && nSteps == nRequested
    && tmsnow < tmsval ) {
    string message="Warning in Pythia::mergeProcessUNLOPS: Les Houches Event";
    message+=" fails merging scale cut. Cut by rejecting event.";
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

  // Discard incomplete histories when reclustering.
  if (  nSteps > 0
    && ( doUNLOPSSubt || doUNLOPSSubtNLO || doUNLOPSLoop )
    && !FullHistory.foundCompleteHistories() ){
    mergingHooksPtr->setWeightCKKWL(0.);
    mergingHooksPtr->setWeightFIRST(0.);
    return -1;
  }

  // Potentially recluster real emission jets for powheg input containing
  // "too many" jets, i.e. real-emission kinematics.
  bool containsRealKin = nSteps > nRequested && nSteps > 0;
  if ( containsRealKin ) nRecluster += nSteps - nRequested;

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
    if ( enforceCutOnLHE && nSteps > 0 && nRequested > 0
      && tnowNew < tmsval ) {
      mergingHooksPtr->setWeightCKKWL(0.);
      mergingHooksPtr->setWeightFIRST(0.);
      return -1;
    }
  }

  // Check reclustering steps to correctly apply MPI.
  int nPerformed = 0;
  if (  doUNLOPSSubt || doUNLOPSSubtNLO || containsRealKin )
    FullHistory.getFirstClusteredEventAboveTMS( RN, nRecluster, newProcess,
      nPerformed, false );
  mergingHooksPtr->nMinMPI(nSteps - nPerformed);

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
  bool doOASSubt  = doUNLOPSSubt && nSteps <= nMaxNLO+1 && nSteps > 0
                 && nMaxNLO != 0;

  // Now begin NLO part for tree-level events
  if ( doOASTree || doOASSubt ) {
    // Calculate terms in expansion of the CKKW-L weight.
    int order = ( nSteps == 1 ) ? 1 : -1;
    if ( nSteps == 2 && nRecluster == 1 && nloTilde ) order = 0;
    wgtFIRST = FullHistory.weight_UNLOPS_CORRECTION( order, 
      trialPartonLevelPtr, mergingHooksPtr->AlphaS_FSR(),
      mergingHooksPtr->AlphaS_ISR(), RN, rndmPtr );
    // Second reclustering term
    if ( nSteps == 1 && doUNLOPSSubt && nloTilde ) wgtFIRST += 1.;

    // If necessary, also dampen the O(\alpha_s)-term
    wgtFIRST *= dampWeight;
    // Set the subtractive weight to the value calculated so far
    mergingHooksPtr->setWeightFIRST(wgtFIRST);
    // Subtract the O(\alpha_s)-term from the CKKW-L weight
    // If PDF contributions have not been included, subtract these later
    wgt = wgt - wgtFIRST;
  }

  // Do not subtract the O(as)-term if more than one reclustering has been
  // performed.
  if ( nPerformed > nRecluster ) mergingHooksPtr->setWeightFIRST(0.);

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

  // Set shower starting scale for non-highest multiplicities to merging
  // scale.
  int nStepsFin = mergingHooksPtr->getNumberOfClusteringSteps( process);

  // If necessary, reattach resonance decay products.
  mergingHooksPtr->reattachResonanceDecays(process); 

  if ( ( doUNLOPSSubt || doUNLOPSSubtNLO || containsRealKin)
    && nStepsFin > 0 )
    mergingHooksPtr->muMI( process.scale() );
  else if ( doUNLOPSSubt || doUNLOPSSubtNLO || containsRealKin )
    mergingHooksPtr->muMI( infoPtr->eCM() );

  // Allow merging hooks to remove emissions from now on.
  mergingHooksPtr->doIgnoreEmissions(false);

  // If no-emission probability is zero.
  if ( wgt == 0. ) return 0;

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

  // Reset the minimal tms value, if necessary.
  tmsNowMin = (nSteps == 0) ? 0. : min(tmsNowMin, tmsnow);

  // Now enfore merging scale cut if the event did not pass the merging scale
  // criterion.
  if ( nSteps > 0 && nSteps == nRequested && tmsnow < tmsval ) {
    string message="Warning in Merging::cutOnProcess: Les Houches Event";
    message+=" fails merging scale cut. Cut by rejecting event.";
    infoPtr->errorMsg(message);
    return true;
  }

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

  // Discard incomplete histories when reclustering.
  if ( nSteps > 0 && !FullHistory.foundCompleteHistories() ) return true;

  // Cut if no history passes the cut on the lowest-multiplicity state.
  double dampWeight = mergingHooksPtr->dampenIfFailCuts( 
           FullHistory.lowestMultProc(RN) );
  if ( dampWeight == 0. ) return true;

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
