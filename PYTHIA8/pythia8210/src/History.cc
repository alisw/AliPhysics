// History.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel.
// Function definitions (not found in the header) for the
// Clustering and History classes.

#include "Pythia8/History.h"

namespace Pythia8 {

//==========================================================================

// The Clustering class.

//--------------------------------------------------------------------------

// Declaration of Clustering class
// This class holds information about one radiator, recoiler,
// emitted system.
// This class is a container class for History class use.

// print for debug
void Clustering::list() const {
  cout << " emt " << emitted
       << " rad " << emittor
       << " rec " << recoiler
       << " partner " << partner
       << " pTscale " << pTscale << endl;
}

//==========================================================================

// The History class.

// A History object represents an event in a given step in the CKKW-L
// clustering procedure. It defines a tree-like recursive structure,
// where the root node represents the state with n jets as given by
// the matrix element generator, and is characterized by the member
// variable mother being null. The leaves on the tree corresponds to a
// fully clustered paths where the original n-jets has been clustered
// down to the Born-level state. Also states which cannot be clustered
// down to the Born-level are possible - these will be called
// incomplete. The leaves are characterized by the vector of children
// being empty.

//--------------------------------------------------------------------------

// Number of trial emission to use for calculating the average number of
// emissions
const int History::NTRIAL = 1;

//--------------------------------------------------------------------------

// Declaration of History class
// The only constructor. Default arguments are used when creating
// the initial history node. The \a depth is the maximum number of
// clusterings requested. \a scalein is the scale at which the \a
// statein was clustered (should be set to the merging scale for the
// initial history node. \a beamAIn and beamBIn are needed to
// calcutate PDF ratios, \a particleDataIn to have access to the
// correct masses of particles. If \a isOrdered is true, the previous
// clusterings has been ordered. \a is the PDF ratio for this
// clustering (=1 for FSR clusterings). \a probin is the accumulated
// probabilities for the previous clusterings, and \ mothin is the
// previous history node (null for the initial node).

History::History( int depth,
         double scalein,
         Event statein,
         Clustering c,
         MergingHooks* mergingHooksPtrIn,
         BeamParticle beamAIn,
         BeamParticle beamBIn,
         ParticleData* particleDataPtrIn,
         Info* infoPtrIn,
         bool isOrdered = true,
         bool isStronglyOrdered = true,
         bool isAllowed = true,
         bool isNextInInput = true,
         double probin = 1.0,
         History * mothin = 0)
    : state(statein),
      mother(mothin),
      sumpath(0.0),
      sumGoodBranches(0.0),
      sumBadBranches(0.0),
      foundOrderedPath(false),
      foundStronglyOrderedPath(false),
      foundAllowedPath(false),
      foundCompletePath(false),
      scale(scalein),
      nextInInput(isNextInInput),
      prob(probin),
      clusterIn(c),
      iReclusteredOld(0),
      doInclude(true),
      mergingHooksPtr(mergingHooksPtrIn),
      beamA(beamAIn),
      beamB(beamBIn),
      particleDataPtr(particleDataPtrIn),
      infoPtr(infoPtrIn)
    {

  // Initialise beam particles
  setupBeams();

  // Update probability with PDF ratio
  if (mother && mergingHooksPtr->includeRedundant()) prob *= pdfForSudakov();

  // Minimal scalar sum of pT used in Herwig to choose history
  // Keep track of scalar PT
  if (mother) {
    double acoll = (mother->state[clusterIn.emittor].isFinal())
                   ? mergingHooksPtr->herwigAcollFSR()
                   : mergingHooksPtr->herwigAcollISR();
    sumScalarPT = mother->sumScalarPT + acoll*scale;
  } else
    sumScalarPT = 0.0;

  // Remember reclustered radiator in lower multiplicity state
  if ( mother ) iReclusteredOld = mother->iReclusteredNew;

  // Check if more steps should be taken.
  int nFinalP = 0;
  int nFinalW = 0;
  for ( int i = 0; i < int(state.size()); ++i )
    if ( state[i].isFinal() ) {
      if ( state[i].colType() != 0 )
        nFinalP++;
      if ( state[i].idAbs() == 24 )
        nFinalW++;
    }
  if ( mergingHooksPtr->doWClustering()
    && nFinalP == 2 && nFinalW == 0 ) depth = 0;

  // If this is not the fully clustered state, try to find possible
  // QCD clusterings.
  vector<Clustering> clusterings;
  if ( depth > 0 ) clusterings = getAllQCDClusterings();

  // If necessary, try to find possible EW clusterings.
  vector<Clustering> clusteringsEW;
  if ( depth > 0 && mergingHooksPtr->doWClustering() )
    clusteringsEW = getAllEWClusterings();
  if ( !clusteringsEW.empty() ) {
    clusterings.insert( clusterings.end(), clusteringsEW.begin(),
                        clusteringsEW.end() );
  }

  // If necessary, try to find possible SQCD clusterings.
  vector<Clustering> clusteringsSQCD;
  if ( depth > 0 && mergingHooksPtr->doSQCDClustering() )
    clusteringsSQCD = getAllSQCDClusterings();
  if ( !clusteringsSQCD.empty() )
    clusterings.insert( clusterings.end(), clusteringsSQCD.begin(),
                        clusteringsSQCD.end() );

  // If no clusterings were found, the recursion is done and we
  // register this node.
  if ( clusterings.empty() ) {
    // Multiply with hard process matrix element.
    prob *= hardProcessME(state);
    registerPath( *this, isOrdered, isStronglyOrdered, isAllowed, depth == 0 );
    return;
  }

  // Now we sort the possible clusterings so that we try the
  // smallest scale first.
  multimap<double, Clustering *> sorted;
  for ( int i = 0, N = clusterings.size(); i < N; ++i ) {
    sorted.insert(make_pair(clusterings[i].pT(), &clusterings[i]));
  }

  for ( multimap<double, Clustering *>::iterator it = sorted.begin();
  it != sorted.end(); ++it ) {

    // If this path is not strongly ordered and we already have found an
    // ordered path, then we don't need to continue along this path.
    bool stronglyOrdered = isStronglyOrdered;
    if ( mergingHooksPtr->enforceStrongOrdering()
      && ( !stronglyOrdered
         || ( mother && ( it->first <
                mergingHooksPtr->scaleSeparationFactor()*scale ) ))) {
      if ( onlyStronglyOrderedPaths()  ) continue;
      stronglyOrdered = false;
    }

    // Check if reclustering follows ordered sequence.
    bool ordered = isOrdered;
    if (  mergingHooksPtr->orderInRapidity()
      && mergingHooksPtr->orderHistories() ) {
      // Get new z value
      double z = getCurrentZ((*it->second).emittor,
                   (*it->second).recoiler,(*it->second).emitted);
      // Get z value of splitting that produced this state
      double zOld = (!mother) ? 0. : mother->getCurrentZ(clusterIn.emittor,
                       clusterIn.recoiler,clusterIn.emitted);
      // If this path is not ordered in pT and y, and we already have found
      // an ordered path, then we don't need to continue along this path.
      if ( !ordered || ( mother && (it->first < scale
         || it->first < pow(1. - z,2) / (z * (1. - zOld ))*scale ))) {
        if ( onlyOrderedPaths()  ) continue;
        ordered = false;
      }

    } else if ( mergingHooksPtr->orderHistories() ) {
      // If this path is not ordered in pT and we already have found an
      // ordered path, then we don't need to continue along this path, unless
      // we have not yet found an allowed path.
      if ( !ordered || ( mother && (it->first < scale) ) ) {
        if ( onlyOrderedPaths() && onlyAllowedPaths() ) continue;
        ordered = false;
      }
    }

    // Check if reclustered state should be disallowed.
    bool doCut = mergingHooksPtr->canCutOnRecState()
              || mergingHooksPtr->allowCutOnRecState();
    bool allowed = isAllowed;
    if (  doCut
      && mergingHooksPtr->doCutOnRecState(cluster(*it->second)) ) {
      if ( onlyAllowedPaths()  ) continue;
      allowed = false;
    }

    // Perform the clustering and recurse and construct the next
    // history node.
    children.push_back(new History(depth - 1,it->first,cluster(*it->second),
           *it->second, mergingHooksPtr, beamA, beamB, particleDataPtr,
           infoPtr, ordered, stronglyOrdered, allowed, true,
           prob*getProb(*it->second), this ));
  }
}

//--------------------------------------------------------------------------

// Function to project all possible paths onto only the desired paths.

bool History::projectOntoDesiredHistories() {
  // At the moment, only trim histories.
  return trimHistories();
}

//--------------------------------------------------------------------------

// In the initial history node, select one of the paths according to
// the probabilities. This function should be called for the initial
// history node.
// IN  trialShower*    : Previously initialised trialShower object,
//                       to perform trial showering and as
//                       repository of pointers to initialise alphaS
//     PartonSystems* : PartonSystems object needed to initialise
//                      shower objects
// OUT double         : (Sukadov) , (alpha_S ratios) , (PDF ratios)

double History::weightTREE(PartonLevel* trial, AlphaStrong * asFSR,
                  AlphaStrong * asISR, double RN) {

  if ( mergingHooksPtr->canCutOnRecState() && !foundAllowedPath ) {
    string message="Warning in History::weightTREE: No allowed history";
    message+=" found. Using disallowed history.";
    infoPtr->errorMsg(message);
  }
  if ( mergingHooksPtr->orderHistories() && !foundOrderedPath ) {
    string message="Warning in History::weightTREE: No ordered history";
    message+=" found. Using unordered history.";
    infoPtr->errorMsg(message);
  }
  if ( mergingHooksPtr->canCutOnRecState()
    && mergingHooksPtr->orderHistories()
    && !foundAllowedPath && !foundOrderedPath ) {
    string message="Warning in History::weightTREE: No allowed or ordered";
    message+=" history found.";
    infoPtr->errorMsg(message);
  }

  // Read alpha_S in ME calculation and maximal scale (eCM)
  double asME     = infoPtr->alphaS();
  double maxScale = (foundCompletePath) ? infoPtr->eCM()
                  : mergingHooksPtr->muFinME();
  // Select a path of clusterings
  History *  selected = select(RN);
  // Set scales in the states to the scales pythia would have set
  selected->setScalesInHistory();

  // Get weight.
  double sudakov   = 1.;
  double asWeight  = 1.;
  double pdfWeight = 1.;

  // Do trial shower, calculation of alpha_S ratios, PDF ratios
  sudakov  = selected->weightTree( trial, asME, maxScale,
               selected->clusterIn.pT(), asFSR, asISR, asWeight,
               pdfWeight );

  // MPI no-emission probability
  int njetsMaxMPI = mergingHooksPtr->nMinMPI();
  double mpiwt = selected->weightTreeEmissions( trial, -1, 0, njetsMaxMPI,
                   maxScale );

  // Set hard process renormalisation scale to default Pythia value.
  bool resetScales = mergingHooksPtr->resetHardQRen();
  // For pure QCD dijet events, evaluate the coupling of the hard process at
  // a more reasonable pT, rather than evaluation \alpha_s at a fixed
  // arbitrary scale.
  if ( resetScales
    && mergingHooksPtr->getProcessString().compare("pp>jj") == 0) {
    // Reset to a running coupling. Here we choose FSR for simplicity.
    double newQ2Ren = pow2( selected->hardRenScale(selected->state) );
    double runningCoupling = (*asFSR).alphaS(newQ2Ren) / asME;
    asWeight *= pow2(runningCoupling);
  }

  // For prompt photon events, evaluate the coupling of the hard process at
  // a more reasonable pT, rather than evaluation \alpha_s at a fixed
  // arbitrary scale.
  if ( resetScales
    && mergingHooksPtr->getProcessString().compare("pp>aj") == 0) {
    // Reset to a running coupling. In prompt photon always ISR.
    double newQ2Ren = pow2( selected->hardRenScale(selected->state) );
    double runningCoupling =
      (*asISR).alphaS( newQ2Ren + pow(mergingHooksPtr->pT0ISR(),2) ) / asME;
    asWeight *= runningCoupling;
  }

  // Done
  return (sudakov*asWeight*pdfWeight*mpiwt);

}

//--------------------------------------------------------------------------

// Function to return weight of virtual correction and subtractive events
// for NL3 merging

double History::weightLOOP(PartonLevel* trial, double RN ) {

  if ( mergingHooksPtr->canCutOnRecState() && !foundAllowedPath ) {
    string message="Warning in History::weightLOOP: No allowed history";
    message+=" found. Using disallowed history.";
    infoPtr->errorMsg(message);
  }

  // Select a path of clusterings
  History *  selected = select(RN);
  // Set scales in the states to the scales pythia would have set
  selected->setScalesInHistory();

  // So far, no reweighting
  double wt = 1.;

  // Only reweighting with MPI no-emission probability
  double maxScale = (foundCompletePath) ? infoPtr->eCM()
                  : mergingHooksPtr->muFinME();
  int njetsMaxMPI = mergingHooksPtr->nMinMPI();
  double mpiwt = selected->weightTreeEmissions( trial, -1, 0, njetsMaxMPI,
                   maxScale );
  wt = mpiwt;
  // Done
  return wt;
}

//--------------------------------------------------------------------------

// Function to calculate O(\alpha_s)-term of CKKWL-weight for NLO merging

double History::weightFIRST(PartonLevel* trial, AlphaStrong* asFSR,
                  AlphaStrong* asISR, double RN, Rndm* rndmPtr ) {

  // Read alpha_S in ME calculation and maximal scale (eCM)
  double asME   = infoPtr->alphaS();
  double muR      = mergingHooksPtr->muRinME();
  double maxScale = (foundCompletePath)
                  ? infoPtr->eCM()
                  : mergingHooksPtr->muFinME();

  // Pick path of clusterings
  History *  selected = select(RN);
  // Set scales in the states to the scales pythia would have set
  selected->setScalesInHistory();

  double nSteps = mergingHooksPtr->getNumberOfClusteringSteps(state);

  // Get the lowest order k-factor and add first two terms in expansion
  double kFactor = asME * mergingHooksPtr->k1Factor(nSteps);

  // If using Bbar, which includes a tree-level part, subtract an
  // additional one, i.e. the O(\as^0) contribution as well
  double wt = 1. + kFactor;

  // Calculate sum of O(alpha) terms
  wt += selected->weightFirst(trial,asME, muR, maxScale, asFSR, asISR,
          rndmPtr );

  // Get starting scale for trial showers.
  double startingScale = (selected->mother) ? state.scale() : infoPtr->eCM();

  // Count emissions: New variant
  // Generate true average, not only one-point
  bool fixpdf = true;
  bool fixas  = true;
  double nWeight1 = 0.;
  for(int i=0; i < NTRIAL; ++i) {
    // Get number of emissions
    vector<double> unresolvedEmissionTerm = countEmissions( trial,
      startingScale, mergingHooksPtr->tms(), 2, asME, asFSR, asISR, 3,
      fixpdf, fixas );
    nWeight1 += unresolvedEmissionTerm[1];
  }

  wt += nWeight1/double(NTRIAL);

  // Done
  return wt;

}

//--------------------------------------------------------------------------

double History::weight_UMEPS_TREE(PartonLevel* trial, AlphaStrong * asFSR,
                  AlphaStrong * asISR, double RN) {
  // No difference to CKKW-L. Recycle CKKW-L function.
  return weightTREE( trial, asFSR, asISR, RN);
}

//--------------------------------------------------------------------------

// Function to return weight of virtual correction events for NLO merging

double History::weight_UMEPS_SUBT(PartonLevel* trial, AlphaStrong * asFSR,
                  AlphaStrong * asISR, double RN ) {

  // Read alpha_S in ME calculation and maximal scale (eCM)
  double asME     = infoPtr->alphaS();
  double maxScale = (foundCompletePath) ? infoPtr->eCM()
                  : mergingHooksPtr->muFinME();
  // Select a path of clusterings
  History *  selected = select(RN);
  // Set scales in the states to the scales pythia would have set
  selected->setScalesInHistory();

  // Get weight.
  double sudakov   = 1.;
  double asWeight  = 1.;
  double pdfWeight = 1.;

  // Do trial shower, calculation of alpha_S ratios, PDF ratios
  sudakov   = selected->weightTree(trial, asME, maxScale,
                selected->clusterIn.pT(), asFSR,asISR, asWeight, pdfWeight);

  // MPI no-emission probability.
  int njetsMaxMPI = mergingHooksPtr->nMinMPI()+1;
  double mpiwt = selected->weightTreeEmissions( trial, -1, 0, njetsMaxMPI,
                   maxScale );

  // Set hard process renormalisation scale to default Pythia value.
  bool resetScales = mergingHooksPtr->resetHardQRen();
  // For pure QCD dijet events, evaluate the coupling of the hard process at
  // a more reasonable pT, rather than evaluation \alpha_s at a fixed
  // arbitrary scale.
  if ( resetScales
    && mergingHooksPtr->getProcessString().compare("pp>jj") == 0) {
    // Reset to a running coupling. Here we choose FSR for simplicity.
    double newQ2Ren = pow2( selected->hardRenScale(selected->state) );
    double runningCoupling = (*asFSR).alphaS(newQ2Ren) / asME;
    asWeight *= pow(runningCoupling,2);
  }

  // For prompt photon events, evaluate the coupling of the hard process at
  // a more reasonable pT, rather than evaluation \alpha_s at a fixed
  // arbitrary scale.
  if ( resetScales
    && mergingHooksPtr->getProcessString().compare("pp>aj") == 0) {
    // Reset to a running coupling. In prompt photon always ISR.
    double newQ2Ren = pow2( selected->hardRenScale(selected->state) );
    double runningCoupling =
      (*asISR).alphaS( newQ2Ren + pow(mergingHooksPtr->pT0ISR(),2) ) / asME;
    asWeight *= runningCoupling;
  }

  // Done
  return (sudakov*asWeight*pdfWeight*mpiwt);

}

//--------------------------------------------------------------------------

double History::weight_UNLOPS_TREE(PartonLevel* trial, AlphaStrong * asFSR,
                  AlphaStrong * asISR, double RN) {

  // Read alpha_S in ME calculation and maximal scale (eCM)
  double asME     = infoPtr->alphaS();
  double maxScale = (foundCompletePath) ? infoPtr->eCM()
                  : mergingHooksPtr->muFinME();
  // Select a path of clusterings
  History *  selected = select(RN);
  // Set scales in the states to the scales pythia would have set
  selected->setScalesInHistory();

  // Get weight.
  double asWeight  = 1.;
  double pdfWeight = 1.;

  // Do trial shower, calculation of alpha_S ratios, PDF ratios
  double wt = selected->weightTree(trial, asME, maxScale,
                selected->clusterIn.pT(), asFSR, asISR, asWeight, pdfWeight);
  // MPI no-emission probability.
  int njetsMaxMPI = mergingHooksPtr->nMinMPI();
  double mpiwt = selected->weightTreeEmissions( trial, -1, 0, njetsMaxMPI,
                   maxScale );

  // Set hard process renormalisation scale to default Pythia value.
  bool resetScales = mergingHooksPtr->resetHardQRen();
  // For pure QCD dijet events, evaluate the coupling of the hard process at
  // a more reasonable pT, rather than evaluation \alpha_s at a fixed
  // arbitrary scale.
  if ( resetScales
    && mergingHooksPtr->getProcessString().compare("pp>jj") == 0) {
    // Reset to a running coupling. Here we choose FSR for simplicity.
    double newQ2Ren = pow2( selected->hardRenScale(selected->state) );
    double runningCoupling = (*asFSR).alphaS(newQ2Ren) / asME;
    asWeight *= pow(runningCoupling,2);
  }

  // For prompt photon events, evaluate the coupling of the hard process at
  // a more reasonable pT, rather than evaluation \alpha_s at a fixed
  // arbitrary scale.
  if ( resetScales
    && mergingHooksPtr->getProcessString().compare("pp>aj") == 0) {
    // Reset to a running coupling. In prompt photon always ISR.
    double newQ2Ren = pow2( selected->hardRenScale(selected->state) );
    double runningCoupling =
      (*asISR).alphaS( newQ2Ren + pow(mergingHooksPtr->pT0ISR(),2) ) / asME;
    asWeight *= runningCoupling;
  }

  // Done
  return (wt*asWeight*pdfWeight*mpiwt);

}

//--------------------------------------------------------------------------

double History::weight_UNLOPS_LOOP(PartonLevel* trial, double RN ) {
  // No difference to default NL3
  return weightLOOP(trial, RN );
}

//--------------------------------------------------------------------------

double History::weight_UNLOPS_SUBT(PartonLevel* trial, AlphaStrong * asFSR,
                  AlphaStrong * asISR, double RN ) {

  // Select a path of clusterings
  History *  selected = select(RN);
  // Set scales in the states to the scales pythia would have set
  selected->setScalesInHistory();
  // So far, no reweighting
  double wt = 1.;

  // Read alpha_S in ME calculation and maximal scale (eCM)
  double asME     = infoPtr->alphaS();
  double maxScale = (foundCompletePath)
                  ? infoPtr->eCM()
                  : mergingHooksPtr->muFinME();

  // Only allow two clusterings if all intermediate states above the
  // merging scale.
  double nSteps = mergingHooksPtr->getNumberOfClusteringSteps(state);
  if ( nSteps == 2 && mergingHooksPtr->nRecluster() == 2
    && ( !foundCompletePath
      || !selected->allIntermediateAboveRhoMS( mergingHooksPtr->tms() )) ) {
    return 0.;
  }

  // Get weights: alpha_S ratios and PDF ratios
  double asWeight  = 1.;
  double pdfWeight = 1.;
  // Do trial shower, calculation of alpha_S ratios, PDF ratios
  double sudakov = selected->weightTree(trial, asME, maxScale,
                     selected->clusterIn.pT(), asFSR, asISR,
                     asWeight, pdfWeight);
  // MPI no-emission probability.
  int njetsMaxMPI = mergingHooksPtr->nMinMPI()+1;
  double mpiwt = selected->weightTreeEmissions( trial, -1, 0, njetsMaxMPI,
                   maxScale );

  // Set weight
  wt = ( mergingHooksPtr->nRecluster() == 2 ) ? 1.
     : asWeight*pdfWeight*sudakov*mpiwt;

  // Done
  return wt;

}

//--------------------------------------------------------------------------

double History::weight_UNLOPS_SUBTNLO(PartonLevel* trial, double RN ) {

  // Select a path of clusterings
  History *  selected = select(RN);
  // Set scales in the states to the scales pythia would have set
  selected->setScalesInHistory();
  // So far, no reweighting
  double wt = 1.;
  // Only reweighting with MPI no-emission probability
  double maxScale = (foundCompletePath) ? infoPtr->eCM()
                  : mergingHooksPtr->muFinME();
  int njetsMaxMPI = mergingHooksPtr->nMinMPI()+1;
  double mpiwt = selected->weightTreeEmissions( trial, -1, 0, njetsMaxMPI,
                   maxScale );
  wt = mpiwt;
  // Done
  return wt;

}

//--------------------------------------------------------------------------

// Function to calculate O(\alpha_s)-term of CKKWL-weight for NLO merging

double History::weight_UNLOPS_CORRECTION( int order, PartonLevel* trial,
                  AlphaStrong* asFSR, AlphaStrong* asISR,
                  double RN, Rndm* rndmPtr ) {

  // Already done if no correction should be calculated
  if ( order < 0 ) return 0.;

  // Read alpha_S in ME calculation and maximal scale (eCM)
  double asME   = infoPtr->alphaS();
  double muR      = mergingHooksPtr->muRinME();
  double maxScale = (foundCompletePath)
                  ? infoPtr->eCM()
                  : mergingHooksPtr->muFinME();

  // Pick path of clusterings
  History *  selected = select(RN);
  // Set scales in the states to the scales pythia would have set
  selected->setScalesInHistory();

  double nSteps = mergingHooksPtr->getNumberOfClusteringSteps(state);

  // Get the lowest order k-factor and add first two terms in expansion
  double kFactor = asME * mergingHooksPtr->k1Factor(nSteps);

  // If using Bbar, which includes a tree-level part, subtract an
  // additional one, i.e. the O(\as^0) contribution as well
  double wt = 1.;

  // If only O(\alpha_s^0)-term is to be calculated, done already.
  if ( order == 0 ) return wt;

  // Start by adding the O(\alpha_s^1)-term of the k-factor.
  wt += kFactor;

  // Calculate sum of O(\alpha_s^1)-terms of the ckkw-l weight WITHOUT
  // the O(\alpha_s^1)-term of the last no-emission probability.
  bool fixpdf = true;
  bool fixas  = true;
  // Get first term in expansion of alpha_s ratios.
  double wA   = selected->weightFirstALPHAS( asME, muR, asFSR, asISR );
  // Add logarithm from \alpha_s expansion to weight.
  wt         += (fixas) ? wA : 0.;
  // Generate true average, not only one-point.
  double nWeight = 0.;
  for ( int i=0; i < NTRIAL; ++i ) {
    // Get average number of emissions.
    double wE   = selected->weightFirstEmissions(trial,asME, maxScale,
      asFSR, asISR, fixpdf, fixas );
    // Add average number of emissions off reconstructed states to weight.
    nWeight    += wE;
    // Get first term in expansion of PDF ratios.
    double pscale = selected->clusterIn.pT();
    double wP   = selected->weightFirstPDFs(asME, maxScale, pscale, rndmPtr);
    // Add integral of DGLAP shifted PDF ratios from \alpha_s expansion to wt.
    nWeight    += (fixpdf) ? wP : 0.;
  }
  wt += nWeight/double(NTRIAL);

  // If O(\alpha_s^1)-term + O(\alpha_s^1)-term is to be calculated, done.
  if ( order == 1 ) return wt;

  // So far, no calculation of  O(\alpha_s^2)-term
  return 0.;

}

//--------------------------------------------------------------------------

// Function to set the state with complete scales for evolution.

void History::getStartingConditions( const double RN, Event& outState ) {

  // Select the history
  History *  selected = select(RN);

  // Set scales in the states to the scales pythia would have set
  selected->setScalesInHistory();

  // Copy the output state
  outState = state;

  // Set the scale of the lowest order process
  if (!selected->mother) {
    int nFinal = 0;
    for(int i=0; i < int(outState.size()); ++i)
      if (outState[i].isFinal()) nFinal++;
    if (nFinal <=2)
      outState.scale(mergingHooksPtr->muF());

      // Save information on last splitting, to allow the next
      // emission in the shower to have smaller rapidity with
      // respect to the last ME splitting
      // For hard process, use dummy values
      if (mergingHooksPtr->getNumberOfClusteringSteps(state) == 0) {
        infoPtr->zNowISR(0.5);
        infoPtr->pT2NowISR(pow(state[0].e(),2));
        infoPtr->hasHistory(true);
      // For incomplete process, try to use real values
      } else {
        infoPtr->zNowISR(selected->zISR());
        infoPtr->pT2NowISR(pow(selected->pTISR(),2));
        infoPtr->hasHistory(true);
      }

  } else {

    // Save information on last splitting, to allow the next
    // emission in the shower to have smaller rapidity with
    // respect to the last ME splitting
    infoPtr->zNowISR(selected->zISR());
    infoPtr->pT2NowISR(pow(selected->pTISR(),2));
    infoPtr->hasHistory(true);

  }

  // Save MPI starting scale
  if (mergingHooksPtr->getNumberOfClusteringSteps(state) == 0)
    mergingHooksPtr->muMI(infoPtr->eCM());
  else
    mergingHooksPtr->muMI(outState.scale());

}

//--------------------------------------------------------------------------

// Function to print the history that would be chosen from the number RN.

void History::printHistory( const double RN ) {
  History *  selected = select(RN);
  selected->printStates();
  // Done
}

//--------------------------------------------------------------------------

// Function to print the states in a history, starting from the hard process.

void History::printStates() {
  if ( !mother ) {
    state.list();
    return;
  }
  // Print.
  state.list();
  // Recurse
  mother->printStates();
  // Done
  return;
}

//--------------------------------------------------------------------------

// Function to set the state with complete scales for evolution.

bool History::getClusteredEvent( const double RN, int nSteps,
                Event& outState) {

  // Select history
  History *  selected = select(RN);
  // Set scales in the states to the scales pythia would have set
  // (Only needed if not done before in calculation of weights or
  //  setting of starting conditions)
  selected->setScalesInHistory();
  // If the history does not allow for nSteps clusterings (e.g. because the
  // history is incomplete), return false
  if (nSteps > selected->nClusterings()) return false;
  // Return event with nSteps-1 additional partons (i.e. recluster the last
  // splitting) and copy the output state
  outState = selected->clusteredState(nSteps-1);
  // Done.
  return true;

}

//--------------------------------------------------------------------------

bool History::getFirstClusteredEventAboveTMS( const double RN, int nDesired,
                Event& process, int& nPerformed, bool doUpdate ) {

  // Do reclustering (looping) steps. Remember process scale.
  int nTried  = nDesired - 1;
  // Get number of clustering steps.
  int nSteps   = select(RN)->nClusterings();
  // Set scales in the states to the scales pythia would have set.
  select(RN)->setScalesInHistory();

  // Recluster until reclustered event is above the merging scale.
  Event dummy = Event();
  do {
    // Initialise temporary output of reclustering.
    dummy.clear();
    dummy.init( "(hard process-modified)", particleDataPtr );
    dummy.clear();
    // Recluster once more.
    nTried++;
    // If reclustered event does not exist, exit.
    if ( !getClusteredEvent( RN, nSteps-nTried+1, dummy ) ) return false;
    if ( nTried >= nSteps ) break;

    // Continue loop if reclustered event has unresolved partons.
  } while ( mergingHooksPtr->getNumberOfClusteringSteps(dummy) > 0
         && mergingHooksPtr->rhoms( dummy, false) < mergingHooksPtr->tms() );

  // Update the hard process.
  if ( doUpdate ) process = dummy;

  // Failed to produce output state.
  if ( nTried > nSteps ) return false;

  nPerformed = nTried;
  if ( doUpdate ) {
    // Update to the actual number of steps.
    mergingHooksPtr->nReclusterSave = nPerformed;
    // Save MPI starting scale
    if (mergingHooksPtr->getNumberOfClusteringSteps(state) == 0)
      mergingHooksPtr->muMI(infoPtr->eCM());
    else
      mergingHooksPtr->muMI(state.scale());
  }

  // Done
  return true;

}

//--------------------------------------------------------------------------

// Calculate and return pdf ratio.

double History::getPDFratio( int side, bool forSudakov, bool useHardPDFs,
                    int flavNum, double xNum, double muNum,
                    int flavDen, double xDen, double muDen) {

  // Do nothing for e+e- beams
  if ( abs(flavNum) > 10 && flavNum != 21 ) return 1.0;
  if ( abs(flavDen) > 10 && flavDen != 21 ) return 1.0;

  // Now calculate PDF ratio if necessary
  double pdfRatio = 1.0;

  // Get mother and daughter pdfs
  double pdfNum = 0.0;
  double pdfDen = 0.0;

  // Use hard process PDFs (i.e. PDFs NOT used in ISR, FSR or MPI).
  if ( useHardPDFs ) {
    if (side == 1) {
      if (forSudakov)
        pdfNum = mother->beamA.xfHard( flavNum, xNum, muNum*muNum);
      else
        pdfNum = beamA.xfHard( flavNum, xNum, muNum*muNum);
      if (forSudakov)
        pdfDen = max(1e-10, beamA.xfHard( flavDen, xDen, muDen*muDen));
      else
        pdfDen = max(1e-10, beamA.xfHard( flavDen, xDen, muDen*muDen));
    } else {
      if (forSudakov)
        pdfNum = mother->beamB.xfHard( flavNum, xNum, muNum*muNum);
      else
        pdfNum = beamB.xfHard( flavNum, xNum, muNum*muNum);

      if (forSudakov)
        pdfDen = max(1e-10,beamB.xfHard( flavDen, xDen, muDen*muDen));
      else
        pdfDen = max(1e-10,beamB.xfHard( flavDen, xDen, muDen*muDen));
    }

  // Use rescaled PDFs in the presence of multiparton interactions
  } else {
    if (side == 1) {
      if (forSudakov)
        pdfNum = mother->beamA.xfISR(0, flavNum, xNum, muNum*muNum);
      else
        pdfNum = beamA.xfISR(0, flavNum, xNum, muNum*muNum);
      if (forSudakov)
        pdfDen = max(1e-10, beamA.xfISR(0, flavDen, xDen, muDen*muDen));
      else
        pdfDen = max(1e-10, beamA.xfISR(0, flavDen, xDen, muDen*muDen));

    } else {
      if (forSudakov)
        pdfNum = mother->beamB.xfISR(0, flavNum, xNum, muNum*muNum);
      else
        pdfNum = beamB.xfISR(0, flavNum, xNum, muNum*muNum);

      if (forSudakov)
        pdfDen = max(1e-10,beamB.xfISR(0, flavDen, xDen, muDen*muDen));
      else
        pdfDen = max(1e-10,beamB.xfISR(0, flavDen, xDen, muDen*muDen));
    }
  }

  // Return ratio of pdfs
  if ( pdfNum > 1e-15 && pdfDen > 1e-10 ) {
    pdfRatio *= pdfNum / pdfDen;
  } else if ( pdfNum < pdfDen ) {
    pdfRatio = 0.;
  } else if ( pdfNum > pdfDen ) {
    pdfRatio = 1.;
  }

  // Done
  return pdfRatio;

}

//--------------------------------------------------------------------------

/*--------------- METHODS USED FOR ONLY ONE PATH OF HISTORY NODES ------- */

// Function to set all scales in the sequence of states. This is a
// wrapper routine for setScales and setEventScales methods

void History::setScalesInHistory() {
  // Find correct links from n+1 to n states (mother --> child), as
  // needed for enforcing ordered scale sequences
  vector<int> ident;
  findPath(ident);
  // Set production scales in the states to the scales pythia would
  // have set and enforce ordering
  setScales(ident,true);
  // Set the overall event scales to the scale of the last branching
  setEventScales();
}

//--------------------------------------------------------------------------

// Function to find the index (in the mother histories) of the
// child history, thus providing a way access the path from both
// initial history (mother == 0) and final history (all children == 0)
// IN vector<int> : The index of each child in the children vector
//                  of the current history node will be saved in
//                  this vector
// NO OUTPUT

void History::findPath(vector<int>& out) {

  // If the initial and final nodes are identical, return
  if (!mother && int(children.size()) < 1) return;
  // Find the child by checking the children vector for the perfomed
  // clustering
  int iChild=-1;
  if ( mother ) {
    int size = int(mother->children.size());
    // Loop through children and identify child chosen
    for ( int i=0; i < size; ++i) {
      if ( mother->children[i]->scale == scale
        && mother->children[i]->prob  == prob
        && equalClustering(mother->children[i]->clusterIn,clusterIn)) {
        iChild = i;
        break;
      }
    }
    // Save the index of the child in the children vector and recurse
    if (iChild >-1)
      out.push_back(iChild);
    mother->findPath(out);
  }
}

//--------------------------------------------------------------------------

// Functions to set the  parton production scales and enforce
// ordering on the scales of the respective clusterings stored in
// the History node:
// Method will start from lowest multiplicity state and move to
// higher states, setting the production scales the shower would
// have used.
// When arriving at the highest multiplicity, the method will switch
// and go back in direction of lower states to check and enforce
// ordering for unordered histories.
// IN vector<int> : Vector of positions of the chosen child
//                  in the mother history to allow to move
//                  in direction initial->final along path
//    bool        : True: Move in direction low->high
//                       multiplicity and set production scales
//                  False: Move in direction high->low
//                       multiplicity and check and enforce
//                       ordering
// NO OUTPUT

void History::setScales( vector<int> index, bool forward) {

  // First, set the scales of the hard process to the kinematial
  // limit (=s)
  if ( children.empty() && forward ) {
    // New "incomplete" configurations showered from mu
    if (!mother) {
      double scaleNew = 1.;

      if (mergingHooksPtr->incompleteScalePrescip()==0) {
        scaleNew = mergingHooksPtr->muF();
      } else if (mergingHooksPtr->incompleteScalePrescip()==1) {
        Vec4 pOut;
        pOut.p(0.,0.,0.,0.);
        for(int i=0; i<int(state.size()); ++i)
          if (state[i].isFinal())
            pOut += state[i].p();
        scaleNew = pOut.mCalc();
      } else if (mergingHooksPtr->incompleteScalePrescip()==2) {
        scaleNew = state[0].e();
      }

      scaleNew = max( mergingHooksPtr->pTcut(), scaleNew);

      state.scale(scaleNew);
      for(int i=3; i < int(state.size());++i)
        if (state[i].colType() != 0)
          state[i].scale(scaleNew);
    } else {
      // 2->2 with non-parton particles showered from eCM
      state.scale( state[0].e() );
      // Count final partons
      bool isLEP = ( state[3].isLepton() && state[4].isLepton() );
      int nFinal = 0;
      int nFinalPartons = 0;
      int nFinalPhotons = 0;
      for ( int i=0; i < int(state.size()); ++i ) {
        if ( state[i].isFinal() ) {
          nFinal++;
          if ( state[i].colType() != 0 ) nFinalPartons++;
          if ( state[i].id() == 22 )     nFinalPhotons++;
        }
      }
      bool isQCD = ( nFinal == 2 && nFinal == nFinalPartons );
      bool isPPh = ( nFinal == 2 && nFinalPartons == 1 && nFinalPhotons == 1);
      // If 2->2, purely partonic, set event scale to kinematic pT
      if ( !isLEP && ( isQCD || isPPh ) )
        state.scale( hardFacScale(state) );

    }
  }
  // Set all particle production scales, starting from lowest
  // multiplicity (final) state
  if (mother && forward) {
    // When choosing splitting scale, beware of unordered splittings:
    double scaleNew = 1.;
    if (mergingHooksPtr->unorderedScalePrescip() == 0) {
      // Use larger scale as common splitting scale for mother and child
      scaleNew = max( mergingHooksPtr->pTcut(), max(scale,mother->scale));
    } else if (mergingHooksPtr->unorderedScalePrescip() == 1) {
      // Use smaller scale as common splitting scale for mother and child
      if (scale < mother->scale)
        scaleNew *= max( mergingHooksPtr->pTcut(), min(scale,mother->scale));
      else
        scaleNew *= max( mergingHooksPtr->pTcut(), max(scale,mother->scale));
    }

    // Rescale the mother state partons to the clustering scales
    // that have been found along the path
    mother->state[clusterIn.emitted].scale(scaleNew);
    mother->state[clusterIn.emittor].scale(scaleNew);
    mother->state[clusterIn.recoiler].scale(scaleNew);

    // Find unchanged copies of partons in higher multiplicity states
    // and rescale those
    mother->scaleCopies(clusterIn.emitted, mother->state, scaleNew);
    mother->scaleCopies(clusterIn.emittor, mother->state, scaleNew);
    mother->scaleCopies(clusterIn.recoiler, mother->state, scaleNew);

    // Recurse
    mother->setScales(index,true);
  }

  // Now, check and correct ordering from the highest multiplicity
  // state backwards to all the clustered states
  if (!mother || !forward) {
    // Get index of child along the path
    int iChild = -1;
    if ( int(index.size()) > 0 ) {
      iChild = index.back();
      index.pop_back();
    }

    // Check that the reclustered scale is above the shower cut
    if (mother) {
      scale = max(mergingHooksPtr->pTcut(), scale);
    }
    // If this is NOT the 2->2 process, check and enforce ordering
    if (iChild != -1 && !children.empty()) {

      if (scale > children[iChild]->scale ) {
        if (mergingHooksPtr->unorderedScalePrescip() == 0) {
          // Use larger scale as common splitting scale for mother and child
          double scaleNew = max( mergingHooksPtr->pTcut(),
                              max(scale,children[iChild]->scale));
          // Enforce ordering in particle production scales
          for( int i = 0; i < int(children[iChild]->state.size()); ++i)
            if (children[iChild]->state[i].scale() == children[iChild]->scale)
              children[iChild]->state[i].scale(scaleNew);
          // Enforce ordering in saved clustering scale
          children[iChild]->scale = scaleNew;

        } else if ( mergingHooksPtr->unorderedScalePrescip() == 1) {
           // Use smaller scale as common splitting scale for mother & child
           double scaleNew = max(mergingHooksPtr->pTcut(),
                               min(scale,children[iChild]->scale));
           // Enforce ordering in particle production scales
           for( int i = 0; i < int(state.size()); ++i)
             if (state[i].scale() == scale)
               state[i].scale(scaleNew);
           // Enforce ordering in saved clustering scale
           scale = scaleNew;
        }
      // Just set the overall event scale to the minimal scale
      } else {
        double scalemin = state[0].e();
        for( int i = 0; i < int(state.size()); ++i)
          if (state[i].colType() != 0)
            scalemin = max(mergingHooksPtr->pTcut(),
                         min(scalemin,state[i].scale()));
        state.scale(scalemin);
        scale = max(mergingHooksPtr->pTcut(), scale);
      }
      //Recurse
      children[iChild]->setScales(index, false);
    }
  }

}

//--------------------------------------------------------------------------

// Function to find a particle in all higher multiplicity events
// along the history path and set its production scale to the input
// scale
// IN  int iPart       : Parton in refEvent to be checked / rescaled
//     Event& refEvent : Reference event for iPart
//     double scale    : Scale to be set as production scale for
//                       unchanged copies of iPart in subsequent steps

void History::scaleCopies(int iPart, const Event& refEvent, double rho) {

  // Check if any parton recently rescaled is found unchanged:
  // Same charge, colours in mother->state
  if ( mother ) {
    for( int i=0; i < mother->state.size(); ++i) {
      if ( ( mother->state[i].id()         == refEvent[iPart].id()
          && mother->state[i].colType()    == refEvent[iPart].colType()
          && mother->state[i].chargeType() == refEvent[iPart].chargeType()
          && mother->state[i].col()        == refEvent[iPart].col()
          && mother->state[i].acol()       == refEvent[iPart].acol() )
         ) {
        // Rescale the unchanged parton
        mother->state[i].scale(rho);
        // Recurse
         if (mother->mother)
          mother->scaleCopies( iPart, refEvent, rho );
       } // end if found unchanged parton case
    } // end loop over particle entries in event
  }
}

//--------------------------------------------------------------------------

// Functions to set the OVERALL EVENT SCALES [=state.scale()] to
// the scale of the last clustering
// NO INPUT
// NO OUTPUT

void History::setEventScales() {
  // Set the event scale to the scale of the last clustering,
  // except for the very lowest multiplicity state
  if (mother) {
    mother->state.scale(scale);
    // Recurse
    mother->setEventScales();
  }
}

//--------------------------------------------------------------------------

// Functions to return the z value of the last ISR splitting
// NO INPUT
// OUTPUT double : z value of last ISR splitting in history

double History::zISR() {

  // Do nothing for ME level state
  if (!mother) return 0.0;
  // Skip FSR splitting
  if (mother->state[clusterIn.emittor].isFinal()) return mother->zISR();
  // Calculate z
  int rad = clusterIn.emittor;
  int rec = clusterIn.recoiler;
  int emt = clusterIn.emitted;
  double z = (mother->state[rad].p() + mother->state[rec].p()
            - mother->state[emt].p()).m2Calc()
    / (mother->state[rad].p() + mother->state[rec].p()).m2Calc();
  // Recurse
  double znew = mother->zISR();
  // Update z
  if (znew > 0.) z = znew;

  return z;
}

//--------------------------------------------------------------------------

// Functions to return the z value of the last FSR splitting
// NO INPUT
// OUTPUT double : z value of last FSR splitting in history

double History::zFSR() {

  // Do nothing for ME level state
  if (!mother) return 0.0;
  // Skip ISR splitting
  if (!mother->state[clusterIn.emittor].isFinal()) return mother->zFSR();
  // Calculate z
  int rad = clusterIn.emittor;
  int rec = clusterIn.recoiler;
  int emt = clusterIn.emitted;
  // Construct 2->3 variables for FSR
  Vec4   sum = mother->state[rad].p() + mother->state[rec].p()
             + mother->state[emt].p();
  double m2Dip = sum.m2Calc();
  double x1 = 2. * (sum * mother->state[rad].p()) / m2Dip;
  double x3 = 2. * (sum * mother->state[emt].p()) / m2Dip;
  // Calculate z of splitting for FSR
  double z = x1/(x1+x3);
  // Recurse
  double znew = mother->zFSR();
  // Update z
  if (znew > 0.) z = znew;

  return z;
}

//--------------------------------------------------------------------------

// Functions to return the pT scale of the last FSR splitting
// NO INPUT
// OUTPUT double : pT scale of last FSR splitting in history

double History::pTISR() {
  // Do nothing for ME level state
  if (!mother) return 0.0;
  // Skip FSR splitting
  if (mother->state[clusterIn.emittor].isFinal()) return mother->pTISR();
  double pT = mother->state.scale();
  // Recurse
  double pTnew = mother->pTISR();
  // Update pT
  if (pTnew > 0.) pT = pTnew;

  return pT;
}

//--------------------------------------------------------------------------

// Functions to return the pT scale of the last FSR splitting
// NO INPUT
// OUTPUT double : pT scale of last FSR splitting in history

double History::pTFSR() {

  // Do nothing for ME level state
  if (!mother) return 0.0;
  // Skip ISR splitting
  if (!mother->state[clusterIn.emittor].isFinal()) return mother->pTFSR();
  double pT = mother->state.scale();
  // Recurse
  double pTnew = mother->pTFSR();
  // Update pT
  if (pTnew > 0.) pT = pTnew;
  return pT;
}

//--------------------------------------------------------------------------

// Function to return the depth of the history (i.e. the number of
// reclustered splittings)
// NO INPUT
// OUTPUT int  : Depth of history

int History::nClusterings() {
  if (!mother) return 0;
  int w = mother->nClusterings();
  w += 1;
  return w;
}

//--------------------------------------------------------------------------

// Functions to return the event after nSteps splittings of the 2->2 process
// Example: nSteps = 1 -> return event with one additional parton
// INPUT  int   : Number of splittings in the event,
//                as counted from core 2->2 process
// OUTPUT Event : event with nSteps additional partons

Event History::clusteredState(int nSteps) {

  // Save state
  Event outState = state;
  // As long as there are steps to do, recursively save state
  if (mother && nSteps > 0)
    outState = mother->clusteredState(nSteps - 1);
  // Done
  return outState;

}

//--------------------------------------------------------------------------

// Function to choose a path from all paths in the tree
// according to their splitting probabilities
// IN double    : Random number
// OUT History* : Leaf of history path chosen

History * History::select(double rnd) {

  // No need to choose if no paths have been constructed.
  if ( goodBranches.empty() && badBranches.empty() ) return this;

  // Choose amongst paths allowed by projections.
  double sum = 0.;
  map<double, History*> selectFrom;
  if ( !goodBranches.empty() ) {
    selectFrom = goodBranches;
    sum        = sumGoodBranches;
  } else {
    selectFrom = badBranches;
    sum        = sumBadBranches;
  }

  if (mergingHooksPtr->pickBySumPT()) {
    // Find index of history with minimal sum of scalar pT
    int nFinal = 0;
    for (int i=0; i < state.size(); ++i)
      if (state[i].isFinal())
        nFinal++;
    double iMin = 0.;
    double sumMin = (nFinal-2)*state[0].e();
    for ( map<double, History*>::iterator it = selectFrom.begin();
      it != selectFrom.end(); ++it ) {

      if (it->second->sumScalarPT < sumMin) {
        sumMin = it->second->sumScalarPT;
        iMin = it->first;
      }
    }
    // Choose history with smallest sum of scalar pT
    return selectFrom.lower_bound(iMin)->second;
  } else {
    // Choose history according to probability, be careful about upper bound
    if ( rnd != 1. ) {
      return selectFrom.upper_bound(sum*rnd)->second;
    } else {
      return selectFrom.lower_bound(sum*rnd)->second;
    }
  }
  // Done
}

//--------------------------------------------------------------------------

// Function to project paths onto desired paths.

bool History::trimHistories() {
  // Do nothing if no paths have been constructed.
  if ( paths.empty() ) return false;
  // Loop through all constructed paths. Check all removal conditions.
  for ( map<double, History*>::iterator it = paths.begin();
    it != paths.end(); ++it ) {
    // Check if history is allowed.
    if ( it->second->keep() && !it->second->keepHistory() )
      it->second->remove();
  }
  // Project onto desired / undesired branches.
  double sumold, sumnew, sumprob, mismatch;
  sumold = sumnew = sumprob = mismatch = 0.;
  // Loop through all constructed paths and store allowed paths.
  // Skip undesired paths.
  for ( map<double, History*>::iterator it = paths.begin();
    it != paths.end(); ++it ) {
    // Update index
    sumnew = it->first;
    if ( it->second->keep() ) {
      // Fill branches with allowed paths.
      goodBranches.insert( make_pair( sumnew - mismatch, it->second) );
      // Add probability of this path.
      sumGoodBranches = sumnew - mismatch;
    } else {
      // Update mismatch in probabilities resulting from not including this
      // path
      double mismatchOld = mismatch;
      mismatch += sumnew - sumold;
      // Fill branches with allowed paths.
      badBranches.insert( make_pair( mismatchOld + sumnew - sumold,
        it->second ) );
      // Add probability of this path.
      sumBadBranches = mismatchOld  + sumnew - sumold;
    }
    // remember index of this path in order to caclulate probability of
    // subsequent path.
    sumold = it->first;
  }
  // Done
  return !goodBranches.empty();
}

//--------------------------------------------------------------------------

// Function implementing checks on a paths, deciding if the path is valid.

bool History::keepHistory() {
  bool keepPath = true;
  // Tag unordered paths for removal.
  if ( mergingHooksPtr->getProcessString().compare("pp>jj") == 0
    || mergingHooksPtr->getProcessString().compare("pp>aj") == 0 ) {
    // Tag unordered paths for removal. Include scale of hard 2->2 process
    // into the ordering definition.
    double maxScale = hardFacScale(state);
    keepPath = isOrderedPath( maxScale );
  }

  //Done
  return keepPath;
}

//--------------------------------------------------------------------------

// Function to check if a path is ordered in evolution pT.

bool History::isOrderedPath( double maxscale ) {
  double newscale = clusterIn.pT();
  if ( !mother ) return true;
  bool ordered = mother->isOrderedPath(newscale);
  if ( !ordered || maxscale < newscale) return false;
  return ordered;
}

//--------------------------------------------------------------------------

// Function to check if all reconstucted states in a path pass the merging
// scale cut.

bool History::allIntermediateAboveRhoMS( double rhoms, bool good ) {
  // If one state below the merging scale has already been found, no need to
  // check further.
  if ( !good ) return false;
  // Check merging scale for states with more than 0 jets
  int nFinal = 0;
  for ( int i = 0; i < state.size(); ++i )
    if ( state[i].isFinal() && state[i].colType() != 0 )
      nFinal++;
  double rhoNew = (nFinal > 0 ) ? mergingHooksPtr->rhoms( state, false )
                : state[0].e();
  // Assume state from ME generator passes merging scale cut.
  if ( !mother ) return good;
  // Recurse.
  return good && mother->allIntermediateAboveRhoMS( rhoms, (rhoNew > rhoms) );
}

//--------------------------------------------------------------------------

// Function to check if any ordered paths were found (and kept).

bool History::foundAnyOrderedPaths() {
  //Do nothing if no paths were found
  if ( paths.empty() ) return false;
  double maxscale = infoPtr->eCM();
  // Loop through paths. Divide probability into ordered and unordered pieces.
  for ( map<double, History*>::iterator it = paths.begin();
    it != paths.end(); ++it )
    if ( it->second->isOrderedPath(maxscale) )
      return true;
  // Done
  return false;
}

//--------------------------------------------------------------------------

// For a full path, find the weight calculated from the ratio of
// couplings, the no-emission probabilities, and possible PDF
// ratios. This function should only be called for the last history
// node of a full path.
// IN  TimeShower : Already initialised shower object to be used as
//                  trial shower
//     double     : alpha_s value used in ME calculation
//     double     : Maximal mass scale of the problem (e.g. E_CM)
//     AlphaStrong: Initialised shower alpha_s object for FSR
//                  alpha_s ratio calculation
//     AlphaStrong: Initialised shower alpha_s object for ISR
//                  alpha_s ratio calculation (can be different from previous)

double History::weightTree(PartonLevel* trial, double as0, double maxscale,
  double pdfScale, AlphaStrong * asFSR, AlphaStrong * asISR,
  double& asWeight, double& pdfWeight) {

  // Use correct scale
  double newScale = scale;

  // For ME state, just multiply by PDF ratios
  if ( !mother ) {

    int sideRad = (state[3].pz() > 0) ? 1 :-1;
    int sideRec = (state[4].pz() > 0) ? 1 :-1;
    // Calculate PDF first leg
    if (state[3].colType() != 0) {
      // Find x value and flavour
      double x = 2.*state[3].e() / state[0].e();
      int flav = state[3].id();

      // Find numerator/denominator scale
      double scaleNum = (children.empty()) ? hardFacScale(state) : maxscale;
      double scaleDen = mergingHooksPtr->muFinME();
      // For initial parton, multiply by PDF ratio
      double ratio = getPDFratio(sideRad, false, false, flav, x, scaleNum,
                       flav, x, scaleDen);
      pdfWeight *= ratio;
    }

    // Calculate PDF ratio for second leg
    if (state[4].colType() != 0) {
      // Find x value and flavour
      double x = 2.*state[4].e() / state[0].e();
      int flav = state[4].id();
      // Find numerator/denominator scale
      double scaleNum = (children.empty()) ? hardFacScale(state) : maxscale;
      double scaleDen = mergingHooksPtr->muFinME();
      // For initial parton, multiply with PDF ratio
      double ratio = getPDFratio(sideRec, false, false, flav, x, scaleNum,
                       flav, x, scaleDen);
      pdfWeight *= ratio;
    }

    return 1.0;
  }

  // Remember new PDF scale n case true sclae should be used for un-ordered
  // splittings.
  double newPDFscale = newScale;
  if (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
    newPDFscale = clusterIn.pT();

  // Recurse
  double w = mother->weightTree(trial, as0, newScale, newPDFscale,
                       asFSR, asISR, asWeight, pdfWeight);

  // Do nothing for empty state
  if (state.size() < 3) return 1.0;
  // If up to now, trial shower was not successful, return zero
  if ( w < 1e-12 ) return 0.0;
  // Do trial shower on current state, return zero if not successful
  w *= doTrialShower(trial, 1, maxscale);
  if ( w < 1e-12 ) return 0.0;

  // Calculate alpha_s ratio for current state
  if ( asFSR && asISR ) {
    double asScale = pow2( newScale );
    if (mergingHooksPtr->unorderedASscalePrescip() == 1)
      asScale = pow2( clusterIn.pT() );
    bool FSR = mother->state[clusterIn.emittor].isFinal();
    double alphaSinPS = (FSR) ? (*asFSR).alphaS(asScale)
                      : (*asISR).alphaS(asScale
                                       + pow2(mergingHooksPtr->pT0ISR()) );
    asWeight *= alphaSinPS / as0;
  }

  // Calculate pdf ratios: Get both sides of event
  int inP = 3;
  int inM = 4;
  int sideP = (mother->state[inP].pz() > 0) ? 1 :-1;
  int sideM = (mother->state[inM].pz() > 0) ? 1 :-1;

  if ( mother->state[inP].colType() != 0 ) {
    // Find x value and flavour
    double x = getCurrentX(sideP);
    int flav = getCurrentFlav(sideP);
    // Find numerator scale
    double scaleNum = (children.empty())
                    ? hardFacScale(state)
                    : ( (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                      ? pdfScale : maxscale );
    double scaleDen = (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                    ? clusterIn.pT() : newScale;
    // Multiply PDF ratio
    double ratio = getPDFratio(sideP, false, false, flav, x, scaleNum,
                     flav, x, scaleDen);
    pdfWeight *= ratio;
  }

  if ( mother->state[inM].colType() != 0 ) {
    // Find x value and flavour
    double x = getCurrentX(sideM);
    int flav = getCurrentFlav(sideM);
    // Find numerator scale
    double scaleNum = (children.empty())
                    ? hardFacScale(state)
                    : ( (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                      ? pdfScale : maxscale );
    double scaleDen = (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                    ? clusterIn.pT() : newScale;
    // Multiply PDF ratio
    double ratio = getPDFratio(sideM, false, false, flav, x, scaleNum,
                     flav, x, scaleDen);
    pdfWeight *= ratio;
  }

  // Done
  return w;
}

//--------------------------------------------------------------------------

// Function to return the \alpha_s-ratio part of the CKKWL weight of a path.

double History::weightTreeALPHAS( double as0, AlphaStrong * asFSR,
  AlphaStrong * asISR ) {

  // For ME state, do nothing.
  if ( !mother ) return 1.;
  // Recurse
  double w = mother->weightTreeALPHAS( as0, asFSR, asISR );
  // Do nothing for empty state
  if (state.size() < 3) return w;

  // Calculate alpha_s ratio for current state
  if ( asFSR && asISR ) {
    double asScale = pow2( scale );
    if (mergingHooksPtr->unorderedASscalePrescip() == 1)
      asScale = pow2( clusterIn.pT() );
    bool FSR = mother->state[clusterIn.emittor].isFinal();
    double alphaSinPS = (FSR)
                ? (*asFSR).alphaS(asScale)
                : (*asISR).alphaS(asScale + pow2(mergingHooksPtr->pT0ISR()) );
    w *= alphaSinPS / as0;
  }

  // Done
  return w;
}

//--------------------------------------------------------------------------

// Function to return the PDF-ratio part of the CKKWL weight of a path.

double History::weightTreePDFs( double maxscale, double pdfScale ) {

  // Use correct scale
  double newScale = scale;

  // For ME state, just multiply by PDF ratios
  if ( !mother ) {

    double wt = 1.;
    int sideRad = (state[3].pz() > 0) ? 1 :-1;
    int sideRec = (state[4].pz() > 0) ? 1 :-1;

    // Calculate PDF first leg
    if (state[3].colType() != 0) {
      // Find x value and flavour
      double x = 2.*state[3].e() / state[0].e();
      int flav = state[3].id();
      // Find numerator/denominator scale
      double scaleNum = (children.empty()) ? hardFacScale(state) : maxscale;
      double scaleDen = mergingHooksPtr->muFinME();
      // For initial parton, multiply by PDF ratio
      wt *= getPDFratio(sideRad, false, false, flav, x, scaleNum, flav, x,
              scaleDen);
    }

    // Calculate PDF ratio for second leg
    if (state[4].colType() != 0) {
      // Find x value and flavour
      double x = 2.*state[4].e() / state[0].e();
      int flav = state[4].id();
      // Find numerator/denominator scale
      double scaleNum = (children.empty()) ? hardFacScale(state) : maxscale;
      double scaleDen = mergingHooksPtr->muFinME();
      // For initial parton, multiply with PDF ratio
      wt *= getPDFratio(sideRec, false, false, flav, x, scaleNum, flav, x,
              scaleDen);
    }

    return wt;
  }

  // Remember new PDF scale n case true sclae should be used for un-ordered
  // splittings.
  double newPDFscale = newScale;
  if ( mergingHooksPtr->unorderedPDFscalePrescip() == 1)
    newPDFscale = clusterIn.pT();

  // Recurse
  double w = mother->weightTreePDFs( newScale, newPDFscale );

  // Do nothing for empty state
  if (state.size() < 3) return w;

  // Calculate pdf ratios: Get both sides of event
  int inP = 3;
  int inM = 4;
  int sideP = (mother->state[inP].pz() > 0) ? 1 :-1;
  int sideM = (mother->state[inM].pz() > 0) ? 1 :-1;

  if ( mother->state[inP].colType() != 0 ) {
    // Find x value and flavour
    double x = getCurrentX(sideP);
    int flav = getCurrentFlav(sideP);
    // Find numerator scale
    double scaleNum = (children.empty())
                    ? hardFacScale(state)
                    : ( (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                      ? pdfScale : maxscale );
    double scaleDen = (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                    ? clusterIn.pT() : newScale;
    // Multiply PDF ratio
    double ratio = getPDFratio(sideP, false, false, flav, x, scaleNum,
                     flav, x, scaleDen);
    w *= ratio;
  }

  if ( mother->state[inM].colType() != 0 ) {
    // Find x value and flavour
    double x = getCurrentX(sideM);
    int flav = getCurrentFlav(sideM);
    // Find numerator scale
    double scaleNum = (children.empty())
                    ? hardFacScale(state)
                    : ( (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                      ? pdfScale : maxscale );
    double scaleDen = (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                    ? clusterIn.pT() : newScale;
    // Multiply PDF ratio
    double ratio = getPDFratio(sideM, false, false, flav, x, scaleNum,
                     flav, x, scaleDen);
    w *= ratio;
  }

  // Done
  return w;
}

//--------------------------------------------------------------------------

// Function to return the no-emission probability part of the CKKWL weight.

double History::weightTreeEmissions( PartonLevel* trial, int type,
  int njetMin, int njetMax, double maxscale ) {

  // Use correct scale
  double newScale = scale;
  // For ME state, just multiply by PDF ratios
  if ( !mother ) return 1.0;
  // Recurse
  double w = mother->weightTreeEmissions(trial,type,njetMin,njetMax,newScale);
  // Do nothing for empty state
  if (state.size() < 3) return 1.0;
  // If up to now, trial shower was not successful, return zero
  if ( w < 1e-12 ) return 0.0;
  // If this node has too many jets, no not calculate no-emission probability.
  int njetNow = mergingHooksPtr->getNumberOfClusteringSteps( state) ;
  if (njetNow >= njetMax) return 1.0;

  if (njetNow < njetMin ) w *= 1.0;
  // Do trial shower on current state, return zero if not successful
  else w *= doTrialShower(trial, type, maxscale);

  if ( w < 1e-12 ) return 0.0;
  // Done
  return w;

}

//--------------------------------------------------------------------------

// Function to generate the O(\alpha_s)-term of the CKKWL-weight.

double History::weightFirst(PartonLevel* trial, double as0, double muR,
  double maxscale, AlphaStrong * asFSR, AlphaStrong * asISR, Rndm* rndmPtr ) {

  // Use correct scale
  double newScale = scale;

  if ( !mother ) {

    double weight = 0.;

    // Calculate PDF first leg
    if (state[3].colType() != 0) {
      // Find x value and flavour
      double x = 2.*state[3].e() / state[0].e();
      int flav = state[3].id();
      // Find numerator/denominator scale
      double scaleNum = (children.empty()) ? hardFacScale(state) : maxscale;
      double scaleDen = mergingHooksPtr->muFinME();
      // Monte Carlo integrand.
      double intPDF4 = monteCarloPDFratios(flav, x, scaleNum, scaleDen,
                         mergingHooksPtr->muFinME(), as0, rndmPtr);
      weight += intPDF4;
    }

    // Calculate PDF ratio for second leg
    if (state[4].colType() != 0) {
      // Find x value and flavour
      double x = 2.*state[4].e() / state[0].e();
      int flav = state[4].id();
      // Find numerator/denominator scale
      double scaleNum = (children.empty()) ? hardFacScale(state) : maxscale;
      double scaleDen = mergingHooksPtr->muFinME();
      // Monte Carlo integrand.
      double intPDF4 = monteCarloPDFratios(flav, x, scaleNum, scaleDen,
                         mergingHooksPtr->muFinME(), as0, rndmPtr);
      weight += intPDF4;
    }

    return weight;
  }

  // Recurse
  double w = mother->weightFirst(trial, as0, muR, newScale, asFSR, asISR,
               rndmPtr );

  // Do nothing for empty state
  if (state.size() < 3) return 0.0;

  // Find right scale
  double b = 1.;
  double asScale2 = newScale*newScale;
  int showerType = (mother->state[clusterIn.emittor].isFinal() ) ? 1 : -1;
  if (showerType == -1) {
    asScale2 += pow(mergingHooksPtr->pT0ISR(),2);
    b = 1.;
  }
  // Find summand beta_0 / 2 * ln(muR^2/t_i) due to as expansion.
  double NF = 4.;
  double BETA0 = 11. - 2./3.* NF;
  // For fixed \alpha_s in matrix element
  w += as0 / (2.*M_PI) * 0.5 * BETA0 * log( (muR*muR) / (b*asScale2) );

  // Count emissions: New variant
  // Generate true average, not only one-point.
  bool fixpdf = true;
  bool fixas  = true;
  double nWeight1 = 0.;
  double nWeight2 = 0.;

  for(int i=0; i < NTRIAL; ++i) {
    // Get number of emissions
    vector<double> unresolvedEmissionTerm = countEmissions(trial, maxscale,
      newScale, 2, as0, asFSR, asISR, 3, fixpdf, fixas);
    nWeight1 += unresolvedEmissionTerm[1];
  }
  w += nWeight1/double(NTRIAL) + nWeight2/double(NTRIAL);

  // Calculate pdf ratios: Get both sides of event
  int inP = 3;
  int inM = 4;
  int sideP = (mother->state[inP].pz() > 0) ? 1 :-1;
  int sideM = (mother->state[inM].pz() > 0) ? 1 :-1;

  if ( mother->state[inP].colType() != 0 ) {
    // Find x value and flavour
    double x = getCurrentX(sideP);
    int flav = getCurrentFlav(sideP);
    // Find numerator scale
    double scaleNum = (children.empty()) ? hardFacScale(state) : maxscale;
    // Monte Carlo integrand.
    double intPDF4 = monteCarloPDFratios(flav, x, scaleNum, newScale,
                       mergingHooksPtr->muFinME(), as0, rndmPtr);
    w += intPDF4;

  }

  if ( mother->state[inM].colType() != 0 ) {
    // Find x value and flavour
    double x = getCurrentX(sideM);
    int flav = getCurrentFlav(sideM);
    // Find numerator scale
    double scaleNum = (children.empty()) ? hardFacScale(state) : maxscale;
    // Monte Carlo integrand.
    double intPDF4 = monteCarloPDFratios(flav, x, scaleNum, newScale,
                       mergingHooksPtr->muFinME(), as0, rndmPtr);
    w += intPDF4;

  }

  // Done
  return w;

}

//--------------------------------------------------------------------------

// Function to generate the O(\alpha_s)-term of the \alpha_s-ratios
// appearing in the CKKWL-weight.

double History::weightFirstALPHAS( double as0, double muR,
  AlphaStrong * asFSR, AlphaStrong * asISR ) {

  // Use correct scale
  double newScale = scale;
  // Done
  if ( !mother ) return 0.;
  // Recurse
  double w = mother->weightFirstALPHAS( as0, muR, asFSR, asISR );
  // Find right scale
  int showerType = (mother->state[clusterIn.emittor].isFinal() ) ? 1 : -1;
  double b = 1.;
  double asScale = pow2( newScale );
  if ( mergingHooksPtr->unorderedASscalePrescip() == 1 )
    asScale = pow2( clusterIn.pT() );
  if (showerType == -1) {
    asScale += pow2( mergingHooksPtr->pT0ISR() );
    b = 1.;
  }
  // Find summand beta_0 / 2 * ln(muR^2/t_i) due to as expansion.
  double NF = 4.;
  double BETA0 = 11. - 2./3.* NF;
  // For fixed \alpha_s in matrix element
  w += as0 / (2.*M_PI) * 0.5 * BETA0 * log( (muR*muR) / (b*asScale) );
  // Done
  return w;

}

//--------------------------------------------------------------------------

// Function to generate the O(\alpha_s)-term of the PDF-ratios
// appearing in the CKKWL-weight.

double History::weightFirstPDFs( double as0, double maxscale, double pdfScale,
  Rndm* rndmPtr ) {

  // Use correct scale
  double newScale = scale;

  if ( !mother ) {

    double wt = 0.;

    // Calculate PDF first leg
    if (state[3].colType() != 0) {
      // Find x value and flavour
      double x        = 2.*state[3].e() / state[0].e();
      int flav        = state[3].id();
      // Find numerator/denominator scale
      double scaleNum = (children.empty()) ? hardFacScale(state) : maxscale;
      double scaleDen = mergingHooksPtr->muFinME();
      // Monte Carlo integrand.
      wt += monteCarloPDFratios(flav, x, scaleNum, scaleDen,
                          mergingHooksPtr->muFinME(), as0, rndmPtr);
    }
    // Calculate PDF ratio for second leg
    if (state[4].colType() != 0) {
      // Find x value and flavour
      double x        = 2.*state[4].e() / state[0].e();
      int flav        = state[4].id();
      // Find numerator/denominator scale
      double scaleNum = (children.empty()) ? hardFacScale(state) : maxscale;
      double scaleDen = mergingHooksPtr->muFinME();
      // Monte Carlo integrand.
      wt += monteCarloPDFratios(flav, x, scaleNum, scaleDen,
                         mergingHooksPtr->muFinME(), as0, rndmPtr);
    }

    // Done
    return wt;
  }

  // Remember new PDF scale n case true sclae should be used for un-ordered
  // splittings.
  double newPDFscale = newScale;
  if (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
    newPDFscale      = clusterIn.pT();

  // Recurse
  double w = mother->weightFirstPDFs( as0, newScale, newPDFscale, rndmPtr);

  // Calculate pdf ratios: Get both sides of event
  int inP   = 3;
  int inM   = 4;
  int sideP = (mother->state[inP].pz() > 0) ? 1 :-1;
  int sideM = (mother->state[inM].pz() > 0) ? 1 :-1;

  if ( mother->state[inP].colType() != 0 ) {
    // Find x value and flavour
    double x        = getCurrentX(sideP);
    int flav        = getCurrentFlav(sideP);
    // Find numerator / denominator scales
    double scaleNum = (children.empty())
                    ? hardFacScale(state)
                    : ( (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                      ? pdfScale : maxscale );
    double scaleDen = (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                    ? clusterIn.pT() : newScale;
    // Monte Carlo integrand.
    w += monteCarloPDFratios(flav, x, scaleNum, scaleDen,
                        mergingHooksPtr->muFinME(), as0, rndmPtr);
  }

  if ( mother->state[inM].colType() != 0 ) {
    // Find x value and flavour
    double x        = getCurrentX(sideM);
    int flav        = getCurrentFlav(sideM);
    // Find numerator / denominator scales
    double scaleNum = (children.empty())
                    ? hardFacScale(state)
                    : ( (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                      ? pdfScale : maxscale );
    double scaleDen = (mergingHooksPtr->unorderedPDFscalePrescip() == 1)
                    ? clusterIn.pT() : newScale;
    // Monte Carlo integrand.
    w += monteCarloPDFratios(flav, x, scaleNum, scaleDen,
                        mergingHooksPtr->muFinME(), as0, rndmPtr);
  }

  // Done
  return w;

}


//--------------------------------------------------------------------------

// Function to generate the O(\alpha_s)-term of the no-emission
// probabilities appearing in the CKKWL-weight.

double History::weightFirstEmissions(PartonLevel* trial, double as0,
  double maxscale, AlphaStrong * asFSR, AlphaStrong * asISR,
  bool fixpdf, bool fixas ) {

  // Use correct scale
  double newScale = scale;
  if ( !mother ) return 0.0;
  // Recurse
  double w = mother->weightFirstEmissions(trial, as0, newScale, asFSR, asISR,
                                          fixpdf, fixas );
  // Do nothing for empty state
  if (state.size() < 3) return 0.0;
  // Generate true average.
  double nWeight1 = 0.;
  double nWeight2 = 0.;
  for(int i=0; i < NTRIAL; ++i) {
    // Get number of emissions
    vector<double> unresolvedEmissionTerm = countEmissions(trial, maxscale,
      newScale, 2, as0, asFSR, asISR, 3, fixpdf, fixas);
    nWeight1 += unresolvedEmissionTerm[1];
  }

  w += nWeight1/double(NTRIAL) + nWeight2/double(NTRIAL);

  // Done
  return w;

}

//--------------------------------------------------------------------------

// Function to return the factorisation scale of the hard process in Pythia.

double History::hardFacScale(const Event& event) {
  // Declare output scale.
  double hardscale = 0.;
  // If scale should not be reset, done.
  if ( !mergingHooksPtr->resetHardQFac() ) return mergingHooksPtr->muF();
  // For pure QCD dijet events, calculate the hadronic cross section
  // of the hard process at the pT of the dijet system, rather than at fixed
  // arbitrary scale.
  if ( mergingHooksPtr->getProcessString().compare("pp>jj") == 0
    || mergingHooksPtr->getProcessString().compare("pp>aj") == 0 ) {
    // Find the mT in the hard sub-process.
    vector <double> mT;
    for ( int i=0; i < event.size(); ++i)
      if ( event[i].isFinal() && event[i].colType() != 0 )
        mT.push_back( abs(event[i].mT2()) );
    if ( int(mT.size()) != 2 )
      hardscale = infoPtr->QFac();
    else
      hardscale = sqrt( min( mT[0], mT[1] ) );
  } else {
    hardscale = mergingHooksPtr->muF();
  }
  // Done
  return hardscale;
}

//--------------------------------------------------------------------------

// Function to return the factorisation scale of the hard process in Pythia.

double History::hardRenScale(const Event& event) {
  // Declare output scale.
  double hardscale = 0.;
  // If scale should not be reset, done.
  if ( !mergingHooksPtr->resetHardQRen() ) return mergingHooksPtr->muR();
  // For pure QCD dijet events, calculate the hadronic cross section
  // of the hard process at the pT of the dijet system, rather than at fixed
  // arbitrary scale.
  if ( mergingHooksPtr->getProcessString().compare("pp>jj") == 0
    || mergingHooksPtr->getProcessString().compare("pp>aj") == 0 ) {
    // Find the mT in the hard sub-process.
    vector <double> mT;
    for ( int i=0; i < event.size(); ++i)
      if ( event[i].isFinal()
        && ( event[i].colType() != 0 || event[i].id() == 22 ) )
        mT.push_back( abs(event[i].mT()) );
    if ( int(mT.size()) != 2 )
      hardscale = infoPtr->QRen();
    else
      hardscale = sqrt( mT[0]*mT[1] );
  } else {
    hardscale = mergingHooksPtr->muR();
  }
  // Done
  return hardscale;
}

//--------------------------------------------------------------------------

// Perform a trial shower using the \a pythia object between
// maxscale down to this scale and return the corresponding Sudakov
// form factor.
// IN  trialShower : Shower object used as trial shower
//     double     : Maximum scale for trial shower branching
// OUT  0.0       : trial shower emission outside allowed pT range
//      1.0       : trial shower successful (any emission was below
//                  the minimal scale )

double History::doTrialShower( PartonLevel* trial, int type,
  double maxscaleIn, double minscaleIn ) {

  // Copy state to local process
  Event process        = state;
  // Set starting scale.
  double startingScale = maxscaleIn;
  // Careful when setting shower starting scale for pure QCD and prompt
  // photon case.
  if ( mergingHooksPtr->getNumberOfClusteringSteps(process) == 0
    && ( mergingHooksPtr->getProcessString().compare("pp>jj") == 0
      || mergingHooksPtr->getProcessString().compare("pp>aj") == 0 ) )
      startingScale = min( startingScale, hardFacScale(process) );

  // Set output.
  bool doVeto          = false;

  while ( true ) {

    // Reset trialShower object
    trial->resetTrial();
    // Construct event to be showered
    Event event = Event();
    event.init("(hard process-modified)", particleDataPtr);
    event.clear();

    // Reset process scale so that shower starting scale is correctly set.
    process.scale(startingScale);
    doVeto = false;

    // Get pT before reclustering
    double minScale = (minscaleIn > 0.) ? minscaleIn : scale;

    // If the maximal scale and the minimal scale coincide (as would
    // be the case for the corrected scales of unordered histories),
    // do not generate Sudakov
    if (minScale >= startingScale) break;

    // Find z and pT values at which the current state was formed, to
    // ensure that the showers can order the next emission correctly in
    // rapidity, if required.
    // NOT CORRECTLY SET FOR HIGHEST MULTIPLICITY STATE!
    double z = ( mergingHooksPtr->getNumberOfClusteringSteps(state) == 0
               || !mother )
             ? 0.5
             : mother->getCurrentZ(clusterIn.emittor,clusterIn.recoiler,
                 clusterIn.emitted);
    // Store z and pT values at which the current state was formed
    infoPtr->zNowISR(z);
    infoPtr->pT2NowISR(pow(startingScale,2));
    infoPtr->hasHistory(true);

    // Perform trial shower emission
    trial->next(process,event);
    // Get trial shower pT.
    double pTtrial   = trial->pTLastInShower();
    int typeTrial    = trial->typeLastInShower();

    // Clear parton systems.
    trial->resetTrial();

    // Get veto (merging) scale value
    double vetoScale  = (mother) ? 0. : mergingHooksPtr->tms();
    // Get merging scale in current event
    double tnow = mergingHooksPtr->tmsNow( event );

    // Done if evolution scale has fallen below minimum
    if ( pTtrial < minScale ) break;
    // Reset starting scale.
    startingScale = pTtrial;

    // Continue if this state is below the veto scale
    if ( tnow < vetoScale && vetoScale > 0. ) continue;

    // Retry if the trial emission was not allowed.
    if ( mergingHooksPtr->canVetoTrialEmission()
      && mergingHooksPtr->doVetoTrialEmission( process, event) ) continue;

    // Only consider allowed emissions for veto:
    // Only allow MPI for MPI no-emission probability.
    if ( type == -1 && typeTrial != 1 ) continue;
    // Only allow ISR or FSR for radiative no-emission probability.
    if ( type ==  1 && !(typeTrial == 2 || typeTrial >= 3) ) continue;

    // Veto event if trial pT was above the next nodal scale.
    if ( pTtrial > minScale ) doVeto = true;

    //// For last no-emission probability, veto event if the trial state
    //// is above the merging scale, i.e. in the matrix element region.
    //if ( !mother && tnow > vetoScale && vetoScale > 0. ) doVeto = true;

    // For 2 -> 2 pure QCD state, do not allow multiparton interactions
    // above the kinematical pT of the 2 -> 2 state.
    if ( type == -1
      && typeTrial == 1
      && mergingHooksPtr->getNumberOfClusteringSteps(process) == 0
      && ( mergingHooksPtr->getProcessString().compare("pp>jj") == 0
        || mergingHooksPtr->getProcessString().compare("pp>aj") == 0 )
      && pTtrial > hardFacScale(process) )
      return 0.0;

    // If pT of trial emission was in suitable range (trial shower
    // successful), return false
    if ( pTtrial < minScale ) doVeto = false;

    // Done
    break;

  }

  // Done
  return ( (doVeto) ? 0. : 1. );
}

//--------------------------------------------------------------------------

// Assume we have a vector of i elements containing indices into
// another vector with N elements. Update the indices so that all
// unique combinations (starting from 0,1,2,3, ...) are
// covered. Return false when all combinations have been ehausted.

bool History::updateind(vector<int> & ind, int i, int N) {
  if ( i < 0 ) return false;
  if ( ++ind[i] < N ) return true;
  if ( !updateind(ind, i - 1, N - 1) ) return false;
  ind[i] = ind[i - 1] + 1;
  return true;
}

//--------------------------------------------------------------------------

// Return the expansion of the no-emission probability up to the Nth
// term. Optionally calculate the the terms using fixed alphaS
// and/or PDF ratios.

vector<double>
History::countEmissions(PartonLevel* trial, double maxscale,
                        double minscale, int showerType, double as0,
                        AlphaStrong * asFSR, AlphaStrong * asISR, int N = 1,
                        bool fixpdf = true, bool fixas = true) {

  if ( N < 0 ) return vector<double>();
  vector<double> result(N+1);
  result[0] = 1.0;
  if ( N < 1 ) return result;

  // Copy state to local process
  Event process = state;

  double startingScale   = maxscale;
  // Careful when setting shower starting scale for pure QCD and prompt
  // photon case.
  if ( mergingHooksPtr->getNumberOfClusteringSteps(process) == 0
    && ( mergingHooksPtr->getProcessString().compare("pp>jj") == 0
      || mergingHooksPtr->getProcessString().compare("pp>aj") == 0 ) )
      startingScale = min( startingScale, hardFacScale(process) );

  vector<double> wts;

  while ( true ) {
    // Reset trialShower object
    trial->resetTrial();
    // Construct event to be showered
    Event event = Event();
    event.init("(hard process-modified)", particleDataPtr);
    event.clear();

    // Reset process scale
    process.scale(startingScale);

    // If the maximal scale and the minimal scale coincide (as would
    // be the case for the corrected scales of unordered histories),
    // do not generate Sudakov
    if (minscale >= startingScale) return result;

    // Find z and pT values at which the current state was formed, to
    // ensure that the showers can order the next emission correctly in
    // rapidity, if required
    if ( mother ) {
      double z = ( mergingHooksPtr->getNumberOfClusteringSteps(state) == 0)
               ? 0.5
               : mother->getCurrentZ(clusterIn.emittor,clusterIn.recoiler,
                   clusterIn.emitted);
      // Store z and pT values at which the current state was formed
      infoPtr->zNowISR(z);
      infoPtr->pT2NowISR(pow(startingScale,2));
      infoPtr->hasHistory(true);
    }

    // Perform trial shower emission
    trial->next(process,event);

    // Get trial shower pT
    double pTtrial = trial->pTLastInShower();
    int typeTrial  = trial->typeLastInShower();

    // Clear parton systems.
    trial->resetTrial();

    // Get veto (merging) scale value
    double vetoScale  = (mother) ? 0. : mergingHooksPtr->tms();
    // Get merging scale in current event
    double tnow = mergingHooksPtr->tmsNow( event );

    // Save scale of current state.
    startingScale   = pTtrial;
    // If the scale of the current state is below the minimal scale, exit.
    if ( pTtrial < minscale ) break;
    // If this state is below the merging scale, do not count emission.
    if ( tnow < vetoScale && vetoScale > 0. ) continue;
    // Retry if the trial emission was not allowed.
    if ( mergingHooksPtr->canVetoTrialEmission()
      && mergingHooksPtr->doVetoTrialEmission( process, event) ) continue;

    // Check if a new emission should be generated, either because
    // the latest emission was not of the desired kind or if the
    // emission was above the minimal scale
    double alphaSinPS = as0;
    double pdfs = 1.0;
    // Initial state splittings.
    if ( (showerType == -1 || showerType == 2) && typeTrial == 2 ) {
      // Get weight to translate to alpha_s at fixed renormalisation scale.
      if ( fixas ) alphaSinPS = (*asISR).alphaS(pTtrial*pTtrial);
      // Get weight to translate to PDFs at fixed factorisation scale.
      if ( fixpdf )
        pdfs = pdfFactor( event, typeTrial, pTtrial,
                          mergingHooksPtr->muFinME() );
    // Final state splittings.
    } else if ( (showerType == 1 || showerType == 2) && typeTrial >= 3 ) {
      // Get weight to translate to alpha_s at fixed renormalisation scale.
      if ( fixas ) alphaSinPS = (*asFSR).alphaS(pTtrial*pTtrial);
      // Get weight to translate to PDFs at fixed factorisation scale. Needed
      // for final state splittings with initial state recoiler.
      if ( fixpdf )
        pdfs = pdfFactor( event, typeTrial, pTtrial,
                          mergingHooksPtr->muFinME() );
    }

    // Save weight correcting to emission generated with fixed scales.
    if ( typeTrial == 2 || typeTrial >= 3 )
      wts.push_back(as0/alphaSinPS*pdfs);

  }

  for ( int n = 1; n <= min(N, int(wts.size())); ++n ) {
    vector<int> ind(N);
    for ( int i = 0; i < N; ++i ) ind[i] = i;
    do {
      double x = 1.0;
      for ( int j = 0; j < n; ++j ) x *= wts[ind[j]];
      result[n] += x;
    }  while ( updateind(ind, n - 1, wts.size()) );
    if ( n%2 ) result[n] *= -1.0;
  }

  // Done
  return result;
}

//--------------------------------------------------------------------------

// Function to integrate PDF ratios between two scales over x and t,
// where the PDFs are always evaluated at the lower t-integration limit

double History::monteCarloPDFratios(int flav, double x, double maxScale,
         double minScale, double pdfScale, double asME, Rndm* rndmPtr) {

  // Perform numerical integration for PDF ratios
  // Prefactor is as/2PI
  double factor = asME / (2.*M_PI);
  // Scale integration just produces a multiplicative logarithm
  factor *= log(maxScale/minScale);

  // For identical scales, done
  if (factor == 0.) return 0.;

  // Declare constants
  double CF = 4./3.;
  double CA = 3.;
  double NF = 4.;
  double TR = 1./2.;

  double integral = 0.;
  double RN = rndmPtr->flat();

  if (flav == 21) {
    double zTrial = pow(x,RN);
    integral  = -log(x) * zTrial *
                integrand(flav, x, pdfScale, zTrial);
    integral += 1./6.*(11.*CA - 4.*NF*TR)
              + 2.*CA*log(1.-x);
  } else {
    double zTrial = x + RN*(1. - x);
    integral  = (1.-x) *
                integrand(flav, x, pdfScale, zTrial);
    integral += 3./2.*CF
              + 2.*CF*log(1.-x);
  }

  // Done
  return (factor*integral);
}

/*--------------- METHODS USED FOR CONTRUCTION OF ALL HISTORIES --------- */

// Check if a ordered (and complete) path has been found in the
// initial node, in which case we will no longer be interested in
// any unordered paths.

bool History::onlyOrderedPaths() {
  if ( !mother || foundOrderedPath ) return foundOrderedPath;
  return  foundOrderedPath = mother->onlyOrderedPaths();
}

//--------------------------------------------------------------------------

// Check if a STRONGLY ordered (and complete) path has been found in the
// initial node, in which case we will no longer be interested in
// any unordered paths.

bool History::onlyStronglyOrderedPaths() {
  if ( !mother || foundStronglyOrderedPath ) return foundStronglyOrderedPath;
  return  foundStronglyOrderedPath = mother->onlyStronglyOrderedPaths();
}

//--------------------------------------------------------------------------

// Check if an allowed (according to user-criterion) path has been found in
// the initial node, in which case we will no longer be interested in
// any forbidden paths.

bool History::onlyAllowedPaths() {
  if ( !mother || foundAllowedPath ) return foundAllowedPath;
  return foundAllowedPath = mother->onlyAllowedPaths();
}

//--------------------------------------------------------------------------

// When a full path has been found, register it with the initial
// history node.
// IN  History : History to be registered as path
//     bool    : Specifying if clusterings so far were ordered
//     bool    : Specifying if path is complete down to 2->2 process
// OUT true if History object forms a plausible path (eg prob>0 ...)

bool History::registerPath(History & l, bool isOrdered,
       bool isStronglyOrdered, bool isAllowed, bool isComplete) {

  // We are not interested in improbable paths.
  if ( l.prob <= 0.0)
    return false;
  // We only register paths in the initial node.
  if ( mother ) return mother->registerPath(l, isOrdered,
                         isStronglyOrdered, isAllowed, isComplete);

  // Again, we are not interested in improbable paths.
  if ( sumpath == sumpath + l.prob )
    return false;
  if ( mergingHooksPtr->canCutOnRecState()
    && foundAllowedPath && !isAllowed )
    return false;
  if ( mergingHooksPtr->enforceStrongOrdering()
    && foundStronglyOrderedPath && !isStronglyOrdered )
    return false;
  if ( mergingHooksPtr->orderHistories()
    && foundOrderedPath && !isOrdered ) {
    // Prefer complete or allowed paths to ordered paths.
    if ( (!foundCompletePath && isComplete)
      || (!foundAllowedPath && isAllowed) ) ;
    else return false;
  }

  if ( foundCompletePath && !isComplete )
    return false;
  if ( !mergingHooksPtr->canCutOnRecState()
    && !mergingHooksPtr->allowCutOnRecState() )
    foundAllowedPath = true;

  if ( mergingHooksPtr->canCutOnRecState() && isAllowed && isComplete) {
    if ( !foundAllowedPath || !foundCompletePath ) {
      // If this is the first complete, allowed path, discard the
      // old, disallowed or incomplete ones.
      paths.clear();
      sumpath = 0.0;
    }
    foundAllowedPath = true;

  }

  if ( mergingHooksPtr->enforceStrongOrdering() && isStronglyOrdered
     && isComplete ) {
    if ( !foundStronglyOrderedPath || !foundCompletePath ) {
      // If this is the first complete, ordered path, discard the
      // old, non-ordered or incomplete ones.
      paths.clear();
      sumpath = 0.0;
    }
    foundStronglyOrderedPath = true;
    foundCompletePath = true;

  }

  if ( mergingHooksPtr->orderHistories() && isOrdered && isComplete ) {
    if ( !foundOrderedPath || !foundCompletePath ) {
      // If this is the first complete, ordered path, discard the
      // old, non-ordered or incomplete ones.
      paths.clear();
      sumpath = 0.0;
    }
    foundOrderedPath = true;
    foundCompletePath = true;

  }

  if ( isComplete ) {
    if ( !foundCompletePath ) {
      // If this is the first complete path, discard the old,
      // incomplete ones.
      paths.clear();
      sumpath = 0.0;
    }
    foundCompletePath = true;
  }

  // Remember, if this path is ordered, even if no ordering is required
  if ( isOrdered ) {
    foundOrderedPath = true;
  }

  // Index path by probability
  sumpath += l.prob;
  paths[sumpath] = &l;

  return true;
}

//--------------------------------------------------------------------------

// For the history-defining state (and if necessary interfering
// states), find all possible clusterings.
// NO INPUT
// OUT vector of all (rad,rec,emt) systems

vector<Clustering> History::getAllQCDClusterings() {
  vector<Clustering> ret;
  // Initialise vectors to keep track of position of partons in the
  // history-defining state
  vector <int> PosFinalPartn;
  vector <int> PosInitPartn;
  vector <int> PosFinalGluon;
  vector <int> PosFinalQuark;
  vector <int> PosFinalAntiq;
  vector <int> PosInitGluon;
  vector <int> PosInitQuark;
  vector <int> PosInitAntiq;

  // Search event record for final state particles and store these in
  // quark, anti-quark and gluon vectors
  for ( int i=0; i < state.size(); ++i )
    if ( state[i].isFinal() && state[i].colType() !=0 ) {
      // Store final partons
      if ( state[i].id() == 21 ) PosFinalGluon.push_back(i);
      else if ( state[i].idAbs() < 10 && state[i].id() > 0)
        PosFinalQuark.push_back(i);
      else if ( state[i].idAbs() < 10 && state[i].id() < 0)
        PosFinalAntiq.push_back(i);
    } else if (state[i].status() == -21 && state[i].colType() != 0 ) {
      // Store initial partons
      if ( state[i].id() == 21 ) PosInitGluon.push_back(i);
      else if ( state[i].idAbs() < 10 && state[i].id() > 0)
        PosInitQuark.push_back(i);
      else if ( state[i].idAbs() < 10 && state[i].id() < 0)
        PosInitAntiq.push_back(i);
    }

  // Get all clusterings for input state
  vector<Clustering> systems;
  systems = getQCDClusterings(state);
  ret.insert(ret.end(), systems.begin(), systems.end());
  systems.resize(0);

  // If valid clusterings were found, return
  if ( !ret.empty() ) return ret;
  // If no clusterings have been found until now, try to find
  // clusterings of diagrams that interfere with the current one
  // (i.e. change the colours of the current event slightly and run
  //  search again)
  else if ( ret.empty()
        && mergingHooksPtr->allowColourShuffling() ) {
    Event NewState = Event(state);
    // Start with changing final state quark colour
    for(int i = 0; i < int(PosFinalQuark.size()); ++i) {
      // Never change the hard process candidates
      if ( mergingHooksPtr->hardProcess.matchesAnyOutgoing(PosFinalQuark[i],
       NewState) )
        continue;
      int col = NewState[PosFinalQuark[i]].col();
      for(int j = 0; j < int(PosInitAntiq.size()); ++j) {
        // Now swap colours
        int acl = NewState[PosInitAntiq[j]].acol();
        if ( col == acl ) continue;
        NewState[PosFinalQuark[i]].col(acl);
        NewState[PosInitAntiq[j]].acol(col);
        systems = getQCDClusterings(NewState);
        if (!systems.empty()) {
          state = NewState;
          NewState.clear();
          ret.insert(ret.end(), systems.begin(), systems.end());
          systems.resize(0);
          return ret;
        }
      }
    }
    // Now change final state antiquark anticolour
    for(int i = 0; i < int(PosFinalAntiq.size()); ++i) {
      // Never change the hard process candidates
      if ( mergingHooksPtr->hardProcess.matchesAnyOutgoing(PosFinalAntiq[i],
       NewState) )
        continue;
      int acl = NewState[PosFinalAntiq[i]].acol();
      for(int j = 0; j < int(PosInitQuark.size()); ++j) {
        // Now swap colours
        int col = NewState[PosInitQuark[j]].col();
        if ( col == acl ) continue;
        NewState[PosFinalAntiq[i]].acol(col);
        NewState[PosInitQuark[j]].col(acl);
        systems = getQCDClusterings(NewState);
        if (!systems.empty()) {
          state = NewState;
          NewState.clear();
          ret.insert(ret.end(), systems.begin(), systems.end());
          systems.resize(0);
          return ret;
        }
      }
    }

    if ( !ret.empty() ) {
      string message="Warning in History::getAllQCDClusterings: Changed";
      message+=" colour structure to allow at least one clustering.";
      infoPtr->errorMsg(message);
    }

  }

  // Done
  return ret;
}

//--------------------------------------------------------------------------

// For one given state, find all possible clusterings.
// IN  Event : state to be investigated
// OUT vector of all (rad,rec,emt) systems in the state

vector<Clustering> History::getQCDClusterings( const Event& event) {
  vector<Clustering> ret;

  // Initialise vectors to keep track of position of partons in the
  // input event
  vector <int> PosFinalPartn;
  vector <int> PosInitPartn;

  vector <int> PosFinalGluon;
  vector <int> PosFinalQuark;
  vector <int> PosFinalAntiq;
  vector <int> PosInitGluon;
  vector <int> PosInitQuark;
  vector <int> PosInitAntiq;

  // Search event record for final state particles and store these in
  // quark, anti-quark and gluon vectors
  for (int i=0; i < event.size(); ++i)
    if ( event[i].isFinal() && event[i].colType() !=0 ) {
      // Store final partons
      PosFinalPartn.push_back(i);
      if ( event[i].id() == 21 ) PosFinalGluon.push_back(i);
      else if ( event[i].idAbs() < 10 && event[i].id() > 0)
        PosFinalQuark.push_back(i);
      else if ( event[i].idAbs() < 10 && event[i].id() < 0)
        PosFinalAntiq.push_back(i);
    } else if ( event[i].status() == -21 && event[i].colType() != 0 ) {
      // Store initial partons
      PosInitPartn.push_back(i);
      if ( event[i].id() == 21 ) PosInitGluon.push_back(i);
      else if ( event[i].idAbs() < 10 && event[i].id() > 0)
        PosInitQuark.push_back(i);
      else if ( event[i].idAbs() < 10 && event[i].id() < 0)
        PosInitAntiq.push_back(i);
    }

  int nFiGluon = int(PosFinalGluon.size());
  int nFiQuark = int(PosFinalQuark.size());
  int nFiAntiq = int(PosFinalAntiq.size());
  int nInGluon = int(PosInitGluon.size());
  int nInQuark = int(PosInitQuark.size());
  int nInAntiq = int(PosInitAntiq.size());

  vector<Clustering> systems;

  // Find rad + emt + rec systems:
  // (1) Start from gluon and find all (rad,rec,emt=gluon) triples
  for (int i = 0; i < nFiGluon; ++i) {
    int EmtGluon = PosFinalGluon[i];
    systems = findQCDTriple( EmtGluon, 2, event, PosFinalPartn, PosInitPartn);
    ret.insert(ret.end(), systems.begin(), systems.end());
    systems.resize(0);
  }

  // For more than one quark-antiquark pair in final state, check for
  // g -> qqbar splittings
  bool check_g2qq = true;
  if ( ( ( nInQuark + nInAntiq == 0 )
          && (nInGluon == 0)
          && (nFiQuark == 1) && (nFiAntiq == 1) )
    || ( ( nFiQuark + nFiAntiq == 0)
          && (nInQuark == 1) && (nInAntiq == 1) ) )
    check_g2qq = false;

  if ( check_g2qq ) {

    // (2) Start from quark and find all (rad,rec,emt=quark) triples
    //     ( when g -> q qbar occured )
    for( int i=0; i < nFiQuark; ++i) {
      int EmtQuark = PosFinalQuark[i];
      systems = findQCDTriple( EmtQuark,1,event, PosFinalPartn, PosInitPartn);
      ret.insert(ret.end(), systems.begin(), systems.end());
      systems.resize(0);
    }

    // (3) Start from anti-quark and find all (rad,rec,emt=anti-quark)
    //     triples ( when g -> q qbar occured )
    for( int i=0; i < nFiAntiq; ++i) {
      int EmtAntiq = PosFinalAntiq[i];
      systems = findQCDTriple( EmtAntiq,1,event, PosFinalPartn, PosInitPartn);
      ret.insert(ret.end(), systems.begin(), systems.end());
      systems.resize(0);
    }
  }

  return ret;
}

//--------------------------------------------------------------------------

// Function to construct (rad,rec,emt) triples from the event
// IN  int   : Position of Emitted in event record for which
//             dipoles should be constructed
//     int   : Colour topogy to be tested
//             1= g -> qqbar, causing 2 -> 2 dipole splitting
//             2= q(bar) -> q(bar) g && g -> gg,
//              causing a 2 -> 3 dipole splitting
//     Event : event record to be checked for ptential partners
// OUT vector of all allowed radiator+recoiler+emitted triples

vector<Clustering> History::findQCDTriple (int EmtTagIn, int colTopIn,
                      const Event& event,
                      vector<int> PosFinalPartn,
                      vector <int> PosInitPartn ) {

  // Copy input parton tag
  int EmtTag = EmtTagIn;
  // Copy input colour topology tag
  // (1: g --> qqbar splitting present, 2:rest)
  int colTop = colTopIn;

  // Initialise FinalSize
  int FinalSize = int(PosFinalPartn.size());
  int InitSize = int(PosInitPartn.size());
  int Size = InitSize + FinalSize;

  vector<Clustering> clus;

  // Search final partons to find partons colour-connected to
  // event[EmtTag], choose radiator, then choose recoiler
  for ( int a = 0; a < Size; ++a ) {
    int i    = (a < FinalSize)? a : (a - FinalSize) ;
    int iRad = (a < FinalSize)? PosFinalPartn[i] : PosInitPartn[i];

    if ( event[iRad].col() == event[EmtTag].col()
      && event[iRad].acol() == event[EmtTag].acol() )
      continue;

    if (iRad != EmtTag ) {
      int pTdef = event[iRad].isFinal() ? 1 : -1;
      int sign = (a < FinalSize)? 1 : -1 ;

      // First colour topology: g --> qqbar. Here, emt & rad should
      // have same flavour (causes problems for gamma->qqbar).
      if (colTop == 1) {

        if ( event[iRad].id() == -sign*event[EmtTag].id() ) {
          int col = -1;
          int acl = -1;
          if (event[iRad].id() < 0) {
            col = event[EmtTag].acol();
            acl = event[iRad].acol();
          } else {
             col = event[EmtTag].col();
             acl = event[iRad].col();
          }
          // Recoiler
          int iRec     = 0;
          // Colour partner
          int iPartner = 0;

          if (col > 0) {
            // Find recoiler by colour
            iRec = FindCol(col,iRad,EmtTag,event,1,true);
            // In initial state splitting has final state colour partner,
            // Save both partner and recoiler
            if ( (sign < 0) && (event[iRec].isFinal()) ) {
              // Save colour recoiler
              iPartner = iRec;
              // Reset kinematic recoiler to initial state parton
              for(int l = 0; l < int(PosInitPartn.size()); ++l)
                if (PosInitPartn[l] != iRad) iRec = PosInitPartn[l];
            // For final state splittings, colour partner and recoiler are
            // identical
            } else {
              iPartner = iRec;
            }
            if ( iRec != 0 && iPartner != 0
             && allowedClustering( iRad, EmtTag, iRec, iPartner, event) ) {
              clus.push_back( Clustering(EmtTag, iRad, iRec, iPartner,
                   pTLund(event[iRad], event[EmtTag], event[iRec], pTdef) ));
              continue;
            }

            // Reset partner
            iPartner = 0;
            // Find recoiler by colour
            iRec = FindCol(col,iRad,EmtTag,event,2,true);
            // In initial state splitting has final state colour partner,
            // Save both partner and recoiler
            if ( (sign < 0) && (event[iRec].isFinal()) ) {
              // Save colour recoiler
              iPartner = iRec;
              // Reset kinematic recoiler to initial state parton
              for(int l = 0; l < int(PosInitPartn.size()); ++l)
                if (PosInitPartn[l] != iRad) iRec = PosInitPartn[l];
            // For final state splittings, colour partner and recoiler are
            // identical
            } else {
              iPartner = iRec;
            }
            if ( iRec != 0 && iPartner != 0
             && allowedClustering( iRad, EmtTag, iRec, iPartner, event) ) {
              clus.push_back( Clustering(EmtTag, iRad, iRec, iPartner,
                   pTLund(event[iRad], event[EmtTag], event[iRec], pTdef) ));
              continue;
            }
          }


          if (acl > 0) {

            // Reset partner
            iPartner = 0;
            // Find recoiler by colour
            iRec = FindCol(acl,iRad,EmtTag,event,1,true);
            // In initial state splitting has final state colour partner,
            // Save both partner and recoiler
            if ( (sign < 0) && (event[iRec].isFinal()) ) {
              // Save colour recoiler
              iPartner = iRec;
              // Reset kinematic recoiler to initial state parton
              for(int l = 0; l < int(PosInitPartn.size()); ++l)
                if (PosInitPartn[l] != iRad) iRec = PosInitPartn[l];
            // For final state splittings, colour partner and recoiler are
            // identical
            } else {
              iPartner = iRec;
            }
            if ( iRec != 0 && iPartner != 0
             && allowedClustering( iRad, EmtTag, iRec, iPartner, event) ) {
              clus.push_back( Clustering(EmtTag, iRad, iRec, iPartner,
                   pTLund(event[iRad], event[EmtTag], event[iRec], pTdef) ));
              continue;
            }

            // Reset partner
            iPartner = 0;
            // Find recoiler by colour
            iRec = FindCol(acl,iRad,EmtTag,event,2,true);
            // In initial state splitting has final state colour partner,
            // Save both partner and recoiler
            if ( (sign < 0) && (event[iRec].isFinal()) ) {
              // Save colour recoiler
              iPartner = iRec;
              // Reset kinematic recoiler to initial state parton
              for(int l = 0; l < int(PosInitPartn.size()); ++l)
                if (PosInitPartn[l] != iRad) iRec = PosInitPartn[l];
            // For final state splittings, colour partner and recoiler are
            // identical
            } else {
              iPartner = iRec;
            }
            if ( iRec != 0 && iPartner != 0
             && allowedClustering( iRad, EmtTag, iRec, iPartner, event) ) {
              clus.push_back( Clustering(EmtTag, iRad, iRec, iPartner,
                   pTLund(event[iRad], event[EmtTag], event[iRec], pTdef) ));
              continue;
            }
          }
        // Initial gluon splitting
        } else if ( event[iRad].id() == 21
                  &&(  event[iRad].col() == event[EmtTag].col()
                    || event[iRad].acol() == event[EmtTag].acol() )) {
          // For an initial state radiator, always set recoiler
          // to the other initial state parton (recoil is taken
          // by full remaining system, so this is just a
          // labelling for such a process)
          int RecInit  = 0;
          for(int l = 0; l < int(PosInitPartn.size()); ++l)
            if (PosInitPartn[l] != iRad) RecInit = PosInitPartn[l];

          // Find the colour connected partner
          // Find colour index of radiator before splitting
          int col = getRadBeforeCol(iRad, EmtTag, event);
          int acl = getRadBeforeAcol(iRad, EmtTag, event);

          // Find the correct partner: If a colour line has split,
          // the partner is connected to the radiator before the splitting
          // by a colour line (same reasoning for anticolour). The colour
          // that split is the colour appearing twice in the
          // radiator + emitted pair.
          // Thus, if we remove a colour index with the clustering,
          // we should look for a colour partner, else look for
          // an anticolour partner
          int colRemove = (event[iRad].col() == event[EmtTag].col())
                  ? event[iRad].col() : 0;

          int iPartner = 0;
          if (colRemove > 0 && col > 0)
            iPartner = FindCol(col,iRad,EmtTag,event,1,true)
                     + FindCol(col,iRad,EmtTag,event,2,true);
          else if (colRemove > 0 && acl > 0)
            iPartner = FindCol(acl,iRad,EmtTag,event,1,true)
                     + FindCol(acl,iRad,EmtTag,event,2,true);

          if ( allowedClustering( iRad, EmtTag, RecInit, iPartner, event ) ) {
            clus.push_back( Clustering(EmtTag, iRad, RecInit, iPartner,
                 pTLund(event[iRad],event[EmtTag],event[RecInit], pTdef) ));
              continue;
          }
        }

      // Second colour topology: Gluon emission

      } else {
        if ( (event[iRad].col() == event[EmtTag].acol())
           || (event[iRad].acol() == event[EmtTag].col())
           || (event[iRad].col() == event[EmtTag].col())
           || (event[iRad].acol() == event[EmtTag].acol()) ) {
          // For the rest, choose recoiler to have a common colour
          // tag with radiator, while not being the "Emitted"

          int col = -1;
          int acl = -1;

          if (event[iRad].isFinal() ) {

            if ( event[iRad].id() < 0) {
              acl = event[EmtTag].acol();
              col = event[iRad].col();
            } else if ( event[iRad].id() > 0 && event[iRad].id() < 10) {
              col = event[EmtTag].col();
              acl = event[iRad].acol();
            } else {
              col = event[iRad].col();
              acl = event[iRad].acol();
            }

            int iRec = 0;
            if (col > 0) {
              iRec = FindCol(col,iRad,EmtTag,event,1,true);
              if ( (sign < 0) && (event[iRec].isFinal()) ) iRec = 0;
              if (iRec != 0
               && allowedClustering( iRad, EmtTag, iRec, iRec, event) ) {
                clus.push_back( Clustering(EmtTag, iRad, iRec, iRec,
                     pTLund(event[iRad],event[EmtTag],event[iRec], pTdef) ));
                continue;
              }

              iRec = FindCol(col,iRad,EmtTag,event,2,true);
              if ( (sign < 0) && (event[iRec].isFinal()) ) iRec = 0;
              if (iRec != 0
               && allowedClustering( iRad, EmtTag, iRec, iRec, event) ) {
                clus.push_back( Clustering(EmtTag, iRad, iRec, iRec,
                     pTLund(event[iRad],event[EmtTag],event[iRec], pTdef) ));
                continue;
              }
            }


            if (acl > 0) {
              iRec = FindCol(acl,iRad,EmtTag,event,1,true);
              if ( (sign < 0) && (event[iRec].isFinal()) ) iRec = 0;
              if (iRec != 0
               && allowedClustering( iRad, EmtTag, iRec, iRec, event) ) {
                clus.push_back( Clustering(EmtTag, iRad, iRec, iRec,
                     pTLund(event[iRad],event[EmtTag],event[iRec], pTdef) ));
                continue;
              }

              iRec = FindCol(acl,iRad,EmtTag,event,2,true);
              if ( (sign < 0) && (event[iRec].isFinal()) ) iRec = 0;
              if (iRec != 0
               && allowedClustering( iRad, EmtTag, iRec, iRec, event) ) {
                clus.push_back( Clustering(EmtTag, iRad, iRec, iRec,
                     pTLund(event[iRad],event[EmtTag],event[iRec], pTdef) ));
                continue;
              }
            }

          } else {

            // For an initial state radiator, always set recoiler
            // to the other initial state parton (recoil is taken

            // by full remaining system, so this is just a
            // labelling for such a process)
            int RecInit = 0;
            int iPartner = 0;
            for(int l = 0; l < int(PosInitPartn.size()); ++l)
              if (PosInitPartn[l] != iRad) RecInit = PosInitPartn[l];

            // Find the colour connected partner
            // Find colour index of radiator before splitting
            col = getRadBeforeCol(iRad, EmtTag, event);
            acl = getRadBeforeAcol(iRad, EmtTag, event);

            // Find the correct partner: If a colour line has split,
            // the partner is connected to the radiator before the splitting
            // by a colour line (same reasoning for anticolour). The colour
            // that split is the colour appearing twice in the
            // radiator + emitted pair.
            // Thus, if we remove a colour index with the clustering,
            // we should look for a colour partner, else look for
            // an anticolour partner
            int colRemove = (event[iRad].col() == event[EmtTag].col())
                    ? event[iRad].col() : 0;
            iPartner = (colRemove > 0)
                     ?   FindCol(col,iRad,EmtTag,event,1,true)
                       + FindCol(col,iRad,EmtTag,event,2,true)
                     :   FindCol(acl,iRad,EmtTag,event,1,true)
                       + FindCol(acl,iRad,EmtTag,event,2,true);

            if ( allowedClustering( iRad, EmtTag, RecInit, iPartner, event)) {
              clus.push_back( Clustering(EmtTag, iRad, RecInit, iPartner,
                   pTLund(event[iRad],event[EmtTag],event[RecInit], pTdef)));

              continue;
            }
          }
        }
      }
    }
  }

  // Done
  return clus;
}

//--------------------------------------------------------------------------

// For the history-defining state (and if necessary interfering
// states), find all possible clusterings.
// NO INPUT
// OUT vector of all (rad,rec,emt) systems

vector<Clustering> History::getAllEWClusterings() {
  vector<Clustering> ret;

  // Get all clusterings for input state
  vector<Clustering> systems;
  systems = getEWClusterings(state);
  ret.insert(ret.end(), systems.begin(), systems.end());
  // Done
  return ret;
}

//--------------------------------------------------------------------------

// For one given state, find all possible clusterings.
// IN  Event : state to be investigated
// OUT vector of all (rad,rec,emt) systems in the state

vector<Clustering> History::getEWClusterings( const Event& event) {
  vector<Clustering> ret;

  // Initialise vectors to keep track of position of partons in the
  // input event
  vector <int> PosFinalPartn;
  vector <int> PosInitPartn;
  vector <int> PosFinalW;

  // Search event record for final state particles and store these in
  // quark, anti-quark and gluon vectors
  for ( int i=0; i < event.size(); ++i )
    if ( event[i].isFinal() && abs(event[i].colType()) == 1 ) {
      // Store final partons
      PosFinalPartn.push_back(i);
    } else if ( event[i].status() == -21 && abs(event[i].colType()) == 1 ) {
      // Store initial partons
      PosInitPartn.push_back(i);
    }
  // Search event record for final W
  for ( int i=0; i < event.size(); ++i )
    if ( event[i].isFinal() && event[i].idAbs() == 24 )
      PosFinalW.push_back( i );

  vector<Clustering> systems;
  // Find rad + emt + rec systems:
  // (1) Start from W boson and find all (rad,rec,emt=W) triples
  for ( int i = 0; i <  int(PosFinalW.size()); ++i ) {
    int EmtW = PosFinalW[i];
    systems = findEWTriple( EmtW, event, PosFinalPartn);
    ret.insert(ret.end(), systems.begin(), systems.end());
    systems.resize(0);
  }

  return ret;
}

//--------------------------------------------------------------------------

// Function to construct (rad,rec,emt) triples from the event
// IN  int   : Position of Emitted in event record for which
//             dipoles should be constructed
//     int   : Colour topogy to be tested
//             1= g -> qqbar, causing 2 -> 2 dipole splitting
//             2= q(bar) -> q(bar) g && g -> gg,
//              causing a 2 -> 3 dipole splitting
//     Event : event record to be checked for ptential partners
// OUT vector of all allowed radiator+recoiler+emitted triples

vector<Clustering> History::findEWTriple ( int EmtTagIn, const Event& event,
                      vector<int> PosFinalPartn ) {
  // Copy input parton tag
  int EmtTag = EmtTagIn;
  // Copy input colour topology tag
  // (1: g --> qqbar splitting present, 2:rest)

  // Initialise FinalSize
  int FinalSize = int(PosFinalPartn.size());

  vector<Clustering> clus;

  // Search final partons to find partons colour-connected to
  // event[EmtTag], choose radiator, then choose recoiler
  for ( int a = 0; a < FinalSize; ++a ) {

    int iRad = PosFinalPartn[a];
    if (iRad != EmtTag ) {
      int pTdef = 1;
      // Find recoiler by flavour.
      int flavRad = event[iRad].id();
      int flavEmt = event[EmtTag].id();

      // Loop through final partons and try to find matching flavours.
      for ( int i = 0; i < FinalSize; ++i ) {
        int iRec = PosFinalPartn[i];
        if ( i != a && flavEmt > 0
          && event[iRec].id() == -flavRad - 1 )
          clus.push_back( Clustering(EmtTag, iRad, iRec, iRec,
               pTLund(event[iRad],event[EmtTag],event[iRec], pTdef) ) );
      }
    }
  }

  // Done
  return clus;
}

//--------------------------------------------------------------------------

// For the history-defining state (and if necessary interfering
// states), find all possible clusterings.
// NO INPUT
// OUT vector of all (rad,rec,emt) systems

vector<Clustering> History::getAllSQCDClusterings() {
  vector<Clustering> ret;

  // Get all clusterings for input state
  vector<Clustering> systems;
  systems = getSQCDClusterings(state);
  ret.insert(ret.end(), systems.begin(), systems.end());
  // Done
  return ret;
}

//--------------------------------------------------------------------------

// For one given state, find all possible clusterings.
// IN  Event : state to be investigated
// OUT vector of all (rad,rec,emt) systems in the state

vector<Clustering> History::getSQCDClusterings( const Event& event) {
  vector<Clustering> ret;

  // Initialise vectors to keep track of position of partons in the
  // input event
  vector <int> PosFinalPartn;
  vector <int> PosInitPartn;

  vector <int> PosFinalGluon;
  vector <int> PosFinalQuark;
  vector <int> PosFinalAntiq;
  vector <int> PosInitGluon;
  vector <int> PosInitQuark;
  vector <int> PosInitAntiq;

  // Search event record for final state particles and store these in
  // quark, anti-quark and gluon vectors
  for (int i=0; i < event.size(); ++i)
    if ( event[i].isFinal() && event[i].colType() !=0 ) {
      // Store final partons
      PosFinalPartn.push_back(i);
      if ( event[i].id() == 21 || event[i].id() == 1000021)
        PosFinalGluon.push_back(i);
      else if ( (event[i].idAbs() < 10 && event[i].id() > 0)
             || (event[i].idAbs() < 1000010 && event[i].idAbs() > 1000000
             && event[i].id() > 0)
             || (event[i].idAbs() < 2000010 && event[i].idAbs() > 2000000
             && event[i].id() > 0))
        PosFinalQuark.push_back(i);
      else if ( (event[i].idAbs() < 10 && event[i].id() < 0)
             || (event[i].idAbs() < 1000010 && event[i].idAbs() > 1000000
             && event[i].id() < 0)
             || (event[i].idAbs() < 2000010 && event[i].idAbs() > 2000000
             && event[i].id() < 0))
        PosFinalAntiq.push_back(i);
    } else if ( event[i].status() == -21 && event[i].colType() != 0 ) {
      // Store initial partons
      PosInitPartn.push_back(i);
      if ( event[i].id() == 21 || event[i].id() == 1000021)
        PosInitGluon.push_back(i);
      else if ( (event[i].idAbs() < 10 && event[i].id() > 0)
             || (event[i].idAbs() < 1000010 && event[i].idAbs() > 1000000
             && event[i].id() > 0)
             || (event[i].idAbs() < 2000010 && event[i].idAbs() > 2000000
             && event[i].id() > 0))
        PosInitQuark.push_back(i);
      else if ( (event[i].idAbs() < 10 && event[i].id() < 0)
             || (event[i].idAbs() < 1000010 && event[i].idAbs() > 1000000
             && event[i].id() < 0)
             || (event[i].idAbs() < 2000010 && event[i].idAbs() > 2000000
             && event[i].id() < 0))
        PosInitAntiq.push_back(i);
    }

  int nFiGluon = int(PosFinalGluon.size());
  int nFiQuark = int(PosFinalQuark.size());
  int nFiAntiq = int(PosFinalAntiq.size());
  int nInGluon = int(PosInitGluon.size());
  int nInQuark = int(PosInitQuark.size());
  int nInAntiq = int(PosInitAntiq.size());
  vector<Clustering> systems;

  // Find rad + emt + rec systems:
  // (1) Start from gluon and find all (rad,rec,emt=gluon) triples
  for (int i = 0; i < nFiGluon; ++i) {
    int EmtGluon = PosFinalGluon[i];
    systems = findSQCDTriple( EmtGluon, 2, event, PosFinalPartn, PosInitPartn);
    ret.insert(ret.end(), systems.begin(), systems.end());
    systems.resize(0);
  }

  // For more than one quark-antiquark pair in final state, check for
  // g -> qqbar splittings
  bool check_g2qq = true;
  if ( ( ( nInQuark + nInAntiq == 0 )
          && (nInGluon == 0)
          && (nFiQuark == 1) && (nFiAntiq == 1) )
    || ( ( nFiQuark + nFiAntiq == 0)
          && (nInQuark == 1) && (nInAntiq == 1) ) )
    check_g2qq = false;

  if ( check_g2qq ) {

    // (2) Start from quark and find all (rad,rec,emt=quark) triples
    //     ( when g -> q qbar occured )
    for( int i=0; i < nFiQuark; ++i) {
      int EmtQuark = PosFinalQuark[i];
      systems = findSQCDTriple( EmtQuark,1,event, PosFinalPartn, PosInitPartn);
      ret.insert(ret.end(), systems.begin(), systems.end());
      systems.resize(0);
    }

    // (3) Start from anti-quark and find all (rad,rec,emt=anti-quark)
    //     triples ( when g -> q qbar occured )
    for( int i=0; i < nFiAntiq; ++i) {
      int EmtAntiq = PosFinalAntiq[i];
      systems = findSQCDTriple( EmtAntiq,1,event, PosFinalPartn, PosInitPartn);
      ret.insert(ret.end(), systems.begin(), systems.end());
      systems.resize(0);
    }

  }

  return ret;
}

//--------------------------------------------------------------------------

// Function to construct (rad,rec,emt) triples from the event
// IN  int   : Position of Emitted in event record for which
//             dipoles should be constructed
//     int   : Colour topogy to be tested
//             1= g -> qqbar, causing 2 -> 2 dipole splitting
//             2= q(bar) -> q(bar) g && g -> gg,
//              causing a 2 -> 3 dipole splitting
//     Event : event record to be checked for ptential partners
// OUT vector of all allowed radiator+recoiler+emitted triples

vector<Clustering> History::findSQCDTriple (int EmtTagIn, int colTopIn,
                      const Event& event,
                      vector<int> PosFinalPartn,
                      vector <int> PosInitPartn ) {

  // Copy input parton tag
  int EmtTag = EmtTagIn;
  // Copy input colour topology tag
  // (1: g --> qqbar splitting present, 2:rest)
  int colTop = colTopIn;

  // PDG numbering offset for squarks
  int offsetL = 1000000;
  int offsetR = 2000000;

  // Initialise FinalSize
  int FinalSize = int(PosFinalPartn.size());
  int InitSize = int(PosInitPartn.size());
  int Size = InitSize + FinalSize;

  vector<Clustering> clus;

  // Search final partons to find partons colour-connected to
  // event[EmtTag], choose radiator, then choose recoiler
  for ( int a = 0; a < Size; ++a ) {
    int i    = (a < FinalSize)? a : (a - FinalSize) ;
    int iRad = (a < FinalSize)? PosFinalPartn[i] : PosInitPartn[i];

    if ( event[iRad].col() == event[EmtTag].col()
      && event[iRad].acol() == event[EmtTag].acol() )
      continue;

    // Save radiator flavour.
    int radID = event[iRad].id();
    // Remember if radiator is BSM.
    bool isSQCDrad = (abs(radID) > offsetL);
    // Remember if emitted is BSM.
    bool isSQCDemt = (event[EmtTag].idAbs() > offsetL );

    if (iRad != EmtTag ) {
      int pTdef = event[iRad].isFinal() ? 1 : -1;
      int sign = (a < FinalSize)? 1 : -1 ;

      // Disalllow clusterings resulting in an initial state sQCD parton!
      int radBefID = getRadBeforeFlav(iRad,EmtTag,event);
      if ( pTdef == -1 && abs(radBefID) > offsetL ) continue;

      // First colour topology: g --> qqbar. Here, emt & rad should
      // have same flavour (causes problems for gamma->qqbar).
      if (colTop == 1) {

        int radSign = (event[iRad].id() < 0) ? -1 : 1;
        int emtSign = (event[EmtTag].id() < 0) ? -1 : 1;

        // Final gluino splitting.
        bool finalSplitting = false;
        if ( abs(radID) < 10
          && radSign*(abs(radID)+offsetL) == -sign*event[EmtTag].id() )
          finalSplitting = true;
        if ( abs(radID) < 10
          && radSign*(abs(radID)+offsetR) == -sign*event[EmtTag].id() )
          finalSplitting = true;
        if ( abs(radID) > offsetL && abs(radID) < offsetL+10
          && radID == -sign*emtSign*( event[EmtTag].idAbs() + offsetL) )
          finalSplitting = true;
        if ( abs(radID) > offsetR && abs(radID) < offsetR+10
          && radID == -sign*emtSign*( event[EmtTag].idAbs() + offsetR) )
          finalSplitting = true;

        // Initial gluon splitting.
        bool initialSplitting = false;
        if ( radID == 21 && ( ( event[EmtTag].idAbs() > offsetL
                             && event[EmtTag].idAbs() < offsetL+10)
                           || ( event[EmtTag].idAbs() > offsetR
                             && event[EmtTag].idAbs() < offsetR+10) )
          && ( event[iRad].col() == event[EmtTag].col()
            || event[iRad].acol() == event[EmtTag].acol() ) )
          initialSplitting = true;

        if ( finalSplitting ) {

          int col = -1;
          int acl = -1;
          if ( radID < 0 && event[iRad].colType() == -1) {
            acl = event[EmtTag].acol();
            col = event[iRad].acol();
          } else if ( event[iRad].colType() == 1 ) {
            col = event[EmtTag].col();
            acl = event[iRad].col();
          }

          // Recoiler
          int iRec     = 0;
          // Colour partner
          int iPartner = 0;

          if (col > 0) {
            // Find recoiler by colour
            iRec = FindCol(col,iRad,EmtTag,event,1,true);
            // In initial state splitting has final state colour partner,
            // Save both partner and recoiler
            if ( (sign < 0) && (event[iRec].isFinal()) ) {
              // Save colour recoiler
              iPartner = iRec;
              // Reset kinematic recoiler to initial state parton
              for(int l = 0; l < int(PosInitPartn.size()); ++l)
                if (PosInitPartn[l] != iRad) iRec = PosInitPartn[l];
            // For final state splittings, colour partner and recoiler are
            // identical
            } else {
              iPartner = iRec;
            }

            // Not interested in pure QCD triples here.
            if ( !isSQCDrad && !isSQCDemt
              && (event[iRec].idAbs() < 10 || event[iRec].id() == 21) )
              iRec = 0;

            if ( iRec != 0 && iPartner != 0
             && allowedClustering( iRad, EmtTag, iRec, iPartner, event) ) {
              clus.push_back( Clustering(EmtTag, iRad, iRec, iPartner,
                   pTLund(event[iRad], event[EmtTag], event[iRec], pTdef) ));
              continue;
            }

            // Reset partner
            iPartner = 0;
            // Find recoiler by colour
            iRec = FindCol(col,iRad,EmtTag,event,2,true);
            // In initial state splitting has final state colour partner,
            // Save both partner and recoiler
            if ( (sign < 0) && (event[iRec].isFinal()) ) {
              // Save colour recoiler
              iPartner = iRec;
              // Reset kinematic recoiler to initial state parton
              for(int l = 0; l < int(PosInitPartn.size()); ++l)
                if (PosInitPartn[l] != iRad) iRec = PosInitPartn[l];
            // For final state splittings, colour partner and recoiler are
            // identical
            } else {
              iPartner = iRec;
            }

            // Not interested in pure QCD triples here.
            if ( !isSQCDrad && !isSQCDemt
              && (event[iRec].idAbs() < 10 || event[iRec].id() == 21) )
              iRec = 0;

            if ( iRec != 0 && iPartner != 0
             && allowedClustering( iRad, EmtTag, iRec, iPartner, event) ) {
              clus.push_back( Clustering(EmtTag, iRad, iRec, iPartner,
                   pTLund(event[iRad], event[EmtTag], event[iRec], pTdef) ));
              continue;
            }
          }

          if (acl > 0) {

            // Reset partner
            iPartner = 0;
            // Find recoiler by colour
            iRec = FindCol(acl,iRad,EmtTag,event,1,true);
            // In initial state splitting has final state colour partner,
            // Save both partner and recoiler
            if ( (sign < 0) && (event[iRec].isFinal()) ) {
              // Save colour recoiler
              iPartner = iRec;
              // Reset kinematic recoiler to initial state parton
              for(int l = 0; l < int(PosInitPartn.size()); ++l)
                if (PosInitPartn[l] != iRad) iRec = PosInitPartn[l];
            // For final state splittings, colour partner and recoiler are
            // identical
            } else {
              iPartner = iRec;
            }

            // Not interested in pure QCD triples here.
            if ( !isSQCDrad && !isSQCDemt
              && (event[iRec].idAbs() < 10 || event[iRec].id() == 21) )
              iRec = 0;

            if ( iRec != 0 && iPartner != 0
             && allowedClustering( iRad, EmtTag, iRec, iPartner, event) ) {
              clus.push_back( Clustering(EmtTag, iRad, iRec, iPartner,
                   pTLund(event[iRad], event[EmtTag], event[iRec], pTdef) ));
              continue;
            }

            // Reset partner
            iPartner = 0;
            // Find recoiler by colour
            iRec = FindCol(acl,iRad,EmtTag,event,2,true);
            // In initial state splitting has final state colour partner,
            // Save both partner and recoiler
            if ( (sign < 0) && (event[iRec].isFinal()) ) {
              // Save colour recoiler
              iPartner = iRec;
              // Reset kinematic recoiler to initial state parton
              for(int l = 0; l < int(PosInitPartn.size()); ++l)
                if (PosInitPartn[l] != iRad) iRec = PosInitPartn[l];
            // For final state splittings, colour partner and recoiler are
            // identical
            } else {
              iPartner = iRec;
            }

            // Not interested in pure QCD triples here.
            if ( !isSQCDrad && !isSQCDemt
              && (event[iRec].idAbs() < 10 || event[iRec].id() == 21) )
              iRec = 0;

            if ( iRec != 0 && iPartner != 0
             && allowedClustering( iRad, EmtTag, iRec, iPartner, event) ) {
              clus.push_back( Clustering(EmtTag, iRad, iRec, iPartner,
                   pTLund(event[iRad], event[EmtTag], event[iRec], pTdef) ));
              continue;
            }
          }
        // Initial gluon splitting
        } else if ( initialSplitting ) {

          // SM splittings already taken care of in findQCDTriple.
          if ( !isSQCDrad && !isSQCDemt ) continue;

          // For an initial state radiator, always set recoiler
          // to the other initial state parton (recoil is taken
          // by full remaining system, so this is just a
          // labelling for such a process)
          int RecInit  = 0;
          for(int l = 0; l < int(PosInitPartn.size()); ++l)
            if (PosInitPartn[l] != iRad) RecInit = PosInitPartn[l];

          // Find the colour connected partner
          // Find colour index of radiator before splitting
          int col = getRadBeforeCol(iRad, EmtTag, event);
          int acl = getRadBeforeAcol(iRad, EmtTag, event);

          // Find the correct partner: If a colour line has split,
          // the partner is connected to the radiator before the splitting
          // by a colour line (same reasoning for anticolour). The colour
          // that split is the colour appearing twice in the
          // radiator + emitted pair.
          // Thus, if we remove a colour index with the clustering,
          // we should look for a colour partner, else look for
          // an anticolour partner
          int colRemove = (event[iRad].col() == event[EmtTag].col())
                  ? event[iRad].col() : 0;

          int iPartner = 0;
          if (colRemove > 0 && col > 0)
            iPartner = FindCol(col,iRad,EmtTag,event,1,true)
                     + FindCol(col,iRad,EmtTag,event,2,true);
          else if (colRemove > 0 && acl > 0)
            iPartner = FindCol(acl,iRad,EmtTag,event,1,true)
                     + FindCol(acl,iRad,EmtTag,event,2,true);

          if ( allowedClustering( iRad, EmtTag, RecInit, iPartner, event ) ) {
            clus.push_back( Clustering(EmtTag, iRad, RecInit, iPartner,
                 pTLund(event[iRad],event[EmtTag],event[RecInit], pTdef) ));
              continue;
          }
        }

      // Second colour topology: Gluino emission

      } else {

        if ( (event[iRad].col() == event[EmtTag].acol())
           || (event[iRad].acol() == event[EmtTag].col())
           || (event[iRad].col() == event[EmtTag].col())
           || (event[iRad].acol() == event[EmtTag].acol()) ) {
          // For the rest, choose recoiler to have a common colour
          // tag with radiator, while not being the "Emitted"

          int col = -1;
          int acl = -1;

          if (event[iRad].isFinal() ) {

            if ( radID < 0 && event[iRad].colType() == -1) {
              acl = event[EmtTag].acol();
              col = event[iRad].col();
            } else if ( radID > 0 && event[iRad].colType() == 1 ) {
              col = event[EmtTag].col();
              acl = event[iRad].acol();
            } else {
              col = event[iRad].col();
              acl = event[iRad].acol();
            }

            int iRec = 0;
            if (col > 0) {
              iRec = FindCol(col,iRad,EmtTag,event,1,true);
              if ( (sign < 0) && (event[iRec].isFinal()) ) iRec = 0;
              // Not interested in pure QCD triples here.
              if ( !isSQCDrad && !isSQCDemt
                && (event[iRec].idAbs() < 10 || event[iRec].id() == 21) )
                iRec = 0;
              if (iRec != 0
               && allowedClustering( iRad, EmtTag, iRec, iRec, event) ) {
                clus.push_back( Clustering(EmtTag, iRad, iRec, iRec,
                     pTLund(event[iRad],event[EmtTag],event[iRec], pTdef) ));
                continue;
              }

              iRec = FindCol(col,iRad,EmtTag,event,2,true);
              if ( (sign < 0) && (event[iRec].isFinal()) ) iRec = 0;
              // Not interested in pure QCD triples here.
              if ( !isSQCDrad && !isSQCDemt
                && (event[iRec].idAbs() < 10 || event[iRec].id() == 21) )
                iRec = 0;
              if (iRec != 0
               && allowedClustering( iRad, EmtTag, iRec, iRec, event) ) {
                clus.push_back( Clustering(EmtTag, iRad, iRec, iRec,
                     pTLund(event[iRad],event[EmtTag],event[iRec], pTdef) ));
                continue;
              }
            }

            if (acl > 0) {
              iRec = FindCol(acl,iRad,EmtTag,event,1,true);
              if ( (sign < 0) && (event[iRec].isFinal()) ) iRec = 0;
              // Not interested in pure QCD triples here.
              if ( !isSQCDrad && !isSQCDemt
                && (event[iRec].idAbs() < 10 || event[iRec].id() == 21) )
                iRec = 0;
              if (iRec != 0
               && allowedClustering( iRad, EmtTag, iRec, iRec, event) ) {
                clus.push_back( Clustering(EmtTag, iRad, iRec, iRec,
                     pTLund(event[iRad],event[EmtTag],event[iRec], pTdef) ));
                continue;
              }

              iRec = FindCol(acl,iRad,EmtTag,event,2,true);
              if ( (sign < 0) && (event[iRec].isFinal()) ) iRec = 0;
              // Not interested in pure QCD triples here.
              if ( !isSQCDrad && !isSQCDemt
                && (event[iRec].idAbs() < 10 || event[iRec].id() == 21) )
                iRec = 0;
              if (iRec != 0
               && allowedClustering( iRad, EmtTag, iRec, iRec, event) ) {
                clus.push_back( Clustering(EmtTag, iRad, iRec, iRec,
                     pTLund(event[iRad],event[EmtTag],event[iRec], pTdef) ));
                continue;
              }
            }

          } else {

            // SM splittings already taken care of in findQCDTriple. Since
            // initial state splittings will not know if the true
            // colour-connected recoiler is a SM particle, any ISR splitting
            // of a SM particle will be included by findQCDTriple. To not
            // include the same splitting twice, continue for SM ISR radiator
            if ( !isSQCDrad || !isSQCDemt ) continue;

            // For an initial state radiator, always set recoiler
            // to the other initial state parton (recoil is taken
            // by full remaining system, so this is just a
            // labelling for such a process)
            int RecInit = 0;
            int iPartner = 0;
            for(int l = 0; l < int(PosInitPartn.size()); ++l)
              if (PosInitPartn[l] != iRad) RecInit = PosInitPartn[l];

            // Find the colour connected partner
            // Find colour index of radiator before splitting
            col = getRadBeforeCol(iRad, EmtTag, event);
            acl = getRadBeforeAcol(iRad, EmtTag, event);

            // Find the correct partner: If a colour line has split,
            // the partner is connected to the radiator before the splitting
            // by a colour line (same reasoning for anticolour). The colour
            // that split is the colour appearing twice in the
            // radiator + emitted pair.
            // Thus, if we remove a colour index with the clustering,
            // we should look for a colour partner, else look for
            // an anticolour partner
            int colRemove = (event[iRad].col() == event[EmtTag].col())
                    ? event[iRad].col() : 0;
            iPartner = (colRemove > 0)
                     ?   FindCol(col,iRad,EmtTag,event,1,true)
                       + FindCol(col,iRad,EmtTag,event,2,true)
                     :   FindCol(acl,iRad,EmtTag,event,1,true)
                       + FindCol(acl,iRad,EmtTag,event,2,true);

            if ( allowedClustering( iRad, EmtTag, RecInit, iPartner, event)) {
              clus.push_back( Clustering(EmtTag, iRad, RecInit, iPartner,
                   pTLund(event[iRad],event[EmtTag],event[RecInit], pTdef)));

              continue;
            }
          }
        }
      }
    }
  }

  // Done
  return clus;

}

//--------------------------------------------------------------------------

// Calculate and return the probability of a clustering.
// IN  Clustering : rad,rec,emt - System for which the splitting
//                  probability should be calcuated
// OUT splitting probability

double History::getProb(const Clustering & SystemIn) {

  // Get local copies of input system
  int Rad = SystemIn.emittor;
  int Rec = SystemIn.recoiler;
  int Emt = SystemIn.emitted;

  // Initialise shower probability
  double showerProb = 0.0;

  // If the splitting resulted in disallowed evolution variable,
  // disallow the splitting
  if (SystemIn.pT() <= 0.) return 0.;

  // Initialise all combinatorical factors
  double CF = 4./3.;
  double NC = 3.;
  // Flavour is known when reclustring, thus n_f=1
  double TR = 1./2.;

  // Split up in FSR and ISR
  bool isFSR = (state[Rad].isFinal() && state[Rec].isFinal());
  bool isFSRinREC = (state[Rad].isFinal() && !state[Rec].isFinal());
  bool isISR = !state[Rad].isFinal();

  // Check if this is the clustering 2->3 to 2->2.
  // If so, use weight for joined evolution
  int nFinal = 0;
  for(int i=0; i < state.size(); ++i)
    if (state[i].isFinal()) nFinal++;
  bool isLast = (nFinal == (mergingHooksPtr->hardProcess.nQuarksOut()
                           +mergingHooksPtr->hardProcess.nLeptonOut()+1));

  if (isISR) {
    // Find incoming particles

    int inP = 0;
    int inM = 0;
    for(int i=0;i< int(state.size()); ++i) {
      if (state[i].mother1() == 1) inP = i;
      if (state[i].mother1() == 2) inM = i;
    }
    // Construct dipole mass, eCM and sHat = x1*x2*s
    Vec4   sum     = state[Rad].p() + state[Rec].p() - state[Emt].p();
    double m2Dip = sum.m2Calc();
    double sHat = (state[inM].p() + state[inP].p()).m2Calc();
    // Energy fraction z=E_q1/E_qi in branch q(i)q(2) -> q(1)g(3)q(2)
    double z1 = m2Dip / sHat;
    // Virtuality of the splittings
    Vec4 Q1( state[Rad].p() - state[Emt].p() );
    Vec4 Q2( state[Rec].p() - state[Emt].p() );
    // Q^2 for emission off radiator line
    double Q1sq = -Q1.m2Calc();
    // pT^2 for emission off radiator line
    double pT1sq = pow(SystemIn.pT(),2);
    // Remember if massive particles involved: Mass corrections for
    // to g->QQ and Q->Qg splittings
    bool g2QQmassive = mergingHooksPtr->includeMassive()
        && state[Rad].id() == 21
        && ( (state[Emt].idAbs() >= 4 && state[Emt].idAbs() < 7)
          || (state[Emt].idAbs() > 1000000 && state[Emt].idAbs() < 1000010 )
          || (state[Emt].idAbs() > 2000000 && state[Emt].idAbs() < 2000010 ));
    bool Q2Qgmassive = mergingHooksPtr->includeMassive()
        && state[Emt].id() == 21
        && ( (state[Rad].idAbs() >= 4 && state[Rad].idAbs() < 7)
          || (state[Rad].idAbs() > 1000000 && state[Rad].idAbs() < 1000010 )
          || (state[Rad].idAbs() > 2000000 && state[Rad].idAbs() < 2000010 ));
    bool isMassive = mergingHooksPtr->includeMassive()
                    && ( g2QQmassive || Q2Qgmassive
                      || state[Emt].id() == 1000021);
    double m2Emt0 = pow(particleDataPtr->m0(state[Emt].id()),2);
    double m2Rad0 = pow(particleDataPtr->m0(state[Rad].id()),2);

    // Correction of virtuality for massive splittings
    if ( g2QQmassive)      Q1sq += m2Emt0;
    else if (Q2Qgmassive)  Q1sq += m2Rad0;

    // pT0 dependence!!!
    double pT0sq = pow(mergingHooksPtr->pT0ISR(),2);
    double Q2sq = -Q2.m2Calc();

    // Correction of virtuality of other splitting
    bool g2QQmassiveRec = mergingHooksPtr->includeMassive()
        && state[Rec].id() == 21
        && ( state[Emt].idAbs() >= 4 && state[Emt].idAbs() < 7);
    bool Q2QgmassiveRec = mergingHooksPtr->includeMassive()
        && state[Emt].id() == 21
        && ( state[Rec].idAbs() >= 4 && state[Rec].idAbs() < 7);
    double m2Rec0 = pow(particleDataPtr->m0(state[Rad].id()),2);
    if ( g2QQmassiveRec)      Q2sq += m2Emt0;
    else if (Q2QgmassiveRec)  Q2sq += m2Rec0;

    bool hasJoinedEvol = (state[Emt].id() == 21
                       || state[Rad].id() == state[Rec].id());

    // Initialise normalization factor multiplying the splitting
    // function numerator
    double fac = 1.;
    if ( mergingHooksPtr->pickByFull() || mergingHooksPtr->pickBySumPT()) {
      double facJoined  = ( Q2sq + pT0sq/(1.-z1) )
                        * 1./(Q1sq*Q2sq + pT0sq*sHat + pow(pT0sq/(1.-z1),2));
      double facSingle = mergingHooksPtr->nonJoinedNorm()*1./( pT1sq + pT0sq);

      fac = (hasJoinedEvol && isLast) ? facJoined : facSingle;

    } else if (mergingHooksPtr->pickByPoPT2()) {
      fac = 1./(pT1sq + pT0sq);
    } else {
      string message="Error in History::getProb: Scheme for calculating";
      message+=" shower splitting probability is undefined.";
      infoPtr->errorMsg(message);
    }

    // Calculate shower splitting probability:
    // Splitting functions*normalization*ME reweighting factors

    // Calculate branching probability for q -> q g
    if ( state[Emt].id() == 21 && state[Rad].id() != 21) {
      // Find splitting kernel
      double num = CF*(1. + pow(z1,2)) / (1.-z1);
      if (isMassive) num -= CF * z1 * (1.-z1) * (m2Rad0/pT1sq);

      // Find ME reweighting factor
      double meReweighting = 1.;
      // Find the number of final state coloured particles, apart
      // from those coming from the hard process
      int nCol = 0;
      for(int i=0; i < state.size(); ++i)
        if (state[i].isFinal() && state[i].colType() != 0
          && !mergingHooksPtr->hardProcess.matchesAnyOutgoing(i,state) )
          nCol++;
      // For first splitting of single vector boson production,
      // apply ME corrections
      if (nCol == 1
       && int(mergingHooksPtr->hardProcess.hardIntermediate.size()) == 1) {
        double sH = m2Dip / z1;
        double tH = -Q1sq;
        double uH = Q1sq - m2Dip * (1. - z1) / z1;
        double misMatch = (uH*tH - (uH + tH)*pT0sq/(1.-z1)
                          + pow(pT0sq/(1.-z1),2) ) / (uH*tH);
        meReweighting *= (tH*tH + uH*uH + 2. * m2Dip * sH)
                       / (sH*sH + m2Dip*m2Dip);
        meReweighting *= misMatch;
      }
      // Multiply factors
      showerProb = num*fac*meReweighting;

    // Calculate branching probability for g -> g g
    } else if ( state[Emt].id() == 21 && state[Rad].id() == 21) {
      // Calculate splitting kernel
      double num = 2.*NC*pow2(1. - z1*(1.-z1)) / (z1*(1.-z1));

      // Include ME reweighting for higgs!!
      // Find ME reweighting factor
      double meReweighting = 1.;
      // Find the number of final state coloured particles, apart
      // from those coming from the hard process
      int nCol = 0;
      for(int i=0; i < state.size(); ++i)
        if (state[i].isFinal() && state[i].colType() != 0
          && !mergingHooksPtr->hardProcess.matchesAnyOutgoing(i,state) )
          nCol++;
      // For first splitting of single vector boson production,
      // apply ME corrections
      if ( nCol == 1
       && mergingHooksPtr->getProcessString().compare("pp>h") == 0
       && int(mergingHooksPtr->hardProcess.hardIntermediate.size()) == 1 ) {
        double sH = m2Dip / z1;
        double tH = -Q1sq;
        double uH = Q1sq - m2Dip * (1. - z1) / z1;
        meReweighting *= 0.5
                       * (pow4(sH) + pow4(tH) + pow4(uH) + pow4(m2Dip))
                       / pow2(sH*sH - m2Dip * (sH - m2Dip));
      }

      // Multiply factors
      showerProb = num*fac*meReweighting;

    // Calculate branching probability for q -> g q
    } else if ( state[Emt].id() != 21 && state[Rad].id() != 21 ) {
      // Calculate splitting kernel
      double num = CF*(1. + pow2(1.-z1)) / z1;

      // Include ME reweighting for higgs!!
      // Find ME reweighting factor
      double meReweighting = 1.;
      // Find the number of final state coloured particles, apart
      // from those coming from the hard process
      int nCol = 0;
      for ( int i=0; i < state.size(); ++i )
        if ( state[i].isFinal() && state[i].colType() != 0
          && !mergingHooksPtr->hardProcess.matchesAnyOutgoing(i,state) )
          nCol++;
      // For first splitting of single vector boson production,
      // apply ME corrections
      if (nCol == 1
       && mergingHooksPtr->getProcessString().compare("pp>h") == 0
       && int(mergingHooksPtr->hardProcess.hardIntermediate.size()) == 1) {
        double sH = m2Dip / z1;
        double uH = Q1sq - m2Dip * (1. - z1) / z1;
        meReweighting *= (sH*sH + uH*uH)
                       / (sH*sH + pow2(sH -m2Dip));
      }

      // Multiply factors
      showerProb = num*fac*meReweighting;

    // Calculate branching probability for g -> q qbar
    } else if ( state[Emt].id() != 21 && state[Rad].id() == 21 ) {
      // Calculate splitting kernel
      double num = TR * ( pow(z1,2) + pow(1.-z1,2) );
      if (isMassive) num += TR * 2.*z1*(1.-z1)*(m2Emt0/pT1sq);
      // Calculate ME reweighting factor
      double meReweighting = 1.;
      // Find the number of final state coloured particles, apart
      // from those coming from the hard process
      int nCol = 0;
      for ( int i=0; i < state.size(); ++i )
        if ( state[i].isFinal() && state[i].colType() != 0
         && !mergingHooksPtr->hardProcess.matchesAnyOutgoing(i,state) )
          nCol++;
      // For first splitting of single vector boson production,
      // apply ME corrections
      if (nCol == 1
        && int(mergingHooksPtr->hardProcess.hardIntermediate.size()) == 1) {
        double sH = m2Dip / z1;
        double tH = -Q1sq;
        double uH = Q1sq - m2Dip * (1. - z1) / z1;
        swap( tH, uH);
        double misMatch = ( uH - pT0sq/(1.-z1) ) / uH;
        double me = (sH*sH + uH*uH + 2. * m2Dip * tH)
                  / (pow2(sH - m2Dip) + m2Dip*m2Dip);
        // Weight with me/overestimate
        meReweighting *= me;
        meReweighting *= misMatch;
      }
      // Multiply factors
      showerProb = num*fac*meReweighting;

    // Print error if no kernel calculated
    } else {
      string message = "Error in History::getProb: Splitting kernel"
        " undefined in ISR in clustering.";
      infoPtr->errorMsg(message);
    }

    // If corrected pT below zero in ISR, put probability to zero
    double m2Sister0 = pow(state[Emt].m0(),2);
    double pT2corr = (Q1sq - z1*(m2Dip + Q1sq)*(Q1sq + m2Sister0)/m2Dip);
    if (pT2corr < 0.) showerProb *= 1e-9;

    // If creating heavy quark by Q -> gQ then next need g -> Q + Qbar.
    // So minimum total mass2 is 4 * m2Sister, but use more to be safe.
    if ( state[Emt].id() == state[Rad].id()
       && ( state[Rad].idAbs() == 4 || state[Rad].idAbs() == 5 )) {
      double m2QQsister =  2.*4.*m2Sister0;
      double pT2QQcorr = Q1sq - z1*(m2Dip + Q1sq)*(Q1sq + m2QQsister)
                       / m2Dip;
      if (pT2QQcorr < 0.0) showerProb *= 1e-9;
    }

    if (mergingHooksPtr->includeRedundant()) {
      // Initialise the spacelike shower alpha_S
      AlphaStrong* asISR = mergingHooksPtr->AlphaS_ISR();
      double as = (*asISR).alphaS(pT1sq + pT0sq) / (2.*M_PI);
      // Multiply with alpha_S
      showerProb *= as;
    }

  // Done for ISR case, begin FSR case

  }  else if (isFSR || isFSRinREC) {

    // Construct dipole mass
    Vec4   sum     = state[Rad].p() + state[Rec].p() + state[Emt].p();
    double m2Dip = sum.m2Calc();
    // Construct 2->3 variables for FSR
    double x1 = 2. * (sum * state[Rad].p()) / m2Dip;
    double x2 = 2. * (sum * state[Rec].p()) / m2Dip;
    double prop1  = max(1e-12, 1. - x1);
    double prop2  = max(1e-12, 1. - x2);
    double x3     = max(1e-12, 2. - x1 - x2);
    // Energy fraction z=E_q1/E_qi in branch q(i)q(2) -> q(1)g(3)q(2)
    double z1 = x1/(x1 + x3);

    // Virtuality of the splittings
    Vec4 Q1( state[Rad].p() + state[Emt].p() );
    Vec4 Q2( state[Rec].p() + state[Emt].p() );

    // Q^2 for emission off radiator line
    double Q1sq = Q1.m2Calc();
    // pT^2 for emission off radiator line
    double pT1sq = pow(SystemIn.pT(),2);
    // Q^2 for emission off recoiler line
    double Q2sq = Q2.m2Calc();

    // Remember if radiator or recoiler is massive.
    bool isMassiveRad = ( state[Rad].idAbs() >= 4
                       && state[Rad].id() != 21 );
    bool isMassiveRec = ( state[Rec].idAbs() >= 4
                       && state[Rec].id() != 21 );
    // Correction of virtuality for massive splittings.
    double m2Rad0 = pow(particleDataPtr->m0(state[Rad].id()),2);
    double m2Rec0 = pow(particleDataPtr->m0(state[Rad].id()),2);
    if ( mergingHooksPtr->includeMassive() && isMassiveRad ) Q1sq -= m2Rad0;
    if ( mergingHooksPtr->includeMassive() && isMassiveRec ) Q2sq -= m2Rec0;

    // Initialise normalization factor multiplying the splitting
    // function numerator
    double fac = 1.;
    if ( mergingHooksPtr->pickByFull() || mergingHooksPtr->pickBySumPT()) {
      double facJoined = (1.-z1)/Q1sq * m2Dip/( Q1sq + Q2sq );
      double facSingle = mergingHooksPtr->fsrInRecNorm() * 1./ pT1sq;
      fac = (!isFSRinREC && isLast) ? facJoined : facSingle;
    } else if (mergingHooksPtr->pickByPoPT2()) {
      fac = 1. / pT1sq;
    } else {
      string message="Error in History::getProb: Scheme for calculating";
      message+=" shower splitting probability is undefined.";
      infoPtr->errorMsg(message);
    }
    // Calculate shower splitting probability:
    // Splitting functions*normalization*ME reweighting factors

    // Calculate branching probability for g -> g_1 g_2
    if ( state[Emt].id() == 21 && state[Rad].id() == 21) {
      // Calculate splitting kernel
      double num = 0.5* NC * (1. + pow3(z1)) / (1.-z1);
      // Multiply factors
      showerProb = num*fac;

    // Calculate branching probability for q -> q g with quark recoiler
    } else if ( state[Emt].id() == 21
             && state[Rad].id() != 21 && state[Rec].id() != 21) {
      // For a qqbar dipole in FSR, ME corrections exist and the
      // splitting function "z-weight" is set to 1.0 (only for 2->2 ??)
      double num = CF * 2./(1.-z1);
      // Find the number of final state coloured particles, apart
      // from those coming from the hard process
      int nCol = 0;
        for(int i=0; i < state.size(); ++i)
          if (state[i].isFinal() && state[i].colType() != 0
          && !mergingHooksPtr->hardProcess.matchesAnyOutgoing(i,state) )
            nCol++;
      // Calculate splitting kernel
      if ( nCol > 3
        || int(mergingHooksPtr->hardProcess.hardIntermediate.size()) > 1)
        num = CF * (1. + pow2(z1)) /(1.-z1);
      // Calculate ME reweighting factor
      double meReweighting = 1.;
      // Correct if this is the process created by the first
      // FSR splitting of a 2->2 process
      if ( nCol == 3
        && int(mergingHooksPtr->hardProcess.hardIntermediate.size()) == 1 ) {
        // Calculate the ME reweighting factor
        double ShowerRate1       = 2./( x3 * prop2 );
        double meDividingFactor1 = prop1 / x3;
        double me                = (pow(x1,2) + pow(x2,2))/(prop1*prop2);
        meReweighting = meDividingFactor1 * me / ShowerRate1;
      }
      // Multiply factors
      showerProb = num*fac*meReweighting;

    // Calculate branching probability for q1 -> q2 w+- with quark recoiler
    } else if ( state[Emt].idAbs() == 24
             && state[Rad].id()    != 21 && state[Rec].id() != 21 ) {
      double m2W = state[Emt].p().m2Calc();
      double num = ( 3.*pow2(m2W / m2Dip)
                   + 2.* (m2W/m2Dip)*(x1 + x2)
                   + pow2(x1) + pow2(x2) ) / ( prop1*prop2 )
                 - (m2W/m2Dip) / pow2(prop1)
                 - (m2W/m2Dip) / pow2(prop2);
      // Multiply factors
      showerProb = num*fac;

    // Calculate branching probability for q -> q g with gluon recoiler
    } else if ( state[Emt].id() == 21 && state[Rad].id() != 21
      && state[Rec].id() == 21 ) {
      // For qg /qbarg dipoles, the splitting function is
      // calculated and not weighted by a ME correction factor
      // Shower splitting function
      double num = CF * (1. + pow2(z1)) / (1.-z1);
      showerProb = num*fac;

    // Calculate branching probability for g -> q qbar
    } else if ( state[Emt].id() != 21 ) {
      // Get flavour of quark / antiquark
      int flavour = state[Emt].id();
      // Get correct masses for the quarks
      // (needed to calculate splitting function?)
      double mFlavour = particleDataPtr->m0(flavour);
      // Get mass of quark/antiquark pair
      double mDipole = m(state[Rad].p(), state[Emt].p());
      // Factor determining if gluon decay was allowed
      double beta = sqrtpos( 1. - 4.*pow2(mFlavour)/pow2(mDipole) );
      // Shower splitting function
      double num = 0.5*TR * ( z1*z1 + (1.-z1)*(1.-z1) );
      if (beta <= 0.) beta = 0.;

      showerProb = num*fac*beta;

    // Print error if no kernel calculated
    } else {
      string message="Error in History::getProb: Splitting kernel undefined";
      message+=" in FSR clustering.";
      infoPtr->errorMsg(message);
    }

    if (mergingHooksPtr->includeRedundant()) {
      // Initialise the spacelike shower alpha_S
      AlphaStrong* asFSR = mergingHooksPtr->AlphaS_FSR();
      double as = (*asFSR).alphaS(pT1sq) / (2.*M_PI);
      // Multiply with alpha_S
      showerProb *= as;
    }

    // Done for FSR
  } else {
    string message="Error in History::getProb: Radiation could not be";
    message+=" interpreted as FSR or ISR.";
    infoPtr->errorMsg(message);
  }

  if (showerProb <= 0.) showerProb = 0.;

//  // Different coupling constants for qcd and ew splittings
//  if ( state[Emt].colType() != 0 ) {
//    AlphaStrong* asFSR = mergingHooksPtr->AlphaS_FSR();
//    double as = (*asFSR).alphaS(91.188*91.188) / (2.*M_PI);
//    showerProb *= as;
//  } else {
//    AlphaEM* aEMFSR = mergingHooksPtr->AlphaEM_FSR();
//    double aEM = (*aEMFSR).alphaEM(91.188*91.188) / (2.*M_PI);
//    showerProb *= aEM;
//  }

  // Done
  return showerProb;
}


//--------------------------------------------------------------------------

// Set up the beams (fill the beam particles with the correct
// current incoming particles) to allow calculation of splitting
// probability.
// For interleaved evolution, set assignments dividing PDFs into
// sea and valence content. This assignment is, until a history path
// is chosen and a first trial shower performed, not fully correct
// (since content is chosen form too high x and too low scale). The
// assignment used for reweighting will be corrected after trial
// showering

void History::setupBeams() {

  // Do nothing for empty event, possible if sequence of
  // clusterings was ill-advised in that it results in
  // colour-disconnected states
  if (state.size() < 4) return;
  // Do nothing for e+e- beams
  if ( state[3].colType() == 0 ) return;
  if ( state[4].colType() == 0 ) return;

  // Incoming partons to hard process are stored in slots 3 and 4.
  int inS = 0;
  int inP = 0;
  int inM = 0;
  for(int i=0;i< int(state.size()); ++i) {
    if (state[i].mother1() == 1) inP = i;
    if (state[i].mother1() == 2) inM = i;
  }

  // Save some info before clearing beams
  // Mothers of incoming partons companion code
  int motherPcompRes = -1;
  int motherMcompRes = -1;

  bool sameFlavP = false;
  bool sameFlavM = false;

  if (mother) {
    int inMotherP = 0;
    int inMotherM = 0;
    for(int i=0;i< int(mother->state.size()); ++i) {
      if (mother->state[i].mother1() == 1) inMotherP = i;
      if (mother->state[i].mother1() == 2) inMotherM = i;
    }
    sameFlavP = (state[inP].id() == mother->state[inMotherP].id());
    sameFlavM = (state[inM].id() == mother->state[inMotherM].id());

    motherPcompRes = (sameFlavP) ? beamA[0].companion() : -2;
    motherMcompRes = (sameFlavM) ? beamB[0].companion() : -2;
  }

  // Append the current incoming particles to the beam
  beamA.clear();
  beamB.clear();

  // Get energy of incoming particles
  double Ep = 2. * state[inP].e();
  double Em = 2. * state[inM].e();

  // If incoming partons are massive then recalculate to put them massless.
  if (state[inP].m() != 0. || state[inM].m() != 0.) {
    Ep = state[inP].pPos() + state[inM].pPos();
    Em = state[inP].pNeg() + state[inM].pNeg();
  }

  // Add incoming hard-scattering partons to list in beam remnants.
  double x1 = Ep / state[inS].m();
  beamA.append( inP, state[inP].id(), x1);
  double x2 = Em / state[inS].m();
  beamB.append( inM, state[inM].id(), x2);

  // Scale. For ME multiplicity history, put scale to mu_F
  // (since sea/valence quark content is chosen from this scale)
  double scalePDF = (mother) ? scale : infoPtr->QFac();
  // Find whether incoming partons are valence or sea. Store.
  // Can I do better, e.g. by setting the scale to the hard process
  // scale (= M_W) or by replacing one of the x values by some x/z??
  beamA.xfISR( 0, state[inP].id(), x1, scalePDF*scalePDF);
  if (!mother) {
    beamA.pickValSeaComp();
  }  else {
    beamA[0].companion(motherPcompRes);
  }
  beamB.xfISR( 0, state[inM].id(), x2, scalePDF*scalePDF);
  if (!mother) {
    beamB.pickValSeaComp();
  } else {
    beamB[0].companion(motherMcompRes);
  }

}

//--------------------------------------------------------------------------

// Calculate the PDF ratio used in the argument of the no-emission
// probability

double History::pdfForSudakov() {

  // Do nothing for e+e- beams
  if ( state[3].colType() == 0 ) return 1.0;
  if ( state[4].colType() == 0 ) return 1.0;

  // Check if splittings was ISR or FSR
  bool FSR = (   mother->state[clusterIn.emittor].isFinal()
             && mother->state[clusterIn.recoiler].isFinal());
  bool FSRinRec = (   mother->state[clusterIn.emittor].isFinal()
                  && !mother->state[clusterIn.recoiler].isFinal());

  // Done for pure FSR
  if (FSR) return 1.0;

  int iInMother = (FSRinRec)? clusterIn.recoiler : clusterIn.emittor;
  //  Find side of event that was reclustered
  int side = ( mother->state[iInMother].pz() > 0 ) ? 1 : -1;

  int inP = 0;
  int inM = 0;
  for(int i=0;i< int(state.size()); ++i) {
    if (state[i].mother1() == 1) inP = i;
    if (state[i].mother1() == 2) inM = i;
  }

  // Save mother id
  int idMother = mother->state[iInMother].id();
  // Find daughter position and id
  int iDau = (side == 1) ? inP : inM;
  int idDaughter = state[iDau].id();
  // Get mother x value
  double xMother = 2. * mother->state[iInMother].e() / mother->state[0].e();
  // Get daughter x value of daughter
  double xDaughter = 2.*state[iDau].e() / state[0].e(); // x1 before isr

  // Calculate pdf ratio
  double ratio = getPDFratio(side, true, false, idMother, xMother, scale,
                   idDaughter, xDaughter, scale);

  // For FSR with incoming recoiler, maximally return 1.0, as
  // is done in Pythia::TimeShower.
  // For ISR, return ratio
  return ( (FSRinRec)? min(1.,ratio) : ratio);
}

//--------------------------------------------------------------------------

// Calculate the hard process matrix element to include in the selection
// probability.

double History::hardProcessME( const Event& event ) {

  // Get hard process.
  string process = mergingHooksPtr->getProcessString();
  double result = 1.;

  if ( process.compare("pp>e+ve") == 0
    || process.compare("pp>e-ve~") == 0
    || process.compare("pp>LEPTONS,NEUTRINOS") == 0 ) {
    // Do nothing for incomplete process.
    int nFinal = 0;
    for ( int i=0; i < int(event.size()); ++i )
      if ( event[i].isFinal() ) nFinal++;
    if ( nFinal != 2 ) return 1.;
    // Get W-boson mass and width.
    double mW = particleDataPtr->m0(24);
    double gW = particleDataPtr->mWidth(24) / mW;
    // Get incoming particles.
    int inP = (event[3].pz() > 0) ? 3 : 4;
    int inM = (event[3].pz() > 0) ? 4 : 3;
    // Get outgoing particles.
    int outP = 0;
    for ( int i=0; i < int(event.size()); ++i ) {
      if ( event[i].isFinal() && event[i].px() > 0 ) outP = i;
    }
    // Get Mandelstam variables.
    double sH = (event[inP].p() + event[inM].p()).m2Calc();
    double tH = (event[inP].p() - event[outP].p()).m2Calc();
    double uH = - sH - tH;

    // Return kinematic part of matrix element.
    result = ( 1. + (tH - uH)/sH ) / ( pow2(sH - mW*mW) + pow2(sH*gW) );
  } else
    result = mergingHooksPtr->hardProcessME(event);

  return result;

}

//--------------------------------------------------------------------------

// Perform the clustering of the current state and return the
// clustered state.
// IN Clustering : rad,rec,emt triple to be clustered to two partons
// OUT clustered state

Event History::cluster( const Clustering & inSystem ) {

  // Initialise tags of particles to be changed
  int Rad = inSystem.emittor;
  int Rec = inSystem.recoiler;
  int Emt = inSystem.emitted;
  // Initialise eCM,mHat
  double eCM = state[0].e();
  // Flags for type of radiation
  int radType = state[Rad].isFinal() ? 1 : -1;
  int recType = state[Rec].isFinal() ? 1 : -1;

  // Construct the clustered event
  Event NewEvent = Event();
  NewEvent.init("(hard process-modified)", particleDataPtr);
  NewEvent.clear();
  // Copy all unchanged particles to NewEvent
  for (int i = 0; i < state.size(); ++i)
    if ( i != Rad && i != Rec && i != Emt )
      NewEvent.append( state[i] );

  // Copy all the junctions one by one
  for (int i = 0; i < state.sizeJunction(); ++i)
    NewEvent.appendJunction( state.getJunction(i) );
  // Find an appropriate scale for the hard process
  double mu = choseHardScale(state);
  // Initialise scales for new event
  NewEvent.saveSize();
  NewEvent.saveJunctionSize();
  NewEvent.scale(mu);
  NewEvent.scaleSecond(mu);

  // Set properties of radiator/recoiler after the clustering
  // Recoiler properties will be unchanged
  Particle RecBefore = Particle( state[Rec] );
  RecBefore.setEvtPtr(&NewEvent);
  RecBefore.daughters(0,0);
  // Find flavour of radiator before splitting
  int radID = getRadBeforeFlav(Rad, Emt, state);
  int recID = state[Rec].id();
  Particle RadBefore = Particle( state[Rad] );
  RadBefore.setEvtPtr(&NewEvent);
  RadBefore.id(radID);
  RadBefore.daughters(0,0);
  // Put dummy values for colours
  RadBefore.cols(RecBefore.acol(),RecBefore.col());

  // Reset status if the reclustered radiator is a resonance.
  if ( particleDataPtr->isResonance(radID) && radType == 1)
    RadBefore.status(state[Rad].status());

  // Put mass for radiator and recoiler
  double radMass = particleDataPtr->m0(radID);
  double recMass = particleDataPtr->m0(recID);
  if (radType == 1 ) RadBefore.m(radMass);
  else RadBefore.m(0.0);
  if (recType == 1 ) RecBefore.m(recMass);
  else RecBefore.m(0.0);

  // Construct momenta and  colours of clustered particles
  // ISR/FSR splittings are treated differently
  if ( radType + recType == 2 && state[Emt].idAbs() != 24 ) {
    // Clustering of final(rad)/final(rec) dipole splitting
    // Get eCM of (rad,rec,emt) triple
    Vec4   sum     = state[Rad].p() + state[Rec].p() + state[Emt].p();
    double eCMME   = sum.mCalc();

    // Define radiator and recoiler back-to-back in the dipole
    // rest frame [=(state[rad]+state[emt])+state[rec] rest frame]
    Vec4 Rad4mom;
    Vec4 Rec4mom;
    double mDsq   = pow2(eCMME);
    // If possible, keep the invariant mass of the radiator.
    double mRsq   = (radID == state[Rad].id() )
                  ? abs(state[Rad].p().m2Calc())
                  : pow2(particleDataPtr->m0(radID));
    double mSsq   = abs(state[Rec].p().m2Calc());
    double a = 0.5*sqrt(mDsq);
    double b = 0.25*(mRsq - mSsq) / a;
    double c = sqrt(pow2(a) + pow2(b) - 2.*a*b - mSsq);

    Rad4mom.p( 0., 0., c, a+b);
    Rec4mom.p( 0., 0.,-c, a-b);

    // Find boost from Rad4mom+Rec4mom rest frame to event cm frame
    Vec4 old1 = Vec4(state[Rad].p() + state[Emt].p());
    Vec4 old2 = Vec4(state[Rec].p());
    RotBstMatrix fromCM;
    fromCM.fromCMframe(old1, old2);
    // Transform momenta
    Rad4mom.rotbst(fromCM);
    Rec4mom.rotbst(fromCM);

    RadBefore.p(Rad4mom);
    RecBefore.p(Rec4mom);
    RadBefore.m(sqrt(mRsq));
    RecBefore.m(sqrt(mSsq));

  } else if ( radType + recType == 2 && state[Emt].idAbs() == 24 ) {
    // Clustering of final(rad)/final(rec) dipole splitting
    // Get eCM of (rad,rec,emt) triple
    Vec4   sum     = state[Rad].p() + state[Rec].p() + state[Emt].p();
    double eCMME   = sum.mCalc();

    // Define radiator and recoiler back-to-back in the dipole
    // rest frame [=(state[rad]+state[emt])+state[rec] rest frame]
    Vec4 Rad4mom;
    Vec4 Rec4mom;
    double mDsq   = pow2(eCMME);
    // If possible, keep the invariant mass of the radiator.
    double mRsq   = (radID == state[Rad].id() )
                  ? abs(state[Rad].p().m2Calc())
                  : pow2(particleDataPtr->m0(radID));
    double mSsq   = abs(state[Rec].p().m2Calc());
    double a = 0.5*sqrt(mDsq);
    double b = 0.25*(mRsq - mSsq) / a;
    double c = sqrt(pow2(a) + pow2(b) - 2.*a*b - mSsq);

    Rad4mom.p( 0., 0., c, a+b);
    Rec4mom.p( 0., 0.,-c, a-b);

    // Find boost from Rad4mom+Rec4mom rest frame to event cm frame
    Vec4 old1 = Vec4(state[Rad].p() + state[Emt].p());
    Vec4 old2 = Vec4(state[Rec].p());
    RotBstMatrix fromCM;
    fromCM.fromCMframe(old1, old2);
    // Transform momenta
    Rad4mom.rotbst(fromCM);
    Rec4mom.rotbst(fromCM);

    RadBefore.p(Rad4mom);
    RecBefore.p(Rec4mom);
    RadBefore.m(sqrt(mRsq));
    RecBefore.m(sqrt(mSsq));

  } else if ( radType + recType == 0 ) {
    // Clustering of final(rad)/initial(rec) dipole splitting
    // Get eCM of (rad,rec,emt) triple
    Vec4   sum     = state[Rad].p() + state[Rec].p() + state[Emt].p();
    double eCMME   = sum.mCalc();
    // Define radiator and recoiler back-to-back in the dipole
    // rest frame [=(state[rad]+state[emt])+state[rec] rest frame]
    Vec4 Rad4mom;
    Vec4 Rec4mom;
    double mDsq   = pow2(eCMME);
    // If possible, keep the invariant mass of the radiator.
    double mRsq   = (radID == state[Rad].id() )
                  ? abs(state[Rad].p().m2Calc())
                  : pow2(particleDataPtr->m0(radID));
    double mSsq   = abs(state[Rec].p().m2Calc());
    double a = 0.5*sqrt(mDsq);
    double b = 0.25*(mRsq - mSsq) / a;
    double c = sqrt(pow2(a) + pow2(b) - 2.*a*b - mSsq);

    Rad4mom.p( 0., 0., c, a+b);
    Rec4mom.p( 0., 0.,-c, a-b);

    // Find boost from Rad4mom+Rec4mom rest frame to event cm frame
    Vec4 old1 = Vec4(state[Rad].p() + state[Emt].p());
    Vec4 old2 = Vec4(state[Rec].p());
    RotBstMatrix fromCM;
    fromCM.fromCMframe(old1, old2);
    // Transform momenta
    Rad4mom.rotbst(fromCM);
    Rec4mom.rotbst(fromCM);

    // Rescale recoiler momentum
    Rec4mom = 2.*state[Rec].p() - Rec4mom;

    RadBefore.p(Rad4mom);
    RecBefore.p(Rec4mom);
    RadBefore.m(sqrt(mRsq));

    // Set mass of initial recoiler to zero
    RecBefore.m( 0.0 );

  } else {
    // Clustering of initial(rad)/initial(rec) dipole splitting
    // We want to cluster: Meaning doing the inverse of a process
    //            ( pDaughter + pRecoiler -> pOut )
    //        ==> ( pMother + pPartner -> pOut' + pSister )
    // produced by an initial state splitting. The matrix element
    // provides us with pMother, pPartner, pSister and pOut'
    Vec4 pMother( state[Rad].p() );
    Vec4 pSister( state[Emt].p() );
    Vec4 pPartner( state[Rec].p() );
    Vec4 pDaughter( 0.,0.,0.,0. );
    Vec4 pRecoiler( 0.,0.,0.,0. );

    // Find side that radiates event (mother moving in sign * p_z direction).
    int sign = (state[Rad].pz() > 0.) ? 1 : -1;

    // Find rotation by phi that would have been done for a
    // splitting daughter -> mother + sister
    double phi = pSister.phi();
    // Find rotation with -phi
    RotBstMatrix rot_by_mphi;
    rot_by_mphi.rot(0.,-phi);
    // Find rotation with +phi
    RotBstMatrix rot_by_pphi;
    rot_by_pphi.rot(0.,phi);

    // Transform pMother and outgoing momenta
    pMother.rotbst( rot_by_mphi );
    pSister.rotbst( rot_by_mphi );
    pPartner.rotbst( rot_by_mphi );
    for(int i=3; i< NewEvent.size(); ++i)
      NewEvent[i].rotbst( rot_by_mphi );

    // Get mother and partner x values
    // x1 after isr
    double x1 = 2. * pMother.e() / eCM;
    // x2 after isr
    double x2 = 2. * pPartner.e() / eCM;

    pDaughter.p( pMother - pSister);
    pRecoiler.p( pPartner );

    // Find boost from event cm frame to rest frame of
    // of-shell daughter + on-shell recoiler
    RotBstMatrix from_CM_to_DR;
    if (sign == 1)
      from_CM_to_DR.toCMframe(pDaughter, pRecoiler);
    else
      from_CM_to_DR.toCMframe(pRecoiler, pDaughter);

    // Transform all momenta
    pMother.rotbst( from_CM_to_DR );
    pPartner.rotbst( from_CM_to_DR );
    pSister.rotbst( from_CM_to_DR );
    for(int i=3; i< NewEvent.size(); ++i)
      NewEvent[i].rotbst( from_CM_to_DR );

    // Find theta angle between pMother and z-axis and undo
    // rotation that would have been done by shower
    double theta = pMother.theta();
    if ( pMother.px() < 0. ) theta *= -1.;
    if (sign == -1) theta += M_PI;

    // Find rotation by +theta
    RotBstMatrix rot_by_ptheta;
    rot_by_ptheta.rot(theta, 0.);

    // Transform all momenta
    pMother.rotbst( rot_by_ptheta );
    pPartner.rotbst( rot_by_ptheta );
    pSister.rotbst( rot_by_ptheta );
    for(int i=3; i< NewEvent.size(); ++i)
      NewEvent[i].rotbst( rot_by_ptheta );

    // Find z of the splitting
    Vec4 qDip( pMother - pSister);
    Vec4 qAfter(pMother + pPartner);
    Vec4 qBefore(qDip + pPartner);
    double z = qBefore.m2Calc() / qAfter.m2Calc();

    // Calculate new e_CM^2
    double x1New = z*x1; // x1 before isr
    double x2New = x2;   // x2 before isr
    double sHat = x1New*x2New*eCM*eCM;

    // Construct daughter and recoiler momenta
    pDaughter.p( 0., 0.,  sign*0.5*sqrt(sHat), 0.5*sqrt(sHat));
    pRecoiler.p( 0., 0., -sign*0.5*sqrt(sHat), 0.5*sqrt(sHat));

    // Find boost from current (daughter+recoiler rest frame)
    // frame to rest frame of daughter+unchanged recoiler to
    // recover the old x2 value
    RotBstMatrix from_DR_to_CM;
    from_DR_to_CM.bst( 0., 0., sign*( x1New - x2New ) / ( x1New + x2New ) );

    // Correct for momentum mismatch by transforming all momenta
    pMother.rotbst( from_DR_to_CM );
    pPartner.rotbst( from_DR_to_CM );
    pSister.rotbst( from_DR_to_CM );
    pDaughter.rotbst( from_DR_to_CM );
    pRecoiler.rotbst( from_DR_to_CM );
    for(int i=3; i< NewEvent.size(); ++i)
      NewEvent[i].rotbst( from_DR_to_CM );

    // Transform pMother and outgoing momenta
    pMother.rotbst( rot_by_pphi );
    pPartner.rotbst( rot_by_pphi );
    pSister.rotbst( rot_by_pphi );
    pDaughter.rotbst( rot_by_pphi );
    pRecoiler.rotbst( rot_by_pphi );
    for(int i=3; i< NewEvent.size(); ++i)
      NewEvent[i].rotbst( rot_by_pphi );

    // Set momenta of particles to be attached to new event record
    RecBefore.p( pRecoiler );
    RadBefore.p( pDaughter );
    if (RecBefore.pz() > 0.) RecBefore.mother1(1);
    else RecBefore.mother1(2);
    if (RadBefore.pz() > 0.) RadBefore.mother1(1);
    else RadBefore.mother1(2);

  }

  // Put some dummy production scales for RecBefore, RadBefore
  RecBefore.scale(mu);
  RadBefore.scale(mu);

  // Append new recoiler and find new radiator colour
  NewEvent.append(RecBefore);

  // Assign the correct colour to re-clustered radiator.
  // Keep old radiator colours for electroweak emission.
  int emtID = state[Emt].id();
  if ( emtID == 22 || emtID == 23 || abs(emtID) == 24 )
    RadBefore.cols( state[Rad].col(), state[Rad].acol() );
  // For QCD, carefully construct colour.
  else if ( !connectRadiator( RadBefore, radType, RecBefore, recType,
               NewEvent ) ) {
    // Could happen if previous clustering produced several colour
    // singlett subsystems in the event
    NewEvent.reset();
    return NewEvent;
  }

  // Build the clustered event
  Event outState = Event();
  outState.init("(hard process-modified)", particleDataPtr);
  outState.clear();

  // Copy system and incoming beam particles to outState
  for (int i = 0; i < 3; ++i)
    outState.append( NewEvent[i] );
  // Copy all the junctions one by one
  for (int i = 0; i < state.sizeJunction(); ++i)
    outState.appendJunction( state.getJunction(i) );
  // Initialise scales for new event
  outState.saveSize();
  outState.saveJunctionSize();
  outState.scale(mu);
  outState.scaleSecond(mu);
  bool radAppended = false;
  bool recAppended = false;
  int size = int(outState.size());
  // Save position of radiator in new event record
  int radPos = 0;

  // Append first incoming particle
  if ( RecBefore.mother1() == 1) {
    outState.append( RecBefore );
    recAppended = true;
  } else if ( RadBefore.mother1() == 1 ) {
    radPos = outState.append( RadBefore );
    radAppended = true;
  } else {
    // Find second incoming in input event
    int in1 = 0;
    for(int i=0; i < int(state.size()); ++i)
      if (state[i].mother1() == 1) in1 =i;
    outState.append( state[in1] );
    size++;
  }
  // Append second incoming particle
  if ( RecBefore.mother1() == 2) {
    outState.append( RecBefore );
    recAppended = true;
  } else if ( RadBefore.mother1() == 2 ) {
    radPos = outState.append( RadBefore );
    radAppended = true;
  } else {
    // Find second incoming in input event
    int in2 = 0;
    for(int i=0; i < int(state.size()); ++i)
      if (state[i].mother1() == 2) in2 =i;

    outState.append( state[in2] );
    size++;
  }

  // Append new recoiler if not done already
  if (!recAppended && !RecBefore.isFinal()) {
    recAppended = true;
    outState.append( RecBefore);
  }
  // Append new radiator if not done already
  if (!radAppended && !RadBefore.isFinal()) {
    radAppended = true;
    radPos = outState.append( RadBefore);
  }

  // Append intermediate particle
  // (careful not to append reclustered recoiler)
  // Append intermediate particle
  // (careful not to append reclustered recoiler)
  for (int i = 0; i < int(NewEvent.size()-1); ++i)
    if (NewEvent[i].status() == -22) outState.append( NewEvent[i] );
  // Append final state particles, resonances first
  for (int i = 0; i < int(NewEvent.size()-1); ++i)
    if (NewEvent[i].status() == 22) outState.append( NewEvent[i] );
  // Then start appending partons
  if (!radAppended && RadBefore.statusAbs() == 22)
    radPos = outState.append(RadBefore);
  if (!recAppended)
    outState.append(RecBefore);
  if (!radAppended && RadBefore.statusAbs() != 22)
    radPos = outState.append(RadBefore);
  // Then partons (not reclustered recoiler)
  for(int i = 0; i < int(NewEvent.size()-1); ++i)
    if ( NewEvent[i].status()  != 22
      && NewEvent[i].colType() != 0
      && NewEvent[i].isFinal())
      outState.append( NewEvent[i] );
  // Then the rest
  for(int i = 0; i < int(NewEvent.size()-1); ++i)
    if ( NewEvent[i].status() != 22
      && NewEvent[i].colType() == 0
      && NewEvent[i].isFinal() )
      outState.append( NewEvent[i]);

  // Find intermediate and respective daughters
  vector<int> PosIntermediate;
  vector<int> PosDaughter1;
  vector<int> PosDaughter2;
  for(int i=0; i < int(outState.size()); ++i)
    if (outState[i].status() == -22) {
      PosIntermediate.push_back(i);
      int d1 = outState[i].daughter1();
      int d2 = outState[i].daughter2();
      // Find daughters in output state
      int daughter1 = FindParticle( state[d1], outState);
      int daughter2 = FindParticle( state[d2], outState);
      // If both daughters found, done
      // Else put first final particle as first daughter
      // and last final particle as second daughter
      if (daughter1 > 0)
        PosDaughter1.push_back( daughter1);
      else {
        daughter1 = 0;
        while(!outState[daughter1].isFinal() ) daughter1++;
        PosDaughter1.push_back( daughter1);
      }
      if (daughter2 > 0)
        PosDaughter2.push_back( daughter2);
      else {
        daughter2 = outState.size()-1;
        while(!outState[daughter2].isFinal() ) daughter2--;
        PosDaughter2.push_back( daughter2);
      }
    }
  // Set daughters and mothers
  for(int i=0; i < int(PosIntermediate.size()); ++i) {
    outState[PosIntermediate[i]].daughters(PosDaughter1[i],PosDaughter2[i]);
    outState[PosDaughter1[i]].mother1(PosIntermediate[i]);
    outState[PosDaughter2[i]].mother1(PosIntermediate[i]);
  }

  // Find range of final state partons
  int minParFinal = int(outState.size());
  int maxParFinal = 0;
  for(int i=0; i < int(outState.size()); ++i)
    if (outState[i].mother1() == 3 && outState[i].mother2() == 4) {
      minParFinal = min(i,minParFinal);
      maxParFinal = max(i,maxParFinal);
    }

  if (minParFinal == maxParFinal) maxParFinal = 0;
  outState[3].daughters(minParFinal,maxParFinal);
  outState[4].daughters(minParFinal,maxParFinal);

  // Update event properties
  outState.saveSize();
  outState.saveJunctionSize();

  // Almost there...
  // If an intermediate coloured parton exists which was directly
  // colour connected to the radiator before the splitting, and the
  // radiator before and after the splitting had only one colour, problems
  // will arise since the colour of the radiator will be changed, whereas
  // the intermediate parton still has the old colour. In effect, this
  // means that when setting up a event for trial showering, one colour will
  // be free.
  // Hence, check for an intermediate coloured triplet resonance has been
  // colour-connected to the "old" radiator.
  // Find resonance
  int iColRes = 0;
  if ( radType == -1 && state[Rad].colType() == 1) {
      // Find resonance connected to initial colour
      for(int i=0; i < int(state.size()); ++i)
        if ( i != Rad && i != Emt && i != Rec
          && state[i].status() == -22
          && state[i].col() == state[Rad].col() )
          iColRes = i;
  } else if ( radType == -1 && state[Rad].colType() == -1) {
      // Find resonance connected to initial anticolour
      for(int i=0; i < int(state.size()); ++i)
        if ( i != Rad && i != Emt && i != Rec
          && state[i].status() == -22
          && state[i].acol() == state[Rad].acol() )
          iColRes = i;
  } else if ( radType == 1 && state[Rad].colType() == 1) {
      // Find resonance connected to final state colour
      for(int i=0; i < int(state.size()); ++i)
        if ( i != Rad && i != Emt && i != Rec
          && state[i].status() == -22
          && state[i].col() == state[Rad].col() )
          iColRes = i;
  } else if ( radType == 1 && state[Rad].colType() == -1) {
      // Find resonance connected to final state anticolour
      for(int i=0; i < int(state.size()); ++i)
        if ( i != Rad && i != Emt && i != Rec
          && state[i].status() == -22
          && state[i].acol() == state[Rad].acol() )
          iColRes = i;
  }

  if (iColRes > 0) {
    // Now find this resonance in the reclustered state
    int iColResNow = FindParticle( state[iColRes], outState);

    // Find reclustered radiator colours
    int radCol = outState[radPos].col();
    int radAcl = outState[radPos].acol();
    // Find resonance radiator colours
    int resCol = outState[iColResNow].col();
    int resAcl = outState[iColResNow].acol();
    // Check if any of the reclustered radiators colours match the resonance
    bool matchesRes =  (radCol > 0
                          && ( radCol == resCol || radCol == resAcl))
                    || (radAcl > 0
                          && ( radAcl == resCol || radAcl == resAcl));

    // If a resonance has been found, but no colours match, change
    // the colour of the resonance
    if (!matchesRes && iColResNow > 0) {
      if ( radType == -1 && outState[radPos].colType() == 1)
        outState[iColResNow].col(radCol);
      else if ( radType ==-1 && outState[radPos].colType() ==-1)
        outState[iColResNow].acol(radAcl);
      else if ( radType == 1 && outState[radPos].colType() == 1)
        outState[iColResNow].col(radCol);
      else if ( radType == 1 && outState[radPos].colType() ==-1)
        outState[iColResNow].acol(radAcl);
    }


    // If a resonance has been found, but no colours match, and the position
    // of the resonance in the event record has been changed, update the
    // radiator mother
    if (!matchesRes && iColResNow > 0 && iColRes != iColResNow)
      outState[radPos].mother1(iColResNow);

  }

  // If event is not constructed properly, return false
  if ( !validEvent(outState) ) {
    outState.reset();
    return outState;
  }

  // Remember position of reclustered radiator in state
  iReclusteredNew = radPos;

  // Done
  return outState;
}

//--------------------------------------------------------------------------

// Function to get the flavour of the radiator before the splitting
// for clustering
// IN int  : Flavour of the radiator after the splitting
//    int  : Flavour of the emitted after the splitting
// OUT int : Flavour of the radiator before the splitting

int History::getRadBeforeFlav(const int RadAfter, const int EmtAfter,
      const Event& event) {

  int type = event[RadAfter].isFinal() ? 1 :-1;
  int emtID  = event[EmtAfter].id();
  int radID  = event[RadAfter].id();
  int emtCOL = event[EmtAfter].col();
  int radCOL = event[RadAfter].col();
  int emtACL = event[EmtAfter].acol();
  int radACL = event[RadAfter].acol();

  bool colConnected = ((type == 1) && ( (emtCOL !=0 && (emtCOL ==radACL))
                                     || (emtACL !=0 && (emtACL ==radCOL)) ))
                    ||((type ==-1) && ( (emtCOL !=0 && (emtCOL ==radCOL))
                                     || (emtACL !=0 && (emtACL ==radACL)) ));
  // QCD splittings
  // Gluon radiation
  if ( emtID == 21 )
    return radID;
  // Final state gluon splitting
  if ( type == 1 && emtID == -radID && !colConnected )
    return 21;
  // Initial state s-channel gluon splitting
  if ( type ==-1 && radID == 21 )
    return -emtID;
  // Initial state t-channel gluon splitting
  if ( type ==-1 && !colConnected
    && emtID != 21 && radID != 21 && abs(emtID) < 10 && abs(radID) < 10)
    return 21;

  // SQCD splittings
  int radSign = (radID < 0) ? -1 : 1;
  int offsetL = 1000000;
  int offsetR = 2000000;
  // Gluino radiation
  if ( emtID == 1000021 ) {
    // Gluino radiation combined with quark yields squark.
    if (abs(radID) < 10 ) {
      int offset = offsetL;
      // Check if righthanded squark present. If so, make the reclustered
      // squark match. Works for squark pair production + gluino emission.
      for (int i=0; i < int(event.size()); ++i)
        if ( event[i].isFinal()
          && event[i].idAbs() < offsetR+10 && event[i].idAbs() > offsetR)
          offset = offsetR;
      return radSign*(abs(radID)+offset);
    }
    // Gluino radiation combined with squark yields quark.
    if (abs(radID) > offsetL && abs(radID) < offsetL+10 )
      return radSign*(abs(radID)-offsetL);
    if (abs(radID) > offsetR && abs(radID) < offsetR+10 )
      return radSign*(abs(radID)-offsetR);
    // Gluino radiation off gluon yields gluino.
    if (radID == 21 ) return emtID;
  }

  int emtSign = (emtID < 0) ? -1 : 1;
  // Get PDG numbering offsets.
  int emtOffset = 0;
  if ( abs(emtID) > offsetL && abs(emtID) < offsetL+10 )
    emtOffset = offsetL;
  if ( abs(emtID) > offsetR && abs(emtID) < offsetR+10 )
    emtOffset = offsetR;
  int radOffset = 0;
  if ( abs(radID) > offsetL && abs(radID) < offsetL+10 )
    radOffset = offsetL;
  if ( abs(radID) > offsetR && abs(radID) < offsetR+10 )
    radOffset = offsetR;

  // Final state gluino splitting
  if ( type == 1 && !colConnected ) {
    // Emitted squark, radiating quark.
    if ( emtOffset > 0 && radOffset == 0
      && emtSign*(abs(emtID) - emtOffset) == -radID )
      return 1000021;
    // Emitted quark, radiating squark.
    if ( emtOffset == 0 && radOffset > 0
      && emtID == -radSign*(abs(radID) - radOffset) )
      return 1000021;
  }

  // Initial state s-channel gluino splitting
  if ( type ==-1 && radID == 1000021 ) {
    // Quark entering underlying hard process.
    if ( emtOffset > 0 ) return -emtSign*(abs(emtID) - emtOffset);
    // Squark entering underlying hard process.
    else return -emtSign*(abs(emtID) + emtOffset);
  }

  // Initial state t-channel gluino splitting.
  if ( type ==-1
    && ( (abs(emtID) > offsetL && abs(emtID) < offsetL+10)
      || (abs(emtID) > offsetR && abs(emtID) < offsetR+10))
    && ( (abs(radID) > offsetL && abs(radID) < offsetL+10)
      || (abs(radID) > offsetR && abs(radID) < offsetR+10))
    && emtSign*(abs(emtID)+emtOffset) == radSign*(abs(radID) - radOffset)
    && !colConnected ) {
    return 1000021;
  }

  // Electroweak splittings splittings
  // Photon / Z radiation: Calculate invariant mass of system
  double m2final = (event[RadAfter].p()+ event[EmtAfter].p()).m2Calc();

  if ( emtID == 22 || emtID == 23 ) return radID;
  // Final state Photon splitting
  if ( type == 1 && emtID == -radID && colConnected && sqrt(m2final) <= 10. )
    return 22;
  // Final state Photon splitting
  if ( type == 1 && emtID == -radID && colConnected && sqrt(m2final)  > 10. )
    return 23;
  // Initial state s-channel photon/ Z splitting
  if ( type ==-1 && (radID == 22 || radID == 23) )
    return -emtID;
  // Initial state t-channel photon / Z splitting: Always bookkeep as photon
  if ( type ==-1 && abs(emtID) < 10 && abs(radID) < 10 && colConnected )
    return 22;

  // W+ radiation
  // Final state W+ splitting
  if ( emtID == 24 && radID < 0 ) return radID + 1;
  if ( emtID == 24 && radID > 0 ) return radID + 1;

  // W- radiation
  // Final state W- splitting
  if ( emtID ==-24 && radID < 0 ) return radID - 1;
  if ( emtID ==-24 && radID > 0 ) return radID - 1;

  // Done.
  return 0;

}

//--------------------------------------------------------------------------

// Function to properly colour-connect the radiator to the rest of
// the event, as needed during clustering
// IN  Particle& : Particle to be connected
//     Particle  : Recoiler forming a dipole with Radiator
//     Event     : event to which Radiator shall be appended
// OUT true               : Radiator could be connected to the event
//     false              : Radiator could not be connected to the
//                          event or the resulting event was
//                          non-valid

bool History::connectRadiator( Particle& Radiator, const int RadType,
                      const Particle& Recoiler, const int RecType,
                      const Event& event ) {

  // Start filling radiator colour indices with dummy values
  Radiator.cols( -1, -1 );

  // Radiator should always be colour-connected to recoiler.
  // Three cases (rad = Anti-Quark, Quark, Gluon) to be considered
  if ( Radiator.colType() == -1 ) {
    // For final state antiquark radiator, the anticolour is fixed
    // by the final / initial state recoiler colour / anticolour
    if ( RadType + RecType == 2 )
      Radiator.cols( 0, Recoiler.col());
    else if ( RadType + RecType == 0 )
      Radiator.cols( 0, Recoiler.acol());
    // For initial state antiquark radiator, the anticolour is fixed
    // by the colour of the emitted gluon (which will be the
    // leftover anticolour of a final state particle or the leftover
    // colour of an initial state particle ( = the recoiler))
    else {
      // Set colour of antiquark radiator to zero
      Radiator.col( 0 );
      for (int i = 0; i < event.size(); ++i) {
        int col = event[i].col();
        int acl = event[i].acol();

        if ( event[i].isFinal()) {
          // Search for leftover anticolour in final / initial state
          if ( acl > 0 && FindCol(acl,i,0,event,1,true) == 0
              && FindCol(acl,i,0,event,2,true) == 0 )
            Radiator.acol(event[i].acol());
        } else {
          // Search for leftover colour in initial / final state
          if ( col > 0 && FindCol(col,i,0,event,1,true) == 0
              && FindCol(col,i,0,event,2,true) == 0 )
            Radiator.acol(event[i].col());
        }
      } // end loop over particles in event record
    }

  } else if ( Radiator.colType() == 1 ) {
    // For final state quark radiator, the colour is fixed
    // by the final / initial state recoiler anticolour / colour
    if ( RadType + RecType == 2 )
      Radiator.cols( Recoiler.acol(), 0);

    else if ( RadType + RecType == 0 )
      Radiator.cols( Recoiler.col(), 0);
    // For initial state quark radiator, the colour is fixed
    // by the anticolour of the emitted gluon (which will be the
    // leftover colour of a final state particle or the leftover
    // anticolour of an initial state particle ( = the recoiler))

    else {
      // Set anticolour of quark radiator to zero
      Radiator.acol( 0 );
      for (int i = 0; i < event.size(); ++i) {
        int col = event[i].col();
        int acl = event[i].acol();

        if ( event[i].isFinal()) {
          // Search for leftover colour in final / initial state
          if ( col > 0 && FindCol(col,i,0,event,1,true) == 0
              && FindCol(col,i,0,event,2,true) == 0)
            Radiator.col(event[i].col());
        } else {
          // Search for leftover anticolour in initial / final state
          if ( acl > 0 && FindCol(acl,i,0,event,1,true) == 0
              && FindCol(acl,i,0,event,2,true) == 0)
            Radiator.col(event[i].acol());
        }
      } // end loop over particles in event record

    } // end distinction between fsr / fsr+initial recoiler / isr

  } else if ( Radiator.colType() == 2 ) {
    // For a gluon radiator, one (anticolour) colour index is defined
    // by the recoiler colour (anticolour).
    // The remaining index is chosen to match the free index in the
    // event
    // Search for leftover colour (anticolour) in the final state
    for (int i = 0; i < event.size(); ++i) {
      int col = event[i].col();
      int acl = event[i].acol();
      int iEx = i;

      if ( event[i].isFinal()) {
        if ( col > 0 && FindCol(col,iEx,0,event,1,true) == 0
          && FindCol(col,iEx,0,event,2,true) == 0) {
          if (Radiator.status() < 0 ) Radiator.col(event[i].col());
          else Radiator.acol(event[i].col());
        }
        if ( acl > 0 && FindCol(acl,iEx,0,event,2,true) == 0
          && FindCol(acl,iEx,0,event,1,true) == 0 ) {
          if (Radiator.status() < 0 )  Radiator.acol(event[i].acol());
          else Radiator.col(event[i].acol());
        }
      } else {
        if ( col > 0 && FindCol(col,iEx,0,event,1,true) == 0
          && FindCol(col,iEx,0,event,2,true) == 0) {
          if (Radiator.status() < 0 ) Radiator.acol(event[i].col());
          else Radiator.col(event[i].col());
        }
        if ( acl > 0 && (FindCol(acl,iEx,0,event,2,true) == 0
          && FindCol(acl,iEx,0,event,1,true) == 0)) {
          if (Radiator.status() < 0 ) Radiator.col(event[i].acol());
          else Radiator.acol(event[i].acol());
        }
      }
    } // end loop over particles in event record
  } // end cases of different radiator colour type

  // If either colour or anticolour has not been set, return false
  if (Radiator.col() < 0 || Radiator.acol() < 0) return false;
  // Done
  return true;
}

//--------------------------------------------------------------------------

// Function to find a colour (anticolour) index in the input event
// IN  int col       : Colour tag to be investigated
//     int iExclude1 : Identifier of first particle to be excluded
//                     from search
//     int iExclude2 : Identifier of second particle to be excluded
//                     from  search
//     Event event   : event to be searched for colour tag
//     int type      : Tag to define if col should be counted as
//                      colour (type = 1) [->find anti-colour index
//                                         contracted with col]
//                      anticolour (type = 2) [->find colour index
//                                         contracted with col]
// OUT int           : Position of particle in event record
//                     contraced with col [0 if col is free tag]

int History::FindCol(int col, int iExclude1, int iExclude2,
            const Event& event, int type, bool isHardIn) {

  bool isHard = isHardIn;
  int index = 0;

  if (isHard) {
    // Search event record for matching colour & anticolour
    for(int n = 0; n < event.size(); ++n) {
      if ( n != iExclude1 && n != iExclude2
        && event[n].colType() != 0
        &&(   event[n].status() > 0          // Check outgoing
           || event[n].status() == -21) ) {  // Check incoming
         if ( event[n].acol() == col ) {
          index = -n;
          break;
        }
        if ( event[n].col()  == col ) {
          index =  n;
          break;
        }
      }
    }
  } else {

    // Search event record for matching colour & anticolour
    for(int n = 0; n < event.size(); ++n) {
      if (  n != iExclude1 && n != iExclude2
        && event[n].colType() != 0
        &&(   event[n].status() == 43        // Check outgoing from ISR
           || event[n].status() == 51        // Check outgoing from FSR
           || event[n].status() == -41       // first initial
           || event[n].status() == -42) ) {  // second initial
        if ( event[n].acol() == col ) {
          index = -n;
          break;
        }
        if ( event[n].col()  == col ) {
          index =  n;
          break;
        }
      }
    }
  }
  // if no matching colour / anticolour has been found, return false
  if ( type == 1 && index < 0) return abs(index);
  if ( type == 2 && index > 0) return abs(index);

  return 0;
}

//--------------------------------------------------------------------------

// Function to in the input event find a particle with quantum
// numbers matching those of the input particle
// IN  Particle : Particle to be searched for
//     Event    : Event to be searched in
// OUT int      : > 0 : Position of matching particle in event
//                < 0 : No match in event

int History::FindParticle( const Particle& particle, const Event& event,
  bool checkStatus ) {

  int index = -1;

  for ( int i = int(event.size()) - 1; i > 0; --i )
    if ( event[i].id()         == particle.id()
      && event[i].colType()    == particle.colType()
      && event[i].chargeType() == particle.chargeType()
      && event[i].col()        == particle.col()
      && event[i].acol()       == particle.acol()
      && event[i].charge()     == particle.charge() ) {
      index = i;
      break;
    }

  if ( checkStatus && event[index].status() != particle.status() )
    index = -1;

  return index;
}

//--------------------------------------------------------------------------

// Function to get the colour of the radiator before the splitting
// for clustering
// IN  int   : Position of the radiator after the splitting, in the event
//     int   : Position of the emitted after the splitting, in the event
//     Event : Reference event
// OUT int   : Colour of the radiator before the splitting

int History::getRadBeforeCol(const int rad, const int emt,
      const Event& event) {

  // Save type of splitting
  int type = (event[rad].isFinal()) ? 1 :-1;
  // Get flavour of radiator after potential clustering
  int radBeforeFlav = getRadBeforeFlav(rad,emt,event);
  // Get colours of the radiator before the potential clustering
  int radBeforeCol = -1;
  // Get reconstructed gluon colours
  if (radBeforeFlav == 21) {

    // Start with quark emissions in FSR
    if (type == 1 && event[emt].id() != 21) {
      radBeforeCol = (event[rad].col()  > 0)
                   ? event[rad].col() : event[emt].col();
    // Quark emissions in ISR
    } else if (type == -1 && event[emt].id() != 21) {
      radBeforeCol = (event[rad].col()  > 0)
                   ? event[rad].col() : event[emt].acol();
    //Gluon emissions in FSR
    } else if (type == 1 && event[emt].id() == 21) {
      // If emitted is a gluon, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].col() == event[emt].acol())
                    ? event[rad].acol() : event[rad].col();
      radBeforeCol  = (event[rad].col()  == colRemove)
                    ? event[emt].col() : event[rad].col();
    //Gluon emissions in ISR
    } else if (type == -1 && event[emt].id() == 21) {
      // If emitted is a gluon, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].col() == event[emt].col())
                    ? event[rad].col() : event[rad].acol();
      radBeforeCol  = (event[rad].col()  == colRemove)
                    ? event[emt].acol() : event[rad].col();
    }

  // Get reconstructed quark colours
  } else if ( radBeforeFlav != 21 && radBeforeFlav > 0) {

    // Quark emission in FSR
    if (type == 1 && event[emt].id() != 21) {
      // If radiating is a quark, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].col() == event[emt].acol())
                    ? event[rad].acol() : 0;
      radBeforeCol  = (event[rad].col()  == colRemove)
                    ? event[emt].col() : event[rad].col();
    //Gluon emissions in FSR
    } else if (type == 1 && event[emt].id() == 21) {
      // If emitted is a gluon, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].col() == event[emt].acol())
                    ? event[rad].col() : 0;
      radBeforeCol  = (event[rad].col()  == colRemove)
                    ? event[emt].col() : event[rad].col();
    //Quark emissions in ISR
    } else if (type == -1 && event[emt].id() != 21) {
      // If emitted is a quark, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].col() == event[emt].col())
                    ? event[rad].col() : 0;
      radBeforeCol  = (event[rad].col()  == colRemove)
                    ? event[emt].acol() : event[rad].col();
    //Gluon emissions in ISR
    } else if (type == -1 && event[emt].id() == 21) {
      // If emitted is a gluon, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].col() == event[emt].col())
                    ? event[rad].col() : 0;
      radBeforeCol  = (event[rad].col()  == colRemove)
                    ? event[emt].acol() : event[rad].col();
    }
  // Other particles are assumed uncoloured
  } else {
    radBeforeCol = 0;
  }

  return radBeforeCol;

}

//--------------------------------------------------------------------------

// Function to get the anticolour of the radiator before the splitting
// for clustering
// IN  int   : Position of the radiator after the splitting, in the event
//     int   : Position of the emitted after the splitting, in the event
//     Event : Reference event
// OUT int   : Anticolour of the radiator before the splitting

int History::getRadBeforeAcol(const int rad, const int emt,
      const Event& event) {

  // Save type of splitting
  int type = (event[rad].isFinal()) ? 1 :-1;
  // Get flavour of radiator after potential clustering

  int radBeforeFlav = getRadBeforeFlav(rad,emt,event);
  // Get colours of the radiator before the potential clustering
  int radBeforeAcl = -1;
  // Get reconstructed gluon colours
  if (radBeforeFlav == 21) {

    // Start with quark emissions in FSR
    if (type == 1 && event[emt].id() != 21) {
      radBeforeAcl = (event[rad].acol() > 0)
                   ? event[rad].acol() : event[emt].acol();
    // Quark emissions in ISR
    } else if (type == -1 && event[emt].id() != 21) {
      radBeforeAcl = (event[rad].acol() > 0)
                   ? event[rad].acol() : event[emt].col();
    //Gluon emissions in FSR
    } else if (type == 1 && event[emt].id() == 21) {
      // If emitted is a gluon, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].col() == event[emt].acol())
                    ? event[rad].acol() : event[rad].col();
      radBeforeAcl  = (event[rad].acol() == colRemove)
                    ? event[emt].acol() : event[rad].acol();
    //Gluon emissions in ISR
    } else if (type == -1 && event[emt].id() == 21) {
      // If emitted is a gluon, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].col() == event[emt].col())
                    ? event[rad].col() : event[rad].acol();
      radBeforeAcl  = (event[rad].acol() == colRemove)
                    ? event[emt].col() : event[rad].acol();
    }

  // Get reconstructed anti-quark colours
  } else if ( radBeforeFlav != 21 && radBeforeFlav < 0) {

    // Antiquark emission in FSR
    if (type == 1 && event[emt].id() != 21) {
      // If radiating is a antiquark, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].col() == event[emt].acol())
                    ? event[rad].acol() : 0;
      radBeforeAcl  = (event[rad].acol()  == colRemove)
                    ? event[emt].acol() : event[rad].acol();
    //Gluon emissions in FSR
    } else if (type == 1 && event[emt].id() == 21) {
      // If emitted is a gluon, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].acol() == event[emt].col())
                    ? event[rad].acol() : 0;
      radBeforeAcl  = (event[rad].acol()  == colRemove)
                    ? event[emt].acol() : event[rad].acol();
    //Antiquark emissions in ISR
    } else if (type == -1 && event[emt].id() != 21) {
      // If emitted is an antiquark, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].acol() == event[emt].acol())
                    ? event[rad].acol() : 0;
      radBeforeAcl  = (event[rad].acol()  == colRemove)
                    ? event[emt].col() : event[rad].acol();
    //Gluon emissions in ISR
    } else if (type == -1 && event[emt].id() == 21) {
      // If emitted is a gluon, remove the repeated index, and take
      // the remaining indices as colour and anticolour
      int colRemove = (event[rad].acol() == event[emt].acol())
                    ? event[rad].acol() : 0;
      radBeforeAcl  = (event[rad].acol()  == colRemove)
                    ? event[emt].col() : event[rad].acol();
    }
  // Other particles are considered uncoloured
  } else {
    radBeforeAcl = 0;
  }

  return radBeforeAcl;

}

//--------------------------------------------------------------------------

  // Function to get the parton connected to in by a colour line
  // IN  int   : Position of parton for which partner should be found
  //     Event : Reference event
  // OUT int   : If a colour line connects the "in" parton with another
  //             parton, return the Position of the partner, else return 0

int History::getColPartner(const int in, const Event& event) {

  if (event[in].col() == 0) return 0;

  int partner = 0;
  // Try to find anticolour index first
  partner = FindCol(event[in].col(),in,0,event,1,true);
  // If no anticolour index has been found, try colour
  if (partner == 0)
   partner = FindCol(event[in].col(),in,0,event,2,true);

  return partner;

}

//--------------------------------------------------------------------------


  // Function to get the parton connected to in by an anticolour line
  // IN  int   : Position of parton for which partner should be found
  //     Event : Reference event
  // OUT int   : If an anticolour line connects the "in" parton with another
  //             parton, return the Position of the partner, else return 0

int History::getAcolPartner(const int in, const Event& event) {

  if (event[in].acol() == 0) return 0;

  int partner = 0;
  // Try to find colour index first
  partner = FindCol(event[in].acol(),in,0,event,2,true);
  // If no colour index has been found, try anticolour
  if (partner == 0)
   partner = FindCol(event[in].acol(),in,0,event,1,true);

  return partner;

}

//--------------------------------------------------------------------------

// Function to get the list of partons connected to the particle
// formed by reclusterinf emt and rad by colour and anticolour lines
// IN  int          : Position of radiator in the clustering
// IN  int          : Position of emitted in the clustering
//     Event        : Reference event
// OUT vector<int>  : List of positions of all partons that are connected
//                    to the parton that will be formed
//                    by clustering emt and rad.

vector<int> History::getReclusteredPartners(const int rad, const int emt,
  const Event& event) {

  // Save type
  int type = event[rad].isFinal() ? 1 : -1;
  // Get reclustered colours
  int radBeforeCol = getRadBeforeCol(rad, emt, event);
  int radBeforeAcl = getRadBeforeAcol(rad, emt, event);
  // Declare output
  vector<int> partners;


  // Start with FSR clusterings
  if (type == 1) {

    for(int i=0; i < int(event.size()); ++i) {
      // Check all initial state partons
      if ( i != emt && i != rad
        && event[i].status() == -21
        && event[i].col() > 0
        && event[i].col() == radBeforeCol)
          partners.push_back(i);
      // Check all final state partons
      if ( i != emt && i != rad
        && event[i].isFinal()
        && event[i].acol() > 0
        && event[i].acol() == radBeforeCol)
          partners.push_back(i);
      // Check all initial state partons
      if ( i != emt && i != rad
        && event[i].status() == -21
        && event[i].acol() > 0
        && event[i].acol() == radBeforeAcl)
          partners.push_back(i);
      // Check all final state partons
      if ( i != emt && i != rad
        && event[i].isFinal()
        && event[i].col() > 0
        && event[i].col() == radBeforeAcl)
          partners.push_back(i);
    }
  // Start with ISR clusterings
  } else {

    for(int i=0; i < int(event.size()); ++i) {
      // Check all initial state partons
      if ( i != emt && i != rad
        && event[i].status() == -21
        && event[i].acol() > 0
        && event[i].acol() == radBeforeCol)
          partners.push_back(i);
      // Check all final state partons
      if ( i != emt && i != rad
        && event[i].isFinal()
        && event[i].col() > 0
        && event[i].col() == radBeforeCol)
          partners.push_back(i);
      // Check all initial state partons
      if ( i != emt && i != rad
        && event[i].status() == -21
        && event[i].col() > 0
        && event[i].col() == radBeforeAcl)
          partners.push_back(i);
      // Check all final state partons
      if ( i != emt && i != rad
        && event[i].isFinal()
        && event[i].acol() > 0
        && event[i].acol() == radBeforeAcl)
          partners.push_back(i);
    }

  }
  // Done
  return partners;
}

//--------------------------------------------------------------------------

// Function to extract a chain of colour-connected partons in
// the event
// IN     int          : Type of parton from which to start extracting a
//                       parton chain. If the starting point is a quark
//                       i.e. flavType = 1, a chain of partons that are
//                       consecutively connected by colour-lines will be
//                       extracted. If the starting point is an antiquark
//                       i.e. flavType =-1, a chain of partons that are
//                       consecutively connected by anticolour-lines
//                       will be extracted.
// IN      int         : Position of the parton from which a
//                       colour-connected chain should be derived
// IN      Event       : Refernence event
// IN/OUT  vector<int> : Partons that should be excluded from the search.
// OUT     vector<int> : Positions of partons along the chain
// OUT     bool        : Found singlet / did not find singlet

bool History::getColSinglet( const int flavType, const int iParton,
  const Event& event, vector<int>& exclude, vector<int>& colSinglet) {

  // If no possible flavour to start from has been found
  if (iParton < 0) return false;

  // If no further partner has been found in a previous iteration,
  // and the whole final state has been excluded, we're done
  if (iParton == 0) {

    // Count number of final state partons
    int nFinal = 0;
    for(int i=0; i < int(event.size()); ++i)
      if ( event[i].isFinal() && event[i].colType() != 0)
        nFinal++;

    // Get number of initial state partons in the list of
    // excluded partons
    int nExclude = int(exclude.size());
    int nInitExclude = 0;
    if (!event[exclude[2]].isFinal())
      nInitExclude++;
    if (!event[exclude[3]].isFinal())
      nInitExclude++;

    // If the whole final state has been considered, return
    if (nFinal == nExclude - nInitExclude)
      return true;
    else
      return false;

  }

  // Declare colour partner
  int colP = 0;
  // Save the colour partner
  colSinglet.push_back(iParton);
  // Remove the partner from the list
  exclude.push_back(iParton);
  // When starting out from a quark line, follow the colour lines
  if (flavType == 1)
    colP = getColPartner(iParton,event);
  // When starting out from an antiquark line, follow the anticolour lines
  else
    colP = getAcolPartner(iParton,event);

  // Do not count excluded partons twice
  for(int i = 0; i < int(exclude.size()); ++i)
    if (colP == exclude[i])
      return true;

  // Recurse
  return getColSinglet(flavType,colP,event,exclude,colSinglet);

}

//--------------------------------------------------------------------------

// Function to check that a set of partons forms a colour singlet
// IN  Event       : Reference event
// IN  vector<int> : Positions of the partons in the set
// OUT bool        : Is a colour singlet / is not

bool History::isColSinglet( const Event& event,
  vector<int> system ) {

  // Check if system forms a colour singlet
  for(int i=0; i < int(system.size()); ++i ) {
    // Match quark and gluon colours
    if ( system[i] > 0
      && (event[system[i]].colType() == 1
       || event[system[i]].colType() == 2) ) {
      for(int j=0; j < int(system.size()); ++j)
        // If flavour matches, remove both partons and continue
        if ( system[i] > 0
          && system[j] > 0
          && event[system[i]].col() == event[system[j]].acol()) {
          // Remove index and break
          system[i] = 0;
          system[j] = 0;
          break;
        }
    }
    // Match antiquark and gluon anticolours
    if ( system[i] > 0
      && (event[system[i]].colType() == -1
       || event[system[i]].colType() == 2) ) {
      for(int j=0; j < int(system.size()); ++j)
        // If flavour matches, remove both partons and continue
        if ( system[i] > 0
          && system[j] > 0
          && event[system[i]].acol() == event[system[j]].col()) {
          // Remove index and break
          system[i] = 0;
          system[j] = 0;
          break;
        }
    }

  }

  // The system is a colour singlet if for all colours,
  // an anticolour was found
  bool isColSing = true;
  for(int i=0; i < int(system.size()); ++i)
    if ( system[i] != 0 )
      isColSing = false;

  // Return
  return isColSing;


}

//--------------------------------------------------------------------------

// Function to check that a set of partons forms a flavour singlet
// IN  Event       : Reference event
// IN  vector<int> : Positions of the partons in the set
// IN  int         : Flavour of all the quarks in the set, if
//                   all quarks in a set should have a fixed flavour
// OUT bool        : Is a flavour singlet / is not

bool History::isFlavSinglet( const Event& event,
  vector<int> system, int flav) {

  // If a decoupled colour singlet has been found, check if this is also
  // a flavour singlet
  // Check that each quark matches an antiquark
  for(int i=0; i < int(system.size()); ++i)
    if ( system[i] > 0 ) {
      for(int j=0; j < int(system.size()); ++j) {
        // If flavour of outgoing partons matches,
        // remove both partons and continue.
        // Skip all bosons
        if ( event[i].idAbs() != 21
          && event[i].idAbs() != 22
          && event[i].idAbs() != 23
          && event[i].idAbs() != 24
          && system[i] > 0
          && system[j] > 0
          && event[system[i]].isFinal()
          && event[system[j]].isFinal()
          && event[system[i]].id() == -1*event[system[j]].id()) {
          // If we want to check if only one flavour of quarks
          // exists
          if (abs(flav) > 0 && event[system[i]].idAbs() != flav)
            return false;
          // Remove index and break
          system[i] = 0;
          system[j] = 0;
          break;
        }
        // If flavour of outgoing and incoming partons match,
        // remove both partons and continue.
        // Skip all bosons
        if ( event[i].idAbs() != 21
          && event[i].idAbs() != 22
          && event[i].idAbs() != 23
          && event[i].idAbs() != 24
          && system[i] > 0
          && system[j] > 0
          && ( ( !event[system[i]].isFinal() && event[system[j]].isFinal())
             ||( !event[system[j]].isFinal() && event[system[i]].isFinal()) )
          && event[system[i]].id() == event[system[j]].id()) {
          // If we want to check if only one flavour of quarks
          // exists
          if (abs(flav) > 0 && event[system[i]].idAbs() != flav)
            return false;
          // Remove index and break
          system[i] = 0;
          system[j] = 0;
          break;
        }

      }
    }

  // The colour singlet is a flavour singlet if for all quarks,
  // an antiquark was found
  bool isFlavSing = true;
  for(int i=0; i < int(system.size()); ++i)
    if ( system[i] != 0 )
      isFlavSing = false;

  // Return
  return isFlavSing;

}

//--------------------------------------------------------------------------

// Function to check if rad,emt,rec triple is allowed for clustering
// IN int rad,emt,rec : Positions (in event record) of the three
//                      particles considered for clustering
//    Event event     : Reference event

bool History::allowedClustering( int rad, int emt, int rec, int partner,
                const Event& event ) {

  // Declare output
  bool allowed = true;

  // CONSTRUCT SOME PROPERTIES FOR LATER INVESTIGATION

  // Check if the triple forms a colour singlett
  bool isSing = isSinglett(rad,emt,partner,event);
  int type = (event[rad].isFinal()) ? 1 :-1;
  // Get flavour of radiator after potential clustering
  int radBeforeFlav = getRadBeforeFlav(rad,emt,event);
  // Get colours of the radiator before the potential clustering
  int radBeforeCol = getRadBeforeCol(rad,emt,event);
  int radBeforeAcl = getRadBeforeAcol(rad,emt,event);
  // Get colour partner of reclustered parton
  vector<int> radBeforeColP = getReclusteredPartners(rad, emt, event);

  // Count coloured partons in hard process
  int nPartonInHard = 0;
  for(int i=0; i < int(event.size()); ++i)
    // Check all final state partons
    if ( event[i].isFinal()
      && event[i].colType() != 0
      && mergingHooksPtr->hardProcess.matchesAnyOutgoing(i, event) )
      nPartonInHard++;

  // Count coloured final state partons in event, excluding
  // rad, rec, emt and hard process
  int nPartons = 0;
  for(int i=0; i < int(event.size()); ++i)
    // Check all final state partons
    if ( i!=emt && i!=rad && i!=rec
      &&  event[i].isFinal()
      &&  event[i].colType() != 0
      && !mergingHooksPtr->hardProcess.matchesAnyOutgoing(i, event) )
      nPartons++;

  // Count number of initial state partons
  int nInitialPartons = 0;
  for(int i=0; i < int(event.size()); ++i)
    if ( event[i].status() == -21
      && event[i].colType() != 0 )
      nInitialPartons++;

  // Get number of non-charged final state particles
  int nFinalEW = 0;
  for(int i=0; i < int(event.size()); ++i)
    if ( event[i].isFinal()
      &&(  event[i].id() == 22
        || event[i].id() == 23
        || event[i].id() == 24
        ||(event[i].idAbs() > 10 && event[i].idAbs() < 20)
        ||(event[i].idAbs() > 10 && event[i].idAbs() < 20)
        ||(event[i].idAbs() > 1000010 && event[i].idAbs() < 1000020)
        ||(event[i].idAbs() > 2000010 && event[i].idAbs() < 2000020) ))
      nFinalEW++;

  // Check if event after potential clustering contains an even
  // number of quarks and/or antiquarks
  // (otherwise no electroweak vertex could be formed!)
  // Get number of final quarks
  int nFinalQuark = 0;
  // Get number of excluded final state quarks as well
  int nFinalQuarkExc = 0;
  for(int i=0; i < int(event.size()); ++i) {
    if (i !=rad && i != emt && i != rec) {
      if (event[i].isFinal() && abs(event[i].colType()) == 1 ) {
        if ( !mergingHooksPtr->hardProcess.matchesAnyOutgoing(i,event) )
          nFinalQuark++;
        else
          nFinalQuarkExc++;
      }
    }
  }

  // Add recoiler to number of final quarks
  if (event[rec].isFinal() && event[rec].isQuark()) nFinalQuark++;
  // Add radiator after clustering to number of final quarks
  if (event[rad].isFinal() && abs(radBeforeFlav) < 10) nFinalQuark++;

  // Get number of initial quarks
  int nInitialQuark = 0;
  if (type == 1) {
    if (event[rec].isFinal()) {
      if (event[3].isQuark()) nInitialQuark++;
      if (event[4].isQuark()) nInitialQuark++;
    } else {
      int iOtherIn = (rec == 3) ? 4 : 3;
      if (event[rec].isQuark()) nInitialQuark++;
      if (event[iOtherIn].isQuark()) nInitialQuark++;
    }
  } else {
    // Add recoiler to number of initial quarks
    if (event[rec].isQuark()) nInitialQuark++;
    // Add radiator after clustering to number of initial quarks
    if (abs(radBeforeFlav) < 10) nInitialQuark++;
  }

  // BEGIN CHECKING THE CLUSTERING

  // Do not allow clusterings that lead to a disallowed proton content.
  int proton[] = {1,2,3,4,5,21,22,23,24};
  bool isInProton = false;
  for(int i=0; i < 9; ++i)
    if (abs(radBeforeFlav) == proton[i]) isInProton = true;
  if (type == -1 && !isInProton) return false;

  // Check if colour is conserved
  vector<int> unmatchedCol;
  vector<int> unmatchedAcl;
  // Check all unmatched colours
  for ( int i = 0; i < event.size(); ++i)
    if ( i != emt && i != rad
      && (event[i].isFinal() || event[i].status() == -21)
      && event[i].colType() != 0 ) {

      int colP = getColPartner(i,event);
      int aclP = getAcolPartner(i,event);

      if (event[i].col() > 0
        && (colP == emt || colP == rad || colP == 0) )
        unmatchedCol.push_back(i);
      if (event[i].acol() > 0
        && (aclP == emt || aclP == rad || aclP == 0) )
        unmatchedAcl.push_back(i);

    }

  // If more than one colour or more than one anticolour are unmatched,
  // there is no way to make this clustering work
  if (int(unmatchedCol.size()) + int(unmatchedAcl.size()) > 2)
    return false;

  // If triple forms colour singlett, check that resulting state
  // matches hard core process
  if (isSing)
    allowed = false;
  if (isSing && (abs(radBeforeFlav)<10 && event[rec].isQuark()) )
    allowed = true;

  // Never recluster any outgoing partons of the core V -> qqbar' splitting!
  if ( mergingHooksPtr->hardProcess.matchesAnyOutgoing(emt,event) ) {
    // Check if any other particle could replace "emt" as part of the candidate
    // core process. If so, replace emt with the new candidate and allow the
    // clustering.
    bool canReplace = mergingHooksPtr->hardProcess.findOtherCandidates(emt,
                        event, true);
    if (canReplace) allowed = true;
    else allowed = false;
  }

  // Never allow clustering of any outgoing partons of the hard process
  // which would change the flavour of one of the hard process partons!
  if ( mergingHooksPtr->hardProcess.matchesAnyOutgoing(rad,event)
      && event[rad].id() != radBeforeFlav )
    allowed = false;

  // If only gluons in initial state and no quarks in final state,
  // reject (no electroweak vertex can be formed)
  if (nFinalEW != 0 && nInitialQuark == 0
    && nFinalQuark == 0 && nFinalQuarkExc == 0)
    allowed = false;

  if ( (nInitialQuark + nFinalQuark + nFinalQuarkExc)%2 != 0 )
    allowed = false;

  // Disallow final state splittings that lead to a purely gluonic final
  // state, while having a completely colour-connected initial state.
  // This means that the clustering is discarded if it does not lead to the
  // t-channel gluon needed to connect the final state to a qq~ initial state.
  // Here, partons excluded from clustering are not counted as possible
  // partners to form a t-channel gluon
  if (event[3].col() == event[4].acol()
    && event[3].acol() == event[4].col()
    && nFinalQuark == 0){
    // Careful if rad and rec are the only quarks in the final state, but
    // were both excluded from the list of final state quarks.
    int nTripletts = abs(event[rec].colType())
                   + abs(particleDataPtr->colType(radBeforeFlav));
    if (event[3].isGluon())                            allowed = false;
    else if (nTripletts != 2 && nFinalQuarkExc%2 == 0) allowed = false;
  }

  // No problems with gluon radiation
  if (event[emt].id() == 21)
    return allowed;

  // No problems with gluino radiation
  if (event[emt].id() == 1000021)
    return allowed;

  // Save all hard process candidates
  vector<int> outgoingParticles;
  int nOut1 = int(mergingHooksPtr->hardProcess.PosOutgoing1.size());
  for ( int i=0; i < nOut1;  ++i ) {
    int iPos = mergingHooksPtr->hardProcess.PosOutgoing1[i];
    outgoingParticles.push_back(
                      mergingHooksPtr->hardProcess.state[iPos].id() );
  }
  int nOut2 = int(mergingHooksPtr->hardProcess.PosOutgoing2.size());
  for ( int i=0; i < nOut2; ++i ) {
    int iPos = mergingHooksPtr->hardProcess.PosOutgoing2[i];
    outgoingParticles.push_back(
                      mergingHooksPtr->hardProcess.state[iPos].id() );
  }

  // Start more involved checks. g -> q_1 qbar_1 splittings are
  // particularly problematic if more than one quark of the emitted
  // flavour is present.
  // Count number of initial quarks of radiator or emitted flavour
  vector<int> iInQuarkFlav;
  for(int i=0; i < int(event.size()); ++i)
    // Check all initial state partons
    if ( i != emt && i != rad
      && event[i].status() == -21
      && event[i].idAbs() == event[emt].idAbs() )
      iInQuarkFlav.push_back(i);

  // Count number of final quarks of radiator or emitted flavour
  vector<int> iOutQuarkFlav;
  for(int i=0; i < int(event.size()); ++i)
  // Check all final state partons
  if ( i != emt && i != rad
    && event[i].isFinal()
    && event[i].idAbs() == event[emt].idAbs() ) {

    // Loop through final state hard particles. If one matches, remove the
    // matching one, and do not count.
    bool matchOut = false;
    for (int j = 0; j < int(outgoingParticles.size()); ++j)
    if ( event[i].idAbs() == abs(outgoingParticles[j])) {
      matchOut = true;
      outgoingParticles[j] = 99;
    }
    if (!matchOut) iOutQuarkFlav.push_back(i);

  }

  // Save number of potentially dangerous quarks
  int nInQuarkFlav  = int(iInQuarkFlav.size());
  int nOutQuarkFlav = int(iOutQuarkFlav.size());

  // Easiest problem 0:
  // Radiator before splitting exactly matches the partner
  // after the splitting
  if ( event[partner].isFinal()
    && event[partner].id()   == 21
    && radBeforeFlav         == 21
    && event[partner].col()  == radBeforeCol
    && event[partner].acol() == radBeforeAcl)
    return false;

  // If there are no ambiguities in qqbar pairs, return
  if (nInQuarkFlav + nOutQuarkFlav == 0)
    return allowed;

  // Save all quarks and gluons that will not change colour
  vector<int> gluon;
  vector<int> quark;
  vector<int> antiq;
  vector<int> partons;
  for(int i=0; i < int(event.size()); ++i)
    // Check initial and final state partons
    if ( i!=emt && i!=rad
      && event[i].colType() != 0
      && (event[i].isFinal() || event[i].status() == -21) ) {
      // Save index
      partons.push_back(i);
      // Split into components
      if (event[i].colType() == 2)
        gluon.push_back(i);
      else if (event[i].colType() ==  1)
        quark.push_back(i);
      else if (event[i].colType() == -1)
        antiq.push_back(i);
    }

  // We split up the test of the g->qq splitting into final state
  // and initial state problems
  bool isFSRg2qq = ((type == 1) && (event[rad].id() == -1*event[emt].id()) );
  bool isISRg2qq = ((type ==-1) && (event[rad].id() ==    event[emt].id()) );

  // First check general things about colour connections
  // Check that clustering does not produce a gluon that is exactly
  // matched in the final state, or does not have any colour connections
  if ( (isFSRg2qq || isISRg2qq)
    && int(quark.size()) + int(antiq.size())
     + int(gluon.size()) > nPartonInHard ) {

      vector<int> colours;
      vector<int> anticolours;
      // Add the colour and anticolour of the gluon before the emission
      // to the list, bookkeep initial colour as final anticolour, and
      // initial anticolour as final colour
      if (type == 1) {
        colours.push_back(radBeforeCol);
        anticolours.push_back(radBeforeAcl);
      } else {
        colours.push_back(radBeforeAcl);
        anticolours.push_back(radBeforeCol);
      }
      // Now store gluon colours and anticolours.
      for(int i=0; i < int(gluon.size()); ++i)
        if (event[gluon[i]].isFinal()) {
          colours.push_back(event[gluon[i]].col());
          anticolours.push_back(event[gluon[i]].acol());
        } else {
          colours.push_back(event[gluon[i]].acol());
          anticolours.push_back(event[gluon[i]].col());
        }

      // Loop through colours and check if any match with
      // anticolours. If colour matches, remove from list
      for(int i=0; i < int(colours.size()); ++i)
        for(int j=0; j < int(anticolours.size()); ++j)
          if (colours[i] > 0 && anticolours[j] > 0
            && colours[i] == anticolours[j]) {
            colours[i] = 0;
            anticolours[j] = 0;
          }


      // If all gluon anticolours and all colours matched, disallow
      // the clustering
      bool allMatched = true;
      for(int i=0; i < int(colours.size()); ++i)
        if (colours[i] != 0)
          allMatched = false;
      for(int i=0; i < int(anticolours.size()); ++i)
        if (anticolours[i] != 0)
          allMatched = false;

      if (allMatched)
        return false;

      // Now add the colours of the hard process, and check if all
      // colours match.
      for(int i=0; i < int(quark.size()); ++i)
        if ( event[quark[i]].isFinal()
        && mergingHooksPtr->hardProcess.matchesAnyOutgoing(quark[i], event) )
          colours.push_back(event[quark[i]].col());

      for(int i=0; i < int(antiq.size()); ++i)
        if ( event[antiq[i]].isFinal()
        && mergingHooksPtr->hardProcess.matchesAnyOutgoing(antiq[i], event) )
          anticolours.push_back(event[antiq[i]].acol());

      // Loop through colours again and check if any match with
      // anticolours. If colour matches, remove from list
      for(int i=0; i < int(colours.size()); ++i)

        for(int j=0; j < int(anticolours.size()); ++j)
          if (colours[i] > 0 && anticolours[j] > 0
            && colours[i] == anticolours[j]) {
            colours[i] = 0;
            anticolours[j] = 0;
          }

      // Check if clustering would produce the hard process
      int nNotInHard = 0;
      for ( int i=0; i < int(quark.size()); ++i )
        if ( !mergingHooksPtr->hardProcess.matchesAnyOutgoing( quark[i],
              event) )
          nNotInHard++;
      for ( int i=0; i < int(antiq.size()); ++i )
        if ( !mergingHooksPtr->hardProcess.matchesAnyOutgoing( antiq[i],
              event) )
          nNotInHard++;
      for(int i=0; i < int(gluon.size()); ++i)
        if ( event[gluon[i]].isFinal() )
          nNotInHard++;
      if ( type == 1 )
          nNotInHard++;

      // If all colours are matched now, and since we have more quarks than
      // present in the hard process, disallow the clustering
      allMatched = true;
      for(int i=0; i < int(colours.size()); ++i)
        if (colours[i] != 0)
          allMatched = false;
      for(int i=0; i < int(anticolours.size()); ++i)
        if (anticolours[i] != 0)
          allMatched = false;

      if (allMatched && nNotInHard > 0)
        return false;

  }

  // FSR PROBLEMS

  if (isFSRg2qq && nInQuarkFlav + nOutQuarkFlav > 0) {

    // Easiest problem 1:
    // RECLUSTERED FINAL STATE GLUON MATCHES INITIAL STATE GLUON
    for(int i=0; i < int(gluon.size()); ++i) {
      if (!event[gluon[i]].isFinal()
        && event[gluon[i]].col()  == radBeforeCol
        && event[gluon[i]].acol() == radBeforeAcl)
        return false;
    }

    // Easiest problem 2:
    // RECLUSTERED FINAL STATE GLUON MATCHES FINAL STATE GLUON
    for(int i=0; i < int(gluon.size()); ++i) {
      if (event[gluon[i]].isFinal()
        && event[gluon[i]].col()  == radBeforeAcl
        && event[gluon[i]].acol() == radBeforeCol)
        return false;
    }

    // Easiest problem 3:
    // RECLUSTERED FINAL STATE GLUON MATCHES FINAL STATE Q-QBAR PAIR
    if ( int(radBeforeColP.size()) == 2
      && event[radBeforeColP[0]].isFinal()
      && event[radBeforeColP[1]].isFinal()
      && event[radBeforeColP[0]].id() == -1*event[radBeforeColP[1]].id() ) {

      // This clustering is allowed if there is no colour in the
      // initial state
      if (nInitialPartons > 0)
        return false;
    }

    // Next-to-easiest problem 1:
    // RECLUSTERED FINAL STATE GLUON MATCHES ONE FINAL STARE Q_1
    // AND ONE INITIAL STATE Q_1
    if ( int(radBeforeColP.size()) == 2
      && ((  event[radBeforeColP[0]].status() == -21
          && event[radBeforeColP[1]].isFinal())
        ||(  event[radBeforeColP[0]].isFinal()
          && event[radBeforeColP[1]].status() == -21))
      && event[radBeforeColP[0]].id() == event[radBeforeColP[1]].id() ) {

      // In principle, clustering this splitting can disconnect
      // the colour lines of a graph. However, the colours can be connected
      // again if a final or initial partons of the correct flavour exists.

      // Check which of the partners are final / initial
      int incoming = (event[radBeforeColP[0]].isFinal())
                   ? radBeforeColP[1] : radBeforeColP[0];
      int outgoing = (event[radBeforeColP[0]].isFinal())
                   ? radBeforeColP[0] : radBeforeColP[1];

      // Loop through event to find "recovery partons"
      bool clusPossible = false;
      for(int i=0; i < int(event.size()); ++i)
        if (  i != emt && i != rad
          &&  i != incoming && i != outgoing
          && !mergingHooksPtr->hardProcess.matchesAnyOutgoing(i,event) ) {
          // Check if an incoming parton matches
          if ( event[i].status() == -21
            && (event[i].id() ==    event[outgoing].id()
              ||event[i].id() == -1*event[incoming].id()) )
          clusPossible = true;
          // Check if a final parton matches
          if ( event[i].isFinal()
            && (event[i].id() == -1*event[outgoing].id()
              ||event[i].id() ==    event[incoming].id()) )
          clusPossible = true;
        }

      // There can be a further complication: If e.g. in
      // t-channel photon exchange topologies, both incoming
      // partons are quarks, and form colour singlets with any
      // number of final state partons, at least try to
      // recluster as much as possible.
      // For this, check if the incoming parton
      // connected to the radiator is connected to a
      // colour and flavour singlet
      vector<int> excludeIn1;
      for(int i=0; i < 4; ++i)
        excludeIn1.push_back(0);
      vector<int> colSingletIn1;
      int flavIn1Type = (event[incoming].id() > 0) ? 1 : -1;
      // Try finding colour singlets
      bool isColSingIn1  = getColSinglet(flavIn1Type,incoming,event,
                             excludeIn1,colSingletIn1);
      // Check if colour singlet also is a flavour singlet
      bool isFlavSingIn1 = isFlavSinglet(event,colSingletIn1);

      // Check if the incoming parton not
      // connected to the radiator is connected to a
      // colour and flavour singlet
      int incoming2 = (incoming == 3) ? 4 : 3;
      vector<int> excludeIn2;
      for(int i=0; i < 4; ++i)
        excludeIn2.push_back(0);
      vector<int> colSingletIn2;
      int flavIn2Type = (event[incoming2].id() > 0) ? 1 : -1;
      // Try finding colour singlets
      bool isColSingIn2  = getColSinglet(flavIn2Type,incoming2,event,
                             excludeIn2,colSingletIn2);
      // Check if colour singlet also is a flavour singlet
      bool isFlavSingIn2 = isFlavSinglet(event,colSingletIn2);

      // If no "recovery clustering" is possible, reject clustering
      if (!clusPossible
        && (!isColSingIn1 || !isFlavSingIn1
         || !isColSingIn2 || !isFlavSingIn2))
        return false;

    }

    // Next-to-easiest problem 2:
    // FINAL STATE Q-QBAR CLUSTERING DISCONNECTS SINGLETT SUBSYSTEM WITH
    // FINAL STATE Q-QBAR PAIR FROM GRAPH

    // Prepare to check for colour singlet combinations of final state quarks
    // Start by building a list of partons to exclude when checking for
    // colour singlet combinations
    int flav = event[emt].id();
    vector<int> exclude;
    exclude.push_back(emt);
    exclude.push_back(rad);
    exclude.push_back(radBeforeColP[0]);
    exclude.push_back(radBeforeColP[1]);
    vector<int> colSinglet;
    // Now find parton from which to start checking colour singlets
    int iOther = -1;
    // Loop through event to find a parton of correct flavour
    for(int i=0; i < int(event.size()); ++i)
      // Check final state for parton equalling emitted flavour.
      // Exclude the colour system coupled to the clustering
      if ( i != emt
        && i != rad
        && i != radBeforeColP[0]
        && i != radBeforeColP[1]
        && event[i].isFinal() ) {
        // Stop if one parton of the correct flavour is found
        if (event[i].id() == flav) {
          iOther = i;
          break;
        }
      }
    // Save the type of flavour
    int flavType = (event[iOther].id() > 0) ? 1 : -1;
    // Try finding colour singlets
    bool isColSing = getColSinglet(flavType,iOther,event,exclude,colSinglet);
    // Check if colour singlet also is a flavour singlet
    bool isFlavSing = isFlavSinglet(event,colSinglet);

    // Check if the colour singlet is precisely contained in the hard process.
    // If so, then we're safe to recluster.
    bool isHardSys = true;
    for(int i=0; i < int(colSinglet.size()); ++i)
      isHardSys =
         mergingHooksPtr->hardProcess.matchesAnyOutgoing(colSinglet[i], event);

    // Nearly there...
    // If the decoupled colour singlet system is NOT contained in the hard
    // process, we need to check the whole final state.
    if (isColSing && isFlavSing && !isHardSys) {

      // In a final check, ensure that the final state does not only
      // consist of colour singlets that are also flavour singlets
      // of the identical (!) flavours
      // Loop through event and save all final state partons
      vector<int> allFinal;
      for(int i=0; i < int(event.size()); ++i)
        if ( event[i].isFinal() )
          allFinal.push_back(i);

      // Check if all final partons form a colour singlet
      bool isFullColSing  = isColSinglet(event,allFinal);
      // Check if all final partons form a flavour singlet
      bool isFullFlavSing = isFlavSinglet(event,allFinal,flav);

      // If all final quarks are of identical flavour,
      // no possible clustering should be discriminated.
      // Otherwise, disallow
      if (!isFullColSing || !isFullFlavSing)
        return false;
    }
  }

  // ISR PROBLEMS

  if (isISRg2qq && nInQuarkFlav + nOutQuarkFlav > 0) {

    // Easiest problem 1:
    // RECLUSTERED INITIAL STATE GLUON MATCHES FINAL STATE GLUON
    for(int i=0; i < int(gluon.size()); ++i) {
      if (event[gluon[i]].isFinal()
        && event[gluon[i]].col()  == radBeforeCol
        && event[gluon[i]].acol() == radBeforeAcl)
        return false;
    }

    // Easiest problem 2:
    // RECLUSTERED INITIAL STATE GLUON MATCHES INITIAL STATE GLUON
    for(int i=0; i < int(gluon.size()); ++i) {
      if (event[gluon[i]].status() == -21
        && event[gluon[i]].acol()  == radBeforeCol
        && event[gluon[i]].col() == radBeforeAcl)
        return false;
    }

    // Next-to-easiest problem 1:
    // RECLUSTERED INITIAL STATE GLUON MATCHES FINAL STATE Q-QBAR PAIR
    if ( int(radBeforeColP.size()) == 2
      && event[radBeforeColP[0]].isFinal()
      && event[radBeforeColP[1]].isFinal()
      && event[radBeforeColP[0]].id() == -1*event[radBeforeColP[1]].id() ) {

      // In principle, clustering this splitting can disconnect
      // the colour lines of a graph. However, the colours can be connected
      // again if final state partons of the correct (anti)flavour, or
      // initial state partons of the correct flavour exist
      // Loop through event to check
      bool clusPossible = false;
      for(int i=0; i < int(event.size()); ++i)
        if ( i != emt && i != rad
          && i != radBeforeColP[0]
          && i != radBeforeColP[1]
          && !mergingHooksPtr->hardProcess.matchesAnyOutgoing(i,event) ) {
          if (event[i].status() == -21
            && ( event[radBeforeColP[0]].id() == event[i].id()
              || event[radBeforeColP[1]].id() == event[i].id() ))

            clusPossible = true;
          if (event[i].isFinal()
            && ( event[radBeforeColP[0]].id() == -1*event[i].id()
              || event[radBeforeColP[1]].id() == -1*event[i].id() ))
            clusPossible = true;
        }

      // There can be a further complication: If e.g. in
      // t-channel photon exchange topologies, both incoming
      // partons are quarks, and form colour singlets with any
      // number of final state partons, at least try to
      // recluster as much as possible.
      // For this, check if the incoming parton
      // connected to the radiator is connected to a
      // colour and flavour singlet
      int incoming1 = 3;
      vector<int> excludeIn1;
      for(int i=0; i < 4; ++i)
        excludeIn1.push_back(0);
      vector<int> colSingletIn1;
      int flavIn1Type = (event[incoming1].id() > 0) ? 1 : -1;
      // Try finding colour singlets
      bool isColSingIn1  = getColSinglet(flavIn1Type,incoming1,event,
                             excludeIn1,colSingletIn1);
      // Check if colour singlet also is a flavour singlet
      bool isFlavSingIn1 = isFlavSinglet(event,colSingletIn1);

      // Check if the incoming parton not
      // connected to the radiator is connected to a
      // colour and flavour singlet
      int incoming2 = 4;
      vector<int> excludeIn2;
      for(int i=0; i < 4; ++i)
        excludeIn2.push_back(0);
      vector<int> colSingletIn2;
      int flavIn2Type = (event[incoming2].id() > 0) ? 1 : -1;
      // Try finding colour singlets
      bool isColSingIn2  = getColSinglet(flavIn2Type,incoming2,event,
                             excludeIn2,colSingletIn2);
      // Check if colour singlet also is a flavour singlet
      bool isFlavSingIn2 = isFlavSinglet(event,colSingletIn2);

      // If no "recovery clustering" is possible, reject clustering
      if (!clusPossible
        && (!isColSingIn1 || !isFlavSingIn1
         || !isColSingIn2 || !isFlavSingIn2))
        return false;

    }

  }

  // Done
  return allowed;
}

//--------------------------------------------------------------------------

// Function to check if rad,emt,rec triple is results in
// colour singlet radBefore+recBefore
// IN int rad,emt,rec : Positions (in event record) of the three
//                      particles considered for clustering
//    Event event     : Reference event

bool History::isSinglett( int rad, int emt, int rec, const Event& event ) {

  int radCol = event[rad].col();
  int emtCol = event[emt].col();
  int recCol = event[rec].col();
  int radAcl = event[rad].acol();
  int emtAcl = event[emt].acol();
  int recAcl = event[rec].acol();
  int recType = event[rec].isFinal() ? 1 : -1;

  bool isSing = false;

  if ( ( recType == -1
       && radCol + emtCol == recCol && radAcl + emtAcl == recAcl)
    ||( recType == 1
       && radCol + emtCol == recAcl && radAcl + emtAcl == recCol) )
    isSing = true;

  return isSing;

}

//--------------------------------------------------------------------------

// Function to check if event is sensibly constructed: Meaning
// that all colour indices are contracted and that the charge in
// initial and final states matches
// IN  event : event to be checked
// OUT TRUE  : event is properly construced
//     FALSE : event not valid

bool History::validEvent( const Event& event ) {

  // Check if event is coloured
  bool validColour = true;
  for ( int i = 0; i < event.size(); ++i)
   // Check colour of quarks
   if ( event[i].isFinal() && event[i].colType() == 1
          // No corresponding anticolour in final state
       && ( FindCol(event[i].col(),i,0,event,1,true) == 0
          // No corresponding colour in initial state
         && FindCol(event[i].col(),i,0,event,2,true) == 0 )) {
     validColour = false;
     break;
   // Check anticolour of antiquarks
   } else if ( event[i].isFinal() && event[i].colType() == -1
          // No corresponding colour in final state
       && ( FindCol(event[i].acol(),i,0,event,2,true) == 0
          // No corresponding anticolour in initial state
         && FindCol(event[i].acol(),i,0,event,1,true) == 0 )) {
     validColour = false;
     break;
   // No uncontracted colour (anticolour) charge of gluons
   } else if ( event[i].isFinal() && event[i].colType() == 2
          // No corresponding anticolour in final state
       && ( FindCol(event[i].col(),i,0,event,1,true) == 0
          // No corresponding colour in initial state
         && FindCol(event[i].col(),i,0,event,2,true) == 0 )
          // No corresponding colour in final state
       && ( FindCol(event[i].acol(),i,0,event,2,true) == 0
          // No corresponding anticolour in initial state
         && FindCol(event[i].acol(),i,0,event,1,true) == 0 )) {
     validColour = false;
     break;
   }

  // Check charge sum in initial and final state
  bool validCharge = true;
  double initCharge = event[3].charge() + event[4].charge();
  double finalCharge = 0.0;
  for(int i = 0; i < event.size(); ++i)
    if (event[i].isFinal()) finalCharge += event[i].charge();
  if (abs(initCharge-finalCharge) > 1e-12) validCharge = false;

  return (validColour && validCharge);

}

//--------------------------------------------------------------------------

// Function to check whether two clusterings are identical, used
// for finding the history path in the mother -> children direction

bool History::equalClustering( Clustering clus1 , Clustering clus2 ) {
  return (  (clus1.emittor  == clus2.emittor)
         && (clus1.emitted  == clus2.emitted)
         && (clus1.recoiler == clus2.recoiler)
         && (clus1.partner  == clus2.partner)
         && (clus1.pT()    == clus2.pT()) );
}

//--------------------------------------------------------------------------

// Chose dummy scale for event construction. By default, choose
//     sHat     for 2->Boson(->2)+ n partons processes and
//     M_Boson  for 2->Boson(->)             processes

double History::choseHardScale( const Event& event ) const {

  // Get sHat
  double mHat = (event[3].p() + event[4].p()).mCalc();

  // Find number of final state particles and bosons
  int nFinal = 0;
  int nFinBos= 0;
  int nBosons= 0;
  double mBos = 0.0;
  for(int i = 0; i < event.size(); ++i)
    if ( event[i].isFinal() ) {
      nFinal++;
      // Remember final state unstable bosons
      if ( event[i].idAbs() == 23
        || event[i].idAbs() == 24 ) {
          nFinBos++;
          nBosons++;
          mBos += event[i].m();
      }
    } else if ( abs(event[i].status()) == 22
             && (  event[i].idAbs() == 23
                || event[i].idAbs() == 24 )) {
      nBosons++;
      mBos += event[i].m(); // Real mass
    }

  // Return averaged boson masses
  if ( nBosons > 0 && (nFinal + nFinBos*2) <= 3)
    return (mBos / double(nBosons));
  else return
    mHat;
}


//--------------------------------------------------------------------------

// If the state has an incoming hadron return the flavour of the
// parton entering the hard interaction. Otherwise return 0

int History::getCurrentFlav(const int side) const {
  int in = (side == 1) ? 3 : 4;
  return state[in].id();
}

//--------------------------------------------------------------------------

double History::getCurrentX(const int side) const {
  int in = (side == 1) ? 3 : 4;
  return ( 2.*state[in].e()/state[0].e() );
}

//--------------------------------------------------------------------------

double History::getCurrentZ(const int rad,
  const int rec, const int emt) const {

  int type = state[rad].isFinal() ? 1 : -1;
  double z = 0.;

  if (type == 1) {
    // Construct 2->3 variables for FSR
    Vec4   sum     = state[rad].p() + state[rec].p()
                   + state[emt].p();
    double m2Dip = sum.m2Calc();
    double x1 = 2. * (sum * state[rad].p()) / m2Dip;
    double x3 = 2. * (sum * state[emt].p()) / m2Dip;
    // Calculate z of splitting, different for FSR
    z = x1/(x1+x3);
  } else {
    // Construct momenta of dipole before/after splitting for ISR
    Vec4 qBR(state[rad].p() - state[emt].p() + state[rec].p());
    Vec4 qAR(state[rad].p() + state[rec].p());
    // Calculate z of splitting, different for ISR
    z = (qBR.m2Calc())/( qAR.m2Calc());
  }

  return z;

}

//--------------------------------------------------------------------------

// Function to compute "pythia pT separation" from Particle input

double History::pTLund(const Particle& RadAfterBranch,
              const Particle& EmtAfterBranch,
              const Particle& RecAfterBranch, int ShowerType) {

  // Save type: 1 = FSR pT definition, else ISR definition
  int Type   = ShowerType;
  // Calculate virtuality of splitting
  int sign = (Type==1) ? 1 : -1;
  Vec4 Q(RadAfterBranch.p() + sign*EmtAfterBranch.p());
  double Qsq = sign * Q.m2Calc();
  // Mass term of radiator: Consider all non-light quark or gluon radiators
  // massive partons.
  bool isMassive = ( RadAfterBranch.idAbs() >= 4
                  && RadAfterBranch.id() != 21 );
  double m2Rad   = ( mergingHooksPtr->includeMassive() && isMassive )
                 ? pow2( particleDataPtr->m0(RadAfterBranch.id()) ) : 0.;
  // Construct 2->3 variables for FSR
  Vec4   sum     = RadAfterBranch.p() + RecAfterBranch.p()
                 + EmtAfterBranch.p();
  double m2Dip = sum.m2Calc();
  double x1 = 2. * (sum * RadAfterBranch.p()) / m2Dip;
  double x3 = 2. * (sum * EmtAfterBranch.p()) / m2Dip;
  // Construct momenta of dipole before/after splitting for ISR
  Vec4 qBR(RadAfterBranch.p() - EmtAfterBranch.p() + RecAfterBranch.p());
  Vec4 qAR(RadAfterBranch.p() + RecAfterBranch.p());
  // Calculate z of splitting, different for FSR and ISR

  double z = (Type==1) ? x1/(x1+x3)
                     : (qBR.m2Calc())/( qAR.m2Calc());
  // Separation of splitting, different for FSR and ISR
  double pTpyth = (Type==1) ? z*(1.-z) : (1.-z);
  // pT^2 = separation*virtuality
  pTpyth *= (Qsq - sign*m2Rad);
  if ( pTpyth < 0. ) pTpyth = 0.;

  // Return pT
  return sqrt(pTpyth);
}

//--------------------------------------------------------------------------

// Function to return the position of the initial line before (or after)
// a single (!) splitting.

int History::posChangedIncoming(const Event& event, bool before) {

  // Check for initial state splittings.
  // Consider a splitting to exist if both mother and sister were found.
  // Find sister
  int iSister = 0;
  for (int i =0; i < event.size(); ++i)
    if (event[i].status() == 43) {
      iSister = i;
      break;
    }
  // Find mother
  int iMother = 0;
  if (iSister > 0) iMother = event[iSister].mother1();

  // Initial state splitting has been found.
  if (iSister > 0 && iMother > 0) {

    // Find flavour, mother flavour
    int flavSister  = event[iSister].id();
    int flavMother  = event[iMother].id();

    // Find splitting flavour
    int flavDaughter = 0;
    if ( abs(flavMother) < 21 && flavSister     == 21)
      flavDaughter = flavMother;
    else if ( flavMother     == 21 && flavSister     == 21)
      flavDaughter = flavMother;
    else if ( flavMother     == 21 && abs(flavSister) < 21)
      flavDaughter = -1*flavSister;
    else if ( abs(flavMother) < 21 && abs(flavSister) < 21)
      flavDaughter = 21;

    // Find initial state (!) daughter
    int iDaughter = 0;
    for (int i =0; i < event.size(); ++i)
      if ( !event[i].isFinal()
        && event[i].mother1() == iMother
        && event[i].id()      == flavDaughter )
        iDaughter = i;

    // Done for initial state splitting.
    if ( !before ) return iMother;
    else return iDaughter;

  }

  // Check for final state splittings with initial state recoiler.
  // Consider a splitting to exist if both mother and daughter were found.
  // Find new mother
  iMother = 0;
  for (int i =0; i < event.size(); ++i)
    if ( abs(event[i].status()) == 53 || abs(event[i].status()) == 54) {
      iMother = i;
      break;
    }
  // Find daughter
  int iDaughter = 0;
  if (iMother > 0) iDaughter = event[iMother].daughter1();

  // Done if final state splitting has been found.
  if (iDaughter > 0 && iMother > 0) {

    // Done for final state splitting.
    if ( !before ) return iMother;
    else return iDaughter;

  }

  // If no splitting has been found, return zero.
  return 0;

}

//--------------------------------------------------------------------------

// Function to give back the ratio of PDFs and PDF * splitting kernels needed
// to convert a splitting at scale pdfScale, chosen with running PDFs, to a
// splitting chosen with PDFs at a fixed scale mu. As needed to properly count
// emissions.

double History::pdfFactor( const Event& event, const int type,
  double pdfScale, double mu ) {

  double weight = 1.;

  // Final state splittings
  if (type >= 3) {

    // Find new mother
    int iMother = 0;
    for (int i =0; i < event.size(); ++i)
      if ( abs(event[i].status()) == 53 || abs(event[i].status()) == 54) {
        iMother = i;
        break;
      }
    int flavMother = event[iMother].id();

    // Done if no initial state recoiler was found
    if ( iMother == 0 ) return 1.;

    // Find daughter
    int iDaughter    = event[iMother].daughter1();
    int flavDaughter = event[iDaughter].id();

    // Find x values
    double xMother = 2.*event[iMother].e() / event[0].e();
    double xDaughter = 2.*event[iDaughter].e() / event[0].e();

    // Calculate PDF ratios

    int sideSplit = ( event[iMother].pz() > 0.) ? 1 : -1;
    double pdfDen1, pdfDen2, pdfNum1, pdfNum2;
    pdfDen1 = pdfDen2 = pdfNum1 = pdfNum2 = 1.;
    if ( sideSplit == 1 ) {
      // Find PDFs
      pdfDen1 = max(1e-15,beamA.xfISR(0, flavDaughter, xDaughter, pow2(mu)) );
      pdfNum1 = beamA.xfISR(0, flavDaughter, xDaughter, pow2(pdfScale) );
      pdfNum2 = beamA.xfISR(0, flavMother, xMother, pow2(mu) );
      pdfDen2 = max(1e-15,beamA.xfISR(0,flavMother, xMother, pow2(pdfScale)) );
    } else {
      // Find PDFs
      pdfDen1 = max(1e-15,beamB.xfISR(0, flavDaughter, xDaughter, pow2(mu)) );
      pdfNum1 = beamB.xfISR(0, flavDaughter, xDaughter, pow2(pdfScale) );
      pdfNum2 = beamB.xfISR(0, flavMother, xMother, pow2(mu) );
      pdfDen2 = max(1e-15,beamB.xfISR(0,flavMother, xMother, pow2(pdfScale)) );
    }

    // The magnitude of the PDF ratio in FSR is limited to one. If that was
    // the case, return one.
    if ( pdfDen2/pdfNum1 > 1. ) return 1.;

    // Calculate PDF weight to reweight emission to emission evaluated at
    // constant factorisation scale. No need to include the splitting kernel in
    // the weight, since it will drop out anyway.
    weight = (pdfNum1/pdfDen1) * (pdfNum2)/(pdfDen2);

  // Initial state splittings
  } else if (type == 2) {

    // Find sister
    int iSister = 0;
    for (int i =0; i < event.size(); ++i)
      if (event[i].status() == 43) {
        iSister = i;
        break;
      }
    int flavSister = event[iSister].id();

    // Find mother
    int iMother    = event[iSister].mother1();
    int flavMother = event[iMother].id();

    // Find splitting flavour
    int flavDaughter = 0;
    if ( abs(flavMother) < 21 && flavSister     == 21)
      flavDaughter = flavMother;
    else if ( flavMother     == 21 && flavSister     == 21)
      flavDaughter = flavMother;
    else if ( flavMother     == 21 && abs(flavSister) < 21)
      flavDaughter = -1*flavSister;
    else if ( abs(flavMother) < 21 && abs(flavSister) < 21)
      flavDaughter = 21;

    // Find x values
    double xMother = 2.*event[iMother].e() / event[0].e();

    // Find initial state (!) daughter
    int iDaughter = 0;
    for (int i =0; i < event.size(); ++i)
      if ( !event[i].isFinal()
        && event[i].mother1() == iMother
        && event[i].id()      == flavDaughter )
        iDaughter = i;
    double xDaughter = 2.*event[iDaughter].e() / event[0].e();

    // Calculate PDF weight to reweight emission to emission evaluated at
    // constant factorisation scale. No need to include the splitting kernel
    // in the weight, since it will drop out anyway.
    int sideSplit = ( event[iMother].pz() > 0.) ? 1 : -1;
    double ratio1 = getPDFratio( sideSplit, false, false, flavDaughter,
                      xDaughter, pdfScale, flavDaughter, xDaughter, mu );
    double ratio2 = getPDFratio( sideSplit, false, false, flavMother,
                      xMother, mu, flavMother, xMother, pdfScale );

    weight = ratio1*ratio2;

  // Do nothing for MPI
  } else {
    weight = 1.;
  }

  // Done
  return weight;
}

//--------------------------------------------------------------------------

// Function giving the product of splitting kernels and PDFs so that the
// resulting flavour is given by flav. This is used as a helper routine
// to dgauss

double History::integrand(int flav, double x, double scaleInt, double z) {

  // Declare constants
  double CF = 4./3.;
  double TR = 1./2.;
  double CA = 3.;

  double result = 0.;

  // Integrate NLL sudakov remainder
  if (flav==0) {

    AlphaStrong* as = mergingHooksPtr->AlphaS_ISR();
    double asNow = (*as).alphaS(z);
    result = 1./z *asNow*asNow* ( log(scaleInt/z) -3./2. );

  // Integrand for PDF ratios. Careful about factors if 1/z, since formulae
  // are expressed in terms if f(x,mu), while Pythia uses x*f(x,mu)!
  } else if (flav==21) {

    double measure1 = 1./(1. - z);
    double measure2 = 1.;

    double integrand1 =
      2.*CA
      * z * beamB.xf( 21,x/z,pow(scaleInt,2))
          / beamB.xf( 21,x,  pow(scaleInt,2))
    - 2.*CA;

    double integrand2 =
      // G -> G terms
      2.*CA  *((1. -z)/z + z*(1.-z))
      * beamB.xf( 21,x/z,pow(scaleInt,2))
      / beamB.xf( 21,x,  pow(scaleInt,2))
      // G -> Q terms
    + CF * ((1+pow(1-z,2))/z)
      *( beamB.xf(  1, x/z,pow(scaleInt,2))
       / beamB.xf( 21, x,  pow(scaleInt,2))
       + beamB.xf( -1, x/z,pow(scaleInt,2))
       / beamB.xf( 21, x,  pow(scaleInt,2))
       + beamB.xf(  2, x/z,pow(scaleInt,2))
       / beamB.xf( 21, x,  pow(scaleInt,2))
       + beamB.xf( -2, x/z,pow(scaleInt,2))
       / beamB.xf( 21, x,  pow(scaleInt,2))
       + beamB.xf(  3, x/z,pow(scaleInt,2))
       / beamB.xf( 21, x,  pow(scaleInt,2))
       + beamB.xf( -3, x/z,pow(scaleInt,2))
       / beamB.xf( 21, x,  pow(scaleInt,2))
       + beamB.xf(  4, x/z,pow(scaleInt,2))
       / beamB.xf( 21, x,  pow(scaleInt,2))
       + beamB.xf( -4, x/z,pow(scaleInt,2))
       / beamB.xf( 21, x,  pow(scaleInt,2)) );

    // Done
    result = integrand1*measure1 + integrand2*measure2;

  } else {

    double measure1 = 1./(1. -z);
    double measure2 = 1.;

    // Q -> Q terms
    double integrand1 =
      CF * (1+pow(z,2))
      * beamB.xf( flav, x/z, pow(scaleInt,2))
      / beamB.xf( flav, x,   pow(scaleInt,2))
    - 2.*CF;

    // Q -> G terms
    double integrand2 =
    + TR * (pow(z,2) + pow(1-z,2))
      * beamB.xf( 21,   x/z, pow(scaleInt,2))
      / beamB.xf( flav, x,   pow(scaleInt,2));

    // Done
    result = measure1*integrand1 + measure2*integrand2;
  }

  return result;

}

//==========================================================================

} // end namespace Pythia8
