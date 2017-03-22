// PartonLevel.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class for parton-level event generation
// PartonLevel: administrates showers, multiparton interactions and remnants.

#ifndef Pythia8_PartonLevel_H
#define Pythia8_PartonLevel_H

#include "Basics.h"
#include "BeamParticle.h"
#include "BeamRemnants.h"
#include "Event.h"
#include "Info.h"
#include "MultipartonInteractions.h"
#include "ParticleData.h"
#include "PartonSystems.h"
#include "PythiaStdlib.h"
#include "RHadrons.h"
#include "Settings.h"
#include "SigmaTotal.h"
#include "SpaceShower.h"
#include "StandardModel.h"
#include "TimeShower.h"
#include "UserHooks.h"
#include "MergingHooks.h"

namespace Pythia8 {
 
//==========================================================================

// The PartonLevel class contains the top-level routines to generate
// the partonic activity of an event.

class PartonLevel {

public:

  // Constructor. 
  PartonLevel() : userHooksPtr(0) {} 
 
  // Initialization of all classes at the parton level.
  bool init( Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn, 
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn, 
    BeamParticle* beamPomAPtrIn, BeamParticle* beamPomBPtrIn, 
    Couplings* couplingsPtrIn, PartonSystems* partonSystemsPtrIn, 
    SigmaTotal* sigmaTotPtr, TimeShower* timesDecPtrIn, 
    TimeShower* timesPtrIn, SpaceShower* spacePtrIn, 
    RHadrons* rHadronsPtrIn, UserHooks* userHooksPtrIn,
    MergingHooks* mergingHooksPtr, bool useAsTrial);
 
  // Generate the next parton-level process.
  bool next( Event& process, Event& event); 

  // Perform showers in resonance decay chains. (For special cases.)
  void setupShowerSys( Event& process, Event& event);
  bool resonanceShowers( Event& process, Event& event, bool skipForR); 

  // Tell whether failure was due to vetoing.
  bool hasVetoed() const {return doVeto;}

  // Accumulate, print and reset statistics.
  void accumulate() {if (isResolved && !isDiff) multiPtr->accumulate();}
  void statistics(bool reset = false) {
    if (doMPI) multiMB.statistics(reset);}
    // For now no separate statistics for diffraction??
    //if (doMPISDA && doDiffraction) multiSDA.statistics(reset); 
    //if (doMPISDB && doDiffraction) multiSDB.statistics(reset);}
  void resetStatistics() { if (doMPI) multiMB.resetStatistics(); }

  // Reset PartonLevel object for trial shower usage.
  void resetTrial();
  // Provide the pT scale of the last branching in the shower.
  double pTLastInShower(){ return pTLastBranch; }
  // Provide the type of the last branching in the shower.
  int typeLastInShower(){ return typeLastBranch; }

private: 

  // Constants: could only be changed in the code itself.
  static const int NTRY;

  // Initialization data, mainly read from Settings.
  bool   doMinBias, doDiffraction, doMPI, doMPIMB, doMPISDA, doMPISDB, 
         doMPICD, doMPIinit, doISR, doFSRduringProcess, doFSRafterProcess,  
         doFSRinResonances, doRemnants, doSecondHard, hasLeptonBeams, 
         hasPointLeptons, canVetoPT, canVetoStep, canVetoMPIStep, 
         canVetoEarly, canSetScale, allowRH;
  double mMinDiff, mWidthDiff, pMaxDiff;

  // Event generation strategy. Number of steps. Maximum pT scales.
  bool   doVeto;
  int    nMPI, nISR, nFSRinProc, nFSRinRes, nISRhard, nFSRhard, 
         typeLatest, nVetoStep, typeVetoStep, nVetoMPIStep, iSysNow;
  double pTsaveMPI, pTsaveISR, pTsaveFSR, pTvetoPT;

  // Current event properties.
  bool   isMinBias, isDiffA, isDiffB, isDiffC, isDiff, isSingleDiff, 
         isDoubleDiff, isCentralDiff, isResolved, isResolvedA, 
         isResolvedB, isResolvedC;
  int    sizeProcess, sizeEvent, nHardDone, nHardDoneRHad, iDS;
  double eCMsave; 
  vector<bool> inRHadDecay;
  vector<int>  iPosBefShow;

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointers to the two incoming beams.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;

  // Spare copies of normal pointers. Pointers to Pomeron beam-inside-beam.
  BeamParticle*  beamHadAPtr;  
  BeamParticle*  beamHadBPtr;  
  BeamParticle*  beamPomAPtr;
  BeamParticle*  beamPomBPtr;

  // Pointers to Standard Model couplings.
  Couplings*     couplingsPtr;
  
  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // Pointer to userHooks object for user interaction with program.
  UserHooks*     userHooksPtr;

  // Pointers to timelike showers for resonance decays and the rest.
  TimeShower*    timesDecPtr;
  TimeShower*    timesPtr;

  // Pointer to spacelike showers.
  SpaceShower*   spacePtr;

  // The generator classes for multiparton interactions.
  MultipartonInteractions  multiMB;
  MultipartonInteractions  multiSDA;
  MultipartonInteractions  multiSDB;
  MultipartonInteractions  multiCD;
  MultipartonInteractions* multiPtr;

  // The generator class to construct beam-remnant kinematics. 
  BeamRemnants remnants;
  // Separate instance for central diffraction.
  BeamRemnants remnantsCD;
  
  // The RHadrons class is used to fragment off and decay R-hadrons.
  RHadrons*    rHadronsPtr;

  // Resolved diffraction: find how many systems should have it.
  int decideResolvedDiff( Event& process);

  // Set up an unresolved process, i.e. elastic or diffractive.
  bool setupUnresolvedSys( Event& process, Event& event);

  // Set up the hard process, excluding subsequent resonance decays.
  void setupHardSys( Event& process, Event& event);

  // Resolved diffraction: pick whether to have it and set up for it.
  void setupResolvedDiff( Event& process);

  // Resolved diffraction: restore normal behaviour.
  void leaveResolvedDiff( int iHardLoop, Event& process, Event& event);

  // Pointer to MergingHooks object for user interaction with the merging.
  MergingHooks* mergingHooksPtr;
  // Parameters to specify trial shower usage
  bool doTrial;
  int nTrialEmissions;
  // Parameters to store to veto trial showers
  double pTLastBranch;
  int typeLastBranch;
  // Parameters to specify merging usage
  bool canRemoveEvent, canRemoveEmission;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_PartonLevel_H
