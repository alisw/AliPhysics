// ProcessLevel.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class for process-level event generation.
// ProcessLevel: administrates the selection of "hard" process.

#ifndef Pythia8_ProcessLevel_H
#define Pythia8_ProcessLevel_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonDistributions.h"
#include "Pythia8/ProcessContainer.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/ResonanceDecays.h"
#include "Pythia8/Settings.h"
#include "Pythia8/SigmaTotal.h"
#include "Pythia8/SusyCouplings.h"
#include "Pythia8/SLHAinterface.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/UserHooks.h"

namespace Pythia8 {

//==========================================================================

// The ProcessLevel class contains the top-level routines to generate
// the characteristic "hard" process of an event.

class ProcessLevel {

public:

  // Constructor.
  ProcessLevel() : iLHACont(-1) {}

  // Destructor to delete processes in containers.
  ~ProcessLevel();

  // Initialization.
  bool init( Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    Couplings* couplingsPtrIn, SigmaTotal* sigmaTotPtrIn, bool doLHAin,
    SLHAinterface* slhaInterfacePtrIn, UserHooks* userHooksPtrIn,
    vector<SigmaProcess*>& sigmaPtrs, ostream& os = cout);

  // Store or replace Les Houches pointer.
  void setLHAPtr( LHAup* lhaUpPtrIn) {lhaUpPtr = lhaUpPtrIn;
    if (iLHACont >= 0) containerPtrs[iLHACont]->setLHAPtr(lhaUpPtr);}

  // Generate the next "hard" process.
  bool next( Event& process);

  // Special case: LHA input of resonance decay only.
  bool nextLHAdec( Event& process);

  // Accumulate and update statistics (after possible user veto).
  void accumulate();

  // Print statistics on cross sections and number of events.
  void statistics(bool reset = false, ostream& os = cout);

  // Reset statistics.
  void resetStatistics();

  // Add any junctions to the process event record list.
  void findJunctions( Event& junEvent);

  // Initialize and call resonance decays separately.
  void initDecays( Info* infoPtrIn, ParticleData* particleDataPtrIn,
    Rndm* rndmPtrIn, LHAup* lhaUpPtrIn) { infoPtr = infoPtrIn;
    resonanceDecays.init( infoPtrIn, particleDataPtrIn, rndmPtrIn);
    containerLHAdec.setLHAPtr(lhaUpPtrIn, particleDataPtrIn); }
  bool nextDecays( Event& process) { return resonanceDecays.next( process);}

private:

  // Constants: could only be changed in the code itself.
  static const int MAXLOOP;

  // Generic info for process generation.
  bool   doSecondHard, doSameCuts, allHardSame, noneHardSame,
         someHardSame, cutsAgree, cutsOverlap, doResDecays;
  int    nImpact, startColTag;
  double mHatMin1, mHatMax1, pTHatMin1, pTHatMax1, mHatMin2, mHatMax2,
         pTHatMin2, pTHatMax2, sigmaND, sumImpactFac, sum2ImpactFac;

  // Vector of containers of internally-generated processes.
  vector<ProcessContainer*> containerPtrs;
  int    iContainer, iLHACont;
  double sigmaMaxSum;

  // Ditto for optional choice of a second hard process.
  vector<ProcessContainer*> container2Ptrs;
  int    i2Container;
  double sigma2MaxSum;

  // Single half-dummy container for LHA input of resonance decay only.
  ProcessContainer containerLHAdec;

  // Pointer to various information on the generation.
  Info*           infoPtr;

  // Pointer to the particle data table.
  ParticleData*   particleDataPtr;

  // Pointer to the random number generator.
  Rndm*           rndmPtr;

  // Pointers to the two incoming beams.
  BeamParticle*   beamAPtr;
  BeamParticle*   beamBPtr;

  // Pointer to Standard Model couplings, including alphaS and alphaEM.
  Couplings*      couplingsPtr;

  // Pointer to SigmaTotal object needed to handle soft QCD processes.
  SigmaTotal*     sigmaTotPtr;

  // Pointer to SusyLesHouches object for interface to SUSY spectra.
  SLHAinterface*  slhaInterfacePtr;

  // Pointer to userHooks object for user interaction with program.
  UserHooks*      userHooksPtr;

  // Pointer to LHAup for generating external events.
  LHAup*          lhaUpPtr;

  // ResonanceDecay object does sequential resonance decays.
  ResonanceDecays resonanceDecays;

  // Generate the next event with one interaction.
  bool nextOne( Event& process);

  // Generate the next event with two hard interactions.
  bool nextTwo( Event& process);

  // Append the second to the first process list.
  void combineProcessRecords( Event& process, Event& process2);

  // Check that colours match up.
  bool checkColours( Event& process);

  // Print statistics when two hard processes allowed.
  void statistics2(bool reset, ostream& os = cout);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ProcessLevel_H
