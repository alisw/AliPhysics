// Merging.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel.
// Merging: Wpapper class to interface matrix element merging schemes with
//          Pythia

#ifndef Pythia8_Merging_H
#define Pythia8_Merging_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/PartonLevel.h"
#include "Pythia8/MergingHooks.h"
#include "Pythia8/History.h"

namespace Pythia8 {

//==========================================================================

// Merging is a wrapper class for the interface of matrix element merging and
// Pythia8.

class Merging {

public:

  // Constructor.
  Merging() { settingsPtr = 0; infoPtr = 0; particleDataPtr = 0;
    rndmPtr = 0; beamAPtr = 0; beamBPtr = 0; trialPartonLevelPtr = 0;
    mergingHooksPtr = 0; }

  // Make Pythia class friend
  friend class Pythia;

  // Destructor.
  ~Merging(){}

protected:

  //----------------------------------------------------------------------//
  // The members
  //----------------------------------------------------------------------//


  // Settings: databases of flags/modes/parms/words to control run.
  Settings*      settingsPtr;

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

  // Pointer to random number generator.
  Rndm* rndmPtr;

  // Pointer to trial PartonLevel object
  PartonLevel* trialPartonLevelPtr;

  // Pointer to trial MergingHooks object
  MergingHooks* mergingHooksPtr;

  // Pointers to beam particles.
  BeamParticle* beamAPtr;
  BeamParticle* beamBPtr;

  // Minimal value found for the merging scale in events.
  double tmsNowMin;
  static const double TMSMISMATCH;

  // Initialisation function for internal use inside Pythia source code
  void init( Settings* settingsPtrIn, Info* infoPtrIn,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    MergingHooks* mergingHooksPtrIn, PartonLevel* trialPartonLevelPtrIn );

  // Function to print statistics.
  void statistics(ostream& os = cout);

  //----------------------------------------------------------------------//
  // Functions that implement matrix element merging.
  //----------------------------------------------------------------------//

  // Function to steer different merging prescriptions.
  int mergeProcess( Event& process);

  // Function to perform CKKW-L merging on the event.
  int mergeProcessCKKWL( Event& process);

  // Function to perform UMEPS merging on the event.
  int mergeProcessUMEPS( Event& process);

  // Function to perform NL3 NLO merging on the event.
  int mergeProcessNL3( Event& process);

  // Function to perform UNLOPS merging on the event.
  int mergeProcessUNLOPS( Event& process);

  // Function to apply the merging scale cut on an input event.
  bool cutOnProcess( Event& process);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Merging_H
