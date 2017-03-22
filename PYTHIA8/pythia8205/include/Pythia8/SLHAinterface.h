// SLHAinterface.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for SUSY Les Houches Accord Interface.
// Handles the communication between PYTHIA and the SusyLesHouches classes.

#ifndef Pythia8_SLHAinterface_H
#define Pythia8_SLHAinterface_H

#include "Pythia8/Basics.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SusyCouplings.h"
#include "Pythia8/SusyLesHouches.h"

namespace Pythia8 {

//==========================================================================

// The SLHAinterface class handles communication between Pythia and
// SusyLesHouches.

class SLHAinterface {

public:

  // Constructor.
  SLHAinterface() {} ;

  // Set pointers
  void setPtr( Info* infoPtrIn ) {infoPtr     = infoPtrIn;}

  // Initialize and switch to SUSY couplings if reading SLHA spectrum
  void init( Settings& settings, Rndm* rndmPtr, Couplings* couplingsPtrIn,
    ParticleData* particleDataPtr, bool& useSHLAcouplings );

  // Initialize SUSY Les Houches Accord data.
  bool initSLHA(Settings& settings, ParticleData* particleDataPtr);

  // Initialize SLHA blocks SMINPUTS and MASS from PYTHIA SM parameter values.
  // E.g., to make sure that there are no important unfilled entries
  void pythia2slha(ParticleData* particleDataPtr);

  // SusyLesHouches - SLHA object for interface to SUSY spectra.
  SusyLesHouches slha;

  // SLHA derived couplings class and pointer to Couplings object
  CoupSUSY       coupSUSY;
  Couplings*     couplingsPtr;

  // Pointers to PYTHIA objects
  Info*          infoPtr;
  Settings*      settingsPtr;

  // Internal data members
  int            meMode;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SLHAinterface_H




