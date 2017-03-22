// BeamRemnants.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for beam-remnants handling.
// BeamRemnants: matches the remnants between the two beams.

#ifndef Pythia8_BeamRemnants_H
#define Pythia8_BeamRemnants_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "FragmentationFlavZpT.h"
#include "Info.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "PartonSystems.h"
#include "PythiaStdlib.h"
#include "Settings.h"

namespace Pythia8 {

//==========================================================================

// This class matches the kinematics of the hard-scattering subsystems
// (with primordial kT added) to that of the two beam remnants.  

class BeamRemnants {

public:

  // Constructor.
  BeamRemnants() : iDS(0) { }  

  // Initialization.
  bool init( Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn, 
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn, 
    PartonSystems* partonSystemsPtrIn);

  // New beams possible for handling of hard diffraction.
  void reassignBeamPtrs( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    int iDSin) {beamAPtr = beamAPtrIn; beamBPtr = beamBPtrIn; iDS = iDSin;}

  // Select the flavours/kinematics/colours of the two beam remnants. 
  bool add( Event& event);

private: 

  // Constants: could only be changed in the code itself.
  static const bool   ALLOWCOLOURTWICE, CORRECTMISMATCH;
  static const int    NTRYCOLMATCH, NTRYKINMATCH;

  // Initialization data, read from Settings.
  bool   doPrimordialKT, allowRescatter, doRescatterRestoreY, doReconnect;
  double primordialKTsoft, primordialKThard, primordialKTremnant,
         halfScaleForKT, halfMassForKT, reconnectRange, 
         pT0Ref, ecmRef, ecmPow;

  // Information set for events.
  int    nSys, oldSize, iDS;
  double eCM, sCM, pT0, pT20Rec;

  // Colour collapses (when one colour is mapped onto another).
  vector<int> colFrom, colTo;

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointers to the two incoming beams.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;

  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // Do the kinematics of the collision subsystems and two beam remnants. 
  bool setKinematics( Event& event);

  // Allow colour reconnections.
  bool reconnectColours( Event&  event);

  // Check that colours are consistent.
  bool checkColours( Event& event);

};
 
//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_BeamRemnants_H
