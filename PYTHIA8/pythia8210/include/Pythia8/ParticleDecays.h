// ParticleDecays.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the classes to perform a particle decay.
// DecayHandler: base class for external handling of decays.
// ParticleDecays: decay a particle.

#ifndef Pythia8_ParticleDecays_H
#define Pythia8_ParticleDecays_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/TimeShower.h"
#include "Pythia8/TauDecays.h"

namespace Pythia8 {

//==========================================================================

// DecayHandler is base class for the external handling of decays.
// There is only one pure virtual method, that should do the decay.

class DecayHandler {

public:

  // Destructor.
  virtual ~DecayHandler() {}

  // A pure virtual method, wherein the derived class method does a decay.
  virtual bool decay(vector<int>& idProd, vector<double>& mProd,
    vector<Vec4>& pProd, int iDec, const Event& event) = 0;

};

//==========================================================================

// The ParticleDecays class contains the routines to decay a particle.

class ParticleDecays {

public:

  // Constructor.
  ParticleDecays() {}

  // Initialize: store pointers and find settings
  void init(Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    Couplings* couplingsPtrIn, TimeShower* timesDecPtrIn,
    StringFlav* flavSelPtrIn, DecayHandler* decayHandlePtrIn,
    vector<int> handledParticles);

  // Perform a decay of a single particle.
  bool decay(int iDec, Event& event);

  // Did decay result in new partons to hadronize?
  bool moreToDo() const {return hasPartons && keepPartons;}

private:

  // Constants: could only be changed in the code itself.
  static const int    NTRYDECAY, NTRYPICK, NTRYMEWT, NTRYDALITZ;
  static const double MSAFEDALITZ, WTCORRECTION[11];

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Pointers to Standard Model couplings.
  Couplings*    couplingsPtr;

  // Pointers to timelike showers, for decays to partons (e.g. Upsilon).
  TimeShower*   timesDecPtr;

  // Pointer to class for flavour generation; needed when to pick hadrons.
  StringFlav*   flavSelPtr;

  // Pointer to a handler of external decays.
  DecayHandler* decayHandlePtr;

  // Initialization data, read from Settings.
  bool   limitTau0, limitTau, limitRadius, limitCylinder, limitDecay,
         mixB, doFSRinDecays, doGammaRad;
  int    tauMode;
  double mSafety, tau0Max, tauMax, rMax, xyMax, zMax, xBdMix, xBsMix,
         sigmaSoft, multIncrease, multIncreaseWeak, multRefMass, multGoffset,
         colRearrange, stopMass, sRhoDal, wRhoDal;

  // Multiplicity. Decay products positions and masses.
  bool   hasPartons, keepPartons;
  int    idDec, meMode, mult;
  double scale;
  vector<int>    iProd, idProd, cols, acols, idPartons;
  vector<double> mProd, mInv, rndmOrd;
  vector<Vec4>   pInv, pProd;
  vector<FlavContainer> flavEnds;

  // Pointer to particle data for currently decaying particle
  ParticleDataEntry* decDataPtr;

  // Tau particle decayer.
  TauDecays tauDecayer;

  // Check whether a decay is allowed, given the upcoming decay vertex.
  bool checkVertex(Particle& decayer);

  // Check for oscillations B0 <-> B0bar or B_s0 <-> B_s0bar.
  bool oscillateB(Particle& decayer);

  // Do a one-body decay.
  bool oneBody(Event& event);

  // Do a two-body decay;
  bool twoBody(Event& event);

  // Do a three-body decay;
  bool threeBody(Event& event);

  // Do a multibody decay using the M-generator algorithm.
  bool mGenerator(Event& event);

  // Select mass of lepton pair in a Dalitz decay.
  bool dalitzMass();

  // Do kinematics of gamma* -> l- l+ in Dalitz decay.
  bool dalitzKinematics(Event& event);

  // Translate a partonic content into a set of actual hadrons.
  bool pickHadrons();

  // Set colour flow and scale in a decay explicitly to partons.
  bool setColours(Event& event);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ParticleDecays_H
