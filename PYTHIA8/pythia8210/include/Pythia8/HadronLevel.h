// HadronLevel.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class for hadron-level generation.
// HadronLevel: handles administration of fragmentation and decay.

#ifndef Pythia8_HadronLevel_H
#define Pythia8_HadronLevel_H

#include "Pythia8/Basics.h"
#include "Pythia8/BoseEinstein.h"
#include "Pythia8/ColourTracing.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/FragmentationSystems.h"
#include "Pythia8/HadronScatter.h"
#include "Pythia8/HiddenValleyFragmentation.h"
#include "Pythia8/Info.h"
#include "Pythia8/JunctionSplitting.h"
#include "Pythia8/MiniStringFragmentation.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/ParticleDecays.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/RHadrons.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StringFragmentation.h"
#include "Pythia8/TimeShower.h"

namespace Pythia8 {

//==========================================================================

// The HadronLevel class contains the top-level routines to generate
// the transition from the partonic to the hadronic stage of an event.

class HadronLevel {

public:

  // Constructor.
  HadronLevel() {}

  // Initialize HadronLevel classes as required.
  bool init(Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    Couplings* couplingsPtrIn, TimeShower* timesDecPtr,
    RHadrons* rHadronsPtrIn, DecayHandler* decayHandlePtr,
    vector<int> handledParticles);

  // Get pointer to StringFlav instance (needed by BeamParticle).
  StringFlav* getStringFlavPtr() {return &flavSel;}

  // Generate the next event.
  bool next(Event& event);

  // Special routine to allow more decays if on/off switches changed.
  bool moreDecays(Event& event);

private:

  // Initialization data, read from Settings.
  bool   doHadronize, doDecay, doBoseEinstein, allowRH;
  double mStringMin, eNormJunction, widthSepBE;

  // Settings for hadron scattering --rjc
  bool   doHadronScatter, hsAfterDecay;

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Pointers to Standard Model couplings.
  Couplings*    couplingsPtr;

  // Configuration of colour-singlet systems.
  ColConfig     colConfig;

  // Colour and mass information.
  vector<int>    iParton, iJunLegA, iJunLegB, iJunLegC,
                 iAntiLegA, iAntiLegB, iAntiLegC, iGluLeg;
  vector<double> m2Pair;

  // The generator class for normal string fragmentation.
  StringFragmentation stringFrag;

  // The generator class for special low-mass string fragmentation.
  MiniStringFragmentation ministringFrag;

  // The generator class for normal decays.
  ParticleDecays decays;

  // The generator class for hadron scattering --rjc
  HadronScatter hadronScatter;

  // The generator class for Bose-Einstein effects.
  BoseEinstein boseEinstein;

  // Classes for flavour, pT and z generation.
  StringFlav flavSel;
  StringPT   pTSel;
  StringZ    zSel;

  // Class for colour tracing.
  ColourTracing colTrace;

  // Junction splitting class.
  JunctionSplitting junctionSplitting;

  // The RHadrons class is used to fragment off and decay R-hadrons.
  RHadrons*  rHadronsPtr;

  // Special class for Hidden-Valley hadronization. Not always used.
  HiddenValleyFragmentation hiddenvalleyFrag;
  bool useHiddenValley;

  // Special case: colour-octet onium decays, to be done initially.
  bool decayOctetOnia(Event& event);

  // Trace colour flow in the event to form colour singlet subsystems.
  bool findSinglets(Event& event);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HadronLevel_H
