// JunctionSplitting.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the class JunctionSplitting.
// JunctionSplitting takes an event and seperate junction chains from
// each other, such that no junctions are colour connected to each other.

#ifndef Pythia8_JunctionSplitting_H
#define Pythia8_JunctionSplitting_H

#include "Pythia8/Basics.h"
#include "Pythia8/ColourTracing.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StringLength.h"

namespace Pythia8 {

//==========================================================================

// JunctionSplitting takes an event and seperate junction chains from
// each other, such that no junctions are colour connected to each other.

class JunctionSplitting {

public:

  // Constructor
  JunctionSplitting() : eNormJunction(), allowDoubleJunRem(), infoPtr(),
    rndmPtr() {}

  // Initialization.
  void init(Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn,
    ParticleData* particleDataPtrIn);

  // Test whether an event has a physical colour configuration.
  // It also splits junction pairs into pieces that PYTHIA can hadronize.
  bool checkColours(Event& event);

private:

  // Constants: could only be changed in the code itself.
  static const int    NTRYJNREST;
  static const double JJSTRINGM2MAX, JJSTRINGM2FRAC, CONVJNREST, MTHAD,
                      MINANGLE;

  double eNormJunction;
  bool allowDoubleJunRem;
  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Classes for flavour, pT and z generation.
  StringFlav flavSel;
  StringPT   pTSel;
  StringZ    zSel;

  // The generator class for normal string fragmentation.
  StringFragmentation stringFrag;

  // Colour tracing class used to find the colour chains.
  ColourTracing colTrace;

  // String Length class used to calculate the string length.
  StringLength stringLength;

  // Split connected junction chains into separated, mainly by splitting
  // gluons into q-qbar pairs.
  bool splitJunGluons(Event& event, vector<vector< int > >& iPartonJun,
    vector<vector< int > >& iPartonAntiJun);

  // Split multiple (> 2) directly connected junctions.
  bool splitJunChains(Event& event);

  // Split junction pairs.
  bool splitJunPairs(Event& event, vector<vector< int > >& iPartonJun,
    vector<vector< int > >& iPartonAntiJun);

  // Get the list of partons connected to the junctions.
  bool getPartonLists(Event& event, vector<vector< int > >& iPartonJun,
    vector<vector<int > >& iPartonAntiJun);

  // Change the anticolour of the particle that has acol to be col.
  bool setAcol(Event& event, int col, int acol);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_JunctionSplitting_H
