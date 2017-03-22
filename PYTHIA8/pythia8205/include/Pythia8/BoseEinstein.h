// Bose-Einstein.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the classes to handle Bose-Einstein effects.
// BoseEinsteinHadron: simple working container for particle momenta.
// BoseEinstein: main class to perform the task.

#ifndef Pythia8_BoseEinstein_H
#define Pythia8_BoseEinstein_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// The BoseEinsteinHadron class is a simple container for studied hadrons.

class BoseEinsteinHadron {

public:

  // Constructors.
  BoseEinsteinHadron() : id(0), iPos(0), p(0.), pShift(0.), pComp(0.),
    m2(0.) {}
  BoseEinsteinHadron(int idIn,  int iPosIn, Vec4 pIn, double mIn) :
    id(idIn), iPos(iPosIn), p(pIn), pShift(0.), pComp(0.) {m2 = mIn*mIn;}

  // Information on hadron - all public.
  int    id, iPos;
  Vec4   p, pShift, pComp;
  double m2;

};

//==========================================================================

// The BoseEinstein class shifts the momenta of identical particles relative
// to each other, to simulate Bose-Einstein effects to some approximation.

class BoseEinstein {

public:

  // Constructor.
  BoseEinstein() {}

  // Find settings. Precalculate table used to find momentum shifts.
  bool init(Info* infoPtrIn, Settings& settings, ParticleData& particleData);

  // Perform Bose-Einstein corrections on an event.
  bool shiftEvent( Event& event);

private:

  // Constants: could only be changed in the code itself.
  static const int    IDHADRON[9], ITABLE[9], NCOMPSTEP;
  static const double STEPSIZE, Q2MIN, COMPRELERR, COMPFACMAX;

  // Initialization data, read from Settings.
  bool   doPion, doKaon, doEta;
  double lambda, QRef;

  // Pointer to various information on the generation.
  Info* infoPtr;

  // Table of momentum shifts for different hadron species.
  int    nStep[4], nStep3[4], nStored[10];
  double QRef2, QRef3, R2Ref, R2Ref2, R2Ref3, mHadron[9],
         mPair[4], m2Pair[4], deltaQ[4], deltaQ3[4], maxQ[4], maxQ3[4];
  double shift[4][200], shift3[4][200];

  // Vector of hadrons to study.
  vector<BoseEinsteinHadron> hadronBE;

  // Calculate shift and (unnormalized) compensation for pair.
  void shiftPair(int i1, int i2, int iHad);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_BoseEinstein_H

