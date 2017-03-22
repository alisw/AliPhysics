// WeakShowerMEs.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the class containg matrix elements needed for
// W/Z emission corrections in both the initial and final state shower.
// WeakShowerMEs: contains the matrix elements.

#ifndef Pythia8_WeakShowerMEs_H
#define Pythia8_WeakShowerMEs_H

#include "Pythia8/Basics.h"

namespace Pythia8 {

//==========================================================================

// The WeakShowerMEs provides ME's needed for W/Z emission in ISR or FSR.

class WeakShowerMEs {

public:

  // Constructor.
  WeakShowerMEs() {}

  // Calculate the 2 to 2 ME uG -> uG, up to a known overall factor.
  double getTchanneluGuGME(double sHat,double tHat,double uHat);

  // Calculate the 2 to 2 ME ud -> ud, up to a known overall factor.
  double getTchannelududME(double sHat,double tHat,double uHat);

  // Calculate the 2 to 2 ME uu -> uu, up to a known overall factor.
  double getTchanneluuuuME(double sHat,double tHat,double uHat);

  // Calculate the 2 to 3 ME uG -> uGZ, up to a known overall factor.
  double getTchanneluGuGZME(Vec4 p1,Vec4 p2,Vec4 p3,Vec4 p4,Vec4 p5);

  // Calculate the 2 to 3 ME ud -> udZ, up to a known overall factor,
  // and with the coupling between Z and d set to zero.
  double getTchannelududZME(Vec4 p1,Vec4 p2,Vec4 p3,Vec4 p4,Vec4 p5);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_WeakShowerMEs_H
