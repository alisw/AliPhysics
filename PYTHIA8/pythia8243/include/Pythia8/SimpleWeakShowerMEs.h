// SimpleWeakShowerMEs.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the class containg matrix elements needed for
// W/Z emission corrections in both the initial and final state shower.
// SimpleWeakShowerMEs: contains the matrix elements.

#ifndef Pythia8_SimpleWeakShowerMEs_H
#define Pythia8_SimpleWeakShowerMEs_H

#include "Pythia8/Basics.h"

namespace Pythia8 {

//==========================================================================

// The SimpleWeakShowerMEs provides ME's needed for W/Z emission in ISR or FSR.
// The 2 -> 2 MEs contain the correct kinematics, but for some of
// the 2 -> 3 MEs some couplings have been switched off.
// This class is used for ME corrections in the weak shower and for merging
// whenever weak reclusterings are allowed.
// Also be aware that no phase-space factors are included.

class SimpleWeakShowerMEs {

public:

  // Constructor.
  SimpleWeakShowerMEs() {}

  // Calculate the 2 to 2 ME qg -> qg, up to a known overall factor.
  double getMEqg2qg(double sH, double tH, double uH);

  // Calculate the 2 to 2 ME qq -> qq, up to a known overall factor.
  double getMEqq2qq(double sH, double tH, double uH, bool sameID);

  // Calculate the 2 to 2 ME gg -> gg, up to a known overall factor.
  double getMEgg2gg(double sH, double tH, double uH);

  // Calculate the 2 to 2 ME gg -> qqbar, up to a known overall factor.
  double getMEgg2qqbar(double sH, double tH, double uH);

  // Calculate the 2 to 2 ME qqbar -> gg, up to a known overall factor.
  double getMEqqbar2gg(double sH, double tH, double uH);

  // Calculate the 2 to 2 ME qqbar -> qqbar, up to a known overall factor.
  double getMEqqbar2qqbar(double sH, double tH, double uH, bool sameID);

  // Calculate the 2 to 3 ME uG -> uGZ, up to a known overall factor.
  double getMEqg2qgZ(Vec4 p1,Vec4 p2,Vec4 p3,Vec4 p4,Vec4 p5);

  // Calculate the 2 to 3 ME ud -> udZ, up to a known overall factor,
  // and with the coupling between Z and d set to zero.
  double getMEqq2qqZ(Vec4 p1,Vec4 p2,Vec4 p3,Vec4 p4,Vec4 p5);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SimpleWeakShowerMEs_H
