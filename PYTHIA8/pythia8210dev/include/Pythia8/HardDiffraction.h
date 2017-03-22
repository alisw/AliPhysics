// HardDiffraction.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Christine Rasmussen.

// Header file for the HardDiffraction class.

#ifndef Pythia8_HardDiffraction_H
#define Pythia8_HardDiffraction_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/BeamRemnants.h"
#include "Pythia8/Info.h"
#include "Pythia8/MultipartonInteractions.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/SpaceShower.h"
#include "Pythia8/TimeShower.h"

namespace Pythia8 {

//==========================================================================

// HardDiffraction class.
// This class handles hard diffraction, together with PartonLevel.

class HardDiffraction {

public:

  // Constructor and destructor.
  HardDiffraction() {};
  ~HardDiffraction() {}

  // Initialise constants
  void init(Info* infoPtrIn, Settings& settingsIn, Rndm* rndmPtrIn,
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    BeamParticle* beamPomAPtrIn, BeamParticle* beamPomBPtrIn);

  // Main routine to check if event is from diffractive PDF.
  bool isDiffractive(int iBeamIn = 1, int partonIn = 0,
    double xIn = 0., double Q2In = 0., double xfIncIn = 0.);

  // Get diffractive values.
  double getXPomeronA()     {return xPomA;}
  double getXPomeronB()     {return xPomB;}
  double getTPomeronA()     {return tPomA;}
  double getTPomeronB()     {return tPomB;}
  double getThetaPomeronA() {return thetaPomA;}
  double getThetaPomeronB() {return thetaPomB;}

private:

  // Constants: could only be changed in the code itself.
  static const double TINYPDF;
  static const double POMERONMASS;
  static const double PROTONMASS;

  // Initialization and event data.
  int    pomSet, pomFlux, iBeam, idA, idB;
  double normPom, a1, a2, a3, A1, A2, A3, a0, ap, b0,
         mA, mB, s, s1, s2, s3, s4,
         xPomA, xPomB, tPomA, tPomB, thetaPomA, thetaPomB;

  // Pointer to various information on the generation.
  Info*           infoPtr;

  // Pointer to the settings database.
  Settings        settings;

  // Pointer to the random number generator.
  Rndm*           rndmPtr;

  // Pointers to incoming beams.
  BeamParticle*   beamAPtr;
  BeamParticle*   beamBPtr;
  BeamParticle*   beamPomAPtr;
  BeamParticle*   beamPomBPtr;

  // Pointer to temporary Pomeron PDF.
  BeamParticle*   tmpPDFPtr;

  // Return Pomeron flux inside proton, integrated over t.
  double xfPom(double xIn = 0.);

  // Pick a t value for a given x.
  double pickTNow(double xIn = 0.);

  // Return Pomeron flux inside proton, differential in t.
  double xfPomWithT(double xIn = 0., double tIn = 0.);

  // Make t range available as a pair.
  pair<double, double> tRange(double xIn = 0.);

  // Calculate scattering angle from  given x and t.
  double getThetaNow(double xIn = 0., double tIn = 0.);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HardDiffraction_H
