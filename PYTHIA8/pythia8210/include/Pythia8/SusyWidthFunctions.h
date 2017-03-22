// SusyResonanceWidths.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand
// Main author of this file: N. Desai
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for resonance properties: dynamical widths etc.
// WidthFunctions: base class for SUSY 3-body decay width functions.

#ifndef Pythia8_SusyWidthFunctions_H
#define Pythia8_SusyWidthFunctions_H

#include "Pythia8/ParticleData.h"
#include "Pythia8/SusyCouplings.h"

namespace Pythia8 {

//==========================================================================

class WidthFunction {

public:

  // Constructor and destructor.
  WidthFunction() { };
  virtual ~WidthFunction() { };

  // Public methods.
  void setPointers( ParticleData* particleDataPtrIn, CoupSUSY* coupSUSYPtrIn,
    Info* infoPtrIn);
  virtual double getWidth( int, int) { return 0.0; };

protected:

  virtual double function(double xin);

  // Gaussian integrator.
  double integrateGauss(double xmin, double xmax, double tol);

  ParticleData* particleDataPtr;
  CoupSUSY* coupSUSYPtr;
  Info* infoPtr;
  int idRes, idInt, id1, id2, id3, id4;
  double mRes, mInt, gammaInt, m1, m2 , m3, m4;

};

//==========================================================================

// Class StauWidths.

class StauWidths : public WidthFunction {

public:

  // Destructor.
  ~StauWidths() { };

  // Public method.
  double getWidth(int idResIn, int idIn);

protected:

  int fnSwitch; // Switch between multiple functions
  void setChannel(int idResIn, int idIn);
  double function(double xin);

  double delm, f0, gf, cons, wparam;
  complex gL, gR;

};

//==========================================================================

}  // end namespace Pythia8

#endif // end Pythia8_SusyResonanceWidths_H
