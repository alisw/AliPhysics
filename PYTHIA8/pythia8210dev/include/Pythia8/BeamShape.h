// BeamShape.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header for classes to set beam momentum and interaction vertex spread.

#ifndef Pythia8_BeamShape_H
#define Pythia8_BeamShape_H

#include "Pythia8/Basics.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// Base class to set beam momentum and interaction spot spread.

class BeamShape {

public:

  // Constructor.
  BeamShape() {}

  // Destructor.
  virtual ~BeamShape() {}

  // Initialize beam parameters.
  virtual void init( Settings& settings, Rndm* rndmPtrIn);

  // Set the two beam momentum deviations and the beam vertex.
  virtual void pick();

  // Methods to read out the choice made with the above method.
  Vec4 deltaPA() const {return Vec4( deltaPxA, deltaPyA, deltaPzA, 0);}
  Vec4 deltaPB() const {return Vec4( deltaPxB, deltaPyB, deltaPzB, 0);}
  Vec4 vertex()  const {return Vec4( vertexX, vertexY, vertexZ, vertexT);}

protected:

  // Values to be set.
  double deltaPxA, deltaPyA, deltaPzA, deltaPxB, deltaPyB, deltaPzB,
         vertexX, vertexY, vertexZ, vertexT;

  // Parameters of Gaussian parametrizations.
  bool   allowMomentumSpread, allowVertexSpread;
  double sigmaPxA, sigmaPyA, sigmaPzA, maxDevA, sigmaPxB, sigmaPyB,
         sigmaPzB, maxDevB, sigmaVertexX, sigmaVertexY, sigmaVertexZ,
         maxDevVertex, sigmaTime, maxDevTime, offsetX, offsetY,
         offsetZ, offsetT;

  // Pointer to the random number generator.
  Rndm*  rndmPtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_BeamShape_H
