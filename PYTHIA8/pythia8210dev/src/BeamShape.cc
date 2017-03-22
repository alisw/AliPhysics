// BeamShape.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the BeamShape class.

#include "Pythia8/BeamShape.h"

namespace Pythia8 {

//==========================================================================

// The BeamShape class.

//--------------------------------------------------------------------------

// Initialize beam parameters.

  void BeamShape::init( Settings& settings, Rndm* rndmPtrIn) {

  // Save pointer.
  rndmPtr             = rndmPtrIn;

  // Main flags.
  allowMomentumSpread = settings.flag("Beams:allowMomentumSpread");
  allowVertexSpread   = settings.flag("Beams:allowVertexSpread");

  // Parameters for beam A momentum spread.
  sigmaPxA            = settings.parm("Beams:sigmaPxA");
  sigmaPyA            = settings.parm("Beams:sigmaPyA");
  sigmaPzA            = settings.parm("Beams:sigmaPzA");
  maxDevA             = settings.parm("Beams:maxDevA");

  // Parameters for beam B momentum spread.
  sigmaPxB            = settings.parm("Beams:sigmaPxB");
  sigmaPyB            = settings.parm("Beams:sigmaPyB");
  sigmaPzB            = settings.parm("Beams:sigmaPzB");
  maxDevB             = settings.parm("Beams:maxDevB");

  // Parameters for beam vertex spread.
  sigmaVertexX        = settings.parm("Beams:sigmaVertexX");
  sigmaVertexY        = settings.parm("Beams:sigmaVertexY");
  sigmaVertexZ        = settings.parm("Beams:sigmaVertexZ");
  maxDevVertex        = settings.parm("Beams:maxDevVertex");
  sigmaTime           = settings.parm("Beams:sigmaTime");
  maxDevTime          = settings.parm("Beams:maxDevTime");

  // Parameters for beam vertex offset.
  offsetX             = settings.parm("Beams:offsetVertexX");
  offsetY             = settings.parm("Beams:offsetVertexY");
  offsetZ             = settings.parm("Beams:offsetVertexZ");
  offsetT             = settings.parm("Beams:offsetTime");

}

//--------------------------------------------------------------------------

// Set the two beam momentum deviations and the beam vertex.

void BeamShape::pick() {

  // Reset all values.
  deltaPxA = deltaPyA = deltaPzA = deltaPxB = deltaPyB = deltaPzB
    = vertexX = vertexY = vertexZ = vertexT = 0.;

  // Set beam A momentum deviation by a three-dimensional Gaussian.
  if (allowMomentumSpread) {
    double totalDev, gauss;
    do {
      totalDev = 0.;
      if (sigmaPxA > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPxA  = sigmaPxA * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPyA > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPyA  = sigmaPyA * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPzA > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPzA  = sigmaPzA * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevA * maxDevA);

    // Set beam B momentum deviation by a three-dimensional Gaussian.
    do {
      totalDev = 0.;
      if (sigmaPxB > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPxB  = sigmaPxB * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPyB > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPyB  = sigmaPyB * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPzB > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPzB  = sigmaPzB * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevB * maxDevB);
  }

  // Set beam vertex location by a three-dimensional Gaussian.
  if (allowVertexSpread) {
    double totalDev, gauss;
    do {
      totalDev = 0.;
      if (sigmaVertexX > 0.) {
        gauss     = rndmPtr->gauss();
        vertexX   = sigmaVertexX * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaVertexY > 0.) {
        gauss     = rndmPtr->gauss();
        vertexY   = sigmaVertexY * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaVertexZ > 0.) {
        gauss     = rndmPtr->gauss();
        vertexZ   = sigmaVertexZ * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevVertex * maxDevVertex);

    // Set beam collision time by a Gaussian.
    if (sigmaTime > 0.) {
      do gauss    = rndmPtr->gauss();
      while (abs(gauss) > maxDevTime);
      vertexT     = sigmaTime * gauss;
    }

    // Add offset to beam vertex.
    vertexX      += offsetX;
    vertexY      += offsetY;
    vertexZ      += offsetZ;
    vertexT      += offsetT;
  }

}

//==========================================================================

} // end namespace Pythia8
