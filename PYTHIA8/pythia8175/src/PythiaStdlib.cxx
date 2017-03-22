// PythiaStdlib.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Gamma function.

#include "PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================
 
// The Gamma function for real arguments, using the Lanczos approximation.
// Code based on http://en.wikipedia.org/wiki/Lanczos_approximation

double GammaCoef[9] = {
     0.99999999999980993,     676.5203681218851,   -1259.1392167224028,
      771.32342877765313,   -176.61502916214059,    12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

double GammaReal(double x) {

  // Reflection formula (recursive!) for x < 0.5.
  if (x < 0.5) return M_PI / (sin(M_PI * x) * GammaReal(1 - x));

  // Iterate through terms.
  double z = x - 1.;
  double gamma = GammaCoef[0];
  for (int i = 1; i < 9; ++i) gamma += GammaCoef[i] / (z + i);

  // Answer.
  double t = z + 7.5;
  gamma *= sqrt(2. * M_PI) * pow(t, z + 0.5) * exp(-t); 
  return gamma;

}

//==========================================================================

} // end namespace Pythia8
