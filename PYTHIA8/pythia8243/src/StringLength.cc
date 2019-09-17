// StringLength.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// StringLength class.

// Calculate the lambda measure for string and junctions.

#include "Pythia8/StringLength.h"

namespace Pythia8 {

//==========================================================================

// The StringLength class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Minimum energy of a parton and minimum angle between two partons.
// This is to avoid problems with infinities.
const double StringLength::TINY     = 1e-20;
const double StringLength::MINANGLE = 1e-7;

//--------------------------------------------------------------------------

void StringLength::init(Info* infoPtrIn, Settings& settings) {

  // Save pointers.
  infoPtr = infoPtrIn;

  // Store variables.
  m0         = settings.parm("ColourReconnection:m0");
  m0sqr      = pow2(m0);
  juncCorr   = settings.parm("ColourReconnection:junctionCorrection");
  sqrt2      = sqrt(2.);
  lambdaForm = settings.mode("ColourReconnection:lambdaForm");
}

//--------------------------------------------------------------------------

// Calculate string length for two indices in the event record.

double StringLength::getStringLength( Event& event, int i, int j) {

  // Find rest frame of particles.
  Vec4 p1 =  event[i].p();
  Vec4 p2 =  event[j].p();

  return getStringLength(p1, p2);
}

//--------------------------------------------------------------------------

// Calculate string length for two particles given their four-momenta.

double StringLength::getStringLength( Vec4 p1, Vec4 p2) {

  // Check for very small energies and angles.
  if (p1.e() < TINY || p2.e() < TINY || theta(p1,p2) < MINANGLE) return 1e9;

  // Boost to restframe.
  Vec4 pSum = p1 + p2;
  p1.bstback(pSum);
  p2.bstback(pSum);

  // Calculate string length.
  Vec4 p0( 0., 0., 0., 1.);

  return getLength(p1,p0) + getLength(p2,p0);
}

//--------------------------------------------------------------------------

// Calculate string length of a single particle.
// The first vector is the 4 vector of the particle.
// The second vector represents (1,0,0,0) in dipole restframe.

double StringLength::getLength(Vec4 p, Vec4 v, bool isJunc) {

  double m = m0;
  if (isJunc) m *= juncCorr;

  if      (lambdaForm == 0) return log (1. + sqrt2 * v * p/ m );
  else if (lambdaForm == 1) return log (1. + 2 * v * p/ m );
  else if (lambdaForm == 2) return log (2. * v * p / m);
  else                      return 1e9;
}

//--------------------------------------------------------------------------

// Calculate the length of a single junction given the 3 entries in the event.

double StringLength::getJuncLength( Event& event, int i, int j, int k) {

  if (i == j || i == k || j == k) return 1e9;

  Vec4 p1 = event[i].p();
  Vec4 p2 = event[j].p();
  Vec4 p3 = event[k].p();

  return getJuncLength(p1, p2, p3);
}
//--------------------------------------------------------------------------

// Calculate the length of a single junction given the 3 four-momenta.

double StringLength::getJuncLength(Vec4 p1, Vec4 p2, Vec4 p3) {

  // Check for very small energies and angles.
  if (p1.e() < TINY || p2.e() < TINY || p3.e() < TINY
    || theta(p1,p2) < MINANGLE || theta(p1,p3) < MINANGLE
    || theta(p2,p3) < MINANGLE) return 1e9;

  // Find the junction rest frame.
  RotBstMatrix MfromJRF1 = stringFragmentation.junctionRestFrame(p1,p2,p3);
  MfromJRF1.invert();
  Vec4 v1( 0., 0., 0., 1.);
  v1.rotbst(MfromJRF1);

  // Possible problem when the right system rest frame system is not found.
  if ( pow2(p1*v1) - p1*p1 < 0. || pow2(p2*v1) - p2*p2 < 0.
    || pow2(p3*v1) - p3*p3 < 0.) return 1e9;

  // Calcualte the junction length.
  return getLength(p1, v1, true) + getLength(p2, v1, true)
       + getLength(p3, v1, true);
}

//--------------------------------------------------------------------------

// Calculate the length of a double junction given the 4 entries in the event.
// The first two are expected to be quarks, the second two to be anti quarks.

double StringLength::getJuncLength( Event& event, int i, int j, int k, int l) {

  if (i == j || i == k || i == l || j == k || j == l || k == l) return 1e9;

  // Simple minimum check of lengths.
  double origLength = getStringLength(event, i, k)
    + getStringLength(event, j, l);
  double minLength  = getStringLength(event, i, j)
    + getStringLength(event, k, l);

  if (origLength < minLength) return minLength;

  Vec4 p1 = event[i].p();
  Vec4 p2 = event[j].p();
  Vec4 p3 = event[k].p();
  Vec4 p4 = event[l].p();

  return getJuncLength(p1, p2, p3, p4);
}

//--------------------------------------------------------------------------

// Calculate the length of a double junction given the 4 four-momenta.
// The first two are expected to be quarks, the second two to be anti quarks.

double StringLength::getJuncLength(Vec4 p1, Vec4 p2, Vec4 p3, Vec4 p4) {

  // Check for very small energies and angles.
  if ( p1.e() < TINY || p2.e() < TINY || p3.e() < TINY || p4.e() < TINY
    || theta(p1,p2) < MINANGLE || theta(p1,p3) < MINANGLE
    || theta(p1,p4) < MINANGLE || theta(p2,p3) < MINANGLE
    || theta(p2,p4) < MINANGLE || theta(p3,p4) < MINANGLE) return 1e9;

  // Calculate velocity of first junction.
  Vec4 pSum1 = p3 +p4;

  RotBstMatrix MfromJRF1 = stringFragmentation.junctionRestFrame(p1,p2,pSum1);
  MfromJRF1.invert();
  Vec4 v1( 0., 0., 0., 1.);
  v1.rotbst(MfromJRF1);

  // Calculate velocity of second junction.
  Vec4 pSum2 = p1 + p2;

  RotBstMatrix MfromJRF2 = stringFragmentation.junctionRestFrame(p3,p4,pSum2);
  MfromJRF2.invert();
  Vec4 v2( 0., 0., 0., 1.);
  v2.rotbst(MfromJRF2);

  // This only happens if it is not possible to find the correct rest frame.
  if ( pow2(p1*v1) - p1*p1 < 0. || pow2(p2*v1) - p2*p2 < 0.
    || pow2(p3*v2) - p3*p3 < 0. || pow2(p4*v2) - p4*p4 < 0.) return 1e9;

  return getLength(p1, v1, true) + getLength(p2, v1, true)
       + getLength(p3, v2, true) + getLength(p4, v2, true)
       + log(v1*v2 + sqrt(pow2(v1*v2)-1));
}

//==========================================================================

} // end namespace Pythia8
