// PartonVertex.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the PartonVertex class.

#include "Pythia8/PartonVertex.h"

namespace Pythia8 {

//==========================================================================

// The PartonVertex class.

//--------------------------------------------------------------------------

// Find relevant settings.

void PartonVertex::init() {

    doVertex      = settingsPtr->flag("PartonVertex:setVertex");
    modeVertex    = settingsPtr->mode("PartonVertex:modeVertex");
    rProton       = settingsPtr->parm("PartonVertex:ProtonRadius");
    pTmin         = settingsPtr->parm("PartonVertex:pTmin");
    widthEmission = settingsPtr->parm("PartonVertex:EmissionWidth");
    bScale        = 2.187 / (2. * rProton);

}

//--------------------------------------------------------------------------

// Select vertex for a beam (remnant) particle.

void PartonVertex::vertexBeam( int iNow, int iBeam, Event& event) {
  double xbA = -bNow/2.;
  double xbB = bNow/2.;
  // Don't do any rotation in phi for these simple models.
  if(iBeam == 0)
    event[iNow].vProd(xbA, 0., 0., 0.);
  else if(iBeam == 1)
    event[iNow].vProd(xbB, 0., 0. ,0.);
  else
    infoPtr->errorMsg("Error in PartonVertex:vertexBeam: Wrong beam index.");
}

//--------------------------------------------------------------------------

// Select vertex or vertices for an MPI.

void PartonVertex::vertexMPI( int iBeg, int nAdd, double bNowIn,
  Event& event) {

  // Skip if not implemented option.
  if (!doVertex || modeVertex < 1 || modeVertex > 2) return;

  // Convert the impact parameter to physical units. Prepare selection.
  bNow = bNowIn / bScale;
  if (modeVertex == 1) {
    xMax = rProton - bNow / 2.;
    yMax = sqrt( 4. * rProton * rProton - bNow * bNow);
  } else if (modeVertex == 2) {
    mux = bNow / 2.0;
  }

  // Partons are given separate vertices for current models.
  for (int iNow = iBeg; iNow < iBeg + nAdd; ++iNow) {
    double x = 0.0, y = 0.0;

    // Sample x and y inside a box, and then require it to be within disks.
    if (modeVertex == 1) {
      bool reject = true;
      while (reject) {
        x = (2. * rndmPtr->flat() - 1.) * xMax;
        y = (2. * rndmPtr->flat() - 1.) * yMax;
        if ( (pow2(x + bNow / 2) + y * y < pow2(rProton))
          && (pow2(x - bNow / 2) + y * y < pow2(rProton)) ) reject = false;
      }

    // Sample x and y according to two-dimensional Gaussian.
    } else {
      pair<double,double> xy = rndmPtr->gauss2();
      x = (rProton / 2.0) * (xy.first + mux);
      y = (rProton / 2.0) * xy.second;
    }

    // Set production vertices.
    event[iNow].vProd( x, y, 0., 0.);
  }

}

//--------------------------------------------------------------------------

// Select vertex for an FSR branching.

void PartonVertex::vertexFSR( int iNow, Event& event) {

  // Skip if not implemented option.
  if (!doVertex || modeVertex < 1 || modeVertex > 2) return;

  // Mother index.
  int iMo = event[iNow].mother1();
  // Start from known vertex, or mother one.
  Vec4 vStart = event[iNow].hasVertex() ? event[iNow].vProd()
              : event[iMo].vProd();

  // Add Gaussian smearing.
  double pT = max( event[iNow].pT(), pTmin);
  pair<double, double> xy = rndmPtr->gauss2();
  Vec4 vSmear = (widthEmission / pT) * Vec4( xy.first, xy.second, 0., 0.);
  event[iNow].vProd( vStart + vSmear);


}

//--------------------------------------------------------------------------

// Select vertex for an ISR branching.

void PartonVertex::vertexISR( int iNow, Event& event) {

  // Skip if not implemented option.
  if (!doVertex || modeVertex < 1 || modeVertex > 2) return;

  // Start from known vertex or mother/daughter one.
  int iMoDa = event[iNow].mother1();
  if (iMoDa == 0) iMoDa = event[iNow].daughter1();
  Vec4 vStart = event[iNow].vProd();
  if (!event[iNow].hasVertex() && iMoDa != 0) vStart = event[iMoDa].vProd();

  // Add Gaussian smearing.
  double pT = max( event[iNow].pT(), pTmin);
  pair<double, double> xy = rndmPtr->gauss2();
  Vec4 vSmear = (widthEmission / pT) * Vec4( xy.first, xy.second, 0., 0.);
  event[iNow].vProd( vStart + vSmear);

}

//==========================================================================

} // end namespace Pythia8
