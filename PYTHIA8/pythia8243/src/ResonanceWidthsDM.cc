// ResonanceWidthsDM.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for DM resonance properties
// in the ResonanceS, ResonanceZp, ResonanceS1, ResonanceCha, ResonanceDM2,
// and ResonanceChaD classes.

#include "Pythia8/ResonanceWidthsDM.h"
#include "Pythia8/PythiaComplex.h"

namespace Pythia8 {

//==========================================================================

// The ResonanceS class.
// Derived class for S properties. (DMmed(s=0), PDG id 54.)

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceS::initConstants() {

  // Locally stored properties and couplings.
  double vq = settingsPtr->parm("Sdm:vf");
  double vX = settingsPtr->parm("Sdm:vX");
  double aq = settingsPtr->parm("Sdm:af");
  double aX = settingsPtr->parm("Sdm:aX");

  gq = abs(aq) > 0 ? aq : vq;
  gX = abs(aX) > 0 ? aX : vX;

  if (abs(aX) > 0) pScalar = true;
  else pScalar = false;

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceS::calcPreFac(bool) {

  // Common coupling factors.
  preFac = 1.0 / (12.0 * M_PI * mRes);
  alpS   = couplingsPtr->alphaS(mHat * mHat);

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceS::calcWidth(bool ) {

  // Check that above threshold.
  if (ps == 0.) return;

  // Caclulate current width.
  double mRat2 = pow2(mf1 / mRes);
  double kinfac = (1 - 4 * mRat2) * (1. + 2 * mRat2);
  widNow = 0.;
  if (id1Abs < 7)   widNow = 3. * pow2(gq * mf1) * preFac * kinfac;
  if (id1Abs == 21) widNow = pow2(gq) * preFac * pow2(alpS / M_PI) * eta2gg();
  if (id1Abs == 52) widNow = pow2(gX * mf1) * preFac * kinfac;

}

//--------------------------------------------------------------------------

// Loop integral for H -> gg coupling.

double ResonanceS::eta2gg() {

  // Initial values.
  complex eta = complex(0., 0.);
  double  mLoop, epsilon, root, rootLog;
  complex phi, etaNow;

  // Loop over s, c, b, t quark flavours.
  for (int idNow = 3; idNow < 7; ++idNow) {
    mLoop   = particleDataPtr->m0(idNow);
    epsilon = pow2(2. * mLoop / mHat);

    // Value of loop integral.
    if (epsilon <= 1.) {
      root    = sqrt(1. - epsilon);
      rootLog = (epsilon < 1e-4) ? log(4. / epsilon - 2.)
              : log( (1. + root) / (1. - root) );
      phi     = complex( -0.25 * (pow2(rootLog) - pow2(M_PI)),
                0.5 * M_PI * rootLog );
    }
    else phi  = complex( pow2( asin(1. / sqrt(epsilon)) ), 0.);

    // Factors that depend on Higgs and flavour type.
    if (!pScalar) etaNow = -0.5 * epsilon
      * (complex(1., 0.) + (1. - epsilon) * phi);
    else etaNow = -0.5 * epsilon * phi;


    // Sum up contribution and return square of absolute value.
    eta += etaNow;
  }
  return (pow2(eta.real()) + pow2(eta.imag()));

}

//==========================================================================

// The ResonanceZp class.
// Derived class for Z'^0 properties. (DMmed(s=1), PDG id 55.)

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceZp::initConstants() {

  // Coupling type and common couplings.
  kinMix = settingsPtr->flag("Zp:kineticMixing");
  gZp    = settingsPtr->parm("Zp:gZp");
  eps    = settingsPtr->parm("Zp:epsilon");
  vX     = settingsPtr->parm("Zp:vX");
  aX     = settingsPtr->parm("Zp:aX");

  // SM fermion couplings for kinetic mixing case.
  if (kinMix) {
    vu = eps * (2./3. + couplingsPtr->vf(2));
    au = eps * couplingsPtr->af(2);
    vd = eps * (-1./3. + couplingsPtr->vf(1));
    ad = eps * couplingsPtr->af(1);
    vl = eps * (-1. + couplingsPtr->vf(11));
    al = eps * couplingsPtr->af(11);
    vv = eps * couplingsPtr->vf(12);
    av = eps * couplingsPtr->af(12);

  // SM fermion couplings set by user.
  } else {
    vu = settingsPtr->parm("Zp:vu");
    vd = settingsPtr->parm("Zp:vd");
    vl = settingsPtr->parm("Zp:vl");
    vv = settingsPtr->parm("Zp:vv");
    au = settingsPtr->parm("Zp:au");
    ad = settingsPtr->parm("Zp:ad");
    al = settingsPtr->parm("Zp:al");
    av = settingsPtr->parm("Zp:av");
  }

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceZp::calcPreFac(bool) {

  // Common coupling factors.
  preFac      = mRes / 12.0 / M_PI;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceZp::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0. || id1 * id2 > 0) return;

  // Couplings to Zp.
  double kinFacA  = pow3(ps);
  double kinFacV  = ps * (1. + 2. * mr1);
  double fac      = 0;
  widNow = 0.;
  if (id1Abs < 7) {
    if (id1Abs%2 == 0) fac = vu * vu * kinFacV + au * au * kinFacA;
    else fac = vd * vd * kinFacV + ad * ad * kinFacA;
    widNow *= 3;
  } else if (id1Abs > 10 && id1Abs < 17) {
    if (id1Abs%2 == 1) fac = vl * vl * kinFacV + al * al * kinFacA;
    else fac = vv * vv * kinFacV + av * av * kinFacA;
  } else if (id1Abs == 52) fac = vX * vX * kinFacV + aX * aX * kinFacA;

  // Set partial width.
  double coup = pow2(gZp);
  if (kinMix && id1Abs != 52)
    coup = couplingsPtr->alphaEM(mRes * mRes) * 4.0 * M_PI;
  widNow = coup * fac * preFac;

}

//==========================================================================

// The ResonanceSl class.
// Derived class for Sl properties. (Using PDG id 56.)

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceSl::initConstants() {

  // Locally stored properties and couplings.
  yuk[0] = 0.0;
  yuk[1] = settingsPtr->parm("DM:yuk1");
  yuk[2] = settingsPtr->parm("DM:yuk2");
  yuk[3] = settingsPtr->parm("DM:yuk3");

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceSl::calcPreFac(bool) {

  // Common coupling factors.
  preFac      =  1. / (mRes * 16.0 * M_PI);

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceSl::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // Set partial width.
  kinFac = (mRes * mRes - mf1 * mf1 - mf2 * mf2);
  double coup = 0.0;
  if(abs(id2) == 11) coup = yuk[1];
  if(abs(id2) == 13) coup = yuk[2];
  if(abs(id2) == 15) coup = yuk[3];
  widNow = pow2(coup) * preFac * kinFac * ps;

  return;

}

//==========================================================================

// The ResonanceCha class.
// Derived class for singly-charged partner properties. (Using PDG id 57.)

//--------------------------------------------------------------------------

void ResonanceCha::setMassMix(){

  // Check if using the Nplet model
  // (do not over-ride DM mass if using s-channel mediators).
  doDY = settingsPtr->flag("DM:qqbar2DY");
  if (!doDY) return;

  // Locally stored properties and couplings.
  double M1     = settingsPtr->parm("DM:M1");
  double M2     = settingsPtr->parm("DM:M2");
  int    type   = settingsPtr->mode("DM:Nplet");
  double Lambda = settingsPtr->parm("DM:Lambda");

  // Mixing parameters.
  double vev    = 174.0;
  mixing        = vev / Lambda;
  if (type > 1) mixing *= sqrt(2) * vev;
  if (type > 2) mixing *= pow2(vev) /pow2(Lambda) / sqrt(12);
  double term1  = sqrt(pow2(M2 - M1) + pow2(mixing));
  double sin2th = 0.5 * (1 - abs(M2 - M1) / term1);
  mixN1         = (M1 > M2) ? sqrt(sin2th) : sqrt(1. - sin2th);
  mixN2         = (M1 > M2) ? sqrt(1. - sin2th) : sqrt(sin2th);

  // Set particle masses.
  double m1m    = 0.5 * (M1 + M2 - term1);
  double m2p    = 0.5 * (M1 + M2 + term1);
  double mplet  = (M1 < M2) ? m2p : m1m;
  // Neutral partners.
  particleDataPtr->m0(52, m1m);
  particleDataPtr->m0(58, m2p);
  // Chargino mass with NLO mass splitting 0.16 GeV. Doubly charged mass.
  particleDataPtr->m0(57, mplet + 0.16);
  particleDataPtr->m0(59, mplet + 0.16 + 0.49);

  return;

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceCha::calcPreFac(bool) {

  // Common coupling factors.
  preFac      =  mRes / (16.0 * M_PI);

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceCha::calcWidth(bool) {

  // Check if using the Nplet model.
  if (!doDY) return;

  // Small mass difference normal for compressed spectrum, so smaller
  // MASSMARGIN required. Check that above threshold.
  if (mHat < mf1 + mf2 + 0.01) return;

  // Mass parameters.
  double dm = 0.0;
  widNow = 0.0;
  double mix = (abs(id2) == 58) ? mixN2 : mixN1;
  double mpion = 0.1396;

  // Two-body decays.
  if (mult == 2) {
    dm = particleDataPtr->m0(57) - particleDataPtr->m0(abs(id2));
    // Long-lived two-body decay via pion.
    if (dm > mpion)  widNow =  2.0 * pow2(mix) * 6.993e-13
      * sqrt(1. - pow2(mpion/dm)) * pow3(dm);
    // Two-body decay via W.
    else if (dm > particleDataPtr->m0(24)) ;
    else return;

  // Three-body decays.
  } else { }

  return;

}

//==========================================================================

// The ResonanceDM2 class.
// Derived class for neutral partner properties. (Using PDG id 58.)
// Not yet implemented.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceDM2::initConstants() {

  setMassMix();
  mHiggs = particleDataPtr->m0(25);
  wHiggs = particleDataPtr->mWidth(25);

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceDM2::calcPreFac(bool) {

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceDM2::calcWidth(bool) {

  // Check if using the Nplet model
  if (!doDY) return;

  // Check that above threshold.
  if (ps == 0.) return;

  return;

}

//==========================================================================

// The ResonanceChaD class.
// Derived class for doubly charged properties. (Using PDG id 59.)

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceChaD::calcPreFac(bool) {

  // Common coupling factors.
  double dm = particleDataPtr->m0(59) - particleDataPtr->m0(57);
  preFac = (dm > 0.) ? 4.0 * 6.993e-13 * sqrtpos(1. - pow2(0.1396/dm))
         * pow3(dm) : 0.;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceChaD::calcWidth(bool) {

  // Check if using the Nplet model.
  if (!doDY) return;

  widNow = preFac;
  return;

}

//==========================================================================

} // end namespace Pythia8
