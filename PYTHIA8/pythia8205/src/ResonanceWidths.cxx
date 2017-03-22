// ResonanceWidths.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for
// the ResonanceWidths class and classes derived from it.

#include "Pythia8/ParticleData.h"
#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/PythiaComplex.h"

namespace Pythia8 {

//==========================================================================

// The ResonanceWidths class.
// Base class for the various resonances.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of points in integration direction for numInt routines.
const int    ResonanceWidths::NPOINT         = 100;

// The mass of a resonance must not be too small.
const double ResonanceWidths::MASSMIN        = 0.4;

// The sum of product masses must not be too close to the resonance mass.
const double ResonanceWidths::MASSMARGIN     = 0.1;

//--------------------------------------------------------------------------

// Initialize data members.
// Calculate and store partial and total widths at the nominal mass.

bool ResonanceWidths::init(Info* infoPtrIn, Settings* settingsPtrIn,
   ParticleData* particleDataPtrIn, Couplings* couplingsPtrIn) {

  // Save pointers.
  infoPtr         = infoPtrIn;
  settingsPtr     = settingsPtrIn;
  particleDataPtr = particleDataPtrIn;
  couplingsPtr    = couplingsPtrIn;

  // Perform any model dependent initialisations (pure dummy in base class).
  bool isInit = initBSM();

  // Minimal decaying-resonance width. Minimal phase space for meMode = 103.
  minWidth     = settingsPtr->parm("ResonanceWidths:minWidth");
  minThreshold = settingsPtr->parm("ResonanceWidths:minThreshold");

  // Pointer to particle species.
  particlePtr  = particleDataPtr->particleDataEntryPtr(idRes);
  if (particlePtr == 0) {
    infoPtr->errorMsg("Error in ResonanceWidths::init:"
      " unknown resonance identity code");
    return false;
  }

  // Generic particles should not have meMode < 100, but allow
  // some exceptions: not used Higgses and not used Technicolor.
  if (idRes == 35 || idRes == 36 || idRes == 37
    || idRes/1000000 == 3) isGeneric = false;

  // Resonance properties: antiparticle, mass, width
  hasAntiRes   = particlePtr->hasAnti();
  mRes         = particlePtr->m0();
  GammaRes     = particlePtr->mWidth();
  m2Res        = mRes*mRes;

  // A resonance cannot be too light, in particular not = 0.
  if (mRes < MASSMIN) {
    ostringstream idCode;
    idCode << idRes;
    infoPtr->errorMsg("Error in ResonanceWidths::init:"
      " resonance mass too small", "for id = " + idCode.str(), true);
    return false;
  }

  // For very narrow resonances assign fictitious small width.
  if (GammaRes < minWidth) GammaRes = 0.1 * minWidth;
  GamMRat      = (mRes == 0.) ? 0. : GammaRes / mRes;

  // Secondary decay chains by default all on.
  openPos      = 1.;
  openNeg      = 1.;

  // Allow option where on-shell width is forced to current value.
  // Disable for mixes gamma*/Z0/Z'0
  doForceWidth = particlePtr->doForceWidth();
  if (idRes == 23 && settingsPtr->mode("WeakZ0:gmZmode") != 2)
    doForceWidth = false;
  if (idRes == 33 && settingsPtr->mode("Zprime:gmZmode") != 3)
    doForceWidth = false;
  forceFactor  = 1.;

  // Check if we are supposed to do the width calculation
  // (can be false e.g. if SLHA decay table should take precedence instead).
  allowCalcWidth = isInit && allowCalc();
  if ( allowCalcWidth ) {
    // Initialize constants used for a resonance.
    initConstants();

    // Calculate various common prefactors for the current mass.
    mHat          = mRes;
    calcPreFac(true);
  }

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.;
  double widPos = 0.;
  double widNeg = 0.;
  int    idNow, idAnti;
  double openSecPos, openSecNeg;

  // Loop over all decay channels. Basic properties of channel.
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    iChannel    = i;
    onMode      = particlePtr->channel(i).onMode();
    meMode      = particlePtr->channel(i).meMode();
    mult        = particlePtr->channel(i).multiplicity();
    widNow      = 0.;

    // Warn if not relevant meMode.
    if ( meMode < 0 || meMode > 103 || (isGeneric && meMode < 100) ) {
      stringstream ssIdRes;
      ssIdRes << "for " << idRes;
      infoPtr->errorMsg("Error in ResonanceWidths::init:"
        " resonance meMode not acceptable", ssIdRes.str());
    }

    // Channels with meMode < 100 must be implemented in derived classes.
    if (meMode < 100 && allowCalcWidth) {

      // Read out information on channel: primarily use first two.
      id1       = particlePtr->channel(i).product(0);
      id2       = particlePtr->channel(i).product(1);
      id1Abs    = abs(id1);
      id2Abs    = abs(id2);

      // Order first two in descending order of absolute values.
      if (id2Abs > id1Abs) {swap( id1, id2); swap( id1Abs, id2Abs);}

      // Allow for third product to be treated in derived classes.
      if (mult > 2) {
        id3     = particlePtr->channel(i).product(2);
        id3Abs  = abs(id3);

        // Also order third into descending order of absolute values.
        if (id3Abs > id2Abs) {swap( id2, id3); swap( id2Abs, id3Abs);}
        if (id2Abs > id1Abs) {swap( id1, id2); swap( id1Abs, id2Abs);}
      }

      // Read out masses. Calculate two-body phase space.
      mf1       = particleDataPtr->m0(id1Abs);
      mf2       = particleDataPtr->m0(id2Abs);
      mr1       = pow2(mf1 / mHat);
      mr2       = pow2(mf2 / mHat);
      ps        = (mHat < mf1 + mf2 + MASSMARGIN) ? 0.
                : sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );
      if (mult > 2) {
        mf3     = particleDataPtr->m0(id3Abs);
        mr3     = pow2(mf3 / mHat);
      }

      // Let derived class calculate width for channel provided.
      calcWidth(true);
    }

    // Channels with meMode >= 100 are calculated based on stored values.
    else widNow = GammaRes * particlePtr->channel(i).bRatio();

    // Find secondary open fractions of partial width.
    openSecPos  = 1.;
    openSecNeg  = 1.;
    if (widNow > 0.) for (int j = 0; j < mult; ++j) {
      idNow     = particlePtr->channel(i).product(j);
      idAnti    = (particleDataPtr->hasAnti(idNow)) ? -idNow : idNow;
      // Secondary widths not yet initialized for heavier states,
      // so have to assume unit open fraction there.
      if (idNow == 23 || abs(idNow) == 24 || idNow == 93 || abs(idNow) == 94
        || particleDataPtr->m0(abs(idNow)) < mRes) {
        openSecPos *= particleDataPtr->resOpenFrac(idNow);
        openSecNeg *= particleDataPtr->resOpenFrac(idAnti);
      }
    }

    // Store partial widths and secondary open fractions.
    particlePtr->channel(i).onShellWidth(widNow);
    particlePtr->channel(i).openSec( idRes, openSecPos);
    particlePtr->channel(i).openSec(-idRes, openSecNeg);

    // Update sum over all channnels and over open channels only.
    widTot     += widNow;
    if (onMode == 1 || onMode == 2) widPos += widNow * openSecPos;
    if (onMode == 1 || onMode == 3) widNeg += widNow * openSecNeg;
  }

  // If no decay channels are open then set particle stable and done.
  if (widTot < minWidth) {
    particlePtr->setMayDecay(false, false);
    particlePtr->setMWidth(0., false);
    for (int i = 0; i < particlePtr->sizeChannels(); ++i)
      particlePtr->channel(i).bRatio( 0., false);
    return true;
  }

  // Normalize branching ratios to unity.
  double bRatio;
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    bRatio      = particlePtr->channel(i).onShellWidth() / widTot;
    particlePtr->channel(i).bRatio( bRatio, false);
  }

  // Optionally force total width by rescaling of all partial ones.
  if (doForceWidth) {
    forceFactor = GammaRes / widTot;
    for (int i = 0; i < particlePtr->sizeChannels(); ++i)
      particlePtr->channel(i).onShellWidthFactor( forceFactor);
  }

  // Else update newly calculated partial width.
  else {
    particlePtr->setMWidth(widTot, false);
    GammaRes    = widTot;
  }

  // Updated width-to-mass ratio. Secondary widths for open.
  GamMRat       = GammaRes / mRes;
  openPos       = widPos / widTot;
  openNeg       = widNeg / widTot;

  // Clip wings of Higgses.
  bool isHiggs = (idRes == 25 || idRes == 35 ||idRes == 36 ||idRes == 37);
  bool clipHiggsWings = settingsPtr->flag("Higgs:clipWings");
  if (isHiggs && clipHiggsWings) {
    double mMinNow  = particlePtr->mMin();
    double mMaxNow  = particlePtr->mMax();
    double wingsFac = settingsPtr->parm("Higgs:wingsFac");
    double mMinWing = mRes - wingsFac * GammaRes;
    double mMaxWing = mRes + wingsFac * GammaRes;
    if (mMinWing > mMinNow) particlePtr->setMMinNoChange(mMinWing);
    if (mMaxWing < mMaxNow || mMaxNow < mMinNow)
      particlePtr->setMMaxNoChange(mMaxWing);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Calculate the total width and store phase-space-weighted coupling sums.

double ResonanceWidths::width(int idSgn, double mHatIn, int idInFlavIn,
  bool openOnly, bool setBR, int idOutFlav1, int idOutFlav2) {

  // Calculate various prefactors for the current mass.
  mHat          = mHatIn;
  idInFlav      = idInFlavIn;
  if (allowCalcWidth) calcPreFac(false);

  // Reset quantities to sum. Declare variables inside loop.
  double widSum = 0.;
  double mfSum, psOnShell;

  // Loop over all decay channels. Basic properties of channel.
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    iChannel    = i;
    onMode      = particlePtr->channel(i).onMode();
    meMode      = particlePtr->channel(i).meMode();
    mult        = particlePtr->channel(i).multiplicity();

    // Initially assume vanishing branching ratio.
    widNow      = 0.;
    if (setBR) particlePtr->channel(i).currentBR(widNow);

    // Optionally only consider specific (two-body) decay channel.
    // Currently only used for Higgs -> q qbar, g g or gamma gamma.
    if (idOutFlav1 > 0 || idOutFlav2 > 0) {
      if (mult > 2) continue;
      if (particlePtr->channel(i).product(0) != idOutFlav1) continue;
      if (particlePtr->channel(i).product(1) != idOutFlav2) continue;
    }

    // Optionally only consider open channels.
    if (openOnly) {
      if (idSgn > 0 && onMode !=1 && onMode != 2) continue;
      if (idSgn < 0 && onMode !=1 && onMode != 3) continue;
    }

    // Channels with meMode < 100 must be implemented in derived classes.
    if (meMode < 100) {

      // Read out information on channel: primarily use first two.
      id1       = particlePtr->channel(i).product(0);
      id2       = particlePtr->channel(i).product(1);
      id1Abs    = abs(id1);
      id2Abs    = abs(id2);

      // Order first two in descending order of absolute values.
      if (id2Abs > id1Abs) {swap( id1, id2); swap( id1Abs, id2Abs);}

      // Allow for third product to be treated in derived classes.
      if (mult > 2) {
        id3     = particlePtr->channel(i).product(2);
        id3Abs  = abs(id3);

        // Also order third into descending order of absolute values.
        if (id3Abs > id2Abs) {swap( id2, id3); swap( id2Abs, id3Abs);}
        if (id2Abs > id1Abs) {swap( id1, id2); swap( id1Abs, id2Abs);}
      }

      // Read out masses. Calculate two-body phase space.
      mf1       = particleDataPtr->m0(id1Abs);
      mf2       = particleDataPtr->m0(id2Abs);
      mr1       = pow2(mf1 / mHat);
      mr2       = pow2(mf2 / mHat);
      ps        = (mHat < mf1 + mf2 + MASSMARGIN) ? 0.
                : sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );
      if (mult > 2) {
        mf3     = particleDataPtr->m0(id3Abs);
        mr3     = pow2(mf3 / mHat);
      }

      // Let derived class calculate width for channel provided.
      calcWidth(false);
    }

    // Now on to meMode >= 100. First case: no correction at all.
    else if (meMode == 100)
      widNow    = GammaRes * particlePtr->channel(i).bRatio();

    // Correction by step at threshold.
    else if (meMode == 101) {
      mfSum     = 0.;
      for (int j = 0; j < mult; ++j) mfSum
               += particleDataPtr->m0( particlePtr->channel(i).product(j) );
      if (mfSum + MASSMARGIN < mHat)
        widNow  = GammaRes * particlePtr->channel(i).bRatio();
    }

    // Correction by a phase space factor for two-body decays.
    else if ( (meMode == 102 || meMode == 103) && mult == 2) {
      mf1       = particleDataPtr->m0( particlePtr->channel(i).product(0) );
      mf2       = particleDataPtr->m0( particlePtr->channel(i).product(1) );
      mr1       = pow2(mf1 / mHat);
      mr2       = pow2(mf2 / mHat);
      ps        = (mHat < mf1 + mf2 + MASSMARGIN) ? 0.
                : sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );
      mr1       = pow2(mf1 / mRes);
      mr2       = pow2(mf2 / mRes);
      psOnShell = (meMode == 102) ? 1. : max( minThreshold,
                  sqrtpos( pow2(1.- mr1 - mr2) - 4. * mr1 * mr2) );
      widNow = GammaRes * particlePtr->channel(i).bRatio() * ps / psOnShell;
    }

    // Correction by simple threshold factor for multibody decay.
    else if (meMode == 102 || meMode == 103) {
      mfSum  = 0.;
      for (int j = 0; j < mult; ++j) mfSum
               += particleDataPtr->m0( particlePtr->channel(i).product(j) );
      ps        = sqrtpos(1. - mfSum / mHat);
      psOnShell = (meMode == 102) ? 1. : max( minThreshold,
                  sqrtpos(1. - mfSum / mRes) );
      widNow = GammaRes * particlePtr->channel(i).bRatio() * ps / psOnShell;
    }

    // Optionally multiply by secondary widths.
    if (openOnly) widNow *= particlePtr->channel(i).openSec(idSgn);

    // Optionally include factor to force to fixed width.
    if (doForceWidth) widNow *= forceFactor;

    // Optionally multiply by current/nominal resonance mass??

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->channel(i).currentBR(widNow);
  }

  // Done.
  return widSum;

}

//--------------------------------------------------------------------------

// Numerical integration of matrix-element in two-body decay,
// where one particle is described by a Breit-Wigner mass distribution.
// Normalization to unit integral if matrix element is unity
// and there are no phase-space restrictions.

double ResonanceWidths::numInt1BW(double mHatIn, double m1, double Gamma1,
  double mMin1, double m2, int psMode) {

  // Check that phase space is open for integration.
  if (mMin1 + m2 > mHatIn) return 0.;

  // Precalculate coefficients for Breit-Wigner selection.
  double s1       = m1 * m1;
  double mG1      = m1 * Gamma1;
  double mMax1    = mHatIn - m2;
  double atanMin1 = atan( (mMin1 * mMin1 - s1) / mG1 );
  double atanMax1 = atan( (mMax1 * mMax1 - s1) / mG1 );
  double atanDif1 = atanMax1 - atanMin1;
  double wtDif1   = atanDif1 / (M_PI * NPOINT);

  // Step size in atan-mapped variable.
  double xStep    = 1. / NPOINT;

  // Variables used in loop over integration points.
  double sum      = 0.;
  double mrNow2   = pow2(m2 / mHatIn);
  double xNow1, sNow1, mNow1, mrNow1, psNow, value;

  // Loop with first-particle mass selection.
  for (int ip1 = 0; ip1 < NPOINT; ++ip1) {
    xNow1         = xStep * (ip1 + 0.5);
    sNow1         = s1 + mG1 * tan(atanMin1 + xNow1 * atanDif1);
    mNow1         = min( mMax1, max( mMin1, sqrtpos(sNow1) ) );
    mrNow1        = pow2(mNow1 / mHatIn);

    // Evaluate value and add to sum. Different matrix elements.
    psNow         = sqrtpos( pow2(1. - mrNow1 - mrNow2)
                    - 4. * mrNow1 * mrNow2);
    value         = 1.;
    if      (psMode == 1) value = psNow;
    else if (psMode == 2) value = psNow * psNow;
    else if (psMode == 3) value = pow3(psNow);
    else if (psMode == 5) value = psNow
      * (pow2(1. - mrNow1 - mrNow2) + 8. * mrNow1 * mrNow2);
    else if (psMode == 6) value = pow3(psNow);
    sum          += value;

  // End of  loop over integration points. Overall normalization.
  }
  sum            *= wtDif1;

  // Done.
  return sum;
}

//--------------------------------------------------------------------------

// Numerical integration of matrix-element in two-body decay,
// where both particles are described by Breit-Wigner mass distributions.
// Normalization to unit integral if matrix element is unity
// and there are no phase-space restrictions.

double ResonanceWidths::numInt2BW(double mHatIn, double m1, double Gamma1,
  double mMin1, double m2, double Gamma2, double mMin2, int psMode) {

  // Check that phase space is open for integration.
  if (mMin1 + mMin2 >= mHatIn) return 0.;

  // Precalculate coefficients for Breit-Wigner selection.
  double s1       = m1 * m1;
  double mG1      = m1 * Gamma1;
  double mMax1    = mHatIn - mMin2;
  double atanMin1 = atan( (mMin1 * mMin1 - s1) / mG1 );
  double atanMax1 = atan( (mMax1 * mMax1 - s1) / mG1 );
  double atanDif1 = atanMax1 - atanMin1;
  double wtDif1   = atanDif1 / (M_PI * NPOINT);
  double s2       = m2 * m2;
  double mG2      = m2 * Gamma2;
  double mMax2    = mHatIn - mMin1;
  double atanMin2 = atan( (mMin2 * mMin2 - s2) / mG2 );
  double atanMax2 = atan( (mMax2 * mMax2 - s2) / mG2 );
  double atanDif2 = atanMax2 - atanMin2;
  double wtDif2   = atanDif2 / (M_PI * NPOINT);

  // If on-shell decay forbidden then split integration range
  // to ensure that low-mass region is not forgotten.
  bool mustDiv    = false;
  double mDiv1    = 0.;
  double atanDiv1 = 0.;
  double atanDLo1 = 0.;
  double atanDHi1 = 0.;
  double wtDLo1   = 0.;
  double wtDHi1   = 0.;
  double mDiv2    = 0.;
  double atanDiv2 = 0.;
  double atanDLo2 = 0.;
  double atanDHi2 = 0.;
  double wtDLo2   = 0.;
  double wtDHi2   = 0.;
  if (m1 + m2 > mHatIn) {
    mustDiv       = true;
    double tmpDiv = (mHatIn - m1 - m2) / (Gamma1 + Gamma2);
    mDiv1         = m1 + Gamma1 * tmpDiv;
    atanDiv1      = atan( (mDiv1 * mDiv1 - s1) / mG1 );
    atanDLo1      = atanDiv1 - atanMin1;
    atanDHi1      = atanMax1 - atanDiv1;
    wtDLo1        = atanDLo1 / (M_PI * NPOINT);
    wtDHi1        = atanDHi1 / (M_PI * NPOINT);
    mDiv2         = m2 + Gamma2 * tmpDiv;
    atanDiv2      = atan( (mDiv2 * mDiv2 - s2) / mG2 );
    atanDLo2      = atanDiv2 - atanMin2;
    atanDHi2      = atanMax2 - atanDiv2;
    wtDLo2        = atanDLo2 / (M_PI * NPOINT);
    wtDHi2        = atanDHi2 / (M_PI * NPOINT);
  }

  // Step size in atan-mapped variable.
  double xStep    = 1. / NPOINT;
  int nIter       = (mustDiv) ? 2 * NPOINT : NPOINT;

  // Variables used in loop over integration points.
  double sum      = 0.;
  double xNow1, sNow1, mNow1, mrNow1, xNow2, sNow2, mNow2, mrNow2, psNow,
         value;
  double wtNow1   = wtDif1;
  double wtNow2   = wtDif2;

  // Outer loop with first-particle mass selection.
  for (int ip1 = 0; ip1 < nIter; ++ip1) {
    if (!mustDiv) {
      xNow1       = xStep * (ip1 + 0.5);
      sNow1       = s1 + mG1 * tan(atanMin1 + xNow1 * atanDif1);
    } else if (ip1 < NPOINT) {
      xNow1       = xStep * (ip1 + 0.5);
      sNow1       = s1 + mG1 * tan(atanMin1 + xNow1 * atanDLo1);
      wtNow1      = wtDLo1;
    } else {
      xNow1       = xStep * (ip1 - NPOINT + 0.5);
      sNow1       = s1 + mG1 * tan(atanDiv1 + xNow1 * atanDHi1);
      wtNow1      = wtDHi1;
    }
    mNow1         = min( mMax1, max( mMin1, sqrtpos(sNow1) ) );
    mrNow1        = pow2(mNow1 / mHatIn);

    // Inner loop with second-particle mass selection.
    for (int ip2 = 0; ip2 < nIter; ++ip2) {
      if (!mustDiv) {
        xNow2     = xStep * (ip2 + 0.5);
        sNow2     = s2 + mG2 * tan(atanMin2 + xNow2 * atanDif2);
      } else if (ip2 < NPOINT) {
        xNow2     = xStep * (ip2 + 0.5);
        sNow2     = s2 + mG2 * tan(atanMin2 + xNow2 * atanDLo2);
        wtNow2    = wtDLo2;
      } else {
        xNow2     = xStep * (ip2 - NPOINT + 0.5);
        sNow2     = s2 + mG2 * tan(atanDiv2 + xNow2 * atanDHi2);
        wtNow2    = wtDHi2;
      }
      mNow2       = min( mMax2, max( mMin2, sqrtpos(sNow2) ) );
      mrNow2      = pow2(mNow2 / mHatIn);

      // Check that point is inside phase space.
      if (mNow1 + mNow2 > mHatIn) break;

      // Evaluate value and add to sum. Different matrix elements.
      psNow       = sqrtpos( pow2(1. - mrNow1 - mrNow2)
                    - 4. * mrNow1 * mrNow2);
      value       = 1.;
      if      (psMode == 1) value = psNow;
      else if (psMode == 2) value = psNow * psNow;
      else if (psMode == 3) value = pow3(psNow);
      else if (psMode == 5) value = psNow
        * (pow2(1. - mrNow1 - mrNow2) + 8. * mrNow1 * mrNow2);
      else if (psMode == 6) value = pow3(psNow);
      sum        += value * wtNow1 * wtNow2;

    // End of second and first loop over integration points.
    }
  }

  // Done.
  return sum;
}

//==========================================================================

// The ResonanceGmZ class.
// Derived class for gamma*/Z0 properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceGmZ::initConstants() {

  // Locally stored properties and couplings.
  gmZmode     = settingsPtr->mode("WeakZ0:gmZmode");
  thetaWRat   = 1. / (16. * couplingsPtr->sin2thetaW()
                * couplingsPtr->cos2thetaW());

  // The Z0copy with id = 93 is a pure Z0.
  if (idRes == 93) gmZmode = 2;

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceGmZ::calcPreFac(bool calledFromInit) {

  // Common coupling factors.
  alpEM       = couplingsPtr->alphaEM(mHat * mHat);
  alpS        = couplingsPtr->alphaS(mHat * mHat);
  colQ        = 3. * (1. + alpS / M_PI);
  preFac      = alpEM * thetaWRat * mHat / 3.;

  // When call for incoming flavour need to consider gamma*/Z0 mix.
  if (!calledFromInit) {

    // Couplings when an incoming fermion is specified; elso only pure Z0.
    ei2       = 0.;
    eivi      = 0.;
    vi2ai2    = 1.;
    int idInFlavAbs = abs(idInFlav);
    if (idInFlavAbs > 0 && idInFlavAbs < 19) {
      ei2     = couplingsPtr->ef2(idInFlavAbs);
      eivi    = couplingsPtr->efvf(idInFlavAbs);
      vi2ai2  = couplingsPtr->vf2af2(idInFlavAbs);
    }

    // Calculate prefactors for gamma/interference/Z0 terms.
    double sH = mHat * mHat;
    gamNorm   = ei2;
    intNorm   = 2. * eivi * thetaWRat * sH * (sH - m2Res)
              / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
    resNorm   = vi2ai2 * pow2(thetaWRat * sH)
              / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

    // Optionally only keep gamma* or Z0 term.
    if (gmZmode == 1) {intNorm = 0.; resNorm = 0.;}
    if (gmZmode == 2) {gamNorm = 0.; intNorm = 0.;}
  }

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceGmZ::calcWidth(bool calledFromInit) {

  // Check that above threshold.
  if (ps == 0.) return;

  // Only contributions from three fermion generations, except top.
  if ( (id1Abs > 5 && id1Abs < 11) || id1Abs > 16 ) return;

  // At initialization only the pure Z0 should be considered.
  if (calledFromInit) {

    // Combine kinematics with colour factor and couplings.
    widNow  = preFac * ps * (couplingsPtr->vf2(id1Abs) * (1. + 2. * mr1)
            + couplingsPtr->af2(id1Abs) * ps*ps);
    if (id1Abs < 6) widNow *= colQ;
  }

  // When call for incoming flavour need to consider gamma*/Z0 mix.
  else {

    // Kinematical factors and couplings.
    double kinFacV  = ps * (1. + 2. * mr1);
    double ef2      = couplingsPtr->ef2(id1Abs) * kinFacV;
    double efvf     = couplingsPtr->efvf(id1Abs) * kinFacV;
    double vf2af2   = couplingsPtr->vf2(id1Abs) * kinFacV
                    + couplingsPtr->af2(id1Abs) * pow3(ps);

    // Relative outwidths: combine instate, propagator and outstate.
    widNow = gamNorm * ef2 + intNorm * efvf + resNorm * vf2af2;

    // Colour factor.
    if (id1Abs < 6) widNow *= colQ;
  }

}

//==========================================================================

// The ResonanceW class.
// Derived class for W+- properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceW::initConstants() {

  // Locally stored properties and couplings.
  thetaWRat = 1. / (12. * couplingsPtr->sin2thetaW());

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceW::calcPreFac(bool) {

  // Common coupling factors.
  alpEM     = couplingsPtr->alphaEM(mHat * mHat);
  alpS      = couplingsPtr->alphaS(mHat * mHat);
  colQ      = 3. * (1. + alpS / M_PI);
  preFac    = alpEM * thetaWRat * mHat;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceW::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // Only contributions from three fermion generations, except top.
  if ( (id1Abs > 5 && id1Abs < 11) || id1Abs > 16 ) return;


  // Combine kinematics with colour factor and couplings.
  widNow    = preFac * ps
            * (1. - 0.5 * (mr1 + mr2) - 0.5 * pow2(mr1 - mr2));
  if (id1Abs < 6) widNow *= colQ * couplingsPtr->V2CKMid(id1Abs, id2Abs);

}

//==========================================================================

// The ResonanceTop class.
// Derived class for top/antitop properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceTop::initConstants() {

  // Locally stored properties and couplings.
  thetaWRat = 1. / (16. * couplingsPtr->sin2thetaW());
  m2W       = pow2(particleDataPtr->m0(24));

  // Extra coupling factors for t -> H+ + b.
  tanBeta   = settingsPtr->parm("HiggsHchg:tanBeta");
  tan2Beta  = tanBeta * tanBeta;
  mbRun     = particleDataPtr->mRun( 5, particleDataPtr->m0(6) );

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceTop::calcPreFac(bool) {

  // Common coupling factors.
  alpEM     = couplingsPtr->alphaEM(mHat * mHat);
  alpS      = couplingsPtr->alphaS(mHat * mHat);
  colQ      = 1. - 2.5 * alpS / M_PI;
  preFac    = alpEM * thetaWRat * pow3(mHat) / m2W;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceTop::calcWidth(bool) {


  // Check that above threshold.
  if (ps == 0.) return;

  // Contributions from W + quark.
  if (id1Abs == 24 && id2Abs < 6) {
    widNow  = preFac * ps
            * ( pow2(1. - mr2) + (1. + mr2) * mr1 - 2. * mr1 * mr1 );

    // Combine with colour factor and CKM couplings.
    widNow *= colQ * couplingsPtr->V2CKMid(6, id2Abs);

  // Contributions from H+ + quark (so far only b).
  } else if (id1Abs == 37 && id2Abs == 5) {
    widNow  = preFac * ps * ( (1. + mr2 - mr1)
            * (pow2(mbRun / mHat) * tan2Beta + 1. / tan2Beta)
            + 4. * mbRun * mf2 / pow2(mHat) );
  }

}

//==========================================================================

// The ResonanceFour class.
// Derived class for fourth-generation properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceFour::initConstants() {

  // Locally stored properties and couplings.
  thetaWRat = 1. / (16. * couplingsPtr->sin2thetaW());
  m2W       = pow2(particleDataPtr->m0(24));

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceFour::calcPreFac(bool) {

  // Common coupling factors.
  alpEM     = couplingsPtr->alphaEM(mHat * mHat);
  alpS      = couplingsPtr->alphaS(mHat * mHat);
  colQ      = (idRes < 9) ? 1. - 2.5 * alpS / M_PI : 1.;
  preFac    = alpEM * thetaWRat * pow3(mHat) / m2W;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceFour::calcWidth(bool) {

  // Only contributions from W + fermion.
  if (id1Abs != 24 || id2Abs > 18) return;

  // Check that above threshold. Kinematical factor.
  if (ps == 0.) return;
  widNow    = preFac * ps
            * ( pow2(1. - mr2) + (1. + mr2) * mr1 - 2. * mr1 * mr1 );

  // Combine with colour factor and CKM couplings.
  if (idRes < 9) widNow *= colQ * couplingsPtr->V2CKMid(idRes, id2Abs);

}

//==========================================================================

// The ResonanceH class.
// Derived class for SM and BSM Higgs properties.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Minimal mass for W, Z, top in integration over respective Breit-Wigner.
// Top constrainted by t -> W b decay, which is not seen in simple top BW.
const double ResonanceH::MASSMINWZ = 10.;
const double ResonanceH::MASSMINT  = 100.;

// Number of widths above threshold where B-W integration not needed.
const double ResonanceH::GAMMAMARGIN = 10.;

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceH::initConstants() {

  // Locally stored properties and couplings.
  useCubicWidth  = settingsPtr->flag("Higgs:cubicWidth");
  useRunLoopMass = settingsPtr->flag("Higgs:runningLoopMass");
  sin2tW         = couplingsPtr->sin2thetaW();
  cos2tW         = 1. - sin2tW;
  mT             = particleDataPtr->m0(6);
  mZ             = particleDataPtr->m0(23);
  mW             = particleDataPtr->m0(24);
  mHchg          = particleDataPtr->m0(37);
  GammaT         = particleDataPtr->mWidth(6);
  GammaZ         = particleDataPtr->mWidth(23);
  GammaW         = particleDataPtr->mWidth(24);

  // NLO corrections to SM Higgs width, rescaled to reference alpha_S value.
  useNLOWidths   = (higgsType == 0) && settingsPtr->flag("HiggsSM:NLOWidths");
  rescAlpS       = 0.12833 / couplingsPtr->alphaS(125. * 125.);
  rescColQ       = 1.;

  // Couplings to fermions, Z and W, depending on Higgs type.
  coup2d         = 1.;
  coup2u         = 1.;
  coup2l         = 1.;
  coup2Z         = 1.;
  coup2W         = 1.;
  coup2Hchg      = 0.;
  coup2H1H1      = 0.;
  coup2A3A3      = 0.;
  coup2H1Z       = 0.;
  coup2A3Z       = 0.;
  coup2A3H1      = 0.;
  coup2HchgW     = 0.;
  if (higgsType == 1) {
    coup2d       = settingsPtr->parm("HiggsH1:coup2d");
    coup2u       = settingsPtr->parm("HiggsH1:coup2u");
    coup2l       = settingsPtr->parm("HiggsH1:coup2l");
    coup2Z       = settingsPtr->parm("HiggsH1:coup2Z");
    coup2W       = settingsPtr->parm("HiggsH1:coup2W");
    coup2Hchg    = settingsPtr->parm("HiggsH1:coup2Hchg");
  } else if (higgsType == 2) {
    coup2d       = settingsPtr->parm("HiggsH2:coup2d");
    coup2u       = settingsPtr->parm("HiggsH2:coup2u");
    coup2l       = settingsPtr->parm("HiggsH2:coup2l");
    coup2Z       = settingsPtr->parm("HiggsH2:coup2Z");
    coup2W       = settingsPtr->parm("HiggsH2:coup2W");
    coup2Hchg    = settingsPtr->parm("HiggsH2:coup2Hchg");
    coup2H1H1    = settingsPtr->parm("HiggsH2:coup2H1H1");
    coup2A3A3    = settingsPtr->parm("HiggsH2:coup2A3A3");
    coup2H1Z     = settingsPtr->parm("HiggsH2:coup2H1Z");
    coup2A3Z     = settingsPtr->parm("HiggsA3:coup2H2Z");
    coup2A3H1    = settingsPtr->parm("HiggsH2:coup2A3H1");
    coup2HchgW   = settingsPtr->parm("HiggsH2:coup2HchgW");
  } else if (higgsType == 3) {
    coup2d       = settingsPtr->parm("HiggsA3:coup2d");
    coup2u       = settingsPtr->parm("HiggsA3:coup2u");
    coup2l       = settingsPtr->parm("HiggsA3:coup2l");
    coup2Z       = settingsPtr->parm("HiggsA3:coup2Z");
    coup2W       = settingsPtr->parm("HiggsA3:coup2W");
    coup2Hchg    = settingsPtr->parm("HiggsA3:coup2Hchg");
    coup2H1H1    = settingsPtr->parm("HiggsA3:coup2H1H1");
    coup2H1Z     = settingsPtr->parm("HiggsA3:coup2H1Z");
    coup2HchgW   = settingsPtr->parm("HiggsA3:coup2Hchg");
  }

  // Initialization of threshold kinematical factor by stepwise
  // numerical integration of H -> t tbar, Z0 Z0 and W+ W-.
  int psModeT  = (higgsType < 3) ? 3 : 1;
  int psModeWZ = (higgsType < 3) ? 5 : 6;
  mLowT        = max( 2.02 * MASSMINT, 0.5 * mT);
  mStepT       = 0.01 * (3. * mT - mLowT);
  mLowZ        = max( 2.02 * MASSMINWZ, 0.5 * mZ);
  mStepZ       = 0.01 * (3. * mZ - mLowZ);
  mLowW        = max( 2.02 * MASSMINWZ, 0.5 * mW);
  mStepW       = 0.01 * (3. * mW - mLowW);
  for (int i = 0; i <= 100; ++i) {
    kinFacT[i] = numInt2BW( mLowT + i * mStepT,
                 mT, GammaT, MASSMINT,  mT, GammaT, MASSMINT,  psModeT);
    kinFacZ[i] = numInt2BW( mLowZ + i * mStepZ,
                 mZ, GammaZ, MASSMINWZ, mZ, GammaZ, MASSMINWZ, psModeWZ);
    kinFacW[i] = numInt2BW( mLowW + i * mStepW,
                 mW, GammaW, MASSMINWZ, mW, GammaW, MASSMINWZ, psModeWZ);
  }

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceH::calcPreFac(bool) {

  // Common coupling factors.
  alpEM     = couplingsPtr->alphaEM(mHat * mHat);
  alpS      = couplingsPtr->alphaS(mHat * mHat);
  colQ      = 3. * (1. + alpS / M_PI);
  preFac    = (alpEM / (8. * sin2tW)) * pow3(mHat) / pow2(mW);
  if (useNLOWidths) rescColQ = 3. * (1. + rescAlpS * alpS / M_PI) / colQ;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceH::calcWidth(bool) {

  // Widths of decays Higgs -> f + fbar.
  if ( id2Abs == id1Abs && ( (id1Abs > 0 && id1Abs < 7)
    || (id1Abs > 10 && id1Abs < 17) ) ) {
    kinFac = 0.;

    // Check that above threshold (well above for top). Kinematical factor.
    if ( (id1Abs != 6 && mHat > 2. * mf1 + MASSMARGIN)
      || (id1Abs == 6 && mHat > 3. * mT ) ) {
      // A0 behaves like beta, h0 and H0 like beta**3.
      kinFac = (higgsType < 3) ? pow3(ps) : ps;
    }

    // Top near or below threshold: interpolate in table.
    else if (id1Abs == 6 && mHat > mLowT) {
      double xTab = (mHat - mLowT) / mStepT;
      int    iTab = max( 0, min( 99, int(xTab) ) );
      kinFac      = kinFacT[iTab]
                  * pow( kinFacT[iTab + 1] / kinFacT[iTab], xTab - iTab);
    }

    // Coupling from mass and from BSM deviation from SM.
    double coupFac = pow2(particleDataPtr->mRun(id1Abs, mHat) / mHat);
    if (id1Abs < 7 && id1Abs%2 == 1) coupFac *= coup2d * coup2d;
    else if (id1Abs < 7)             coupFac *= coup2u * coup2u;
    else                             coupFac *= coup2l * coup2l;

    // Combine couplings and phase space with colour factor.
    widNow = preFac * coupFac * kinFac;
    if (id1Abs < 7) widNow *= colQ;
  }

  // Widths of decays Higgs -> g + g.
  else if (id1Abs == 21 && id2Abs == 21)
    widNow = preFac * pow2(alpS / M_PI) * eta2gg();

  // Widths of decays Higgs -> gamma + gamma.
  else if (id1Abs == 22 && id2Abs == 22)
    widNow = preFac * pow2(alpEM / M_PI) * 0.5 * eta2gaga();

  // Widths of decays Higgs -> Z0 + gamma0.
  else if (id1Abs == 23 && id2Abs == 22)
    widNow = preFac * pow2(alpEM / M_PI) * pow3(ps) * eta2gaZ();

  // Widths of decays Higgs (h0, H0) -> Z0 + Z0.
  else if (id1Abs == 23 && id2Abs == 23) {
    // If Higgs heavy use on-shell expression, else interpolation in table
    if (mHat > 3. * mZ) kinFac = (1.  - 4. * mr1 + 12. * mr1 * mr1) * ps;
    else if (mHat > mLowZ) {
      double xTab = (mHat - mLowZ) / mStepZ;
      int    iTab = max( 0, min( 99, int(xTab) ) );
      kinFac      = kinFacZ[iTab]
                  * pow( kinFacZ[iTab + 1] / kinFacZ[iTab], xTab - iTab );
    }
    else kinFac   = 0.;
    // Prefactor, normally rescaled to mRes^2 * mHat rather than mHat^3.
    widNow        = 0.25 * preFac * pow2(coup2Z) * kinFac;
    if (!useCubicWidth) widNow *= pow2(mRes / mHat);
  }

  // Widths of decays Higgs (h0, H0) -> W+ + W-.
  else if (id1Abs == 24 && id2Abs == 24) {
    // If Higgs heavy use on-shell expression, else interpolation in table.
    if (mHat > 3. * mW) kinFac = (1.  - 4. * mr1 + 12. * mr1 * mr1) * ps;
    else if (mHat > mLowW) {
      double xTab = (mHat - mLowW) / mStepW;
      int    iTab = max( 0, min( 99, int(xTab) ) );
      kinFac      = kinFacW[iTab]
                  * pow( kinFacW[iTab + 1] / kinFacW[iTab], xTab - iTab);
    }
    else kinFac   = 0.;
    // Prefactor, normally rescaled to mRes^2 * mHat rather than mHat^3.
    widNow        = 0.5 * preFac * pow2(coup2W) * kinFac;
    if (!useCubicWidth) widNow *= pow2(mRes / mHat);
  }

  // Widths of decays Higgs (H0) -> h0 + h0.
  else if (id1Abs == 25 && id2Abs == 25)
    widNow = 0.25 * preFac * pow4(mZ / mHat) * ps * pow2(coup2H1H1);

  // Widths of decays Higgs (H0) -> A0 + A0.
  else if (id1Abs == 36 && id2Abs == 36)
    widNow = 0.5 * preFac * pow4(mZ / mHat) * ps * pow2(coup2A3A3);

  // Widths of decays Higgs (A0) -> h0 + Z0.
  else if (id1Abs == 25 && id2Abs == 23)
    widNow = 0.5 * preFac * pow3(ps) * pow2(coup2H1Z);

  // Widths of decays Higgs (H0) -> A0 + Z0.
  else if (id1Abs == 36 && id2Abs == 23)
    widNow = 0.5 * preFac * pow3(ps) * pow2(coup2A3Z);

  // Widths of decays Higgs (H0) -> A0 + h0.
  else if (id1Abs == 36 && id2Abs == 25)
    widNow = 0.25 * preFac * pow4(mZ / mHat) * ps * pow2(coup2A3H1);

  // Widths of decays Higgs -> H+- + W-+.
  else if (id1Abs == 37 && id2Abs == 24)
    widNow = 0.5 * preFac * pow3(ps) * pow2(coup2HchgW);

  // NLO multiplicative factors for SM h0 (125 GeV) based on LHCXSWG
  // recommendations.
  if (useNLOWidths) {
    if      (id1Abs == 21 && id2Abs == 21) widNow *= 1.47 * pow2(rescAlpS);
    else if (id1Abs == 22 && id2Abs == 22) widNow *= 0.88;
    else if (id1Abs == 22 && id2Abs == 23) widNow *= 0.95;
    else if (id1Abs == 23 && id2Abs == 23) widNow *= 1.10;
    else if (id1Abs == 24 && id2Abs == 24) widNow *= 1.09;
    else if (id1Abs ==  5 && id2Abs ==  5) widNow *= 1.07  * rescColQ;
    else if (id1Abs ==  4 && id2Abs ==  4) widNow *= 0.937 * rescColQ;
    else if (id1Abs == 13 && id2Abs == 13) widNow *= 0.974;
    else if (id1Abs == 15 && id2Abs == 15) widNow *= 0.992;
  }

}

//--------------------------------------------------------------------------

// Sum up quark loop contributions in Higgs -> g + g.
// Note: running quark masses are used, unlike Pythia6 (not negligible shift).

double ResonanceH::eta2gg() {

  // Initial values.
  complex eta = complex(0., 0.);
  double  mLoop, epsilon, root, rootLog;
  complex phi, etaNow;

  // Loop over s, c, b, t quark flavours.
  for (int idNow = 3; idNow < 7; ++idNow) {
    mLoop   = (useRunLoopMass) ? particleDataPtr->mRun(idNow, mHat)
                               : particleDataPtr->m0(idNow);
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
    if (higgsType < 3) etaNow = -0.5 * epsilon
      * (complex(1., 0.) + (1. - epsilon) * phi);
    else etaNow = -0.5 * epsilon * phi;
    if (idNow%2 == 1) etaNow *= coup2d;
    else              etaNow *= coup2u;

    // Sum up contribution and return square of absolute value.
    eta += etaNow;
  }
  return (pow2(eta.real()) + pow2(eta.imag()));

}

//--------------------------------------------------------------------------

// Sum up quark, lepton, W+- and (for BSM) H+- loop contributions
// in Higgs -> gamma + gamma.

double ResonanceH::eta2gaga() {

  // Initial values.
  complex eta = complex(0., 0.);
  int     idNow;
  double  ef, mLoop, epsilon, root, rootLog;
  complex phi, etaNow;

  // Loop over s, c, b, t, mu, tau, W+-, H+- flavours.
  for (int idLoop = 0; idLoop < 8; ++idLoop) {
    if      (idLoop < 4) idNow = idLoop + 3;
    else if (idLoop < 6) idNow = 2 * idLoop + 5;
    else if (idLoop < 7) idNow = 24;
    else                 idNow = 37;
    if (idNow == 37 && higgsType == 0) continue;

    // Charge and loop integral parameter.
    ef      = (idNow < 20) ? couplingsPtr->ef(idNow) : 1.;
    mLoop   = (useRunLoopMass) ? particleDataPtr->mRun(idNow, mHat)
                               : particleDataPtr->m0(idNow);
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

    // Expressions for quarks and leptons that depend on Higgs type.
    if (idNow < 17) {
      if (higgsType < 3) etaNow = -0.5 * epsilon
        * (complex(1., 0.) + (1. - epsilon) * phi);
      else etaNow = -0.5 * epsilon * phi;
      if (idNow < 7 && idNow%2 == 1) etaNow *= 3. * pow2(ef) * coup2d;
      else if (idNow < 7 )           etaNow *= 3. * pow2(ef) * coup2u;
      else                           etaNow *=      pow2(ef) * coup2l;
    }

    // Expression for W+-.
    else if (idNow == 24) etaNow = (complex(0.5 + 0.75 * epsilon, 0.)
      + 0.75 * epsilon * (2. - epsilon) * phi) * coup2W;

    // Expression for H+-.
   else etaNow = (complex(epsilon, 0.) - epsilon * epsilon * phi)
     * pow2(mW / mHchg) * coup2Hchg;

    // Sum up contribution and return square of absolute value.
    eta       += etaNow;
  }
  return (pow2(eta.real()) + pow2(eta.imag()));

}

//--------------------------------------------------------------------------

// Sum up quark, lepton, W+- and (for BSM) H+- loop contributions
// in Higgs -> gamma + Z0.

double ResonanceH::eta2gaZ() {

  // Initial values.
  complex eta = complex(0., 0.);
  int     idNow;
  double  ef, vf, mLoop, epsilon, epsPrime, root, rootLog, asinEps;
  complex phi, psi, phiPrime, psiPrime, fXY, f1, etaNow;

  // Loop over s, c, b, t, mu , tau, W+-, H+- flavours.
  for (int idLoop = 0; idLoop < 7; ++idLoop) {
    if      (idLoop < 4) idNow = idLoop + 3;
    else if (idLoop < 6) idNow = 2 * idLoop + 5;
    else if (idLoop < 7) idNow = 24;
    else                 idNow = 37;

    // Electroweak charges and loop integral parameters.
    ef        = (idNow < 20) ? couplingsPtr->ef(idNow) : 1.;
    vf        = (idNow < 20) ? couplingsPtr->vf(idNow) : 0.;
    mLoop     = (useRunLoopMass) ? particleDataPtr->mRun(idNow, mHat)
                                 : particleDataPtr->m0(idNow);
    epsilon   = pow2(2. * mLoop / mHat);
    epsPrime  = pow2(2. * mLoop / mZ);

    // Value of loop integral for epsilon = 4 m^2 / sHat.
    if (epsilon <= 1.) {
      root    = sqrt(1. - epsilon);
      rootLog = (epsilon < 1e-4) ? log(4. / epsilon - 2.)
                : log( (1. + root) / (1. - root) );
      phi     = complex( -0.25 * (pow2(rootLog) - pow2(M_PI)),
                0.5 * M_PI * rootLog );
      psi     = 0.5 * root * complex( rootLog, -M_PI);
    } else {
      asinEps = asin(1. / sqrt(epsilon));
      phi     = complex( pow2(asinEps), 0.);
      psi     = complex( sqrt(epsilon - 1.) * asinEps, 0.);
    }

    // Value of loop integral for epsilonPrime = 4 m^2 / m_Z^2.
    if (epsPrime <= 1.) {
      root     = sqrt(1. - epsPrime);
      rootLog  = (epsPrime < 1e-4) ? log(4. / epsPrime - 2.)
                 : log( (1. + root) / (1. - root) );
      phiPrime = complex( -0.25 * (pow2(rootLog) - pow2(M_PI)),
                          0.5 * M_PI * rootLog );
      psiPrime = 0.5 * root * complex( rootLog, -M_PI);
    } else {
      asinEps  = asin(1. / sqrt(epsPrime));
      phiPrime = complex( pow2(asinEps), 0.);
      psiPrime = complex( sqrt(epsPrime - 1.) * asinEps, 0.);
    }

    // Combine the two loop integrals.
    fXY = (epsilon * epsPrime / (8. * pow2(epsilon - epsPrime)))
      * ( complex(epsilon - epsPrime, 0)
      + epsilon * epsPrime * (phi - phiPrime)
      + 2. * epsilon * (psi - psiPrime) );
    f1 = - (epsilon * epsPrime / (2. * (epsilon - epsPrime)))
      * (phi - phiPrime);

    // Expressions for quarks and leptons that depend on Higgs type.
    if (idNow < 17) {
      etaNow = (higgsType < 3) ? -fXY + 0.25 * f1 : 0.25 * f1;
      if (idNow < 7 && idNow%2 == 1) etaNow *= 3. * ef * vf * coup2d;
      else if (idNow < 7)         etaNow *= 3. * ef * vf * coup2u;
      else                     etaNow *=      ef * vf * coup2l;

    // Expression for W+-.
    } else if (idNow == 24) {
      double coef1  = 3. - sin2tW / cos2tW;
      double coefXY = (1. + 2. / epsilon) * sin2tW / cos2tW
        - (5. + 2. / epsilon);
      etaNow = -cos2tW * (coef1 * f1 + coefXY * fXY) * coup2W;

    // Expression for H+-.
    } else etaNow = (1. - 2. * sin2tW) * fXY * pow2(mW / mHchg)
      * coup2Hchg;

    // Sum up contribution and return square of absolute value.
    eta += etaNow;
  }
  return ( (pow2(eta.real()) + pow2(eta.imag())) / (sin2tW * cos2tW) );

}

//==========================================================================

// The ResonanceHchg class.
// Derived class for H+- properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceHchg::initConstants() {

  // Locally stored properties and couplings.
  useCubicWidth = settingsPtr->flag("Higgs:cubicWidth");
  thetaWRat     = 1. / (8. * couplingsPtr->sin2thetaW());
  mW            = particleDataPtr->m0(24);
  tanBeta       = settingsPtr->parm("HiggsHchg:tanBeta");
  tan2Beta      = tanBeta * tanBeta;
  coup2H1W      = settingsPtr->parm("HiggsHchg:coup2H1W");

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceHchg::calcPreFac(bool) {

  // Common coupling factors.
  alpEM         = couplingsPtr->alphaEM(mHat * mHat);
  alpS          = couplingsPtr->alphaS(mHat * mHat);
  colQ          = 3. * (1. + alpS / M_PI);
  preFac        = alpEM * thetaWRat * pow3(mHat) / pow2(mW);

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceHchg::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // H+- decay to fermions involves running masses.
  if (id1Abs < 17 && (id1Abs < 7 || id1Abs > 10)) {
    double mRun1   = particleDataPtr->mRun(id1Abs, mHat);
    double mRun2   = particleDataPtr->mRun(id2Abs, mHat);
    double mrRunDn = pow2(mRun1 / mHat);
    double mrRunUp = pow2(mRun2 / mHat);
    if (id1Abs%2 == 0) swap( mrRunDn, mrRunUp);

    // Width to fermions: couplings, kinematics, colour factor.
    widNow = preFac * max( 0., (mrRunDn * tan2Beta + mrRunUp / tan2Beta)
           * (1. - mrRunDn - mrRunUp) - 4. *mrRunDn * mrRunUp ) * ps;
    if (id1Abs < 7) widNow *= colQ;
  }

  // H+- decay to h0 + W+-.
  else if (id1Abs == 25 && id2Abs == 24)
    widNow    = 0.5 * preFac * pow3(ps) * pow2(coup2H1W);

}

//==========================================================================

// The ResonanceZprime class.
// Derived class for gamma*/Z0/Z'^0 properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceZprime::initConstants() {

  // Locally stored properties and couplings.
  gmZmode     = settingsPtr->mode("Zprime:gmZmode");
  sin2tW      = couplingsPtr->sin2thetaW();
  cos2tW      = 1. - sin2tW;
  thetaWRat   = 1. / (16. * sin2tW * cos2tW);

  // Properties of Z resonance.
  mZ          = particleDataPtr->m0(23);
  GammaZ      = particleDataPtr->mWidth(23);
  m2Z         = mZ*mZ;
  GamMRatZ    = GammaZ / mZ;

  // Ensure that arrays initially empty.
  for (int i = 0; i < 20; ++i) afZp[i] = 0.;
  for (int i = 0; i < 20; ++i) vfZp[i] = 0.;

  // Store first-generation axial and vector couplings.
  afZp[1]     = settingsPtr->parm("Zprime:ad");
  afZp[2]     = settingsPtr->parm("Zprime:au");
  afZp[11]    = settingsPtr->parm("Zprime:ae");
  afZp[12]    = settingsPtr->parm("Zprime:anue");
  vfZp[1]     = settingsPtr->parm("Zprime:vd");
  vfZp[2]     = settingsPtr->parm("Zprime:vu");
  vfZp[11]    = settingsPtr->parm("Zprime:ve");
  vfZp[12]    = settingsPtr->parm("Zprime:vnue");

  // Determine if the 4th generation should be included
  bool coupZp2gen4    = settingsPtr->flag("Zprime:coup2gen4");
  maxZpGen = (coupZp2gen4) ? 8 : 6;

  // Second and third generation could be carbon copy of this...
  if (settingsPtr->flag("Zprime:universality")) {
    for (int i = 3; i <= maxZpGen; ++i) {
      afZp[i]    = afZp[i-2];
      vfZp[i]    = vfZp[i-2];
      afZp[i+10] = afZp[i+8];
      vfZp[i+10] = vfZp[i+8];
    }

  // ... or could have different couplings.
  } else {
    afZp[3]   = settingsPtr->parm("Zprime:as");
    afZp[4]   = settingsPtr->parm("Zprime:ac");
    afZp[5]   = settingsPtr->parm("Zprime:ab");
    afZp[6]   = settingsPtr->parm("Zprime:at");
    afZp[13]  = settingsPtr->parm("Zprime:amu");
    afZp[14]  = settingsPtr->parm("Zprime:anumu");
    afZp[15]  = settingsPtr->parm("Zprime:atau");
    afZp[16]  = settingsPtr->parm("Zprime:anutau");
    vfZp[3]   = settingsPtr->parm("Zprime:vs");
    vfZp[4]   = settingsPtr->parm("Zprime:vc");
    vfZp[5]   = settingsPtr->parm("Zprime:vb");
    vfZp[6]   = settingsPtr->parm("Zprime:vt");
    vfZp[13]  = settingsPtr->parm("Zprime:vmu");
    vfZp[14]  = settingsPtr->parm("Zprime:vnumu");
    vfZp[15]  = settingsPtr->parm("Zprime:vtau");
    vfZp[16]  = settingsPtr->parm("Zprime:vnutau");
    if( coupZp2gen4 ) {
      afZp[7]   = settingsPtr->parm("Zprime:abPrime");
      afZp[8]   = settingsPtr->parm("Zprime:atPrime");
      vfZp[7]   = settingsPtr->parm("Zprime:vbPrime");
      vfZp[8]   = settingsPtr->parm("Zprime:vtPrime");
      afZp[17]  = settingsPtr->parm("Zprime:atauPrime");
      afZp[18]  = settingsPtr->parm("Zprime:anutauPrime");
      vfZp[17]  = settingsPtr->parm("Zprime:vtauPrime");
      vfZp[18]  = settingsPtr->parm("Zprime:vnutauPrime");
    }
  }

  // Coupling for Z' -> W+ W-.
  coupZpWW    = settingsPtr->parm("Zprime:coup2WW");

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceZprime::calcPreFac(bool calledFromInit) {

  // Common coupling factors.
  alpEM       = couplingsPtr->alphaEM(mHat * mHat);
  alpS        = couplingsPtr->alphaS(mHat * mHat);
  colQ        = 3. * (1. + alpS / M_PI);
  preFac      = alpEM * thetaWRat * mHat / 3.;

  // When call for incoming flavour need to consider gamma*/Z0 mix.
  if (!calledFromInit) {

    // Couplings when an incoming fermion is specified; elso only pure Z'0.
    ei2       = 0.;
    eivi      = 0.;
    vai2      = 0.;
    eivpi     = 0.;
    vaivapi   = 0.,
    vapi2     = 1.;
    int idInFlavAbs = abs(idInFlav);
    if ( (idInFlavAbs >  0 && idInFlavAbs <= maxZpGen)
      || (idInFlavAbs > 10 && idInFlavAbs <= maxZpGen + 10) ) {
      double ei  = couplingsPtr->ef(idInFlavAbs);
      double ai  = couplingsPtr->af(idInFlavAbs);
      double vi  = couplingsPtr->vf(idInFlavAbs);
      double api = afZp[idInFlavAbs];
      double vpi = vfZp[idInFlavAbs];
      ei2     = ei * ei;
      eivi    = ei * vi;
      vai2    = vi * vi + ai * ai;
      eivpi   = ei * vpi;
      vaivapi = vi * vpi + ai * api;;
      vapi2   = vpi * vpi + api * api;
    }

    // Calculate prefactors for gamma/interference/Z0 terms.
    double sH     = mHat * mHat;
    double propZ  = sH / ( pow2(sH - m2Z) + pow2(sH * GamMRatZ) );
    double propZp = sH / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
    gamNorm   = ei2;
    gamZNorm  = 2. * eivi * thetaWRat * (sH - m2Z) * propZ;
    ZNorm     = vai2 * pow2(thetaWRat) * sH * propZ;
    gamZpNorm = 2. * eivpi * thetaWRat * (sH - m2Res) * propZp;
    ZZpNorm   = 2. * vaivapi * pow2(thetaWRat) * ((sH - m2Res) * (sH - m2Z)
              + sH * GamMRat * sH * GamMRatZ) * propZ * propZp;
    ZpNorm    = vapi2 * pow2(thetaWRat) * sH * propZp;

    // Optionally only keep some of gamma*, Z0 and Z' terms.
    if (gmZmode == 1) {gamZNorm = 0; ZNorm = 0.; gamZpNorm = 0.;
      ZZpNorm = 0.; ZpNorm = 0.;}
    if (gmZmode == 2) {gamNorm = 0.; gamZNorm = 0.; gamZpNorm = 0.;
      ZZpNorm = 0.; ZpNorm = 0.;}
    if (gmZmode == 3) {gamNorm = 0.; gamZNorm = 0.; ZNorm = 0.;
      gamZpNorm = 0.; ZZpNorm = 0.;}
    if (gmZmode == 4) {gamZpNorm = 0.; ZZpNorm = 0.; ZpNorm = 0.;}
    if (gmZmode == 5) {gamZNorm = 0.; ZNorm = 0.; ZZpNorm = 0.;}
    if (gmZmode == 6) {gamNorm = 0.; gamZNorm = 0.; gamZpNorm = 0.;}
  }

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceZprime::calcWidth(bool calledFromInit) {

  // Check that above threshold.
  if (ps == 0.) return;

  // At initialization only the pure Z'0 should be considered.
  if (calledFromInit) {

    // Contributions from three (4?) fermion generations.
    if ( id1Abs <= maxZpGen || (id1Abs > 10 && id1Abs <= maxZpGen + 10) ) {
      double apf = afZp[id1Abs];
      double vpf = vfZp[id1Abs];
      widNow = preFac * ps * (vpf*vpf * (1. + 2. * mr1)
             + apf*apf * ps*ps);
      if (id1Abs < 9) widNow *= colQ;

    // Contribution from Z'0 -> W^+ W^-.
    } else if (id1Abs == 24) {
      widNow = preFac * pow2(coupZpWW * cos2tW) * pow3(ps)
        * (1. + mr1*mr1 + mr2*mr2 + 10. * (mr1 + mr2 + mr1 * mr2));
    }
  }

  // When call for incoming flavour need to consider full mix.
  else {

    // Contributions from three (4?) fermion generations.
    if ( id1Abs <= maxZpGen || (id1Abs > 10 && id1Abs <= maxZpGen + 10) ) {

      // Couplings of gamma^*/Z^0/Z'^0  to final flavour
      double ef  = couplingsPtr->ef(id1Abs);
      double af  = couplingsPtr->af(id1Abs);
      double vf  = couplingsPtr->vf(id1Abs);
      double apf = afZp[id1Abs];
      double vpf = vfZp[id1Abs];

      // Combine couplings with kinematical factors.
      double kinFacA  = pow3(ps);
      double kinFacV  = ps * (1. + 2. * mr1);
      double ef2      = ef * ef * kinFacV;
      double efvf     = ef * vf * kinFacV;
      double vaf2     = vf * vf * kinFacV + af * af * kinFacA;
      double efvpf    = ef * vpf * kinFacV;
      double vafvapf  = vf * vpf * kinFacV + af * apf * kinFacA;
      double vapf2    = vpf * vpf * kinFacV + apf * apf * kinFacA;

      // Relative outwidths: combine instate, propagator and outstate.
      widNow = gamNorm * ef2 + gamZNorm * efvf + ZNorm * vaf2
             + gamZpNorm * efvpf + ZZpNorm * vafvapf + ZpNorm * vapf2;
      if (id1Abs < 9) widNow *= colQ;

    // Contribution from Z'0 -> W^+ W^-.
    } else if (id1Abs == 24) {
      widNow = ZpNorm * pow2(coupZpWW * cos2tW) * pow3(ps)
        * (1. + mr1*mr1 + mr2*mr2 + 10. * (mr1 + mr2 + mr1 * mr2));
    }
  }

}

//==========================================================================

// The ResonanceWprime class.
// Derived class for W'+- properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceWprime::initConstants() {

  // Locally stored properties and couplings.
  thetaWRat = 1. / (12. * couplingsPtr->sin2thetaW());
  cos2tW    = couplingsPtr->cos2thetaW();

  // Axial and vector couplings of fermions.
  aqWp      = settingsPtr->parm("Wprime:aq");
  vqWp      = settingsPtr->parm("Wprime:vq");
  alWp      = settingsPtr->parm("Wprime:al");
  vlWp      = settingsPtr->parm("Wprime:vl");

  // Coupling for W' -> W Z.
  coupWpWZ    = settingsPtr->parm("Wprime:coup2WZ");

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceWprime::calcPreFac(bool) {

  // Common coupling factors.
  alpEM     = couplingsPtr->alphaEM(mHat * mHat);
  alpS      = couplingsPtr->alphaS(mHat * mHat);
  colQ      = 3. * (1. + alpS / M_PI);
  preFac    = alpEM * thetaWRat * mHat;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceWprime::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // Decay to quarks involves colour factor and CKM matrix.
  if (id1Abs > 0 && id1Abs < 9) widNow
    = preFac * ps * 0.5 * ((vqWp*vqWp + aqWp * aqWp)
    + 6. * (vqWp*vqWp - aqWp * aqWp) * sqrt(mr1 *mr2))
    * (1. - 0.5 * (mr1 + mr2) - 0.5 * pow2(mr1 - mr2))
    * colQ * couplingsPtr->V2CKMid(id1Abs, id2Abs);

  // Decay to leptons simpler.
  else if (id1Abs > 10 && id1Abs < 19) widNow
    = preFac * ps * 0.5 * ((vlWp*vqWp + alWp * aqWp)
    + 6. * (vlWp*vqWp - alWp * aqWp) * sqrt(mr1 *mr2))
    * (1. - 0.5 * (mr1 + mr2) - 0.5 * pow2(mr1 - mr2));

  // Decay to W^+- Z^0.
  else if (id1Abs == 24 && id2Abs == 23) widNow
    = preFac * 0.25 * pow2(coupWpWZ) * cos2tW * (mr1 / mr2) * pow3(ps)
    * (1. + mr1*mr1 + mr2*mr2 + 10. * (mr1 + mr2 + mr1 * mr2));

}

//==========================================================================

// The ResonanceRhorizontal class.
// Derived class for R^0 (horizontal gauge boson) properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceRhorizontal::initConstants() {

  // Locally stored properties and couplings.
  thetaWRat = 1. / (12. * couplingsPtr->sin2thetaW());

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceRhorizontal::calcPreFac(bool) {

  // Common coupling factors.
  alpEM     = couplingsPtr->alphaEM(mHat * mHat);
  alpS      = couplingsPtr->alphaS(mHat * mHat);
  colQ      = 3. * (1. + alpS / M_PI);
  preFac    = alpEM * thetaWRat * mHat;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceRhorizontal::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // R -> f fbar. Colour factor for quarks.
  widNow    = preFac * ps * (2. - mr1 - mr2 - pow2(mr1 - mr2));
  if (id1Abs < 9) widNow *= colQ;

}

//==========================================================================

// The ResonanceExcited class.
// Derived class for excited-fermion properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceExcited::initConstants() {

  // Locally stored properties and couplings.
  Lambda        = settingsPtr->parm("ExcitedFermion:Lambda");
  coupF         = settingsPtr->parm("ExcitedFermion:coupF");
  coupFprime    = settingsPtr->parm("ExcitedFermion:coupFprime");
  coupFcol      = settingsPtr->parm("ExcitedFermion:coupFcol");
  contactDec    = settingsPtr->parm("ExcitedFermion:contactDec");
  sin2tW        = couplingsPtr->sin2thetaW();
  cos2tW        = 1. - sin2tW;

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceExcited::calcPreFac(bool) {

  // Common coupling factors.
  alpEM         = couplingsPtr->alphaEM(mHat * mHat);
  alpS          = couplingsPtr->alphaS(mHat * mHat);
  preFac        = pow3(mHat) / pow2(Lambda);

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceExcited::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // f^* -> f g.
  if (id1Abs == 21) widNow = preFac * alpS * pow2(coupFcol) / 3.;

  // f^* -> f gamma.
  else if (id1Abs == 22) {
    double chgI3 = (id2Abs%2 == 0) ? 0.5 : -0.5;
    double chgY  = (id2Abs < 9) ? 1. / 6. : -0.5;
    double chg   = chgI3 * coupF + chgY * coupFprime;
    widNow       = preFac * alpEM * pow2(chg) / 4.;
  }

  // f^* -> f Z^0.
  else if (id1Abs == 23) {
    double chgI3 = (id2Abs%2 == 0) ? 0.5 : -0.5;
    double chgY  = (id2Abs < 9) ? 1. / 6. : -0.5;
    double chg   = chgI3 * cos2tW * coupF - chgY * sin2tW * coupFprime;
    widNow       = preFac * (alpEM * pow2(chg) / (8. * sin2tW * cos2tW))
                 * ps*ps * (2. + mr1);
  }

  // f^* -> f' W^+-.
  else if (id1Abs == 24) widNow  = preFac * (alpEM * pow2(coupF)
                 / (16. * sin2tW)) * ps*ps * (2. + mr1);

  // f^* -> f f' fbar' contact interaction decays (code by Olga Igonkina).
  else {
    if (id1Abs < 17 && id2Abs < 17 && id3Abs > 0 && id3Abs < 17 ) {
      widNow = preFac * pow2(contactDec * mHat) / (pow2(Lambda) * 96. * M_PI);
      if (id3Abs < 10) widNow *= 3.;
      if (id1Abs == id2Abs && id1Abs == id3Abs) {
        if (idRes - 4000000 < 10) widNow *= 4./3.;
        else                      widNow *= 2.;
      }
    }
  }

}

//==========================================================================

// The ResonanceGraviton class.
// Derived class for excited Graviton properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceGraviton::initConstants() {

  // SMinBulk = off/on, use universal coupling (kappaMG)
  // or individual (Gxx) between graviton and SM particles.
  eDsmbulk   = settingsPtr->flag("ExtraDimensionsG*:SMinBulk");
  eDvlvl = false;
  if (eDsmbulk) eDvlvl = settingsPtr->flag("ExtraDimensionsG*:VLVL");
  kappaMG    = settingsPtr->parm("ExtraDimensionsG*:kappaMG");
  for (int i = 0; i < 27; ++i) eDcoupling[i] = 0.;
  double tmp_coup = settingsPtr->parm("ExtraDimensionsG*:Gqq");
  for (int i = 1; i <= 4; ++i)  eDcoupling[i] = tmp_coup;
  eDcoupling[5] = settingsPtr->parm("ExtraDimensionsG*:Gbb");
  eDcoupling[6] = settingsPtr->parm("ExtraDimensionsG*:Gtt");
  tmp_coup = settingsPtr->parm("ExtraDimensionsG*:Gll");
  for (int i = 11; i <= 16; ++i) eDcoupling[i] = tmp_coup;
  eDcoupling[21] = settingsPtr->parm("ExtraDimensionsG*:Ggg");
  eDcoupling[22] = settingsPtr->parm("ExtraDimensionsG*:Ggmgm");
  eDcoupling[23] = settingsPtr->parm("ExtraDimensionsG*:GZZ");
  eDcoupling[24] = settingsPtr->parm("ExtraDimensionsG*:GWW");
  eDcoupling[25] = settingsPtr->parm("ExtraDimensionsG*:Ghh");

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceGraviton::calcPreFac(bool) {

  // Common coupling factors.
  alpS          = couplingsPtr->alphaS(mHat * mHat);
  colQ          = 3. * (1. + alpS / M_PI);
  preFac        = mHat / M_PI;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceGraviton::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // Widths to fermion pairs.
  if (id1Abs < 19) {
     widNow  = preFac * pow3(ps) * (1. + 8. * mr1 / 3.) / 320.;
     if (id1Abs < 9) widNow *= colQ;

  // Widths to gluon and photon pair.
  } else if (id1Abs == 21) {
    widNow = preFac / 20.;
  } else if (id1Abs == 22) {
    widNow = preFac / 160.;

  // Widths to Z0 Z0 and W+ W- pair.
  } else if (id1Abs == 23 || id1Abs == 24) {
    // Longitudinal W/Z only.
    if (eDvlvl) {
      widNow = preFac * pow(ps,5) / 480.;
    // Transverse W/Z contributions as well.
    } else {
      widNow  = preFac * ps * (13. / 12. + 14. * mr1 / 3. + 4. * mr1 * mr1)
              / 80.;
    }
    if (id1Abs == 23) widNow *= 0.5;

  // Widths to h h pair.
  } else if (id1Abs == 25) {
    widNow = preFac * pow(ps,5) / 960.;
  }

  // RS graviton coupling
  if (eDsmbulk) widNow *= 2. * pow2(eDcoupling[min( id1Abs, 26)] * mHat);
  else          widNow *= pow2(kappaMG * mHat / mRes);

}

//==========================================================================

// The ResonanceKKgluon class.
// Derived class for excited kk-gluon properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceKKgluon::initConstants() {

  // KK-gluon gv/ga couplings and interference.
  for (int i = 0; i < 10; ++i) { eDgv[i] = 0.; eDga[i] = 0.; }
  double tmp_gL = settingsPtr->parm("ExtraDimensionsG*:KKgqL");
  double tmp_gR = settingsPtr->parm("ExtraDimensionsG*:KKgqR");
  for (int i = 1; i <= 4; ++i) {
    eDgv[i] = 0.5 * (tmp_gL + tmp_gR);
    eDga[i] = 0.5 * (tmp_gL - tmp_gR);
  }
  tmp_gL = settingsPtr->parm("ExtraDimensionsG*:KKgbL");
  tmp_gR = settingsPtr->parm("ExtraDimensionsG*:KKgbR");
  eDgv[5] = 0.5 * (tmp_gL + tmp_gR); eDga[5] = 0.5 * (tmp_gL - tmp_gR);
  tmp_gL = settingsPtr->parm("ExtraDimensionsG*:KKgtL");
  tmp_gR = settingsPtr->parm("ExtraDimensionsG*:KKgtR");
  eDgv[6] = 0.5 * (tmp_gL + tmp_gR); eDga[6] = 0.5 * (tmp_gL - tmp_gR);
  interfMode    = settingsPtr->mode("ExtraDimensionsG*:KKintMode");

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceKKgluon::calcPreFac(bool calledFromInit) {

  // Common coupling factors.
  alpS          = couplingsPtr->alphaS(mHat * mHat);
  preFac        = alpS * mHat / 6;

  // When call for incoming flavour need to consider g*/gKK mix.
  if (!calledFromInit) {
    // Calculate prefactors for g/interference/gKK terms.
    int idInFlavAbs = abs(idInFlav);
    double sH = mHat * mHat;
    normSM   = 1;
    normInt  = 2. * eDgv[min(idInFlavAbs, 9)] * sH * (sH - m2Res)
              / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
    normKK   = ( pow2(eDgv[min(idInFlavAbs, 9)])
               + pow2(eDga[min(idInFlavAbs, 9)]) ) * sH * sH
              / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

    // Optionally only keep g* or gKK term.
    if (interfMode == 1) {normInt = 0.; normKK = 0.;}
    if (interfMode == 2) {normSM = 0.; normInt = 0.; normKK = 1.;}
  }

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceKKgluon::calcWidth(bool calledFromInit) {

  // Check that above threshold.
  if (ps == 0.) return;

  // Widths to quark pairs.
  if (id1Abs > 9) return;

  if (calledFromInit) {
    widNow = preFac * ps * (pow2(eDgv[min(id1Abs, 9)]) * (1. + 2.*mr1)
                         +  pow2(eDga[min(id1Abs, 9)]) * (1. - 4.*mr1) );
  } else {
    // Relative outwidths: combine instate, propagator and outstate.
    widNow = normSM  * ps * (1. + 2. * mr1)
           + normInt * ps * eDgv[min(id1Abs, 9)] * (1. + 2. * mr1)
           + normKK  * ps * (pow2(eDgv[min(id1Abs, 9)]) * (1. + 2.*mr1)
                          +  pow2(eDga[min(id1Abs, 9)]) * (1. - 4.*mr1) );
    widNow *= preFac;
  }

}

//==========================================================================

// The ResonanceLeptoquark class.
// Derived class for leptoquark properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceLeptoquark::initConstants() {

  // Locally stored properties and couplings.
  kCoup      = settingsPtr->parm("LeptoQuark:kCoup");

  // Check that flavour info in decay channel is correctly set.
  int id1Now = particlePtr->channel(0).product(0);
  int id2Now = particlePtr->channel(0).product(1);
  if (id1Now < 1 || id1Now > 6) {
    infoPtr->errorMsg("Error in ResonanceLeptoquark::init:"
      " unallowed input quark flavour reset to u");
    id1Now   = 2;
    particlePtr->channel(0).product(0, id1Now);
  }
  if (abs(id2Now) < 11 || abs(id2Now) > 16) {
    infoPtr->errorMsg("Error in ResonanceLeptoquark::init:"
      " unallowed input lepton flavour reset to e-");
    id2Now   = 11;
    particlePtr->channel(0).product(1, id2Now);
  }

  // Set/overwrite charge and name of particle.
  bool changed  = particlePtr->hasChanged();
  int chargeLQ  = particleDataPtr->chargeType(id1Now)
                + particleDataPtr->chargeType(id2Now);
  particlePtr->setChargeType(chargeLQ);
  string nameLQ = "LQ_" + particleDataPtr->name(id1Now) + ","
                + particleDataPtr->name(id2Now);
  particlePtr->setNames(nameLQ, nameLQ + "bar");
  if (!changed) particlePtr->setHasChanged(false);

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceLeptoquark::calcPreFac(bool) {

  // Common coupling factors.
  alpEM         = couplingsPtr->alphaEM(mHat * mHat);
  preFac        = 0.25 * alpEM * kCoup * mHat;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceLeptoquark::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // Width into lepton plus quark.
  if (id1Abs > 10 && id1Abs < 17 && id2Abs < 7) widNow = preFac * pow3(ps);

}

//==========================================================================

// The ResonanceNuRight class.
// Derived class for righthanded Majorana neutrino properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceNuRight::initConstants() {

  // Locally stored properties and couplings: righthanded W mass.
  thetaWRat = 1. / (768. * M_PI * pow2(couplingsPtr->sin2thetaW()));
  mWR       = particleDataPtr->m0(9900024);

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceNuRight::calcPreFac(bool) {

  // Common coupling factors.
  alpEM     = couplingsPtr->alphaEM(mHat * mHat);
  alpS      = couplingsPtr->alphaS(mHat * mHat);
  colQ      = 3. * (1. + alpS / M_PI);
  preFac    = pow2(alpEM) * thetaWRat * pow5(mHat) / pow4(max(mHat, mWR));

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceNuRight::calcWidth(bool) {

  // Check that above threshold.
  if (mHat < mf1 + mf2 + mf3 + MASSMARGIN) return;

  // Coupling part of widths to l- q qbar', l- l'+ nu_lR' and c.c.
  widNow    = (id2Abs < 9 && id3Abs < 9)
            ? preFac * colQ * couplingsPtr->V2CKMid(id2, id3) : preFac;

  // Phase space corrections in decay. Must have y < 1.
  double x  = (mf1 + mf2 + mf3) / mHat;
  double x2 = x * x;
  double fx = 1. - 8. * x2 + 8. * pow3(x2) - pow4(x2)
            - 24. * pow2(x2) * log(x);
  double y  = min( 0.999, pow2(mHat / mWR) );
  double fy = ( 12. * (1. - y) * log(1. - y) + 12. * y - 6. * y*y
            - 2.* pow3(y) ) / pow4(y);
  widNow   *= fx * fy;

}

//==========================================================================

// The ResonanceZRight class.
// Derived class for Z_R^0 properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceZRight::initConstants() {

  // Locally stored properties and couplings: righthanded W mass.
  sin2tW    = couplingsPtr->sin2thetaW();
  thetaWRat = 1. / (48. * sin2tW  * (1. - sin2tW) * (1. - 2. * sin2tW) );

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceZRight::calcPreFac(bool) {

  // Common coupling factors.
  alpEM     = couplingsPtr->alphaEM(mHat * mHat);
  alpS      = couplingsPtr->alphaS(mHat * mHat);
  colQ      = 3. * (1. + alpS / M_PI);
  preFac    = alpEM * thetaWRat * mHat;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceZRight::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // Couplings to q qbar and l+ l-.
  double vf = 0.;
  double af = 0.;
  double symMaj = 1.;
  if (id1Abs < 9 && id1Abs%2 == 1) {
    af      = -1. + 2. * sin2tW;
    vf      = -1. + 4. * sin2tW / 3.;
  } else if (id1Abs < 9) {
    af      = 1. - 2. * sin2tW;
    vf      = 1. - 8. * sin2tW / 3.;
  } else if (id1Abs < 19 && id1Abs%2 == 1) {
    af      = -1. + 2. * sin2tW;
    vf      = -1. + 4. * sin2tW;

  // Couplings to nu_L nu_Lbar and nu_R nu_Rbar, both assumed Majoranas.
  } else if (id1Abs < 19) {
    af      = -2. * sin2tW;
    vf      = 0.;
    symMaj  = 0.5;
  } else {
    af      = 2. * (1. - sin2tW);
    vf      = 0.;
    symMaj  = 0.5;
  }

  // Width expression, including phase space and colour factor.
  widNow    = preFac * (vf*vf * (1. + 2. * mr1) + af*af * ps*ps) * ps
            * symMaj;
  if (id1Abs < 9) widNow *= colQ;

}

//==========================================================================

// The ResonanceWRight class.
// Derived class for W_R+- properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceWRight::initConstants() {

  // Locally stored properties and couplings.
  thetaWRat     = 1. / (12. * couplingsPtr->sin2thetaW());

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceWRight::calcPreFac(bool) {

  // Common coupling factors.
  alpEM     = couplingsPtr->alphaEM(mHat * mHat);
  alpS      = couplingsPtr->alphaS(mHat * mHat);
  colQ      = 3. * (1. + alpS / M_PI);
  preFac    = alpEM * thetaWRat * mHat;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceWRight::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // Combine kinematics with colour factor and CKM couplings.
  widNow    = preFac * (1. - 0.5 * (mr1 + mr2) - 0.5 * pow2(mr1 - mr2))
            * ps;
  if (id1Abs < 9) widNow *= colQ * couplingsPtr->V2CKMid(id1Abs, id2Abs);

}

//==========================================================================

// The ResonanceHchgchgLeft class.
// Derived class for H++/H-- (left) properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceHchgchgLeft::initConstants() {

  // Read in Yukawa matrix for couplings to a lepton pair.
  yukawa[1][1]  = settingsPtr->parm("LeftRightSymmmetry:coupHee");
  yukawa[2][1]  = settingsPtr->parm("LeftRightSymmmetry:coupHmue");
  yukawa[2][2]  = settingsPtr->parm("LeftRightSymmmetry:coupHmumu");
  yukawa[3][1]  = settingsPtr->parm("LeftRightSymmmetry:coupHtaue");
  yukawa[3][2]  = settingsPtr->parm("LeftRightSymmmetry:coupHtaumu");
  yukawa[3][3]  = settingsPtr->parm("LeftRightSymmmetry:coupHtautau");

  // Locally stored properties and couplings.
  gL            = settingsPtr->parm("LeftRightSymmmetry:gL");
  vL            = settingsPtr->parm("LeftRightSymmmetry:vL");
  mW            = particleDataPtr->m0(24);

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceHchgchgLeft::calcPreFac(bool) {

  // Common coupling factors.
  preFac        = mHat / (8. * M_PI);

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceHchgchgLeft::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // H++-- width to a pair of leptons. Combinatorial factor of 2.
  if (id1Abs < 17 && id2Abs < 17) {
    widNow    = preFac * pow2(yukawa[(id1Abs-9)/2][(id2Abs-9)/2]) * ps;
    if (id2Abs != id1Abs) widNow *= 2.;
  }

  // H++-- width to a pair of lefthanded W's.
  else if (id1Abs == 24 && id2Abs == 24)
    widNow    = preFac * 0.5 * pow2(gL*gL * vL / mW)
              * (3. * mr1 + 0.25 / mr1 - 1.) * ps;

}

//==========================================================================

// The ResonanceHchgchgRight class.
// Derived class for H++/H-- (right) properties.

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceHchgchgRight::initConstants() {

  // Read in Yukawa matrix for couplings to a lepton pair.
  yukawa[1][1]  = settingsPtr->parm("LeftRightSymmmetry:coupHee");
  yukawa[2][1]  = settingsPtr->parm("LeftRightSymmmetry:coupHmue");
  yukawa[2][2]  = settingsPtr->parm("LeftRightSymmmetry:coupHmumu");
  yukawa[3][1]  = settingsPtr->parm("LeftRightSymmmetry:coupHtaue");
  yukawa[3][2]  = settingsPtr->parm("LeftRightSymmmetry:coupHtaumu");
  yukawa[3][3]  = settingsPtr->parm("LeftRightSymmmetry:coupHtautau");

  // Locally stored properties and couplings.
  idWR          = 9000024;
  gR            = settingsPtr->parm("LeftRightSymmmetry:gR");

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceHchgchgRight::calcPreFac(bool) {

  // Common coupling factors.
  preFac        = mHat / (8. * M_PI);

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceHchgchgRight::calcWidth(bool) {

  // Check that above threshold.
  if (ps == 0.) return;

  // H++-- width to a pair of leptons. Combinatorial factor of 2.
  if (id1Abs < 17 && id2Abs < 17) {
    widNow    = preFac * pow2(yukawa[(id1Abs-9)/2][(id2Abs-9)/2]) * ps;
    if (id2Abs != id1Abs) widNow *= 2.;
  }

  // H++-- width to a pair of lefthanded W's.
  else if (id1Abs == idWR && id2Abs == idWR)
    widNow    = preFac * pow2(yukawa[(id1Abs-9)/2][(id2Abs-9)/2]) * ps;

}

//==========================================================================

} // end namespace Pythia8
