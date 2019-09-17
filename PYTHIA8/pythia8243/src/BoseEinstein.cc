// BoseEinstein.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the BoseEinsten class.

#include "Pythia8/BoseEinstein.h"

namespace Pythia8 {

//==========================================================================

// The BoseEinstein class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Enumeration of id codes and table for particle species considered.
const int    BoseEinstein::IDHADRON[9] = { 211, -211, 111, 321, -321,
                                           130,  310, 221, 331 };
const int    BoseEinstein::ITABLE[9] = { 0, 0, 0, 1, 1, 1, 1, 2, 3 };

// Distance between table entries, normalized to min( 2*mass, QRef).
const double BoseEinstein::STEPSIZE  = 0.05;

// Skip shift for two extremely close particles, to avoid instabilities.
const double BoseEinstein::Q2MIN     = 1e-8;

// Parameters of energy compensation procedure: maximally allowed
// relative energy error, iterative stepsize, and number of iterations.
const double BoseEinstein::COMPRELERR = 1e-10;
const double BoseEinstein::COMPFACMAX = 1000.;
const int    BoseEinstein::NCOMPSTEP  = 10;

//--------------------------------------------------------------------------

// Find settings. Precalculate table used to find momentum shifts.

bool BoseEinstein::init(Info* infoPtrIn, Settings& settings,
  ParticleData& particleData) {

  // Save pointer.
  infoPtr         = infoPtrIn;

  // Main flags.
  doPion   = settings.flag("BoseEinstein:Pion");
  doKaon   = settings.flag("BoseEinstein:Kaon");
  doEta    = settings.flag("BoseEinstein:Eta");

  // Shape of Bose-Einstein enhancement/suppression.
  lambda   = settings.parm("BoseEinstein:lambda");
  QRef     = settings.parm("BoseEinstein:QRef");

  // Multiples and inverses (= "radii") of distance parameters in Q-space.
  QRef2    = 2. * QRef;
  QRef3    = 3. * QRef;
  R2Ref    = 1. / (QRef * QRef);
  R2Ref2   = 1. / (QRef2 * QRef2);
  R2Ref3   = 1. / (QRef3 * QRef3);

  // Masses of particles with Bose-Einstein implemented.
  for (int iSpecies = 0; iSpecies < 9; ++iSpecies)
    mHadron[iSpecies] = particleData.m0( IDHADRON[iSpecies] );

  // Pair pi, K, eta and eta' masses for use in tables.
  mPair[0] = 2. * mHadron[0];
  mPair[1] = 2. * mHadron[3];
  mPair[2] = 2. * mHadron[7];
  mPair[3] = 2. * mHadron[8];

  // Loop over the four required tables. Local variables.
  double Qnow, Q2now, centerCorr;
  for (int iTab = 0; iTab < 4; ++iTab) {
    m2Pair[iTab]      = mPair[iTab] * mPair[iTab];

    // Step size and number of steps in normal table.
    deltaQ[iTab]      = STEPSIZE * min(mPair[iTab], QRef);
    nStep[iTab]       = min( 199, 1 + int(3. * QRef / deltaQ[iTab]) );
    maxQ[iTab]        = (nStep[iTab] - 0.1) * deltaQ[iTab];
    centerCorr        = deltaQ[iTab] * deltaQ[iTab] / 12.;

    // Construct normal table recursively in Q space.
    shift[iTab][0]    = 0.;
    for (int i = 1; i <= nStep[iTab]; ++i) {
      Qnow            = deltaQ[iTab] * (i - 0.5);
      Q2now           = Qnow * Qnow;
      shift[iTab][i]  = shift[iTab][i - 1] + exp(-Q2now * R2Ref)
        * deltaQ[iTab] * (Q2now + centerCorr) / sqrt(Q2now + m2Pair[iTab]);
    }

    // Step size and number of steps in compensation table.
    deltaQ3[iTab]     = STEPSIZE * min(mPair[iTab], QRef3);
    nStep3[iTab]      = min( 199, 1 + int(9. * QRef / deltaQ3[iTab]) );
    maxQ3[iTab]       = (nStep3[iTab] - 0.1) * deltaQ3[iTab];
    centerCorr        = deltaQ3[iTab] * deltaQ3[iTab] / 12.;

    // Construct compensation table recursively in Q space.
    shift3[iTab][0]   = 0.;
    for (int i = 1; i <= nStep3[iTab]; ++i) {
      Qnow            = deltaQ3[iTab] * (i - 0.5);
      Q2now           = Qnow * Qnow;
      shift3[iTab][i] = shift3[iTab][i - 1] + exp(-Q2now * R2Ref3)
        * deltaQ3[iTab] * (Q2now + centerCorr) / sqrt(Q2now + m2Pair[iTab]);
    }

  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Perform Bose-Einstein corrections on an event.

bool BoseEinstein::shiftEvent( Event& event) {

  // Reset list of identical particles.
  hadronBE.resize(0);

  // Loop over all hadron species with BE effects.
  nStored[0] = 0;
  for (int iSpecies = 0; iSpecies < 9; ++iSpecies) {
    nStored[iSpecies + 1] = nStored[iSpecies];
    if (!doPion && iSpecies <= 2) continue;
    if (!doKaon && iSpecies >= 3 && iSpecies <= 6) continue;
    if (!doEta  && iSpecies >= 7) continue;

    // Properties of current hadron species.
    int idNow = IDHADRON[ iSpecies ];
    int iTab  = ITABLE[ iSpecies ];

    // Loop through event record to store copies of current species.
    for (int i = 0; i < event.size(); ++i)
      if ( event[i].id() == idNow && event[i].isFinal() )
        hadronBE.push_back(
          BoseEinsteinHadron( idNow, i, event[i].p(), event[i].m() ) );
    nStored[iSpecies + 1] = hadronBE.size();

    // Loop through pairs of identical particles and find shifts.
    for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
    for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
      shiftPair( i1, i2, iTab);
  }

  // Must have at least two pairs to carry out compensation.
  if (nStored[9] < 2) return true;

  // Shift momenta and recalculate energies.
  double eSumOriginal = 0.;
  double eSumShifted  = 0.;
  double eDiffByComp  = 0.;
  for (int i = 0; i < nStored[9]; ++i) {
    eSumOriginal  += hadronBE[i].p.e();
    hadronBE[i].p += hadronBE[i].pShift;
    hadronBE[i].p.e( sqrt( hadronBE[i].p.pAbs2() + hadronBE[i].m2 ) );
    eSumShifted   += hadronBE[i].p.e();
    eDiffByComp   += dot3( hadronBE[i].pComp, hadronBE[i].p)
                     / hadronBE[i].p.e();
  }

  // Iterate compensation shift until convergence.
  int iStep = 0;
  while ( abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal
    && abs(eSumShifted - eSumOriginal) < COMPFACMAX * abs(eDiffByComp)
    && iStep < NCOMPSTEP ) {
    ++iStep;
    double compFac   = (eSumOriginal - eSumShifted) / eDiffByComp;
    eSumShifted      = 0.;
    eDiffByComp      = 0.;
    for (int i = 0; i < nStored[9]; ++i) {
      hadronBE[i].p += compFac * hadronBE[i].pComp;
      hadronBE[i].p.e( sqrt( hadronBE[i].p.pAbs2() + hadronBE[i].m2 ) );
      eSumShifted   += hadronBE[i].p.e();
      eDiffByComp   += dot3( hadronBE[i].pComp, hadronBE[i].p)
                       / hadronBE[i].p.e();
    }
  }

  // Error if no convergence, and then return without doing BE shift.
  // However, not grave enough to kill event, so return true.
  if ( abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal ) {
    infoPtr->errorMsg("Warning in BoseEinstein::shiftEvent: "
      "no consistent BE shift topology found, so skip BE");
    return true;
  }

  // Store new particle copies with shifted momenta.
  for (int i = 0; i < nStored[9]; ++i) {
    int iNew = event.copy( hadronBE[i].iPos, 99);
    event[ iNew ].p( hadronBE[i].p );
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Calculate shift and (unnormalized) compensation for pair.

void BoseEinstein::shiftPair( int i1, int i2, int iTab) {

  // Calculate old relative momentum.
  double Q2old = m2(hadronBE[i1].p, hadronBE[i2].p) - m2Pair[iTab];
  if (Q2old < Q2MIN) return;
  double Qold  = sqrt(Q2old);
  double psFac = sqrt(Q2old + m2Pair[iTab]) / Q2old;

  // Calculate new relative momentum for normal shift.
  double Qmove = 0.;
  if (Qold < deltaQ[iTab]) Qmove = Qold / 3.;
  else if (Qold < maxQ[iTab]) {
    double realQbin = Qold / deltaQ[iTab];
    int    intQbin  = int( realQbin );
    double inter    = (pow3(realQbin) - pow3(intQbin))
      / (3 * intQbin * (intQbin + 1) + 1);
    Qmove = ( shift[iTab][intQbin] + inter * (shift[iTab][intQbin + 1]
      - shift[iTab][intQbin]) ) * psFac;
  }
  else Qmove = shift[iTab][nStep[iTab]] * psFac;
  double Q2new = Q2old * pow( Qold / (Qold + 3. * lambda * Qmove), 2. / 3.);

  // Calculate corresponding three-momentum shift.
  double Q2Diff    = Q2new - Q2old;
  double p2DiffAbs = (hadronBE[i1].p - hadronBE[i2].p).pAbs2();
  double p2AbsDiff = hadronBE[i1].p.pAbs2() - hadronBE[i2].p.pAbs2();
  double eSum      = hadronBE[i1].p.e() + hadronBE[i2].p.e();
  double eDiff     = hadronBE[i1].p.e() - hadronBE[i2].p.e();
  double sumQ2E    = Q2Diff + eSum * eSum;
  double rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
  double rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
  double factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
    + Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

  // Add shifts to sum. (Energy component dummy.)
  Vec4   pDiff     = factor * (hadronBE[i1].p - hadronBE[i2].p);
  hadronBE[i1].pShift += pDiff;
  hadronBE[i2].pShift -= pDiff;

  // Calculate new relative momentum for compensation shift.
  double Qmove3 = 0.;
  if (Qold < deltaQ3[iTab]) Qmove3 = Qold / 3.;
  else if (Qold < maxQ3[iTab]) {
    double realQbin = Qold / deltaQ3[iTab];
    int    intQbin  = int( realQbin );
    double inter    = (pow3(realQbin) - pow3(intQbin))
      / (3 * intQbin * (intQbin + 1) + 1);
    Qmove3 = ( shift3[iTab][intQbin] + inter * (shift3[iTab][intQbin + 1]
      - shift3[iTab][intQbin]) ) * psFac;
  }
  else Qmove3 = shift3[iTab][nStep3[iTab]] *psFac;
  double Q2new3 = Q2old * pow( Qold / (Qold + 3. * lambda * Qmove3), 2. / 3.);

  // Calculate corresponding three-momentum shift.
  Q2Diff    = Q2new3 - Q2old;
  sumQ2E    = Q2Diff + eSum * eSum;
  rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
  rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
  factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
    + Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

  // Extra dampening factor to go from BE_3 to BE_32.
  factor   *= 1. - exp(-Q2old * R2Ref2);

  // Add shifts to sum. (Energy component dummy.)
  pDiff     = factor * (hadronBE[i1].p - hadronBE[i2].p);
  hadronBE[i1].pComp += pDiff;
  hadronBE[i2].pComp -= pDiff;

}

//==========================================================================

} // end namespace Pythia8
