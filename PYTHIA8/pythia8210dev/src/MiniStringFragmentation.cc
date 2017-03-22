// MiniStringFragmentation.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the .
// MiniStringFragmentation class

#include "Pythia8/MiniStringFragmentation.h"

namespace Pythia8 {

//==========================================================================

// The MiniStringFragmentation class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Since diffractive by definition is > 1 particle, try hard.
const int MiniStringFragmentation::NTRYDIFFRACTIVE = 200;

// After one-body fragmentation failed, try two-body once more.
const int MiniStringFragmentation::NTRYLASTRESORT  = 100;

// Loop try to combine available endquarks to valid hadron.
const int MiniStringFragmentation::NTRYFLAV        = 10;

//--------------------------------------------------------------------------

// Initialize and save pointers.

void MiniStringFragmentation::init(Info* infoPtrIn, Settings& settings,
   ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
   StringFlav* flavSelPtrIn, StringPT* pTSelPtrIn, StringZ* zSelPtrIn) {

  // Save pointers.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  flavSelPtr      = flavSelPtrIn;
  pTSelPtr        = pTSelPtrIn;
  zSelPtr         = zSelPtrIn;

  // Initialize the MiniStringFragmentation class proper.
  nTryMass        = settings.mode("MiniStringFragmentation:nTry");

  // Initialize the b parameter of the z spectrum, used when joining jets.
  bLund           = zSelPtr->bAreaLund();

}

//--------------------------------------------------------------------------

// Do the fragmentation: driver routine.

bool MiniStringFragmentation::fragment(int iSub, ColConfig& colConfig,
  Event& event, bool isDiff) {

  // Junction topologies not yet considered - is very rare.
  iParton  = colConfig[iSub].iParton;
  if (iParton.front() < 0) {
    infoPtr->errorMsg("Error in MiniStringFragmentation::fragment: "
      "very low-mass junction topologies not yet handled");
    return false;
  }

  // Read in info on system to be treated.
  flav1    = FlavContainer( event[ iParton.front() ].id() );
  flav2    = FlavContainer( event[ iParton.back() ].id() );
  pSum     = colConfig[iSub].pSum;
  mSum     = colConfig[iSub].mass;
  m2Sum    = mSum*mSum;
  isClosed = colConfig[iSub].isClosed;

  // Do not want diffractive systems to easily collapse to one particle.
  int nTryFirst = (isDiff) ? NTRYDIFFRACTIVE : nTryMass;

  // First try to produce two particles from the system.
  if (ministring2two( nTryFirst, event)) return true;

  // If this fails, then form one hadron and shuffle momentum.
  if (ministring2one( iSub, colConfig, event)) return true;

  // If also this fails, then try harder to produce two particles.
  if (ministring2two( NTRYLASTRESORT, event)) return true;

  // Else complete failure.
  infoPtr->errorMsg("Error in MiniStringFragmentation::fragment: "
      "no 1- or 2-body state found above mass threshold");
  return false;

}

//--------------------------------------------------------------------------

  // Attempt to produce two particles from the ministring.

bool MiniStringFragmentation::ministring2two( int nTry, Event& event) {

  // Properties of the produced hadrons.
  int    idHad1  = 0;
  int    idHad2  = 0;
  double mHad1   = 0.;
  double mHad2   = 0.;
  double mHadSum = 0.;

  // Allow a few attempts to find a particle pair with low enough masses.
  for (int iTry = 0; iTry < nTry; ++iTry) {

    // For closed gluon loop need to pick an initial flavour.
    if (isClosed) do {
      int idStart = flavSelPtr->pickLightQ();
      FlavContainer flavStart(idStart, 1);
      flavStart = flavSelPtr->pick( flavStart);
      flav1 = flavSelPtr->pick( flavStart);
      flav2.anti(flav1);
    } while (flav1.id == 0 || flav1.nPop > 0);

    // Create a new q qbar flavour to form two hadrons.
    // Start from a diquark, if any.
    do {
      FlavContainer flav3 =
        (flav1.isDiquark() || (!flav2.isDiquark() && rndmPtr->flat() < 0.5) )
        ? flavSelPtr->pick( flav1) : flavSelPtr->pick( flav2).anti();
      idHad1 = flavSelPtr->combine( flav1, flav3);
      idHad2 = flavSelPtr->combine( flav2, flav3.anti());
    } while (idHad1 == 0 || idHad2 == 0);

    // Check whether the mass sum fits inside the available phase space.
    mHad1 = particleDataPtr->mSel(idHad1);
    mHad2 = particleDataPtr->mSel(idHad2);
    mHadSum = mHad1 + mHad2;
    if (mHadSum < mSum) break;
  }
  if (mHadSum >= mSum) return false;

  // Define an effective two-parton string, by splitting intermediate
  // gluon momenta in proportion to their closeness to either endpoint.
  Vec4 pSum1 = event[ iParton.front() ].p();
  Vec4 pSum2 = event[ iParton.back() ].p();
  if (iParton.size() > 2) {
    Vec4 pEnd1 = pSum1;
    Vec4 pEnd2 = pSum2;
    Vec4 pEndSum = pEnd1 + pEnd2;
    for (int i = 1; i < int(iParton.size()) - 1 ; ++i) {
      Vec4 pNow = event[ iParton[i] ].p();
      double ratio = (pEnd2 * pNow) / (pEndSum * pNow);
      pSum1 += ratio * pNow;
      pSum2 += (1. - ratio) * pNow;
    }
  }

  // If split did not provide an axis then pick random axis to break tie.
  // (Needed for low-mass q-g-qbar with q-qbar perfectly parallel.)
  if (pSum1.mCalc() + pSum2.mCalc() > 0.999999 * mSum) {
    double cthe = 2. * rndmPtr->flat() - 1.;
    double sthe = sqrtpos(1. - cthe * cthe);
    double phi  = 2. * M_PI * rndmPtr->flat();
    Vec4 delta  = 0.5 * min( pSum1.e(), pSum2.e())
        * Vec4( sthe * sin(phi), sthe * cos(phi), cthe, 0.);
    pSum1 += delta;
    pSum2 -= delta;
    infoPtr->errorMsg("Warning in MiniStringFragmentation::ministring2two: "
      "random axis needed to break tie");
  }

  // Set up a string region based on the two effective endpoints.
  StringRegion region;
  region.setUp( pSum1, pSum2);

  // Generate an isotropic decay in the ministring rest frame,
  // suppressed at large pT by a fragmentation pT Gaussian.
  double pAbs2 = 0.25 * ( pow2(m2Sum - mHad1*mHad1 - mHad2*mHad2)
    - pow2(2. * mHad1 * mHad2) ) / m2Sum;
  double pT2 = 0.;
  do {
    double cosTheta = rndmPtr->flat();
    pT2 = (1. - pow2(cosTheta)) * pAbs2;
  } while (pTSelPtr->suppressPT2(pT2) < rndmPtr->flat() );

  // Construct the forward-backward asymmetry of the two particles.
  double mT21 = mHad1*mHad1 + pT2;
  double mT22 = mHad2*mHad2 + pT2;
  double lambda = sqrtpos( pow2(m2Sum  - mT21 - mT22) - 4. * mT21 * mT22 );
  double probReverse = 1. / (1. + exp( min( 50., bLund * lambda) ) );

  // Construct kinematics, as viewed in the transverse rest frame.
  double xpz1 = 0.5 * lambda/ m2Sum;
  if (probReverse > rndmPtr->flat()) xpz1 = -xpz1;
  double xmDiff = (mT21 - mT22) / m2Sum;
  double xe1 = 0.5 * (1. + xmDiff);
  double xe2 = 0.5 * (1. - xmDiff );

  // Distribute pT isotropically in angle.
  double phi = 2. * M_PI * rndmPtr->flat();
  double pT  = sqrt(pT2);
  double px  = pT * cos(phi);
  double py  = pT * sin(phi);

  // Translate this into kinematics in the string frame.
  Vec4 pHad1 = region.pHad( xe1 + xpz1, xe1 - xpz1,  px,  py);
  Vec4 pHad2 = region.pHad( xe2 - xpz1, xe2 + xpz1, -px, -py);

  // Mark hadrons from junction fragmentation with different status.
  int statusHadPos = 82, statusHadNeg = 82;
  if (abs(idHad1) > 1000 && abs(idHad1) < 10000 &&
      abs(idHad2) > 1000 && abs(idHad2) < 10000) {
    if (event[ iParton.front() ].statusAbs() == 74) statusHadPos = 89;
    if (event[ iParton.back() ].statusAbs() == 74)  statusHadNeg = 89;
  }
  else if (abs(idHad1) > 1000 && abs(idHad1) < 10000) {
    if (event[ iParton.front() ].statusAbs() == 74 ||
        event[ iParton.back() ].statusAbs() == 74) statusHadPos = 89;
  }
  else if (abs(idHad2) > 1000 && abs(idHad2) < 10000) {
    if (event[ iParton.front() ].statusAbs() == 74 ||
        event[ iParton.back() ].statusAbs() == 74) statusHadNeg = 89;
  }
  // Add produced particles to the event record.
  int iFirst = event.append( idHad1, statusHadPos, iParton.front(),
    iParton.back(), 0, 0, 0, 0, pHad1, mHad1);
  int iLast = event.append( idHad2, statusHadNeg, iParton.front(),
    iParton.back(), 0, 0, 0, 0, pHad2, mHad2);

  // Set decay vertex when this is displaced.
  if (event[iParton.front()].hasVertex()) {
    Vec4 vDec = event[iParton.front()].vDec();
    event[iFirst].vProd( vDec );
    event[iLast].vProd( vDec );
  }

  // Set lifetime of hadrons.
  event[iFirst].tau( event[iFirst].tau0() * rndmPtr->exp() );
  event[iLast].tau( event[iLast].tau0() * rndmPtr->exp() );

  // Mark original partons as hadronized and set their daughter range.
  for (int i = 0; i < int(iParton.size()); ++i) {
    event[ iParton[i] ].statusNeg();
    event[ iParton[i] ].daughters(iFirst, iLast);
  }

  // Successfully done.
  return true;

}

//--------------------------------------------------------------------------

// Attempt to produce one particle from a ministring.
// Current algorithm: find the system with largest invariant mass
// relative to the existing one, and boost that system appropriately.
// Try more sophisticated alternatives later?? (Z0 mass shifted??)
// Also, if problems, attempt several times to obtain closer mass match??

bool MiniStringFragmentation::ministring2one( int iSub,
  ColConfig& colConfig, Event& event) {

  // Cannot handle qq + qbarqbar system.
  if (abs(flav1.id) > 100 && abs(flav2.id) > 100) return false;

  // For closed gluon loop need to pick an initial flavour.
  if (isClosed) do {
    int idStart = flavSelPtr->pickLightQ();
    FlavContainer flavStart(idStart, 1);
    flav1 = flavSelPtr->pick( flavStart);
    flav2 = flav1.anti();
  } while (abs(flav1.id) > 100);

  // Select hadron flavour from available quark flavours.
  int idHad = 0;
  for (int iTryFlav = 0; iTryFlav < NTRYFLAV; ++iTryFlav) {
    idHad = flavSelPtr->combine( flav1, flav2);
    if (idHad != 0) break;
  }
  if (idHad == 0) return false;

  // Find mass.
  double mHad = particleDataPtr->mSel(idHad);

  // Find the untreated parton system which combines to the largest
  // squared mass above mimimum required.
  int iMax = -1;
  double deltaM2 = mHad*mHad - mSum*mSum;
  double delta2Max = 0.;
  for (int iRec = iSub + 1; iRec < colConfig.size(); ++iRec) {
    double delta2Rec = 2. * (pSum * colConfig[iRec].pSum) - deltaM2
      - 2. * mHad * colConfig[iRec].mass;
    if (delta2Rec > delta2Max) { iMax = iRec; delta2Max = delta2Rec;}
  }
  if (iMax == -1) return false;

  // Construct kinematics of the hadron and recoiling system.
  Vec4& pRec     = colConfig[iMax].pSum;
  double mRec    = colConfig[iMax].mass;
  double vecProd = pSum * pRec;
  double coefOld = mSum*mSum + vecProd;
  double coefNew = mHad*mHad + vecProd;
  double coefRec = mRec*mRec + vecProd;
  double coefSum = coefOld + coefNew;
  double sHat    = coefOld + coefRec;
  double root    = sqrtpos( (pow2(coefSum) - 4. * sHat * mHad*mHad)
    / (pow2(vecProd) - pow2(mSum * mRec)) );
  double k2      = 0.5 * (coefOld * root - coefSum) / sHat;
  double k1      = (coefRec * k2 + 0.5 * deltaM2) / coefOld;
  Vec4 pHad      = (1. + k1) * pSum - k2 * pRec;
  Vec4 pRecNew   = (1. + k2) * pRec - k1 * pSum;

  // Mark hadrons from junction split off with status 89.
  int statusHad = 81;
  if (abs(idHad) > 1000 && abs(idHad) < 10000 &&
      (event[ iParton.front() ].statusAbs() == 74 ||
       event[ iParton.back() ].statusAbs() == 74)) statusHad = 89;

  // Add the produced particle to the event record.
  int iHad = event.append( idHad, statusHad, iParton.front(), iParton.back(),
    0, 0, 0, 0, pHad, mHad);

  // Set decay vertex when this is displaced.
  if (event[iParton.front()].hasVertex()) {
    Vec4 vDec = event[iParton.front()].vDec();
    event[iHad].vProd( vDec );
  }

  // Set lifetime of hadron.
  event[iHad].tau( event[iHad].tau0() * rndmPtr->exp() );

  // Mark original partons as hadronized and set their daughter range.
  for (int i = 0; i < int(iParton.size()); ++i) {
    event[ iParton[i] ].statusNeg();
    event[ iParton[i] ].daughters(iHad, iHad);
  }

  // Copy down recoiling system, with boosted momentum. Update current partons.
  RotBstMatrix M;
  M.bst(pRec, pRecNew);
  for (int i = 0; i < colConfig[iMax].size(); ++i) {
    int iOld = colConfig[iMax].iParton[i];
    // Do not touch negative iOld = beginning of new junction leg.
    if (iOld >= 0) {
      int iNew;
      // Keep track of 74 throughout the event.
      if (event[iOld].status() == 74) iNew = event.copy(iOld, 74);
      else iNew = event.copy(iOld, 72);
      event[iNew].rotbst(M);
      colConfig[iMax].iParton[i] = iNew;
    }
  }
  colConfig[iMax].pSum = pRecNew;
  colConfig[iMax].isCollected = true;

  // Successfully done.
  return true;

}

//==========================================================================

} // end namespace Pythia8
