// FragmentationSystems.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// ColConfig, StringRegion and StringSystem classes.

#include "Pythia8/FragmentationSystems.h"

namespace Pythia8 {

//==========================================================================

// The ColConfig class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// A typical u/d constituent mass.
const double ColConfig::CONSTITUENTMASS = 0.325;

//--------------------------------------------------------------------------

// Initialize and save pointers.

void ColConfig::init(Info* infoPtrIn, Settings& settings,
  StringFlav* flavSelPtrIn) {

  // Save pointers.
  infoPtr       = infoPtrIn;
  flavSelPtr    = flavSelPtrIn;

  // Joining of nearby partons along the string.
  mJoin         = settings.parm("FragmentationSystems:mJoin");

  // For consistency ensure that mJoin is bigger than in StringRegion.
  mJoin         = max( mJoin, 2. * StringRegion::MJOIN);

  // Simplification of q q q junction topology to quark - diquark one.
  mJoinJunction = settings.parm("FragmentationSystems:mJoinJunction");
  mStringMin    = settings.parm("HadronLevel:mStringMin");

}

//--------------------------------------------------------------------------

// Insert a new colour singlet system in ascending mass order.
// Calculate its properties. Join nearby partons.

bool ColConfig::insert( vector<int>& iPartonIn, Event& event) {

  // Find momentum and invariant mass of system, minus endpoint masses.
  Vec4 pSumIn;
  double mSumIn = 0.;
  bool hasJunctionIn = false;
  int  nJunctionLegs = 0;
  for (int i = 0; i < int(iPartonIn.size()); ++i) {
    if (iPartonIn[i] < 0) {
      hasJunctionIn = true;
      ++nJunctionLegs;
    } else {
      pSumIn += event[ iPartonIn[i] ].p();
      if (!event[ iPartonIn[i] ].isGluon())
        mSumIn += event[ iPartonIn[i] ].constituentMass();
    }
  }
  double massIn = pSumIn.mCalc();
  double massExcessIn = massIn - mSumIn;

  // Check for rare triple- and higher junction systems (like J-Jbar-J)
  if (nJunctionLegs >= 5) {
    infoPtr->errorMsg("Error in ColConfig::insert: "
      "junction topology too complicated; too many junction legs");
    return false;
  }
  // Check that junction systems have at least three legs.
  else if (nJunctionLegs > 0 && nJunctionLegs <= 2) {
    infoPtr->errorMsg("Error in ColConfig::insert: "
      "junction topology inconsistent; too few junction legs");
    return false;
  }

  // Check that momenta do not contain not-a-number.
  if (abs(massExcessIn) >= 0.);
  else {
    infoPtr->errorMsg("Error in ColConfig::insert: "
      "not-a-number system mass");
    return false;
  }

  // Identify closed gluon loop. Assign "endpoint" masses as light quarks.
  bool isClosedIn = (iPartonIn[0] >= 0 && event[ iPartonIn[0] ].col() != 0
    && event[ iPartonIn[0] ].acol() != 0 );
  if (isClosedIn) massExcessIn -= 2. * CONSTITUENTMASS;

  // For junction topology: join two nearby legs into a diquark.
  if (hasJunctionIn && joinJunction( iPartonIn, event, massExcessIn))
    hasJunctionIn = false;

  // Loop while > 2 partons left and hope of finding joining pair.
  bool hasJoined = true;
  while (hasJoined && iPartonIn.size() > 2) {

    // Look for the pair of neighbour partons (along string) with
    // the smallest invariant mass (subtracting quark masses).
    int iJoinMin = -1;
    double mJoinMin = 2. * mJoin;
    int nSize = iPartonIn.size();
    int nPair = (isClosedIn) ? nSize : nSize - 1;
    for (int i = 0; i < nPair; ++i) {
      // Keep three legs of junction separate.
      if (iPartonIn[i] < 0 || iPartonIn[(i + 1)%nSize] < 0) continue;
      Particle& parton1 = event[ iPartonIn[i] ];
      Particle& parton2 = event[ iPartonIn[(i + 1)%nSize] ];
      // Avoid joining non-partons, e.g. gluino/squark for R-hadron.
      if (!parton1.isParton() || !parton2.isParton()) continue;
      Vec4 pSumNow;
      pSumNow += (parton1.isGluon()) ? 0.5 * parton1.p() : parton1.p();
      pSumNow += (parton2.isGluon()) ? 0.5 * parton2.p() : parton2.p();
      double mJoinNow = pSumNow.mCalc();
      if (!parton1.isGluon()) mJoinNow -= parton1.m();
      if (!parton2.isGluon()) mJoinNow -= parton2.m();
      if (mJoinNow < mJoinMin) { iJoinMin = i; mJoinMin = mJoinNow; }
    }

    // If sufficiently nearby then join into one new parton.
    // Note: error sensitivity to mJoin indicates unstable precedure??
    hasJoined = false;
    if (mJoinMin < mJoin) {
      int iJoin1  = iPartonIn[iJoinMin];
      int iJoin2  = iPartonIn[(iJoinMin + 1)%nSize];
      int idNew   = (event[iJoin1].isGluon()) ? event[iJoin2].id()
                                              : event[iJoin1].id();
      int iMoth1  = min(iJoin1, iJoin2);
      int iMoth2  = max(iJoin1, iJoin2);
      // When g + q -> q flip to ensure that mother1 = q.
      if (event[iMoth1].id() == 21 && event[iMoth2].id() != 21)
        swap( iMoth1, iMoth2);
      int colNew  = event[iJoin1].col();
      int acolNew = event[iJoin2].acol();
      if (colNew == acolNew) {
        colNew    = event[iJoin2].col();
        acolNew   = event[iJoin1].acol();
      }
      Vec4 pNew   = event[iJoin1].p() + event[iJoin2].p();

      int statusHad = 73;
      // Need to keep status as 74 for junctions in order to keep track
      // of them.
      if (event[iMoth1].statusAbs() == 74) statusHad = 74;

      // Append joined parton to event record.
      int iNew = event.append( idNew, statusHad, iMoth1, iMoth2, 0, 0,
        colNew, acolNew, pNew, pNew.mCalc() );

      // Displaced lifetime/vertex; mothers should be same but prefer quark.
      int iVtx = (event[iJoin1].isGluon()) ? iJoin2 : iJoin1;
      event[iNew].tau( event[iVtx].tau() );
      if (event[iVtx].hasVertex()) event[iNew].vProd( event[iVtx].vProd() );

      // Mark joined partons and reduce remaining system.
      event[iJoin1].statusNeg();
      event[iJoin2].statusNeg();
      event[iJoin1].daughter1(iNew);
      event[iJoin2].daughter1(iNew);
      if (iJoinMin == nSize - 1) iPartonIn[0] = iNew;
      else {
        iPartonIn[iJoinMin] = iNew;
        for (int i = iJoinMin + 1; i < nSize - 1; ++i)
          iPartonIn[i] = iPartonIn[i + 1];
      }
      iPartonIn.pop_back();

      // If joined,then loopback to look for more.
      hasJoined = true;
    }
  }

  // Store new colour singlet system at the end.
  singlets.push_back( ColSinglet(iPartonIn, pSumIn, massIn,
    massExcessIn, hasJunctionIn, isClosedIn) );

  // Now move around, so that smallest mass excesses come first.
  int iInsert = singlets.size() - 1;
  for (int iSub = singlets.size() - 2; iSub >= 0; --iSub) {
    if (massExcessIn > singlets[iSub].massExcess) break;
    singlets[iSub + 1] = singlets[iSub];
    iInsert = iSub;
  }
  if (iInsert < int(singlets.size()) - 1) singlets[iInsert] =
    ColSinglet(iPartonIn, pSumIn, massIn, massExcessIn,
    hasJunctionIn, isClosedIn);

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Join two legs of junction to a diquark for small invariant masses.
// Note: for junction system, iPartonIn points to structure
// (-code0) g...g.q0 (-code1) g...g.q1 (-code2) g...g.q2

bool ColConfig::joinJunction( vector<int>& iPartonIn, Event& event,
  double massExcessIn) {

  // Find four-momentum and endpoint quarks and masses on the three legs.
  Vec4   pLeg[3];
  double mLeg[3] = { 0., 0., 0.};
  int    idAbsLeg[3];
  int leg = -1;
  for (int i = 0; i < int(iPartonIn.size()); ++ i) {
    if (iPartonIn[i] < 0) ++leg;
    else {
      pLeg[leg]    += event[ iPartonIn[i] ].p();
      mLeg[leg]     = event[ iPartonIn[i] ].m();
      idAbsLeg[leg] = event[ iPartonIn[i] ].idAbs();
    }
  }

  // Calculate invariant mass of three pairs, minus endpoint masses.
  double m01  = (pLeg[0] + pLeg[1]).mCalc() - mLeg[0] - mLeg[1];
  double m02  = (pLeg[0] + pLeg[2]).mCalc() - mLeg[0] - mLeg[2];
  double m12  = (pLeg[1] + pLeg[2]).mCalc() - mLeg[1] - mLeg[2];

  // Find lowest-mass pair not involving diquark.
  double mMin = mJoinJunction + 1.;
  int    legA = -1;
  int    legB = -1;
  if (m01 < mMin && idAbsLeg[0] < 9 && idAbsLeg[1] < 9) {
    mMin = m01;
    legA = 0;
    legB = 1;
  }
  if (m02 < mMin && idAbsLeg[0] < 9 && idAbsLeg[2] < 9) {
    mMin = m02;
    legA = 0;
    legB = 2;
  }
  if (m12 < mMin && idAbsLeg[1] < 9 && idAbsLeg[2] < 9) {
    mMin = m12;
    legA = 1;
    legB = 2;
  }
  int legC = 3 - legA - legB;

  // Nothing to do if no two legs have small invariant mass, and
  // system as a whole is above MiniStringFragmentation threshold.
  if (legA == -1 || (mMin > mJoinJunction && massExcessIn > mStringMin))
    return false;

  // Construct separate index arrays for the three legs.
  vector<int> iLegA, iLegB, iLegC;
  leg = -1;
  for (int i = 0; i < int(iPartonIn.size()); ++ i) {
    if (iPartonIn[i] < 0) ++leg;
    else if( leg == legA) iLegA.push_back( iPartonIn[i] );
    else if( leg == legB) iLegB.push_back( iPartonIn[i] );
    else if( leg == legC) iLegC.push_back( iPartonIn[i] );
  }

  // First step: successively combine any gluons on the two legs.
  // (Presumably overkill; not likely to be (m)any extra gluons.)
  // (Do as successive binary joinings, so only need two mothers.)
  for (leg = 0; leg < 2; ++leg) {
    vector<int>& iLegNow = (leg == 0) ? iLegA : iLegB;
    int sizeNow = iLegNow.size();
    for (int i = sizeNow - 2; i >= 0; --i) {
      int iQ = iLegNow.back();
      int iG = iLegNow[i];
      int colNew = (event[iQ].id() > 0) ? event[iG].col() : 0;
      int acolNew = (event[iQ].id() < 0) ? event[iG].acol() : 0;
      Vec4 pNew = event[iQ].p() + event[iG].p();
      int iNew = event.append( event[iQ].id(), 74, iQ, iG, 0, 0,
        colNew, acolNew, pNew, pNew.mCalc() );

      // Mark joined partons and update iLeg end.
      event[iQ].statusNeg();
      event[iG].statusNeg();
      event[iQ].daughter1(iNew);
      event[iG].daughter1(iNew);
      iLegNow.back() = iNew;
    }
  }

  // Second step: combine two quarks into a diquark.
  int iQA     = iLegA.back();
  int iQB     = iLegB.back();
  int idQA    = event[iQA].id();
  int idQB    = event[iQB].id();
  int idNew   = flavSelPtr->makeDiquark( idQA, idQB );
  // Diquark colour is opposite to parton closest to junction on third leg.
  int colNew  = (idNew > 0) ? 0 : event[ iLegC[0] ].acol();
  int acolNew = (idNew > 0) ? event[ iLegC[0] ].col() : 0;
  Vec4 pNew   = pLeg[legA] + pLeg[legB];
  int iNew    = event.append( idNew, 74, min(iQA, iQB), max( iQA, iQB),
     0, 0, colNew, acolNew, pNew, pNew.mCalc() );

  // Mark joined partons and reduce remaining system.
  event[iQA].statusNeg();
  event[iQB].statusNeg();
  event[iQA].daughter1(iNew);
  event[iQB].daughter1(iNew);
  iPartonIn.resize(0);
  iPartonIn.push_back( iNew);
  for (int i = 0; i < int(iLegC.size()) ; ++i)
    iPartonIn.push_back( iLegC[i]);

  // Remove junction from event record list, identifying by colour.
  int iJun = -1;
  for (int i = 0; i < event.sizeJunction(); ++i)
    for (int j = 0; j < 3; ++ j)
      if ( event.colJunction(i,j) == max(colNew, acolNew) ) iJun = i;
  if (iJun >= 0) event.eraseJunction(iJun);

  // Done, having eliminated junction.
  return true;

}

//--------------------------------------------------------------------------

// Collect all partons of singlet to be consecutively ordered.

void ColConfig::collect(int iSub, Event& event, bool skipTrivial) {

  // Check that all partons have positive energy.
  for (int i = 0; i < singlets[iSub].size(); ++i) {
    int iNow = singlets[iSub].iParton[i];
    if (iNow > 0 && event[iNow].e() < 0.)
    infoPtr->errorMsg("Warning in ColConfig::collect: "
      "negative-energy parton encountered");
  }

  // Partons may already have been collected, e.g. at ministring collapse.
  if (singlets[iSub].isCollected) return;
  singlets[iSub].isCollected = true;

  // Check if partons already "by chance" happen to be ordered.
  bool inOrder = true;
  for (int i = 0; i < singlets[iSub].size() - 1; ++i) {
    int iFirst = singlets[iSub].iParton[i];
    if (iFirst < 0) continue;
    int iSecond = singlets[iSub].iParton[i + 1];
    if (iSecond < 0) iSecond = singlets[iSub].iParton[i + 2];
    if (iSecond != iFirst + 1) { inOrder = false; break;}
  }

  // Normally done if in order, but sometimes may need to copy anyway.
  if (inOrder && skipTrivial) return;

  // Copy down system. Update current partons.
  for (int i = 0; i < singlets[iSub].size(); ++i) {
    int iOld = singlets[iSub].iParton[i];
    if (iOld < 0) continue;
    int iNew;
    if (event[iOld].status() == 74) iNew = event.copy(iOld, 74);
    else iNew = event.copy(iOld, 71);
    singlets[iSub].iParton[i] = iNew;
  }

  // Done.
}

//--------------------------------------------------------------------------

// Find to which singlet system a particle belongs.

int ColConfig::findSinglet(int i) {

  // Loop through all systems and all members in them.
  for (int iSub = 0; iSub < int(singlets.size()); ++iSub)
  for (int iMem = 0; iMem < singlets[iSub].size(); ++iMem)
    if (singlets[iSub].iParton[iMem] == i) return iSub;

  // Done without having found particle; return -1 = error code.
  return -1;
}

//--------------------------------------------------------------------------

// List all currently identified singlets.

void ColConfig::list(ostream& os) const {

  // Header. Loop over all individual singlets.
  os << "\n --------  Colour Singlet Systems Listing -------------------\n";
  for (int iSub = 0; iSub < int(singlets.size()); ++iSub) {

    // List all partons belonging to each singlet.
    os << " singlet " << iSub << " contains " ;
    for (int i = 0; i < singlets[iSub].size(); ++i)
      os << singlets[iSub].iParton[i] << " ";
    os << "\n";

  // Done.
  }
}

//==========================================================================

// The StringRegion class.

// Currently a number of simplifications, in particular ??
// 1) No popcorn baryon production.
// 2) Simplified treatment of pT in stepping and joining.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// If a string region is smaller thsan this it is assumed empty.
const double StringRegion::MJOIN = 0.1;

// Avoid division by zero.
const double StringRegion::TINY  = 1e-20;

//--------------------------------------------------------------------------

// Set up four-vectors for longitudinal and transverse directions.

void StringRegion::setUp(Vec4 p1, Vec4 p2, bool isMassless) {

  // Simple case: the two incoming four-vectors guaranteed massless.
  if (isMassless) {

    // Calculate w2, minimum value. Lightcone directions = input.
    w2 = 2. * (p1 * p2);
    if (w2 < MJOIN*MJOIN) {isSetUp = true; isEmpty = true; return;}
    pPos = p1;
    pNeg = p2;

  // Else allow possibility of masses for incoming partons (also gluons!).
  } else {

    // Generic four-momentum combinations.
    double m1Sq = p1 * p1;
    double m2Sq = p2 * p2;
    double p1p2 = p1 * p2;
    w2 = m1Sq + 2. * p1p2 + m2Sq;
    double rootSq = pow2(p1p2) - m1Sq * m2Sq;

    // If crazy kinematics (should not happen!) modify energies.
    if (w2 <= 0. || rootSq <= 0.) {
      if (m1Sq < 0.) m1Sq = 0.;
      p1.e( sqrt(m1Sq + p1.pAbs2()) );
      if (m2Sq < 0.) m2Sq = 0.;
      p2.e( sqrt(m2Sq + p2.pAbs2()) );
      p1p2 = p1 * p2;
      w2 = m1Sq + 2. * p1p2 + m2Sq;
      rootSq = pow2(p1p2) - m1Sq * m2Sq;
    }

    // If still small invariant mass then empty region (e.g. in gg system).
    if (w2 < MJOIN*MJOIN) {isSetUp = true; isEmpty = true; return;}

    // Find two lightconelike longitudinal four-vector directions.
    double root = sqrt( max(TINY, rootSq) );
    double k1 = 0.5 * ( (m2Sq + p1p2) / root - 1.);
    double k2 = 0.5 * ( (m1Sq + p1p2) / root - 1.);
    pPos = (1. + k1) * p1 - k2 * p2;
    pNeg = (1. + k2) * p2 - k1 * p1;
  }

  // Find two spacelike transverse four-vector directions.
  // Begin by picking two sensible trial directions.
  Vec4 eDiff = pPos / pPos.e() - pNeg / pNeg.e();
  double eDx = pow2( eDiff.px() );
  double eDy = pow2( eDiff.py() );
  double eDz = pow2( eDiff.pz() );
  if (eDx < min(eDy, eDz)) {
    eX = Vec4( 1., 0., 0., 0.);
    eY = (eDy < eDz) ? Vec4( 0., 1., 0., 0.) : Vec4( 0., 0., 1., 0.);
  } else if (eDy < eDz) {
    eX = Vec4( 0., 1., 0., 0.);
    eY = (eDx < eDz) ? Vec4( 1., 0., 0., 0.) : Vec4( 0., 0., 1., 0.);
  } else {
    eX = Vec4( 0., 0., 1., 0.);
    eY = (eDx < eDy) ? Vec4( 1., 0., 0., 0.) : Vec4( 0., 1., 0., 0.);
  }

  // Then construct orthogonal linear combinations.
  double pPosNeg = pPos * pNeg;
  double kXPos = eX * pPos / pPosNeg;
  double kXNeg = eX * pNeg / pPosNeg;
  double kXX = 1. / sqrt( 1. + 2. * kXPos * kXNeg * pPosNeg );
  double kYPos = eY * pPos / pPosNeg;
  double kYNeg = eY * pNeg / pPosNeg;
  double kYX = kXX * (kXPos * kYNeg + kXNeg * kYPos) * pPosNeg;
  double kYY = 1. / sqrt(1. + 2. * kYPos * kYNeg * pPosNeg - pow2(kYX));
  eX = kXX * (eX - kXNeg * pPos - kXPos * pNeg);
  eY = kYY * (eY - kYNeg * pPos - kYPos * pNeg - kYX * eX);

  // Done.
  isSetUp = true;
  isEmpty = false;

}

//--------------------------------------------------------------------------

// Project a four-momentum onto (x+, x-, px, py).

void StringRegion::project(Vec4 pIn) {

  // Perform projections by four-vector multiplication.
  xPosProj = 2. * (pIn * pNeg) / w2;
  xNegProj = 2. * (pIn * pPos) / w2;
  pxProj = - (pIn * eX);
  pyProj = - (pIn * eY);

}

//==========================================================================

// The StringSystem class.

//--------------------------------------------------------------------------

// Set up system from parton list.

void StringSystem::setUp(vector<int>& iSys, Event& event) {

  // Figure out how big the system is. (Closed gluon loops?)
  sizePartons = iSys.size();
  sizeStrings = sizePartons - 1;
  sizeRegions = (sizeStrings * (sizeStrings + 1)) / 2;
  indxReg = 2 * sizeStrings + 1;
  iMax = sizeStrings - 1;

  // Reserve space for the required number of regions.
  system.clear();
  system.resize(sizeRegions);

  // Set up the lowest-lying regions.
  for (int i = 0; i < sizeStrings; ++i) {
    Vec4 p1 = event[ iSys[i] ].p();
    if ( event[ iSys[i] ].isGluon() ) p1 *= 0.5;
    Vec4 p2 = event[ iSys[i+1] ].p();
    if ( event[ iSys[i+1] ].isGluon() ) p2 *= 0.5;
    system[ iReg(i, iMax - i) ].setUp( p1, p2, false);
  }

}

//==========================================================================

} // end namespace Pythia8
