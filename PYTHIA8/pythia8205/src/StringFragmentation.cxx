// StringFragmentation.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the StringEnd and
// StringFragmentation classes.

#include "Pythia8/StringFragmentation.h"

namespace Pythia8 {

//==========================================================================

// The StringEnd class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.

// Avoid unphysical solutions to equation system.
const double StringEnd::TINY = 1e-6;

// Assume two (eX, eY) regions are related if pT2 differs by less.
const double StringEnd::PT2SAME = 0.01;

//--------------------------------------------------------------------------

// Set up initial endpoint values from input.

void StringEnd::setUp(bool fromPosIn, int iEndIn, int idOldIn, int iMaxIn,
  double pxIn, double pyIn, double GammaIn, double xPosIn, double xNegIn) {

  // Simple transcription from input.
  fromPos  = fromPosIn;
  iEnd     = iEndIn;
  iMax     = iMaxIn;
  flavOld  = FlavContainer(idOldIn);
  pxOld    = pxIn;
  pyOld    = pyIn;
  GammaOld = GammaIn;
  iPosOld  = (fromPos) ? 0 : iMax;
  iNegOld  = (fromPos) ? iMax : 0;
  xPosOld  = xPosIn;
  xNegOld  = xNegIn;

}

//--------------------------------------------------------------------------

// Fragment off one hadron from the string system, in flavour and pT.

void StringEnd::newHadron() {

  // Pick new flavour and form a new hadron.
  do {
    flavNew = flavSelPtr->pick( flavOld);
    idHad   = flavSelPtr->combine( flavOld, flavNew);
  } while (idHad == 0);

  // Pick its transverse momentum.
  pair<double, double> pxy = pTSelPtr->pxy();
  pxNew = pxy.first;
  pyNew = pxy.second;
  pxHad = pxOld + pxNew;
  pyHad = pyOld + pyNew;

  // Pick its mass and thereby define its transverse mass.
  mHad   = particleDataPtr->mSel(idHad);
  mT2Had = pow2(mHad) + pow2(pxHad) + pow2(pyHad);

}

//--------------------------------------------------------------------------

// Fragment off one hadron from the string system, in momentum space,
// by taking steps from positive end.

Vec4 StringEnd::kinematicsHadron( StringSystem& system) {

  // Pick fragmentation step z and calculate new Gamma.
  zHad = zSelPtr->zFrag( flavOld.id, flavNew.id, mT2Had);
  GammaNew = (1. - zHad) * (GammaOld + mT2Had / zHad);

  // Set up references that are direction-neutral;
  // ...Dir for direction of iteration and ...Inv for its inverse.
  int&    iDirOld = (fromPos) ? iPosOld : iNegOld;
  int&    iInvOld = (fromPos) ? iNegOld : iPosOld;
  int&    iDirNew = (fromPos) ? iPosNew : iNegNew;
  int&    iInvNew = (fromPos) ? iNegNew : iPosNew;
  double& xDirOld = (fromPos) ? xPosOld : xNegOld;
  double& xInvOld = (fromPos) ? xNegOld : xPosOld;
  double& xDirNew = (fromPos) ? xPosNew : xNegNew;
  double& xInvNew = (fromPos) ? xNegNew : xPosNew;
  double& xDirHad = (fromPos) ? xPosHad : xNegHad;
  double& xInvHad = (fromPos) ? xNegHad : xPosHad;

  // Start search for new breakup in the old region.
  iDirNew = iDirOld;
  iInvNew = iInvOld;
  Vec4 pTNew;

  // Each step corresponds to trying a new string region.
  for (int iStep = 0; ; ++iStep) {

    // Referance to current string region.
    StringRegion& region = system.region( iPosNew, iNegNew);

    // Now begin special section for rapid processing of low region.
    if (iStep == 0 && iPosOld + iNegOld == iMax) {

      // A first step within a low region is easy.
      if (mT2Had < zHad * xDirOld * (1. - xInvOld) * region.w2) {

        // Translate into x coordinates.
        xDirHad = zHad * xDirOld;
        xInvHad = mT2Had / (xDirHad * region.w2);
        xDirNew = xDirOld - xDirHad;
        xInvNew = xInvOld + xInvHad;

        // Find and return four-momentum of the produced particle.
        return region.pHad( xPosHad, xNegHad, pxHad, pyHad);

      // A first step out of a low region also OK, if there are more regions.
      // Negative energy signals failure, i.e. in last region.
      } else {
        --iInvNew;
        if (iInvNew < 0) return Vec4(0., 0., 0., -1.);

        // Momentum taken by stepping out of region. Continue to next region.
        xInvHad = 1. - xInvOld;
        xDirHad = 0.;
        pSoFar  = region.pHad( xPosHad, xNegHad, pxOld, pyOld);
        continue;
      }

    // Else, for first step, take into account starting pT.
    } else if (iStep == 0) {
      pSoFar = region.pHad( 0., 0., pxOld, pyOld);
      pTNew  = region.pHad( 0., 0., pxNew, pyNew);
    }

    // Now begin normal treatment of nontrivial regions.
    // Set up four-vectors in a region not visited before.
    if (!region.isSetUp) region.setUp(
      system.regionLowPos(iPosNew).pPos,
      system.regionLowNeg(iNegNew).pNeg, true);

    // If new region is vanishingly small, continue immediately to next.
    // Negative energy signals failure to do this, i.e. moved too low.
    if (region.isEmpty) {
      xDirHad = (iDirNew == iDirOld) ? xDirOld : 1.;
      xInvHad = 0.;
      pSoFar += region.pHad( xPosHad, xNegHad, 0., 0.);
      ++iDirNew;
      if (iDirNew + iInvNew > iMax) return Vec4(0., 0., 0., -1.);
      continue;
    }

    // Reexpress pTNew w.r.t. base vectors in new region, if possible.
    // Recall minus sign from normalization e_x^2 = e_y^2 = -1.
    double pxNewTemp = -pTNew * region.eX;
    double pyNewTemp = -pTNew * region.eY;
    if (abs( pxNewTemp * pxNewTemp + pyNewTemp * pyNewTemp
      - pxNew * pxNew - pyNew * pyNew) < PT2SAME) {
      pxNew = pxNewTemp;
      pyNew = pyNewTemp;
    }

    // Four-momentum taken so far, including new pT.
    Vec4 pTemp = pSoFar + region.pHad( 0., 0., pxNew, pyNew);

    // Derive coefficients for m2 expression.
    // cM2 * x+ + cM3 * x- + cM4 * x+ * x- = m^2 - cM1;
    double cM1 = pTemp.m2Calc();
    double cM2 = 2. * (pTemp * region.pPos);
    double cM3 = 2. * (pTemp * region.pNeg);
    double cM4 = region.w2;
    if (!fromPos) swap( cM2, cM3);

    // Derive coefficients for Gamma expression.
    // cGam2 * x+ + cGam3 * x- + cGam4 * x+ * x- = Gamma_new - cGam1;
    double cGam1 = 0.;
    double cGam2 = 0.;
    double cGam3 = 0.;
    double cGam4 = 0.;
    for (int iInv = iInvNew; iInv <= iMax - iDirNew; ++iInv) {
      double xInv = 1.;
      if (iInv == iInvNew) xInv = (iInvNew == iInvOld) ? xInvOld : 0.;
      for (int iDir = iDirNew; iDir <= iMax - iInv; ++iDir) {
        double xDir = (iDir == iDirOld) ? xDirOld : 1.;
        int iPos = (fromPos) ? iDir : iInv;
        int iNeg = (fromPos) ? iInv : iDir;
        StringRegion& regionGam =  system.region( iPos, iNeg);
        if (!regionGam.isSetUp) regionGam.setUp(
          system.regionLowPos(iPos).pPos,
          system.regionLowNeg(iNeg).pNeg, true);
        double w2 = regionGam.w2;
        cGam1 += xDir * xInv * w2;
        if (iDir == iDirNew) cGam2 -= xInv * w2;
        if (iInv == iInvNew) cGam3 += xDir * w2;
        if (iDir == iDirNew && iInv == iInvNew) cGam4 -= w2;
      }
    }

    // Solve (m2, Gamma) equation system => r2 * x-^2 + r1 * x- + r0 = 0.
    double cM0   = pow2(mHad) - cM1;
    double cGam0 = GammaNew - cGam1;
    double r2    = cM3 * cGam4 - cM4 * cGam3;
    double r1    = cM4 * cGam0 - cM0 * cGam4 + cM3 * cGam2 - cM2 * cGam3;
    double r0    = cM2 * cGam0 - cM0 * cGam2;
    double root  = sqrtpos( r1*r1 - 4. * r2 * r0 );
    if (abs(r2) < TINY || root < TINY) return Vec4(0., 0., 0., -1.);
    xInvHad      = 0.5 * (root / abs(r2) - r1 / r2);
    xDirHad      = (cM0 - cM3 * xInvHad) / (cM2 + cM4 * xInvHad);

    // Define position of new trial vertex.
    xDirNew = (iDirNew == iDirOld) ? xDirOld - xDirHad : 1. - xDirHad;
    xInvNew = (iInvNew == iInvOld) ? xInvOld + xInvHad : xInvHad;

    // Step up to new region if new x- > 1.
    if (xInvNew > 1.) {
      xInvHad = (iInvNew == iInvOld) ? 1. - xInvOld : 1.;
      xDirHad = 0.;
      pSoFar += region.pHad( xPosHad, xNegHad, 0., 0.);
      --iInvNew;
      if (iInvNew < 0) return Vec4(0., 0., 0., -1.);
      continue;

    // Step down to new region if new x+ < 0.
    } else if (xDirNew < 0.) {
      xDirHad = (iDirNew == iDirOld) ? xDirOld : 1.;
      xInvHad = 0.;
      pSoFar += region.pHad( xPosHad, xNegHad, 0., 0.);
      ++iDirNew;
      if (iDirNew + iInvNew > iMax) return Vec4(0., 0., 0., -1.);
      continue;
    }

    // Else we have found the correct region, and can return four-momentum.
    return pSoFar + region.pHad( xPosHad, xNegHad, pxNew, pyNew);

  // End of "infinite" loop of stepping to new region.
  }

}

//--------------------------------------------------------------------------

// Update string end information after a hadron has been removed.

void StringEnd::update() {

  flavOld.anti(flavNew);
  iPosOld  = iPosNew;
  iNegOld  = iNegNew;
  pxOld    = -pxNew;
  pyOld    = -pyNew;
  GammaOld = GammaNew;
  xPosOld  = xPosNew;
  xNegOld  = xNegNew;

}

//==========================================================================

// The StringFragmentation class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to (flavour-, energy) join the two string ends.
const int    StringFragmentation::NTRYFLAV      = 10;
const int    StringFragmentation::NTRYJOIN      = 30;

// The last few times gradually increase the stop mass to make it easier.
const int    StringFragmentation::NSTOPMASS     = 15;
const double StringFragmentation::FACSTOPMASS   = 1.05;

// For closed string, pick a Gamma by taking a step with fictitious mass.
const double StringFragmentation::CLOSEDM2MAX   = 25.;
const double StringFragmentation::CLOSEDM2FRAC  = 0.1;

// Do not allow too large argument to exp function.
const double StringFragmentation::EXPMAX        = 50.;

// Matching criterion that p+ and p- not the same (can happen in gg loop).
const double StringFragmentation::MATCHPOSNEG   = 1e-4;

// For pull on junction, do not trace too far down each leg.
const double StringFragmentation::EJNWEIGHTMAX  = 10.;

// Iterate junction rest frame boost until convergence or too many tries.
const double StringFragmentation::CONVJNREST   = 1e-5;
const int    StringFragmentation::NTRYJNREST   = 20;

// Fail and try again when two legs combined to diquark (3 loops).
const int    StringFragmentation::NTRYJNMATCH   = 20;
const double StringFragmentation::EEXTRAJNMATCH = 0.5;
const double StringFragmentation::MDIQUARKMIN   = -2.0;

// Consider junction-leg parton as massless if m2 tiny.
const double StringFragmentation::M2MAXJRF      = 1e-4;

// Iterate junction rest frame equation until convergence or too many tries.
const double StringFragmentation::CONVJRFEQ     = 1e-12;
const int    StringFragmentation::NTRYJRFEQ     = 40;

//--------------------------------------------------------------------------

// Initialize and save pointers.

void StringFragmentation::init(Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn, StringFlav* flavSelPtrIn,
  StringPT* pTSelPtrIn, StringZ* zSelPtrIn) {

  // Save pointers.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  flavSelPtr      = flavSelPtrIn;
  pTSelPtr        = pTSelPtrIn;
  zSelPtr         = zSelPtrIn;

  // Initialize the StringFragmentation class.
  stopMass        = zSelPtr->stopMass();
  stopNewFlav     = zSelPtr->stopNewFlav();
  stopSmear       = zSelPtr->stopSmear();
  eNormJunction   = settings.parm("StringFragmentation:eNormJunction");
  eBothLeftJunction
     = settings.parm("StringFragmentation:eBothLeftJunction");
  eMaxLeftJunction
    = settings.parm("StringFragmentation:eMaxLeftJunction");
  eMinLeftJunction
    = settings.parm("StringFragmentation:eMinLeftJunction");

  // Joining of nearby partons along the string.
  mJoin           = settings.parm("FragmentationSystems:mJoin");

  // Initialize the b parameter of the z spectrum, used when joining jets.
  bLund           = zSelPtr->bAreaLund();

  // Initialize the hadrons instance of an event record.
  hadrons.init( "(string fragmentation)", particleDataPtr);

  // Send on pointers to the two StringEnd instances.
  posEnd.init( particleDataPtr, flavSelPtr, pTSelPtr, zSelPtr);
  negEnd.init( particleDataPtr, flavSelPtr, pTSelPtr, zSelPtr);

}

//--------------------------------------------------------------------------

// Perform the fragmentation.

bool StringFragmentation::fragment( int iSub, ColConfig& colConfig,
  Event& event) {

  // Find partons and their total four-momentum.
  iParton            = colConfig[iSub].iParton;
  iPos               = iParton[0];
  if (iPos < 0) iPos = iParton[1];
  int idPos          = event[iPos].id();
  iNeg               = iParton.back();
  int idNeg          = event[iNeg].id();
  pSum               = colConfig[iSub].pSum;

  // Reset the local event record.
  hadrons.clear();

  // For closed gluon string: pick first breakup region.
  isClosed = colConfig[iSub].isClosed;
  if (isClosed) iParton = findFirstRegion(iParton, event);

  // For junction topology: fragment off two of the string legs.
  // Then iParton overwritten to remaining leg + leftover diquark.
  pJunctionHadrons = 0.;
  hasJunction = colConfig[iSub].hasJunction;
  if (hasJunction && !fragmentToJunction(event)) return false;
  int junctionHadrons = hadrons.size();
  if (hasJunction) {
    idPos  = event[ iParton[0] ].id();
    idNeg  = event.back().id();
    pSum  -= pJunctionHadrons;
  }

  // Set up kinematics of string evolution ( = motion).
  system.setUp(iParton, event);
  stopMassNow = stopMass;
  int nExtraJoin = 0;

  // Fallback loop, when joining in the middle fails.  Bailout if stuck.
  for ( int iTry = 0; ; ++iTry) {
    if (iTry > NTRYJOIN) {
      infoPtr->errorMsg("Error in StringFragmentation::fragment: "
        "stuck in joining");
      if (hasJunction) ++nExtraJoin;
      if (nExtraJoin > 0) event.popBack(nExtraJoin);
      return false;
    }

    // After several failed tries join some (extra) nearby partons.
    if (iTry == NTRYJOIN / 3) nExtraJoin = extraJoin( 2., event);
    if (iTry == 2 * NTRYJOIN / 3) nExtraJoin += extraJoin( 4., event);

    // After several failed tries gradually allow larger stop mass.
    if (iTry > NTRYJOIN - NSTOPMASS) stopMassNow *= FACSTOPMASS;

    // Set up flavours of two string ends, and reset other info.
    setStartEnds(idPos, idNeg, system);
    pRem = pSum;

    // Begin fragmentation loop, interleaved from the two ends.
    bool fromPos;

    // Variables used to help identifying baryons from junction splittings.
    bool usedPosJun = false, usedNegJun = false;

    for ( ; ; ) {

      // Take a step either from the positive or the negative end.
      fromPos           = (rndmPtr->flat() < 0.5);
      StringEnd& nowEnd = (fromPos) ? posEnd : negEnd;

      // Construct trial hadron and check that energy remains.
      nowEnd.newHadron();
      if ( energyUsedUp(fromPos) ) break;

      // Construct kinematics of the new hadron and store it.
      Vec4 pHad = nowEnd.kinematicsHadron(system);
      int statusHad = (fromPos) ? 83 : 84;

      // Change status code if hadron from junction.
      if (abs(nowEnd.idHad) > 1000 && abs(nowEnd.idHad) < 10000) {
        if (fromPos && event[iPos].statusAbs() == 74 && !usedPosJun)  {
          statusHad = 87;
          usedPosJun = true;
        }
        if (!fromPos && event[iNeg].statusAbs() == 74 && !usedNegJun)  {
          statusHad = 88;
          usedNegJun = true;
        }
        if (!fromPos && hasJunction && !usedNegJun) {
          statusHad = 88;
          usedNegJun = true;
        }
      }
      hadrons.append( nowEnd.idHad, statusHad, iPos, iNeg,
        0, 0, 0, 0, pHad, nowEnd.mHad);
      if (pHad.e() < 0.) break;

      // Update string end and remaining momentum.
      nowEnd.update();
      pRem -= pHad;

    // End of fragmentation loop.
    }

    // When done, join in the middle. If this works, then really done.
    if ( finalTwo(fromPos, event, usedPosJun, usedNegJun) ) break;

    // Else remove produced particles (except from first two junction legs)
    // and start all over.
    int newHadrons = hadrons.size() - junctionHadrons;
    hadrons.popBack(newHadrons);
  }


  // Junctions & extra joins: remove fictitious end, restore original partons.
  if (hasJunction) ++nExtraJoin;
  if (nExtraJoin > 0) {
    event.popBack(nExtraJoin);
    iParton = colConfig[iSub].iParton;
  }

  // Store the hadrons in the normal event record, ordered from one end.
  store(event);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Find region where to put first string break for closed gluon loop.

vector<int> StringFragmentation::findFirstRegion(vector<int>& iPartonIn,
  Event& event) {

  // Evaluate mass-squared for all adjacent gluon pairs.
  vector<double> m2Pair;
  double m2Sum = 0.;
  int size = iPartonIn.size();
  for (int i = 0; i < size; ++i) {
    double m2Now = 0.5 * event[ iPartonIn[i] ].p()
      * event[ iPartonIn[(i + 1)%size] ].p();
    m2Pair.push_back(m2Now);
    m2Sum += m2Now;
  }

  // Pick breakup region with probability proportional to mass-squared.
  double m2Reg = m2Sum * rndmPtr->flat();
  int iReg = -1;
  do m2Reg -= m2Pair[++iReg];
  while (m2Reg > 0. && iReg < size - 1);

  // Create reordered parton list, with breakup string region duplicated.
  vector<int> iPartonOut;
  for (int i = 0; i < size + 2; ++i)
    iPartonOut.push_back( iPartonIn[(i + iReg)%size] );

  // Done.
  return iPartonOut;

}

//--------------------------------------------------------------------------

// Set flavours and momentum position for initial string endpoints.

void StringFragmentation::setStartEnds( int idPos, int idNeg,
  StringSystem systemNow) {

  // Variables characterizing string endpoints: defaults for open string.
  double px          = 0.;
  double py          = 0.;
  double Gamma       = 0.;
  double xPosFromPos = 1.;
  double xNegFromPos = 0.;
  double xPosFromNeg = 0.;
  double xNegFromNeg = 1.;

  // For closed gluon loop need to pick an initial flavour.
  if (isClosed) {
    do {
      int idTry = flavSelPtr->pickLightQ();
      FlavContainer flavTry(idTry, 1);
      flavTry = flavSelPtr->pick( flavTry);
      flavTry = flavSelPtr->pick( flavTry);
      idPos   = flavTry.id;
      idNeg   = -idPos;
    } while (idPos == 0);

    // Also need pT and breakup vertex position in region.
   pair<double, double> pxy = pTSelPtr->pxy();
   px = pxy.first;
   py = pxy.second;
    double m2Region = systemNow.regionLowPos(0).w2;
    double m2Temp   = min( CLOSEDM2MAX, CLOSEDM2FRAC * m2Region);
    do {
      double zTemp = zSelPtr->zFrag( idPos, idNeg, m2Temp);
      xPosFromPos  = 1. - zTemp;
      xNegFromPos  = m2Temp / (zTemp * m2Region);
    } while (xNegFromPos > 1.);
    Gamma = xPosFromPos * xNegFromPos * m2Region;
    xPosFromNeg = xPosFromPos;
    xNegFromNeg = xNegFromPos;
  }

  // Initialize two string endpoints.
  posEnd.setUp(  true, iPos, idPos, systemNow.iMax,  px,  py,
    Gamma, xPosFromPos, xNegFromPos);
  negEnd.setUp( false, iNeg, idNeg, systemNow.iMax, -px, -py,
    Gamma, xPosFromNeg, xNegFromNeg);

  // For closed gluon loop can allow popcorn on one side but not both.
  if (isClosed) {
    flavSelPtr->assignPopQ(posEnd.flavOld);
    flavSelPtr->assignPopQ(negEnd.flavOld);
    if (rndmPtr->flat() < 0.5) posEnd.flavOld.nPop = 0;
    else                    negEnd.flavOld.nPop = 0;
    posEnd.flavOld.rank = 1;
    negEnd.flavOld.rank = 1;
  }

  // Done.

}

//--------------------------------------------------------------------------

// Check remaining energy-momentum whether it is OK to continue.

bool StringFragmentation::energyUsedUp(bool fromPos) {

  // If remaining negative energy then abort right away.
  if (pRem.e() < 0.) return true;

  // Calculate W2_minimum and done if remaining W2 is below it.
  double wMin = stopMassNow
    + particleDataPtr->constituentMass(posEnd.flavOld.id)
    + particleDataPtr->constituentMass(negEnd.flavOld.id);
  if (fromPos) wMin += stopNewFlav
    * particleDataPtr->constituentMass(posEnd.flavNew.id);
  else         wMin += stopNewFlav
    * particleDataPtr->constituentMass(negEnd.flavNew.id);
  wMin *= 1. + (2. * rndmPtr->flat() - 1.) * stopSmear;
  w2Rem = pRem.m2Calc();
  if (w2Rem < pow2(wMin)) return true;

  // Else still enough energy left to continue iteration.
  return false;

}

//--------------------------------------------------------------------------

// Produce the final two partons to complete the system.

bool StringFragmentation::finalTwo(bool fromPos, Event& event, bool usedPosJun,
  bool usedNegJun) {

  // Check whether we went too far in p+-.
  if (pRem.e() < 0.  || w2Rem < 0. || (hadrons.size() > 0
    && hadrons.back().e() < 0.) ) return false;
  if ( posEnd.iPosOld > negEnd.iPosOld || negEnd.iNegOld > posEnd.iNegOld)
    return false;
  if ( posEnd.iPosOld == negEnd.iPosOld && posEnd.xPosOld < negEnd.xPosOld)
    return false;
  if ( posEnd.iNegOld == negEnd.iNegOld && posEnd.xNegOld > negEnd.xNegOld)
    return false;

  // Construct the final hadron from the leftover flavours.
  // Impossible to join two diquarks. Also break if stuck for other reason.
  FlavContainer flav1 = (fromPos) ? posEnd.flavNew.anti() : posEnd.flavOld;
  FlavContainer flav2 = (fromPos) ? negEnd.flavOld : negEnd.flavNew.anti();
  if (flav1.isDiquark() && flav2.isDiquark()) return false;
  int idHad;
  for (int iTry = 0; iTry < NTRYFLAV; ++iTry) {
    idHad = flavSelPtr->combine( flav1, flav2);
    if (idHad != 0) break;
  }
  if (idHad == 0) return false;

  // Store the final particle and its new pT, and construct its mass.
  if (fromPos) {
    negEnd.idHad = idHad;
    negEnd.pxNew = -posEnd.pxNew;
    negEnd.pyNew = -posEnd.pyNew;
    negEnd.mHad  = particleDataPtr->mSel(idHad);
  } else {
    posEnd.idHad = idHad;
    posEnd.pxNew = -negEnd.pxNew;
    posEnd.pyNew = -negEnd.pyNew;
    posEnd.mHad  = particleDataPtr->mSel(idHad);
  }

  // String region in which to do the joining.
  StringRegion region = finalRegion();
  if (region.isEmpty) return false;

  // Project remaining momentum along longitudinal and transverse directions.
  region.project( pRem);
  double pxRem   = region.px() - posEnd.pxOld - negEnd.pxOld;
  double pyRem   = region.py() - posEnd.pyOld - negEnd.pyOld;
  double xPosRem = region.xPos();
  double xNegRem = region.xNeg();

  // Share extra pT kick evenly between final two hadrons.
  posEnd.pxOld += 0.5 * pxRem;
  posEnd.pyOld += 0.5 * pyRem;
  negEnd.pxOld += 0.5 * pxRem;
  negEnd.pyOld += 0.5 * pyRem;

  // Construct new pT and mT of the final two particles.
  posEnd.pxHad  = posEnd.pxOld + posEnd.pxNew;
  posEnd.pyHad  = posEnd.pyOld + posEnd.pyNew;
  posEnd.mT2Had = pow2(posEnd.mHad) + pow2(posEnd.pxHad)
    + pow2(posEnd.pyHad);
  negEnd.pxHad  = negEnd.pxOld + negEnd.pxNew;
  negEnd.pyHad  = negEnd.pyOld + negEnd.pyNew;
  negEnd.mT2Had = pow2(negEnd.mHad) + pow2(negEnd.pxHad)
    + pow2(negEnd.pyHad);

  // Construct remaining system transverse mass.
  double wT2Rem = w2Rem + pow2( posEnd.pxHad + negEnd.pxHad)
    + pow2( posEnd.pyHad + negEnd.pyHad);

  // Check that kinematics possible.
  if ( sqrt(wT2Rem) < sqrt(posEnd.mT2Had) + sqrt(negEnd.mT2Had) )
    return false;
  double lambda2 = pow2( wT2Rem - posEnd.mT2Had - negEnd.mT2Had)
    - 4. * posEnd.mT2Had * negEnd.mT2Had;
  if (lambda2 <= 0.) return false;

  // Construct kinematics, as viewed in the transverse rest frame.
  double lambda = sqrt(lambda2);
  double probReverse = 1. / (1. + exp( min( EXPMAX, bLund * lambda) ) );
  double xpzPos = 0.5 * lambda/ wT2Rem;
  if (probReverse > rndmPtr->flat()) xpzPos = -xpzPos;
  double xmDiff = (posEnd.mT2Had - negEnd.mT2Had) / wT2Rem;
  double xePos  = 0.5 * (1. + xmDiff);
  double xeNeg  = 0.5 * (1. - xmDiff );

  // Translate this into kinematics in the string frame.
  Vec4 pHadPos = region.pHad( (xePos + xpzPos) *  xPosRem,
    (xePos - xpzPos) *  xNegRem, posEnd.pxHad, posEnd.pyHad);
  Vec4 pHadNeg = region.pHad( (xeNeg - xpzPos) *  xPosRem,
    (xeNeg + xpzPos) *  xNegRem, negEnd.pxHad, negEnd.pyHad);

  // Update status codes for junction baryons.
  int statusHadPos = 83;
  int statusHadNeg = 84;

  if (fromPos) {
    if (abs(posEnd.idHad) > 1000 && abs(posEnd.idHad) < 10000) {
      if (event[iPos].statusAbs() == 74 && !usedPosJun)  {
        statusHadPos = 87;
        usedPosJun = true;
      }
    }
    if (abs(idHad) > 1000 && abs(idHad) < 10000) {
      if ((!usedNegJun && (event[iNeg].statusAbs() == 74 || hasJunction))
          || (!usedPosJun && event[iPos].statusAbs() == 74)) {
        statusHadNeg = 88;
      }
    }
  } else {
    if (abs(negEnd.idHad) > 1000 && abs(negEnd.idHad) < 10000) {
        if (!usedNegJun && (event[iNeg].statusAbs() == 74 || hasJunction)) {
          statusHadNeg = 87;
          usedNegJun = true;
        }
      }
    if (abs(idHad) > 1000 && abs(idHad) < 10000) {
      if ((!usedNegJun && (event[iNeg].statusAbs() == 74 || hasJunction))
          || (!usedPosJun && event[iPos].statusAbs() == 74)) {
        statusHadPos = 88;
      }
    }
  }

  // Add produced particles to the event record.
  hadrons.append( posEnd.idHad, statusHadPos, posEnd.iEnd, negEnd.iEnd,
    0, 0, 0, 0, pHadPos, posEnd.mHad);
  hadrons.append( negEnd.idHad, statusHadNeg, posEnd.iEnd, negEnd.iEnd,
    0, 0, 0, 0, pHadNeg, negEnd.mHad);

  // It worked.
  return true;

}

//--------------------------------------------------------------------------

//  Construct a special joining region for the final two hadrons.

StringRegion StringFragmentation::finalRegion() {

  // Simple case when both string ends are in the same region.
  if (posEnd.iPosOld == negEnd.iPosOld && posEnd.iNegOld == negEnd.iNegOld)
    return system.region( posEnd.iPosOld, posEnd.iNegOld);

  // Start out with empty region. (Empty used for error returns.)
  StringRegion region;

  // Add up all remaining p+.
  Vec4 pPosJoin;
  if ( posEnd.iPosOld == negEnd.iPosOld) {
    double xPosJoin = posEnd.xPosOld - negEnd.xPosOld;
    if (xPosJoin < 0.) return region;
    pPosJoin = system.regionLowPos(posEnd.iPosOld).pHad( xPosJoin, 0., 0., 0.);
  } else {
    for (int iPosNow = posEnd.iPosOld; iPosNow <= negEnd.iPosOld; ++iPosNow) {
      if (iPosNow == posEnd.iPosOld) pPosJoin
        += system.regionLowPos(iPosNow).pHad( posEnd.xPosOld, 0., 0., 0.);
      else if (iPosNow == negEnd.iPosOld) pPosJoin
        += system.regionLowPos(iPosNow).pHad( 1. - negEnd.xPosOld, 0., 0., 0.);
      else pPosJoin += system.regionLowPos(iPosNow).pHad( 1., 0., 0., 0.);
    }
  }

  // Add up all remaining p-.
  Vec4 pNegJoin;
  if ( negEnd.iNegOld == posEnd.iNegOld) {
    double xNegJoin = negEnd.xNegOld - posEnd.xNegOld;
    if (xNegJoin < 0.) return region;
    pNegJoin = system.regionLowNeg(negEnd.iNegOld).pHad( 0., xNegJoin, 0., 0.);
  } else {
    for (int iNegNow = negEnd.iNegOld; iNegNow <= posEnd.iNegOld; ++iNegNow) {
      if (iNegNow == negEnd.iNegOld) pNegJoin
        += system.regionLowNeg(iNegNow).pHad( 0., negEnd.xNegOld, 0., 0.);
      else if (iNegNow == posEnd.iNegOld) pNegJoin
        += system.regionLowNeg(iNegNow).pHad( 0., 1. - posEnd.xNegOld, 0., 0.);
      else pNegJoin += system.regionLowNeg(iNegNow).pHad( 0., 1., 0., 0.);
    }
  }

  // For a closed gluon loop pPosJoin == pNegJoin and the above does not work.
  // So reshuffle; "perfect" for g g systems, OK in general.
  Vec4 pTest = pPosJoin - pNegJoin;
  if ( abs(pTest.px()) + abs(pTest.py()) + abs(pTest.pz()) + abs(pTest.e())
    < MATCHPOSNEG * (pPosJoin.e() + pNegJoin.e()) ) {
    Vec4 delta
      = system.regionLowPos(posEnd.iPosOld + 1).pHad( 1., 0., 0., 0.)
      - system.regionLowNeg(negEnd.iNegOld + 1).pHad( 0., 1., 0., 0.);
    // If reshuffle did not help then pick random axis to break tie.
    // (Needed for low-mass q-g-qbar with q-qbar perfectly parallel.)
    if ( abs(delta.px()) + abs(delta.py()) + abs(delta.pz()) + abs(delta.e())
      < MATCHPOSNEG * (pPosJoin.e() + pNegJoin.e()) ) {
      double cthe = 2. * rndmPtr->flat() - 1.;
      double sthe = sqrtpos(1. - cthe * cthe);
      double phi  = 2. * M_PI * rndmPtr->flat();
      delta = 0.5 * min( pPosJoin.e(), pNegJoin.e()) 
        * Vec4( sthe * sin(phi), sthe * cos(phi), cthe, 0.);
      infoPtr->errorMsg("Warning in StringFragmentation::finalRegion: "
        "random axis needed to break tie");
    }  
    pPosJoin -= delta;
    pNegJoin += delta;
  }

  // Construct a new region from remaining p+ and p-.
  region.setUp( pPosJoin, pNegJoin);
  if (region.isEmpty) return region;

  // Project the existing pTold vectors onto the new directions.
  Vec4 pTposOld = system.region( posEnd.iPosOld, posEnd.iNegOld).pHad(
    0., 0.,  posEnd.pxOld, posEnd.pyOld);
  region.project( pTposOld);
  posEnd.pxOld = region.px();
  posEnd.pyOld = region.py();
  Vec4 pTnegOld = system.region( negEnd.iPosOld, negEnd.iNegOld).pHad(
      0., 0.,  negEnd.pxOld, negEnd.pyOld);
  region.project( pTnegOld);
  negEnd.pxOld = region.px();
  negEnd.pyOld = region.py();

  // Done.
  return region;

}

//--------------------------------------------------------------------------

// Store the hadrons in the normal event record, ordered from one end.

void StringFragmentation::store(Event& event) {

  // Starting position.
  int iFirst = event.size();

  // Copy straight over from first two junction legs.
  if (hasJunction) {
    for (int i = 0; i < hadrons.size(); ++i)
    if (hadrons[i].status() == 85 || hadrons[i].status() == 86)
      event.append( hadrons[i] );
  }

  // Loop downwards, copying all from the positive end.
  for (int i = 0; i < hadrons.size(); ++i)
    if (hadrons[i].status() == 83 || hadrons[i].status() == 87)
      event.append( hadrons[i] );

  // Loop upwards, copying all from the negative end.
  for (int i = hadrons.size() - 1; i >= 0 ; --i)
    if (hadrons[i].status() == 84 || hadrons[i].status() == 88)
      event.append( hadrons[i] );

  int iLast = event.size() - 1;

  // Set decay vertex when this is displaced.
  if (event[posEnd.iEnd].hasVertex()) {
    Vec4 vDec = event[posEnd.iEnd].vDec();
    for (int i = iFirst; i <= iLast; ++i) event[i].vProd( vDec );
  }

  // Set lifetime of hadrons.
  for (int i = iFirst; i <= iLast; ++i)
    event[i].tau( event[i].tau0() * rndmPtr->exp() );

  // Mark original partons as hadronized and set their daughter range.
  for (int i = 0; i < int(iParton.size()); ++i)
  if (iParton[i] >= 0) {
    event[ iParton[i] ].statusNeg();
    event[ iParton[i] ].daughters(iFirst, iLast);
  }

}

//--------------------------------------------------------------------------

// Fragment off two of the string legs in to a junction.

bool StringFragmentation::fragmentToJunction(Event& event) {

  // Identify range of partons on the three legs.
  // (Each leg begins with an iParton[i] = -(10 + 10*junctionNumber + leg),
  // and partons then appear ordered from the junction outwards.)
  int legBeg[3] = { 0, 0, 0};
  int legEnd[3] = { 0, 0, 0};
  int leg = -1;
  // PS (4/10/2011) Protect against invalid systems
  if (iParton[0] > 0) {
    infoPtr->errorMsg("Error in StringFragmentation::fragment"
                      "ToJunction: iParton[0] not a valid junctionNumber");
    return false;
  }
  for (int i = 0; i < int(iParton.size()); ++i) {
    if (iParton[i] < 0) {
      if (leg == 2) {
        infoPtr->errorMsg("Error in StringFragmentation::fragment"
                          "ToJunction: unprocessed multi-junction system");
        return false;
      }
      legBeg[++leg] = i + 1;
    }
    else legEnd[leg] = i;
  }

  // Iterate from system rest frame towards the junction rest frame (JRF).
  RotBstMatrix MtoJRF, Mstep;
  MtoJRF.bstback(pSum);
  Vec4 pWTinJRF[3];
  int iter = 0;
  double errInCM = 0.;
  do {
    ++iter;

    // Find weighted sum of momenta on the three sides of the junction.
    for (leg = 0; leg < 3; ++ leg) {
      pWTinJRF[leg] = 0.;
      double eWeight = 0.;
      for (int i = legBeg[leg]; i <= legEnd[leg]; ++i) {
        Vec4 pTemp = event[ iParton[i] ].p();
        pTemp.rotbst(MtoJRF);
        pWTinJRF[leg] += pTemp * exp(-eWeight);
        eWeight += pTemp.e() / eNormJunction;
        if (eWeight > EJNWEIGHTMAX) break;
      }
    }

    // Store original deviation from 120 degree topology.
    if (iter == 1) errInCM = pow2(costheta(pWTinJRF[0], pWTinJRF[1]) + 0.5)
      + pow2(costheta(pWTinJRF[0], pWTinJRF[2]) + 0.5)
      + pow2(costheta(pWTinJRF[1], pWTinJRF[2]) + 0.5);

    // Find new JRF from the set of weighted momenta.
    Mstep = junctionRestFrame( pWTinJRF[0], pWTinJRF[1], pWTinJRF[2]);
    // Fortran code will not take full step after the first few
    // iterations. How implement this in terms of an M matrix??
    MtoJRF.rotbst( Mstep );
  } while (iter < 3 || (Mstep.deviation() > CONVJNREST && iter < NTRYJNREST) );

  // If final deviation from 120 degrees is bigger than in CM then revert.
  double errInJRF = pow2(costheta(pWTinJRF[0], pWTinJRF[1]) + 0.5)
    + pow2(costheta(pWTinJRF[0], pWTinJRF[2]) + 0.5)
    + pow2(costheta(pWTinJRF[1], pWTinJRF[2]) + 0.5);
  if (errInJRF > errInCM + CONVJNREST) {
    infoPtr->errorMsg("Warning in StringFragmentation::fragmentTo"
      "Junction: bad convergence junction rest frame");
    MtoJRF.reset();
    MtoJRF.bstback(pSum);
  }

  // Opposite operation: boost from JRF to original system.
  RotBstMatrix MfromJRF = MtoJRF;
  MfromJRF.invert();

  // Sum leg four-momenta in original frame and in JRF.
  Vec4 pInLeg[3], pInJRF[3];
  for (leg = 0; leg < 3; ++leg) {
    pInLeg[leg] = 0.;
    for (int i = legBeg[leg]; i <= legEnd[leg]; ++i)
      pInLeg[leg] += event[ iParton[i] ].p();
    pInJRF[leg] = pInLeg[leg];
    pInJRF[leg].rotbst(MtoJRF);
  }

  // Pick the two legs with lowest energy in JRF.
  int legMin = (pInJRF[0].e() < pInJRF[1].e()) ? 0 : 1;
  int legMax = 1 - legMin;
  if (pInJRF[2].e() < min(pInJRF[0].e(), pInJRF[1].e()) ) legMin = 2;
  else if (pInJRF[2].e() > max(pInJRF[0].e(), pInJRF[1].e()) ) legMax = 2;
  int legMid = 3 - legMin - legMax;

  // Save info on which status codes belong with the three legs.
  int iJunction = (-iParton[0]) / 10 - 1;
  event.statusJunction( iJunction, legMin, 85);
  event.statusJunction( iJunction, legMid, 86);
  event.statusJunction( iJunction, legMax, 83);

  // Temporarily copy the partons on the low-energy legs, into the JRF,
  // in reverse order, so (anti)quark leg end first.
  vector<int> iPartonMin;
  for (int i = legEnd[legMin]; i >= legBeg[legMin]; --i) {
    int iNew = event.append( event[ iParton[i] ] );
    event[iNew].rotbst(MtoJRF);
    iPartonMin.push_back( iNew );
  }
  vector<int> iPartonMid;
  for (int i = legEnd[legMid]; i >= legBeg[legMid]; --i) {
    int iNew = event.append( event[ iParton[i] ] );
    event[iNew].rotbst(MtoJRF);
    iPartonMid.push_back( iNew );
  }

  // Find final weighted sum of momenta on each of the two legs.
  double eWeight = 0.;
  pWTinJRF[legMin] = 0.;
  for (int i = iPartonMin.size() - 1; i >= 0; --i) {
    pWTinJRF[legMin] += event[ iPartonMin[i] ].p() * exp(-eWeight);
    eWeight += event[ iPartonMin[i] ].e() / eNormJunction;
    if (eWeight > EJNWEIGHTMAX) break;
  }
  eWeight = 0.;
  pWTinJRF[legMid] = 0.;
  for (int i = iPartonMid.size() - 1; i >= 0; --i) {
    pWTinJRF[legMid] += event[ iPartonMid[i] ].p() * exp(-eWeight);
    eWeight += event[ iPartonMid[i] ].e() / eNormJunction;
    if (eWeight > EJNWEIGHTMAX) break;
  }

  // Define fictitious opposing partons in JRF and store as string ends.
  Vec4 pOppose = pWTinJRF[legMin];
  pOppose.flip3();
  int idOppose = (rndmPtr->flat() > 0.5) ? 2 : 1;
  if (event[ iPartonMin[0] ].col() > 0) idOppose = -idOppose;
  int iOppose = event.append( idOppose, 77, 0, 0, 0, 0, 0, 0,
    pOppose, 0.);
  iPartonMin.push_back( iOppose);
  pOppose = pWTinJRF[legMid];
  pOppose.flip3();
  idOppose = (rndmPtr->flat() > 0.5) ? 2 : 1;
  if (event[ iPartonMid[0] ].col() > 0) idOppose = -idOppose;
  iOppose = event.append( idOppose, 77, 0, 0, 0, 0, 0, 0,
    pOppose, 0.);
  iPartonMid.push_back( iOppose);

  // Set up kinematics of string evolution in low-energy temporary systems.
  systemMin.setUp(iPartonMin, event);
  systemMid.setUp(iPartonMid, event);

  // Outer fallback loop, when too little energy left for third leg.
  int idMin = 0;
  int idMid = 0;
  Vec4 pDiquark;
  for ( int iTryOuter = 0; ; ++iTryOuter) {

    // Middle fallback loop, when much unused energy in leg remnants.
    double eLeftMin = 0.;
    double eLeftMid = 0.;
    for ( int iTryMiddle = 0; ; ++iTryMiddle) {

      // Loop over the two lowest-energy legs.
      for (int legLoop = 0; legLoop < 2; ++ legLoop) {
        int legNow = (legLoop == 0) ? legMin : legMid;

        // Read in properties specific to this leg.
        StringSystem& systemNow = (legLoop == 0) ? systemMin : systemMid;
        int idPos = (legLoop == 0) ? event[ iPartonMin[0] ].id()
          : event[ iPartonMid[0] ].id();
        idOppose = (legLoop == 0) ? event[ iPartonMin.back() ].id()
          : event[ iPartonMid.back() ].id();
        double eInJRF = pInJRF[legNow].e();
        int statusHad = (legLoop == 0) ? 85 : 86;

        // Inner fallback loop, when a diquark comes in to junction.
        double eUsed = 0.;
        for ( int iTryInner = 0; ; ++iTryInner) {
          if (iTryInner > 2 * NTRYJNMATCH) {
            infoPtr->errorMsg("Error in StringFragmentation::fragment"
              "ToJunction: caught in junction flavour loop");
            event.popBack( iPartonMin.size() + iPartonMid.size() );
            return false;
          }
          bool needBaryon = (abs(idPos) > 10 && iTryInner > NTRYJNMATCH);
          double eExtra   = (iTryInner > NTRYJNMATCH) ? EEXTRAJNMATCH : 0.;

          // Set up two string ends, and begin fragmentation loop.
          setStartEnds(idPos, idOppose, systemNow);
          eUsed = 0.;
          int nHadrons = 0;
          bool noNegE = true;
          for ( ; ; ++nHadrons) {

            // Construct trial hadron from positive end.
            posEnd.newHadron();
            Vec4 pHad = posEnd.kinematicsHadron(systemNow);

            // Negative energy signals failure in construction.
            if (pHad.e() < 0. ) { noNegE = false; break; }

            // Break if passed system midpoint ( = junction) in energy.
            // Exceptions: small systems, and/or with diquark end.
            bool delayedBreak = false;
            if (eUsed + pHad.e() + eExtra > eInJRF) {
              if (nHadrons > 0 || !needBaryon) break;
              delayedBreak = true;
            }

            // Else construct kinematics of the new hadron and store it.
            hadrons.append( posEnd.idHad, statusHad, iPos, iNeg,
              0, 0, 0, 0, pHad, posEnd.mHad);

            // Update string end and remaining momentum.
            posEnd.update();
            eUsed += pHad.e();

            // Delayed break in small systems, and/or with diquark end.
            if (delayedBreak) {
              ++nHadrons;
              break;
            }
          }

          // Possible to produce zero hadrons if the endpoint is not a diquark.
          if (iTryInner > NTRYJNMATCH && !noNegE && nHadrons == 0 &&
            abs(idPos) < 10) break;

          // End of fragmentation loop. Inner loopback if ends on a diquark.
          if ( noNegE && abs(posEnd.flavOld.id) < 10 ) break;
          hadrons.popBack(nHadrons);
        }

        // End of one-leg fragmentation. Store end quark and remnant energy.
        if (legNow == legMin) {
          idMin = posEnd.flavOld.id;
          eLeftMin = eInJRF - eUsed;
        } else {
          idMid = posEnd.flavOld.id;
          eLeftMid = eInJRF - eUsed;
        }
      }

      // End of both-leg fragmentation.
      // Middle loopback if too much energy left.
      double eTrial = eBothLeftJunction + rndmPtr->flat() * eMaxLeftJunction;
      if (iTryMiddle > NTRYJNMATCH
        || ( min( eLeftMin, eLeftMid) < eBothLeftJunction
        && max( eLeftMin, eLeftMid) < eTrial ) ) break;
      hadrons.clear();
    }

    // Boost hadrons away from the JRF to the original frame.
    for (int i = 0; i < hadrons.size(); ++i) {
      hadrons[i].rotbst(MfromJRF);
      // Recalculate energy to compensate for numerical precision loss
      // in iterative calculation of MfromJRF.
      hadrons[i].e( hadrons[i].eCalc() );
      pJunctionHadrons += hadrons[i].p();
    }

    // Outer loopback if combined diquark mass too negative
    // or too little energy left in third leg.
    pDiquark = pInLeg[legMin] + pInLeg[legMid] - pJunctionHadrons;
    double m2Left = m2( pInLeg[legMax], pDiquark);
    if (iTryOuter >  NTRYJNMATCH || (pDiquark.mCalc() > MDIQUARKMIN
      && m2Left > eMinLeftJunction * pInLeg[legMax].e()) ) break;
    hadrons.clear();
    pJunctionHadrons = 0.;
  }

  // Now found solution; no more loopback. Remove temporary parton copies.
  event.popBack( iPartonMin.size() + iPartonMid.size() );

  // Construct and store an effective diquark string end from the
  // two remnant quark ends, for temporary usage.
  int    idDiquark = flavSelPtr->makeDiquark( idMin, idMid);
  double mDiquark  = pDiquark.mCalc();
  int    iDiquark  = event.append( idDiquark, 78, 0, 0, 0, 0, 0, 0,
    pDiquark, mDiquark);

  // Find the partons on the last leg, again in reverse order.
  vector<int> iPartonMax;
  for (int i = legEnd[legMax]; i >= legBeg[legMax]; --i)
    iPartonMax.push_back( iParton[i] );
  iPartonMax.push_back( iDiquark );

  // Recluster gluons nearby to diquark end when taken too much energy.
  int iPsize       = iPartonMax.size();
  double m0Diquark = event[iDiquark].m0();
  while (iPsize > 2) {
    Vec4 pGluNear = event[ iPartonMax[iPsize - 2] ].p();
    if ( (pDiquark + 0.5 * pGluNear).mCalc() > m0Diquark + mJoin ) break;
    pDiquark += pGluNear;
    event[iDiquark].p( pDiquark );
    event[iDiquark].m( pDiquark.mCalc() );
    iPartonMax.pop_back();
    --iPsize;
    iPartonMax[iPsize - 1] = iDiquark;
  }

  // Modify parton list to remaining leg + remnant of the first two.
  iParton = iPartonMax;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Find the boost matrix to the rest frame of a junction,
// given the three respective endpoint four-momenta.

RotBstMatrix StringFragmentation::junctionRestFrame(Vec4& p0, Vec4& p1,
  Vec4& p2) {

  // Calculate masses and other invariants.
  Vec4 pSumJun  = p0 + p1 + p2;
  double sHat   = pSumJun.m2Calc();
  double pp[3][3];
  pp[0][0]      = p0.m2Calc();
  pp[1][1]      = p1.m2Calc();
  pp[2][2]      = p2.m2Calc();
  pp[0][1] = pp[1][0] = p0 * p1;
  pp[0][2] = pp[2][0] = p0 * p2;
  pp[1][2] = pp[2][1] = p1 * p2;

  // Requirement (eiMax)_j = pi*pj/mj < (eiMax)_k = pi*pk/mk, used below,
  // here rewritten as pi*pj * mk < pi*pk * mj and squared.
  double eMax01 = pow2(pp[0][1]) * pp[2][2];
  double eMax02 = pow2(pp[0][2]) * pp[1][1];
  double eMax12 = pow2(pp[1][2]) * pp[0][0];

  // Initially pick i to be the most massive parton. but allow other tries.
  int i = (pp[1][1] > pp[0][0]) ? 1 : 0;
  if (pp[2][2] > max(pp[0][0], pp[1][1])) i = 2;
  int j, k;
  double ei     = 0.;
  double ej     = 0.;
  double ek     = 0.;
  for (int iTry = 0; iTry < 3; ++iTry) {

    // Pick j to give minimal eiMax, and k the third vector.
    if (i == 0) j = (eMax02 < eMax01) ? 2 : 1;
    else if (i == 1) j = (eMax12 < eMax01) ? 2 : 0;
    else j = (eMax12 < eMax02) ? 1 : 0;
    k = 3 - i - j;

    // Alternative names according to i, j, k conventions.
    double m2i  = pp[i][i];
    double m2j  = pp[j][j];
    double m2k  = pp[k][k];
    double pipj = pp[i][j];
    double pipk = pp[i][k];
    double pjpk = pp[j][k];

    // Trivial to find new parton energies if all three partons are massless.
    if (m2i < M2MAXJRF) {
      ei        = sqrt( 2. * pipk * pipj / (3. * pjpk) );
      ej        = sqrt( 2. * pjpk * pipj / (3. * pipk) );
      ek        = sqrt( 2. * pipk * pjpk / (3. * pipj) );

    // Else find three-momentum range for parton i and values at extremes.
    } else {
      // Minimum when i is at rest.
      double piMin = 0.;
      double eiMin = sqrt(m2i);
      double ejMin = pipj / eiMin;
      double ekMin = pipk / eiMin;
      double pjMin = sqrtpos( ejMin*ejMin - m2j );
      double pkMin = sqrtpos( ekMin*ekMin - m2k );
      double fMin  = ejMin * ekMin + 0.5 * pjMin * pkMin - pjpk;
      // Maximum estimated when j + k is at rest, alternatively j at rest.
      double eiMax = (pipj + pipk) / sqrt(m2j + m2k + 2. * pjpk);
      if (m2j > M2MAXJRF) eiMax = min( eiMax, pipj / sqrt(m2j) );
      double piMax = sqrtpos( eiMax*eiMax - m2i );
      double temp  = eiMax*eiMax - 0.25 *piMax*piMax;
      double pjMax = (eiMax * sqrtpos( pipj*pipj - m2j * temp )
        - 0.5 * piMax * pipj) / temp;
      double pkMax = (eiMax * sqrtpos( pipk*pipk - m2k * temp )
        - 0.5 * piMax * pipk) / temp;
      double ejMax = sqrt(pjMax*pjMax + m2j);
      double ekMax = sqrt(pkMax*pkMax + m2k);
      double fMax  = ejMax * ekMax + 0.5 * pjMax * pkMax - pjpk;

      // If unexpected values at upper endpoint then pick another parton.
      if (fMax > 0.) {
        int iPrel = (i + 1)%3;
        if (pp[iPrel][iPrel] > M2MAXJRF) {i = iPrel; continue;}
        ++iTry;
        iPrel = (i + 2)%3;
        if (iTry < 3 && pp[iPrel][iPrel] > M2MAXJRF) {i = iPrel; continue;}
      }

      // Start binary + linear search to find solution inside range.
      int iterMin = 0;
      int iterMax = 0;
      double pi   = 0.5 * (piMin + piMax);
      for (int iter = 0; iter < NTRYJRFEQ; ++iter) {

        // Derive momentum of other two partons and distance to root.
        ei = sqrt(pi*pi + m2i);
        temp = ei*ei - 0.25 * pi*pi;
        double pj = (ei * sqrtpos( pipj*pipj - m2j * temp )
          - 0.5 * pi * pipj) / temp;
        double pk = (ei * sqrtpos( pipk*pipk - m2k * temp )
          - 0.5 * pi * pipk) / temp;
        ej = sqrt(pj*pj + m2j);
        ek = sqrt(pk*pk + m2k);
        double fNow = ej * ek + 0.5 * pj * pk - pjpk;

        // Replace lower or upper bound by new value.
        if (fNow > 0.) { ++iterMin; piMin = pi; fMin = fNow;}
        else {++iterMax; piMax = pi; fMax = fNow;}

        // Pick next i momentum to explore, hopefully closer to root.
        if (2 * iter < NTRYJRFEQ
          && (iterMin < 2 || iterMax < 2 || 4 * iter < NTRYJRFEQ))
          { pi = 0.5 * (piMin + piMax); continue;}
        if (fMin < 0. || fMax > 0. || abs(fNow) < CONVJRFEQ * sHat) break;
        pi = piMin + (piMax - piMin) * fMin / (fMin - fMax);
      }

    // If arrived here then either succeeded or exhausted possibilities.
    } break;
  }

  // Now we know the energies in the junction rest frame.
  double eNew[3] = { 0., 0., 0.};
  eNew[i] = ei;
  eNew[j] = ej;
  eNew[k] = ek;

  // Boost (copy of) partons to their rest frame.
  RotBstMatrix Mmove;
  Vec4 p0cm = p0;
  Vec4 p1cm = p1;
  Vec4 p2cm = p2;
  Mmove.bstback(pSumJun);
  p0cm.rotbst(Mmove);
  p1cm.rotbst(Mmove);
  p2cm.rotbst(Mmove);

  // Construct difference vectors and the boost to junction rest frame.
  Vec4 pDir01      = p0cm / p0cm.e() - p1cm / p1cm.e();
  Vec4 pDir02      = p0cm / p0cm.e() - p2cm / p2cm.e();
  double pDiff01   = pDir01.pAbs2();
  double pDiff02   = pDir02.pAbs2();
  double pDiff0102 = dot3(pDir01, pDir02);
  double eDiff01   = eNew[0] / p0cm.e() - eNew[1] / p1cm.e();
  double eDiff02   = eNew[0] / p0cm.e() - eNew[2] / p2cm.e();
  double denom     = pDiff01 * pDiff02 - pDiff0102*pDiff0102;
  double coef01    = (eDiff01 * pDiff02 - eDiff02 * pDiff0102) / denom;
  double coef02    = (eDiff02 * pDiff01 - eDiff01 * pDiff0102) / denom;
  Vec4 vJunction   = coef01 * pDir01 + coef02 * pDir02;
  vJunction.e( sqrt(1. + vJunction.pAbs2()) );

  // Add two boosts, giving final result.
  Mmove.bst(vJunction);
  return Mmove;

}

//--------------------------------------------------------------------------

// When string fragmentation has failed several times,
// try to join some more nearby partons.

int StringFragmentation::extraJoin(double facExtra, Event& event) {

  // Keep on looping while pairs found below joining threshold.
  int nJoin  = 0;
  int iPsize = iParton.size();
  while (iPsize > 2) {

    // Look for the pair of neighbour partons (along string) with
    // the smallest invariant mass (subtracting quark masses).
    int iJoinMin    = -1;
    double mJoinMin = 2. * facExtra * mJoin;
    for (int i = 0; i < iPsize - 1; ++i) {
      Particle& parton1 = event[ iParton[i] ];
      Particle& parton2 = event[ iParton[i + 1] ];
      Vec4 pSumNow;
      pSumNow += (parton2.isGluon()) ? 0.5 * parton1.p() : parton1.p();
      pSumNow += (parton2.isGluon()) ? 0.5 * parton2.p() : parton2.p();
      double mJoinNow = pSumNow.mCalc();
      if (!parton1.isGluon()) mJoinNow -= parton1.m0();
      if (!parton2.isGluon()) mJoinNow -= parton2.m0();
      if (mJoinNow < mJoinMin) { iJoinMin = i; mJoinMin = mJoinNow; }
    }

    // Decide whether to join, if not finished.
    if (iJoinMin < 0 || mJoinMin > facExtra * mJoin) return nJoin;
    ++nJoin;

    // Create new joined parton.
    int iJoin1  = iParton[iJoinMin];
    int iJoin2  = iParton[iJoinMin + 1];
    int idNew   = (event[iJoin1].isGluon()) ? event[iJoin2].id()
                                            : event[iJoin1].id();
    int colNew  = event[iJoin1].col();
    int acolNew = event[iJoin2].acol();
    if (colNew == acolNew) {
      colNew    = event[iJoin2].col();
      acolNew   = event[iJoin1].acol();
    }
    Vec4 pNew   = event[iJoin1].p() + event[iJoin2].p();

    // Append joined parton to event record and reduce parton list.
    int iNew = event.append( idNew, 73, min(iJoin1, iJoin2),
      max(iJoin1, iJoin2), 0, 0, colNew, acolNew, pNew, pNew.mCalc() );
    iParton[iJoinMin] = iNew;
    for (int i = iJoinMin + 1; i < iPsize - 1; ++i)
      iParton[i] = iParton[i + 1];
    iParton.pop_back();
    --iPsize;

  // Done.
  }
  return nJoin;
}

//==========================================================================

} // end namespace Pythia8
