// StringFragmentation.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the classes for string fragmentation.
// StringEnd: keeps track of the fragmentation step.
// StringFragmentation: is the top-level class.

#ifndef Pythia8_StringFragmentation_H
#define Pythia8_StringFragmentation_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/FragmentationSystems.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Ropewalk.h"
#include "Pythia8/Settings.h"
#include "Pythia8/UserHooks.h"

namespace Pythia8 {

//==========================================================================

// The StringEnd class contains the information related to
// one of the current endpoints of the string system.
// Only to be used inside StringFragmentation, so no private members.

class StringEnd {

public:

  // Constructor.
  StringEnd() : particleDataPtr(), flavSelPtr(), pTSelPtr(), zSelPtr(),
    fromPos(), thermalModel(), mT2suppression(), iEnd(), iMax(), idHad(),
    iPosOld(), iNegOld(), iPosNew(), iNegNew(), hadSoFar(), colOld(), colNew(),
    pxOld(), pyOld(), pxNew(), pyNew(), pxHad(), pyHad(), mHad(), mT2Had(),
    zHad(), GammaOld(), GammaNew(), xPosOld(), xPosNew(), xPosHad(), xNegOld(),
    xNegNew(), xNegHad(), aLund(), bLund() {}

  // Save pointers.
  void init( ParticleData* particleDataPtrIn, StringFlav* flavSelPtrIn,
    StringPT* pTSelPtrIn, StringZ* zSelPtrIn, Settings& settings) {
    particleDataPtr = particleDataPtrIn; flavSelPtr = flavSelPtrIn;
    pTSelPtr = pTSelPtrIn; zSelPtr = zSelPtrIn;
    bLund = zSelPtr->bAreaLund(); aLund = zSelPtr->aAreaLund();
    thermalModel   = settings.flag("StringPT:thermalModel");
    mT2suppression = settings.flag("StringPT:mT2suppression"); }

  // Set up initial endpoint values from input.
  void setUp(bool fromPosIn, int iEndIn, int idOldIn, int iMaxIn,
    double pxIn, double pyIn, double GammaIn, double xPosIn,
    double xNegIn, int colIn);

  // Fragment off one hadron from the string system, in flavour and pT.
  void newHadron(double nNSP = 0.0);

  // Fragment off one hadron from the string system, in momentum space,
  // by taking steps either from positive or from negative end.
  Vec4 kinematicsHadron(StringSystem& system,
    vector<StringVertex>& stringVertices, bool useInputZ = false,
    double zHadIn = 0.);

  // Generate momentum for some possible next hadron, based on mean values
  // to get an estimate for rapidity and pT.
  Vec4 kinematicsHadronTmp(StringSystem system, Vec4 pRem, double phi,
    double mult);

  // Update string end information after a hadron has been removed.
  void update();

  // Constants: could only be changed in the code itself.
  static const double TINY, PT2SAME, MEANMMIN, MEANM, MEANPT;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointers to classes for flavour, pT and z generation.
  StringFlav*   flavSelPtr;
  StringPT*     pTSelPtr;
  StringZ*      zSelPtr;

  // Data members.
  bool   fromPos, thermalModel, mT2suppression;
  int    iEnd, iMax, idHad, iPosOld, iNegOld, iPosNew, iNegNew, hadSoFar,
         colOld, colNew;
  double pxOld, pyOld, pxNew, pyNew, pxHad, pyHad, mHad, mT2Had, zHad,
         GammaOld, GammaNew, xPosOld, xPosNew, xPosHad, xNegOld, xNegNew,
         xNegHad, aLund, bLund;
  FlavContainer flavOld, flavNew;
  Vec4   pHad, pSoFar;

};

//==========================================================================

// The StringFragmentation class contains the top-level routines
// to fragment a colour singlet partonic system.

class StringFragmentation {

public:

  // Constructor.
  StringFragmentation() : infoPtr(), particleDataPtr(), rndmPtr(),
    flavSelPtr(), pTSelPtr(), zSelPtr(), flavRopePtr(), userHooksPtr(),
    closePacking(), doFlavRope(), setVertices(), constantTau(), smearOn(),
    traceColours(false), hadronVertex(), stopMass(), stopNewFlav(),
    stopSmear(), eNormJunction(), eBothLeftJunction(), eMaxLeftJunction(),
    eMinLeftJunction(), mJoin(), bLund(), pT20(), xySmear(), kappaVtx(),
    mc(), mb(), hasJunction(), isClosed(), iPos(), iNeg(), w2Rem(),
    stopMassNow(), idDiquark(), legMin(), legMid() {}

  // Initialize and save pointers.
  void init(Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    StringFlav* flavSelPtrIn, StringPT* pTSelPtrIn, StringZ* zSelPtrIn,
    FlavourRope* flavRopePtrIn = NULL, UserHooks* userHooksPtrIn = NULL);

  // Do the fragmentation: driver routine.
  bool fragment( int iSub, ColConfig& colConfig, Event& event);

  // Find the boost matrix to the rest frame of a junction.
  RotBstMatrix junctionRestFrame(Vec4& p0, Vec4& p1, Vec4& p2);

private:

  // Constants: could only be changed in the code itself.
  static const int    NTRYFLAV, NTRYJOIN, NSTOPMASS, NTRYJNREST,
                      NTRYJNMATCH, NTRYJRFEQ;
  static const double FACSTOPMASS, CLOSEDM2MAX, CLOSEDM2FRAC, EXPMAX,
                      MATCHPOSNEG, EJNWEIGHTMAX, CONVJNREST, M2MAXJRF,
                      M2MINJRF, EEXTRAJNMATCH, MDIQUARKMIN, CONVJRFEQ,
                      CHECKPOS;

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Pointers to classes for flavour, pT and z generation.
  StringFlav*   flavSelPtr;
  StringPT*     pTSelPtr;
  StringZ*      zSelPtr;

  // Pointer to flavour-composition-changing ropes.
  FlavourRope*  flavRopePtr;

  // Pointer to the User Hooks class for user intervention
  UserHooks*    userHooksPtr;

  // Initialization data, read from Settings.
  bool   closePacking, doFlavRope, setVertices, constantTau, smearOn,
         traceColours;
  int    hadronVertex;
  double stopMass, stopNewFlav, stopSmear, eNormJunction,
         eBothLeftJunction, eMaxLeftJunction, eMinLeftJunction,
         mJoin, bLund, pT20, xySmear, kappaVtx, mc, mb;

  // Data members.
  bool   hasJunction, isClosed;
  int    iPos, iNeg;
  double w2Rem, stopMassNow;
  Vec4   pSum, pRem, pJunctionHadrons;

  // List of partons in string system.
  vector<int> iParton, iPartonMinLeg, iPartonMidLeg, iPartonMax;

  // Vertex information from the fragmentation process.
  vector<StringVertex> stringVertices, legMinVertices, legMidVertices;

  // Boost from/to rest frame of a junction to original frame.
  RotBstMatrix MfromJRF, MtoJRF;

  // Information on diquark created at the junction.
  int    idDiquark;

  // Fictitious opposing partons in JRF: string ends for vertex location.
  Vec4 pMinEnd, pMidEnd;

  // Temporary event record for the produced particles.
  Event hadrons;

  // Information on the system of string regions.
  StringSystem system, systemMin, systemMid;

  // Information on the two current endpoints of the fragmenting system.
  StringEnd posEnd, negEnd;

  // Find region where to put first string break for closed gluon loop.
  vector<int> findFirstRegion(int iSub, ColConfig& colConfig, Event& event);

  // Set flavours and momentum position for initial string endpoints.
  void setStartEnds(int idPos, int idNeg, StringSystem systemNow,
    int legNow = 3);

  // Check remaining energy-momentum whether it is OK to continue.
  bool energyUsedUp(bool fromPos);

  // Produce the final two partons to complete the system.
  bool finalTwo(bool fromPos, Event& event, bool usedPosJun, bool usedNegJun,
  double nNSP);

  // Final region information.
  Vec4 pPosFinalReg, pNegFinalReg, eXFinalReg, eYFinalReg;

  // Set hadron production points in space-time picture.
  void setHadronVertices(Event& event);

  // Construct a special joining region for the final two hadrons.
  StringRegion finalRegion();

  // Store the hadrons in the normal event record, ordered from one end.
  void store(Event& event);

  // Fragment off two of the string legs in to a junction.
  bool fragmentToJunction(Event& event);

  // Initially considered legs from the junction.
  int legMin, legMid;

  // Join extra nearby partons when stuck.
  int extraJoin(double facExtra, Event& event);

  // Get the number of nearby strings given the energies.
  double nearStringPieces(StringEnd end,
    vector< vector< pair<double,double> > >& rapPairs);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_StringFragmentation_H
