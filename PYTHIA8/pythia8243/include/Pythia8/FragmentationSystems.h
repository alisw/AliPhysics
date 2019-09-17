// FragmentationSystems.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains auxiliary classes in the fragmentation process.
// ColSinglet contains info on an individual singlet.
// ColConfig describes the colour configuration of the whole event.
// StringRegion keeps track on string momenta and directions.
// StringSystem contains all the StringRegions of the colour singlet.
// StringVertex contains information on space-time location of string breaks.

#ifndef Pythia8_FragmentationSystems_H
#define Pythia8_FragmentationSystems_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// The ColSinglet class contains info on an individual singlet.
// Only to be used inside ColConfig, so no private members.

class ColSinglet {

public:

  // Constructors.
  ColSinglet() : pSum(0., 0., 0., 0.), mass(0.), massExcess(0.),
    hasJunction(false), isClosed(false), isCollected(false) {}
  ColSinglet(vector<int>& iPartonIn, Vec4 pSumIn, double massIn,
    double massExcessIn, bool hasJunctionIn = false,
    bool isClosedIn = false, bool isCollectedIn = false)
    : iParton(iPartonIn), pSum(pSumIn), mass(massIn),
    massExcess(massExcessIn), hasJunction(hasJunctionIn),
    isClosed(isClosedIn), isCollected(isCollectedIn) {}

  // Size of iParton array.
  int size() const { return iParton.size();}

  // Stored quantities.
  vector<int> iParton;
  Vec4   pSum;
  double mass, massExcess;
  bool   hasJunction, isClosed, isCollected;

};

//==========================================================================

// The ColConfig class describes the colour configuration of the whole event.

class ColConfig {

public:

  // Constructor.
  ColConfig() : infoPtr(), flavSelPtr(), mJoin(), mJoinJunction(),
    mStringMin() {singlets.resize(0);}

  // Initialize and save pointers.
  void init(Info* infoPtrIn, Settings& settings, StringFlav* flavSelPtrIn);

  // Number of colour singlets.
  int size() const {return singlets.size();}

  // Overload index operator to access separate colour singlets.
  ColSinglet& operator[](int iSub) {return singlets[iSub];}

  // Clear contents.
  void clear() {singlets.resize(0);}

  // Insert a new colour singlet system in ascending mass order.
  // Calculate its properties. Join nearby partons.
  bool insert( vector<int>& iPartonIn, Event& event);

  // Insert a new qqbar colour singlet system in ascending mass order.
  // Calculate its properties.
  bool simpleInsert( vector<int>& iPartonIn, Event& event);

  // Erase a colour singlet system. (Rare operation.)
  void erase(int iSub) {singlets.erase(singlets.begin() + iSub);}

  // Collect all partons of singlet to be consecutively ordered.
  void collect(int iSub, Event& event, bool skipTrivial = true);

  // Find to which singlet system a particle belongs.
  int findSinglet(int i);

  // List all currently identified singlets.
  void list() const;

  // Rapidity range [y_min, y_max] of all string pieces in all singlets.
  // Only used when stringPT:closePacking is on.
  vector< vector< pair<double,double> > > rapPairs;

private:

  // Constants: could only be changed in the code itself.
  static const double CONSTITUENTMASS;

  // Pointer to various information on the generation.
  Info*       infoPtr;

  // Pointer to class for flavour generation.
  StringFlav* flavSelPtr;

  // Initialization data, to be read from Settings.
  double mJoin, mJoinJunction, mStringMin;

  // List of all separate colour singlets.
  vector<ColSinglet> singlets;

  // Join two legs of junction to a diquark for small invariant masses.
  bool joinJunction( vector<int>& iPartonIn, Event& event,
    double massExcessIn);

};

//==========================================================================

// The StringRegion class contains the information related to
// one string section in the evolution of a multiparton system.
// Only to be used inside StringFragmentation and MiniStringFragmentation,
// so no private members.

class StringRegion {

public:

  // Constructor.
  StringRegion() : isSetUp(false), isEmpty(true), w2(0.), xPosProj(0.),
    xNegProj(0.), pxProj(0.), pyProj(0.), colPos(0), colNeg(0) {}

  // Constants: could only be changed in the code itself.
  static const double MJOIN, TINY;

  // Data members.
  bool   isSetUp, isEmpty;
  Vec4   pPos, pNeg, eX, eY, pPosMass, pNegMass, massOffset;
  double w2, xPosProj, xNegProj, pxProj, pyProj;
  int    colPos, colNeg;

  // Calculate offset of the region from parton list. Special junction case.
  Vec4 gluonOffset(vector<int>& iSys, Event& event, int iPos, int iNeg);
  Vec4 gluonOffsetJRF(vector<int>& iSys, Event& event, int iPos, int iNeg,
    RotBstMatrix MtoJRF);

  // If massive case, the offset of the initial regions is calculated.
  bool massiveOffset(int iPos, int iNeg, int iMax, int id1, int id2,
    double mc, double mb);

  // Set up four-vectors for longitudinal and transverse directions.
  void setUp(Vec4 p1, Vec4 p2, int col1, int col2, bool isMassless = false);

  // Construct a four-momentum from (x+, x-, px, py).
  Vec4 pHad( double xPosIn, double xNegIn, double pxIn, double pyIn)
    { return xPosIn * pPos + xNegIn * pNeg + pxIn * eX + pyIn * eY; }

  // Project a four-momentum onto (x+, x-, px, py). Read out projection.
  void project(Vec4 pIn);
  void project( double pxIn, double pyIn, double pzIn, double eIn)
    { project( Vec4( pxIn, pyIn, pzIn, eIn) ); }
  double xPos() const {return xPosProj;}
  double xNeg() const {return xNegProj;}
  double px() const {return pxProj;}
  double py() const {return pyProj;}

};

//==========================================================================

// The StringSystem class contains the complete set of all string regions.
// Only to be used inside StringFragmentation, so no private members.

class StringSystem {

public:

  // Constructor.
  StringSystem() : sizePartons(), sizeStrings(), sizeRegions(), indxReg(),
    iMax(), mJoin(), m2Join() {}

  // Set up system from parton list.
  void setUp(vector<int>& iSys, Event& event);

  // Calculate string region from (iPos, iNeg) pair.
  int iReg( int iPos, int iNeg) const
    {return (iPos * (indxReg - iPos)) / 2 + iNeg;}

  // Reference to string region specified by (iPos, iNeg) pair.
  StringRegion& region(int iPos, int iNeg) {return system[iReg(iPos, iNeg)];}

  // Reference to low string region specified either by iPos or iNeg.
  StringRegion& regionLowPos(int iPos) {
    return system[iReg(iPos, iMax - iPos)]; }
  StringRegion& regionLowNeg(int iNeg) {
    return system[iReg(iMax - iNeg, iNeg)]; }

  // Main content: a vector with all the string regions of the system.
  vector<StringRegion> system;

  // Other data members.
  int    sizePartons, sizeStrings, sizeRegions, indxReg, iMax;
  double mJoin, m2Join;

};

//==========================================================================

// The StringVertex class contains the space-time vertex location information
// stored during the fragmentation process. No private members.

class StringVertex {

public:

  // Constructors.
  StringVertex(bool fromPosIn = true, int iRegPosIn = 0,
    int iRegNegIn = 0, double xRegPosIn = 0., double xRegNegIn = 0.)
    : fromPos(fromPosIn), iRegPos(iRegPosIn), iRegNeg(iRegNegIn),
    xRegPos(xRegPosIn), xRegNeg(xRegNegIn) { }

  StringVertex(const StringVertex& v): fromPos(v.fromPos),
    iRegPos(v.iRegPos), iRegNeg(v.iRegNeg),
    xRegPos(v.xRegPos), xRegNeg(v.xRegNeg) { }

  StringVertex& operator = (const StringVertex& v) {if (this != &v)
    {fromPos = v.fromPos; iRegPos = v.iRegPos; iRegNeg = v.iRegNeg;
    xRegPos = v.xRegPos; xRegNeg = v.xRegNeg;} return *this; }

  // Variable members.
  bool fromPos;
  int iRegPos, iRegNeg;
  double xRegPos, xRegNeg;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_FragmentationSystems_H
