// Ropewalk.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Rope Hadronization. The Ropewalk takes care of setting
// up string geometry and calculating overlaps etc. The RopeDipole classes
// take care of the dynamics of string shoving. Flavour-composition-changing
// ropes are handled by FlavourRope which changes parameters, and RopeFragPars
// which calculates parameters for the rope.
//
// The file contains the following classes: RopeDipoleEnd,
// OverlappingRopeDipole, RopeDipole, Ropewalk, RopeFragPars and FlavourRope.

#ifndef Pythia8_Ropewalk_H
#define Pythia8_Ropewalk_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationSystems.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/Settings.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==================================================================

// Define the end of a dipole, containing a pointer to the particle,
// and its index in the event record.
// Includes some methods for kinematics output.

class RopeDipoleEnd {

public:

  // Constructors sets event pointer and event record index.
  RopeDipoleEnd() : e(NULL), ne(-1) { }
  RopeDipoleEnd(Event* eIn, int neIn) : e(eIn), ne(neIn) { }

  // Get a pointer to the particle.
  Particle* getParticlePtr() { if (!e) return NULL; return &(*e)[ne]; }

  // Get the particle index in event record.
  int getNe() {return ne;}

  // Output methods for (modified) rapidity.
  double labrap() { return getParticlePtr()->y(); }
  double rap(double m0){ return getParticlePtr()->y(m0); }
  double rap(double m0, RotBstMatrix& r) { return getParticlePtr()->y(m0, r); }

private:

  // Pointer to the event and the particle index in event record.
  Event* e;
  int ne;

};

//==================================================================

// A dipole has many dipoles overlapping with it. The OverlappingRopeDipole
// class does bookkeeping of this. Holds a pointer to the original dipole.

// Forward declaration of RopeDipole class.
class RopeDipole;

class OverlappingRopeDipole {

public:

  // Constructor sets up coordinates in the rest frame of other dipole.
  OverlappingRopeDipole(RopeDipole* d, double m0, RotBstMatrix& r);

  // Calculate the overlap at given y and b.
  bool overlap(double y, Vec4 ba, double r0);

  // Has the dipole been hadronized?
  bool hadronized();

private:

  // Pointer to the original dipole.
  RopeDipole* dipole;

// Data members are made public rather than making
// the RopeDipole class a friend.
public:

  int dir;
  double y1, y2;
  Vec4 b1, b2;

};

//==================================================================

// The RopeDipole class holds information about a colour dipole, as well as
// functionality to do shoving and to calculate effective string tension.

class RopeDipole {

public:

  // The RopeDipole constructor makes sure that d1 is always the colored
  // end and d2 the anti-colored.
  RopeDipole(RopeDipoleEnd d1In, RopeDipoleEnd d2In, int iSubIn,
    Info* infoPtrIn);

  // Insert an excitation on dipole, if not already there.
  void addExcitation(double ylab, Particle* ex);

  // Deep access to the RopeDipoleEnds (needed by OverlappingRopeDipole).
  RopeDipoleEnd* d1Ptr() { return &d1; }
  RopeDipoleEnd* d2Ptr() { return &d2; }

  // Get the rotation matrix to go to dipole rest frame.
  RotBstMatrix getDipoleRestFrame();
  RotBstMatrix getDipoleLabFrame();

  // Get the dipole momentum four-vector.
  Vec4 dipoleMomentum();
  // Get the spatial point interpolated to given rapidity.
  // In the dipole rest frame.
  Vec4 bInterpolateDip(double y, double m0);
  // In the lab frame.
  Vec4 bInterpolateLab(double y, double m0);
  // Given a Lorentz matrix.
  Vec4 bInterpolate(double y, RotBstMatrix rb, double m0);

  // Get the quantum numbers m,n characterizing all dipole overlaps
  // at a given rapidity value.
  pair<int, int> getOverlaps(double yfrac, double m0, double r0);

  // Add an overlapping dipole.
  void addOverlappingDipole(OverlappingRopeDipole& d) {
    overlaps.push_back(d); }

  // Get the maximal and minimal rapidity of the dipole.
  double maxRapidity(double m0) { return (max(d1.rap(m0), d2.rap(m0))); }
  double minRapidity(double m0) { return (min(d1.rap(m0), d2.rap(m0))); }

  // Get the maximal and minimal boosted rapidity of the dipole.
  double maxRapidity(double m0, RotBstMatrix& r) { return (max(d1.rap(m0,r),
    d2.rap(m0,r))); }
  double minRapidity(double m0, RotBstMatrix& r) { return (min(d1.rap(m0,r),
    d2.rap(m0,r))); }

  // Propagate the dipole itself.
  void propagateInit(double deltat);

  // Propagate both dipole ends as well as all excitations.
  void propagate(double deltat, double m0);

  // Redistribute momentum to two particles.
  void splitMomentum(Vec4 mom, Particle* p1, Particle* p2, double frac = 0.5);

  // Put gluon excitations on the dipole.
  void excitationsToString(double m0, Event& event);

  // Test if the dipole is hadronized.
  bool hadronized() { return isHadronized; }

  // Get the (event colconfig) index.
  int index() { return iSub; }

  // Recoil the dipole from adding a gluon. If the "dummy" option is set,
  // the recoil will not be added, but only checked.
  // Note: the gluon will not actually be added, only the recoil (if possible).
  bool recoil(Vec4& pg, bool dummy = false);

  // Set dipole hadronized flag.
  void hadronized(bool h) { isHadronized = h; }

  // The number of excitations on the dipole.
  int nExcitations() { return int(excitations.size()); }

private:

  // The ends (ie. particles) of the dipole.
  RopeDipoleEnd d1, d2;

  // The propagated positions in the lab frame.
  Vec4 b1, b2;

  // The string index (internal to the event).
  int iSub;

  // Lorentz matrices to go to and from dipole rest frame.
  RotBstMatrix rotFrom, rotTo;
  bool hasRotFrom, hasRotTo;

  // The dipoles overlapping with this one.
  vector<OverlappingRopeDipole> overlaps;

  // All excitations belonging to this dipole ordered in rapidity in lab frame.
  map<double, Particle*> excitations;

  bool isHadronized;
  Info* infoPtr;

};

//==================================================================

// The Ropewalk class keeps track of all the strings making up ropes
// for shoving as well as flavour enhancement.

class Ropewalk {

public:

  // Constructor.
  Ropewalk() : r0(), m0(), pTcut(), doShoving(), shoveJunctionStrings(),
    shoveMiniStrings(), shoveGluonLoops(), mStringMin(), limitMom(), rCutOff(),
    gAmplitude(), gExponent(), deltay(), deltat(), tShove(), tInit(),
    showerCut(), alwaysHighest(), infoPtr(), rndmPtr() {}

  // The Ropewalk init function sets parameters and pointers.
  bool init(Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn);

  // Extract all dipoles from an event.
  bool extractDipoles(Event& event, ColConfig& colConfig);

  // Calculate all overlaps of all dipoles and store as OverlappingRopeDipoles.
  bool calculateOverlaps();

  // Calculate the effective string tension a fraction yfrac in on the dipole
  // given by indices e1 and e2.
  double getKappaHere(int e1, int e2, double yfrac);

  // The multiplicity of a colour state given its quantum numbers.
  double multiplicity(double p, double q) {
    return ( p < 0 || q < 0 || p + q == 0 )
      ? 0.0 : 0.5 * (p + 1) * (q + 1) * (p + q + 2);
  }

  // Calculate the average string tension of the event, in units of the default
  // string tension (ie. 1 GeV/fm), using random walk in colour space.
  double averageKappa();

  // Invoke the random walk and select a state.
  pair<int, int> select(int m, int n, Rndm* rndm);

  // Shove all dipoles in the event.
  void shoveTheDipoles(Event& event);

private:

  // Parameters of the ropewalk.
  double r0, m0, pTcut;
  // Do shoving flag.
  bool doShoving;
  // Include junction strings in shoving.
  bool shoveJunctionStrings;
  // Include ministrings in shoving.
  bool shoveMiniStrings;
  // Include gluon loops in shoving.
  bool shoveGluonLoops;
  // Limit between string frag and ministring.
  double mStringMin;
  // Limit participating dipoles by their momentum.
  bool limitMom;
  // Radius cutoff in multiples of r0.
  double rCutOff;
  // Parameters of the shoving.
  double gAmplitude, gExponent;
  // Rapidity slicing.
  double deltay;
  // Time steps.
  double deltat;
  // Shove time.
  double tShove;
  // Intial propagation time.
  double tInit;
  // Final state shower cut-off.
  double showerCut;
  // Assume we are always in highest multiplet.
  bool alwaysHighest;

  // Info pointer.
  Info* infoPtr;

  // Pointer the the random number generator.
  Rndm* rndmPtr;

  // All dipoles in the event sorted by event record.
  // Index of the two partons.
  typedef multimap<pair<int,int>, RopeDipole> DMap;
  DMap dipoles;

  // All excitations.
  vector< vector<Particle> > eParticles;

  // Random walk states and their weights.
  vector<pair<int, int> > states;
  vector<double> weights;

  // Private assignment operator.
  Ropewalk& operator=(const Ropewalk);

};

//==================================================================

// RopeFragPars recalculates fragmentation parameters according to a
// changed string tension. Helper class to FlavourRope.

class RopeFragPars {

public:

  // Constructor.
  RopeFragPars() : infoPtr(), aIn(), adiqIn(), bIn(), rhoIn(), xIn(),
    yIn(), xiIn(), sigmaIn(), kappaIn(), aEff(), adiqEff(), bEff(),
    rhoEff(), xEff(), yEff(), xiEff(), sigmaEff(), kappaEff(),
    beta() {}

  // The init function sets up initial parameters from settings.
  void init(Info* infoPtrIn, Settings& settings);

  // Return parameters at given string tension, ordered by their
  // name for easy insertion in settings.
  map<string,double> getEffectiveParameters(double h);

private:

  // Constants: can only be changed in the code itself.
  static const double DELTAA, ACONV, ZCUT;

  // Get the Fragmentation function a parameter from cache or calculate it.
  double getEffectiveA(double thisb, double mT2, bool isDiquark);

  // Calculate the effective parameters.
  bool calculateEffectiveParameters(double h);

  // Insert calculated parameters in cache for later (re-)use.
  bool insertEffectiveParameters(double h);

  // Calculate the a parameter.
  double aEffective(double aOrig, double thisb, double mT2);

  // The Lund fragmentation function.
  double fragf(double z, double a, double b, double mT2);

  // Integral of the Lund fragmentation function with given parameter values.
  double integrateFragFun(double a, double b, double mT2);

  // Helper function for integration.
  double trapIntegrate(double a, double b, double mT2, double sOld, int n);

  // The info pointer.
  Info* infoPtr;

  // Parameter caches to re-use calculations. Sets of parameters, ordered in h.
  map<double, map<string, double> > parameters;

  // Values of the a-parameter ordered in b*mT2 grid.
  map<double, double> aMap;

  // The same for aDiquark.
  map<double, double> aDiqMap;

  // Initial values of parameters.
  double aIn, adiqIn, bIn, rhoIn, xIn, yIn, xiIn, sigmaIn, kappaIn;

  // Effective values of parameters.
  double aEff, adiqEff, bEff, rhoEff, xEff, yEff, xiEff, sigmaEff, kappaEff;

  // Junction parameter.
  double beta;

};

//==================================================================

// The FlavourRope class takes care of placing a string breakup in
// the event, and assigning the string breakup effective parameters.
// It is a UserHooks derived class, and one must make sure to add it
// to the UserHooksVector in the main program or somewhere else.

class FlavourRope {

public:

  // Constructor.
  FlavourRope() : settingsPtr(), rndmPtr(), particleDataPtr(), infoPtr(),
    rwPtr(), ePtr(), doBuffon(), rapiditySpan(), stringProtonRatio(),
    fixedKappa(), h() {}

  // Initialize. Set pointers.
  void init(Settings* settingsPtrIn, Rndm* rndmPtrIn, ParticleData*
    particleDataPtrIn, Info* infoPtrIn, Ropewalk* rwPtrIn) {
    settingsPtr = settingsPtrIn, rndmPtr = rndmPtrIn,
    particleDataPtr = particleDataPtrIn, infoPtr = infoPtrIn,
    rwPtr = rwPtrIn;
    // Initialize event pointer such that it can be tested.
    ePtr = NULL;
    h = settingsPtr->parm("Ropewalk:presetKappa");
    fixedKappa = settingsPtr->flag("Ropewalk:setFixedKappa");
    doBuffon = settingsPtr->flag("Ropewalk:doBuffon");
    rapiditySpan = settingsPtr->parm("Ropewalk:rapiditySpan");
    stringProtonRatio = settingsPtr->parm("Ropewalk:stringProtonRatio");
    // Initialize FragPar.
    fp.init(infoPtr, *settingsPtr);
  }

  // Change the fragmentation parameters.
  bool doChangeFragPar(StringFlav* flavPtr, StringZ* zPtr,
   StringPT * pTPtr, double m2Had, vector<int> iParton, int endId);

  // Set enhancement manually.
  void setEnhancement(double hIn) { h = hIn;}

  // Set pointer to the event.
  void setEventPtr(Event& event) { ePtr = &event;}

private:

  // Find breakup placement and fetch effective parameters.
  // For model depending on vertex information.
  map<string, double> fetchParameters(double m2Had, vector<int> iParton,
    int endId);
  // For simple Buffon model.
  map<string, double> fetchParametersBuffon(double m2Had, vector<int> iParton,
    int endId);

  // Pointer to settings.
  Settings* settingsPtr;

  // Random number generator needed for reinitialization.
  Rndm* rndmPtr;

  // Particle data object needed for reinitialization.
  ParticleData* particleDataPtr;

  // Pythia info pointer.
  Info* infoPtr;

  // Pointer to the ropewalk object.
  Ropewalk* rwPtr;

  // Pointer to the event.
  Event* ePtr;

  // The object which handles change in parameters.
  RopeFragPars fp;

  // Use Buffon prescription without vertex information.
  bool doBuffon;

  // Parameters of Buffon prescription
  double rapiditySpan, stringProtonRatio;

  // Keep track of hadronized dipoles in Buffon prescription.
  vector<int> hadronized;

  // Use preset kappa from settings.
  bool fixedKappa;

  // Locally stored string tension.
  double h;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Ropewalk_H
