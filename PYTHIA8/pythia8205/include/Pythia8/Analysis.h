// Analysis.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the Sphericity, Thrust, ClusterJet and CellJet classes.
// Sphericity: sphericity analysis of the event.
// Thrust: thrust analysis of the event.
// ClusterJet: clustering jet finder.
// CellJet: calorimetric cone jet finder.
// SlowJet: recombination algorithm; lightweight version of FastJet.

#ifndef Pythia8_Analysis_H
#define Pythia8_Analysis_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Sphericity class.
// This class performs (optionally modified) sphericity analysis on an event.

class Sphericity {

public:

  // Constructor.
  Sphericity(double powerIn = 2., int selectIn = 2) : power(powerIn),
    select(selectIn), nFew(0) {powerInt = 0;
    if (abs(power - 1.) < 0.01) powerInt = 1;
    if (abs(power - 2.) < 0.01) powerInt = 2;
    powerMod = 0.5 * power - 1.;}

  // Analyze event.
  bool analyze(const Event& event, ostream& os = cout);

  // Return info on results of analysis.
  double sphericity()      const {return 1.5 * (eVal2 + eVal3);}
  double aplanarity()      const {return 1.5 * eVal3;}
  double eigenValue(int i) const {return (i < 2) ? eVal1 :
    ( (i < 3) ? eVal2 : eVal3 ) ;}
  Vec4 eventAxis(int i)    const {return (i < 2) ? eVec1 :
    ( (i < 3) ? eVec2 : eVec3 ) ;}

  // Provide a listing of the info.
  void list(ostream& os = cout) const;

  // Tell how many events could not be analyzed.
  int nError() const {return nFew;}

private:

  // Constants: could only be changed in the code itself.
  static const int    NSTUDYMIN, TIMESTOPRINT;
  static const double P2MIN, EIGENVALUEMIN;

  // Properties of analysis.
  double power;
  int    select, powerInt;
  double powerMod;

  // Outcome of analysis.
  double eVal1, eVal2, eVal3;
  Vec4   eVec1, eVec2, eVec3;

  // Error statistics;
  int    nFew;

};

//==========================================================================

// Thrust class.
// This class performs thrust analysis on an event.

class Thrust {

public:

  // Constructor.
  Thrust(int selectIn = 2) : select(selectIn), nFew(0) {}

  // Analyze event.
  bool analyze(const Event& event, ostream& os = cout);

  // Return info on results of analysis.
  double thrust()       const {return eVal1;}
  double tMajor()       const {return eVal2;}
  double tMinor()       const {return eVal3;}
  double oblateness()   const {return eVal2 - eVal3;}
  Vec4 eventAxis(int i) const {return (i < 2) ? eVec1 :
    ( (i < 3) ? eVec2 : eVec3 ) ;}

  // Provide a listing of the info.
  void list(ostream& os = cout) const;

  // Tell how many events could not be analyzed.
  int nError() const {return nFew;}

private:

  // Constants: could only be changed in the code itself.
  static const int    NSTUDYMIN, TIMESTOPRINT;
  static const double MAJORMIN;

  // Properties of analysis.
  int    select;

  // Outcome of analysis.
  double eVal1, eVal2, eVal3;
  Vec4   eVec1, eVec2, eVec3;

  // Error statistics;
  int    nFew;

};

//==========================================================================

// SingleClusterJet class.
// Simple helper class to ClusterJet for a jet and its contents.

class SingleClusterJet {

public:

  // Constructors.
  SingleClusterJet(Vec4 pJetIn = 0., int motherIn = 0) :
    pJet(pJetIn), mother(motherIn), daughter(0), multiplicity(1),
    isAssigned(false) {pAbs = max( PABSMIN, pJet.pAbs());}
  SingleClusterJet& operator=(const SingleClusterJet& j) { if (this != &j)
    { pJet = j.pJet;  mother = j.mother; daughter = j.daughter;
    multiplicity = j.multiplicity; pAbs = j.pAbs;
    isAssigned = j.isAssigned;} return *this; }

  // Properties of jet.
  // Note: mother, daughter and isAssigned only used for original
  // particles, multiplicity and pTemp only for reconstructed jets.
  Vec4   pJet;
  int    mother, daughter, multiplicity;
  bool   isAssigned;
  double pAbs;
  Vec4   pTemp;

  // Distance measures (Lund, JADE, Durham) with friend.
  friend double dist2Fun(int measure, const SingleClusterJet& j1,
    const SingleClusterJet& j2);

private:

  // Constants: could only be changed in the code itself.
  static const double PABSMIN;

} ;

//--------------------------------------------------------------------------

// Namespace function declarations; friend of SingleClusterJet.

// Distance measures (Lund, JADE, Durham) with friend.
double dist2Fun(int measure, const SingleClusterJet& j1,
  const SingleClusterJet& j2);

//==========================================================================

// ClusterJet class.
// This class performs a jet clustering according to different
// distance measures: Lund, JADE or Durham.

class ClusterJet {

public:

  // Constructor.
  ClusterJet(string measureIn = "Lund", int selectIn = 2, int massSetIn = 2,
    bool preclusterIn = false, bool reassignIn = false) : measure(1),
    select(selectIn), massSet(massSetIn), doPrecluster(preclusterIn),
    doReassign(reassignIn), nFew(0) {
    char firstChar = toupper(measureIn[0]);
    if (firstChar == 'J') measure = 2;
    if (firstChar == 'D') measure = 3;
  }

  // Analyze event.
  bool analyze(const Event& event, double yScaleIn, double pTscaleIn,
    int nJetMinIn = 1, int nJetMaxIn = 0, ostream& os = cout);

  // Return info on jets produced.
  int  size()      const {return jets.size();}
  Vec4 p(int i)    const {return jets[i].pJet;}
  int  mult(int i) const {return jets[i].multiplicity;}

  // Return belonging of particle to one of the jets (-1 if none).
  int jetAssignment(int i) const {
    for (int iP = 0; iP < int(particles.size()); ++iP)
    if (particles[iP].mother == i) return particles[iP].daughter;
    return -1;}

  // Provide a listing of the info.
  void list(ostream& os = cout) const;

  // Return info on clustering values.
  int    distanceSize() const {return distances.size();}
  double distance(int i) const {
    return (i < distanceSize()) ? distances[i] : 0.; }

  // Tell how many events could not be analyzed.
  int nError() const {return nFew;}

private:

  // Constants: could only be changed in the code itself.
  static const int    TIMESTOPRINT;
  static const double PIMASS, PABSMIN, PRECLUSTERFRAC, PRECLUSTERSTEP;

  // Properties of analysis.
  int    measure, select, massSet;
  bool   doPrecluster, doReassign;
  double yScale, pTscale;
  int    nJetMin, nJetMax;

  // Temporary results.
  double dist2Join, dist2BigMin, distPre, dist2Pre;
  vector<SingleClusterJet> particles;
  int    nParticles;

  // Error statistics;
  int    nFew;

  // Member functions for some operations (for clarity).
  void precluster();
  void reassign();

  // Outcome of analysis: ET-ordered list of jets.
  vector<SingleClusterJet> jets;

  // Outcome of analysis: the distance values where the jets were merged.
  deque<double> distances;

};

//==========================================================================

// SingleCell class.
// Simple helper class to CellJet for a cell and its contents.

class SingleCell {

public:

  // Constructor.
  SingleCell(int iCellIn = 0, double etaCellIn = 0., double phiCellIn = 0.,
    double eTcellIn = 0., int multiplicityIn = 0) : iCell(iCellIn),
    etaCell(etaCellIn), phiCell(phiCellIn), eTcell(eTcellIn),
    multiplicity(multiplicityIn), canBeSeed(true), isUsed(false),
    isAssigned(false) {}

  // Properties of cell.
  int    iCell;
  double etaCell, phiCell, eTcell;
  int    multiplicity;
  bool   canBeSeed, isUsed, isAssigned;

} ;

//==========================================================================

// SingleCellJet class.
// Simple helper class to CellJet for a jet and its contents.

class SingleCellJet {

public:

  // Constructor.
  SingleCellJet(double eTjetIn = 0., double etaCenterIn = 0.,
    double phiCenterIn = 0., double etaWeightedIn = 0.,
    double phiWeightedIn = 0., int multiplicityIn = 0,
    Vec4 pMassiveIn = 0.) : eTjet(eTjetIn), etaCenter(etaCenterIn),
    phiCenter(phiCenterIn), etaWeighted(etaWeightedIn),
    phiWeighted(phiWeightedIn), multiplicity(multiplicityIn),
    pMassive(pMassiveIn) {}

  // Properties of jet.
  double eTjet, etaCenter, phiCenter, etaWeighted, phiWeighted;
  int    multiplicity;
  Vec4   pMassive;

} ;

//==========================================================================

// CellJet class.
// This class performs a cone jet search in (eta, phi, E_T) space.

class CellJet {

public:

  // Constructor.
  CellJet(double etaMaxIn = 5., int nEtaIn = 50, int nPhiIn = 32,
    int selectIn = 2, int smearIn = 0, double resolutionIn = 0.5,
    double upperCutIn = 2., double thresholdIn = 0., Rndm* rndmPtrIn = 0)
    : etaMax(etaMaxIn), nEta(nEtaIn), nPhi(nPhiIn), select(selectIn),
    smear(smearIn), resolution(resolutionIn), upperCut(upperCutIn),
    threshold(thresholdIn), nFew(0), rndmPtr(rndmPtrIn) { }

  // Analyze event.
  bool analyze(const Event& event, double eTjetMinIn = 20.,
    double coneRadiusIn = 0.7, double eTseedIn = 1.5, ostream& os = cout);

  // Return info on results of analysis.
  int    size()              const {return jets.size();}
  double eT(int i)           const {return jets[i].eTjet;}
  double etaCenter(int i)    const {return jets[i].etaCenter;}
  double phiCenter(int i)    const {return jets[i].phiCenter;}
  double etaWeighted(int i)  const {return jets[i].etaWeighted;}
  double phiWeighted(int i)  const {return jets[i].phiWeighted;}
  int    multiplicity(int i) const {return jets[i].multiplicity;}
  Vec4   pMassless(int i)    const {return jets[i].eTjet * Vec4(
           cos(jets[i].phiWeighted),  sin(jets[i].phiWeighted),
          sinh(jets[i].etaWeighted), cosh(jets[i].etaWeighted) );}
  Vec4   pMassive(int i)     const {return jets[i].pMassive;}
  double m(int i)            const {return jets[i].pMassive.mCalc();}

  // Provide a listing of the info.
  void list(ostream& os = cout) const;

  // Tell how many events could not be analyzed: so far never.
  int nError() const {return nFew;}

private:

  // Constants: could only be changed in the code itself.
  static const int    TIMESTOPRINT;

  // Properties of analysis.
  double etaMax;
  int    nEta, nPhi, select, smear;
  double resolution, upperCut, threshold;
  double eTjetMin, coneRadius, eTseed;

  // Error statistics;
  int    nFew;

  // Outcome of analysis: ET-ordered list of jets.
  vector<SingleCellJet> jets;

  // Pointer to the random number generator (needed for energy smearing).
  Rndm* rndmPtr;

};

//==========================================================================

// SlowJetHook class.
// Base class, used to derive your own class with your selection criteria.

class SlowJetHook {

public:

  // Destructor.
  virtual ~SlowJetHook() { }

  // Method to be overloaded.
  // It will be called for all final-state particles, one at a time, and
  // should return true if the particle should be analyzed, false if not.
  // The particle is in location iSel of the event record.
  // If you wish you can also modify the four-momentum and mass that will
  //  be used in the analysis, without affecting the event record itself,
  // by changing pSel and mSel. Remember to respect E^2 - p^2 = m^2.
  virtual bool include(int iSel, const Event& event, Vec4& pSel,
    double& mSel) = 0;

};

//==========================================================================

// SingleSlowJet class.
// Simple helper class to SlowJet for a jet and its contents.

class SingleSlowJet {

public:

  // Constructors.
  SingleSlowJet( Vec4 pIn = 0., double pT2In = 0., double yIn = 0.,
      double phiIn = 0., int idxIn = 0) : p(pIn), pT2(pT2In), y(yIn),
      phi(phiIn), mult(1) { idx.insert(idxIn); }
  SingleSlowJet(const SingleSlowJet& ssj) : p(ssj.p), pT2(ssj.pT2),
    y(ssj.y), phi(ssj.phi), mult(ssj.mult), idx(ssj.idx) { }
  SingleSlowJet& operator=(const SingleSlowJet& ssj) { if (this != &ssj)
    { p = ssj.p; pT2 = ssj.pT2; y = ssj.y; phi = ssj.phi;
    mult = ssj.mult; idx = ssj.idx; } return *this; }

  // Properties of jet.
  Vec4     p;
  double   pT2, y, phi;
  int      mult;
  set<int> idx;

};

//==========================================================================

// SlowJet class.
// This class performs a recombination jet search in (y, phi, pT) space.

class SlowJet {

public:

  // Constructor.
  SlowJet(int powerIn, double Rin, double pTjetMinIn = 0.,
    double etaMaxIn = 25., int selectIn = 2, int massSetIn = 2,
    SlowJetHook* sjHookPtrIn = 0, bool useFJcoreIn = true,
    bool useStandardRin = true) : power(powerIn), R(Rin),
    pTjetMin(pTjetMinIn), etaMax(etaMaxIn), select(selectIn),
    massSet(massSetIn), sjHookPtr(sjHookPtrIn), useFJcore(useFJcoreIn),
    useStandardR(useStandardRin) { isAnti = (power < 0); isKT = (power > 0);
    R2 = R*R; pT2jetMin = pTjetMin*pTjetMin; cutInEta = (etaMax <= 20.);
    chargedOnly = (select > 2); visibleOnly = (select == 2);
    modifyMass = (massSet < 2); noHook = (sjHookPtr == 0);}

  // Analyze event, all in one go.
  bool analyze(const Event& event) {
    if ( !setup(event) ) return false;
    if (useFJcore) return clusterFJ();
    while (clSize > 0) doStep(); return true; }

  // Set up list of particles to analyze, and initial distances.
  bool setup(const Event& event);

  // Do one recombination step, possibly giving a jet.
  bool doStep();

  // Do several recombinations steps, if possible.
  bool doNSteps(int nStep) { if (useFJcore) return false;
    while(nStep > 0 && clSize > 0) { doStep(); --nStep;}
    return (nStep == 0); }

  // Do recombinations until fixed numbers of clusters and jets remain.
  bool stopAtN(int nStop) { if (useFJcore) return false;
    while (clSize + jtSize > nStop && clSize > 0) doStep();
    return (clSize + jtSize == nStop); }

  // Return info on jet (+cluster) results of analysis.
  int    sizeOrig()          const {return origSize;}
  int    sizeJet()           const {return jtSize;}
  int    sizeAll()           const {return jtSize + clSize;}
  double pT(int i)           const {return (i < jtSize)
    ? sqrt(jets[i].pT2) : sqrt(clusters[i - jtSize].pT2);}
  double y(int i)            const {return (i < jtSize)
    ? jets[i].y : clusters[i - jtSize].y;}
  double phi(int i)          const {return (i < jtSize)
    ? jets[i].phi : clusters[i - jtSize].phi;}
  Vec4   p(int i)            const {return (i < jtSize)
    ? jets[i].p : clusters[i - jtSize].p;}
  double m(int i)            const {return (i < jtSize)
    ? jets[i].p.mCalc() : clusters[i - jtSize].p.mCalc();}
  int    multiplicity(int i) const {return (i < jtSize)
    ? jets[i].mult : clusters[i - jtSize].mult;}

  // Return info on next step to be taken.
  int    iNext() const {return (iMin == -1) ? -1 : iMin + jtSize;}
  int    jNext() const {return (jMin == -1) ? -1 : jMin + jtSize;}
  double dNext() const {return dMin;}

  // Provide a listing of the info.
  void list(bool listAll = false, ostream& os = cout) const;

  // Give a list of all particles in the jet.
  vector<int> constituents(int j) { vector<int> cTemp;
    for (set<int>::iterator idxTmp = jets[j].idx.begin();
       idxTmp != jets[j].idx.end(); ++idxTmp)
       cTemp.push_back( *idxTmp); return cTemp;
  }

  // Give a list of all particles in the cluster.
  vector<int> clusConstituents(int j) { vector<int> cTemp;
    for (set<int>::iterator idxTmp = clusters[j].idx.begin();
       idxTmp != clusters[j].idx.end(); ++idxTmp)
       cTemp.push_back( *idxTmp); return cTemp;
  }

  // Give the index of the jet that the particle i of the event record
  // belongs to. Returns -1 if particle i is not found in a jet.
  int jetAssignment(int i) {
    for (int j = 0; j < sizeJet(); ++j)
      if (jets[j].idx.find(i) != jets[j].idx.end()) return j;
    return -1;
  }

  // Remove a jet.
  void removeJet(int i) {
    if (i < 0 || i >= jtSize) return;
    jets.erase(jets.begin() + i);
    jtSize--;
  }

private:

  // Constants: could only be changed in the code itself.
  static const int    TIMESTOPRINT;
  static const double PIMASS, TINY;

  // Properties of analysis as such.
  int    power;
  double R, pTjetMin, etaMax, R2, pT2jetMin;
  int    select, massSet;
  // SlowJetHook can be used to tailor particle selection step.
  SlowJetHook* sjHookPtr;
  bool   useFJcore, useStandardR, isAnti, isKT, cutInEta, chargedOnly,
         visibleOnly, modifyMass, noHook;

  // Intermediate clustering objects and final jet objects.
  vector<SingleSlowJet> clusters;
  vector<SingleSlowJet> jets;

  // Intermediate clustering distances.
  vector<double> diB;
  vector<double> dij;

  // Other intermediate variables.
  int    origSize, clSize, clLast, jtSize, iMin, jMin;
  double dPhi, dijTemp, dMin;

  // Find next cluster pair to join.
  void findNext();

  // Use FJcore interface to perform clustering.
  bool clusterFJ();

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_Analysis_H

