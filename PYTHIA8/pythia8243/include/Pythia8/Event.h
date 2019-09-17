// Event.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the Particle and Event classes.
// Particle: information on an instance of a particle.
// Junction: information on a junction between three colours.
// Event: list of particles in the current event.

#ifndef Pythia8_Event_H
#define Pythia8_Event_H

#include "Pythia8/Basics.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Forward references to ParticleDataEntry and ResonanceWidths classes.
class ParticleDataEntry;
class ResonanceWidths;
class Event;

//==========================================================================

// Particle class.
// This class holds info on a particle in general.

class Particle {

public:

  // Constructors.
  Particle() : idSave(0), statusSave(0), mother1Save(0), mother2Save(0),
    daughter1Save(0), daughter2Save(0), colSave(0), acolSave(0),
    pSave(Vec4(0.,0.,0.,0.)), mSave(0.), scaleSave(0.), polSave(9.),
    hasVertexSave(false), vProdSave(Vec4(0.,0.,0.,0.)), tauSave(0.),
    pdePtr(0), evtPtr(0) { }
  Particle(int idIn, int statusIn = 0, int mother1In = 0,
    int mother2In = 0, int daughter1In = 0, int daughter2In = 0,
    int colIn = 0, int acolIn = 0, double pxIn = 0., double pyIn = 0.,
    double pzIn = 0., double eIn = 0., double mIn = 0.,
    double scaleIn = 0., double polIn = 9.)
    : idSave(idIn), statusSave(statusIn), mother1Save(mother1In),
    mother2Save(mother2In), daughter1Save(daughter1In),
    daughter2Save(daughter2In), colSave(colIn), acolSave(acolIn),
    pSave(Vec4(pxIn, pyIn, pzIn, eIn)), mSave(mIn), scaleSave(scaleIn),
    polSave(polIn), hasVertexSave(false), vProdSave(Vec4(0.,0.,0.,0.)),
    tauSave(0.), pdePtr(0), evtPtr(0) { }
  Particle(int idIn, int statusIn, int mother1In, int mother2In,
    int daughter1In, int daughter2In, int colIn, int acolIn,
    Vec4 pIn, double mIn = 0., double scaleIn = 0., double polIn = 9.)
    : idSave(idIn), statusSave(statusIn), mother1Save(mother1In),
    mother2Save(mother2In), daughter1Save(daughter1In),
    daughter2Save(daughter2In), colSave(colIn), acolSave(acolIn),
    pSave(pIn), mSave(mIn), scaleSave(scaleIn), polSave(polIn),
    hasVertexSave(false), vProdSave(Vec4(0.,0.,0.,0.)), tauSave(0.),
    pdePtr(0), evtPtr(0) { }
  Particle(const Particle& pt) : idSave(pt.idSave),
    statusSave(pt.statusSave), mother1Save(pt.mother1Save),
    mother2Save(pt.mother2Save), daughter1Save(pt.daughter1Save),
    daughter2Save(pt.daughter2Save), colSave(pt.colSave),
    acolSave(pt.acolSave), pSave(pt.pSave), mSave(pt.mSave),
    scaleSave(pt.scaleSave), polSave(pt.polSave),
    hasVertexSave(pt.hasVertexSave), vProdSave(pt.vProdSave),
    tauSave(pt.tauSave), pdePtr(pt.pdePtr), evtPtr(pt.evtPtr) { }
  Particle& operator=(const Particle& pt) {if (this != &pt) {
    idSave = pt.idSave; statusSave = pt.statusSave;
    mother1Save = pt.mother1Save; mother2Save = pt.mother2Save;
    daughter1Save = pt.daughter1Save; daughter2Save = pt.daughter2Save;
    colSave = pt.colSave; acolSave = pt.acolSave; pSave = pt.pSave;
    mSave = pt.mSave; scaleSave = pt.scaleSave; polSave = pt.polSave;
    hasVertexSave = pt.hasVertexSave; vProdSave = pt.vProdSave;
    tauSave = pt.tauSave; pdePtr = pt.pdePtr; evtPtr = pt.evtPtr; }
    return *this; }

  // Destructor.
  virtual ~Particle() {}

  // Member functions to set the Event and ParticleDataEntry pointers.
  void setEvtPtr(Event* evtPtrIn) { evtPtr = evtPtrIn; setPDEPtr();}
  void setPDEPtr(ParticleDataEntry* pdePtrIn = 0);

  // Member functions for input.
  void id(int idIn) {idSave = idIn; setPDEPtr();}
  void status(int statusIn) {statusSave = statusIn;}
  void statusPos() {statusSave = abs(statusSave);}
  void statusNeg() {statusSave = -abs(statusSave);}
  void statusCode(int statusIn) {statusSave =
    (statusSave > 0) ? abs(statusIn) : -abs(statusIn);}
  void mother1(int mother1In) {mother1Save = mother1In;}
  void mother2(int mother2In) {mother2Save = mother2In;}
  void mothers(int mother1In = 0, int mother2In = 0)
    {mother1Save = mother1In; mother2Save = mother2In;}
  void daughter1(int daughter1In) {daughter1Save = daughter1In;}
  void daughter2(int daughter2In) {daughter2Save = daughter2In;}
  void daughters(int daughter1In = 0, int daughter2In = 0)
    {daughter1Save = daughter1In; daughter2Save = daughter2In;}
  void col(int colIn) {colSave = colIn;}
  void acol(int acolIn) {acolSave = acolIn;}
  void cols(int colIn = 0,int acolIn = 0) {colSave = colIn;
    acolSave = acolIn;}
  void p(Vec4 pIn) {pSave = pIn;}
  void p(double pxIn, double pyIn, double pzIn, double eIn)
    {pSave.p(pxIn, pyIn, pzIn, eIn);}
  void px(double pxIn) {pSave.px(pxIn);}
  void py(double pyIn) {pSave.py(pyIn);}
  void pz(double pzIn) {pSave.pz(pzIn);}
  void e(double eIn) {pSave.e(eIn);}
  void m(double mIn) {mSave = mIn;}
  void scale(double scaleIn) {scaleSave = scaleIn;}
  void pol(double polIn) {polSave = polIn;}
  void vProd(Vec4 vProdIn) {vProdSave = vProdIn; hasVertexSave = true;}
  void vProd(double xProdIn, double yProdIn, double zProdIn, double tProdIn)
    {vProdSave.p(xProdIn, yProdIn, zProdIn, tProdIn); hasVertexSave = true;}
  void xProd(double xProdIn) {vProdSave.px(xProdIn); hasVertexSave = true;}
  void yProd(double yProdIn) {vProdSave.py(yProdIn); hasVertexSave = true;}
  void zProd(double zProdIn) {vProdSave.pz(zProdIn); hasVertexSave = true;}
  void tProd(double tProdIn) {vProdSave.e(tProdIn); hasVertexSave = true;}
  void vProdAdd(Vec4 vProdIn) {vProdSave += vProdIn; hasVertexSave = true;}
  void tau(double tauIn) {tauSave = tauIn;}

  // Member functions for output.
  int    id()        const {return idSave;}
  int    status()    const {return statusSave;}
  int    mother1()   const {return mother1Save;}
  int    mother2()   const {return mother2Save;}
  int    daughter1() const {return daughter1Save;}
  int    daughter2() const {return daughter2Save;}
  int    col()       const {return colSave;}
  int    acol()      const {return acolSave;}
  Vec4   p()         const {return pSave;}
  double px()        const {return pSave.px();}
  double py()        const {return pSave.py();}
  double pz()        const {return pSave.pz();}
  double e()         const {return pSave.e();}
  double m()         const {return mSave;}
  double scale()     const {return scaleSave;}
  double pol()       const {return polSave;}
  bool   hasVertex() const {return hasVertexSave;}
  Vec4   vProd()     const {return vProdSave;}
  double xProd()     const {return vProdSave.px();}
  double yProd()     const {return vProdSave.py();}
  double zProd()     const {return vProdSave.pz();}
  double tProd()     const {return vProdSave.e();}
  double tau()       const {return tauSave;}

  // Member functions for output; derived int and bool quantities.
  int    idAbs()     const {return abs(idSave);}
  int    statusAbs() const {return abs(statusSave);}
  bool   isFinal()   const {return (statusSave > 0);}
  int    intPol()    const;
  bool   isRescatteredIncoming() const {return
    (statusSave == -34 || statusSave == -45 ||
     statusSave == -46 || statusSave == -54);}

  // Member functions for output; derived double quantities.
  double m2()        const {return (mSave >= 0.) ?  mSave*mSave
                                                 : -mSave*mSave;}
  double mCalc()     const {return pSave.mCalc();}
  double m2Calc()    const {return pSave.m2Calc();}
  double eCalc()     const {return sqrt(abs(m2() + pSave.pAbs2()));}
  double pT()        const {return pSave.pT();}
  double pT2()       const {return pSave.pT2();}
  double mT()        const {double temp = m2() + pSave.pT2();
    return (temp >= 0.) ? sqrt(temp) : -sqrt(-temp);}
  double mT2()       const {return m2() + pSave.pT2();}
  double pAbs()      const {return pSave.pAbs();}
  double pAbs2()     const {return pSave.pAbs2();}
  double eT()        const {return pSave.eT();}
  double eT2()       const {return pSave.eT2();}
  double theta()     const {return pSave.theta();}
  double phi()       const {return pSave.phi();}
  double thetaXZ()   const {return pSave.thetaXZ();}
  double pPos()      const {return pSave.pPos();}
  double pNeg()      const {return pSave.pNeg();}
  double y()         const;
  double eta()       const;
  double y(double mCut) const;
  double y(double mCut, RotBstMatrix& M) const;
  Vec4   vDec()      const {return (tauSave > 0. && mSave > 0.)
    ? vProdSave + tauSave * pSave / mSave : vProdSave;}
  double xDec()      const {return (tauSave > 0. && mSave > 0.)
    ? vProdSave.px() + tauSave * pSave.px() / mSave : vProdSave.px();}
  double yDec()      const {return (tauSave > 0. && mSave > 0.)
    ? vProdSave.py() + tauSave * pSave.py() / mSave : vProdSave.py();}
  double zDec()      const {return (tauSave > 0. && mSave > 0.)
    ? vProdSave.pz() + tauSave * pSave.pz() / mSave : vProdSave.pz();}
  double tDec()      const {return (tauSave > 0. && mSave > 0.)
    ? vProdSave.e()  + tauSave * pSave.e()  / mSave : vProdSave.e();}

  // Methods that can refer back to the event the particle belongs to.
  virtual int index() const;
  int iTopCopy()      const;
  int iBotCopy()      const;
  int iTopCopyId(bool simplify = false) const;
  int iBotCopyId(bool simplify = false) const;
  vector<int> motherList()   const;
  vector<int> daughterList() const;
  vector<int> daughterListRecursive() const;
  vector<int> sisterList(bool traceTopBot = false) const;
  bool isAncestor(int iAncestor) const;
  int statusHepMC()  const;
  bool isFinalPartonLevel() const;
  bool undoDecay();

  // Further output, based on a pointer to a ParticleDataEntry object.
  string name()      const {
    return (pdePtr != 0) ? pdePtr->name(idSave) : " ";}
  string nameWithStatus(int maxLen = 20) const;
  int    spinType()  const {
    return (pdePtr != 0) ? pdePtr->spinType() : 0;}
  int    chargeType() const {
    return (pdePtr != 0) ? pdePtr->chargeType(idSave) : 0;}
  double charge()    const {
    return (pdePtr != 0) ?  pdePtr->charge(idSave) : 0;}
  bool   isCharged() const {
    return (pdePtr != 0) ? (pdePtr->chargeType(idSave) != 0) : false;}
  bool   isNeutral() const {
    return (pdePtr != 0) ? (pdePtr->chargeType(idSave) == 0) : false;}
  int    colType()   const {
    return (pdePtr != 0) ? pdePtr->colType(idSave) : 0;}
  double m0()        const {
    return (pdePtr != 0) ? pdePtr->m0() : 0.;}
  double mWidth()    const {
    return (pdePtr != 0) ? pdePtr->mWidth() : 0.;}
  double mMin()      const {
    return (pdePtr != 0) ? pdePtr->mMin() : 0.;}
  double mMax()      const {
    return (pdePtr != 0) ? pdePtr->mMax() : 0.;}
  double mSel()      const {
    return (pdePtr != 0) ? pdePtr->mSel() : 0.;}
  double constituentMass() const {
    return (pdePtr != 0) ? pdePtr->constituentMass() : 0.;}
  double tau0()      const {
    return (pdePtr != 0) ? pdePtr->tau0() : 0.;}
  bool   mayDecay()  const {
    return (pdePtr != 0) ? pdePtr->mayDecay() : false;}
  bool   canDecay()  const {
    return (pdePtr != 0) ? pdePtr->canDecay() : false;}
  bool   doExternalDecay() const {
    return (pdePtr != 0) ? pdePtr->doExternalDecay() : false;}
  bool   isResonance() const {
    return (pdePtr != 0) ? pdePtr->isResonance() : false;}
  bool   isVisible() const {
    return (pdePtr != 0) ? pdePtr->isVisible() : false;}
  bool   isLepton()  const {
    return (pdePtr != 0) ? pdePtr->isLepton() : false;}
  bool   isQuark()   const {
    return (pdePtr != 0) ? pdePtr->isQuark() : false;}
  bool   isGluon()   const {
    return (pdePtr != 0) ? pdePtr->isGluon() : false;}
  bool   isDiquark()   const {
    return (pdePtr != 0) ? pdePtr->isDiquark() : false;}
  bool   isParton()   const {
    return (pdePtr != 0) ? pdePtr->isParton() : false;}
  bool   isHadron()  const {
    return (pdePtr != 0) ? pdePtr->isHadron() : false;}
  ParticleDataEntry& particleDataEntry() const {return *pdePtr;}

  // Member functions that perform operations.
  void rescale3(double fac) {pSave.rescale3(fac);}
  void rescale4(double fac) {pSave.rescale4(fac);}
  void rescale5(double fac) {pSave.rescale4(fac); mSave *= fac;}
  void rot(double thetaIn, double phiIn) {pSave.rot(thetaIn, phiIn);
    if (hasVertexSave) vProdSave.rot(thetaIn, phiIn);}
  void bst(double betaX, double betaY, double betaZ) {
    pSave.bst(betaX, betaY, betaZ);
    if (hasVertexSave) vProdSave.bst(betaX, betaY, betaZ);}
  void bst(double betaX, double betaY, double betaZ, double gamma) {
    pSave.bst(betaX, betaY, betaZ, gamma);
    if (hasVertexSave) vProdSave.bst(betaX, betaY, betaZ, gamma);}
  void bst(const Vec4& pBst) {pSave.bst(pBst);
    if (hasVertexSave) vProdSave.bst(pBst);}
  void bst(const Vec4& pBst, double mBst) {pSave.bst(pBst, mBst);
    if (hasVertexSave) vProdSave.bst(pBst, mBst);}
  void bstback(const Vec4& pBst) {pSave.bstback(pBst);
    if (hasVertexSave) vProdSave.bstback(pBst);}
  void bstback(const Vec4& pBst, double mBst) {pSave.bstback(pBst, mBst);
    if (hasVertexSave) vProdSave.bstback(pBst, mBst);}
  void rotbst(const RotBstMatrix& M, bool boostVertex = true) {pSave.rotbst(M);
    if (hasVertexSave && boostVertex) vProdSave.rotbst(M);}
  void offsetHistory( int minMother, int addMother, int minDaughter,
    int addDaughter);
  void offsetCol( int addCol);

protected:

  // Constants: could only be changed in the code itself.
  static const double TINY;

  // Properties of the current particle.
  int    idSave, statusSave, mother1Save, mother2Save, daughter1Save,
         daughter2Save, colSave, acolSave;
  Vec4   pSave;
  double mSave, scaleSave, polSave;
  bool   hasVertexSave;
  Vec4   vProdSave;
  double tauSave;

  // Pointer to properties of the particle species.
  // Should no be saved in a persistent copy of the event record.
  // The //! below is ROOT notation that this member should not be saved.
  // Event::restorePtrs() can be called to restore the missing information.
  ParticleDataEntry* pdePtr;  //!

  // Pointer to the whole event record to which the particle belongs (if any).
  // As above it should not be saved.
  Event*             evtPtr;  //!

};

// Invariant mass of a pair and its square.
// (Not part of class proper, but tightly linked.)
double m(const Particle&, const Particle&);
double m2(const Particle&, const Particle&);

//==========================================================================

// The junction class stores what kind of junction it is, the colour indices
// of the legs at the junction and as far out as legs have been traced,
// and the status codes assigned for fragmentation of each leg.

class Junction {

public:

  // Constructors.
  Junction() : remainsSave(true), kindSave(0), colSave(), endColSave(),
    statusSave() { }

  Junction( int kindIn, int col0In, int col1In, int col2In)
    : remainsSave(true), kindSave(kindIn), colSave(), endColSave(),
    statusSave() {
      colSave[0] = col0In; colSave[1] = col1In; colSave[2] = col2In;
      for (int j = 0; j < 3; ++j) {
      endColSave[j] = colSave[j];  } }
  Junction(const Junction& ju) : remainsSave(ju.remainsSave),
    kindSave(ju.kindSave), colSave(), endColSave(), statusSave() {
    for (int j = 0; j < 3; ++j) {
    colSave[j] = ju.colSave[j]; endColSave[j] = ju.endColSave[j];
    statusSave[j] = ju.statusSave[j]; } }
  Junction& operator=(const Junction& ju) {if (this != &ju) {
    remainsSave = ju.remainsSave; kindSave =  ju.kindSave;
    for (int j = 0; j < 3; ++j) { colSave[j] = ju.colSave[j];
    endColSave[j] = ju.endColSave[j]; statusSave[j] = ju.statusSave[j]; } }
    return *this; }

  // Set values.
  void remains(bool remainsIn) {remainsSave = remainsIn;}
  void col(int j, int colIn) {colSave[j] = colIn; endColSave[j] = colIn;}
  void cols(int j, int colIn, int endColIn) {colSave[j] = colIn;
    endColSave[j] = endColIn;}
  void endCol(int j, int endColIn) {endColSave[j] = endColIn;}
  void status(int j, int statusIn) {statusSave[j] = statusIn;}

  // Read out value.
  bool   remains()     const {return remainsSave;}
  int    kind()        const {return kindSave;}
  int    col(int j)    const {return colSave[j];}
  int    endCol(int j) const {return endColSave[j];}
  int    status(int j) const {return statusSave[j];}

private:

  // Kind, positions of the three ends and their status codes.
  bool remainsSave;
  int kindSave, colSave[3], endColSave[3], statusSave[3];

};

//==========================================================================

// The Event class holds all info on the generated event.

class Event {

public:

  // Constructors.
  Event(int capacity = 100) : startColTag(100), maxColTag(100),
    savedSize(0), savedJunctionSize(0), savedPartonLevelSize(0),
    scaleSave(0.), scaleSecondSave(0.),
    headerList("----------------------------------------"),
    particleDataPtr(0) { entry.reserve(capacity); }
  Event& operator=(const Event& oldEvent);
  Event(const Event& oldEvent) {*this = oldEvent;}

  // Initialize header for event listing, particle data table, and colour.
  void init( string headerIn = "", ParticleData* particleDataPtrIn = 0,
    int startColTagIn = 100) {
    headerList.replace(0, headerIn.length() + 2, headerIn + "  ");
     particleDataPtr = particleDataPtrIn; startColTag = startColTagIn;}

  // Clear event record.
  void clear() {entry.resize(0); maxColTag = startColTag;
    savedPartonLevelSize = 0; scaleSave = 0.; scaleSecondSave = 0.;
    clearJunctions();}
  void free() {vector<Particle>().swap(entry); maxColTag = startColTag;
    savedPartonLevelSize = 0; scaleSave = 0.; scaleSecondSave = 0.;
    clearJunctions();}

  // Clear event record, and set first particle empty.
  void reset() {clear(); append(90, -11, 0, 0, 0., 0., 0., 0., 0.);}

  // Overload index operator to access element of event record.
  Particle& operator[](int i) {return entry.at(i);}
  const Particle& operator[](int i) const {return entry.at(i);}

  // Implement standard references to elements in the particle array.
  Particle& front()   {return entry.front();}
  Particle& at(int i) {return entry.at(i);}
  Particle& back()    {return entry.back();}

  // Event record size.
  int size() const {return entry.size();}

  // Put a new particle at the end of the event record; return index.
  int append(Particle entryIn) {
    entry.push_back(entryIn); setEvtPtr();
    if (entryIn.col() > maxColTag) maxColTag = entryIn.col();
    if (entryIn.acol() > maxColTag) maxColTag = entryIn.acol();
    return entry.size() - 1;
  }
  int append(int id, int status, int mother1, int mother2, int daughter1,
    int daughter2, int col, int acol, double px, double py, double pz,
    double e, double m = 0., double scaleIn = 0., double polIn = 9.) {
    entry.push_back( Particle(id, status, mother1, mother2, daughter1,
    daughter2, col, acol, px, py, pz, e, m, scaleIn, polIn) ); setEvtPtr();
    if (col > maxColTag) maxColTag = col;
    if (acol > maxColTag) maxColTag = acol;
    return entry.size() - 1;
  }
  int append(int id, int status, int mother1, int mother2, int daughter1,
    int daughter2, int col, int acol, Vec4 p, double m = 0.,
    double scaleIn = 0., double polIn = 9.) {
    entry.push_back( Particle(id, status, mother1, mother2, daughter1,
    daughter2, col, acol, p, m, scaleIn, polIn) ); setEvtPtr();
    if (col > maxColTag) maxColTag = col;
    if (acol > maxColTag) maxColTag = acol;
    return entry.size() - 1;
  }

  // Brief versions of append: no mothers and no daughters.
  int append(int id, int status, int col, int acol, double px, double py,
    double pz, double e, double m = 0., double scaleIn = 0.,
    double polIn = 9.) { entry.push_back( Particle(id, status, 0, 0, 0, 0,
    col, acol, px, py, pz, e, m, scaleIn, polIn) ); setEvtPtr();
    if (col > maxColTag) maxColTag = col;
    if (acol > maxColTag) maxColTag = acol;
    return entry.size() - 1;
  }
  int append(int id, int status, int col, int acol, Vec4 p, double m = 0.,
    double scaleIn = 0., double polIn = 9.) {entry.push_back( Particle(id,
    status, 0, 0, 0, 0, col, acol, p, m, scaleIn, polIn) ); setEvtPtr();
    if (col > maxColTag) maxColTag = col;
    if (acol > maxColTag) maxColTag = acol;
    return entry.size() - 1;
  }

  // Set pointer to the event for a particle, by default latest one.
  void setEvtPtr(int iSet = -1) {if (iSet < 0) iSet = entry.size() - 1;
    entry[iSet].setEvtPtr( this);}

  // Add a copy of an existing particle at the end of the event record.
  int copy(int iCopy, int newStatus = 0);

  // List the particles in an event.
  void list(bool showScaleAndVertex = false,
    bool showMothersAndDaughters = false, int precision = 3) const;

  // Remove last n entries.
  void popBack(int nRemove = 1) { if (nRemove ==1) entry.pop_back();
    else {int newSize = max( 0, size() - nRemove);
    entry.resize(newSize);} }

  // Remove entries from iFirst to iLast, including endpoints, anf fix history.
  // (To the extent possible; history pointers in removed range are zeroed.)
  void remove(int iFirst, int iLast, bool shiftHistory = true);

  // Restore all ParticleDataEntry* pointers in the Particle vector.
  // Useful when a persistent copy of the event record is read back in.
  void restorePtrs() { for (int i = 0; i < size(); ++i) setEvtPtr(i); }

  // Save or restore the size of the event record (throwing at the end).
  void saveSize() {savedSize = entry.size();}
  void restoreSize() {entry.resize(savedSize);}
  int  savedSizeValue() {return savedSize;}

  // Initialize and access colour tag information.
  void initColTag(int colTag = 0) {maxColTag = max( colTag,startColTag);}
  int lastColTag() const {return maxColTag;}
  int nextColTag() {return ++maxColTag;}

  // Access scale for which event as a whole is defined.
  void scale( double scaleIn) {scaleSave = scaleIn;}
  double scale() const {return scaleSave;}

  // Need a second scale if two hard interactions in event.
  void scaleSecond( double scaleSecondIn) {scaleSecondSave = scaleSecondIn;}
  double scaleSecond() const {return scaleSecondSave;}

  // Find complete list of daughters.
  // Note: temporarily retained for CMS compatibility. Do not use!
  vector<int> daughterList(int i) const {return entry[i].daughterList();}

  // Return number of final-state particles, optionally charged only.
  int nFinal(bool chargedOnly = false) const {
    int nFin = 0;
    for (int i = 0; i < size(); ++i)
      if (entry[i].isFinal() && (!chargedOnly || entry[i].isCharged()))
        ++nFin;
    return nFin; }

  // Find separation in y, eta, phi or R between two particles.
  double dyAbs(int i1, int i2) const {
    return abs( entry[i1].y() - entry[i2].y() ); }
  double detaAbs(int i1, int i2) const {
    return abs( entry[i1].eta() - entry[i2].eta() ); }
  double dphiAbs(int i1, int i2) const {
    double dPhiTmp = abs( entry[i1].phi() - entry[i2].phi() );
    if (dPhiTmp > M_PI)
      dPhiTmp = 2. * M_PI - dPhiTmp;
    return dPhiTmp; }
  double RRapPhi(int i1, int i2) const {
    return sqrt( pow2(dyAbs(i1, i2)) + pow2(dphiAbs(i1, i2)) ); }
  double REtaPhi(int i1, int i2) const {
    return sqrt( pow2(detaAbs(i1, i2)) + pow2(dphiAbs(i1, i2)) ); }

  // Member functions for rotations and boosts of an event.
  void rot(double theta, double phi)
    {for (int i = 0; i < size(); ++i) entry[i].rot(theta, phi);}
  void bst(double betaX, double betaY, double betaZ)
    {for (int i = 0; i < size(); ++i) entry[i].bst(betaX, betaY, betaZ);}
  void bst(double betaX, double betaY, double betaZ, double gamma)
    {for (int i = 0; i < size(); ++i) entry[i].bst(betaX, betaY, betaZ,
    gamma);}
  void bst(const Vec4& vec)
    {for (int i = 0; i < size(); ++i) entry[i].bst(vec);}
  void rotbst(const RotBstMatrix& M, bool boostVertices = true)
    {for (int i = 0; i < size(); ++i) entry[i].rotbst(M, boostVertices);}

  // Clear the list of junctions.
  void clearJunctions() {junction.resize(0);}

  // Add a junction to the list, study it or extra input.
  int appendJunction( int kind, int col0, int col1, int col2)
    { junction.push_back( Junction( kind, col0, col1, col2) );
    return junction.size() - 1;}
  int appendJunction(Junction junctionIn) {junction.push_back(junctionIn);
    return junction.size() - 1;}
  int sizeJunction() const {return junction.size();}
  bool remainsJunction(int i) const {return junction[i].remains();}
  void remainsJunction(int i, bool remainsIn) {junction[i].remains(remainsIn);}
  int kindJunction(int i) const {return junction[i].kind();}
  int colJunction( int i, int j) const {return junction[i].col(j);}
  void colJunction( int i, int j, int colIn) {junction[i].col(j, colIn);}
  int endColJunction( int i, int j) const {return junction[i].endCol(j);}
  void endColJunction( int i, int j, int endColIn)
    {junction[i].endCol(j, endColIn);}
  int statusJunction( int i, int j) const {return junction[i].status(j);}
  void statusJunction( int i, int j, int statusIn)
    {junction[i].status(j, statusIn);}
  Junction& getJunction(int i) {return junction[i];}
  const Junction& getJunction(int i) const {return junction[i];}
  void eraseJunction(int i);

  // Save or restore the size of the junction list (throwing at the end).
  void saveJunctionSize() {savedJunctionSize = junction.size();}
  void restoreJunctionSize() {junction.resize(savedJunctionSize);}

  // List any junctions in the event; for debug mainly.
  void listJunctions() const;

  // Save event record size at Parton Level, i.e. before hadronization.
  void savePartonLevelSize() {savedPartonLevelSize = entry.size();}

  // Operator overloading allows to append one event to an existing one.
  // Warning: particles should be OK, but some other information unreliable.
  Event& operator+=(const Event& addEvent);

private:

  // The Particle class needs to access particle data.
  friend class Particle;

  // Constants: could only be changed in the code itself.
  static const int IPERLINE;

  // Initialization data, normally only set once.
  int startColTag;

  // The event: a vector containing all particles (entries).
  // The explicit use of Pythia8:: qualifier patches a limitation in ROOT.
  vector<Pythia8::Particle> entry;

  // The list of junctions.
  // The explicit use of Pythia8:: qualifier patches a limitation in ROOT.
  vector<Pythia8::Junction> junction;

  // The maximum colour tag of the event so far.
  int maxColTag;

  // Saved entry and junction list sizes, for simple restoration.
  int savedSize, savedJunctionSize, savedPartonLevelSize;

  // The scale of the event; linear quantity in GeV.
  double scaleSave, scaleSecondSave;

  // Header specification in event listing (at most 40 characters wide).
  string headerList;

  // Pointer to the particle data table.
  // The //! below is ROOT notation that this member should not be saved.
  ParticleData* particleDataPtr;  //!

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_Event_H
