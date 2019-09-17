// ParticleData.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the classes containing particle data.
// DecayChannel contains info on a single decay channel.
// ParticleDataEntry contains info on a single particle species.
// ParticleData collects info on all particles as a map.

#ifndef Pythia8_ParticleData_H
#define Pythia8_ParticleData_H

#include "Pythia8/Basics.h"
#include "Pythia8/Info.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"

namespace Pythia8 {

//==========================================================================

// Forward reference to some classes.
class ParticleData;
class ResonanceWidths;
class Couplings;
class CoupSUSY;
class SUSYResonanceWidths;

//==========================================================================

// This class holds info on a single decay channel.

class DecayChannel {

public:

  // Constructor.
  DecayChannel(int onModeIn = 0, double bRatioIn = 0., int meModeIn = 0,
    int prod0 = 0, int prod1 = 0, int prod2 = 0, int prod3 = 0,
    int prod4 = 0, int prod5 = 0, int prod6 = 0, int prod7 = 0)
    : onModeSave(onModeIn), bRatioSave(bRatioIn), currentBRSave(0.),
    onShellWidthSave(0.), openSecPos(1.), openSecNeg(1.),
    meModeSave(meModeIn), nProd(0), prod(), hasChangedSave(true) {
    prod[0] = prod0; prod[1] = prod1; prod[2] = prod2; prod[3] = prod3;
    prod[4] = prod4; prod[5] = prod5; prod[6] = prod6; prod[7] = prod7;
    for (int j = 0; j < 8; ++j) if (prod[j] != 0 && j == nProd) ++nProd; }

  // Copy constructor.
  DecayChannel& operator=( const DecayChannel& oldDC) { if (this != &oldDC) {
    onModeSave = oldDC.onModeSave; bRatioSave = oldDC.bRatioSave;
    currentBRSave = oldDC.currentBRSave;
    onShellWidthSave = oldDC.onShellWidthSave; openSecPos = oldDC.openSecPos;
    openSecNeg = oldDC.openSecNeg; meModeSave = oldDC.meModeSave;
    nProd = oldDC.nProd; for (int j = 0; j < 8; ++j) prod[j] = oldDC.prod[j];
    hasChangedSave = oldDC.hasChangedSave; } return *this; }

  // Member functions for input.
  void onMode(int onModeIn) {onModeSave = onModeIn; hasChangedSave = true;}
  void bRatio(double bRatioIn, bool countAsChanged = true) {
    bRatioSave = bRatioIn; if (countAsChanged) hasChangedSave = true;}
  void rescaleBR(double fac) {bRatioSave *= fac; hasChangedSave = true;}
  void meMode(int meModeIn) {meModeSave = meModeIn; hasChangedSave = true;}
  void multiplicity(int multIn)  {nProd = multIn; hasChangedSave = true;}
  void product(int i, int prodIn) {prod[i] = prodIn; nProd = 0;
    for (int j = 0; j < 8; ++j) if (prod[j] != 0 && j == nProd) ++nProd;
    hasChangedSave = true;}
  void setHasChanged(bool hasChangedIn) {hasChangedSave = hasChangedIn;}

  // Member functions for output.
  int    onMode()       const {return onModeSave;}
  double bRatio()       const {return bRatioSave;}
  int    meMode()       const {return meModeSave;}
  int    multiplicity() const {return nProd;}
  int    product(int i) const {return (i >= 0 && i < nProd) ? prod[i] : 0;}
  bool   hasChanged()   const { return hasChangedSave;}

  // Check for presence of particles anywhere in decay list.
  bool   contains(int id1) const;
  bool   contains(int id1, int id2) const;
  bool   contains(int id1, int id2, int id3) const;

  // Input/output for current selection of decay modes.
  // Takes into account on/off switches and dynamic width for resonances.
  void   currentBR(double currentBRIn) {currentBRSave = currentBRIn;}
  double currentBR() const {return currentBRSave;}

  // Input/output for nominal partial width; used by resonances.
  void   onShellWidth(double onShellWidthIn) {
         onShellWidthSave = onShellWidthIn;}
  double onShellWidth() const {return onShellWidthSave;}
  void   onShellWidthFactor(double factor) {onShellWidthSave *= factor;}

  // Input/output for fraction of secondary open widths; used by resonances.
  void   openSec(int idSgn, double openSecIn) {
    if (idSgn > 0) openSecPos = openSecIn; else openSecNeg = openSecIn;}
  double openSec(int idSgn) const {
    return (idSgn > 0) ? openSecPos : openSecNeg;}

private:

  // Decay channel info.
  int    onModeSave;
  double bRatioSave, currentBRSave, onShellWidthSave, openSecPos,
         openSecNeg;
  int    meModeSave, nProd, prod[8];
  bool   hasChangedSave;

};

//==========================================================================

// This class holds info on a single particle species.

class ParticleDataEntry {

public:

  // Constructors: for antiparticle exists or not.
  ParticleDataEntry(int idIn = 0, string nameIn = " ",
    int spinTypeIn = 0, int chargeTypeIn = 0, int colTypeIn = 0,
    double m0In = 0., double mWidthIn = 0., double mMinIn = 0.,
    double mMaxIn = 0., double tau0In = 0.) : idSave(abs(idIn)),
    nameSave(nameIn), antiNameSave("void"),  spinTypeSave(spinTypeIn),
    chargeTypeSave(chargeTypeIn), colTypeSave(colTypeIn), m0Save(m0In),
    mWidthSave (mWidthIn), mMinSave(mMinIn), mMaxSave(mMaxIn),
    tau0Save(tau0In), constituentMassSave(), hasAntiSave(false),
    isResonanceSave(), mayDecaySave(), doExternalDecaySave(), isVisibleSave(),
    doForceWidthSave(), hasChangedSave(true), hasChangedMMinSave(false),
    hasChangedMMaxSave(false), modeBWnow(), modeTau0now(), atanLow(),
    atanDif(), mThr(), currentBRSum(), resonancePtr(0), particleDataPtr() {
    setDefaults();}
  ParticleDataEntry(int idIn, string nameIn, string antiNameIn,
    int spinTypeIn = 0, int chargeTypeIn = 0, int colTypeIn = 0,
    double m0In = 0., double mWidthIn = 0., double mMinIn = 0.,
    double mMaxIn = 0., double tau0In = 0.) : idSave(abs(idIn)),
    nameSave(nameIn), antiNameSave(antiNameIn), spinTypeSave(spinTypeIn),
    chargeTypeSave(chargeTypeIn), colTypeSave(colTypeIn), m0Save(m0In),
    mWidthSave (mWidthIn), mMinSave(mMinIn), mMaxSave(mMaxIn),
    tau0Save(tau0In), constituentMassSave(), hasAntiSave(true),
    isResonanceSave(), mayDecaySave(), doExternalDecaySave(), isVisibleSave(),
    doForceWidthSave(), hasChangedSave(true), hasChangedMMinSave(false),
    hasChangedMMaxSave(false), modeBWnow(), modeTau0now(), atanLow(),
    atanDif(), mThr(), currentBRSum(), resonancePtr(0), particleDataPtr() {
    setDefaults(); if (toLower(antiNameIn) == "void") hasAntiSave = false;}

  // Copy constructor.
  ParticleDataEntry& operator=( const ParticleDataEntry& oldPDE) {
    if (this != &oldPDE) { idSave = oldPDE.idSave;
    nameSave = oldPDE.nameSave; antiNameSave = oldPDE.antiNameSave;
    spinTypeSave = oldPDE.spinTypeSave; chargeTypeSave = oldPDE.chargeTypeSave;
    colTypeSave = oldPDE.colTypeSave; m0Save = oldPDE.m0Save;
    mWidthSave = oldPDE.mWidthSave;  mMinSave = oldPDE.mMinSave;
    mMaxSave = oldPDE.mMaxSave;  tau0Save = oldPDE.tau0Save;
    constituentMassSave = oldPDE.constituentMassSave;
    hasAntiSave = oldPDE.hasAntiSave; isResonanceSave = oldPDE.isResonanceSave;
    mayDecaySave = oldPDE.mayDecaySave; doExternalDecaySave
    = oldPDE.doExternalDecaySave; isVisibleSave = oldPDE.isVisibleSave;
    doForceWidthSave = oldPDE.doForceWidthSave; hasChangedSave
    = oldPDE.hasChangedSave; hasChangedMMinSave = oldPDE.hasChangedMMinSave;
    hasChangedMMaxSave = oldPDE.hasChangedMMaxSave;
    modeBWnow = oldPDE.modeBWnow; atanLow = oldPDE.atanLow;
    atanDif = oldPDE.atanDif; mThr = oldPDE.mThr;
    for (int i = 0; i < int(oldPDE.channels.size()); ++i) {
      DecayChannel oldDC = oldPDE.channels[i]; channels.push_back(oldDC); }
    currentBRSum = oldPDE.currentBRSum; resonancePtr = 0;
    particleDataPtr = 0; } return *this; }

  // Destructor: delete any ResonanceWidths object.
  ~ParticleDataEntry();

  // Initialization of some particle flags.
  void setDefaults();

  // Store pointer to whole particle data table/database.
  void initPtr( ParticleData* particleDataPtrIn) {
    particleDataPtr = particleDataPtrIn;}

  // Reset all the properties of an existing particle.
  void setAll(string nameIn, string antiNameIn, int spinTypeIn = 0,
    int chargeTypeIn = 0, int colTypeIn = 0, double m0In = 0.,
    double mWidthIn = 0., double mMinIn = 0., double mMaxIn = 0.,
    double tau0In = 0.)
    {nameSave = nameIn; antiNameSave = antiNameIn; hasAntiSave = true;
    if (toLower(antiNameIn) == "void") hasAntiSave = false;
    spinTypeSave = spinTypeIn; chargeTypeSave = chargeTypeIn;
    colTypeSave = colTypeIn; m0Save = m0In; mWidthSave = mWidthIn;
    setMMin(mMinIn); setMMax(mMaxIn); tau0Save = tau0In;
    setDefaults(); hasChangedSave = true;}

  // Change current values one at a time (or set if not set before).
  // (Must use set here since else name+signature clash with get methods.)
  void setName(string nameIn) {nameSave = nameIn; hasChangedSave = true;}
  void setAntiName(string antiNameIn) {antiNameSave = antiNameIn;
    hasAntiSave = (toLower(antiNameIn) != "void"); hasChangedSave = true;}
  void setNames(string nameIn, string antiNameIn) {nameSave = nameIn;
    antiNameSave = antiNameIn; hasAntiSave = (toLower(antiNameIn) != "void");
    hasChangedSave = true;}
  void setSpinType(int spinTypeIn) {spinTypeSave = spinTypeIn;
    hasChangedSave = true;}
  void setChargeType(int chargeTypeIn) {chargeTypeSave = chargeTypeIn;
    hasChangedSave = true;}
  void setColType(int colTypeIn) {colTypeSave = colTypeIn;
    hasChangedSave = true;}
  void setM0(double m0In) {m0Save = m0In; setConstituentMass();
    hasChangedSave = true;}
  void setMWidth(double mWidthIn, bool countAsChanged = true) {
    mWidthSave = mWidthIn; if (countAsChanged) hasChangedSave = true;}
  void setMMin(double mMinIn) {mMinSave = mMinIn; hasChangedSave = true;
    hasChangedMMinSave=true;}
  void setMMax(double mMaxIn) {mMaxSave = mMaxIn; hasChangedSave = true;
    hasChangedMMaxSave=true;}
  // Special options specifically when cutting wings of Breit-Wigners.
  void setMMinNoChange(double mMinIn) {mMinSave = mMinIn;}
  void setMMaxNoChange(double mMaxIn) {mMaxSave = mMaxIn;}
  void setTau0(double tau0In, bool countAsChanged = true)
    {tau0Save = tau0In; if (countAsChanged) hasChangedSave = true;}
  void setIsResonance(bool isResonanceIn) {isResonanceSave = isResonanceIn;
    hasChangedSave = true;}
  void setMayDecay(bool mayDecayIn, bool countAsChanged = true) {
    mayDecaySave = mayDecayIn; if (countAsChanged) hasChangedSave = true;}
  void setDoExternalDecay(bool doExternalDecayIn)
    {doExternalDecaySave = doExternalDecayIn; hasChangedSave = true;}
  void setIsVisible(bool isVisibleIn) {isVisibleSave = isVisibleIn;
    hasChangedSave = true;}
  void setDoForceWidth(bool doForceWidthIn) {doForceWidthSave = doForceWidthIn;
    hasChangedSave = true;}
  void setHasChanged(bool hasChangedIn) {hasChangedSave = hasChangedIn;
    for (int i = 0; i < int(channels.size()); ++i)
      channels[i].setHasChanged(hasChangedIn);
    if (!hasChangedIn) {hasChangedMMinSave=false; hasChangedMMaxSave=false;}}

  // Give back current values.
  int    id()                     const { return idSave; }
  bool   hasAnti()                const { return hasAntiSave; }
  string name(int idIn = 1)       const {
         return (idIn > 0) ? nameSave : antiNameSave; }
  int    spinType() const {return spinTypeSave; }
  int    chargeType(int idIn = 1) const {
         return (idIn > 0) ? chargeTypeSave : -chargeTypeSave; }
  double charge(int idIn = 1)     const {
         return (idIn > 0) ? chargeTypeSave / 3. : -chargeTypeSave / 3.; }
  int    colType(int idIn = 1)    const {
         if (colTypeSave == 2) return colTypeSave;
         return (idIn > 0) ? colTypeSave : -colTypeSave; }
  double m0()                     const { return m0Save; }
  double mWidth()                 const { return mWidthSave; }
  double mMin()                   const { return mMinSave; }
  double mMax()                   const { return mMaxSave; }
  double m0Min()                  const {
         return (modeBWnow == 0) ? m0Save : mMinSave; }
  double m0Max()                  const {
         return (modeBWnow == 0) ? m0Save : mMaxSave; }
  double tau0()                   const { return tau0Save; }
  bool   isResonance()            const { return isResonanceSave; }
  bool   mayDecay()               const { return mayDecaySave; }
  bool   doExternalDecay()        const { return doExternalDecaySave; }
  bool   isVisible()              const { return isVisibleSave; }
  bool   doForceWidth()           const { return doForceWidthSave; }
  bool   hasChanged()     const { if (hasChangedSave) return true;
         for (int i = 0; i < int(channels.size()); ++i)
           if (channels[i].hasChanged()) return true;
         return false;}
  bool   hasChangedMMin()         const { return hasChangedMMinSave; }
  bool   hasChangedMMax()         const { return hasChangedMMaxSave; }

  // Set and give back several mass-related quantities.
  void   initBWmass();
  double constituentMass()        const { return constituentMassSave; }
  double mSel();
  double mRun(double mH);

  // Give back other quantities.
  bool   useBreitWigner() const { return (modeBWnow > 0); }
  bool   canDecay()       const { return (channels.size() > 0);}
  bool   isLepton()       const { return (idSave > 10 && idSave < 19);}
  bool   isQuark()        const { return (idSave != 0 && idSave < 9);}
  bool   isGluon()        const { return (idSave == 21);}
  bool   isDiquark()      const { return (idSave > 1000 && idSave < 10000
         && (idSave/10)%10 == 0);}
  bool   isParton()       const { return ( idSave == 21
    || (idSave != 0 && idSave < 6)
    || (idSave > 1000 && idSave < 5510 && (idSave/10)%10 == 0) );}
  bool   isHadron()       const;
  bool   isMeson()        const;
  bool   isBaryon()       const;
  bool   isOnium()        const;

  // Intermediate octet ccbar or bbar states in colour-octet model.
  bool   isOctetHadron()  const {return idSave >= 9940000
      && idSave < 9960000; }
  int    heaviestQuark(int idIn = 1)    const;
  int    baryonNumberType(int idIn = 1) const;
  int    nQuarksInCode(int idQIn)       const;

  // Reset to empty decay table.
  void clearChannels() {channels.resize(0);}

  // Add a decay channel to the decay table.
  void addChannel(int onMode = 0, double bRatio = 0., int meMode = 0,
    int prod0 = 0, int prod1 = 0, int prod2 = 0, int prod3 = 0,
    int prod4 = 0, int prod5 = 0, int prod6 = 0, int prod7 = 0) {
    channels.push_back( DecayChannel( onMode, bRatio, meMode, prod0,
    prod1, prod2, prod3, prod4, prod5, prod6, prod7) ); }

  // Decay table size.
  int sizeChannels() const {return channels.size();}

  // Gain access to a channel in the decay table.
  DecayChannel& channel(int i){return channels[i];}
  const DecayChannel& channel(int i) const {return channels[i];}

  // Rescale sum of branching ratios to unity.
  void rescaleBR(double newSumBR = 1.);

  // Random choice of decay channel according to branching ratios.
  bool preparePick(int idSgn, double mHat = 0., int idInFlav = 0);
  DecayChannel& pickChannel();

  // Access methods stored in ResonanceWidths.
  void   setResonancePtr(ResonanceWidths* resonancePtrIn);
  ResonanceWidths* getResonancePtr() {return resonancePtr;}
  void   resInit(Info* infoPtrIn, Settings* settingsPtrIn,
    ParticleData* particleDataPtrIn, Couplings* couplingsPtrIn);
  double resWidth(int idSgn, double mHat, int idIn = 0,
    bool openOnly = false, bool setBR = false);
  double resWidthOpen(int idSgn, double mHat, int idIn = 0);
  double resWidthStore(int idSgn, double mHat, int idIn = 0);
  double resOpenFrac(int idSgn);
  double resWidthRescaleFactor();
  double resWidthChan(double mHat, int idAbs1 = 0, int idAbs2 = 0);

private:

  // Constants: could only be changed in the code itself.
  static const int    INVISIBLENUMBER, INVISIBLETABLE[80], KNOWNNOWIDTH[3];
  static const double MAXTAU0FORDECAY,MINMASSRESONANCE, NARROWMASS,
                      CONSTITUENTMASSTABLE[10];

  // Particle data.
  int    idSave;
  string nameSave, antiNameSave;
  int    spinTypeSave, chargeTypeSave, colTypeSave;
  double m0Save, mWidthSave, mMinSave, mMaxSave, tau0Save,
         constituentMassSave;
  bool   hasAntiSave, isResonanceSave, mayDecaySave, doExternalDecaySave,
         isVisibleSave, doForceWidthSave, hasChangedSave, hasChangedMMinSave,
         hasChangedMMaxSave;

  // Extra data for mass selection according to a Breit-Wigner and lifetime.
  int    modeBWnow, modeTau0now;
  double atanLow, atanDif, mThr;

  // A vector containing all the decay channels of the particle.
  vector<DecayChannel> channels;

  // Summed branching ratio of currently open channels.
  double currentBRSum;

  // Pointer to ResonanceWidths object; only used for some particles.
  ResonanceWidths* resonancePtr;

  // Pointer to the full particle data table.
  ParticleData* particleDataPtr;

  // Set constituent mass.
  void setConstituentMass();

};

//==========================================================================

// This class holds a map of all ParticleDataEntries.

class ParticleData {

public:

  // Constructor.
  ParticleData() : setRapidDecayVertex(), modeBreitWigner(), maxEnhanceBW(),
    mQRun(), Lambda5Run(), intermediateTau0(), infoPtr(0), settingsPtr(0),
    rndmPtr(0), couplingsPtr(0), particlePtr(0), isInit(false),
    readingFailedSave(false) {}

  // Copy constructors.
  ParticleData& operator=( const ParticleData& oldPD) { if (this != &oldPD) {
    modeBreitWigner = oldPD.modeBreitWigner; maxEnhanceBW = oldPD.maxEnhanceBW;
    for (int i = 0; i < 7; ++i) mQRun[i] = oldPD.mQRun[i];
    Lambda5Run = oldPD.Lambda5Run;
    infoPtr = 0; settingsPtr = 0; rndmPtr = 0; couplingsPtr = 0;
    for ( map<int, ParticleDataEntry>::const_iterator pde = oldPD.pdt.begin();
      pde != oldPD.pdt.end(); pde++) { int idTmp = pde->first;
      pdt[idTmp] = pde->second; pdt[idTmp].initPtr(this); }
    particlePtr = 0; isInit = oldPD.isInit;
    readingFailedSave = oldPD.readingFailedSave; } return *this; }

  // Initialize pointers.
  void initPtr(Info* infoPtrIn, Settings* settingsPtrIn, Rndm* rndmPtrIn,
    Couplings* couplingsPtrIn) {infoPtr = infoPtrIn;
    settingsPtr = settingsPtrIn; rndmPtr = rndmPtrIn;
    couplingsPtr = couplingsPtrIn;}

  // Read in database from specific file.
  bool init(string startFile = "../share/Pythia8/xmldoc/ParticleData.xml") {
    initCommon(); return readXML(startFile);}

  // Read in database from saved file stored in memory.
  bool init(const ParticleData& particleDataIn) {
    initCommon(); return copyXML(particleDataIn);}

  // Read in database from an istream.
  bool init(istream& is) { initCommon(); return readXML(is);}

  // Overwrite existing database by reading from specific file.
  bool reInit(string startFile, bool xmlFormat = true) { initCommon();
    return (xmlFormat) ? readXML(startFile) : readFF(startFile);}

  // Initialize pointers, normal Breit-Wigners and special resonances.
  void initWidths(vector<ResonanceWidths*> resonancePtrs);

  // Read and process or list whole (or part of) database from/to an XML file.
  bool readXML(string inFile, bool reset = true) ;
  void listXML(string outFile);
  bool readXML(istream& is, bool reset=true);

  // Copy and process XML information from another particleData object.
  bool copyXML(const ParticleData &particleDataIn);

  // Auxiliary functions to readXML() and copyXML().
  bool loadXML(string inFile, bool reset = true) ;
  bool loadXML(istream& is, bool reset=true);
  bool processXML(bool reset = true) ;

  // Read or list whole (or part of) database from/to a free format file.
  bool readFF(string inFile, bool reset = true) ;
  bool readFF(istream& is, bool reset = true);
  void listFF(string outFile);

  // Read in one update from a single line.
  bool readString(string lineIn, bool warn = true) ;

  // Keep track whether any readings have failed, invalidating run setup.
  bool readingFailed() {return readingFailedSave;}

  // Print out table of whole database, or of only part of it.
  void listAll() {list(false, true);}
  void listChanged(bool changedRes = false) {list(true, changedRes);}
  void list(bool changedOnly = false, bool changedRes = true);

  // Print out specified particles.
  void list(int idList) {vector<int> idListTemp;
    idListTemp.push_back(idList); list( idListTemp);}
  void list(vector<int> idList);

  // Retrieve readString history (e.g., for inspection). Everything
  // (subrun=-999), up to first subrun (=-1), or subrun-specific (>=0).
  vector<string> getReadHistory(int subrun=-999) {
    if (subrun == -999) return readStringHistory;
    else if (readStringSubrun.find(subrun) != readStringSubrun.end())
      return readStringSubrun[subrun];
    else return vector<string>();
  }

  // Check that table makes sense, especially for decays.
  void checkTable(int verbosity = 1) ;

  // Add new entry.
  void addParticle(int idIn, string nameIn = " ", int spinTypeIn = 0,
    int chargeTypeIn = 0, int colTypeIn = 0, double m0In = 0.,
    double mWidthIn = 0., double mMinIn = 0., double mMaxIn = 0.,
    double tau0In = 0.) { pdt[abs(idIn)] = ParticleDataEntry(idIn,
    nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn,
    mMinIn, mMaxIn, tau0In); pdt[abs(idIn)].initPtr(this); }
  void addParticle(int idIn, string nameIn, string antiNameIn,
    int spinTypeIn = 0, int chargeTypeIn = 0, int colTypeIn = 0,
    double m0In = 0., double mWidthIn = 0., double mMinIn = 0.,
    double mMaxIn = 0., double tau0In = 0.) { pdt[abs(idIn)]
    = ParticleDataEntry(idIn, nameIn, antiNameIn, spinTypeIn,
    chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn, tau0In);
    pdt[abs(idIn)].initPtr(this); }

  // Reset all the properties of an entry in one go.
  void setAll(int idIn, string nameIn, string antiNameIn,
    int spinTypeIn = 0, int chargeTypeIn = 0, int colTypeIn = 0,
    double m0In = 0., double mWidthIn = 0., double mMinIn = 0.,
    double mMaxIn = 0.,double tau0In = 0.) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setAll( nameIn, antiNameIn, spinTypeIn, chargeTypeIn,
    colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn, tau0In); }

  // Query existence of an entry.
  bool isParticle(int idIn) const {
    map<int,ParticleDataEntry>::const_iterator found = pdt.find( abs(idIn) );
    if ( found == pdt.end() ) return false;
    if ( idIn > 0 || found->second.hasAnti() ) return true;
    return false;
  }

  // Query existence of an entry and return an iterator.
  ParticleDataEntry* findParticle(int idIn) {
    map<int,ParticleDataEntry>::iterator found = pdt.find( abs(idIn) );
    if( found == pdt.end() ) return NULL;
    if ( idIn > 0 || found->second.hasAnti() ) return &((*found).second);
    return NULL;
  }

  // Query existence of an entry and return a const iterator.
  const ParticleDataEntry* findParticle(int idIn) const {
    map<int,ParticleDataEntry>::const_iterator found = pdt.find( abs(idIn) );
    if( found == pdt.end() ) return NULL;
    if ( idIn > 0 || found->second.hasAnti() ) return &((*found).second);
    return NULL;
  }

  // Return the id of the sequentially next particle stored in table.
  int nextId(int idIn) ;

  // Change current values one at a time (or set if not set before).
  void name(int idIn, string nameIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setName(nameIn); }
  void antiName(int idIn, string antiNameIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setAntiName(antiNameIn); }
  void names(int idIn, string nameIn, string antiNameIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setNames(nameIn, antiNameIn); }
  void spinType(int idIn, int spinTypeIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setSpinType(spinTypeIn); }
  void chargeType(int idIn, int chargeTypeIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setChargeType(chargeTypeIn); }
  void colType(int idIn, int colTypeIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setColType(colTypeIn); }
  void m0(int idIn, double m0In) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setM0(m0In); }
  void mWidth(int idIn, double mWidthIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setMWidth(mWidthIn); }
  void mMin(int idIn, double mMinIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setMMin(mMinIn); }
  void mMax(int idIn, double mMaxIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setMMax(mMaxIn); }
  void tau0(int idIn, double tau0In) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setTau0(tau0In); }
  void isResonance(int idIn, bool isResonanceIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setIsResonance(isResonanceIn); }
  void mayDecay(int idIn, bool mayDecayIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setMayDecay(mayDecayIn); }
  void doExternalDecay(int idIn, bool doExternalDecayIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setDoExternalDecay(doExternalDecayIn); }
  void isVisible(int idIn, bool isVisibleIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setIsVisible(isVisibleIn); }
  void doForceWidth(int idIn, bool doForceWidthIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setDoForceWidth(doForceWidthIn); }
  void hasChanged(int idIn, bool hasChangedIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setHasChanged(hasChangedIn); }

  // Give back current values.
  bool hasAnti(int idIn) const {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->hasAnti() : false; }
  string name(int idIn) const {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->name(idIn) : " "; }
  int spinType(int idIn) const {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->spinType() : 0; }
  int chargeType(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->chargeType(idIn) : 0; }
  double charge(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->charge(idIn) : 0; }
  int colType(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->colType(idIn) : 0 ; }
  double m0(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->m0() : 0. ; }
  double mWidth(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->mWidth() : 0. ; }
  double mMin(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->mMin() : 0. ; }
  double m0Min(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->m0Min() : 0. ; }
  double mMax(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->mMax() : 0. ; }
  double m0Max(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->m0Max() : 0. ; }
  double tau0(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->tau0() : 0. ; }
  bool isResonance(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isResonance() : false ; }
  bool mayDecay(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->mayDecay() : false ; }
  bool doExternalDecay(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->doExternalDecay() : false ; }
  bool isVisible(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isVisible() : false ; }
  bool doForceWidth(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->doForceWidth() : false ; }
  bool hasChanged(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->hasChanged() : false ; }
  bool hasChangedMMin(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->hasChangedMMin() : false ; }
  bool hasChangedMMax(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->hasChangedMMax() : false ; }

  // Give back special mass-related quantities.
  bool useBreitWigner(int idIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->useBreitWigner() : false ; }
  double constituentMass(int idIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->constituentMass() : 0. ; }
  double mSel(int idIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->mSel() : 0. ; }
  double mRun(int idIn, double mH) {
    ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->mRun(mH) : 0. ; }

  // Give back other quantities.
  bool canDecay(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->canDecay() : false ; }
  bool isLepton(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isLepton() : false ; }
  bool isQuark(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isQuark() : false ; }
  bool isGluon(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isGluon() : false ; }
  bool isDiquark(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isDiquark() : false ; }
  bool isParton(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isParton() : false ; }
  bool isHadron(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isHadron() : false ; }
  bool isMeson(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isMeson() : false ; }
  bool isBaryon(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isBaryon() : false ; }
  bool isOnium(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isOnium() : false ; }
  bool isOctetHadron(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->isOctetHadron() : false ; }
  int heaviestQuark(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->heaviestQuark(idIn) : 0 ; }
  int baryonNumberType(int idIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->baryonNumberType(idIn) : 0 ; }
  int nQuarksInCode(int idIn, int idQIn) {
    const ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->nQuarksInCode(idQIn) : 0 ; }

  // Change branching ratios.
  void rescaleBR(int idIn, double newSumBR = 1.) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->rescaleBR(newSumBR); }

  // Access methods stored in ResonanceWidths.
  void setResonancePtr(int idIn, ResonanceWidths* resonancePtrIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->setResonancePtr( resonancePtrIn);}
  void resInit(int idIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    if ( ptr ) ptr->resInit(infoPtr, settingsPtr, this, couplingsPtr);}
  double resWidth(int idIn, double mHat, int idInFlav = 0,
    bool openOnly = false, bool setBR = false) {
    ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->resWidth(idIn, mHat,
    idInFlav, openOnly, setBR) : 0.;}
  double resWidthOpen(int idIn, double mHat, int idInFlav = 0) {
    ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->resWidthOpen(idIn, mHat, idInFlav) : 0.;}
  double resWidthStore(int idIn, double mHat, int idInFlav = 0) {
    ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->resWidthStore(idIn, mHat, idInFlav) : 0.;}
  double resOpenFrac(int id1In, int id2In = 0, int id3In = 0);
  double resWidthRescaleFactor(int idIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->resWidthRescaleFactor() : 0.;}
  double resWidthChan(int idIn, double mHat, int idAbs1 = 0,
    int idAbs2 = 0) {
    ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr->resWidthChan( mHat, idAbs1, idAbs2) : 0.;}

  // Return pointer to entry.
  ParticleDataEntry* particleDataEntryPtr(int idIn) {
    ParticleDataEntry* ptr = findParticle(idIn);
    return ( ptr ) ? ptr : &pdt[0]; }

  // Check initialisation status.
  bool getIsInit() {return isInit;}

private:

  // Common data, accessible for the individual particles.
  bool   setRapidDecayVertex;
  int    modeBreitWigner;
  double maxEnhanceBW, mQRun[7], Lambda5Run, intermediateTau0;

  // The individual particle need access to the full database.
  friend class ParticleDataEntry;

  // Pointer to various information on the generation.
  Info*     infoPtr;

  // Pointer to the settings database.
  Settings* settingsPtr;

  // Pointer to the random number generator.
  Rndm*     rndmPtr;

  // Pointer to Standard Model couplings.
  Couplings*   couplingsPtr;

  // All particle data stored in a map.
  map<int, ParticleDataEntry> pdt;

  // Pointer to current particle (e.g. when reading decay channels).
  ParticleDataEntry* particlePtr;

  // Flag that initialization has been performed; whether any failures.
  bool   isInit, readingFailedSave;

  // Method for common setting of particle-specific info.
  void   initCommon();

  // Useful functions for string handling.
  bool   boolString(string tag) { string tagLow = toLower(tag);
    return ( tagLow == "true" || tagLow == "1" || tagLow == "on"
    || tagLow == "yes" || tagLow == "ok" ); }

  // Extract XML value following XML attribute.
  string attributeValue(string line, string attribute);
  bool   boolAttributeValue(string line, string attribute);
  int    intAttributeValue(string line, string attribute);
  double doubleAttributeValue(string line, string attribute);

  // Vector of strings containing the readable lines of the XML file.
  vector<string> xmlFileSav;

  // Stored history of readString statements (common and by subrun).
  vector<string> readStringHistory;
  map<int, vector<string> > readStringSubrun;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ParticleData_H
