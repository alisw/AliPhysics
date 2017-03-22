// ParticleData.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the classes containing particle data.
// DecayChannel contains info on a single decay channel.
// ParticleDataEntry contains info on a single particle species.
// ParticleData collects info on all particles as a map.

#ifndef Pythia8_ParticleData_H
#define Pythia8_ParticleData_H

#include "Basics.h"
#include "Info.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "StandardModel.h"

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
    meModeSave(meModeIn), nProd(0), hasChangedSave(true) {
    prod[0] = prod0; prod[1] = prod1; prod[2] = prod2; prod[3] = prod3; 
    prod[4] = prod4; prod[5] = prod5; prod[6] = prod6; prod[7] = prod7; 
    for (int j = 0; j < 8; ++j) if (prod[j] != 0 && j == nProd) ++nProd; }  

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
    tau0Save(tau0In), hasAntiSave(false), hasChangedSave(true), 
    resonancePtr(0) {setDefaults();} 
  ParticleDataEntry(int idIn, string nameIn, string antiNameIn, 
    int spinTypeIn = 0, int chargeTypeIn = 0, int colTypeIn = 0, 
    double m0In = 0., double mWidthIn = 0., double mMinIn = 0., 
    double mMaxIn = 0., double tau0In = 0.) : idSave(abs(idIn)), 
    nameSave(nameIn), antiNameSave(antiNameIn), spinTypeSave(spinTypeIn), 
    chargeTypeSave(chargeTypeIn), colTypeSave(colTypeIn), m0Save(m0In), 
    mWidthSave (mWidthIn), mMinSave(mMinIn), mMaxSave(mMaxIn), 
    tau0Save(tau0In), hasAntiSave(true), hasChangedSave(true),
    resonancePtr(0) {setDefaults(); 
    if (toLower(antiNameIn) == "void") hasAntiSave = false;} 

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
    mMinSave = mMinIn; mMaxSave = mMaxIn; tau0Save = tau0In; 
    setDefaults(); hasChangedSave = true;}

  // Change current values one at a time (or set if not set before).
  // (Must use set here since else name+signature clash with get methods.)
  void setName(string nameIn) {nameSave = nameIn; hasChangedSave = true;}
  void setAntiName(string antiNameIn) {antiNameSave = antiNameIn; 
    hasChangedSave = true;}
  void setNames(string nameIn, string antiNameIn) {nameSave = nameIn; 
    antiNameSave = antiNameIn; hasAntiSave = true; if (toLower(antiNameIn) 
    == "void") hasAntiSave = false; hasChangedSave = true;}
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
  void setMMin(double mMinIn) {mMinSave = mMinIn; hasChangedSave = true;}
  void setMMax(double mMaxIn) {mMaxSave = mMaxIn; hasChangedSave = true;}
  // Special options specifically when cutting wings of Breit-Wigners.
  void setMMinNoChange(double mMinIn) {mMinSave = mMinIn;}
  void setMMaxNoChange(double mMaxIn) {mMaxSave = mMaxIn;}
  void setTau0(double tau0In) {tau0Save = tau0In; hasChangedSave = true;}
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
      channels[i].setHasChanged(hasChangedIn);}

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
         if (channels[i].hasChanged()) return true; return false;}

  // Set and give back several mass-related quantities.
  void   initBWmass(); 
  double constituentMass()        const { return constituentMassSave; } 
  double mass(); 
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

  // Intermediate octet ccbar or bbar states in colour-octet model. 
  bool   isOctetHadron()  const {return (idSave == 9900441
         || idSave == 9900443 || idSave == 9900551 || idSave == 9900553 
         || idSave == 9910441 || idSave == 9910551); }
  int    heaviestQuark(int idIn = 1)    const; 
  int    baryonNumberType(int idIn = 1) const;

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
  static const int    INVISIBLENUMBER, INVISIBLETABLE[50], KNOWNNOWIDTH[3];
  static const double MAXTAU0FORDECAY,MINMASSRESONANCE, NARROWMASS,
                      CONSTITUENTMASSTABLE[10];

  // Particle data.
  int    idSave;
  string nameSave, antiNameSave;
  int    spinTypeSave, chargeTypeSave, colTypeSave;
  double m0Save, mWidthSave, mMinSave, mMaxSave, tau0Save, 
         constituentMassSave;
  bool   hasAntiSave, isResonanceSave, mayDecaySave, doExternalDecaySave, 
         isVisibleSave, doForceWidthSave, hasChangedSave;

  // Extra data for mass selection according to a Breit-Wigner.
  int    modeBWnow;
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

  // Useful functions for string handling.
  string toLower(const string& nameConv) { string temp(nameConv);
    for (int i = 0; i < int(temp.length()); ++i) temp[i] = tolower(temp[i]); 
    return temp; }

};

//==========================================================================

// This class holds a map of all ParticleDataEntries.

class ParticleData {

public:

  // Constructor.
  ParticleData() : infoPtr(0), settingsPtr(0), rndmPtr(0), couplingsPtr(0), 
    particlePtr(0), isInit(false), readingFailedSave(false) {}

  // Initialize pointers.
  void initPtr(Info* infoPtrIn, Settings* settingsPtrIn, Rndm* rndmPtrIn, 
    Couplings* couplingsPtrIn) {infoPtr = infoPtrIn; 
    settingsPtr = settingsPtrIn; rndmPtr = rndmPtrIn; 
    couplingsPtr = couplingsPtrIn;}
 
  // Read in database from specific file.
  bool init(string startFile = "../xmldoc/ParticleData.xml") {
    initCommon(); return readXML(startFile);}

  // Overwrite existing database by reading from specific file.
  bool reInit(string startFile, bool xmlFormat = true) { initCommon();
    return (xmlFormat) ? readXML(startFile) : readFF(startFile);}

  // Initialize pointers, normal Breit-Wigners and special resonances.
  void initWidths(vector<ResonanceWidths*> resonancePtrs);

  // Read or list whole (or part of) database from/to an XML file.
  bool readXML(string inFile, bool reset = true) ; 
  void listXML(string outFile); 

  // Read or list whole (or part of) database from/to a free format file.
  bool readFF(string inFile, bool reset = true) ; 
  void listFF(string outFile); 

  // Read in one update from a single line.
  bool readString(string lineIn, bool warn = true, ostream& os = cout) ; 

  // Keep track whether any readings have failed, invalidating run setup.
  bool readingFailed() {return readingFailedSave;} 

  // Print out table of whole database, or of only part of it.
  void listAll(ostream& os = cout) {list(false, true, os);} 
  void listChanged(ostream& os = cout) {list(true, false, os);} 
  void listChanged(bool changedRes, ostream& os = cout) {
    list(true, changedRes, os);} 
  void list(bool changedOnly = false, bool changedRes = true, 
    ostream& os = cout); 

  // Print out specified particles.
  void list(int idList, ostream& os = cout) {vector<int> idListTemp; 
    idListTemp.push_back(idList); list( idListTemp, os);} 
  void list(vector<int> idList, ostream& os = cout); 

  // Check that table makes sense, especially for decays.
  void checkTable(ostream& os = cout) {checkTable(1, os);};
  void checkTable(int verbosity, ostream& os = cout) ;
 
  // Add new entry.
  void addParticle(int idIn, string nameIn = " ", int spinTypeIn = 0, 
    int chargeTypeIn = 0, int colTypeIn = 0, double m0In = 0., 
    double mWidthIn = 0., double mMinIn = 0., double mMaxIn = 0., 
    double tau0In = 0.) { pdt[abs(idIn)] = ParticleDataEntry(idIn, 
    nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, 
    mMinIn, mMaxIn, tau0In); }  
  void addParticle(int idIn, string nameIn, string antiNameIn, 
    int spinTypeIn = 0, int chargeTypeIn = 0, int colTypeIn = 0, 
    double m0In = 0., double mWidthIn = 0., double mMinIn = 0., 
    double mMaxIn = 0., double tau0In = 0.) { pdt[abs(idIn)] 
    = ParticleDataEntry(idIn, nameIn, antiNameIn, spinTypeIn, 
    chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn, tau0In); }  

  // Reset all the properties of an entry in one go.
  void setAll(int idIn, string nameIn, string antiNameIn, 
    int spinTypeIn = 0, int chargeTypeIn = 0, int colTypeIn = 0, 
    double m0In = 0., double mWidthIn = 0., double mMinIn = 0., 
    double mMaxIn = 0.,double tau0In = 0.) { if (isParticle(idIn)) 
    pdt[abs(idIn)].setAll( nameIn, antiNameIn, spinTypeIn, chargeTypeIn, 
    colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn, tau0In); }  

  // Query existence of an entry.
  bool isParticle(int idIn) {
    if (pdt.find(abs(idIn)) == pdt.end()) return false;
    if (idIn > 0 || pdt[abs(idIn)].hasAnti()) return true;
    return false; }

  // Return the id of the sequentially next particle stored in table.
  int nextId(int idIn) ;

  // Change current values one at a time (or set if not set before).
  void name(int idIn, string nameIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setName(nameIn); }
  void antiName(int idIn, string antiNameIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setAntiName(antiNameIn); }
  void names(int idIn, string nameIn, string antiNameIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setNames(nameIn, antiNameIn); }
  void spinType(int idIn, int spinTypeIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setSpinType(spinTypeIn); }
  void chargeType(int idIn, int chargeTypeIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setChargeType(chargeTypeIn); }
  void colType(int idIn, int colTypeIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setColType(colTypeIn); }
  void m0(int idIn, double m0In) {
    if (isParticle(idIn)) pdt[abs(idIn)].setM0(m0In); }
  void mWidth(int idIn, double mWidthIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setMWidth(mWidthIn); }
  void mMin(int idIn, double mMinIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setMMin(mMinIn); }
  void mMax(int idIn, double mMaxIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setMMax(mMaxIn); }
  void tau0(int idIn, double tau0In) {
    if (isParticle(idIn)) pdt[abs(idIn)].setTau0(tau0In); }
  void isResonance(int idIn, bool isResonanceIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setIsResonance(isResonanceIn); }
  void mayDecay(int idIn, bool mayDecayIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setMayDecay(mayDecayIn); }
  void doExternalDecay(int idIn, bool doExternalDecayIn) {
    if (isParticle(idIn)) 
    pdt[abs(idIn)].setDoExternalDecay(doExternalDecayIn); }
  void isVisible(int idIn, bool isVisibleIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setIsVisible(isVisibleIn); }
  void doForceWidth(int idIn, bool doForceWidthIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setDoForceWidth(doForceWidthIn); }
  void hasChanged(int idIn, bool hasChangedIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setHasChanged(hasChangedIn); }
 
  // Give back current values. 
  bool hasAnti(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].hasAnti() : false ; } 
  string name(int idIn) {
    return (isParticle(abs(idIn))) ? pdt[abs(idIn)].name(idIn) : " "; }
  int spinType(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].spinType() : 0 ; } 
  int chargeType(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].chargeType(idIn) : 0 ; } 
  double charge(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].charge(idIn) : 0 ; } 
  int colType(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].colType(idIn) : 0 ; } 
  double m0(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].m0() : 0. ; } 
  double mWidth(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].mWidth() : 0. ; } 
  double mMin(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].mMin() : 0. ; } 
  double m0Min(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].m0Min() : 0. ; } 
  double mMax(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].mMax() : 0. ; } 
  double m0Max(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].m0Max() : 0. ; } 
  double tau0(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].tau0() : 0. ; } 
  bool isResonance(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isResonance() : false ; } 
  bool mayDecay(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].mayDecay() : false ; } 
  bool doExternalDecay(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].doExternalDecay() : false ; }
  bool isVisible(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isVisible() : false ; } 
  bool doForceWidth(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].doForceWidth() : false ; } 
  bool hasChanged(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].hasChanged() : false ; } 

  // Give back special mass-related quantities.
  bool useBreitWigner(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].useBreitWigner() : false ; } 
  double constituentMass(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].constituentMass() : 0. ; } 
  double mass(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].mass() : 0. ; } 
  double mRun(int idIn, double mH) {
    return isParticle(idIn) ? pdt[abs(idIn)].mRun(mH) : 0. ; } 

  // Give back other quantities.
  bool canDecay(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].canDecay() : false ; }
  bool isLepton(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isLepton() : false ; } 
  bool isQuark(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isQuark() : false ; } 
  bool isGluon(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isGluon() : false ; } 
  bool isDiquark(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isDiquark() : false ; } 
  bool isParton(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isParton() : false ; } 
  bool isHadron(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isHadron() : false ; } 
  bool isMeson(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isMeson() : false ; } 
  bool isBaryon(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isBaryon() : false ; } 
  bool isOctetHadron(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isOctetHadron() : false ; } 
  int heaviestQuark(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].heaviestQuark(idIn) : 0 ; }  
  int baryonNumberType(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].baryonNumberType(idIn) : 0 ; }  

  // Change branching ratios.
  void rescaleBR(int idIn, double newSumBR = 1.) {
    if (isParticle(idIn)) pdt[abs(idIn)].rescaleBR(newSumBR); }

  // Access methods stored in ResonanceWidths.
  void setResonancePtr(int idIn, ResonanceWidths* resonancePtrIn) { 
    if (isParticle(idIn)) pdt[abs(idIn)].setResonancePtr( resonancePtrIn);} 
  void resInit(int idIn) { if (isParticle(idIn)) 
    pdt[abs(idIn)].resInit(infoPtr, settingsPtr, this, couplingsPtr);} 
  double resWidth(int idIn, double mHat, int idInFlav = 0, 
    bool openOnly = false, bool setBR = false) {
    return isParticle(idIn) ? pdt[abs(idIn)].resWidth(idIn, mHat,
    idInFlav, openOnly, setBR) : 0.;}
  double resWidthOpen(int idIn, double mHat, int idInFlav = 0) {
    return isParticle(idIn) ? pdt[abs(idIn)].resWidthOpen(idIn, mHat, 
    idInFlav) : 0.;}
  double resWidthStore(int idIn, double mHat, int idInFlav = 0) {
    return isParticle(idIn) ? pdt[abs(idIn)].resWidthStore(idIn, mHat, 
    idInFlav) : 0.;}
  double resOpenFrac(int id1In, int id2In = 0, int id3In = 0);
  double resWidthRescaleFactor(int idIn) { return isParticle(idIn) 
    ? pdt[abs(idIn)].resWidthRescaleFactor() : 0.;}
  double resWidthChan(int idIn, double mHat, int idAbs1 = 0, 
    int idAbs2 = 0) { return isParticle(idIn) 
    ? pdt[abs(idIn)].resWidthChan( mHat, idAbs1, idAbs2) : 0.;}
  
  // Return pointer to entry.
  ParticleDataEntry* particleDataEntryPtr(int idIn) {
    return (isParticle(idIn)) ? &pdt[abs(idIn)] : &pdt[0]; }

private:

  // Common data, accessible for the individual particles.
  int    modeBreitWigner;
  double maxEnhanceBW, mQRun[7], Lambda5Run;

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
  string toLower(const string& nameConv) { string temp(nameConv);
    for (int i = 0; i < int(temp.length()); ++i) temp[i] = tolower(temp[i]); 
    return temp; }
  bool   boolString(string tag) { string tagLow = toLower(tag);
    return ( tagLow == "true" || tagLow == "1" || tagLow == "on" 
    || tagLow == "yes" || tagLow == "ok" ); }  

  // Extract XML value following XML attribute.
  string attributeValue(string line, string attribute);
  bool   boolAttributeValue(string line, string attribute);
  int    intAttributeValue(string line, string attribute);
  double doubleAttributeValue(string line, string attribute);

};
 
//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ParticleData_H
