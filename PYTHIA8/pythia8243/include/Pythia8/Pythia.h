// Pythia.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class for event generation.
// Pythia: provide the main user interface to everything else.

#ifndef Pythia8_Pythia_H
#define Pythia8_Pythia_H

// Version number defined for use in macros and for consistency checks.
#define PYTHIA_VERSION 8.243
#define PYTHIA_VERSION_INTEGER 8243

// Header files for the Pythia class and for what else the user may need.
#include "Pythia8/Analysis.h"
#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/BeamShape.h"
#include "Pythia8/ColourReconnection.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/HadronLevel.h"
#include "Pythia8/History.h"
#include "Pythia8/Info.h"
#include "Pythia8/JunctionSplitting.h"
#include "Pythia8/LesHouches.h"
#include "Pythia8/Merging.h"
#include "Pythia8/MergingHooks.h"
#include "Pythia8/PartonLevel.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonDistributions.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/PartonVertex.h"
#include "Pythia8/ProcessLevel.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/RHadrons.h"
#include "Pythia8/Ropewalk.h"
#include "Pythia8/Settings.h"
#include "Pythia8/SigmaTotal.h"
#include "Pythia8/SimpleSpaceShower.h"
#include "Pythia8/SimpleTimeShower.h"
#include "Pythia8/SpaceShower.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SLHAinterface.h"
#include "Pythia8/TimeShower.h"
#include "Pythia8/UserHooks.h"

namespace Pythia8 {

//==========================================================================

// Forward declaration of HeavyIons and HIUserHooks classes.
class HeavyIons;
class HIUserHooks;

// The Pythia class contains the top-level routines to generate an event.

class Pythia {

public:

  // Constructor. (See Pythia.cc file.)
  Pythia(string xmlDir = "../share/Pythia8/xmldoc", bool printBanner = true);

  // Constructor to copy settings and particle database from another Pythia
  // object instead of XML files (to speed up multiple initialisations).
  Pythia(Settings& settingsIn, ParticleData& particleDataIn,
    bool printBanner = true);

  // Constructor taking input from streams instead of files.
  Pythia( istream& settingsStrings, istream& particleDataStrings,
    bool printBanner = true);

  // Destructor. (See Pythia.cc file.)
  ~Pythia();

  // Initialise new Pythia object (called by constructors).
  void initPtrs();

  // Check consistency of version numbers (called by constructors).
  bool checkVersion();

  // Read in one update for a setting or particle data from a single line.
  bool readString(string, bool warn = true);

  // Read in updates for settings or particle data from user-defined file.
  bool readFile(string fileName, bool warn = true,
    int subrun = SUBRUNDEFAULT);
  bool readFile(string fileName, int subrun) {
    return readFile(fileName, true, subrun);}
  bool readFile(istream& is = cin, bool warn = true,
    int subrun = SUBRUNDEFAULT);
  bool readFile(istream& is, int subrun) {
    return readFile(is, true, subrun);}

  // Possibility to pass in pointers to PDF's.
  bool setPDFPtr( PDF* pdfAPtrIn, PDF* pdfBPtrIn, PDF* pdfHardAPtrIn = 0,
    PDF* pdfHardBPtrIn = 0, PDF* pdfPomAPtrIn = 0, PDF* pdfPomBPtrIn = 0,
    PDF* pdfGamAPtrIn = 0, PDF* pdfGamBPtrIn = 0, PDF* pdfHardGamAPtrIn = 0,
    PDF* pdfHardGamBPtrIn = 0, PDF* pdfUnresAPtrIn = 0,
    PDF* pdfUnresBPtrIn = 0, PDF* pdfUnresGamAPtrIn = 0,
    PDF* pdfUnresGamBPtrIn = 0, PDF* pdfVMDAPtrIn = 0, PDF* pdfVMDBPtrIn = 0);
  bool setPDFAPtr( PDF* pdfAPtrIn );
  bool setPDFBPtr( PDF* pdfBPtrIn );

  // Set photon fluxes externally. Used with option "PDF:lepton2gammaSet = 2".
  bool setPhotonFluxPtr( PDF* photonFluxAIn, PDF* photonFluxBIn) {
    if ( photonFluxAIn != 0 ) pdfGamFluxAPtr = photonFluxAIn;
    if ( photonFluxBIn != 0 ) pdfGamFluxBPtr = photonFluxBIn;
    return true;}

  // Possibility to pass in pointer to external LHA-interfaced generator.
  bool setLHAupPtr( LHAup* lhaUpPtrIn) {lhaUpPtr = lhaUpPtrIn; return true;}

  // Possibility to pass in pointer for external handling of some decays.
  bool setDecayPtr( DecayHandler* decayHandlePtrIn,
    vector<int> handledParticlesIn) {decayHandlePtr = decayHandlePtrIn;
    handledParticles.resize(0);
    for (int i = 0; i < int(handledParticlesIn.size()); ++i)
      handledParticles.push_back( handledParticlesIn[i] );
    return true;}

  // Possibility to pass in pointer for external random number generation.
  bool setRndmEnginePtr( RndmEngine* rndmEnginePtrIn)
    { return rndm.rndmEnginePtr( rndmEnginePtrIn);}

  // Possibility to pass in pointer for user hooks.
  bool setUserHooksPtr( UserHooks* userHooksPtrIn) {
    if (hasUserHooksVector) delete userHooksPtr;
    hasUserHooksVector = false;
    userHooksPtr = userHooksPtrIn; return true;}

  // Possibility to add further pointers to allow multiple user hooks.
  bool addUserHooksPtr( UserHooks* userHooksPtrIn) {
    if ( !userHooksPtr ) return setUserHooksPtr(userHooksPtrIn);
    UserHooksVector* uhv = dynamic_cast<UserHooksVector*>(userHooksPtr);
    if ( !uhv ) { uhv = new UserHooksVector();
      uhv->hooks.push_back(userHooksPtr); userHooksPtr = uhv; }
    uhv->hooks.push_back(userHooksPtrIn);
    hasUserHooksVector = true; return true;}

  // Possibility to pass in pointer for full merging class.
  bool setMergingPtr( Merging* mergingPtrIn)
    { mergingPtr = mergingPtrIn; return true;}

  // Possibility to pass in pointer for merging hooks.
  bool setMergingHooksPtr( MergingHooks* mergingHooksPtrIn)
    { mergingHooksPtr = mergingHooksPtrIn; return true;}

  // Possibility to pass in pointer for beam shape.
  bool setBeamShapePtr( BeamShape* beamShapePtrIn)
    { beamShapePtr = beamShapePtrIn; return true;}

  // Possibility to pass in pointer(s) for external cross section,
  // with option to include external phase-space generator(s).
  bool setSigmaPtr( SigmaProcess* sigmaPtrIn, PhaseSpace* phaseSpacePtrIn = 0)
    { sigmaPtrs.push_back( sigmaPtrIn);
      phaseSpacePtrs.push_back(phaseSpacePtrIn); return true;}

  // Possibility to pass in pointer(s) for external resonance.
  bool setResonancePtr( ResonanceWidths* resonancePtrIn)
    { resonancePtrs.push_back( resonancePtrIn); return true;}

  // Possibility to pass in pointer for external showers.
  bool setShowerPtr( TimeShower* timesDecPtrIn,
    TimeShower* timesPtrIn = 0, SpaceShower* spacePtrIn = 0)
    { timesDecPtr = timesDecPtrIn; timesPtr = timesPtrIn;
    spacePtr = spacePtrIn; return true;}

  // Possibility to pass in pointer for modelling of heavy ion collisions.
  bool setHeavyIonsPtr( HeavyIons* heavyIonsPtrIn)
    { heavyIonsPtr = heavyIonsPtrIn; return true;}

  // Possibility to pass a HIUserHooks pointer for modifying the
  // behavior of the heavy ion modelling.
  bool setHIHooks(HIUserHooks* hiHooksPtrIn)
    { hiHooksPtr = hiHooksPtrIn; return true; }

  // Possibility to get the pointer to a object modelling heavy ion
  // collisions.
  HeavyIons* getHeavyIonsPtr() { return heavyIonsPtr;}

  // Possibility to pass in pointer for setting of parton space-time vertices.
  bool setPartonVertexPtr( PartonVertex* partonVertexPtrIn)
    { partonVertexPtr = partonVertexPtrIn; return true;}

  // Initialize.
  bool init();

  // Generate the next event.
  bool next();

  // Generate the next event, either with new energies or new beam momenta.
  bool next(double eCMin);
  bool next(double eAin, double eBin);
  bool next(double pxAin, double pyAin, double pzAin,
            double pxBin, double pyBin, double pzBin);

  // Generate only a single timelike shower as in a decay.
  int forceTimeShower( int iBeg, int iEnd, double pTmax, int nBranchMax = 0)
    {  partonSystems.clear(); info.setScalup( 0, pTmax);
    return timesDecPtr->shower( iBeg, iEnd, event, pTmax, nBranchMax); }

  // Generate only the hadronization/decay stage.
  bool forceHadronLevel( bool findJunctions = true);

  // Special routine to allow more decays if on/off switches changed.
  bool moreDecays() {return hadronLevel.moreDecays(event);}

  // Special routine to force R-hadron decay when not done before.
  bool forceRHadronDecays() {return doRHadronDecays();}

  // List the current Les Houches event.
  void LHAeventList() { if (lhaUpPtr != 0) lhaUpPtr->listEvent();}

  // Skip a number of Les Houches events at input.
  bool LHAeventSkip(int nSkip) {
    if (lhaUpPtr != 0) return lhaUpPtr->skipEvent(nSkip);
    return false;}

  // Main routine to provide final statistics on generation.
  void stat();

  // Read in settings values: shorthand, not new functionality.
  bool   flag(string key) {return settings.flag(key);}
  int    mode(string key) {return settings.mode(key);}
  double parm(string key) {return settings.parm(key);}
  string word(string key) {return settings.word(key);}

  // Auxiliary to set parton densities among list of possibilities.
  PDF* getPDFPtr(int idIn, int sequence = 1, string beam = "",
    bool resolved = true);

  // The event record for the parton-level central process.
  Event          process;

  // The event record for the complete event history.
  Event          event;

  // Information on the generation: current subprocess and error statistics.
  Info           info;

  // Settings: databases of flags/modes/parms/words to control run.
  Settings       settings;

  // ParticleData: the particle data table/database.
  ParticleData   particleData;

  // Random number generator.
  Rndm           rndm;

  // Standard Model couplings, including alphaS and alphaEM.
  Couplings      couplings;
  Couplings*     couplingsPtr;

  // SLHA Interface
  SLHAinterface  slhaInterface;

  // The partonic content of each subcollision system (auxiliary to event).
  PartonSystems  partonSystems;

  // Merging object as wrapper for matrix element merging routines.
  Merging*       mergingPtr;

  // Pointer to MergingHooks object for user interaction with the merging.
  // MergingHooks also more generally steers the matrix element merging.
  MergingHooks*  mergingHooksPtr;

  // Pointer to a HeavyIons object for generating heavy ion collisions
  HeavyIons* heavyIonsPtr;

  // Pointer to a HIUserHooks object to modify heavy ion modelling.
  HIUserHooks* hiHooksPtr;

  // The two incoming beams.
  BeamParticle beamA;
  BeamParticle beamB;

private:

  // Copy and = constructors are made private so they cannot be used.
  Pythia(const Pythia&);
  Pythia& operator=(const Pythia&);

  // Constants: could only be changed in the code itself.
  static const double VERSIONNUMBERHEAD, VERSIONNUMBERCODE;
  static const int    NTRY, SUBRUNDEFAULT;

  // Initialization data, extracted from database.
  string xmlPath;
  bool   doProcessLevel, doPartonLevel, doHadronLevel, doSoftQCDall,
         doSoftQCDinel, doCentralDiff, doDiffraction,
         doSoftQCD, doVMDsideA, doVMDsideB, doHardDiff, doResDec,
         doFSRinRes, decayRHadrons, abortIfVeto, checkEvent, checkHistory;
  int    nErrList;
  double epTolErr, epTolWarn, mTolErr, mTolWarn;

  // Initialization data related to photon-photon interactions.
  bool   beamHasGamma, beamAisResGamma, beamBisResGamma, beamAhasResGamma,
         beamBhasResGamma;
  int    gammaMode;

  // Initialization data, extracted from init(...) call.
  bool   isConstructed, isInit, isUnresolvedA, isUnresolvedB, showSaV,
         showMaD, doReconnect, forceHadronLevelCR;
  int    idA, idB, frameType, boostType, nCount, nShowLHA, nShowInfo,
         nShowProc, nShowEvt, reconnectMode;
  double mA, mB, pxA, pxB, pyA, pyB, pzA, pzB, eA, eB,
         pzAcm, pzBcm, eCM, betaZ, gammaZ;
  Vec4   pAinit, pBinit, pAnow, pBnow;
  RotBstMatrix MfromCM, MtoCM;

  // information for error checkout.
  int    nErrEvent;
  vector<int> iErrId, iErrCol, iErrEpm, iErrNan, iErrNanVtx;

  // Pointers to the parton distributions of the two incoming beams.
  PDF* pdfAPtr;
  PDF* pdfBPtr;

  // Extra PDF pointers to be used in hard processes only.
  PDF* pdfHardAPtr;
  PDF* pdfHardBPtr;

  // Extra Pomeron PDF pointers to be used in diffractive processes only.
  PDF* pdfPomAPtr;
  PDF* pdfPomBPtr;

  // Extra Photon PDF pointers to be used in lepton -> gamma processes.
  PDF* pdfGamAPtr;
  PDF* pdfGamBPtr;

  // Extra PDF pointers to be used in hard lepton -> gamma processes.
  PDF* pdfHardGamAPtr;
  PDF* pdfHardGamBPtr;

  // Alternative unresolved PDFs when mixing resolved and unresolved processes.
  PDF* pdfUnresAPtr;
  PDF* pdfUnresBPtr;
  PDF* pdfUnresGamAPtr;
  PDF* pdfUnresGamBPtr;

  // PDF pointers to externally provided photon fluxes.
  PDF* pdfGamFluxAPtr;
  PDF* pdfGamFluxBPtr;

  // Extra VMD PDF pointers to be used in SoftQCD with gammas.
  PDF* pdfVMDAPtr;
  PDF* pdfVMDBPtr;

  // Keep track when "new" has been used and needs a "delete" for PDF's etc.
  bool useNewPdfA, useNewPdfB, useNewPdfHard, useNewPdfPomA, useNewPdfPomB,
    useNewPdfGamA, useNewPdfGamB, useNewPdfHardGamA, useNewPdfHardGamB,
    useNewPdfUnresA, useNewPdfUnresB, useNewPdfUnresGamA, useNewPdfUnresGamB,
    useNewPdfVMDA, useNewPdfVMDB, hasUserHooksVector;

  // Alternative Pomeron beam-inside-beam.
  BeamParticle beamPomA;
  BeamParticle beamPomB;

  // Alternative photon beam-inside-beam.
  BeamParticle beamGamA;
  BeamParticle beamGamB;

  // Alternative VMD beam-inside-beam.
  BeamParticle beamVMDA;
  BeamParticle beamVMDB;

  // LHAup object for generating external events.
  bool   doLHA, useNewLHA;
  LHAup* lhaUpPtr;

  // Pointer to external decay handler and list of particles it handles.
  DecayHandler* decayHandlePtr;
  vector<int>   handledParticles;

  // Pointer to UserHooks object for user interaction with program.
  UserHooks* userHooksPtr;
  bool       hasUserHooks, doVetoProcess, doVetoPartons, retryPartonLevel;

  // Pointer to BeamShape object for beam momentum and interaction vertex.
  BeamShape* beamShapePtr;
  bool       useNewBeamShape, doMomentumSpread, doVertexSpread, doVarEcm;
  double     eMinPert, eWidthPert;

  // Pointers to external processes derived from the Pythia base classes.
  vector<SigmaProcess*> sigmaPtrs;

  // Pointers to external phase-space generators derived from Pythia
  // base classes.
  vector<PhaseSpace*> phaseSpacePtrs;

  // Pointers to external calculation of resonance widths.
  vector<ResonanceWidths*> resonancePtrs;

  // Pointers to timelike and spacelike showers.
  TimeShower*  timesDecPtr;
  TimeShower*  timesPtr;
  SpaceShower* spacePtr;
  bool         useNewTimesDec, useNewTimes, useNewSpace;

  // Pointer to assign space-time vertices during parton evolution.
  PartonVertex* partonVertexPtr;
  bool          useNewPartonVertex;

  // The main generator class to define the core process of the event.
  ProcessLevel processLevel;

  // The main generator class to produce the parton level of the event.
  PartonLevel partonLevel;

  // The main generator class to perform trial showers of the event.
  PartonLevel trialPartonLevel;

  // Flags for defining the merging scheme.
  bool        hasMerging, hasOwnMerging;
  bool        hasMergingHooks, hasOwnMergingHooks, doMerging;

  // The Colour reconnection class.
  ColourReconnection colourReconnection;

  // The junction spltiting class.
  JunctionSplitting junctionSplitting;

  // The main generator class to produce the hadron level of the event.
  HadronLevel hadronLevel;

  // The total cross section class is used both on process and parton level.
  SigmaTotal sigmaTot;

  // The RHadrons class is used both at PartonLevel and HadronLevel.
  RHadrons   rHadrons;

  // Flags for handling generation of heavy ion collisons.
  bool        hasHeavyIons, hasOwnHeavyIons, doHeavyIons;

  // Write the Pythia banner, with symbol and version information.
  void banner();

  // Check for lines in file that mark the beginning of new subrun.
  int readSubrun(string line, bool warn = true);

  // Check for lines that mark the beginning or end of commented section.
  int readCommented(string line);

  // Check that combinations of settings are allowed; change if not.
  void checkSettings();

  // Check that beams and beam combination can be handled.
  bool checkBeams();

  // Calculate kinematics at initialization.
  bool initKinematics();

  // Set up pointers to PDFs.
  bool initPDFs();

  // Recalculate kinematics for each event when beam momentum has a spread.
  void nextKinematics();

  // Boost from CM frame to lab frame, or inverse. Set production vertex.
  void boostAndVertex(bool toLab, bool setVertex);

  // Perform R-hadron decays.
  bool doRHadronDecays();

  // Check that the final event makes sense.
  bool check();

  // Initialization of SLHA data.
  bool initSLHA ();
  stringstream particleDataBuffer;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Pythia_H
