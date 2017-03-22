// Pythia.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class for event generation.
// Pythia: provide the main user interface to everything else.

#ifndef Pythia8_Pythia_H
#define Pythia8_Pythia_H

#include "Analysis.h"
#include "Basics.h"
#include "BeamParticle.h"
#include "BeamShape.h"
#include "Event.h"
#include "FragmentationFlavZpT.h"
#include "HadronLevel.h"
#include "History.h"
#include "Info.h"
#include "LesHouches.h"
#include "PartonLevel.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "PartonSystems.h"
#include "ProcessLevel.h"
#include "PythiaStdlib.h"
#include "ResonanceWidths.h"
#include "RHadrons.h"
#include "Settings.h"
#include "SigmaTotal.h"
#include "SpaceShower.h"
#include "StandardModel.h"
#include "SusyLesHouches.h"
#include "TimeShower.h"
#include "UserHooks.h"
#include "MergingHooks.h"
#include "Merging.h"

namespace Pythia8 {
 
//==========================================================================

// The Pythia class contains the top-level routines to generate an event.

class Pythia {

public:

  // Constructor. (See Pythia.cc file.)
  Pythia(string xmlDir = "../xmldoc", bool printBanner = true);

  // Destructor. (See Pythia.cc file.)
  ~Pythia();

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
    PDF* pdfHardBPtrIn = 0, PDF* pdfPomAPtrIn = 0, PDF* pdfPomBPtrIn = 0);

  // Possibility to pass in pointer to external LHA-interfaced generator.
  bool setLHAupPtr( LHAup* lhaUpPtrIn) {lhaUpPtr = lhaUpPtrIn; return true;}

  // Possibility to pass in pointer for external handling of some decays.
  bool setDecayPtr( DecayHandler* decayHandlePtrIn, 
    vector<int> handledParticlesIn) {decayHandlePtr = decayHandlePtrIn; 
    handledParticles.resize(0); 
    for(int i = 0; i < int(handledParticlesIn.size()); ++i)
    handledParticles.push_back( handledParticlesIn[i] ); return true;}  

  // Possibility to pass in pointer for external random number generation.
  bool setRndmEnginePtr( RndmEngine* rndmEnginePtrIn) 
    { return rndm.rndmEnginePtr( rndmEnginePtrIn);}  

  // Possibility to pass in pointer for user hooks. 
  bool setUserHooksPtr( UserHooks* userHooksPtrIn) 
    { userHooksPtr = userHooksPtrIn; return true;} 

  // Possibility to pass in pointer for merging hooks. 
  bool setMergingHooksPtr( MergingHooks* mergingHooksPtrIn) 
    { mergingHooksPtr = mergingHooksPtrIn; return true;} 

  // Possibility to pass in pointer for beam shape. 
  bool setBeamShapePtr( BeamShape* beamShapePtrIn) 
    { beamShapePtr = beamShapePtrIn; return true;} 

  // Possibility to pass in pointer(s) for external cross section.
  bool setSigmaPtr( SigmaProcess* sigmaPtrIn) 
    { sigmaPtrs.push_back( sigmaPtrIn); return true;} 

  // Possibility to pass in pointer(s) for external resonance.
  bool setResonancePtr( ResonanceWidths* resonancePtrIn) 
    { resonancePtrs.push_back( resonancePtrIn); return true;} 

  // Possibility to pass in pointer for external showers.
  bool setShowerPtr( TimeShower* timesDecPtrIn, 
    TimeShower* timesPtrIn = 0, SpaceShower* spacePtrIn = 0) 
    { timesDecPtr = timesDecPtrIn; timesPtr = timesPtrIn;
    spacePtr = spacePtrIn; return true;} 

  // Initialization using the Beams variables.
  bool init();

  // Deprecated: initialization in the CM frame.
  bool init( int idAin, int idBin, double eCMin);

  // Deprecated: initialization with two collinear beams, e.g. fixed target.
  bool init( int idAin, int idBin, double eAin, double eBin);

  // Deprecated: initialization with two acollinear beams.
  bool init( int idAin, int idBin, double pxAin, double pyAin, 
    double pzAin, double pxBin, double pyBin, double pzBin);

  // Deprecated: initialization by a Les Houches Event File.
  bool init( string LesHouchesEventFile, bool skipInit = false);

  // Deprecated: initialization according to the Les Houches Accord.
  bool init( LHAup* lhaUpPtrIn);

  // Generate the next event.
  bool next(); 

  // Generate only a single timelike shower as in a decay.
  int forceTimeShower( int iBeg, int iEnd, double pTmax, int nBranchMax = 0)
    { return timesDecPtr->shower( iBeg, iEnd, event, pTmax, nBranchMax); }

  // Generate only the hadronization/decay stage.
  bool forceHadronLevel( bool findJunctions = true);

  // Special routine to allow more decays if on/off switches changed.
  bool moreDecays() {return hadronLevel.moreDecays(event);}

  // Special routine to force R-hadron decay when not done before.
  bool forceRHadronDecays() {return doRHadronDecays();}

  // List the current Les Houches event.
  void LHAeventList(ostream& os = cout) {
    if (lhaUpPtr != 0) lhaUpPtr->listEvent(os);}

  // Skip a number of Les Houches events at input.
  bool LHAeventSkip(int nSkip) {
    if (lhaUpPtr != 0) return lhaUpPtr->skipEvent(nSkip); return false;}

  // Main routine to provide final statistics on generation.
  void stat();

  // Deprecated: alternative to provide final statistics on generation.
  void statistics(bool all = false, bool reset = false);

  // Read in settings values: shorthand, not new functionality.
  bool   flag(string key) {return settings.flag(key);}
  int    mode(string key) {return settings.mode(key);} 
  double parm(string key) {return settings.parm(key);}
  string word(string key) {return settings.word(key);}

  // Auxiliary to set parton densities among list of possibilities.
  PDF* getPDFPtr(int idIn, int sequence = 1);

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
  CoupSUSY       coupSUSY;
  Couplings*     couplingsPtr;

  // SusyLesHouches - SLHA object for interface to SUSY spectra.
  SusyLesHouches slha;

  // The partonic content of each subcollision system (auxiliary to event).
  PartonSystems  partonSystems; 

  // Merging object as wrapper for matrix element merging routines.
  Merging        merging;

  // Pointer to MergingHooks object for user interaction with the merging.
  // MergingHooks also more generally steers the matrix element merging.
  MergingHooks*  mergingHooksPtr;

private: 

  // Copy and = constructors are made private so they cannot be used.
  Pythia(const Pythia&);
  Pythia& operator=(const Pythia&);

  // Constants: could only be changed in the code itself.
  static const double VERSIONNUMBERCODE;
  static const int    NTRY, SUBRUNDEFAULT;

  // Initialization data, extracted from database.
  string xmlPath;
  bool   doProcessLevel, doPartonLevel, doHadronLevel, doDiffraction, 
         doResDec, doFSRinRes, decayRHadrons, abortIfVeto, checkEvent, 
         checkHistory;
  int    nErrList;
  double epTolErr, epTolWarn;

  // Initialization data, extracted from init(...) call.
  bool   isConstructed, isInit, isUnresolvedA, isUnresolvedB, showSaV, 
         showMaD;
  int    idA, idB, frameType, boostType, nCount, nShowLHA, nShowInfo,
         nShowProc, nShowEvt;  
  double mA, mB, pxA, pxB, pyA, pyB, pzA, pzB, eA, eB, 
         pzAcm, pzBcm, eCM, betaZ, gammaZ;
  Vec4   pAinit, pBinit, pAnow, pBnow;
  RotBstMatrix MfromCM, MtoCM;

  // information for error checkout.
  int    nErrEvent;
  vector<int> iErrId, iErrCol, iErrNan, iErrNanVtx;

  // Pointers to the parton distributions of the two incoming beams.
  PDF* pdfAPtr;  
  PDF* pdfBPtr; 

  // Extra PDF pointers to be used in hard processes only. 
  PDF* pdfHardAPtr;  
  PDF* pdfHardBPtr; 

  // Extra Pomeron PDF pointers to be used in diffractive processes only. 
  PDF* pdfPomAPtr;  
  PDF* pdfPomBPtr; 

  // Keep track when "new" has been used and needs a "delete" for PDF's.  
  bool useNewPdfA, useNewPdfB, useNewPdfHard, useNewPdfPomA, useNewPdfPomB;

  // The two incoming beams.
  BeamParticle beamA;
  BeamParticle beamB;

  // Alternative Pomeron beam-inside-beam.
  BeamParticle beamPomA;
  BeamParticle beamPomB;

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
  bool       useNewBeamShape, doMomentumSpread, doVertexSpread;

  // Pointers to external processes derived from the Pythia base classes.
  vector<SigmaProcess*> sigmaPtrs;  

  // Pointers to external calculation of resonance widths.
  vector<ResonanceWidths*> resonancePtrs;

  // Pointers to timelike and spacelike showers.
  TimeShower*  timesDecPtr;
  TimeShower*  timesPtr;
  SpaceShower* spacePtr;
  bool         useNewTimes, useNewSpace;

  // The main generator class to define the core process of the event.
  ProcessLevel processLevel;

  // The main generator class to produce the parton level of the event.
  PartonLevel partonLevel;

  // The main generator class to perform trial showers of the event.
  PartonLevel trialPartonLevel;

  // Flags for defining the merging scheme.
  bool        hasMergingHooks, hasOwnMergingHooks, doMerging;

  // The main generator class to produce the hadron level of the event.
  HadronLevel hadronLevel;

  // The total cross section class is used both on process and parton level.
  SigmaTotal sigmaTot; 

  // The RHadrons class is used both at PartonLevel and HadronLevel.
  RHadrons   rHadrons;

  // Write the Pythia banner, with symbol and version information.
  void banner(ostream& os = cout);

  // Check for lines in file that mark the beginning of new subrun.
  int readSubrun(string line, bool warn = true, ostream& os = cout);

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
  bool check(ostream& os = cout);

  // Initialization of SLHA data.
  bool initSLHA ();

};
 
//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Pythia_H
