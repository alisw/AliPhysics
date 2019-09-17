// TimeShower.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the base class of timelike final-state showers.
// TimeShower: handles the showering description.

#ifndef Pythia8_TimeShower_H
#define Pythia8_TimeShower_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/PartonVertex.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/UserHooks.h"
#include "Pythia8/MergingHooks.h"

namespace Pythia8 {

//==========================================================================

// The TimeShower class does timelike showers.

class TimeShower {

public:

  // Constructor.
  TimeShower() : mergingHooksPtr(), infoPtr(), settingsPtr(),
    particleDataPtr(), rndmPtr(), coupSMPtr(), beamAPtr(), beamBPtr(),
    partonSystemsPtr(), userHooksPtr(), partonVertexPtr(),
    doUncertainties(), uVarMuSoftCorr(), uVarMPIshowers(), noResVariations(),
    noProcVariations(), nUncertaintyVariations(), nVarQCD(), uVarNflavQ(),
    dASmax(), cNSpTmin(), uVarpTmin2(), overFactor(), varG2GGmuRfac(),
    varQ2QGmuRfac(), varG2QQmuRfac(), varX2XGmuRfac(), varG2GGcNS(),
    varQ2QGcNS(), varG2QQcNS(), varX2XGcNS(), varPDFplus(), varPDFminus(),
    varPDFmember() {}

  // Destructor.
  virtual ~TimeShower() {}

  // Initialize various pointers.
  // (Separated from rest of init since not virtual.)
  void initPtr(Info* infoPtrIn, Settings* settingsPtrIn,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    CoupSM* coupSMPtrIn, PartonSystems* partonSystemsPtrIn,
    UserHooks* userHooksPtrIn, MergingHooks* mergingHooksPtrIn,
    PartonVertex* partonVertexPtrIn) { infoPtr = infoPtrIn;
    settingsPtr = settingsPtrIn; particleDataPtr = particleDataPtrIn;
    rndmPtr = rndmPtrIn; coupSMPtr = coupSMPtrIn;
    partonSystemsPtr = partonSystemsPtrIn; userHooksPtr = userHooksPtrIn;
    mergingHooksPtr = mergingHooksPtrIn; partonVertexPtr = partonVertexPtrIn;
  }

  // New beams possible for handling of hard diffraction. (Not virtual.)
  void reassignBeamPtrs( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    int beamOffsetIn = 0) {beamAPtr = beamAPtrIn; beamBPtr = beamBPtrIn;
    beamOffset = beamOffsetIn;}

  // Initialize alphaStrong and related pTmin parameters.
  // Usage: init( beamAPtr, beamBPtr).
  virtual void init( BeamParticle* = 0, BeamParticle* = 0) {}

  // Find whether to limit maximum scale of emissions, and whether to dampen.
  // Usage: limitPTmax( event, Q2Fac, double Q2Ren).
  virtual bool limitPTmax( Event& , double = 0., double = 0.) {return true;}

  // Top-level routine to do a full time-like shower in resonance decay.
  // Usage: shower( iBeg, iEnd, event, pTmax, nBranchMax).
  virtual int shower( int , int , Event& , double , int = 0) {return 0;}

  // Top-level routine for QED radiation in hadronic decay to two leptons.
  // Usage: showerQED( i1, i2, event, pTmax).
  virtual int showerQED( int , int , Event& , double ) {return 0;}

  // Optional method to add QED showers after remnants have been added
  // but before hadronisation. (Called from PartonLevel.)
  virtual int showerQEDafterRemnants(Event&) { return 0; }

  // Global recoil: reset counters and store locations of outgoing partons.
  // Usage: prepareGlobal( event).
  virtual void prepareGlobal( Event& ) {}

  // Prepare system for evolution after each new interaction; identify ME.
  // Usage: prepare( iSys, event, limitPTmax).
  virtual void prepare( int , Event& , bool = true) {}

  // Update dipole list after a multiparton interactions rescattering.
  // Usage: rescatterUpdate( iSys, event).
  virtual void rescatterUpdate( int , Event& ) {}

  // Update dipole list after each ISR emission.
  // Usage: update( iSys, event, hasWeakRad).
  virtual void update( int , Event& , bool = false) {}

  // Select next pT in downwards evolution.
  // Usage: pTnext( event, pTbegAll, pTendAll, isFirstTrial, doTrialIn).
  virtual double pTnext( Event& , double , double , bool = false, bool = false)
    { return 0.;}

  // ME corrections and kinematics that may give failure.
  // Usage: branch( event, isInterleaved).
  virtual bool branch( Event& , bool = false) {return true;}

  // Print dipole list; for debug mainly.
  virtual void list() const {}

  // Initialize data members for calculation of uncertainty bands.
  virtual bool initUncertainties() {return false;}

  // Tell whether FSR has done a weak emission.
  virtual bool getHasWeaklyRadiated() {return false;}

  // Tell which system was the last processed one.
  virtual int system() const {return 0;}

  // Potential enhancement factor of pTmax scale for hardest emission.
  virtual double enhancePTmax() {return 1.;}

  // Provide the pT scale of the last branching in the above shower.
  virtual double pTLastInShower() {return 0.;}

  // Functions to allow usage of shower kinematics, evolution variables,
  // and splitting probabilities outside of shower.
  // Virtual so that shower plugins can overwrite these functions.
  // This makes it possible for another piece of the code to request
  // these - which is very convenient for merging.
  // Function variable names are not included to avoid compiler warnings.
  // Please see the documentation under "Implement New Showers" for details.

  // Return clustering kinematics - as needed for merging.
  virtual Event clustered( const Event& , int , int , int , string )
    { return Event();}

  // Return the evolution variable(s).
  // Important note: this map must contain the following entries
  // - a key "t" for the value of the shower evolution variable;
  // - a key "tRS" for the value of the shower evolution variable
  //   from which the shower would be restarted after a branching;
  // - a key "scaleAS" for the argument of alpha_s used for the branching;
  // - a key "scalePDF" for the argument of the PDFs used for the branching.
  // Usage: getStateVariables( event, iRad, iEmt, iRec,  name)
  virtual map<string, double> getStateVariables (const Event& , int , int ,
    int , string ) { return map<string,double>();}

  // Check if attempted clustering is handled by timelike shower
  // Usage: isTimelike( event, iRad, iEmt, iRec, name)
  virtual bool isTimelike(const Event& , int , int , int , string )
    { return false; }

  // Return a string identifier of a splitting.
  // Usage: getSplittingName( event, iRad, iEmt, iRec)
  virtual vector<string> getSplittingName( const Event& , int, int , int)
    { return vector<string>();}

  // Return the splitting probability.
  // Usage: getSplittingProb( event, iRad, iEmt, iRec)
  virtual double getSplittingProb( const Event& , int , int , int , string )
    { return 0.;}
  virtual bool allowedSplitting( const Event& , int , int)
    { return true; }
  virtual vector<int> getRecoilers( const Event&, int, int, string)
    { return vector<int>(); }

  // Pointer to MergingHooks object for NLO merging.
  MergingHooks*  mergingHooksPtr;

protected:

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the settings database.
  Settings*      settingsPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointer to Standard Model couplings.
  CoupSM*        coupSMPtr;

  // Pointers to the two incoming beams. Offset their location in event.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;
  int            beamOffset;

  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // Pointer to userHooks object for user interaction with program.
  UserHooks*     userHooksPtr;

  // Pointer to assign space-time vertices during parton evolution.
  PartonVertex*  partonVertexPtr;

  // Store uncertainty variations relevant to TimeShower.
  bool   doUncertainties, uVarMuSoftCorr, uVarMPIshowers, noResVariations,
         noProcVariations;
  int    nUncertaintyVariations, nVarQCD, uVarNflavQ;
  double dASmax, cNSpTmin, uVarpTmin2, overFactor;
  map<int,double> varG2GGmuRfac, varQ2QGmuRfac, varG2QQmuRfac, varX2XGmuRfac,
                  varG2GGcNS, varQ2QGcNS, varG2QQcNS, varX2XGcNS;
  map<int,double>* varPDFplus;
  map<int,double>* varPDFminus;
  map<int,double>* varPDFmember;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_TimeShower_H
