// UserHooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file to allow user access to program at different stages.
// UserHooks: almost empty base class, with user to write the rela code.
// MyUserHooks: derived class, only intended as an example.

#ifndef Pythia8_UserHooks_H
#define Pythia8_UserHooks_H

#include "Pythia8/Event.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// Forward references to the PhaseSpace and StringEnd classes.
class PhaseSpace;
class StringEnd;

//==========================================================================

// UserHooks is base class for user access to program execution.

class UserHooks {

public:

  // Destructor.
  virtual ~UserHooks() {}

  // Initialize pointers and workEvent. Note: not virtual.
  void initPtr( Info* infoPtrIn, Settings* settingsPtrIn,
    ParticleData* particleDataPtrIn,  Rndm* rndmPtrIn,
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    BeamParticle* beamPomAPtrIn, BeamParticle* beamPomBPtrIn,
    CoupSM* coupSMPtrIn, PartonSystems* partonSystemsPtrIn,
    SigmaTotal* sigmaTotPtrIn) { infoPtr = infoPtrIn;
    settingsPtr = settingsPtrIn; particleDataPtr = particleDataPtrIn;
    rndmPtr = rndmPtrIn; beamAPtr = beamAPtrIn; beamBPtr = beamBPtrIn;
    beamPomAPtr = beamPomAPtrIn; beamPomBPtr = beamPomBPtrIn;
    coupSMPtr = coupSMPtrIn; partonSystemsPtr = partonSystemsPtrIn;
    sigmaTotPtr = sigmaTotPtrIn;
    workEvent.init("(work event)", particleDataPtr);}

  // Initialisation after beams have been set by Pythia::init().
  virtual bool initAfterBeams() { return true; }

  // Possibility to modify cross section of process.
  virtual bool canModifySigma() {return false;}

  // Multiplicative factor modifying the cross section of a hard process.
  virtual double multiplySigmaBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool inEvent);

  // Possibility to bias selection of events, compensated by a weight.
  virtual bool canBiasSelection() {return false;}

  // Multiplicative factor in the phase space selection of a hard process.
  virtual double biasSelectionBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool inEvent);

  // Event weight to compensate for selection weight above.
  virtual double biasedSelectionWeight() {return 1./selBias;}

  // Possibility to veto event after process-level selection.
  virtual bool canVetoProcessLevel() {return false;}

  // Decide whether to veto current process or not, based on process record.
  // Usage: doVetoProcessLevel( process).
  virtual bool doVetoProcessLevel(Event& ) {return false;}

  // Possibility to veto resonance decay chain.
  virtual bool canVetoResonanceDecays() {return false;}

  // Decide whether to veto current resonance decay chain or not, based on
  // process record. Usage: doVetoProcessLevel( process).
  virtual bool doVetoResonanceDecays(Event& ) {return false;}

  // Possibility to veto MPI + ISR + FSR evolution and kill event,
  // making decision at a fixed pT scale. Useful for MLM-style matching.
  virtual bool canVetoPT() {return false;}

  // Transverse-momentum scale for veto test.
  virtual double scaleVetoPT() {return 0.;}

  // Decide whether to veto current event or not, based on event record.
  // Usage: doVetoPT( iPos, event), where iPos = 0: no emissions so far;
  // iPos = 1/2/3 joint evolution, latest step was MPI/ISR/FSR;
  // iPos = 4: FSR only afterwards; iPos = 5: FSR in resonance decay.
  virtual bool doVetoPT( int , const Event& ) {return false;}

  // Possibility to veto MPI + ISR + FSR evolution and kill event,
  // making decision after fixed number of ISR or FSR steps.
  virtual bool canVetoStep() {return false;}

  // Up to how many ISR + FSR steps of hardest interaction should be checked.
  virtual int numberVetoStep() {return 1;}

  // Decide whether to veto current event or not, based on event record.
  // Usage: doVetoStep( iPos, nISR, nFSR, event), where iPos as above,
  // nISR and nFSR number of emissions so far for hard interaction only.
  virtual bool doVetoStep( int , int , int , const Event& ) {return false;}

  // Possibility to veto MPI + ISR + FSR evolution and kill event,
  // making decision after fixed number of MPI steps.
  virtual bool canVetoMPIStep() {return false;}

  // Up to how many MPI steps should be checked.
  virtual int numberVetoMPIStep() {return 1;}

  // Decide whether to veto current event or not, based on event record.
  // Usage: doVetoMPIStep( nMPI, event), where nMPI is number of MPI's so far.
  virtual bool doVetoMPIStep( int , const Event& ) {return false;}

  // Possibility to veto event after ISR + FSR + MPI in parton level,
  // but before beam remnants and resonance decays.
  virtual bool canVetoPartonLevelEarly() {return false;}

  // Decide whether to veto current partons or not, based on event record.
  // Usage: doVetoPartonLevelEarly( event).
  virtual bool doVetoPartonLevelEarly( const Event& ) {return false;}

  // Retry same ProcessLevel with a new PartonLevel after a veto in
  // doVetoPT, doVetoStep, doVetoMPIStep or doVetoPartonLevelEarly
  // if you overload this method to return true.
  virtual bool retryPartonLevel() {return false;}

  // Possibility to veto event after parton-level selection.
  virtual bool canVetoPartonLevel() {return false;}

  // Decide whether to veto current partons or not, based on event record.
  // Usage: doVetoPartonLevel( event).
  virtual bool doVetoPartonLevel( const Event& ) {return false;}

  // Possibility to set initial scale in TimeShower for resonance decay.
  virtual bool canSetResonanceScale() {return false;}

  // Initial scale for TimeShower evolution.
  // Usage: scaleResonance( iRes, event), where iRes is location
  // of decaying resonance in the event record.
  virtual double scaleResonance( int, const Event& ) {return 0.;}

  // Possibility to veto an emission in the ISR machinery.
  virtual bool canVetoISREmission() {return false;}

  // Decide whether to veto current emission or not, based on event record.
  // Usage: doVetoISREmission( sizeOld, event, iSys) where sizeOld is size
  // of event record before current emission-to-be-scrutinized was added,
  // and iSys is the system of the radiation (according to PartonSystems).
  virtual bool doVetoISREmission( int, const Event&, int ) {return false;}

  // Possibility to veto an emission in the FSR machinery.
  virtual bool canVetoFSREmission() {return false;}

  // Decide whether to veto current emission or not, based on event record.
  // Usage: doVetoFSREmission( sizeOld, event, iSys, inResonance) where
  // sizeOld is size of event record before current emission-to-be-scrutinized
  // was added, iSys is the system of the radiation (according to
  // PartonSystems), and inResonance is true if the emission takes place in a
  // resonance decay.
  virtual bool doVetoFSREmission( int, const Event&, int, bool = false )
      {return false;}

  // Possibility to veto an MPI.
  virtual bool canVetoMPIEmission() { return false; }

  // Decide whether to veto an MPI based on event record.
  // Usage: doVetoMPIEmission( sizeOld, event) where sizeOld
  // is size of event record before the current MPI.
  virtual bool doVetoMPIEmission( int, const Event &) { return false; }

  // Possibility to reconnect colours from resonance decay systems.
  virtual bool canReconnectResonanceSystems() { return false; }

  // Do reconnect colours from resonance decay systems.
  // Usage: doVetoFSREmission( oldSizeEvt, event)
  // where oldSizeEvent is the event size before resonance decays.
  // Should normally return true, while false means serious failure.
  // Value of PartonLevel:earlyResDec determines where method is called.
  virtual bool doReconnectResonanceSystems( int, Event &) {return true;}

  // Enhance emission rates (sec. 4 in EPJC (2013) 73).
  virtual bool canEnhanceEmission() {return false;}
  virtual double enhanceFactor( string ) {return 1.;}
  virtual double vetoProbability( string ) {return 0.;}
  void setEnhancedEventWeight(double wt) { enhancedEventWeight = wt;}
  double getEnhancedEventWeight() { return enhancedEventWeight;}

  // Bookkeeping of weights for enhanced actual or trial emissions
  // (sec. 3 in EPJC (2013) 73).
  virtual bool canEnhanceTrial() {return false;}
  void setEnhancedTrial( double pTIn, double wtIn) { pTEnhanced = pTIn;
    wtEnhanced = wtIn; }
  double getEnhancedTrialPT() { return pTEnhanced;}
  double getEnhancedTrialWeight() { return wtEnhanced;}

  // Can change fragmentation parameters.
  virtual bool canChangeFragPar() { return false;}

  // Set initial ends of a string to be fragmented. This is done once
  // for each string. Note that the second string end may be zero in case
  // we are hadronising a string piece leading to a junction.
  virtual void setStringEnds( const StringEnd*, const StringEnd*,
    vector<int>) {}

  // Do change fragmentation parameters.
  // Input: flavPtr, zPtr, pTPtr, idEnd, m2Had, iParton and posEnd (or
  // negEnd).
  virtual bool doChangeFragPar( StringFlav*, StringZ*, StringPT*, int,
    double, vector<int>, const StringEnd* ) { return false;}

  // Do a veto on a hadron just before it is added to the final state.
  // The StringEnd from which the the hadron was produced is included
  // for information.
  virtual bool doVetoFragmentation( Particle, const StringEnd *)
    { return false;}

  // Do a veto on a hadron just before it is added to the final state
  // (final two hadron case).
  virtual bool doVetoFragmentation(Particle, Particle,
    const StringEnd*, const StringEnd* ) { return false;}

  // Can set the overall impact parameter for the MPI treatment.
  virtual bool canSetImpactParameter() const { return false; }

  // Set the overall impact parameter for the MPI treatment.
  virtual double doSetImpactParameter() { return 0.0; }

protected:

  // Constructor.
  UserHooks() : infoPtr(0), settingsPtr(0), particleDataPtr(0), rndmPtr(0),
    beamAPtr(0), beamBPtr(0), beamPomAPtr(0), beamPomBPtr(0), coupSMPtr(0),
    partonSystemsPtr(0), sigmaTotPtr(0), selBias(1.), enhancedEventWeight(),
    pTEnhanced(), wtEnhanced() {}

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the settings database.
  Settings*      settingsPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

 // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointers to the two incoming beams and to Pomeron beam-inside-beam.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;
  BeamParticle*  beamPomAPtr;
  BeamParticle*  beamPomBPtr;

  // Pointers to Standard Model couplings.
  CoupSM*        coupSMPtr;

  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // Pointer to the total/elastic/diffractive cross sections.
  SigmaTotal*    sigmaTotPtr;

  // omitResonanceDecays omits resonance decay chains from process record.
  void omitResonanceDecays(const Event& process, bool finalOnly = false);

  // subEvent extracts currently resolved partons in the hard process.
  void subEvent(const Event& event, bool isHardest = true);

  // Have one event object around as work area.
  Event workEvent;

  // User-imposed selection bias.
  double selBias;

  // Bookkept quantities for boosted event weights.
  double enhancedEventWeight, pTEnhanced, wtEnhanced;

};

//==========================================================================

// SuppressSmallPT is a derived class for user access to program execution.
// It is a simple example, illustrating how to suppress the cross section
// of 2 -> 2 processes by a factor pT^4 / (pT0^2 + pT^2)^2, with pT0 input,
// and also modify alpha_strong scale similarly.

class SuppressSmallPT : public UserHooks {

public:

  // Constructor.
  SuppressSmallPT( double pT0timesMPIIn = 1., int numberAlphaSIn = 0,
    bool useSameAlphaSasMPIIn = true) : pT20(0.) {isInit = false;
    pT0timesMPI = pT0timesMPIIn; numberAlphaS = numberAlphaSIn;
    useSameAlphaSasMPI = useSameAlphaSasMPIIn;}

  // Possibility to modify cross section of process.
  virtual bool canModifySigma() {return true;}

  // Multiplicative factor modifying the cross section of a hard process.
  // Usage: inEvent is true for event generation, false for initialization.
  virtual double multiplySigmaBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool );

private:

  // Save input properties and the squared pT0 scale.
  bool   isInit, useSameAlphaSasMPI;
  int    numberAlphaS;
  double pT0timesMPI, pT20;

  // Alpha_strong calculation.
  AlphaStrong alphaS;

};

//==========================================================================

// UserHooksVector implements a vector of UserHooks and is itself a UserHooks.

class UserHooksVector: public UserHooks {

private:

  // The default constructor is private, and should only be used
  // internally in Pythia.
  UserHooksVector() {}
  friend class Pythia;

public:

  // Destructor.
  virtual ~UserHooksVector() {}

  // Initialisation after beams have been set by Pythia::init().
  // Check that there are no (obvious) clashes.
  virtual bool initAfterBeams() {
    int nCanSetResonanceScale  = 0;
    int nCanChangeFragPar      = 0;
    int nCanSetImpactParameter = 0;
    for ( int i = 0, N = hooks.size(); i < N; ++i ) {
      hooks[i]->initPtr(infoPtr, settingsPtr, particleDataPtr, rndmPtr,
                        beamAPtr, beamBPtr, beamPomAPtr, beamPomBPtr,
                        coupSMPtr, partonSystemsPtr, sigmaTotPtr);
      if ( !hooks[i]->initAfterBeams() ) return false;
      if (hooks[i]->canSetResonanceScale()) ++nCanSetResonanceScale;
      if (hooks[i]->canChangeFragPar()) ++nCanChangeFragPar;
      if (hooks[i]->canSetImpactParameter()) ++nCanSetImpactParameter;
    }
    if (nCanSetResonanceScale > 1) {
      infoPtr->errorMsg("Error in UserHooksVector::initAfterBeams "
        "multiple UserHooks with canSetResonanceScale() not allowed");
      return false;
    }
    if (nCanChangeFragPar > 1) {
      infoPtr->errorMsg("Error in UserHooksVector::initAfterBeams "
        "multiple UserHooks with canChangeFragPar() not allowed");
      return false;
    }
    if (nCanSetImpactParameter > 1) {
      infoPtr->errorMsg("Error in UserHooksVector::initAfterBeams "
        "multiple UserHooks with canSetImpactParameter() not allowed");
      return false;
    }
    return true;
  }

  // Possibility to modify cross section of process.
  virtual bool canModifySigma() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canModifySigma() ) return true;
    return false;
  }

  // Multiplicative factor modifying the cross section of a hard process.
  virtual double multiplySigmaBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool inEvent) {
    double f = 1.0;
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canModifySigma() )
        f *= hooks[i]->multiplySigmaBy(sigmaProcessPtr, phaseSpacePtr,inEvent);
    return f;
  }

  // Possibility to bias selection of events, compensated by a weight.
  virtual bool canBiasSelection() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canBiasSelection() ) return true;
    return false;
  }

  // Multiplicative factor in the phase space selection of a hard process.
  virtual double biasSelectionBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool inEvent) {
    double f = 1.0;
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canBiasSelection() )
        f *= hooks[i]->biasSelectionBy(sigmaProcessPtr, phaseSpacePtr,
             inEvent);
    return f;
  }

  // Event weight to compensate for selection weight above.
  virtual double biasedSelectionWeight() {
    double f = 1.0;
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canBiasSelection() )
        f *= hooks[i]->biasedSelectionWeight();
    return f;
  }

  // Possibility to veto event after process-level selection.
  virtual bool canVetoProcessLevel() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoProcessLevel() ) return true;
    return false;
  }

  // Decide whether to veto current process or not, based on process record.
  // Usage: doVetoProcessLevel( process).
  virtual bool doVetoProcessLevel(Event& e) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoProcessLevel() &&
           hooks[i]->doVetoProcessLevel(e) ) return true;
    return false;
  }

  // Possibility to veto resonance decay chain.
  virtual bool canVetoResonanceDecays() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoResonanceDecays() ) return true;
    return false;
  }

  // Decide whether to veto current resonance decay chain or not, based on
  // process record. Usage: doVetoProcessLevel( process).
  virtual bool doVetoResonanceDecays(Event& e) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoResonanceDecays() &&
           hooks[i]->doVetoResonanceDecays(e) ) return true;
    return false;
  }

  // Possibility to veto MPI + ISR + FSR evolution and kill event,
  // making decision at a fixed pT scale. Useful for MLM-style matching.
  virtual bool canVetoPT() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoPT() ) return true;
    return false;
  }

  // Transverse-momentum scale for veto test.
  virtual double scaleVetoPT() {
    double s = 0.0;
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoPT() ) s = max(s, hooks[i]->scaleVetoPT());
    return s;
  }

  // Decide whether to veto current event or not, based on event record.
  // Usage: doVetoPT( iPos, event), where iPos = 0: no emissions so far;
  // iPos = 1/2/3 joint evolution, latest step was MPI/ISR/FSR;
  // iPos = 4: FSR only afterwards; iPos = 5: FSR in resonance decay.
  virtual bool doVetoPT( int iPos, const Event& e) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoPT() && hooks[i]->doVetoPT(iPos, e) ) return true;
    return false;
  }

  // Possibility to veto MPI + ISR + FSR evolution and kill event,
  // making decision after fixed number of ISR or FSR steps.
  virtual bool canVetoStep() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoStep() ) return true;
    return false;
  }

  // Up to how many ISR + FSR steps of hardest interaction should be checked.
  virtual int numberVetoStep() {
    int n = 1;
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoStep() ) n = max(n, hooks[i]->numberVetoStep());
    return n;
  }

  // Decide whether to veto current event or not, based on event record.
  // Usage: doVetoStep( iPos, nISR, nFSR, event), where iPos as above,
  // nISR and nFSR number of emissions so far for hard interaction only.
  virtual bool doVetoStep( int iPos, int nISR, int nFSR, const Event& e) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoStep()
        && hooks[i]->doVetoStep(iPos, nISR, nFSR, e) ) return true;
    return false;
  }

  // Possibility to veto MPI + ISR + FSR evolution and kill event,
  // making decision after fixed number of MPI steps.
  virtual bool canVetoMPIStep() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoMPIStep() ) return true;
    return false;
  }

  // Up to how many MPI steps should be checked.
  virtual int numberVetoMPIStep() {
    int n = 1;
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoMPIStep() )
        n = max(n, hooks[i]->numberVetoMPIStep());
    return n;
  }

  // Decide whether to veto current event or not, based on event record.
  // Usage: doVetoMPIStep( nMPI, event), where nMPI is number of MPI's so far.
  virtual bool doVetoMPIStep( int nMPI, const Event& e) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoMPIStep() && hooks[i]->doVetoMPIStep(nMPI, e) )
        return true;
    return false;
  }

  // Possibility to veto event after ISR + FSR + MPI in parton level,
  // but before beam remnants and resonance decays.
  virtual bool canVetoPartonLevelEarly() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoPartonLevelEarly() ) return true;
    return false;
  }

  // Decide whether to veto current partons or not, based on event record.
  // Usage: doVetoPartonLevelEarly( event).
  virtual bool doVetoPartonLevelEarly( const Event& e) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoPartonLevelEarly()
        && hooks[i]->doVetoPartonLevelEarly(e) ) return true;
    return false;
  }

  // Retry same ProcessLevel with a new PartonLevel after a veto in
  // doVetoPT, doVetoStep, doVetoMPIStep or doVetoPartonLevelEarly
  // if you overload this method to return true.
  virtual bool retryPartonLevel() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->retryPartonLevel() ) return true;
    return false;
  }

  // Possibility to veto event after parton-level selection.
  virtual bool canVetoPartonLevel() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoPartonLevel() ) return true;
    return false;
  }

  // Decide whether to veto current partons or not, based on event record.
  // Usage: doVetoPartonLevel( event).
  virtual bool doVetoPartonLevel( const Event& e) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoPartonLevel()
        && hooks[i]->doVetoPartonLevel(e) ) return true;
   return false;
  }

  // Possibility to set initial scale in TimeShower for resonance decay.
  virtual bool canSetResonanceScale() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canSetResonanceScale() ) return true;
    return false;
  }

  // Initial scale for TimeShower evolution.
  // Usage: scaleResonance( iRes, event), where iRes is location
  // of decaying resonance in the event record.
  virtual double scaleResonance( int iRes, const Event& e) {
    double s = 0.0;
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canSetResonanceScale() )
        s = max(s, hooks[i]->scaleResonance(iRes, e));
    return s;
  }

  // Possibility to veto an emission in the ISR machinery.
  virtual bool canVetoISREmission() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoISREmission() ) return true;
    return false;
  }

  // Decide whether to veto current emission or not, based on event record.
  // Usage: doVetoISREmission( sizeOld, event, iSys) where sizeOld is size
  // of event record before current emission-to-be-scrutinized was added,
  // and iSys is the system of the radiation (according to PartonSystems).
  virtual bool doVetoISREmission( int sizeOld, const Event& e, int iSys) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoISREmission()
        && hooks[i]->doVetoISREmission(sizeOld, e, iSys) ) return true;
    return false;
  }

  // Possibility to veto an emission in the FSR machinery.
  virtual bool canVetoFSREmission() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoFSREmission() ) return true;
    return false;
  }

  // Decide whether to veto current emission or not, based on event record.
  // Usage: doVetoFSREmission( sizeOld, event, iSys, inResonance) where
  // sizeOld is size of event record before current emission-to-be-scrutinized
  // was added, iSys is the system of the radiation (according to
  // PartonSystems), and inResonance is true if the emission takes place in a
  // resonance decay.
  virtual bool doVetoFSREmission(int sizeOld, const Event& e,
    int iSys, bool inResonance = false ) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoFSREmission()
        && hooks[i]->doVetoFSREmission(sizeOld, e, iSys, inResonance) )
        return true;
    return false;
  }


  // Possibility to veto an MPI.
  virtual bool canVetoMPIEmission() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoMPIEmission() ) return true;
    return false;
  }

  // Decide whether to veto an MPI based on event record.
  // Usage: doVetoMPIEmission( sizeOld, event) where sizeOld
  // is size of event record before the current MPI.
  virtual bool doVetoMPIEmission( int sizeOld, const Event & e) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canVetoMPIEmission()
        && hooks[i]->doVetoMPIEmission(sizeOld, e) )
        return true;
    return false;
  }

  // Possibility to reconnect colours from resonance decay systems.
  virtual bool canReconnectResonanceSystems() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canReconnectResonanceSystems() ) return true;
    return false;
  }

  // Do reconnect colours from resonance decay systems.
  // Usage: doVetoFSREmission( oldSizeEvt, event)
  // where oldSizeEvent is the event size before resonance decays.
  // Should normally return true, while false means serious failure.
  // Value of PartonLevel:earlyResDec determines where method is called.
  virtual bool doReconnectResonanceSystems( int j, Event & e) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canReconnectResonanceSystems()
        && hooks[i]->doReconnectResonanceSystems(j, e) ) return true;
    return false;
  }

  // Enhance emission rates (sec. 4 in EPJC (2013) 73).
  virtual bool canEnhanceEmission() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canEnhanceEmission() ) return true;
    return false;
  }
  virtual double enhanceFactor( string s) {
    double f = 1.0;
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canEnhanceEmission() ) f *= hooks[i]->enhanceFactor(s);
    return f;
  }
  virtual double vetoProbability( string s) {
    double keep = 1.0;
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canEnhanceEmission() )
        keep *= 1.0 - hooks[i]->vetoProbability(s);
    return 1.0 - keep;
  }

  // Bookkeeping of weights for enhanced actual or trial emissions
  // (sec. 3 in EPJC (2013) 73).
  virtual bool canEnhanceTrial() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canEnhanceTrial() ) return true;
    return false;
  }

  // Can change fragmentation parameters.
  virtual bool canChangeFragPar() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canChangeFragPar() ) return true;
    return false;
  }

  // Do a veto on a hadron just before it is added to the final state.
  virtual bool doVetoFragmentation(Particle p, const StringEnd* nowEnd) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canChangeFragPar()
        && hooks[i]->doVetoFragmentation(p, nowEnd) ) return true;
    return false;
  }

  virtual bool doVetoFragmentation(Particle p1, Particle p2,
    const StringEnd* e1, const StringEnd* e2) {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canChangeFragPar()
        && hooks[i]->doVetoFragmentation(p1, p2, e1, e2) ) return true;
    return false;
  }

  // Can set the overall impact parameter for the MPI treatment.
  virtual bool canSetImpactParameter() const {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canSetImpactParameter() ) return true;
    return false;
  }

  // Set the overall impact parameter for the MPI treatment.
  virtual double doSetImpactParameter() {
    for ( int i = 0, N = hooks.size(); i < N; ++i )
      if ( hooks[i]->canSetImpactParameter() )
        return hooks[i]->doSetImpactParameter();
    return 0.0;
  }

public:

  vector<UserHooks*> hooks;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_UserHooks_H
