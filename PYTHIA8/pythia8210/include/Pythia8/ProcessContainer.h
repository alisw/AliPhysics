// ProcessContainer.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the collected machinery of a process.
// ProcessContainer: contains information on a particular process.
// SetupContainers: administrates the selection/creation of processes.

#ifndef Pythia8_ProcessContainer_H
#define Pythia8_ProcessContainer_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonDistributions.h"
#include "Pythia8/PhaseSpace.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/ResonanceDecays.h"
#include "Pythia8/Settings.h"
#include "Pythia8/SigmaProcess.h"
#include "Pythia8/SigmaTotal.h"
#include "Pythia8/SigmaOnia.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SusyCouplings.h"
#include "Pythia8/SLHAinterface.h"
#include "Pythia8/UserHooks.h"

namespace Pythia8 {

//==========================================================================

// The ProcessContainer class combines pointers to matrix element and
// phase space generator with general generation info.

class ProcessContainer {

public:

  // Constructor.
  ProcessContainer(SigmaProcess* sigmaProcessPtrIn = 0,
    bool externalPtrIn = false, PhaseSpace* phaseSpacePtrIn = 0) :
      sigmaProcessPtr(sigmaProcessPtrIn),
      externalPtr(externalPtrIn), phaseSpacePtr(phaseSpacePtrIn) {}

  // Destructor. Do not destroy external sigmaProcessPtr.
  ~ProcessContainer() {delete phaseSpacePtr;
    if (!externalPtr) delete sigmaProcessPtr;}

  // Initialize phase space and counters.
  bool init(bool isFirst, Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn, BeamParticle* beamAPtr,
    BeamParticle* beamBPtr, Couplings* couplings, SigmaTotal* sigmaTotPtr,
    ResonanceDecays* resDecaysPtrIn, SLHAinterface* slhaInterfacePtr,
    UserHooks* userHooksPtr);

  // Store or replace Les Houches pointer.
  void setLHAPtr( LHAup* lhaUpPtrIn,  ParticleData* particleDataPtrIn = 0)
    {lhaUpPtr = lhaUpPtrIn;
    if (particleDataPtrIn != 0) particleDataPtr = particleDataPtrIn;
    if (sigmaProcessPtr != 0) sigmaProcessPtr->setLHAPtr(lhaUpPtr);
    if (phaseSpacePtr != 0) phaseSpacePtr->setLHAPtr(lhaUpPtr);}

  // Update the CM energy of the event.
  void newECM(double eCM) {phaseSpacePtr->newECM(eCM);}

  // Generate a trial event; accepted or not.
  bool trialProcess();

  // Pick flavours and colour flow of process.
  bool constructState();

  // Give the hard subprocess (with option for a second hard subprocess).
  bool constructProcess( Event& process, bool isHardest = true);

  // Give the Les Houches decay chain for external resonance.
  bool constructDecays( Event& process);

  // Do resonance decays.
  bool decayResonances( Event& process);

  // Accumulate statistics after user veto.
  void accumulate();

  // Reset statistics on events generated so far.
  void reset();

  // Process name and code, and the number of final-state particles.
  string name()        const {return sigmaProcessPtr->name();}
  int    code()        const {return sigmaProcessPtr->code();}
  int    nFinal()      const {return sigmaProcessPtr->nFinal();}
  bool   isSUSY()      const {return sigmaProcessPtr->isSUSY();}

  // Member functions for info on generation process.
  bool   newSigmaMax() const {return newSigmaMx;}
  double sigmaMax()    const {return sigmaMx;}
  long   nTried()      const {return nTry;}
  long   nSelected()   const {return nSel;}
  long   nAccepted()   const {return nAcc;}
  double weightSum()   const {return wtAccSum;}
  double sigmaSelMC()  {if (nTry > nTryStat) sigmaDelta(); return sigmaAvg;}
  double sigmaMC()     {if (nTry > nTryStat) sigmaDelta(); return sigmaFin;}
  double deltaMC()     {if (nTry > nTryStat) sigmaDelta(); return deltaFin;}

  // Some kinematics quantities.
  int    id1()         const {return sigmaProcessPtr->id(1);}
  int    id2()         const {return sigmaProcessPtr->id(2);}
  double x1()          const {return phaseSpacePtr->x1();}
  double x2()          const {return phaseSpacePtr->x2();}
  double Q2Fac()       const {return sigmaProcessPtr->Q2Fac();}
  double mHat()        const {return sqrtpos(phaseSpacePtr->sHat());}
  double pTHat()       const {return phaseSpacePtr->pTHat();}

  // Tell whether container is for Les Houches events.
  bool   isLHAContainer() const {return isLHA;}
  int    lhaStrategy()    const {return lhaStrat;}

  // Info on Les Houches events.
  int    codeLHASize()       const {return codeLHA.size();}
  int    subCodeLHA(int i)   const {return codeLHA[i];}
  long   nTriedLHA(int i)    const {return nTryLHA[i];}
  long   nSelectedLHA(int i) const {return nSelLHA[i];}
  long   nAcceptedLHA(int i) const {return nAccLHA[i];}

  // When two hard processes set or get info whether process is matched.
  void   isSame( bool isSameIn) { isSameSave = isSameIn;}
  bool   isSame()      const {return isSameSave;}

private:

  // Constants: could only be changed in the code itself.
  static const int N12SAMPLE, N3SAMPLE;

  // Pointer to the subprocess matrix element. Mark if external.
  SigmaProcess*    sigmaProcessPtr;
  bool             externalPtr;

  // Pointer to the phase space generator.
  PhaseSpace*      phaseSpacePtr;

  // Pointer to various information on the generation.
  Info*            infoPtr;

  // Pointer to the particle data table.
  ParticleData*    particleDataPtr;

  // Pointer to the random number generator.
  Rndm*            rndmPtr;

  // Pointer to ResonanceDecays object for sequential resonance decays.
  ResonanceDecays* resDecaysPtr;

  // Pointer to userHooks object for user interaction with program.
  UserHooks*       userHooksPtr;

  // Pointer to LHAup for generating external events.
  LHAup*           lhaUpPtr;

  // Possibility to modify Les Houches input.
  bool   matchInOut;
  int    idRenameBeams, setLifetime, setLeptonMass, idLep[3];
  double mRecalculate, mLep[3];

  // Info on process.
  bool   isLHA, isNonDiff, isResolved, isDiffA, isDiffB, isDiffC, isQCD3body,
         allowNegSig, isSameSave, increaseMaximum, canVetoResDecay;
  int    lhaStrat, lhaStratAbs;
  bool   useStrictLHEFscales;

  // Statistics on generation process. (Long integers just in case.)
  bool   newSigmaMx;
  long   nTry, nSel, nAcc, nTryStat;
  double sigmaMx, sigmaSgn, sigmaSum, sigma2Sum, sigmaNeg, sigmaAvg,
         sigmaFin, deltaFin, weightNow, wtAccSum;

  // Statistics for Les Houches event classification.
  vector<int> codeLHA;
  vector<long> nTryLHA, nSelLHA, nAccLHA;

  // Estimate integrated cross section and its uncertainty.
  void sigmaDelta();

};

//==========================================================================

// The SetupContainers class turns the list of user-requested processes
// into a vector of ProcessContainer objects, each with a process.

class SetupContainers {

public:

  // Constructor.
  SetupContainers() {}

  // Initialization assuming all necessary data already read.
  bool init(vector<ProcessContainer*>& containerPtrs, Info* infoPtr,
    Settings& settings, ParticleData* particleDataPtr, Couplings* couplings);

  // Initialization of a second hard process.
  bool init2(vector<ProcessContainer*>& container2Ptrs, Settings& settings);

private:

  // Methods to check that outgoing SUSY particles are allowed ones.
  void setupIdVecs( Settings& settings);
  bool allowIdVals( int idCheck1, int idCheck2);

  // Arrays of allowed outgoing SUSY particles and their lengths.
  vector<int> idVecA, idVecB;
  int nVecA, nVecB;

  // Helper class to setup onia production.
  SigmaOniaSetup charmonium, bottomonium;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ProcessContainer_H
