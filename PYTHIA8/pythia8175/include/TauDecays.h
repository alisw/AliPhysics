// TauDecays.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the TauDecays class.

#ifndef Pythia8_TauDecays_H
#define Pythia8_TauDecays_H

#include "Basics.h"
#include "Event.h"
#include "HelicityBasics.h"
#include "HelicityMatrixElements.h"
#include "PythiaComplex.h"
#include "PythiaStdlib.h"
#include "Settings.h"

namespace Pythia8 {

//==========================================================================

// TauDecays class.
// This class decays tau leptons, with helicity information.

class TauDecays {
public:

  // Constructor and destructor.
  TauDecays() :
  correlated(0), 
    tauMode(0), tauMother(0), tauModeSave(0), tauMotherSave(0),
    polarization(0), polSave(0),
    hardME(0), decayME(0), 
    infoPtr(0), settingsPtr(0), particleDataPtr(0), rndmPtr(0), couplingsPtr(0),
    tau0Max(0), tauMax(0), rMax(0), xyMax(0), zMax(0),
    limitTau0(0), limitTau(0), limitRadius(0), limitCylinder(0), limitDecay(0)
  {};

  ~TauDecays() {}

  // Initializer.
  void init(Info* infoPtrIn, Settings* settingsPtrIn, 
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    Couplings* couplingsPtrIn);

  // Decay a tau or correlated tau pair.
  bool decay(int iDec, Event& event);

  // Choose a decay channel for a particle.
  vector<HelicityParticle> createChildren(HelicityParticle parent);

  // Perform an N-body isotropic decay.
  void isotropicDecay(vector<HelicityParticle>& p);

  // Write the decay to event record.
  void writeDecay(Event& event, vector<HelicityParticle>& p);

private: 

  // Constants: could only be changed in the code itself.
  static const int    NTRYCHANNEL, NTRYDECAY;
  static const double WTCORRECTION[11];

  // Flag whether a correlated tau decay should be performed.
  bool   correlated;

  // User selected mode and mother for tau decays.
  int    tauMode, tauMother, tauModeSave, tauMotherSave;

  // User selected polarization for tau decays.
  double polarization, polSave;

  // Helicity matrix element pointers.
  HelicityMatrixElement* hardME;
  HelicityMatrixElement* decayME;

  // Hard process helicity matrix elements.
  HMETwoFermions2W2TwoFermions      hmeTwoFermions2W2TwoFermions;
  HMETwoFermions2Z2TwoFermions      hmeTwoFermions2Z2TwoFermions;
  HMETwoFermions2Gamma2TwoFermions  hmeTwoFermions2Gamma2TwoFermions;
  HMETwoFermions2GammaZ2TwoFermions hmeTwoFermions2GammaZ2TwoFermions;
  HMEZ2TwoFermions                  hmeZ2TwoFermions;
  HMEHiggsEven2TwoFermions          hmeHiggsEven2TwoFermions;
  HMEHiggsOdd2TwoFermions           hmeHiggsOdd2TwoFermions;
  HMEHiggsCharged2TwoFermions       hmeHiggsCharged2TwoFermions;
  HMEUnpolarized                    hmeUnpolarized;

  // Tau decay helicity matrix elements.
  HMETau2Meson                    hmeTau2Meson;
  HMETau2TwoLeptons               hmeTau2TwoLeptons;
  HMETau2TwoMesonsViaVector       hmeTau2TwoMesonsViaVector;
  HMETau2TwoMesonsViaVectorScalar hmeTau2TwoMesonsViaVectorScalar;
  HMETau2ThreePions               hmeTau2ThreePions;
  HMETau2ThreeMesonsWithKaons     hmeTau2ThreeMesonsWithKaons;
  HMETau2ThreeMesonsGeneric       hmeTau2ThreeMesonsGeneric;
  HMETau2TwoPionsGamma            hmeTau2TwoPionsGamma;
  HMETau2FourPions                hmeTau2FourPions;
  HMETau2FivePions                hmeTau2FivePions;
  HMETau2PhaseSpace               hmeTau2PhaseSpace;

  // Particles of the hard process.
  HelicityParticle in1, in2, mediator, out1, out2;
  vector<HelicityParticle> particles;

  // The info pointer for the Pythia class.
  Info*         infoPtr;

  // Pointer to settings database.
  Settings*     settingsPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Pointer to SM coupling data.
  Couplings*    couplingsPtr;

  // Parameters to determine if correlated partner should decay.
  double tau0Max, tauMax, rMax, xyMax, zMax;
  bool limitTau0, limitTau, limitRadius, limitCylinder, limitDecay;
};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_TauDecays_H
