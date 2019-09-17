// TauDecays.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the TauDecays class.

#ifndef Pythia8_TauDecays_H
#define Pythia8_TauDecays_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/HelicityBasics.h"
#include "Pythia8/HelicityMatrixElements.h"
#include "Pythia8/PythiaComplex.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// TauDecays class.
// This class decays tau leptons, with helicity information.

class TauDecays {

public:

  // Constructor and destructor.
  TauDecays() : correlated(), tauExt(), tauMode(), tauMother(), tauPol(),
    hardME(), decayME(), infoPtr(), settingsPtr(), particleDataPtr(),
    rndmPtr(), couplingsPtr(), tau0Max(), tauMax(), rMax(), xyMax(), zMax(),
    limitTau0(), limitTau(), limitRadius(), limitCylinder(), limitDecay() {};
  ~TauDecays() {}

  // Initializer.
  void init(Info* infoPtrIn, Settings* settingsPtrIn,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    Couplings* couplingsPtrIn);

  // Decay a tau or correlated tau pair.
  bool decay(int iDec, Event& event);

  // Determine internal or external polarization and correlation mechanism.
  bool internalMechanism(Event &event);
  bool externalMechanism(Event &event);

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
  bool correlated;

  // User selected mode and mother for tau decays.
  int tauExt, tauMode, tauMother, tauPol;

  // Helicity matrix element pointers.
  HelicityMatrixElement* hardME;
  HelicityMatrixElement* decayME;

  // Hard process helicity matrix elements.
  HMETwoFermions2W2TwoFermions      hmeTwoFermions2W2TwoFermions;
  HMETwoFermions2GammaZ2TwoFermions hmeTwoFermions2GammaZ2TwoFermions;
  HMEW2TwoFermions                  hmeW2TwoFermions;
  HMEZ2TwoFermions                  hmeZ2TwoFermions;
  HMEGamma2TwoFermions              hmeGamma2TwoFermions;
  HMEHiggs2TwoFermions              hmeHiggs2TwoFermions;

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
