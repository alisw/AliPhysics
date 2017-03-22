// TauDecays.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the TauDecays class.

#include "TauDecays.h"

namespace Pythia8 {

//==========================================================================

// The TauDecays class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.
  
// Number of times to try a decay channel.
const int    TauDecays::NTRYCHANNEL = 10;
  
// Number of times to try a decay sampling.
  const int    TauDecays::NTRYDECAY   = 10000;
  //const int    TauDecays::NTRYDECAY   = 100000;

// These numbers are hardwired empirical parameters, 
// intended to speed up the M-generator.
const double TauDecays::WTCORRECTION[11] = { 1., 1., 1., 
  2., 5., 15., 60., 250., 1250., 7000., 50000. };

//--------------------------------------------------------------------------

// Initialize the TauDecays class with the necessary pointers to info,
// particle data, random numbers, and Standard Model couplings.
// Additionally, the necessary matrix elements are initialized with the
// Standard Model couplings, and particle data pointers.

void TauDecays::init(Info* infoPtrIn, Settings* settingsPtrIn, 
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn, 
  Couplings* couplingsPtrIn) {

  // Set the pointers.
  infoPtr         = infoPtrIn; 
  settingsPtr     = settingsPtrIn; 
  particleDataPtr = particleDataPtrIn; 
  rndmPtr         = rndmPtrIn; 
  couplingsPtr    = couplingsPtrIn;

  // Initialize the hard matrix elements.
  hmeTwoFermions2W2TwoFermions   .initPointers(particleDataPtr, couplingsPtr);
  hmeTwoFermions2Z2TwoFermions   .initPointers(particleDataPtr, couplingsPtr);
  hmeTwoFermions2Gamma2TwoFermions .initPointers(particleDataPtr, 
                                                 couplingsPtr);
  hmeTwoFermions2GammaZ2TwoFermions.initPointers(particleDataPtr, 
                                                 couplingsPtr);
  hmeZ2TwoFermions               .initPointers(particleDataPtr, couplingsPtr);
  hmeHiggsEven2TwoFermions       .initPointers(particleDataPtr, couplingsPtr);
  hmeHiggsOdd2TwoFermions        .initPointers(particleDataPtr, couplingsPtr);
  hmeHiggsCharged2TwoFermions    .initPointers(particleDataPtr, couplingsPtr);
  hmeUnpolarized                 .initPointers(particleDataPtr, couplingsPtr);

  // Initialize the tau decay matrix elements.
  hmeTau2Meson                   .initPointers(particleDataPtr, couplingsPtr);
  hmeTau2TwoLeptons              .initPointers(particleDataPtr, couplingsPtr);
  hmeTau2TwoMesonsViaVector      .initPointers(particleDataPtr, couplingsPtr);
  hmeTau2TwoMesonsViaVectorScalar.initPointers(particleDataPtr, couplingsPtr);
  hmeTau2ThreePions              .initPointers(particleDataPtr, couplingsPtr);
  hmeTau2ThreeMesonsWithKaons    .initPointers(particleDataPtr, couplingsPtr);
  hmeTau2ThreeMesonsGeneric      .initPointers(particleDataPtr, couplingsPtr);
  hmeTau2TwoPionsGamma           .initPointers(particleDataPtr, couplingsPtr);
  hmeTau2FourPions               .initPointers(particleDataPtr, couplingsPtr);
  hmeTau2FivePions               .initPointers(particleDataPtr, couplingsPtr);
  hmeTau2PhaseSpace              .initPointers(particleDataPtr, couplingsPtr);

  // User selected tau decay mode.
  tauModeSave   = settingsPtr->mode("ParticleDecays:sophisticatedTau");
    
  // User selected tau decay mother.
  tauMotherSave = settingsPtr->mode("ParticleDecays:tauMother");

  // User selected tau polarization.
  polSave       = settingsPtr->parm("ParticleDecays:tauPolarization");

  // Parameters to determine if correlated partner should decay.
  limitTau0     = settingsPtr->flag("ParticleDecays:limitTau0");
  tau0Max       = settingsPtr->parm("ParticleDecays:tau0Max");
  limitTau      = settingsPtr->flag("ParticleDecays:limitTau");
  tauMax        = settingsPtr->parm("ParticleDecays:tauMax");
  limitRadius   = settingsPtr->flag("ParticleDecays:limitRadius");
  rMax          = settingsPtr->parm("ParticleDecays:rMax");
  limitCylinder = settingsPtr->flag("ParticleDecays:limitCylinder");
  xyMax         = settingsPtr->parm("ParticleDecays:xyMax");
  zMax          = settingsPtr->parm("ParticleDecays:zMax");
  limitDecay    = limitTau0 || limitTau || limitRadius || limitCylinder;
}

//--------------------------------------------------------------------------

// Main method of the TauDecays class. Pass the index of the tau requested
// to be decayed along with the event record in which the tau exists. The
// tau is then decayed with proper spin correlations as well any partner. 
// The children of the decays are written to the event record, and if the 
// decays were succesful, a return value of true is supplied.

bool TauDecays::decay(int idxOut1, Event& event) {

  // User selected tau decay mode, mother, and polarization.
  tauMode      = tauModeSave;
  tauMother    = tauMotherSave;
  polarization = polSave;

  // Set the first outgoing particle of the hard process.
  out1 = HelicityParticle(event[idxOut1]); 
  out1.idx = idxOut1;

  // Begin PS April 2012.
  // Check if this tau already has helicity information (eg from LHEF).
  bool   hasHelicity = false;
  double helicityNow = 0.;
  if (tauMode >= 1 && abs(out1.pol()) <= 1.001) {
    hasHelicity = true;
    helicityNow = out1.pol();
  }
  // End PS April 2012.
  
  // Find the mediator of the hard process. Create temporary copy.
  int idxMediator  = out1.mother1();
  int idxFirstOut1 = idxOut1;
  while(idxMediator > 0 && event[idxMediator].id() == out1.id()) {
    idxFirstOut1   = idxMediator;
    idxMediator    = event[idxMediator].mother1();
  }
  Particle medTmp  = event[idxMediator];

  // Find and set up the incoming particles of the hard process.
  int idxIn1 = medTmp.mother1();
  int idxIn2 = medTmp.mother2();
  while(idxIn1 > 0 && event[idxIn1].id() == medTmp.id()) {
    idxIn1   = event[idxIn1].mother1();
    idxIn2   = event[idxIn2].mother2();
  }
  in1           = HelicityParticle(event[idxIn1]); 
  in1.idx       = idxIn1; 
  in1.direction = -1;
  in2           = HelicityParticle(event[idxIn2]); 
  in2.idx       = idxIn2; 
  in2.direction = -1;
  
  // Find and set up the second outgoing particle of the hard process.
  int idxOut2 = (medTmp.daughter1() == idxFirstOut1)
    ? medTmp.daughter2() : medTmp.daughter1();
  while (idxOut2 > 0 && event[idxOut2].daughter1() != 0 
    && event[event[idxOut2].daughter1()].id() == event[idxOut2].id()) {
    idxOut2 = event[idxOut2].daughter1();
  }
  out2     = HelicityParticle(event[idxOut2]); 
  out2.idx = idxOut2;

  // Set up the mediator. Special case for dipole shower, 
  // where a massless photon can branch to a tau pair.
  if (medTmp.id() == 22 && out2.idAbs() == 15 
    && medTmp.m() < out1.m() + out2.m()) {
    Vec4 pTmp        = out1.p() + out2.p();
    medTmp.p( pTmp);
    medTmp.m( pTmp.mCalc() );
  } 
  mediator           = HelicityParticle(medTmp);
  mediator.idx       = idxMediator; 
  mediator.direction = -1;

  // Set the particles vector.
  particles.clear();
  particles.push_back(in1);
  particles.push_back(in2);
  particles.push_back(out1);
  particles.push_back(out2);
  
  // Set the hard matrix element.
  // Polarized tau (decayed one by one).
  if (hasHelicity) {
    correlated = false;

  // Produced from a W.
  } else if (abs(mediator.id()) == 24) {
    // Produced from quarks: s-channel.
    if (abs(in1.id()) <= 18 && abs(in2.id()) <= 18)
      hardME = hmeTwoFermions2W2TwoFermions.initChannel(particles);
    // Produced from quarks: t-channel. 
    else if (abs(in1.id()) <= 18 || abs(in2.id()) <= 18) {
      bool fermion = (abs(in1.id()) <= 18) ? 0 : 1;
      particles[!fermion] 
        = (event[particles[fermion].daughter1()].id() == mediator.id())
        ? HelicityParticle(event[particles[fermion].daughter2()]) 
        : HelicityParticle(event[particles[fermion].daughter1()]);
      particles[!fermion].direction = 1;
      if (abs(particles[!fermion].id()) <= 18)
        hardME = hmeTwoFermions2W2TwoFermions.initChannel(particles);
      else {
        infoPtr->errorMsg("Warning in TauDecays::decay: unknown "
          "tau production, assuming unpolarized and uncorrelated");
	hardME = hmeUnpolarized.initChannel(particles);
      }
    // Unknown W production: assume negative helicity.
    } else if (tauMode == 1) {
      tauMode      = 3;
      polarization = -1;
    }
    correlated = false;

  // Produced from a photon.
  } else if (abs(mediator.id()) == 22 && abs(in1.id()) <= 18) {
    particles.push_back(mediator);
    hardME = hmeTwoFermions2Gamma2TwoFermions.initChannel(particles);
    correlated = true;

  // Produced from a photon/Z.
  } else if (abs(mediator.id()) == 23 && abs(in1.id()) <= 18) {
    particles.push_back(mediator);
    if (settingsPtr->mode("WeakZ0:gmZmode") == 0)
      hardME = hmeTwoFermions2GammaZ2TwoFermions.initChannel(particles);
    else if (settingsPtr->mode("WeakZ0:gmZmode") == 1)
      hardME = hmeTwoFermions2Gamma2TwoFermions.initChannel(particles);
    else if (settingsPtr->mode("WeakZ0:gmZmode") == 2)
      hardME = hmeTwoFermions2Z2TwoFermions.initChannel(particles);
    correlated = true;

  // Unkown Z production: assume unpolarized Z.
  } else if (abs(mediator.id()) == 23) {
    particles[1] = mediator;
    hardME = hmeZ2TwoFermions.initChannel(particles);
    correlated = true;

  // Produced from a CP even Higgs.
  } else if (abs(mediator.id()) == 25 || abs(mediator.id()) == 35) {
    hardME = hmeHiggsEven2TwoFermions.initChannel(particles);
    correlated = true;

  // Produced from a CP odd Higgs.
  } else if (abs(mediator.id()) == 36) {
    hardME = hmeHiggsOdd2TwoFermions.initChannel(particles);
    correlated = true;

  // Produced from a charged Higgs.
  } else if (abs(mediator.id()) == 37) {
    hardME = hmeHiggsCharged2TwoFermions.initChannel(particles);
    correlated = false;

  // Produced from a D or B hadron decay with a single tau.
  } else if ((abs(mediator.id()) == 411 || abs(mediator.id()) == 431 
           || abs(mediator.id()) == 511 || abs(mediator.id()) == 521 
           || abs(mediator.id()) == 531 || abs(mediator.id()) == 541
           || (abs(mediator.id()) > 5100 && abs(mediator.id()) < 5600) ) 
           && abs(out2.id()) == 16) {
    int idBmother = (mediator.id() > 0) ? -5 : 5;
    if (abs(mediator.id()) > 5100) idBmother = -idBmother; 
    particles[0] = HelicityParticle(  idBmother, 0, 0, 0, 0, 0, 0, 0, 
      0., 0., 0., 0., 0., 0., particleDataPtr);
    particles[1] = HelicityParticle( -idBmother, 0, 0, 0, 0, 0, 0, 0, 
      0., 0., 0., 0., 0., 0., particleDataPtr);
    particles[0].idx = 0; 
    particles[1].idx = 1;

    // D or B meson decays into neutrino + tau + meson.
    if (mediator.daughter1() + 2 == mediator.daughter2()) {	
      particles[0].p(mediator.p());
      particles[1].direction = 1; 
      particles[1].id(-particles[1].id());
      particles[1].p(particles[0].p() - particles[2].p() - particles[3].p());
    }

    // D or B meson decays into neutrino + tau.
    else {
      particles[0].p(mediator.p()/2);
      particles[1].p(mediator.p()/2);
    }
    hardME = hmeTwoFermions2W2TwoFermions.initChannel(particles);
    correlated = false;

  // Produced from a virtual photon with correlated taus.
  } else if (abs(out1.id()) == 15 && abs(out2.id()) == 15) {
    particles.push_back(mediator);
    particles[0] = HelicityParticle(-1, 0, 0, 0, 0, 0, 0, 0, 
      mediator.p()/2, 0., 0., particleDataPtr);
    particles[1] = HelicityParticle(1, 0, 0, 0, 0, 0, 0, 0, 
      mediator.p()/2, 0., 0., particleDataPtr);
    particles[0].direction = -1; 
    particles[1].direction = -1;
    particles[0].idx = 0; 
    particles[1].idx = 0;
    hardME = hmeTwoFermions2Gamma2TwoFermions.initChannel(particles);
    correlated = true;

  // Produced from an unknown process, assume unpolarized and uncorrelated.
  } else {
    if (tauMode <= 1) 
    infoPtr->errorMsg("Warning in TauDecays::decay: unknown "
      "tau production, assuming unpolarized and uncorrelated");
    hardME = hmeUnpolarized.initChannel(particles);
    correlated = false;
  }

  // Check if correlated partner should decay.
  if (correlated) {
    // Check vertex is within limits.
    if (limitTau0 && out2.tau0() > tau0Max) correlated = false;
    else if (limitTau && out2.tau() > tauMax) correlated = false;
    else if (limitRadius && pow2(out2.xDec()) + pow2(out2.yDec())
      + pow2(out2.zDec()) > pow2(rMax)) correlated = false;
    else if (limitCylinder && (pow2(out2.xDec()) + pow2(out2.yDec())
      > pow2(xyMax) || abs(out2.zDec()) > zMax)) correlated = false;
    // Check partner can decay.
    else if (!out2.canDecay()) correlated = false;
    else if (!out2.mayDecay()) correlated = false;
    // Check partner is compatible with hard matrix element (only leptons).
    else if (out2.idAbs() < 11 || out2.idAbs() > 16) {
      infoPtr->errorMsg("Warning in TauDecays::decay: incompatible "
        "correlated partner in tau decay");
      correlated = false;
    }
    // Undecay correlated partner if already decayed.
    else if (!out2.isFinal()) event.undoDecay(out2.idx);
  }

  // Pick the first tau to decay.
  HelicityParticle* tau;
  int idx;
  if (correlated) idx = (rndmPtr->flat() < 0.5) ? 2 : 3;
  else idx = (abs(particles[2].id()) == 15) ? 2 : 3;
  tau = &particles[idx];

  // Calculate the density matrix and decay the tau.
  if ( (tauMode == 2 && abs(mediator.id()) == tauMother) || tauMode == 3 ) {
    tau->rho[0][0] = (tau->id() > 0) 
      ? (1 - polarization) / 2 : (1 + polarization) / 2;
    tau->rho[1][1] = (tau->id() > 0) 
      ? (1 + polarization) / 2 : (1 - polarization) / 2;
    correlated = false;
  }

  // Begin PS April 2012.
  // Else use tau helicity provided by event record (LHEF).
  else if (hasHelicity) {
    tau->rho[0][0] = (1. - helicityNow) / 2.;
    tau->rho[1][1] = (1. + helicityNow) / 2.;
  }
  // End PS April 2012.

  // Else compute density matrix according to matrix element.
  else
    hardME->calculateRho(idx, particles);
  vector<HelicityParticle> children = createChildren(*tau);
  if (children.size() == 0) return false;

  // Decay the first tau.
  bool accepted = false;
  int  tries    = 0;
  while (!accepted) {
    isotropicDecay(children);
    double decayWeight    = decayME->decayWeight(children);
    double decayWeightMax = decayME->decayWeightMax(children);
    accepted = (rndmPtr->flat() < decayWeight / decayWeightMax);
    if (decayWeight > decayWeightMax)
      infoPtr->errorMsg("Warning in TauDecays::decay: maximum "
	"decay weight exceeded in tau decay");
    if (tries > NTRYDECAY) {
      infoPtr->errorMsg("Warning in TauDecays::decay: maximum "
	"number of decay attempts exceeded");
      break;
    }
    ++tries;
  }
  writeDecay(event,children);

  // If a correlated second tau exists, decay that tau as well.
  if (correlated) {
    idx = (idx == 2) ? 3 : 2;
    // Calculate the first tau decay matrix.
    decayME->calculateD(children);
    // Update the decay matrix for the tau.
    (*tau).D = children[0].D;
    // Switch the taus.
    tau = &particles[idx];
    // Calculate second tau's density matrix.
    hardME->calculateRho(idx, particles);

    // Decay the second tau.
    children.clear();
    children = createChildren(*tau);
    if (children.size() == 0) return false;
    accepted = false;
    tries    = 0;
    while (!accepted) {
      isotropicDecay(children);
      double decayWeight    = decayME->decayWeight(children);
      double decayWeightMax = decayME->decayWeightMax(children);
      accepted = (rndmPtr->flat() < decayWeight / decayWeightMax);
      if (decayWeight > decayWeightMax)
	infoPtr->errorMsg("Warning in TauDecays::decay: maximum "
	  "decay weight exceeded in correlated tau decay"); 
      if (tries > NTRYDECAY) {
	infoPtr->errorMsg("Warning in TauDecays::decay: maximum "
	  "number of decay attempts exceeded");
	break;
      }
      ++tries;
    }
    writeDecay(event,children);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Given a HelicityParticle parent, select the decay channel and return
// a vector of HelicityParticles containing the children, with the parent
// particle duplicated in the first entry of the vector.

vector<HelicityParticle> TauDecays::createChildren(HelicityParticle parent) {

  // Initial values.
  int meMode = 0;
  vector<HelicityParticle> children;

  // Set the parent as incoming.
  parent.direction = -1;

  // Setup decay data for the decaying particle.
  ParticleDataEntry decayData = parent.particleDataEntry();

  // Initialize the decay data.
  if (!decayData.preparePick(parent.id())) return children;

  // Try to pick a decay channel.
  bool decayed = false;
  int decayTries = 0;
  while (!decayed && decayTries < NTRYCHANNEL) {

    // Pick a decay channel.
    DecayChannel decayChannel = decayData.pickChannel();
    meMode = decayChannel.meMode();
    int decayMult = decayChannel.multiplicity();

    // Select children masses.
    bool allowed = false;
    int channelTries = 0;
    while (!allowed && channelTries < NTRYCHANNEL) {
      children.resize(0);
      children.push_back(parent);
      for (int i = 0; i < decayMult; ++i) {
        // Grab child ID.
        int childId = decayChannel.product(i);  
	// Flip sign for anti-particle decay.
        if (parent.id() < 0 && particleDataPtr->hasAnti(childId))
	  childId = -childId;
	// Grab child mass.
        double childMass = particleDataPtr->mass(childId);  
	// Push back the child into the children vector.
	children.push_back( HelicityParticle(childId, 91, parent.idx, 0, 0, 
	  0, 0, 0, 0., 0., 0., 0.,childMass, 0., particleDataPtr) );
      }
	
      // Check there is enough phase space for decay.
      if (decayMult > 1) {
        double massDiff = parent.m();
        for (int i = 0; i < decayMult; ++i) massDiff -= children[i].m();
	// For now we just check kinematically available.
        if (massDiff > 0) { 
          allowed = true; 
          decayed = true; 
        }
      }	

      // End pick a channel.
      ++channelTries;
    }
    ++decayTries;
  }

  // Swap the children ordering for muons.
  if (parent.idAbs() == 13 && children.size() == 4 && meMode == 22)
    swap(children[1], children[3]);

  // Set the decay matrix element.  
  // Two body decays.
  if (children.size() == 3) {
    if (meMode == 1521) 
      decayME = hmeTau2Meson.initChannel(children);
    else decayME = hmeTau2PhaseSpace.initChannel(children);
  } 

  // Three body decays.
  else if (children.size() == 4) {
    // Leptonic decay.
    if (meMode == 1531 || (parent.idAbs() == 13 && meMode == 22))	
      decayME = hmeTau2TwoLeptons.initChannel(children);
    // Two meson decay via vector meson.
    else if (meMode == 1532)
      decayME = hmeTau2TwoMesonsViaVector.initChannel(children);
    // Two meson decay via vector or scalar meson.
    else if (meMode == 1533)
      decayME = hmeTau2TwoMesonsViaVectorScalar.initChannel(children);
    // Flat phase space.
    else decayME = hmeTau2PhaseSpace.initChannel(children);
  }

  // Four body decays.
  else if (children.size() == 5) {
    // Three pion CLEO decay.
    if (meMode == 1541)
      decayME = hmeTau2ThreePions.initChannel(children);
    // Three meson decay with one or more kaons decay.
    else if (meMode == 1542)
      decayME = hmeTau2ThreeMesonsWithKaons.initChannel(children);
    // Generic three meson decay.
    else if (meMode == 1543)
      decayME = hmeTau2ThreeMesonsGeneric.initChannel(children);
    // Two pions and photon decay.
    else if (meMode == 1544)
      decayME = hmeTau2TwoPionsGamma.initChannel(children);
    // Flat phase space.
    else decayME = hmeTau2PhaseSpace.initChannel(children);
  }

  // Five body decays.
  else if (children.size() == 6) {
    // Four pion Novosibirsk current.
    if (meMode == 1551)
      decayME = hmeTau2FourPions.initChannel(children);
    // Flat phase space.
    else decayME = hmeTau2PhaseSpace.initChannel(children);
  }

  // Six body decays.
  else if (children.size() == 7) {
    // Four pion Novosibirsk current.
    if (meMode == 1561)
      decayME = hmeTau2FivePions.initChannel(children);
    // Flat phase space.
    else decayME = hmeTau2PhaseSpace.initChannel(children);
  }

  // Flat phase space.
  else decayME = hmeTau2PhaseSpace.initChannel(children);

  // Done.
  return children;
}

//--------------------------------------------------------------------------

// N-body decay using the M-generator algorithm described in
// "Monte Carlo Phase Space" by F. James in CERN 68-15, May 1968. Taken
// from ParticleDecays::mGenerator but modified to handle spin particles.
// Given a vector of HelicityParticles where the first particle is 
// the mother, the remaining particles are decayed isotropically.

void TauDecays::isotropicDecay(vector<HelicityParticle>& children) {

  // Mother and sum daughter masses.
  int decayMult = children.size() - 1;
  double m0      = children[0].m();
  double mSum    = children[1].m();
  for (int i = 2; i <= decayMult; ++i) mSum += children[i].m(); 
  double mDiff   = m0 - mSum;
    
  // Begin setup of intermediate invariant masses.
  vector<double> mInv;
  for (int i = 0; i <= decayMult; ++i) mInv.push_back( children[i].m());
    
  // Calculate the maximum weight in the decay.
  double wtPS;
  double wtPSmax = 1. / WTCORRECTION[decayMult];
  double mMax    = mDiff + children[decayMult].m();
  double mMin    = 0.; 
  for (int i = decayMult - 1; i > 0; --i) {
    mMax        += children[i].m();
    mMin        += children[i+1].m();
    double mNow  = children[i].m();
    wtPSmax     *= 0.5 * sqrtpos( (mMax - mMin - mNow) * (mMax + mMin + mNow)
		 * (mMax + mMin - mNow) * (mMax - mMin + mNow) ) / mMax;  
  }
      
  // Begin loop to find the set of intermediate invariant masses.
  vector<double> rndmOrd;
  do {
    wtPS  = 1.;
	  
    // Find and order random numbers in descending order.
    rndmOrd.clear();
    rndmOrd.push_back(1.);
    for (int i = 1; i < decayMult - 1; ++i) {
      double random = rndmPtr->flat();
      rndmOrd.push_back(random);
      for (int j = i - 1; j > 0; --j) {
        if (random > rndmOrd[j]) swap( rndmOrd[j], rndmOrd[j+1] );
        else break;
      } 
    }
    rndmOrd.push_back(0.);
	  
    // Translate into intermediate masses and find weight.
    for (int i = decayMult - 1; i > 0; --i) {
      mInv[i] = mInv[i+1] + children[i].m() 
              + (rndmOrd[i-1] - rndmOrd[i]) * mDiff; 
      wtPS   *= 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - children[i].m()) 
	      * (mInv[i] + mInv[i+1] + children[i].m()) 
              * (mInv[i] + mInv[i+1] - children[i].m()) 
	      * (mInv[i] - mInv[i+1] + children[i].m()) ) / mInv[i];  
    }
	  
    // If rejected, try again with new invariant masses.
  } while ( wtPS < rndmPtr->flat() * wtPSmax ); 
	
  // Perform two-particle decays in the respective rest frame.
  vector<Vec4> pInv(decayMult + 1);
  for (int i = 1; i < decayMult; ++i) {
    double pAbs = 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - children[i].m()) 
                * (mInv[i] + mInv[i+1] + children[i].m()) 
                * (mInv[i] + mInv[i+1] - children[i].m())
                * (mInv[i] - mInv[i+1] + children[i].m()) ) / mInv[i]; 

    // Isotropic angles give three-momentum.
    double cosTheta = 2. * rndmPtr->flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi      = 2. * M_PI * rndmPtr->flat();
    double pX       = pAbs * sinTheta * cos(phi);  
    double pY       = pAbs * sinTheta * sin(phi);  
    double pZ       = pAbs * cosTheta;  

    // Calculate energies, fill four-momenta.
    double eHad     = sqrt( children[i].m()*children[i].m() + pAbs*pAbs);
    double eInv     = sqrt( mInv[i+1]*mInv[i+1] + pAbs*pAbs);
    children[i].p( pX, pY, pZ, eHad);
    pInv[i+1].p( -pX, -pY, -pZ, eInv);
  }       

  // Boost decay products to the mother rest frame.
  children[decayMult].p( pInv[decayMult] );
  for (int iFrame = decayMult - 1; iFrame > 1; --iFrame) 
    for (int i = iFrame; i <= decayMult; ++i) 
      children[i].bst( pInv[iFrame], mInv[iFrame]);

  // Boost decay products to the current frame. 
  pInv[1].p( children[0].p() );
  for (int i = 1; i <= decayMult; ++i) children[i].bst( pInv[1], mInv[1] );

  // Done.
  return;
}

//--------------------------------------------------------------------------

// Write the vector of HelicityParticles to the event record, excluding the
// first particle. Set the lifetime and production vertex of the particles
// and mark the first particle of the vector as decayed.
 
void TauDecays::writeDecay(Event& event, vector<HelicityParticle>& children) {

  // Set additional information and append children to event.
  int  decayMult   = children.size() - 1;
  Vec4 decayVertex = children[0].vDec();  
  for (int i = 1; i <= decayMult; i++) {
    // Set child lifetime.
    children[i].tau(children[i].tau0() * rndmPtr->exp());
    // Set child production vertex.
    children[i].vProd(decayVertex);
    // Append child to record.
    children[i].idx = event.append(children[i]);
  }
  
  // Mark the parent as decayed and set children.
  event[children[0].idx].statusNeg(); 
  event[children[0].idx].daughters(children[1].idx, children[decayMult].idx);

}

//==========================================================================

} // end namespace Pythia8
