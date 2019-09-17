// TauDecays.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the TauDecays class.

#include "Pythia8/TauDecays.h"

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
  hmeTwoFermions2W2TwoFermions
    .initPointers(particleDataPtr, couplingsPtr, settingsPtr);
  hmeTwoFermions2GammaZ2TwoFermions
    .initPointers(particleDataPtr, couplingsPtr, settingsPtr);
  hmeW2TwoFermions
    .initPointers(particleDataPtr, couplingsPtr, settingsPtr);
  hmeZ2TwoFermions
    .initPointers(particleDataPtr, couplingsPtr, settingsPtr);
  hmeGamma2TwoFermions
    .initPointers(particleDataPtr, couplingsPtr);
  hmeHiggs2TwoFermions
    .initPointers(particleDataPtr, couplingsPtr, settingsPtr);

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

  // User selected tau settings.
  tauExt    = settingsPtr->mode("TauDecays:externalMode");
  tauMode   = settingsPtr->mode("TauDecays:mode");
  tauMother = settingsPtr->mode("TauDecays:tauMother");
  tauPol    = settingsPtr->parm("TauDecays:tauPolarization");

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

  // Set the outgoing particles of the hard process.
  out1                    = HelicityParticle(event[idxOut1]);
  int         idxOut1Top  = out1.iTopCopyId();
  vector<int> sistersOut1 = event[idxOut1Top].sisterList();
  int         idxOut2Top  = idxOut1Top;
  if (sistersOut1.size() == 1) idxOut2Top = sistersOut1[0];
  else {
    // If more then one sister, select by preference tau, nu_tau, lep, nu_lep.
    int tau(-1), tnu(-1), lep(-1), lnu(-1);
    for (int i = 0; i < int(sistersOut1.size()); ++i) {
      int sn = out1.id() == 15 ? -1 : 1;
      int id = event[sistersOut1[i]].id();
      if      (id == sn * 15 && tau == -1) tau = sistersOut1[i];
      else if (id == sn * 16 && tnu == -1) tnu = sistersOut1[i];
      else if ((id == sn * 11 || (id == sn * 13)) && lep == -1)
        lep = sistersOut1[i];
      else if ((id == sn * 12 || (id == sn * 14)) && lnu == -1)
        lnu = sistersOut1[i];
    }
    if      (tau > 0) idxOut2Top = tau;
    else if (tnu > 0) idxOut2Top = tnu;
    else if (lep > 0) idxOut2Top = lep;
    else if (lnu > 0) idxOut2Top = lnu;
  }
  int idxOut2 = event[idxOut2Top].iBotCopyId();
  out2        = HelicityParticle(event[idxOut2]);

  // Set the mediator of the hard process.
  int idxMediator    = event[idxOut1Top].mother1();
  mediator           = HelicityParticle(event[idxMediator]);
  mediator.direction = -1;
  if (mediator.m() < out1.m() + out2.m()) {
    Vec4 p = out1.p() + out2.p();
    mediator.p(p);
    mediator.m(p.mCalc());
  }

  // Set the incoming particles of the hard process.
  int idxMediatorTop = mediator.iTopCopyId();
  int idxIn1         = event[idxMediatorTop].mother1();
  int idxIn2         = event[idxMediatorTop].mother2();
  in1                = HelicityParticle(event[idxIn1]);
  in1.direction      = -1;
  in2                = HelicityParticle(event[idxIn2]);
  in2.direction      = -1;

  // Set the particles vector.
  particles.clear();
  particles.push_back(in1);
  particles.push_back(in2);
  particles.push_back(out1);
  particles.push_back(out2);

  // Determine if correlated (allow lepton flavor violating partner).
  correlated = false;
  if (idxOut1 != idxOut2 &&
      (abs(out2.id()) == 11 || abs(out2.id()) == 13 || abs(out2.id()) == 15)) {
    correlated = true;
    // Check partner vertex is within limits.
    if (limitTau0 && out2.tau0() > tau0Max) correlated = false;
    else if (limitTau && out2.tau() > tauMax) correlated = false;
    else if (limitRadius && pow2(out2.xDec()) + pow2(out2.yDec())
             + pow2(out2.zDec()) > pow2(rMax)) correlated = false;
    else if (limitCylinder && (pow2(out2.xDec()) + pow2(out2.yDec())
                               > pow2(xyMax) || abs(out2.zDec()) > zMax))
      correlated = false;
    // Check partner can decay.
    else if (!out2.canDecay()) correlated = false;
    else if (!out2.mayDecay()) correlated = false;
  }

  // Set the production mechanism.
  bool known = false;
  hardME = 0;
  if      (tauMode == 4) known = internalMechanism(event);
  else if (tauMode == 5) known = externalMechanism(event);
  else {
    if ((tauMode == 2 && abs(mediator.id()) == tauMother) || tauMode == 3) {
      known       = true;
      correlated  = false;
      double sign = out1.id() == -15 ? -1 : 1;
      particles[2].rho[0][0] = (1 - sign * tauPol) / 2;
      particles[2].rho[1][1] = (1 + sign * tauPol) / 2;
    } else {
      if (!externalMechanism(event)) known = internalMechanism(event);
      else known = true;
    }
  }

  // Catch unknown production mechanims.
  if (!known) {
    particles[1] = mediator;
    if (abs(mediator.id()) == 22)
      hardME = hmeGamma2TwoFermions.initChannel(particles);
    else if (abs(mediator.id()) == 23 || abs(mediator.id()) == 32)
      hardME = hmeZ2TwoFermions.initChannel(particles);
    else if (abs(mediator.id()) == 24 || abs(mediator.id()) == 34)
      hardME = hmeW2TwoFermions.initChannel(particles);
    else if (correlated) {
      Vec4 p = out1.p() + out2.p();
      particles[1] = HelicityParticle(22, -22, idxIn1, idxIn2, idxOut1,
        idxOut2, 0, 0, p, p.mCalc(), 0, particleDataPtr);
      hardME = hmeGamma2TwoFermions.initChannel(particles);
      infoPtr->errorMsg("Warning in TauDecays::decay: unknown correlated "
                        "tau production, assuming from unpolarized photon");
    } else {
      infoPtr->errorMsg("Warning in TauDecays::decay: unknown uncorrelated "
                        "tau production, assuming unpolarized");
    }
  }

  // Undecay correlated partner if already decayed.
  if (correlated && !out2.isFinal()) event[out2.index()].undoDecay();

  // Pick the first tau to decay.
  HelicityParticle* tau;
  int idx = 2;
  if (correlated) idx = (rndmPtr->flat() < 0.5) ? 2 : 3;
  tau = &particles[idx];

  // Calculate the density matrix (if needed) and select channel.
  if (hardME) hardME->calculateRho(idx, particles);
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
    tau->D = children[0].D;
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

// Determine the tau polarization and tau decay correlation using the internal
// helicity matrix elements.

bool TauDecays::internalMechanism(Event&) {

  // Flag if process is known.
  bool known = true;

  // Produced from a photon, Z, or Z'.
  if (abs(mediator.id()) == 22 || abs(mediator.id()) == 23 ||
      abs(mediator.id()) == 32) {
    // Produced from fermions: s-channel.
    if (abs(in1.id()) <= 18 && abs(in2.id()) <= 18 && in1.daughter2() == 0 &&
        in2.daughter2() == 0 && in1.daughter1() == in2.daughter1()) {
      particles.push_back(mediator);
      hardME = hmeTwoFermions2GammaZ2TwoFermions.initChannel(particles);
    // Unknown photon production.
    } else known = false;

  // Produced from a W or W'.
  } else if (abs(mediator.id()) == 24 || abs(mediator.id()) == 34) {
    // Produced from fermions: s-channel.
    if (abs(in1.id()) <= 18 && abs(in2.id()) <= 18 && in1.daughter2() == 0 &&
        in2.daughter2() == 0 && in1.daughter1() == in2.daughter1()) {
      particles.push_back(mediator);
      hardME = hmeTwoFermions2W2TwoFermions.initChannel(particles);
    // Unknown W production.
    } else known = false;

  // Produced from a Higgs.
  } else if (abs(mediator.id()) == 25 || abs(mediator.id()) == 35 ||
             abs(mediator.id()) == 36 || abs(mediator.id()) == 37) {
    particles[1] = mediator;
    hardME = hmeHiggs2TwoFermions.initChannel(particles);

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
    particles[0].index(-1);
    particles[1].index(-1);

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

  // Unknown production.
  } else known = false;
  return known;

}

//--------------------------------------------------------------------------

// Determine the tau polarization and tau decay correlation using the provided
// SPINUP digits.

bool TauDecays::externalMechanism(Event &event) {

  // Flag if process is known.
  bool known = true;

  // Uncorrelated, take directly from SPINUP if valid.
  if (tauExt == 0) correlated = false;
  if (!correlated) {
    double spinup = particles[2].pol();
    if (abs(spinup) > 1.001) spinup = event[particles[2].iTopCopyId()].pol();
    if (abs(spinup) > 1.001) known = false;
    else {
      particles[2].rho[0][0] = (1 - spinup) / 2;
      particles[2].rho[1][1] = (1 + spinup) / 2;
    }

  // Correlated, try mother.
  } else if (tauExt == 1) {
    double spinup = mediator.pol();
    if (abs(spinup) > 1.001) spinup = event[mediator.iTopCopyId()].pol();
    if (abs(spinup) > 1.001) spinup = 0;
    if (mediator.rho.size() > 1) {
      mediator.rho[0][0] = (1 - spinup) / mediator.spinStates();
      mediator.rho[1][1] = (1 + spinup) / mediator.spinStates();
    }
    particles[1] = mediator;
    if (abs(mediator.id()) == 22)
      hardME = hmeGamma2TwoFermions.initChannel(particles);
    else if (abs(mediator.id()) == 23 || abs(mediator.id()) == 32)
      hardME = hmeZ2TwoFermions.initChannel(particles);
    else if (abs(mediator.id()) == 24 || abs(mediator.id()) == 34)
      hardME = hmeZ2TwoFermions.initChannel(particles);
    else if (abs(mediator.id()) == 25 || abs(mediator.id()) == 35 ||
             abs(mediator.id()) == 36 || abs(mediator.id()) == 37)
      hardME = hmeHiggs2TwoFermions.initChannel(particles);
    else known = false;

    // Unknown mechanism.
  } else known = false;
  return known;

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
        double childMass = particleDataPtr->mSel(childId);
        // Push back the child into the children vector.
        children.push_back( HelicityParticle(childId, 91, parent.index(), 0, 0,
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
    children[i].index(event.append(children[i]));
  }

  // Mark the parent as decayed and set children.
  event[children[0].index()].statusNeg();
  event[children[0].index()].daughters(children[1].index(),
    children[decayMult].index());

}

//==========================================================================

} // end namespace Pythia8
