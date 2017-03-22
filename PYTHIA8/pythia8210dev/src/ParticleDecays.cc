// ParticleDecays.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// ParticleDecays class.

#include "Pythia8/ParticleDecays.h"

namespace Pythia8 {

//==========================================================================

// The ParticleDecays class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of times one tries to let decay happen (for 2 nested loops).
const int    ParticleDecays::NTRYDECAY   = 10;

// Number of times one tries to pick valid hadronic content in decay.
const int    ParticleDecays::NTRYPICK    = 100;

// Number of times one tries to pick decay topology.
const int    ParticleDecays::NTRYMEWT    = 1000;

// Maximal loop count in Dalitz decay treatment.
const int    ParticleDecays::NTRYDALITZ  = 1000;

// Minimal Dalitz pair mass this factor above threshold.
const double ParticleDecays::MSAFEDALITZ = 1.000001;

// These numbers are hardwired empirical parameters,
// intended to speed up the M-generator.
const double ParticleDecays::WTCORRECTION[11] = { 1., 1., 1.,
  2., 5., 15., 60., 250., 1250., 7000., 50000. };

//--------------------------------------------------------------------------

// Initialize and save pointers.

void ParticleDecays::init(Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
  Couplings* couplingsPtrIn, TimeShower* timesDecPtrIn,
  StringFlav* flavSelPtrIn, DecayHandler* decayHandlePtrIn,
  vector<int> handledParticles) {

  // Save pointers to error messages handling and flavour generation.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  couplingsPtr    = couplingsPtrIn;
  flavSelPtr      = flavSelPtrIn;

  // Save pointer to timelike shower, as needed in some few decays.
  timesDecPtr     = timesDecPtrIn;

  // Save pointer for external handling of some decays.
  decayHandlePtr  = decayHandlePtrIn;

  // Set which particles should be handled externally.
  if (decayHandlePtr != 0)
  for (int i = 0; i < int(handledParticles.size()); ++i)
    particleDataPtr->doExternalDecay(handledParticles[i], true);

  // Safety margin in mass to avoid troubles.
  mSafety       = settings.parm("ParticleDecays:mSafety");

  // Lifetime and vertex rules for determining whether decay allowed.
  limitTau0     = settings.flag("ParticleDecays:limitTau0");
  tau0Max       = settings.parm("ParticleDecays:tau0Max");
  limitTau      = settings.flag("ParticleDecays:limitTau");
  tauMax        = settings.parm("ParticleDecays:tauMax");
  limitRadius   = settings.flag("ParticleDecays:limitRadius");
  rMax          = settings.parm("ParticleDecays:rMax");
  limitCylinder = settings.flag("ParticleDecays:limitCylinder");
  xyMax         = settings.parm("ParticleDecays:xyMax");
  zMax          = settings.parm("ParticleDecays:zMax");
  limitDecay    = limitTau0 || limitTau || limitRadius || limitCylinder;

  // B-Bbar mixing parameters.
  mixB          = settings.flag("ParticleDecays:mixB");
  xBdMix        = settings.parm("ParticleDecays:xBdMix");
  xBsMix        = settings.parm("ParticleDecays:xBsMix");

  // Suppression of extra-hadron momenta in semileptonic decays.
  sigmaSoft     = settings.parm("ParticleDecays:sigmaSoft");

  // Selection of multiplicity and colours in "phase space" model.
  multIncrease     = settings.parm("ParticleDecays:multIncrease");
  multIncreaseWeak = settings.parm("ParticleDecays:multIncreaseWeak");
  multRefMass      = settings.parm("ParticleDecays:multRefMass");
  multGoffset      = settings.parm("ParticleDecays:multGoffset");
  colRearrange     = settings.parm("ParticleDecays:colRearrange");

  // Minimum energy in system (+ m_q) from StringFragmentation.
  stopMass      = settings.parm("StringFragmentation:stopMass");

  // Parameters for Dalitz decay virtual gamma mass spectrum.
  sRhoDal       = pow2(particleDataPtr->m0(113));
  wRhoDal       = pow2(particleDataPtr->mWidth(113));

  // Allow showers in decays to qqbar/gg/ggg/gammagg.
  doFSRinDecays = settings.flag("ParticleDecays:FSRinDecays");
  doGammaRad    = settings.flag("ParticleDecays:allowPhotonRadiation");

  // Use standard decays or dedicated tau decay package
  tauMode       = settings.mode("TauDecays:mode");

  // Initialize the dedicated tau decay handler.
  if (tauMode) tauDecayer.init(infoPtr, &settings,
    particleDataPtr, rndmPtr, couplingsPtr);

}

//--------------------------------------------------------------------------

// Decay a particle; main method.

bool ParticleDecays::decay( int iDec, Event& event) {

  // Check whether a decay is allowed, given the upcoming decay vertex.
  Particle& decayer = event[iDec];
  hasPartons  = false;
  keepPartons = false;
  if (limitDecay && !checkVertex(decayer)) return true;

  // Do not allow resonance decays (beyond handling capability).
  if (decayer.isResonance()) {
    infoPtr->errorMsg("Warning in ParticleDecays::decay: "
      "resonance left undecayed");
    return true;
  }

  // Fill the decaying particle in slot 0 of arrays.
  idDec = decayer.id();
  iProd.resize(0);
  idProd.resize(0);
  mProd.resize(0);
  iProd.push_back( iDec );
  idProd.push_back( idDec );
  mProd.push_back( decayer.m() );

  // Check for oscillations B0 <-> B0bar or B_s0 <-> B_s0bar.
  bool hasOscillated = (abs(idDec) == 511 || abs(idDec) == 531)
    ? oscillateB(decayer) : false;
  if (hasOscillated) {idDec = - idDec; idProd[0] = idDec;}

  // Particle data for decaying particle.
  decDataPtr = &decayer.particleDataEntry();

  // Optionally send on to external decay program.
  bool doneExternally = false;
  if (decDataPtr->doExternalDecay()) {
    pProd.resize(0);
    pProd.push_back(decayer.p());
    doneExternally = decayHandlePtr->decay(idProd, mProd, pProd,
      iDec, event);

    // If it worked, then store the decay products in the event record.
    if (doneExternally) {
      mult = idProd.size() - 1;
      int status = (hasOscillated) ? 94 : 93;
      for (int i = 1; i <= mult; ++i) {
        int iPos = event.append( idProd[i], status, iDec, 0, 0, 0,
        0, 0, pProd[i], mProd[i]);
        iProd.push_back( iPos);
      }

      // Also mark mother decayed and store daughters.
      event[iDec].statusNeg();
      event[iDec].daughters( iProd[1], iProd[mult]);
    }
  }

  // Check if the particle is tau and let the special tau decayer handle it.
  if (decayer.idAbs() == 15 && !doneExternally && tauMode) {
    doneExternally = tauDecayer.decay(iDec, event);
    if (doneExternally) return true;
  }

  // Now begin normal internal decay treatment.
  if (!doneExternally) {

    // Allow up to ten tries to pick a channel.
    if (!decDataPtr->preparePick(idDec, decayer.m())) return false;
    bool foundChannel = false;
    bool hasStored    = false;
    for (int iTryChannel = 0; iTryChannel < NTRYDECAY; ++iTryChannel) {

      // Remove previous failed channel.
      if (hasStored) event.popBack(mult);
      hasStored = false;

      // Pick new channel. Read out basics.
      DecayChannel& channel = decDataPtr->pickChannel();
      meMode = channel.meMode();
      keepPartons = (meMode > 90 && meMode <= 100);
      mult = channel.multiplicity();

      // Allow up to ten tries for each channel (e.g with different masses).
      bool foundMode = false;
      for (int iTryMode = 0; iTryMode < NTRYDECAY; ++iTryMode) {
        idProd.resize(1);
        mProd.resize(1);
        scale = 0.;

        // Extract and store the decay products in local arrays.
        hasPartons = false;
        for (int i = 0; i < mult; ++i) {
          int idNow = channel.product(i);
          int idAbs = abs(idNow);
          if ( idAbs < 10 || idAbs == 21 || idAbs == 81 || idAbs == 82
            || idAbs == 83 || (idAbs > 1000 && idAbs < 10000
            && (idAbs/10)%10 == 0) ) hasPartons = true;
          if (idDec < 0 && particleDataPtr->hasAnti(idNow)) idNow = -idNow;
          double mNow = particleDataPtr->mSel(idNow);
          idProd.push_back( idNow);
          mProd.push_back( mNow);
        }

        // Decays into partons usually translate into hadrons.
        if (hasPartons && !keepPartons && !pickHadrons()) continue;

        // Need to set colour flow if explicit decay to partons.
        cols.resize(0);
        acols.resize(0);
        for (int i = 0; i <= mult; ++i) {
          cols.push_back(0);
          acols.push_back(0);
        }
        if (hasPartons && keepPartons && !setColours(event)) continue;

        // Check that enough phase space for decay.
        if (mult > 1) {
          double mDiff = mProd[0];
          for (int i = 1; i <= mult; ++i) mDiff -= mProd[i];
          if (mDiff < mSafety) continue;
        }

        // End of inner trial loops. Check if succeeded or not.
        foundMode = true;
        break;
      }
      if (!foundMode) continue;

      // Store decay products in the event record.
      int status = (hasOscillated) ? 92 : 91;
      for (int i = 1; i <= mult; ++i) {
        int iPos = event.append( idProd[i], status, iDec, 0, 0, 0,
          cols[i], acols[i], Vec4(0., 0., 0., 0.), mProd[i], scale);
        iProd.push_back( iPos);
      }
      hasStored = true;

      // Pick mass of Dalitz decay. Temporarily change multiplicity.
      if ( (meMode == 11 || meMode == 12 || meMode == 13)
        && !dalitzMass() ) continue;

      // Do a decay, split by multiplicity.
      bool decayed = false;
      if      (mult == 1) decayed = oneBody(event);
      else if (mult == 2) decayed = twoBody(event);
      else if (mult == 3) decayed = threeBody(event);
      else                decayed = mGenerator(event);
      if (!decayed) continue;

      // Kinematics of gamma* -> l- l+ in Dalitz decay. Restore multiplicity.
      if (meMode == 11 || meMode == 12 || meMode == 13)
        dalitzKinematics(event);

      // End of outer trial loops.
      foundChannel = true;
      break;
    }

    // If the decay worked, then mark mother decayed and store daughters.
    if (foundChannel) {
      event[iDec].statusNeg();
      event[iDec].daughters( iProd[1], iProd[mult]);

    // Else remove unused daughters and return failure.
    } else {
      if (hasStored) event.popBack(mult);
      infoPtr->errorMsg("Error in ParticleDecays::decay: "
        "failed to find workable decay channel");
      return false;
    }

  // Now finished normal internal decay treatment.
  }

  // Set decay vertex when this is displaced.
  if (event[iDec].hasVertex() || event[iDec].tau() > 0.) {
    Vec4 vDec = event[iDec].vDec();
    for (int i = 1; i <= mult; ++i) event[iProd[i]].vProd( vDec );
  }

  // Set lifetime of daughters.
  for (int i = 1; i <= mult; ++i)
    event[iProd[i]].tau( event[iProd[i]].tau0() * rndmPtr->exp() );

  // In a decay explicitly to partons then optionally do a shower,
  // and always flag that partonic system should be fragmented.
  if (hasPartons && keepPartons && doFSRinDecays)
    timesDecPtr->shower( iProd[1], iProd.back(), event, mProd[0]);

  // Photon radiation implemented only for two-body decay to leptons.
  else if (doGammaRad && mult == 2 && event[iProd[1]].isLepton()
  && event[iProd[2]].isLepton())
    timesDecPtr->showerQED( iProd[1], iProd[2], event, mProd[0]);

  // For Hidden Valley particles also allow leptons to shower.
  else if (event[iDec].idAbs() > 4900000 && event[iDec].idAbs() < 5000000
  && doFSRinDecays && mult == 2 && event[iProd[1]].isLepton()) {
    event[iProd[1]].scale(mProd[0]);
    event[iProd[2]].scale(mProd[0]);
    timesDecPtr->shower( iProd[1], iProd.back(), event, mProd[0]);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Check whether a decay is allowed, given the upcoming decay vertex.

bool ParticleDecays::checkVertex(Particle& decayer) {

  // Check whether any of the conditions are not fulfilled.
  if (limitTau0 && decayer.tau0() > tau0Max) return false;
  if (limitTau && decayer.tau() > tauMax) return false;
  if (limitRadius && pow2(decayer.xDec()) + pow2(decayer.yDec())
    + pow2(decayer.zDec()) > pow2(rMax)) return false;
  if (limitCylinder && (pow2(decayer.xDec()) + pow2(decayer.yDec())
    > pow2(xyMax) || abs(decayer.zDec()) > zMax) ) return false;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Check for oscillations B0 <-> B0bar or B_s0 <-> B_s0bar.

bool ParticleDecays::oscillateB(Particle& decayer) {

  // Extract relevant information and decide.
  if (!mixB) return false;
  double xBmix   = (abs(decayer.id()) == 511) ? xBdMix : xBsMix;
  double tau     = decayer.tau();
  double tau0    = decayer.tau0();
  double probosc = pow2(sin(0.5 * xBmix * tau / tau0));
  return (probosc > rndmPtr->flat());

}

//--------------------------------------------------------------------------

// Do a one-body decay. (Rare; e.g. for K0 -> K0_short.)

bool ParticleDecays::oneBody(Event& event) {

  // References to the particles involved.
  Particle& decayer = event[iProd[0]];
  Particle& prod    = event[iProd[1]];

  // Set momentum and expand mother information.
  prod.p( decayer.p() );
  prod.m( decayer.m() );
  prod.mother2( iProd[0] );

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do a two-body decay.

bool ParticleDecays::twoBody(Event& event) {

  // References to the particles involved.
  Particle& decayer = event[iProd[0]];
  Particle& prod1   = event[iProd[1]];
  Particle& prod2   = event[iProd[2]];

  // Masses.
  double m0   = mProd[0];
  double m1   = mProd[1];
  double m2   = mProd[2];

  // Energies and absolute momentum in the rest frame.
  if (m1 + m2 + mSafety > m0) return false;
  double e1   = 0.5 * (m0*m0 + m1*m1 - m2*m2) / m0;
  double e2   = 0.5 * (m0*m0 + m2*m2 - m1*m1) / m0;
  double pAbs = 0.5 * sqrtpos( (m0 - m1 - m2) * (m0 + m1 + m2)
    * (m0 + m1 - m2) * (m0 - m1 + m2) ) / m0;

  // When meMode = 2, for V -> PS2 + PS3 (V = vector, pseudoscalar),
  // need to check if production is PS0 -> PS1/gamma + V.
  int iMother = event[iProd[0]].mother1();
  int idSister = 0;
  if (meMode == 2) {
    if (iMother <= 0 || iMother >= iProd[0]) meMode = 0;
    else {
      int iDaughter1 = event[iMother].daughter1();
      int iDaughter2 = event[iMother].daughter2();
      if (iDaughter2 != iDaughter1 + 1) meMode = 0;
      else {
        int idMother = abs( event[iMother].id() );
        if (idMother <= 100 || idMother%10 !=1
          || (idMother/1000)%10 != 0) meMode = 0;
        else {
          int iSister = (iProd[0] == iDaughter1) ? iDaughter2 : iDaughter1;
          idSister = abs( event[iSister].id() );
          if ( (idSister <= 100 || idSister%10 !=1
            || (idSister/1000)%10 != 0) && idSister != 22) meMode = 0;
        }
      }
    }
  }

  // Begin loop over matrix-element corrections.
  double wtME, wtMEmax;
  int loop = 0;
  do {
    wtME = 1.;
    wtMEmax = 1.;
    ++loop;

    // Isotropic angles give three-momentum.
    double cosTheta = 2. * rndmPtr->flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi      = 2. * M_PI * rndmPtr->flat();
    double pX       = pAbs * sinTheta * cos(phi);
    double pY       = pAbs * sinTheta * sin(phi);
    double pZ       = pAbs * cosTheta;

    // Fill four-momenta and boost them away from mother rest frame.
    prod1.p(  pX,  pY,  pZ, e1);
    prod2.p( -pX, -pY, -pZ, e2);
    prod1.bst( decayer.p(), decayer.m() );
    prod2.bst( decayer.p(), decayer.m() );

    // Matrix element for PS0 -> PS1 + V1 -> PS1 + PS2 + PS3 of form
    // cos**2(theta02) in V1 rest frame, and for PS0 -> gamma + V1
    // -> gamma + PS2 + PS3 of form sin**2(theta02).
    if (meMode == 2) {
      double p10 = decayer.p() * event[iMother].p();
      double p12 = decayer.p() * prod1.p();
      double p02 = event[iMother].p() * prod1.p();
      double s0  = pow2(event[iMother].m());
      double s1  = pow2(decayer.m());
      double s2  =  pow2(prod1.m());
      if (idSister != 22) wtME = pow2(p10 * p12 - s1 * p02);
      else wtME = s1 * (2. * p10 * p12 * p02 - s1 * p02*p02
        - s0 * p12*p12 - s2 * p10*p10 + s1 * s0 * s2);
      wtME = max( wtME, 1e-6 * s1*s1 * s0 * s2);
      wtMEmax = (p10*p10 - s1 * s0) * (p12*p12 - s1 * s2);
    }

    // Break out of loop if no sensible ME weight.
    if(loop > NTRYMEWT) {
      infoPtr->errorMsg("ParticleDecays::twoBody: "
        "caught in infinite ME weight loop");
      wtME = abs(wtMEmax);
    }

  // If rejected, try again with new invariant masses.
  } while ( wtME < rndmPtr->flat() * wtMEmax );

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do a three-body decay (except Dalitz decays).

bool ParticleDecays::threeBody(Event& event) {

  // References to the particles involved.
  Particle& decayer = event[iProd[0]];
  Particle& prod1   = event[iProd[1]];
  Particle& prod2   = event[iProd[2]];
  Particle& prod3   = event[iProd[3]];

  // Mother and sum daughter masses. Fail if too close.
  double m0      = mProd[0];
  double m1      = mProd[1];
  double m2      = mProd[2];
  double m3      = mProd[3];
  double mSum    = m1 + m2 + m3;
  double mDiff   = m0 - mSum;
  if (mDiff < mSafety) return false;

  // Kinematical limits for 2+3 mass. Maximum phase-space weight.
  double m23Min  = m2 + m3;
  double m23Max  = m0 - m1;
  double p1Max   = 0.5 * sqrtpos( (m0 - m1 - m23Min) * (m0 + m1 + m23Min)
    * (m0 + m1 - m23Min) * (m0 - m1 + m23Min) ) / m0;
  double p23Max  = 0.5 * sqrtpos( (m23Max - m2 - m3) * (m23Max + m2 + m3)
    * (m23Max + m2 - m3) * (m23Max - m2 + m3) ) / m23Max;
  double wtPSmax = 0.5 * p1Max * p23Max;

  // Begin loop over matrix-element corrections.
  double wtME, wtMEmax, wtPS, m23, p1Abs, p23Abs;
  do {
    wtME     = 1.;
    wtMEmax  = 1.;

    // Pick an intermediate mass m23 flat in the allowed range.
    do {
      m23    = m23Min + rndmPtr->flat() * mDiff;

      // Translate into relative momenta and find phase-space weight.
      p1Abs  = 0.5 * sqrtpos( (m0 - m1 - m23) * (m0 + m1 + m23)
        * (m0 + m1 - m23) * (m0 - m1 + m23) ) / m0;
      p23Abs = 0.5 * sqrtpos( (m23 - m2 - m3) * (m23 + m2 + m3)
        * (m23 + m2 - m3) * (m23 - m2 + m3) ) / m23;
      wtPS   = p1Abs * p23Abs;

    // If rejected, try again with new invariant masses.
    } while ( wtPS < rndmPtr->flat() * wtPSmax );

    // Set up m23 -> m2 + m3 isotropic in its rest frame.
    double cosTheta = 2. * rndmPtr->flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi      = 2. * M_PI * rndmPtr->flat();
    double pX       = p23Abs * sinTheta * cos(phi);
    double pY       = p23Abs * sinTheta * sin(phi);
    double pZ       = p23Abs * cosTheta;
    double e2       = sqrt( m2*m2 + p23Abs*p23Abs);
    double e3       = sqrt( m3*m3 + p23Abs*p23Abs);
    prod2.p(  pX,  pY,  pZ, e2);
    prod3.p( -pX, -pY, -pZ, e3);

    // Set up m0 -> m1 + m23 isotropic in its rest frame.
    cosTheta        = 2. * rndmPtr->flat() - 1.;
    sinTheta        = sqrt(1. - cosTheta*cosTheta);
    phi             = 2. * M_PI * rndmPtr->flat();
    pX              = p1Abs * sinTheta * cos(phi);
    pY              = p1Abs * sinTheta * sin(phi);
    pZ              = p1Abs * cosTheta;
    double e1       = sqrt( m1*m1 + p1Abs*p1Abs);
    double e23      = sqrt( m23*m23 + p1Abs*p1Abs);
    prod1.p( pX, pY, pZ, e1);

    // Boost 2 + 3 to the 0 rest frame.
    Vec4 p23( -pX, -pY, -pZ, e23);
    prod2.bst( p23, m23 );
    prod3.bst( p23, m23 );

    // Matrix-element weight for omega/phi -> pi+ pi- pi0.
    if (meMode == 1) {
      double p1p2 = prod1.p() * prod2.p();
      double p1p3 = prod1.p() * prod3.p();
      double p2p3 = prod2.p() * prod3.p();
      wtME = pow2(m1 * m2 * m3) - pow2(m1 * p2p3) - pow2(m2 * p1p3)
        - pow2(m3 * p1p2) + 2. * p1p2 * p1p3 * p2p3;
      wtMEmax = pow3(m0 * m0) / 150.;

    // Effective matrix element for nu spectrum in tau -> nu + hadrons.
    } else if (meMode == 21) {
      double x1 = 2. *  prod1.e() / m0;
      wtME = x1 * (3. - 2. * x1);
      double xMax = min( 0.75, 2. * (1. - mSum / m0) );
      wtMEmax = xMax * (3. - 2. * xMax);

    // Matrix element for weak decay (only semileptonic for c and b).
    } else if ((meMode == 22 || meMode == 23) && prod1.isLepton()) {
      wtME = m0 * prod1.e() * (prod2.p() * prod3.p());
      wtMEmax = min( pow4(m0) / 16., m0 * (m0 - m1 - m2) * (m0 - m1 - m3)
        * (m0 - m2 - m3) );

    // Effective matrix element for weak decay to hadrons (B -> D, D -> K).
    } else if (meMode == 22 || meMode == 23) {
      double x1 = 2. * prod1.pAbs() / m0;
      wtME = x1 * (3. - 2. * x1);
      double xMax = min( 0.75, 2. * (1. - mSum / m0) );
      wtMEmax = xMax * (3. - 2. * xMax);

    // Effective matrix element for gamma spectrum in B -> gamma + hadrons.
    } else if (meMode == 31) {
      double x1 = 2. * prod1.e() / m0;
      wtME = pow3(x1);
      double x1Max = 1. - pow2(mSum / m0);
      wtMEmax = pow3(x1Max);

    // Matrix-element weight for "onium" -> g + g + g or gamma + g + g.
    } else if (meMode == 92) {
      double x1 = 2. * prod1.e() / m0;
      double x2 = 2. * prod2.e() / m0;
      double x3 = 2. * prod3.e() / m0;
      wtME = pow2( (1. - x1) / (x2 * x3) ) + pow2( (1. - x2) / (x1 * x3) )
        + pow2( (1. - x3) / (x1 * x2) );
      wtMEmax = 2.;
      // For gamma + g + g require minimum mass for g + g system.
      if (prod1.id() == 22 && sqrt(1. - x1) * m0 < 2. * stopMass) wtME = 0.;
      if (prod2.id() == 22 && sqrt(1. - x2) * m0 < 2. * stopMass) wtME = 0.;
      if (prod3.id() == 22 && sqrt(1. - x3) * m0 < 2. * stopMass) wtME = 0.;
    }

  // If rejected, try again with new invariant masses.
  } while ( wtME < rndmPtr->flat() * wtMEmax );

  // Boost 1 + 2 + 3 to the current frame.
  prod1.bst( decayer.p(), decayer.m() );
  prod2.bst( decayer.p(), decayer.m() );
  prod3.bst( decayer.p(), decayer.m() );

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do a multibody decay using the M-generator algorithm.

bool ParticleDecays::mGenerator(Event& event) {

  // Mother and sum daughter masses. Fail if too close or inconsistent.
  double m0      = mProd[0];
  double mSum    = mProd[1];
  for (int i = 2; i <= mult; ++i) mSum += mProd[i];
  double mDiff   = m0 - mSum;
  if (mDiff < mSafety) return false;

  // Begin setup of intermediate invariant masses.
  mInv.resize(0);
  for (int i = 0; i <= mult; ++i) mInv.push_back( mProd[i]);

  // Calculate the maximum weight in the decay.
  double wtPS, wtME, wtMEmax;
  double wtPSmax = 1. / WTCORRECTION[mult];
  double mMax    = mDiff + mProd[mult];
  double mMin    = 0.;
  for (int i = mult - 1; i > 0; --i) {
    mMax        += mProd[i];
    mMin        += mProd[i+1];
    double mNow  = mProd[i];
    wtPSmax     *= 0.5 * sqrtpos( (mMax - mMin - mNow) * (mMax + mMin + mNow)
                 * (mMax + mMin - mNow) * (mMax - mMin + mNow) ) / mMax;
  }

  // Begin loop over matrix-element corrections.
  do {
    wtME    = 1.;
    wtMEmax = 1.;

    // Begin loop to find the set of intermediate invariant masses.
    do {
      wtPS  = 1.;

      // Find and order random numbers in descending order.
      rndmOrd.resize(0);
      rndmOrd.push_back(1.);
      for (int i = 1; i < mult - 1; ++i) {
        double rndm = rndmPtr->flat();
        rndmOrd.push_back(rndm);
        for (int j = i - 1; j > 0; --j) {
          if (rndm > rndmOrd[j]) swap( rndmOrd[j], rndmOrd[j+1] );
          else break;
        }
      }
      rndmOrd.push_back(0.);

      // Translate into intermediate masses and find weight.
      for (int i = mult - 1; i > 0; --i) {
        mInv[i] = mInv[i+1] + mProd[i] + (rndmOrd[i-1] - rndmOrd[i]) * mDiff;
        wtPS   *= 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - mProd[i])
          * (mInv[i] + mInv[i+1] + mProd[i]) * (mInv[i] + mInv[i+1] - mProd[i])
          * (mInv[i] - mInv[i+1] + mProd[i]) ) / mInv[i];
      }

    // If rejected, try again with new invariant masses.
    } while ( wtPS < rndmPtr->flat() * wtPSmax );

    // Perform two-particle decays in the respective rest frame.
    pInv.resize(mult + 1);
    for (int i = 1; i < mult; ++i) {
      double pAbs = 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - mProd[i])
        * (mInv[i] + mInv[i+1] + mProd[i]) * (mInv[i] + mInv[i+1] - mProd[i])
        * (mInv[i] - mInv[i+1] + mProd[i]) ) / mInv[i];

      // Isotropic angles give three-momentum.
      double cosTheta = 2. * rndmPtr->flat() - 1.;
      double sinTheta = sqrt(1. - cosTheta*cosTheta);
      double phi      = 2. * M_PI * rndmPtr->flat();
      double pX       = pAbs * sinTheta * cos(phi);
      double pY       = pAbs * sinTheta * sin(phi);
      double pZ       = pAbs * cosTheta;

      // Calculate energies, fill four-momenta.
      double eHad     = sqrt( mProd[i]*mProd[i] + pAbs*pAbs);
      double eInv     = sqrt( mInv[i+1]*mInv[i+1] + pAbs*pAbs);
      event[iProd[i]].p( pX, pY, pZ, eHad);
      pInv[i+1].p( -pX, -pY, -pZ, eInv);
    }

    // Boost decay products to the mother rest frame.
    event[iProd[mult]].p( pInv[mult] );
    for (int iFrame = mult - 1; iFrame > 1; --iFrame)
      for (int i = iFrame; i <= mult; ++i)
        event[iProd[i]].bst( pInv[iFrame], mInv[iFrame]);

    // Effective matrix element for nu spectrum in tau -> nu + hadrons.
    if (meMode == 21 && event[iProd[1]].isLepton()) {
      double x1 = 2. * event[iProd[1]].e() / m0;
      wtME = x1 * (3. - 2. * x1);
      double xMax = min( 0.75, 2. * (1. - mSum / m0) );
      wtMEmax = xMax * (3. - 2. * xMax);

    // Effective matrix element for weak decay (only semileptonic for c and b).
    // Particles 4 onwards should be made softer explicitly?
    } else if ((meMode == 22 || meMode == 23) && event[iProd[1]].isLepton()) {
      Vec4 pRest = event[iProd[3]].p();
      for (int i = 4; i <= mult; ++i) pRest += event[iProd[i]].p();
      wtME = m0 * event[iProd[1]].e() * (event[iProd[2]].p() * pRest);
      for (int i = 4; i <= mult; ++i) wtME
        *= exp(- event[iProd[i]].pAbs2() / pow2(sigmaSoft) );
      wtMEmax = pow4(m0) / 16.;

    // Effective matrix element for weak decay to hadrons (B -> D, D -> K).
    } else if (meMode == 22 || meMode == 23) {
      double x1 = 2. * event[iProd[1]].pAbs() / m0;
      wtME = x1 * (3. - 2. * x1);
      double xMax = min( 0.75, 2. * (1. - mSum / m0) );
      wtMEmax = xMax * (3. - 2. * xMax);

    // Effective matrix element for gamma spectrum in B -> gamma + hadrons.
    } else if (meMode == 31) {
      double x1 = 2. * event[iProd[1]].e() / m0;
      wtME = pow3(x1);
      double x1Max = 1. - pow2(mSum / m0);
      wtMEmax = pow3(x1Max);
    }

  // If rejected, try again with new invariant masses.
  } while ( wtME < rndmPtr->flat() * wtMEmax );

  // Boost decay products to the current frame.
  pInv[1].p( event[iProd[0]].p() );
  for (int i = 1; i <= mult; ++i) event[iProd[i]].bst( pInv[1], mInv[1] );

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select mass of lepton pair in a Dalitz decay.

bool ParticleDecays::dalitzMass() {

  // Mother and sum daughter masses.
  double mSum1 = 0;
  for (int i = 1; i <= mult - 2; ++i) mSum1 += mProd[i];
  if (meMode == 13) mSum1 *= MSAFEDALITZ;
  double mSum2 = MSAFEDALITZ * (mProd[mult -1] + mProd[mult]);
  double mDiff = mProd[0] - mSum1 - mSum2;

  // Fail if too close or inconsistent.
  if (mDiff < mSafety) return false;
  if (idProd[mult - 1] + idProd[mult] != 0
    || mProd[mult - 1] != mProd[mult]) {
    infoPtr->errorMsg("Error in ParticleDecays::dalitzMass:"
    " inconsistent flavour/mass assignments");
    return false;
  }
  if ( meMode == 13 && (idProd[1] + idProd[2] != 0
    || mProd[1] != mProd[2]) ) {
    infoPtr->errorMsg("Error in ParticleDecays::dalitzMass:"
    " inconsistent flavour/mass assignments");
    return false;
  }

  // Case 1: one Dalitz pair.
  if (meMode == 11 || meMode == 12) {

    // Kinematical limits for gamma* squared mass.
    double sGamMin = pow2(mSum2);
    double sGamMax = pow2(mProd[0] - mSum1);
    // Select virtual gamma squared mass. Guessed form for meMode == 12.
    double sGam, wtGam;
    int loop = 0;
    do {
      if (++loop > NTRYDALITZ) return false;
      sGam = sGamMin * pow( sGamMax / sGamMin, rndmPtr->flat() );
      wtGam = (1. + 0.5 * sGamMin / sGam) *  sqrt(1. - sGamMin / sGam)
        * pow3(1. - sGam / sGamMax) * sRhoDal * (sRhoDal + wRhoDal)
        / ( pow2(sGam - sRhoDal) + sRhoDal * wRhoDal );
    } while ( wtGam < rndmPtr->flat() );

    // Store results in preparation for doing a one-less-body decay.
    --mult;
    mProd[mult] = sqrt(sGam);

  // Case 2: two Dalitz pairs.
  } else {

    // Kinematical limits for 1 + 2 and 3 + 4 gamma* masses.
    double s0 = pow2(mProd[0]);
    double s12Min = pow2(mSum1);
    double s12Max = pow2(mProd[0] - mSum2);
    double s34Min = pow2(mSum2);
    double s34Max = pow2(mProd[0] - mSum1);

    // Select virtual gamma squared masses. Guessed form for meMode == 13.
    double s12, s34, wt12, wt34, wtPAbs, wtAll;
    int loop = 0;
    do {
      if (++loop > NTRYDALITZ) return false;
      s12 = s12Min * pow( s12Max / s12Min, rndmPtr->flat() );
      wt12 = (1. + 0.5 * s12Min / s12) *  sqrt(1. - s12Min / s12)
        * sRhoDal * (sRhoDal + wRhoDal)
        / ( pow2(s12 - sRhoDal) + sRhoDal * wRhoDal );
      s34 = s34Min * pow( s34Max / s34Min, rndmPtr->flat() );
      wt34 = (1. + 0.5 * s34Min / s34) *  sqrt(1. - s34Min / s34)
        * sRhoDal * (sRhoDal + wRhoDal)
        / ( pow2(s34 - sRhoDal) + sRhoDal * wRhoDal );
      wtPAbs = sqrtpos( pow2(1. - (s12 + s34)/ s0)
        - 4. * s12 * s34 / (s0 * s0) );
      wtAll = wt12 * wt34 * pow3(wtPAbs);
      if (wtAll > 1.) infoPtr->errorMsg(
        "Error in ParticleDecays::dalitzMass: weight > 1");
    } while (wtAll < rndmPtr->flat());

    // Store results in preparation for doing a two-body decay.
    mult = 2;
    mProd[1] = sqrt(s12);
    mProd[2] = sqrt(s34);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do kinematics of gamma* -> l- l+ in Dalitz decay.

bool ParticleDecays::dalitzKinematics(Event& event) {

  // Restore multiplicity.
  int nDal = (meMode < 13) ? 1 : 2;
  mult += nDal;

  // Loop over one or two lepton pairs.
  for (int iDal = 0; iDal < nDal; ++iDal) {

    // References to the particles involved.
    Particle& decayer = event[iProd[0]];
    Particle& prodA = (iDal == 0) ? event[iProd[mult - 1]]
      : event[iProd[1]];
    Particle& prodB = (iDal == 0) ? event[iProd[mult]]
      : event[iProd[2]];

    // Reconstruct required rotations and boosts backwards.
    Vec4 pDec    = decayer.p();
    int  iGam    = (meMode < 13) ? mult - 1 : 2 - iDal;
    Vec4 pGam    = event[iProd[iGam]].p();
    pGam.bstback( pDec, decayer.m() );
    double phiGam = pGam.phi();
    pGam.rot( 0., -phiGam);
    double thetaGam = pGam.theta();
    pGam.rot( -thetaGam, 0.);

    // Masses and phase space in gamma* rest frame.
    double mGam     = (meMode < 13) ? mProd[mult - 1] : mProd[2 - iDal];
    double mA       = prodA.m();
    double mB       = prodB.m();
    double mGamMin  = MSAFEDALITZ * (mA + mB);
    double mGamRat  = pow2(mGamMin / mGam);
    double pGamAbs  = 0.5 * sqrtpos( (mGam - mA - mB) * (mGam + mA + mB) );

    // Set up decay in gamma* rest frame, reference along +z axis.
    double cosTheta, cos2Theta;
    do {
      cosTheta      = 2. * rndmPtr->flat() - 1.;
      cos2Theta     = cosTheta * cosTheta;
    } while ( 1. + cos2Theta + mGamRat * (1. - cos2Theta)
      < 2. * rndmPtr->flat() );
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi      = 2. * M_PI * rndmPtr->flat();
    double pX       = pGamAbs * sinTheta * cos(phi);
    double pY       = pGamAbs * sinTheta * sin(phi);
    double pZ       = pGamAbs * cosTheta;
    double eA       = sqrt( mA*mA + pGamAbs*pGamAbs);
    double eB       = sqrt( mB*mB + pGamAbs*pGamAbs);
    prodA.p(  pX,  pY,  pZ, eA);
    prodB.p( -pX, -pY, -pZ, eB);

    // Boost to lab frame.
    prodA.bst( pGam, mGam);
    prodB.bst( pGam, mGam);
    prodA.rot( thetaGam, phiGam);
    prodB.rot( thetaGam, phiGam);
    prodA.bst( pDec, decayer.m() );
    prodB.bst( pDec, decayer.m() );
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Translate a partonic content into a set of actual hadrons.

bool ParticleDecays::pickHadrons() {

  // Find partonic decay products. Rest are known id's, mainly hadrons,
  // when necessary shuffled to beginning of idProd list.
  idPartons.resize(0);
  int nPartons = 0;
  int nKnown = 0;
  bool closedGLoop = false;
  for (int i = 1; i <= mult; ++i) {
    int idAbs = abs(idProd[i]);
    if ( idAbs < 9 || (idAbs > 1000 && idAbs < 10000 && (idAbs/10)%10 == 0)
      || idAbs == 81 || idAbs == 82 || idAbs == 83) {
      ++nPartons;
      idPartons.push_back(idProd[i]);
      if (idAbs == 83) closedGLoop = true;
    } else {
      ++nKnown;
      if (nPartons > 0) {
        idProd[nKnown] = idProd[i];
        mProd[nKnown] = mProd[i];
      }
    }
  }

  // Replace generic spectator flavour code by the actual one.
  for (int i = 0; i < nPartons; ++i) {
    int idPart = idPartons[i];
    int idNew = idPart;
    if (idPart == 81) {
      int idAbs = abs(idDec);
      if ( (idAbs/1000)%10 == 0 ) {
        idNew = -(idAbs/10)%10;
        if ((idAbs/100)%2 == 1) idNew = -idNew;
      } else if ( (idAbs/100)%10 >= (idAbs/10)%10 )
        idNew = 100 * ((idAbs/10)%100) + 3;
      else idNew = 1000 * ((idAbs/10)%10) + 100 * ((idAbs/100)%10) + 1;
      if (idDec < 0) idNew = -idNew;

    // Replace generic random flavour by a randomly selected one.
    } else if (idPart == 82 || idPart == 83) {
      double mFlav;
      do {
        int idDummy = -flavSelPtr->pickLightQ();
        FlavContainer flavDummy(idDummy, idPart - 82);
        do idNew = flavSelPtr->pick(flavDummy).id;
        while (idNew == 0);
        mFlav = particleDataPtr->constituentMass(idNew);
      } while (2. * mFlav + stopMass > mProd[0]);
    } else if (idPart == -82 || idPart == -83) {
      idNew = -idPartons[i-1];
    }
    idPartons[i] = idNew;
  }

  // Determine whether fixed multiplicity or to be selected at random.
  int nMin = max( 2, nKnown + nPartons / 2);
  if (meMode == 23) nMin = 3;
  if (meMode > 41 && meMode <= 50) nMin = meMode - 40;
  if (meMode > 51 && meMode <= 60) nMin = meMode - 50;
  int nFix = 0;
  if (meMode == 0) nFix = nMin;
  if (meMode == 11) nFix = 3;
  if (meMode == 12) nFix = 4;
  if (meMode > 61 && meMode <= 70) nFix = meMode - 60;
  if (meMode > 71 && meMode <= 80) nFix = meMode - 70;
  if (nFix > 0 && nKnown + nPartons/2 > nFix) return false;

  // Initial values for loop to set new hadronic content.
  int nFilled, nTotal, nNew, nSpec, nLeft;
  double mDiff;
  int nTry = 0;
  bool diquarkClash = false;
  bool usedChannel  = false;

  // Begin loop; interrupt if multiple tries fail.
  do {
    ++nTry;
    if (nTry > NTRYPICK) return false;

    // Initialize variables inside new try.
    nFilled = nKnown + 1;
    idProd.resize(nFilled);
    mProd.resize(nFilled);
    nTotal = nKnown;
    nSpec = 0;
    nLeft = nPartons;
    mDiff = mProd[0];
    for (int i = 1; i < nFilled; ++i) mDiff -= mProd[i];
    diquarkClash = false;
    usedChannel = false;

    // For weak decays collapse spectators to one particle.
    if ( (meMode == 22 || meMode == 23) && nLeft > 1) {
      FlavContainer flav1( idPartons[nPartons - 2] );
      FlavContainer flav2( idPartons[nPartons - 1] );
      int idHad;
      do idHad = flavSelPtr->combine( flav1, flav2);
      while (idHad == 0);
      double mHad = particleDataPtr->mSel(idHad);
      mDiff -= mHad;
      idProd.push_back( idHad);
      mProd.push_back( mHad);
      ++nFilled;
      nSpec = 1;
      nLeft -= 2;
    }

    // If there are partons left, then determine how many hadrons to pick.
    if (nLeft > 0) {

      // For B -> gamma + X use geometrical distribution.
      if (meMode == 31) {
        double geom = rndmPtr->flat();
        nTotal = 1;
        do {
          ++nTotal;
          geom *= 2.;
        } while (geom < 1. && nTotal < 10);

      // Calculate mass excess and from there average multiplicity.
      } else if (nFix == 0) {
        double multIncreaseNow = (meMode == 23)
          ? multIncreaseWeak : multIncrease;
        double mDiffPS = mDiff;
        for (int i = 0; i < nLeft; ++i)
          mDiffPS -= particleDataPtr->constituentMass( idPartons[i] );
        double average = 0.5 * (nKnown + nSpec) + 0.25 * nPartons
          + multIncreaseNow * log( max( 1.1, mDiffPS / multRefMass ) );
        if (closedGLoop) average += multGoffset;

        // Pick multiplicity according to Poissonian.
        double value = 1.;
        double sum = 1.;
        for (int nNow = nMin + 1; nNow <= 10; ++nNow) {
          value *= average / nNow;
          sum += value;
        }
        nTotal = nMin;
        value = 1.;
        sum *= rndmPtr->flat();
        sum -= value;
        if (sum > 0.) do {
          ++nTotal;
          value *= average / nTotal;
          sum -= value;
        } while (sum > 0. && nTotal < 10);

      // Alternatively predetermined multiplicity.
      } else {
        nTotal = nFix;
      }
      nNew = nTotal - nKnown - nSpec;

      // Set up ends of fragmented system, as copy of idPartons.
      flavEnds.resize(0);
      for (int i = 0; i < nLeft; ++i) {
        flavEnds.push_back( FlavContainer(idPartons[i]) );
        if (abs(idPartons[i]) > 100) flavSelPtr->assignPopQ( flavEnds[i] );
      }

      // Fragment off at random, but save nLeft/2 for final recombination.
      if (nNew > nLeft/2) {
        FlavContainer flavNew;
        int idHad;
        for (int i = 0; i < nNew - nLeft/2; ++i) {
          // When four quarks consider last one to be spectator.
          int iEnd = int( (nLeft - 1.) * rndmPtr->flat() );
          // Pick new flavour and form a new hadron.
          do {
            flavNew = flavSelPtr->pick( flavEnds[iEnd] );
            idHad = flavSelPtr->combine( flavEnds[iEnd], flavNew);
          } while (idHad == 0);
          // Store new hadron and endpoint flavour.
          idProd.push_back( idHad);
          flavEnds[iEnd].anti(flavNew);
        }
      }

      // When only two quarks left, combine to form final hadron.
      if (nLeft == 2) {
        int idHad;
        if ( abs(flavEnds[0].id) > 8 && abs(flavEnds[1].id) > 8)
          diquarkClash = true;
        else {
          do idHad = flavSelPtr->combine( flavEnds[0], flavEnds[1]);
          while (idHad == 0);
          idProd.push_back( idHad);
        }

      // If four quarks, decide how to pair them up.
      } else {
        int iEnd1 = 0;
        int iEnd2 = 1;
        int iEnd3 = 2;
        int iEnd4 = 3;
        if ( rndmPtr->flat() < colRearrange) iEnd2 = 3;
        int relColSign =
          ( (flavEnds[iEnd1].id > 0 && flavEnds[iEnd1].id < 9)
          || flavEnds[iEnd1].id < -10 ) ? 1 : -1;
        if ( (flavEnds[iEnd2].id < 0 && flavEnds[iEnd2].id > -9)
          || flavEnds[iEnd2].id > 10 ) relColSign *= -1;
        if (relColSign == 1) iEnd2 = 2;
        if (iEnd2 == 2) iEnd3 = 1;
        if (iEnd2 == 3) iEnd4 = 1;

        // Then combine to get final two hadrons.
        int idHad;
        if ( abs(flavEnds[iEnd1].id) > 8 && abs(flavEnds[iEnd2].id) > 8)
          diquarkClash = true;
        else {
          do idHad = flavSelPtr->combine( flavEnds[iEnd1], flavEnds[iEnd2]);
          while (idHad == 0);
          idProd.push_back( idHad);
        }
        if ( abs(flavEnds[iEnd3].id) > 8 && abs(flavEnds[iEnd4].id) > 8)
          diquarkClash = true;
        else {
          do idHad = flavSelPtr->combine( flavEnds[iEnd3], flavEnds[iEnd4]);
          while (idHad == 0);
          idProd.push_back( idHad);
        }
      }

      // Find masses of the new hadrons.
      for (int i = nFilled; i < int(idProd.size()) ; ++i) {
        double mHad = particleDataPtr->mSel(idProd[i]);
        mProd.push_back( mHad);
        mDiff -= mHad;
      }
    }

    // Optional: check that this decay mode is not explicitly defined.
    if ( (meMode > 61 && meMode <= 80) && mDiff > mSafety && !diquarkClash ) {
      int idMatch[10], idPNow;
      usedChannel = false;
      bool matched = false;
      // Loop through all channels. Done if not same multiplicity.
      for (int i = 0; i < decDataPtr->sizeChannels(); ++i) {
        DecayChannel& channel = decDataPtr->channel(i);
        if (channel.multiplicity() != nTotal) continue;
        for (int k = 0; k < nTotal; ++k) idMatch[k] = channel.product(k);
        // Match particles one by one until fail.
        // Do not distinguish K0/K0bar/K0short/K0long at this stage.
        for (int j = 0; j < nTotal; ++j) {
          matched = false;
          idPNow = idProd[j + 1];
          if (idPNow == -311 || idPNow == 130 || idPNow == 310) idPNow = 311;
          for (int k = 0; k < nTotal - j; ++k)
          if (idMatch[k] == idPNow || (idMatch[k] == -311 && idPNow == 311)) {
            // Compress list of unmatched when matching worked.
            idMatch[k] = idMatch[nTotal - 1 - j];
            matched = true;
            break;
          }
          if (!matched) break;
        }
        // If matching worked, then chosen channel to be rejected.
        if (matched) {usedChannel = true; break;}
      }
    }

  // Keep on trying until enough phase space and no clash.
  } while (mDiff < mSafety || diquarkClash || usedChannel);

  // Update particle multiplicity.
  mult = idProd.size() - 1;

  // For Dalitz decays shuffle Dalitz pair to the end of the list.
  if (meMode == 11 || meMode == 12) {
    int iL1 = 0;
    int iL2 = 0;
    for (int i = 1; i <= mult; ++i) {
      if (idProd[i] ==  11 || idProd[i] ==  13 || idProd[i] ==  15) iL1 = i;
      if (idProd[i] == -11 || idProd[i] == -13 || idProd[i] == -15) iL2 = i;
    }
    if (iL1 > 0 && iL2 > 0) {
      int idL1 = idProd[iL1];
      int idL2 = idProd[iL2];
      double mL1 = mProd[iL1];
      double mL2 = mProd[iL2];
      int iMove = 0;
      for (int i = 1; i <= mult; ++i) if (i != iL1 && i != iL2) {
        ++iMove;
        idProd[iMove] = idProd[i];
        mProd[iMove] = mProd[i];
      }
      idProd[mult - 1] = idL1;
      idProd[mult] = idL2;
      mProd[mult - 1] = mL1;
      mProd[mult] = mL2;
    }
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Set colour flow and scale in a decay explicitly to partons.

bool ParticleDecays::setColours(Event& event) {

  // Decay to q qbar (or qbar q).
  if (meMode == 91 && idProd[1] > 0 && idProd[1] < 9) {
    int newCol = event.nextColTag();
    cols[1] = newCol;
    acols[2] = newCol;
  } else if (meMode == 91 && idProd[1] < 0 && idProd[1] > -9) {
    int newCol = event.nextColTag();
    cols[2] = newCol;
    acols[1] = newCol;

  // Decay to g g.
  } else if (meMode == 91 && idProd[1] == 21) {
    int newCol1 = event.nextColTag();
    int newCol2 = event.nextColTag();
    cols[1] = newCol1;
    acols[1] = newCol2;
    cols[2] = newCol2;
    acols[2] = newCol1;

  // Decay to g g g.
  } else if (meMode == 92 && idProd[1] == 21 && idProd[2] == 21
    &&  idProd[3] == 21) {
    int newCol1 = event.nextColTag();
    int newCol2 = event.nextColTag();
    int newCol3 = event.nextColTag();
    cols[1] = newCol1;
    acols[1] = newCol2;
    cols[2] = newCol2;
    acols[2] = newCol3;
    cols[3] = newCol3;
    acols[3] = newCol1;

  // Decay to g g gamma: locate which is gamma.
  } else if (meMode == 92) {
    int iGlu1 = (idProd[1] == 21) ? 1 : 3;
    int iGlu2 = (idProd[2] == 21) ? 2 : 3;
    int newCol1 = event.nextColTag();
    int newCol2 = event.nextColTag();
    cols[iGlu1] = newCol1;
    acols[iGlu1] = newCol2;
    cols[iGlu2] = newCol2;
    acols[iGlu2] = newCol1;

  // Unknown decay mode means failure.
  } else return false;

  // Set maximum scale to be mass of decaying particle.
  scale = mProd[0];

  // Done.
  return true;

}

//==========================================================================

} // end namespace Pythia8
