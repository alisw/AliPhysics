// SigmaProcess.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// SigmaProcess class, and classes derived from it.

#include "Pythia8/SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// The SigmaProcess class.
// Base class for cross sections.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Conversion of GeV^{-2} to mb for cross section.
const double SigmaProcess::CONVERT2MB    = 0.389380;

// The sum of outgoing masses must not be too close to the cm energy.
const double SigmaProcess::MASSMARGIN    = 0.1;

// Parameters of momentum rescaling procedure: maximally allowed
// relative energy error and number of iterations.
const double SigmaProcess::COMPRELERR = 1e-10;
const int    SigmaProcess::NCOMPSTEP  = 10;

//--------------------------------------------------------------------------

// Perform simple initialization and store pointers.

void SigmaProcess::init(Info* infoPtrIn, Settings* settingsPtrIn,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn, BeamParticle* beamAPtrIn,
  BeamParticle* beamBPtrIn, Couplings* couplingsPtrIn,
  SigmaTotal* sigmaTotPtrIn, SLHAinterface* slhaInterfacePtrIn) {

  // Store pointers.
  infoPtr         = infoPtrIn;
  settingsPtr     = settingsPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  beamAPtr        = beamAPtrIn;
  beamBPtr        = beamBPtrIn;
  couplingsPtr    = couplingsPtrIn;
  sigmaTotPtr     = sigmaTotPtrIn;
  // Pointer to SLHA object allows semi-internal processes to access
  // SLHA blocks via getEntry() methods, see arXiv:1109.5852
  slhaPtr         = (slhaInterfacePtrIn != 0) ? &slhaInterfacePtrIn->slha : 0;

  // Read out some properties of beams to allow shorthand.
  idA             = (beamAPtr != 0) ? beamAPtr->id() : 0;
  idB             = (beamBPtr != 0) ? beamBPtr->id() : 0;
  mA              = (beamAPtr != 0) ? beamAPtr->m() : 0.;
  mB              = (beamBPtr != 0) ? beamBPtr->m() : 0.;
  isLeptonA       = (beamAPtr != 0) ? beamAPtr->isLepton() : false;
  isLeptonB       = (beamBPtr != 0) ? beamBPtr->isLepton() : false;
  hasLeptonBeams  = isLeptonA || isLeptonB;

  // Photon beams from leptons.
  bool isLepton2gamma = settingsPtr->flag("PDF:lepton2gamma");
  lepton2gammaA   = (beamAPtr != 0) ?
    beamAPtr->isLepton() && isLepton2gamma : false;
  lepton2gammaB   = (beamBPtr != 0) ?
    beamBPtr->isLepton() && isLepton2gamma : false;

  // K factor, multiplying resolved processes. (But not here for MPI.)
  Kfactor         = settingsPtr->parm("SigmaProcess:Kfactor");

  // Maximum incoming quark flavour.
  nQuarkIn        = settingsPtr->mode("PDFinProcess:nQuarkIn");

  // Medium heavy fermion masses set massless or not in ME expressions.
  mcME            = (settingsPtr->flag("SigmaProcess:cMassiveME"))
                  ? particleDataPtr->m0(4)  : 0.;
  mbME            = (settingsPtr->flag("SigmaProcess:bMassiveME"))
                  ? particleDataPtr->m0(5)  : 0.;
  mmuME           = (settingsPtr->flag("SigmaProcess:muMassiveME"))
                  ? particleDataPtr->m0(13) : 0.;
  mtauME          = (settingsPtr->flag("SigmaProcess:tauMassiveME"))
                  ? particleDataPtr->m0(15) : 0.;

  // Renormalization scale choice.
  renormScale1    = settingsPtr->mode("SigmaProcess:renormScale1");
  renormScale2    = settingsPtr->mode("SigmaProcess:renormScale2");
  renormScale3    = settingsPtr->mode("SigmaProcess:renormScale3");
  renormScale3VV  = settingsPtr->mode("SigmaProcess:renormScale3VV");
  renormMultFac   = settingsPtr->parm("SigmaProcess:renormMultFac");
  renormFixScale  = settingsPtr->parm("SigmaProcess:renormFixScale");

  // Factorization scale choice.
  factorScale1    = settingsPtr->mode("SigmaProcess:factorScale1");
  factorScale2    = settingsPtr->mode("SigmaProcess:factorScale2");
  factorScale3    = settingsPtr->mode("SigmaProcess:factorScale3");
  factorScale3VV  = settingsPtr->mode("SigmaProcess:factorScale3VV");
  factorMultFac   = settingsPtr->parm("SigmaProcess:factorMultFac");
  factorFixScale  = settingsPtr->parm("SigmaProcess:factorFixScale");

  // CP violation parameters for the BSM Higgs sector.
  higgsH1parity   = settingsPtr->mode("HiggsH1:parity");
  higgsH1eta      = settingsPtr->parm("HiggsH1:etaParity");
  higgsH1phi      = settingsPtr->parm("HiggsH1:phiParity");
  higgsH2parity   = settingsPtr->mode("HiggsH2:parity");
  higgsH2eta      = settingsPtr->parm("HiggsH2:etaParity");
  higgsH2phi      = settingsPtr->parm("HiggsH2:phiParity");
  higgsA3parity   = settingsPtr->mode("HiggsA3:parity");
  higgsA3eta      = settingsPtr->parm("HiggsA3:etaParity");
  higgsA3phi      = settingsPtr->parm("HiggsA3:phiParity");

  // If BSM not switched on then H1 should have SM properties.
  if (!settingsPtr->flag("Higgs:useBSM")){
    higgsH1parity = 1;
    higgsH1eta    = 0.;
    higgsH1phi    = M_PI / 2.;
  }

}

//--------------------------------------------------------------------------

// Set up allowed flux of incoming partons.
// addBeam: set up PDF's that need to be evaluated for the two beams.
// addPair: set up pairs of incoming partons from the two beams.

bool SigmaProcess::initFlux() {

  // Reset arrays (in case of several init's in same run).
  inBeamA.clear();
  inBeamB.clear();
  inPair.clear();

  // Read in process-specific channel information.
  string fluxType = inFlux();

  // Case with g g incoming state.
  if (fluxType == "gg") {
    addBeamA(21);
    addBeamB(21);
    addPair(21, 21);
  }

  // Case with q g incoming state.
  else if (fluxType == "qg") {
    for (int i = -nQuarkIn; i <= nQuarkIn; ++i) {
      int idNow = (i == 0) ? 21 : i;
      addBeamA(idNow);
      addBeamB(idNow);
    }
    for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
    if (idNow != 0) {
      addPair(idNow, 21);
      addPair(21, idNow);
    }
  }

  // Case with q q', q qbar' or qbar qbar' incoming state.
  else if (fluxType == "qq") {
    for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
    if (idNow != 0) {
      addBeamA(idNow);
      addBeamB(idNow);
    }
    for (int id1Now = -nQuarkIn; id1Now <= nQuarkIn; ++id1Now)
    if (id1Now != 0)
    for (int id2Now = -nQuarkIn; id2Now <= nQuarkIn; ++id2Now)
    if (id2Now != 0)
      addPair(id1Now, id2Now);
  }

  // Case with q qbar' incoming state.
  else if (fluxType == "qqbar") {
    for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
    if (idNow != 0) {
      addBeamA(idNow);
      addBeamB(idNow);
    }
    for (int id1Now = -nQuarkIn; id1Now <= nQuarkIn; ++id1Now)
    if (id1Now != 0)
    for (int id2Now = -nQuarkIn; id2Now <= nQuarkIn; ++id2Now)
    if (id2Now != 0 && id1Now * id2Now < 0)
      addPair(id1Now, id2Now);
  }

  // Case with q qbar incoming state.
  else if (fluxType == "qqbarSame") {
    for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
    if (idNow != 0) {
      addBeamA(idNow);
      addBeamB(idNow);
    }
    for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
    if (idNow != 0)
      addPair(idNow, -idNow);
  }

  // Case with f f', f fbar', fbar fbar' incoming state.
  else if (fluxType == "ff") {
    // If beams are leptons then they are also the colliding partons
    // unless lepton includes a photon beam.
    if ( isLeptonA && isLeptonB && !lepton2gammaA && !lepton2gammaB ) {
      addBeamA(idA);
      addBeamB(idB);
      addPair(idA, idB);
    // First beam is lepton and second is hadron.
    } else if ( isLeptonA && !lepton2gammaA ) {
      addBeamA(idA);
      for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamB(idNow);
        addPair(idA, idNow);
      }
    // First beam is hadron and second is lepton.
    } else if ( isLeptonB && !lepton2gammaB ) {
      addBeamB(idB);
      for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamA(idNow);
        addPair(idNow, idB);
      }
    // Hadron beams gives quarks.
    } else {
      for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamA(idNow);
        addBeamB(idNow);
      }
      for (int id1Now = -nQuarkIn; id1Now <= nQuarkIn; ++id1Now)
      if (id1Now != 0)
      for (int id2Now = -nQuarkIn; id2Now <= nQuarkIn; ++id2Now)
      if (id2Now != 0)
        addPair(id1Now, id2Now);
    }
  }

  // Case with f fbar' generic incoming state.
  else if (fluxType == "ffbar") {
    // If beams are leptons then also colliding partons
    // unless lepton includes a photon beam.
    if (isLeptonA && isLeptonB && idA * idB < 0
        && !lepton2gammaA && !lepton2gammaB) {
      addBeamA(idA);
      addBeamB(idB);
      addPair(idA, idB);
    // Hadron beams gives quarks.
    } else {
      for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamA(idNow);
        addBeamB(idNow);
      }
      for (int id1Now = -nQuarkIn; id1Now <= nQuarkIn; ++id1Now)
      if (id1Now != 0)
      for (int id2Now = -nQuarkIn; id2Now <= nQuarkIn; ++id2Now)
      if (id2Now != 0 && id1Now * id2Now < 0)
        addPair(id1Now, id2Now);
    }
  }

  // Case with f fbar incoming state.
  else if (fluxType == "ffbarSame") {
    // If beams are antiparticle pair and leptons then also colliding partons
    // unless lepton includes a photon beam.
    if ( idA + idB == 0 && isLeptonA && !lepton2gammaA && !lepton2gammaB) {
      addBeamA(idA);
      addBeamB(idB);
      addPair(idA, idB);
    // Else assume both to be hadrons, for better or worse.
    } else {
      for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamA(idNow);
        addBeamB(idNow);
      }
      for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0)
        addPair(idNow, -idNow);
    }
  }

  // Case with f fbar' charged(+-1) incoming state.
  else if (fluxType == "ffbarChg") {
    // If beams are leptons then also colliding partons
    // unless lepton includes a photon beam.
    if ( isLeptonA && isLeptonB && !lepton2gammaA && !lepton2gammaB
         && abs( particleDataPtr->chargeType(idA)
           + particleDataPtr->chargeType(idB) ) == 3 ) {
      addBeamA(idA);
      addBeamB(idB);
      addPair(idA, idB);
    // Hadron beams gives quarks.
    } else {
      for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamA(idNow);
        addBeamB(idNow);
      }
      for (int id1Now = -nQuarkIn; id1Now <= nQuarkIn; ++id1Now)
      if (id1Now != 0)
      for (int id2Now = -nQuarkIn; id2Now <= nQuarkIn; ++id2Now)
      if (id2Now != 0 && id1Now * id2Now < 0
        && (abs(id1Now) + abs(id2Now))%2 == 1) addPair(id1Now, id2Now);
    }
  }

  // Case with f gamma incoming state.
  else if (fluxType == "fgm") {
    // Fermion from incoming side A if no photon beam inside.
    if ( isLeptonA && !lepton2gammaA ) {
      addBeamA( idA);
      addPair(idA, 22);
    } else {
      for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamA(idNow);
        addPair(idNow, 22);
      }
    }
    // Fermion from incoming side B if no photon beam inside.
    if ( isLeptonB && !lepton2gammaB ) {
      addBeamB( idB);
      addPair(22, idB);
    } else {
      for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamB(idNow);
        addPair(22, idNow);
      }
    }
    // Photons in the beams.
    addBeamA(22);
    addBeamB(22);
  }

  // Case with quark gamma incoming state.
  else if (fluxType == "qgm") {
    for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamA(idNow);
        addPair(idNow, 22);
      }
    for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamB(idNow);
        addPair(22, idNow);
      }

    // Photons in the beams.
    addBeamA(22);
    addBeamB(22);
  }

  // Case with gamma quark incoming state.
  // Need this when both resolved and unresolved photon beams.
  else if (fluxType == "gmq") {
    for (int idNow = -nQuarkIn; idNow <= nQuarkIn; ++idNow)
      if (idNow != 0) {
        addBeamB(idNow);
        addPair(22, idNow);
      }
    // Photon in the beam.
    addBeamA(22);
  }

  // Case with gluon gamma incoming state.
  else if (fluxType == "ggm") {
    addBeamA(21);
    addBeamA(22);
    addBeamB(21);
    addBeamB(22);
    addPair(21, 22);
    addPair(22, 21);
  }

  // Case with gamma gluon incoming state.
  // Need this when both resolved and unresolved photon beams.
  else if (fluxType == "gmg") {
    addBeamA(22);
    addBeamB(21);
    addPair(22, 21);
  }

  // Case with gamma gamma incoming state.
  else if (fluxType == "gmgm") {
    addBeamA(22);
    addBeamB(22);
    addPair(22, 22);
  }

  // Unrecognized fluxType is bad sign. Else done.
  else {
    infoPtr->errorMsg("Error in SigmaProcess::initFlux: "
    "unrecognized inFlux type", fluxType);
    return false;
  }
  return true;

}

//--------------------------------------------------------------------------

// Convolute matrix-element expression(s) with parton flux and K factor.
// Possibly different PDFs for the phase-space initialization.
// Can also take new values for x's to correct for oversampling, as
// needed with external photon flux.

double SigmaProcess::sigmaPDF(bool initPS, bool samexGamma,
    bool useNewXvalues, double x1New, double x2New) {

  // Evaluate and store the required parton densities.
  for (int j = 0; j < sizeBeamA(); ++j) {
    if ( initPS)
      inBeamA[j].pdf = beamAPtr->xfMax( inBeamA[j].id, x1Save, Q2FacSave);
    else if ( samexGamma)
      inBeamA[j].pdf = beamAPtr->xfSame( inBeamA[j].id, x1Save, Q2FacSave);
    else if ( useNewXvalues && x1New > 0.)
      inBeamA[j].pdf = beamAPtr->xfGamma( inBeamA[j].id, x1New, Q2FacSave);
    else
      inBeamA[j].pdf = beamAPtr->xfHard( inBeamA[j].id, x1Save, Q2FacSave);
  }
  for (int j = 0; j < sizeBeamB(); ++j){
    if ( initPS)
      inBeamB[j].pdf = beamBPtr->xfMax( inBeamB[j].id, x2Save, Q2FacSave);
    else if ( samexGamma)
      inBeamB[j].pdf = beamBPtr->xfSame( inBeamB[j].id, x2Save, Q2FacSave);
    else if ( useNewXvalues && x2New > 0.)
      inBeamB[j].pdf = beamBPtr->xfGamma( inBeamB[j].id, x2New, Q2FacSave);
    else
      inBeamB[j].pdf = beamBPtr->xfHard( inBeamB[j].id, x2Save, Q2FacSave);
  }

  // Save the x_gamma values after PDFs are called if new value is sampled
  // if using internal photon flux from leptons.
  if ( !useNewXvalues && !samexGamma && beamAPtr->hasResGamma() )
    beamAPtr->xGammaPDF();
  if ( !useNewXvalues && !samexGamma && beamBPtr->hasResGamma() )
    beamBPtr->xGammaPDF();

  // Loop over allowed incoming channels.
  sigmaSumSave = 0.;
  for (int i = 0; i < sizePair(); ++i) {

    // Evaluate hard-scattering cross section. Include K factor.
    inPair[i].pdfSigma = Kfactor
                       * sigmaHatWrap(inPair[i].idA, inPair[i].idB);

    // Multiply by respective parton densities.
    for (int j = 0; j < sizeBeamA(); ++j)
    if (inPair[i].idA == inBeamA[j].id) {
      inPair[i].pdfA      = inBeamA[j].pdf;
      inPair[i].pdfSigma *= inBeamA[j].pdf;
      break;
    }
    for (int j = 0; j < sizeBeamB(); ++j)
    if (inPair[i].idB == inBeamB[j].id) {
      inPair[i].pdfB      = inBeamB[j].pdf;
      inPair[i].pdfSigma *= inBeamB[j].pdf;
      break;
    }

    // Sum for all channels.
    sigmaSumSave += inPair[i].pdfSigma;
  }

  // Done.
  return sigmaSumSave;

}

//--------------------------------------------------------------------------

// Select incoming parton channel and extract parton densities (resolved).

void SigmaProcess::pickInState(int id1in, int id2in) {

  // Multiparton interactions: partons already selected.
  if (id1in != 0 && id2in != 0) {
    id1 = id1in;
    id2 = id2in;
    return;
  }

  // Pick channel. Extract channel flavours and pdf's.
  double sigmaRand =  sigmaSumSave * rndmPtr->flat();
  for (int i = 0; i < sizePair(); ++i) {
    sigmaRand -= inPair[i].pdfSigma;
    if (sigmaRand <= 0.) {
      id1      = inPair[i].idA;
      id2      = inPair[i].idB;
      pdf1Save = inPair[i].pdfA;
      pdf2Save = inPair[i].pdfB;
      break;
    }
  }

}

//--------------------------------------------------------------------------

// Calculate incoming modified masses and four-vectors for matrix elements.

bool SigmaProcess::setupForMEin() {

  // Initially assume it will work out to set up modified kinematics.
  bool allowME = true;

  // Correct incoming c, b, mu and tau to be massive or not.
  mME[0] = 0.;
  int id1Tmp = abs(id1);
  if (id1Tmp ==  4) mME[0] = mcME;
  if (id1Tmp ==  5) mME[0] = mbME;
  if (id1Tmp == 13) mME[0] = mmuME;
  if (id1Tmp == 15) mME[0] = mtauME;
  mME[1] = 0.;
  int id2Tmp = abs(id2);
  if (id2Tmp ==  4) mME[1] = mcME;
  if (id2Tmp ==  5) mME[1] = mbME;
  if (id2Tmp == 13) mME[1] = mmuME;
  if (id2Tmp == 15) mME[1] = mtauME;

  // If kinematically impossible return to massless case, but set error.
  if (mME[0] + mME[1] >= mH) {
    mME[0] = 0.;
    mME[1] = 0.;
    allowME = false;
  }

  // Do incoming two-body kinematics for massless or massive cases.
  if (mME[0] == 0. && mME[1] == 0.) {
  pME[0] = 0.5 * mH * Vec4( 0., 0.,  1., 1.);
  pME[1] = 0.5 * mH * Vec4( 0., 0., -1., 1.);
  } else {
    double e0   = 0.5 * (mH * mH + mME[0] * mME[0] - mME[1] * mME[1]) / mH;
    double pz0  = sqrtpos(e0 * e0 - mME[0] * mME[0]);
    pME[0] = Vec4( 0., 0.,  pz0, e0);
    pME[1] = Vec4( 0., 0., -pz0, mH - e0);
  }

  // Done.
  return allowME;

}

//--------------------------------------------------------------------------

// Evaluate weight for W decay distribution in t -> W b -> f fbar b.

double SigmaProcess::weightTopDecay( Event& process, int iResBeg,
  int iResEnd) {

  // If not pair W d/s/b and mother t then return unit weight.
  if (iResEnd - iResBeg != 1) return 1.;
  int iW1  = iResBeg;
  int iB2  = iResBeg + 1;
  int idW1 = process[iW1].idAbs();
  int idB2 = process[iB2].idAbs();
  if (idW1 != 24) {
    swap(iW1, iB2);
    swap(idW1, idB2);
  }
  if (idW1 != 24 || (idB2 != 1 && idB2 != 3 && idB2 != 5)) return 1.;
  int iT   = process[iW1].mother1();
  if (iT <= 0 || process[iT].idAbs() != 6) return 1.;

  // Find sign-matched order of W decay products.
  int iF    = process[iW1].daughter1();
  int iFbar = process[iW1].daughter2();
  if (iFbar - iF != 1) return 1.;
  if (process[iT].id() * process[iF].id() < 0) swap(iF, iFbar);

  // Weight and maximum weight.
  double wt    = (process[iT].p() * process[iFbar].p())
               * (process[iF].p() * process[iB2].p());
  double wtMax = ( pow4(process[iT].m()) - pow4(process[iW1].m()) ) / 8.;

  // Done.
  return wt / wtMax;

}

//--------------------------------------------------------------------------

// Evaluate weight for Z0/W+- decay distributions in H -> Z0/W+ Z0/W- -> 4f
// and H -> gamma Z0 -> gamma f fbar.

double SigmaProcess::weightHiggsDecay( Event& process, int iResBeg,
  int iResEnd) {

  // If not pair Z0 Z0, W+ W- or gamma Z0 then return unit weight.
  if (iResEnd - iResBeg != 1) return 1.;
  int iZW1  = iResBeg;
  int iZW2  = iResBeg + 1;
  int idZW1 = process[iZW1].id();
  int idZW2 = process[iZW2].id();
  if (idZW1 < 0 || idZW2 == 22) {
    swap(iZW1, iZW2);
    swap(idZW1, idZW2);
  }
  if ( (idZW1 != 23 || idZW2 != 23) && (idZW1 != 24 || idZW2 != -24)
    && (idZW1 != 22 || idZW2 != 23) ) return 1.;

  // If mother is not Higgs then return unit weight.
  int iH  = process[iZW1].mother1();
  if (iH <= 0) return 1.;
  int idH = process[iH].id();
  if (idH != 25 && idH != 35 && idH !=36) return 1.;

  // H -> gamma Z0 -> gamma f fbar is 1 + cos^2(theta) in Z rest frame.
  if (idZW1 == 22) {
    int i5 = process[iZW2].daughter1();
    int i6 = process[iZW2].daughter2();
    double pgmZ = process[iZW1].p() * process[iZW2].p();
    double pgm5 = process[iZW1].p() * process[i5].p();
    double pgm6 = process[iZW1].p() * process[i6].p();
    return (pow2(pgm5) + pow2(pgm6)) / pow2(pgmZ);
  }

  // Parameters depend on Higgs type: H0(H_1), H^0(H_2) or A^0(H_3).
  int    higgsParity = higgsH1parity;
  double higgsEta    = higgsH1eta;
  if (idH == 35) {
    higgsParity      = higgsH2parity;
    higgsEta         = higgsH2eta;
  } else if (idH == 36) {
    higgsParity      = higgsA3parity;
    higgsEta         = higgsA3eta;
  }

  // Option with isotropic decays (also for pseudoscalar fermion couplings).
  if (higgsParity == 0 || higgsParity > 3) return 1.;

  // Maximum and initial weight.
  double wtMax = pow4(process[iH].m());
  double wt    = wtMax;

  // Find sign-matched order of Z0/W+- decay products.
  int i3 = process[iZW1].daughter1();
  int i4 = process[iZW1].daughter2();
  if (process[i3].id() < 0) swap( i3, i4);
  int i5 = process[iZW2].daughter1();
  int i6 = process[iZW2].daughter2();
  if (process[i5].id() < 0) swap( i5, i6);

  // Evaluate four-vector products and find masses..
  double p35  = 2. * process[i3].p() * process[i5].p();
  double p36  = 2. * process[i3].p() * process[i6].p();
  double p45  = 2. * process[i4].p() * process[i5].p();
  double p46  = 2. * process[i4].p() * process[i6].p();
  double p34  = 2. * process[i3].p() * process[i4].p();
  double p56  = 2. * process[i5].p() * process[i6].p();
  double mZW1 = process[iZW1].m();
  double mZW2 = process[iZW2].m();

  // For mixed CP states need epsilon product and gauge boson masses.
  double epsilonProd = 0.;
  if (higgsParity == 3) {
    double p[4][4];
    for (int i = 0; i < 4; ++i) {
      int         ii = i3;
      if (i == 1) ii = i4;
      if (i == 2) ii = i5;
      if (i == 3) ii = i6;
      p[i][0] = process[ii].e();
      p[i][1] = process[ii].px();
      p[i][2] = process[ii].py();
      p[i][3] = process[ii].pz();
    }
    epsilonProd
      = p[0][0]*p[1][1]*p[2][2]*p[3][3] - p[0][0]*p[1][1]*p[2][3]*p[3][2]
      - p[0][0]*p[1][2]*p[2][1]*p[3][3] + p[0][0]*p[1][2]*p[2][3]*p[3][1]
      + p[0][0]*p[1][3]*p[2][1]*p[3][2] - p[0][0]*p[1][3]*p[2][2]*p[3][1]
      - p[0][1]*p[1][0]*p[2][2]*p[3][3] + p[0][1]*p[1][0]*p[2][3]*p[3][2]
      + p[0][1]*p[1][2]*p[2][0]*p[3][3] - p[0][1]*p[1][2]*p[2][3]*p[3][0]
      - p[0][1]*p[1][3]*p[2][0]*p[3][2] + p[0][1]*p[1][3]*p[2][2]*p[3][0]
      + p[0][2]*p[1][0]*p[2][1]*p[3][3] - p[0][2]*p[1][0]*p[2][3]*p[3][1]
      - p[0][2]*p[1][1]*p[2][0]*p[3][3] + p[0][2]*p[1][1]*p[2][3]*p[3][0]
      + p[0][2]*p[1][3]*p[2][0]*p[3][1] - p[0][2]*p[1][3]*p[2][1]*p[3][0]
      - p[0][3]*p[1][0]*p[2][1]*p[3][2] + p[0][3]*p[1][0]*p[2][2]*p[3][1]
      + p[0][3]*p[1][1]*p[2][0]*p[3][2] - p[0][3]*p[1][1]*p[2][2]*p[3][0]
      - p[0][3]*p[1][2]*p[2][0]*p[3][1] + p[0][3]*p[1][2]*p[2][1]*p[3][0];
  }

  // Z0 Z0 decay: vector and axial couplings of two fermion pairs.
  if (idZW1 == 23) {
    double vf1 = couplingsPtr->vf(process[i3].idAbs());
    double af1 = couplingsPtr->af(process[i3].idAbs());
    double vf2 = couplingsPtr->vf(process[i5].idAbs());
    double af2 = couplingsPtr->af(process[i5].idAbs());
    double va12asym = 4. * vf1 * af1 * vf2 * af2
      / ( (vf1*vf1 + af1*af1) * (vf2*vf2 + af2*af2) );
    double vh = 1;
    double ah = higgsEta / pow2( particleDataPtr->m0(23) );

    // Normal CP-even decay.
    if (higgsParity == 1) wt = 8. * (1. + va12asym) * p35 * p46
      + 8. * (1. - va12asym) * p36 * p45;

    // CP-odd decay (normal for A0(H_3)).
    else if (higgsParity == 2) wt = ( pow2(p35 + p46)
      + pow2(p36 + p45) - 2. * p34 * p56
      - 2. * pow2(p35 * p46 - p36 * p45) / (p34 * p56)
      + va12asym * (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) )
      / (1. +  va12asym);

    // Mixed CP states.
    else wt = 32. * ( 0.25 * pow2(vh) * ( (1. + va12asym) * p35 * p46
      + (1. - va12asym) * p36 * p45 ) - 0.5 * vh * ah * epsilonProd
      * ( (1. + va12asym) * (p35 + p46) - (1. - va12asym) * (p36 + p45) )
      + 0.0625 * pow2(ah) * (-2. * pow2(p34 * p56)
      - 2. * pow2(p35 * p46 - p36 * p45)
      + p34 * p56 * (pow2(p35 + p46) + pow2(p36 + p45))
      + va12asym * p34 * p56 * (p35 + p36 - p45 - p46)
      * (p35 + p45 - p36 - p46) ) )
      / ( pow2(vh) + 2. * abs(vh * ah) * mZW1 * mZW2
      + 2. * pow2(ah * mZW1 * mZW2) * (1. + va12asym) );

  // W+ W- decay.
  } else if (idZW1 == 24) {
    double vh = 1;
    double ah = higgsEta / pow2( particleDataPtr->m0(24) );

    // Normal CP-even decay.
    if (higgsParity == 1) wt = 16. * p35 * p46;

    // CP-odd decay (normal for A0(H_3)).
    else if (higgsParity == 2) wt = 0.5 * ( pow2(p35 + p46)
      + pow2(p36 + p45) - 2. * p34 * p56
      - 2. * pow2(p35 * p46 - p36 * p45) / (p34 * p56)
      + (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) );

    // Mixed CP states.
    else wt = 32. * ( 0.25 * pow2(vh) * 2. * p35 * p46
      - 0.5 * vh * ah * epsilonProd * 2. * (p35 + p46)
      + 0.0625 * pow2(ah) * (-2. * pow2(p34 * p56)
      - 2. * pow2(p35 * p46 - p36 * p45)
      + p34 * p56 * (pow2(p35 + p46) + pow2(p36 + p45))
      + p34 * p56 * (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) ) )
      / ( pow2(vh) + 2. * abs(vh * ah) * mZW1 * mZW2
      + 2. * pow2(ah * mZW1 * mZW2) );
  }

  // Done.
  return wt / wtMax;

}

//==========================================================================

// The Sigma1Process class.
// Base class for resolved 2 -> 1 cross sections; derived from SigmaProcess.

//--------------------------------------------------------------------------

// Wrapper to sigmaHat, to (a) store current incoming flavours,
// (b) convert from GeV^-2 to mb where required, and
// (c) convert from |M|^2 to d(sigmaHat)/d(tHat) where required.

double Sigma1Process::sigmaHatWrap(int id1in, int id2in) {

  id1 = id1in;
  id2 = id2in;
  double sigmaTmp = sigmaHat();
  if (convertM2()) {
    sigmaTmp /= 2. * sH;
    // Convert 2 * pi * delta(p^2 - m^2) to Breit-Wigner with same area.
    int idTmp     = resonanceA();
    double mTmp   = particleDataPtr->m0(idTmp);
    double GamTmp = particleDataPtr->mWidth(idTmp);
    sigmaTmp *= 2. * mTmp * GamTmp / ( pow2(sH - mTmp * mTmp)
      + pow2(mTmp * GamTmp) );
  }
  if (convert2mb()) sigmaTmp *= CONVERT2MB;
  return sigmaTmp;

}

//--------------------------------------------------------------------------

// Input and complement kinematics for resolved 2 -> 1 process.

void Sigma1Process::store1Kin( double x1in, double x2in, double sHin) {

  // Default value only sensible for these processes.
  swapTU = false;

  // Incoming parton momentum fractions and sHat.
  x1Save = x1in;
  x2Save = x2in;
  sH     = sHin;
  mH     = sqrt(sH);
  sH2    = sH * sH;

  // Different options for renormalization scale, but normally sHat.
  Q2RenSave                        = renormMultFac * sH;
  if (renormScale1 == 2) Q2RenSave = renormFixScale;

  // Different options for factorization scale, but normally sHat.
  Q2FacSave                        = factorMultFac * sH;
  if (factorScale1 == 2) Q2FacSave = factorFixScale;

  // Evaluate alpha_strong and alpha_EM.
  alpS   = couplingsPtr->alphaS(Q2RenSave);
  alpEM  = couplingsPtr->alphaEM(Q2RenSave);

}

//--------------------------------------------------------------------------

// Calculate modified masses and four-vectors for matrix elements.

bool Sigma1Process::setupForME() {

  // Common initial-state handling.
  bool allowME = setupForMEin();

  // Final state trivial here.
  mME[2] = mH;
  pME[2] = Vec4( 0., 0., 0., mH);

  // Done.
  return allowME;

}

//==========================================================================

// The Sigma2Process class.
// Base class for resolved 2 -> 2 cross sections; derived from SigmaProcess.

//--------------------------------------------------------------------------

// Input and complement kinematics for resolved 2 -> 2 process.

void Sigma2Process::store2Kin( double x1in, double x2in, double sHin,
  double tHin, double m3in, double m4in, double runBW3in, double runBW4in) {

  // Default ordering of particles 3 and 4.
  swapTU   = false;

  // Incoming parton momentum fractions.
  x1Save   = x1in;
  x2Save   = x2in;

  // Incoming masses and their squares.
  bool masslessKin = (id3Mass() == 0) && (id4Mass() == 0);
  if (masslessKin) {
    m3     = 0.;
    m4     = 0.;
  } else {
    m3     = m3in;
    m4     = m4in;
  }
  mSave[3] = m3;
  mSave[4] = m4;
  s3       = m3 * m3;
  s4       = m4 * m4;

  // Standard Mandelstam variables and their squares.
  sH       = sHin;
  tH       = tHin;
  uH       = (masslessKin) ? -(sH + tH) : s3 + s4 - (sH + tH);
  mH       = sqrt(sH);
  sH2      = sH * sH;
  tH2      = tH * tH;
  uH2      = uH * uH;

  // The nominal Breit-Wigner factors with running width.
  runBW3   = runBW3in;
  runBW4   = runBW4in;

  // Calculate squared transverse momentum.
  pT2 = (masslessKin) ?  tH * uH / sH : (tH * uH - s3 * s4) / sH;

  // Special case: pick scale as if 2 -> 1 process in disguise.
  if (isSChannel()) {

    // Different options for renormalization scale, but normally sHat.
    Q2RenSave                        = renormMultFac * sH;
    if (renormScale1 == 2) Q2RenSave = renormFixScale;

    // Different options for factorization scale, but normally sHat.
    Q2FacSave                        = factorMultFac * sH;
    if (factorScale1 == 2) Q2FacSave = factorFixScale;

  // Normal case with "true" 2 -> 2.
  } else {

    // Different options for renormalization scale.
    if (masslessKin)            Q2RenSave = (renormScale2 < 4) ? pT2 : sH;
    else if (renormScale2 == 1) Q2RenSave = pT2 + min(s3, s4);
    else if (renormScale2 == 2) Q2RenSave = sqrt((pT2 + s3) * (pT2 + s4));
    else if (renormScale2 == 3) Q2RenSave = pT2 + 0.5 * (s3 + s4);
    else                        Q2RenSave = sH;
    Q2RenSave                            *= renormMultFac;
    if      (renormScale2 == 5) Q2RenSave = renormFixScale;
    if      (renormScale2 == 6) Q2RenSave = -tH * renormMultFac;

    // Different options for factorization scale.
    if (masslessKin)            Q2FacSave = (factorScale2 < 4) ? pT2 : sH;
    else if (factorScale2 == 1) Q2FacSave = pT2 + min(s3, s4);
    else if (factorScale2 == 2) Q2FacSave = sqrt((pT2 + s3) * (pT2 + s4));
    else if (factorScale2 == 3) Q2FacSave = pT2 + 0.5 * (s3 + s4);
    else                        Q2FacSave = sH;
    Q2FacSave                            *= factorMultFac;
    if      (factorScale2 == 5) Q2FacSave = factorFixScale;
    if      (factorScale2 == 6) Q2FacSave = -tH * factorMultFac;
  }

  // Evaluate alpha_strong and alpha_EM.
  alpS  = couplingsPtr->alphaS(Q2RenSave);
  alpEM = couplingsPtr->alphaEM(Q2RenSave);

}

//--------------------------------------------------------------------------

// As above, special kinematics for multiparton interactions.

void Sigma2Process::store2KinMPI( double x1in, double x2in,
  double sHin, double tHin, double uHin, double alpSin, double alpEMin,
  bool needMasses, double m3in, double m4in) {

  // Default ordering of particles 3 and 4.
  swapTU    = false;

  // Incoming x values.
  x1Save    = x1in;
  x2Save    = x2in;

  // Standard Mandelstam variables and their squares.
  sH        = sHin;
  tH        = tHin;
  uH        = uHin;
  mH        = sqrt(sH);
  sH2       = sH * sH;
  tH2       = tH * tH;
  uH2       = uH * uH;

  // Strong and electroweak couplings.
  alpS      = alpSin;
  alpEM     = alpEMin;

  // Assume vanishing masses. (Will be modified in final kinematics.)
  m3        = 0.;
  s3        = 0.;
  m4        = 0.;
  s4        = 0.;
  sHBeta    = sH;

  // Scattering angle.
  cosTheta  = (tH - uH) / sH;
  sinTheta  = 2. * sqrtpos( tH * uH ) / sH;

  // In some cases must use masses and redefine meaning of tHat and uHat.
  if (needMasses) {
    m3      = m3in;
    s3      = m3 * m3;
    m4      = m4in;
    s4      = m4 * m4;
    sHMass  = sH - s3 - s4;
    sHBeta  = sqrtpos(sHMass*sHMass - 4. * s3 * s4);
    tH      = -0.5 * (sHMass - sHBeta * cosTheta);
    uH      = -0.5 * (sHMass + sHBeta * cosTheta);
    tH2     = tH * tH;
    uH2     = uH * uH;
  }

  // pT2 with masses (at this stage) included.
  pT2Mass   = 0.25 * sHBeta * pow2(sinTheta);

}

//--------------------------------------------------------------------------

// Perform kinematics for a multiparton interaction, including a rescattering.

bool Sigma2Process::final2KinMPI( int i1Res, int i2Res, Vec4 p1Res, Vec4 p2Res,
  double m1Res, double m2Res) {

  // Have to set flavours and colours.
  setIdColAcol();

  // Check that masses of outgoing particles not too big.
  if (m3 == 0.) m3  = particleDataPtr->m0(idSave[3]);
  if (m4 == 0.) m4  = particleDataPtr->m0(idSave[4]);
  mH           = sqrt(sH);
  if (m3 + m4 + MASSMARGIN > mH) return false;
  s3           = m3 * m3;
  s4           = m4 * m4;

  // Do kinematics of the production; without or with masses.
  double e1In  = 0.5 * mH;
  double e2In  = e1In;
  double pzIn  = e1In;
  if (i1Res > 0 || i2Res > 0) {
    double s1  = m1Res * m1Res;
    double s2  = m2Res * m2Res;
    e1In       = 0.5 * (sH + s1 - s2) / mH;
    e2In       = 0.5 * (sH + s2 - s1) / mH;
    pzIn       = sqrtpos( e1In*e1In - s1 );
  }

  // Do kinematics of the decay.
  double e3    = 0.5 * (sH + s3 - s4) / mH;
  double e4    = 0.5 * (sH + s4 - s3) / mH;
  double pAbs  = sqrtpos( e3*e3 - s3 );
  phi          = 2. * M_PI * rndmPtr->flat();
  double pZ    = pAbs * cosTheta;
  pTFin        = pAbs * sinTheta;
  double pX    = pTFin * sin(phi);
  double pY    = pTFin * cos(phi);
  double scale = 0.5 * mH * sinTheta;
  if (swappedTU()) pZ = -pZ;

  // Fill particle info.
  int status1  = (i1Res == 0) ? -31 : -34;
  int status2  = (i2Res == 0) ? -31 : -34;
  parton[1]    = Particle( idSave[1], status1, 0, 0, 3, 4,
    colSave[1], acolSave[1],  0.,  0.,  pzIn, e1In, m1Res, scale);
  parton[2]    = Particle( idSave[2], status2, 0, 0, 3, 4,
    colSave[2], acolSave[2],  0.,  0., -pzIn, e2In, m2Res, scale);
  parton[3]    = Particle( idSave[3],      33, 1, 2, 0, 0,
    colSave[3], acolSave[3],  pX,  pY,    pZ,   e3,    m3, scale);
  parton[4]    = Particle( idSave[4],      33, 1, 2, 0, 0,
    colSave[4], acolSave[4], -pX, -pY,   -pZ,   e4,    m4, scale);

  // Boost particles from subprocess rest frame to event rest frame.
  // Normal multiparton interaction: only longitudinal boost.
  if (i1Res == 0 && i2Res == 0) {
    double betaZ = (x1Save - x2Save) / (x1Save + x2Save);
    for (int i = 1; i <= 4; ++i) parton[i].bst(0., 0., betaZ);
  // Rescattering: generic rotation and boost required.
  } else {
    RotBstMatrix M;
    M.fromCMframe( p1Res, p2Res);
    for (int i = 1; i <= 4; ++i) parton[i].rotbst(M);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Calculate modified masses and four-vectors for matrix elements.

bool Sigma2Process::setupForME() {

  // Common initial-state handling.
  bool allowME = setupForMEin();

  // Correct outgoing c, b, mu and tau to be massive or not.
  mME[2] = m3;
  int id3Tmp = abs(id3Mass());
  if (id3Tmp ==  4) mME[2] = mcME;
  if (id3Tmp ==  5) mME[2] = mbME;
  if (id3Tmp == 13) mME[2] = mmuME;
  if (id3Tmp == 15) mME[2] = mtauME;
  mME[3] = m4;
  int id4Tmp = abs(id4Mass());
  if (id4Tmp ==  4) mME[3] = mcME;
  if (id4Tmp ==  5) mME[3] = mbME;
  if (id4Tmp == 13) mME[3] = mmuME;
  if (id4Tmp == 15) mME[3] = mtauME;

  // If kinematically impossible turn to massless case, but set error.
  if (mME[2] + mME[3] >= mH) {
    mME[2] = 0.;
    mME[3] = 0.;
    allowME = false;
  }

  // Calculate scattering angle in subsystem rest frame.
  double sH34 = sqrtpos( pow2(sH - s3 - s4) - 4. * s3 * s4);
  double cThe = (tH - uH) / sH34;
  double sThe = sqrtpos(1. - cThe * cThe);

  // Setup massive kinematics with preserved scattering angle.
  double s3ME   = pow2(mME[2]);
  double s4ME   = pow2(mME[3]);
  double sH34ME = sqrtpos( pow2(sH - s3ME - s4ME) - 4. * s3ME * s4ME);
  double pAbsME = 0.5 * sH34ME / mH;

  // Normally allowed with unequal (or vanishing) masses.
  if (id3Tmp == 0 || id3Tmp != id4Tmp) {
    pME[2] = Vec4(  pAbsME * sThe, 0.,  pAbsME * cThe,
             0.5 * (sH + s3ME - s4ME) / mH);
    pME[3] = Vec4( -pAbsME * sThe, 0., -pAbsME * cThe,
             0.5 * (sH + s4ME - s3ME) / mH);

  // For equal (anti)particles (e.g. W+ W-) use averaged mass.
  } else {
    mME[2] = sqrtpos(0.5 * (s3ME + s4ME) - 0.25 * pow2(s3ME - s4ME) / sH);
    mME[3] = mME[2];
    pME[2] = Vec4(  pAbsME * sThe, 0.,  pAbsME * cThe, 0.5 * mH);
    pME[3] = Vec4( -pAbsME * sThe, 0., -pAbsME * cThe, 0.5 * mH);
  }

  // Done.
  return allowME;

}

//==========================================================================

// The Sigma3Process class.
// Base class for resolved 2 -> 3 cross sections; derived from SigmaProcess.

//--------------------------------------------------------------------------

// Input and complement kinematics for resolved 2 -> 3 process.

void Sigma3Process::store3Kin( double x1in, double x2in, double sHin,
  Vec4 p3cmIn, Vec4 p4cmIn, Vec4 p5cmIn, double m3in, double m4in,
  double m5in, double runBW3in, double runBW4in, double runBW5in) {

  // Default ordering of particles 3 and 4 - not relevant here.
  swapTU   = false;

  // Incoming parton momentum fractions.
  x1Save   = x1in;
  x2Save   = x2in;

  // Incoming masses and their squares.
  if (id3Mass() == 0 && id4Mass() == 0 && id5Mass() == 0) {
    m3     = 0.;
    m4     = 0.;
    m5     = 0.;
  } else {
    m3     = m3in;
    m4     = m4in;
    m5     = m5in;
  }
  mSave[3] = m3;
  mSave[4] = m4;
  mSave[5] = m5;
  s3       = m3 * m3;
  s4       = m4 * m4;
  s5       = m5 * m5;

  // Standard Mandelstam variables and four-momenta in rest frame.
  sH       = sHin;
  mH       = sqrt(sH);
  sH2      = sH * sH;
  p3cm     = p3cmIn;
  p4cm     = p4cmIn;
  p5cm     = p5cmIn;

  // The nominal Breit-Wigner factors with running width.
  runBW3   = runBW3in;
  runBW4   = runBW4in;
  runBW5   = runBW5in;

  // Special case: pick scale as if 2 -> 1 process in disguise.
  if (isSChannel()) {

    // Different options for renormalization scale, but normally sHat.
    Q2RenSave = renormMultFac * sH;
    if (renormScale1 == 2) Q2RenSave = renormFixScale;

    // Different options for factorization scale, but normally sHat.
    Q2FacSave = factorMultFac * sH;
    if (factorScale1 == 2) Q2RenSave = factorFixScale;

  // "Normal" 2 -> 3 processes, i.e. not vector boson fusion.
  } else if ( idTchan1() != 23 && idTchan1() != 24 && idTchan2() != 23
    && idTchan2() != 24 ) {
    double mT3S = s3 + p3cm.pT2();
    double mT4S = s4 + p4cm.pT2();
    double mT5S = s5 + p5cm.pT2();

    // Different options for renormalization scale.
    if      (renormScale3 == 1) Q2RenSave = min( mT3S, min(mT4S, mT5S) );
    else if (renormScale3 == 2) Q2RenSave = sqrt( mT3S * mT4S * mT5S
      / max( mT3S, max(mT4S, mT5S) ) );
    else if (renormScale3 == 3) Q2RenSave = pow( mT3S * mT4S * mT5S,
                                            1./3. );
    else if (renormScale3 == 4) Q2RenSave = (mT3S + mT4S + mT5S) / 3.;
    else                        Q2RenSave = sH;
    Q2RenSave                            *= renormMultFac;
    if      (renormScale3 == 6) Q2RenSave = renormFixScale;

    // Different options for factorization scale.
    if      (factorScale3 == 1) Q2FacSave = min( mT3S, min(mT4S, mT5S) );
    else if (factorScale3 == 2) Q2FacSave = sqrt( mT3S * mT4S * mT5S
      / max( mT3S, max(mT4S, mT5S) ) );
    else if (factorScale3 == 3) Q2FacSave = pow( mT3S * mT4S * mT5S,
                                            1./3. );
    else if (factorScale3 == 4) Q2FacSave = (mT3S + mT4S + mT5S) / 3.;
    else                        Q2FacSave = sH;
    Q2FacSave                            *= factorMultFac;
    if      (factorScale3 == 6) Q2FacSave = factorFixScale;

  // Vector boson fusion 2 -> 3 processes; recoils in positions 4 and 5.
  } else {
    double sV4   = pow2( particleDataPtr->m0(idTchan1()) );
    double sV5   = pow2( particleDataPtr->m0(idTchan2()) );
    double mT3S  = s3  + p3cm.pT2();
    double mTV4S = sV4 + p4cm.pT2();
    double mTV5S = sV5 + p5cm.pT2();

    // Different options for renormalization scale.
    if      (renormScale3VV == 1) Q2RenSave = max( sV4, sV5);
    else if (renormScale3VV == 2) Q2RenSave = sqrt( mTV4S * mTV5S );
    else if (renormScale3VV == 3) Q2RenSave = pow( mT3S * mTV4S * mTV5S,
                                              1./3. );
    else if (renormScale3VV == 4) Q2RenSave = (mT3S * mTV4S * mTV5S) / 3.;
    else                          Q2RenSave = sH;
    Q2RenSave                              *= renormMultFac;
    if      (renormScale3VV == 6) Q2RenSave = renormFixScale;

    // Different options for factorization scale.
    if      (factorScale3VV == 1) Q2FacSave = max( sV4, sV5);
    else if (factorScale3VV == 2) Q2FacSave = sqrt( mTV4S * mTV5S );
    else if (factorScale3VV == 3) Q2FacSave = pow( mT3S * mTV4S * mTV5S,
                                              1./3. );
    else if (factorScale3VV == 4) Q2FacSave = (mT3S * mTV4S * mTV5S) / 3.;
    else                          Q2FacSave = sH;
    Q2FacSave                              *= factorMultFac;
    if      (factorScale3VV == 6) Q2FacSave = factorFixScale;
  }

  // Evaluate alpha_strong and alpha_EM.
  alpS  = couplingsPtr->alphaS(Q2RenSave);
  alpEM = couplingsPtr->alphaEM(Q2RenSave);

}

//--------------------------------------------------------------------------

// Calculate modified masses and four-vectors for matrix elements.

bool Sigma3Process::setupForME() {

  // Common initial-state handling.
  bool allowME = setupForMEin();

  // Correct outgoing c, b, mu and tau to be massive or not.
  mME[2] = m3;
  int id3Tmp = abs(id3Mass());
  if (id3Tmp ==  4) mME[2] = mcME;
  if (id3Tmp ==  5) mME[2] = mbME;
  if (id3Tmp == 13) mME[2] = mmuME;
  if (id3Tmp == 15) mME[2] = mtauME;
  mME[3] = m4;
  int id4Tmp = abs(id4Mass());
  if (id4Tmp ==  4) mME[3] = mcME;
  if (id4Tmp ==  5) mME[3] = mbME;
  if (id4Tmp == 13) mME[3] = mmuME;
  if (id4Tmp == 15) mME[3] = mtauME;
  mME[4] = m5;
  int id5Tmp = abs(id5Mass());
  if (id5Tmp ==  4) mME[4] = mcME;
  if (id5Tmp ==  5) mME[4] = mbME;
  if (id5Tmp == 13) mME[4] = mmuME;
  if (id5Tmp == 15) mME[4] = mtauME;

  // If kinematically impossible turn to massless case, but set error.
  if (mME[2] + mME[3] + mME[4] >= mH) {
    mME[2] = 0.;
    mME[3] = 0.;
    mME[4] = 0.;
    allowME = false;
  }

  // Form new average masses if identical particles.
  if (id3Tmp != 0 && id4Tmp == id3Tmp && id5Tmp == id3Tmp) {
    double mAvg = (mME[2] + mME[3] + mME[4]) / 3.;
    mME[2] = mAvg;
    mME[3] = mAvg;
    mME[4] = mAvg;
  } else if (id3Tmp != 0 && id4Tmp == id3Tmp) {
    mME[2] = sqrtpos(0.5 * (pow2(mME[2]) + pow2(mME[3]))
           - 0.25 * pow2(pow2(mME[2]) - pow2(mME[3])) / sH);
    mME[3] = mME[2];
  } else if (id3Tmp != 0 && id5Tmp == id3Tmp) {
    mME[2] = sqrtpos(0.5 * (pow2(mME[2]) + pow2(mME[4]))
           - 0.25 * pow2(pow2(mME[2]) - pow2(mME[4])) / sH);
    mME[4] = mME[2];
  } else if (id4Tmp != 0 && id5Tmp == id4Tmp) {
    mME[3] = sqrtpos(0.5 * (pow2(mME[3]) + pow2(mME[4]))
           - 0.25 * pow2(pow2(mME[3]) - pow2(mME[4])) / sH);
    mME[4] = mME[2];
  }

  // Iterate rescaled three-momenta until convergence.
  double m2ME3 = pow2(mME[2]);
  double m2ME4 = pow2(mME[3]);
  double m2ME5 = pow2(mME[4]);
  double p2ME3 = p3cm.pAbs2();
  double p2ME4 = p4cm.pAbs2();
  double p2ME5 = p5cm.pAbs2();
  double p2sum = p2ME3 + p2ME4 + p2ME5;
  double eME3  = sqrt(m2ME3 + p2ME3);
  double eME4  = sqrt(m2ME4 + p2ME4);
  double eME5  = sqrt(m2ME5 + p2ME5);
  double esum  = eME3 + eME4 + eME5;
  double p2rat = p2ME3 / eME3 + p2ME4 / eME4 + p2ME5 / eME5;
  int iStep = 0;
  while ( abs(esum - mH) > COMPRELERR * mH && iStep < NCOMPSTEP ) {
    ++iStep;
    double compFac = 1. + 2. * (mH - esum) / p2rat;
    p2ME3 *= compFac;
    p2ME4 *= compFac;
    p2ME5 *= compFac;
    eME3   = sqrt(m2ME3 + p2ME3);
    eME4   = sqrt(m2ME4 + p2ME4);
    eME5   = sqrt(m2ME5 + p2ME5);
    esum   = eME3 + eME4 + eME5;
    p2rat  = p2ME3 / eME3 + p2ME4 / eME4 + p2ME5 / eME5;
  }

  // If failed convergence set error flag.
  if (abs(esum - mH) > COMPRELERR * mH) allowME = false;

  // Set up accepted kinematics.
  double totFac = sqrt( (p2ME3 + p2ME4 + p2ME5) / p2sum);
  pME[2] = totFac * p3cm;
  pME[2].e( eME3);
  pME[3] = totFac * p4cm;
  pME[3].e( eME4);
  pME[4] = totFac * p5cm;
  pME[4].e( eME5);

  // Done.
  return allowME;

}

//==========================================================================

// The SigmaLHAProcess class.
// Wrapper for Les Houches Accord external input; derived from SigmaProcess.
// Note: arbitrary subdivision into PhaseSpaceLHA and SigmaLHAProcess tasks.

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double SigmaLHAProcess::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Do nothing if decays present already at input.
  if (iResBeg < process.savedSizeValue()) return 1.;

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For Higgs decay hand over to standard routine.
  if (idMother == 25 || idMother == 35 || idMother == 36)
    return weightHiggsDecay( process, iResBeg, iResEnd);

  // For top decay hand over to standard routine.
  if (idMother == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // Else done.
  return 1.;

}

//--------------------------------------------------------------------------

// Set scale, alpha_strong and alpha_EM when not set.

void SigmaLHAProcess::setScale() {

  // If scale has not been set, then to set.
  double scaleLHA = lhaUpPtr->scale();
  if (scaleLHA < 0.) {

    // Final-state partons and their invariant mass.
    vector<int> iFin;
    Vec4 pFinSum;
    for (int i = 3; i < lhaUpPtr->sizePart(); ++i)
    if (lhaUpPtr->mother1(i) == 1) {
      iFin.push_back(i);
      pFinSum += Vec4( lhaUpPtr->px(i), lhaUpPtr->py(i),
        lhaUpPtr->pz(i), lhaUpPtr->e(i) );
    }
    int nFin = iFin.size();
    sH       = pFinSum * pFinSum;
    mH       = sqrt(sH);
    sH2      = sH * sH;

    // If 1 final-state particle then use Sigma1Process logic.
    if (nFin == 1) {
      Q2RenSave                             = renormMultFac * sH;
      if (renormScale1 == 2) Q2RenSave      = renormFixScale;
      Q2FacSave                             = factorMultFac * sH;
      if (factorScale1 == 2) Q2FacSave      = factorFixScale;

    // If 2 final-state particles then use Sigma2Process logic.
    } else if (nFin == 2) {
      double s3  = pow2(lhaUpPtr->m(iFin[0]));
      double s4  = pow2(lhaUpPtr->m(iFin[1]));
      double pT2 = pow2(lhaUpPtr->px(iFin[0])) + pow2(lhaUpPtr->py(iFin[0]));
      if      (renormScale2 == 1) Q2RenSave = pT2 + min(s3, s4);
      else if (renormScale2 == 2) Q2RenSave = sqrt((pT2 + s3) * (pT2 + s4));
      else if (renormScale2 == 3) Q2RenSave = pT2 + 0.5 * (s3 + s4);
      else                        Q2RenSave = sH;
      Q2RenSave                            *= renormMultFac;
      if      (renormScale2 == 5) Q2RenSave = renormFixScale;
      if      (factorScale2 == 1) Q2FacSave = pT2 + min(s3, s4);
      else if (factorScale2 == 2) Q2FacSave = sqrt((pT2 + s3) * (pT2 + s4));
      else if (factorScale2 == 3) Q2FacSave = pT2 + 0.5 * (s3 + s4);
      else                        Q2FacSave = sH;
      Q2FacSave                            *= factorMultFac;
      if      (factorScale2 == 5) Q2FacSave = factorFixScale;

    // If 3 or more final-state particles then use Sigma3Process logic.
    } else {
      double mTSlow  = sH;
      double mTSmed  = sH;
      double mTSprod = 1.;
      double mTSsum  = 0.;
      for (int i = 0; i < nFin; ++i) {
        double mTSnow = pow2(lhaUpPtr->m(iFin[i]))
          + pow2(lhaUpPtr->px(iFin[i])) + pow2(lhaUpPtr->py(iFin[i]));
        if      (mTSnow < mTSlow) {mTSmed = mTSlow; mTSlow = mTSnow;}
        else if (mTSnow < mTSmed) mTSmed = mTSnow;
        mTSprod *= mTSnow;
        mTSsum  += mTSnow;
      }
      if      (renormScale3 == 1) Q2RenSave = mTSlow;
      else if (renormScale3 == 2) Q2RenSave = sqrt(mTSlow * mTSmed);
      else if (renormScale3 == 3) Q2RenSave = pow(mTSprod, 1. / nFin);
      else if (renormScale3 == 4) Q2RenSave = mTSsum / nFin;
      else                        Q2RenSave = sH;
      Q2RenSave                            *= renormMultFac;
      if      (renormScale3 == 6) Q2RenSave = renormFixScale;
      if      (factorScale3 == 1) Q2FacSave = mTSlow;
      else if (factorScale3 == 2) Q2FacSave = sqrt(mTSlow * mTSmed);
      else if (factorScale3 == 3) Q2FacSave = pow(mTSprod, 1. / nFin);
      else if (factorScale3 == 4) Q2FacSave = mTSsum / nFin;
      else                        Q2FacSave = sH;
      Q2FacSave                            *= factorMultFac;
      if      (factorScale3 == 6) Q2FacSave = factorFixScale;
    }
  }

  // If alpha_strong and alpha_EM have not been set, then set them.
  if (lhaUpPtr->alphaQCD() < 0.001) {
    double Q2RenNow = (scaleLHA < 0.) ? Q2RenSave : pow2(scaleLHA);
    alpS = couplingsPtr->alphaS(Q2RenNow);
  }
  if (lhaUpPtr->alphaQED() < 0.001) {
    double Q2RenNow = (scaleLHA < 0.) ? Q2RenSave : pow2(scaleLHA);
    alpEM = couplingsPtr->alphaEM(Q2RenNow);
  }

}

//--------------------------------------------------------------------------

// Obtain number of final-state partons from LHA object.

int SigmaLHAProcess::nFinal() const {

  // At initialization size unknown, so return 0.
  if (lhaUpPtr->sizePart() <= 0) return 0;

  // Sum up all particles that has first mother = 1.
  int nFin = 0;
  for (int i = 3; i < lhaUpPtr->sizePart(); ++i)
    if (lhaUpPtr->mother1(i) == 1) ++nFin;
  return nFin;

}

//==========================================================================

} // end namespace Pythia8
