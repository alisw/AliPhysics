// SigmaHiggs.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// Part of code written by Marc Montull, CERN summer student 2007.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Function definitions (not found in the header) for the
// Higgs simulation classes.

#include "Pythia8/SigmaHiggs.h"

namespace Pythia8 {

//==========================================================================

// Sigma1ffbar2H class.
// Cross section for f fbar -> H0 , H1, H2 or A3.
// (f is quark or lepton, H0 SM Higgs and H1, H2, A3 BSM Higgses ).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ffbar2H::initProc() {

  // Properties specific to Higgs state.
  if (higgsType == 0) {
    nameSave = "f fbar -> H (SM)";
    codeSave = 901;
    idRes    = 25;
  }
  else if (higgsType == 1) {
    nameSave = "f fbar -> h0(H1)";
    codeSave = 1001;
    idRes    = 25;
  }
  else if (higgsType == 2) {
    nameSave = "f fbar -> H0(H2)";
    codeSave = 1021;
    idRes    = 35;
  }
  else if (higgsType == 3) {
    nameSave = "f fbar -> A0(A3)";
    codeSave = 1041;
    idRes    = 36;
  }

  // Find pointer to H0, H1, H2 or A3 depending on the value of idRes.
  HResPtr = particleDataPtr->particleDataEntryPtr(idRes);

  // Store H0, H1, H2 or A3 mass and width for propagator.
  mRes     = HResPtr->m0();
  GammaRes = HResPtr->mWidth();
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

}


//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma1ffbar2H::sigmaKin() {

  // Set up Breit-Wigner.
  double width = HResPtr->resWidth(idRes, mH);
  sigBW        = 4. * M_PI/ ( pow2(sH - m2Res) + pow2(mH * width) );

  // Width out only includes open channels.
  widthOut     = width * HResPtr->resOpenFrac(idRes);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma1ffbar2H::sigmaHat() {

  // Calculate mass-dependent incoming width, including colour factor.
  int idAbs      = abs(id1);
  double widthIn = HResPtr->resWidthChan( mH, idAbs, -idAbs);
  if (idAbs < 9) widthIn /= 9.;

  // Done.
  return widthIn * sigBW * widthOut;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2H::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idRes);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma1ffbar2H::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma1gg2H class.
// Cross section for g g -> H0, H1, H2 or A3 (H0 SM Higgs, H1, H2, A3 BSM).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1gg2H::initProc() {

  // Properties specific to Higgs state.
  if (higgsType == 0) {
    nameSave = "g g -> H (SM)";
    codeSave = 902;
    idRes    = 25;
  }
  else if (higgsType == 1) {
    nameSave = "g g -> h0(H1)";
    codeSave = 1002;
    idRes    = 25;
  }
  else if (higgsType == 2) {
    nameSave = "g g -> H0(H2)";
    codeSave = 1022;
    idRes    = 35;
  }
  else if (higgsType == 3) {
    nameSave = "g g -> A0(A3)";
    codeSave = 1042;
    idRes    = 36;
  }

  // Find pointer to H0, H1, H2 or A3 depending on idRes.
  HResPtr = particleDataPtr->particleDataEntryPtr(idRes);

  // Store H0, H1, H2 or A3 mass and width for propagator.
  mRes     = HResPtr->m0();
  GammaRes = HResPtr->mWidth();
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1gg2H::sigmaKin() {

  // Incoming width for gluons, gives colour factor of 1/8 * 1/8.
  double widthIn  = HResPtr->resWidthChan( mH, 21, 21) / 64.;

  // Set up Breit-Wigner.
  double width    = HResPtr->resWidth(idRes, mH);
  double sigBW    = 8. * M_PI/ ( pow2(sH - m2Res) + pow2(mH * width) );

  // Width out only includes open channels.
  double widthOut = width * HResPtr->resOpenFrac(idRes);

  // Done.
  sigma = widthIn * sigBW * widthOut;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1gg2H::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, idRes);

  // Colour flow topology.
  setColAcol( 1, 2, 2, 1, 0, 0);

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma1gg2H::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma1gmgm2H class.
// Cross section for gamma gamma -> H0, H1, H2 or H3.
// (H0 SM Higgs, H1, H2 and A3 BSM Higgses).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1gmgm2H::initProc() {

  // Properties specific to Higgs state.
  if (higgsType == 0) {
    nameSave = "gamma gamma -> H (SM)";
    codeSave = 903;
    idRes    = 25;
  }
  else if (higgsType == 1) {
    nameSave = "gamma gamma -> h0(H1)";
    codeSave = 1003;
    idRes    = 25;
  }
  else if (higgsType == 2) {
    nameSave = "gamma gamma -> H0(H2)";
    codeSave = 1023;
    idRes    = 35;
  }
  else if (higgsType == 3) {
    nameSave = "gamma gamma -> A0(A3)";
    codeSave = 1043;
    idRes    = 36;
  }

  // Find pointer to H0, H1, H2 or A3.
  HResPtr = particleDataPtr->particleDataEntryPtr(idRes);

  // Store H0, H1, H2 or A3 mass and width for propagator.
  mRes     = HResPtr->m0();
  GammaRes = HResPtr->mWidth();
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1gmgm2H::sigmaKin() {

  // Incoming width for photons.
  double widthIn  = HResPtr->resWidthChan( mH, 22, 22);

  // Set up Breit-Wigner.
  double width    = HResPtr->resWidth(idRes, mH);
  double sigBW    = 8. * M_PI/ ( pow2(sH - m2Res) + pow2(mH * width) );

  // Width out only includes open channels.
  double widthOut = width * HResPtr->resOpenFrac(idRes);

  // Done.
  sigma = widthIn * sigBW * widthOut;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1gmgm2H::setIdColAcol() {

  // Flavours trivial.
  setId( 22, 22, idRes);

  // Colour flow trivial.
  setColAcol( 0, 0, 0, 0, 0, 0);

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma1gmgm2H::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma2ffbar2HZ class.
// Cross section for f fbar -> H0 Z0, H1 Z0, H2 Z0 or A3 Z0.
// (H0 SM Higgs, H1, H2 and A3 BSM Higgses).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2HZ::initProc() {

  // Properties specific to Higgs state.
  if (higgsType == 0) {
     nameSave = "f fbar -> H0 Z0 (SM)";
     codeSave = 904;
     idRes    = 25;
     coup2Z   = 1.;
  }
  else if (higgsType == 1) {
    nameSave = "f fbar -> h0(H1) Z0";
    codeSave = 1004;
    idRes    = 25;
    coup2Z   = settingsPtr->parm("HiggsH1:coup2Z");
  }
  else if (higgsType == 2) {
    nameSave = "f fbar -> H0(H2) Z0";
    codeSave = 1024;
    idRes    = 35;
    coup2Z   = settingsPtr->parm("HiggsH2:coup2Z");
  }
  else if (higgsType == 3) {
    nameSave = "f fbar -> A0(A3) ZO";
    codeSave = 1044;
    idRes    = 36;
    coup2Z   = settingsPtr->parm("HiggsA3:coup2Z");
  }

  // Store Z0 mass and width for propagator. Common coupling factor.
  mZ           = particleDataPtr->m0(23);
  widZ         = particleDataPtr->mWidth(23);
  mZS          = mZ*mZ;
  mwZS         = pow2(mZ * widZ);
  thetaWRat    = 1. / (16. * couplingsPtr->sin2thetaW()
                 * couplingsPtr->cos2thetaW());

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(idRes, 23);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2HZ::sigmaKin() {

  // Evaluate differential cross section.
  sigma0 = (M_PI / sH2) * 8. * pow2(alpEM * thetaWRat * coup2Z)
    * (tH * uH - s3 * s4 + 2. * sH * s4) / (pow2(sH - mZS) + mwZS);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2HZ::sigmaHat() {

  // Coupling a_f^2 + v_f^2 to s-channel Z0 and colour factor.
  int idAbs    = abs(id1);
  double sigma = sigma0 * couplingsPtr->vf2af2(idAbs);
  if (idAbs < 9) sigma /= 3.;

  // Secondary width for H0 and Z0 or H1 and Z0 or H2 and Z0 or A3 and Z0.
  sigma       *= openFracPair;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2HZ::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idRes, 23);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2ffbar2HZ::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For Higgs decay hand over to standard routine.
  if (idMother == 25 || idMother == 35 || idMother == 36)
    return weightHiggsDecay( process, iResBeg, iResEnd);

  // For top decay hand over to standard routine.
  if (idMother == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // If not decay of Z0 created along with Higgs then done.
  if (iResBeg != 5 || iResEnd != 6) return 1.;

  // Order so that fbar(1) f(2) -> H() f'(3) fbar'(4).
  int i1       = (process[3].id() < 0) ? 3 : 4;
  int i2       = 7 - i1;
  int i3       = process[6].daughter1();
  int i4       = process[6].daughter2();
  if (process[i3].id() < 0) swap( i3, i4);

  // Find left- and righthanded couplings of fermion pairs.
  int    idAbs = process[i1].idAbs();
  double liS   = pow2( couplingsPtr->lf(idAbs) );
  double riS   = pow2( couplingsPtr->rf(idAbs) );
  idAbs        = process[i3].idAbs();
  double lfS   = pow2( couplingsPtr->lf(idAbs) );
  double rfS   = pow2( couplingsPtr->rf(idAbs) );

  // Evaluate relevant four-products.
  double pp13  = process[i1].p() * process[i3].p();
  double pp14  = process[i1].p() * process[i4].p();
  double pp23  = process[i2].p() * process[i3].p();
  double pp24  = process[i2].p() * process[i4].p();

  // Weight and maximum.
  double wt    = (liS * lfS + riS * rfS) * pp13 * pp24
               + (liS * rfS + riS * lfS) * pp14 * pp23;
  double wtMax = (liS + riS) * (lfS + rfS) * (pp13 + pp14) * (pp23 + pp24);

  // Done.
  return wt / wtMax;

}

//==========================================================================

// Sigma2ffbar2HW class.
// Cross section for f fbar -> H0 W+-, H1 W+-, H2 W+- or A3 W+-.
// (H0 SM Higgs, H1, H2 and A3 BSM Higgses).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2HW::initProc() {

  // Properties specific to Higgs state.
  if (higgsType == 0) {
    nameSave = "f fbar -> H0 W+- (SM)";
    codeSave = 905;
    idRes    = 25;
    coup2W   = 1.;
  }
  else if (higgsType == 1) {
    nameSave = "f fbar -> h0(H1) W+-";
    codeSave = 1005;
    idRes    = 25;
    coup2W   = settingsPtr->parm("HiggsH1:coup2W");
  }
  else if (higgsType == 2) {
    nameSave = "f fbar -> H0(H2) W+-";
    codeSave = 1025;
    idRes    = 35;
    coup2W   = settingsPtr->parm("HiggsH2:coup2W");
  }
  else if (higgsType == 3) {
    nameSave = "f fbar -> A0(A3) W+-";
    codeSave = 1045;
    idRes    = 36;
    coup2W   = settingsPtr->parm("HiggsA3:coup2W");
  }

  // Store W+- mass and width for propagator. Common coupling factor.
  mW              = particleDataPtr->m0(24);
  widW            = particleDataPtr->mWidth(24);
  mWS             = mW*mW;
  mwWS            = pow2(mW * widW);
  thetaWRat       = 1. / (4. * couplingsPtr->sin2thetaW());

  // Secondary open width fractions.
  openFracPairPos = particleDataPtr->resOpenFrac(idRes,  24);
  openFracPairNeg = particleDataPtr->resOpenFrac(idRes, -24);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2HW::sigmaKin() {

  //  Evaluate differential cross section.
  sigma0 = (M_PI / sH2) * 2. * pow2(alpEM * thetaWRat * coup2W)
    * (tH * uH - s3 * s4 + 2. * sH * s4) / (pow2(sH - mWS) + mwWS);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2HW::sigmaHat() {

  // CKM and colour factors.
  double sigma = sigma0;
  if (abs(id1) < 9) sigma *= couplingsPtr->V2CKMid(abs(id1), abs(id2)) / 3.;

  // Secondary width for H0 and W+-.
  int idUp     = (abs(id1)%2 == 0) ? id1 : id2;
  sigma       *= (idUp > 0) ? openFracPairPos : openFracPairNeg;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2HW::setIdColAcol() {

  // Sign of outgoing W.
  int sign = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, idRes, 24 * sign);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2ffbar2HW::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For Higgs decay hand over to standard routine.
  if (idMother == 25 || idMother == 35 || idMother == 36)
    return weightHiggsDecay( process, iResBeg, iResEnd);

  // For top decay hand over to standard routine.
  if (idMother == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // If not decay of W+- created along with Higgs then done.
  if (iResBeg != 5 || iResEnd != 6) return 1.;

  // Order so that fbar(1) f(2) -> H() f'(3) fbar'(4).
  int i1       = (process[3].id() < 0) ? 3 : 4;
  int i2       = 7 - i1;
  int i3       = process[6].daughter1();
  int i4       = process[6].daughter2();
  if (process[i3].id() < 0) swap( i3, i4);

  // Evaluate relevant four-products.
  double pp13  = process[i1].p() * process[i3].p();
  double pp14  = process[i1].p() * process[i4].p();
  double pp23  = process[i2].p() * process[i3].p();
  double pp24  = process[i2].p() * process[i4].p();

  // Weight and maximum.
  double wt    = pp13 * pp24;
  double wtMax = (pp13 + pp14) * (pp23 + pp24);

  // Done.
  return wt / wtMax;

}

//==========================================================================

// Sigma3ff2HfftZZ class.
// Cross section for f f' -> H f f' (Z0 Z0 fusion of SM or BSM Higgs).
// (H can be H0 SM or H1, H2, A3 from BSM).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma3ff2HfftZZ::initProc() {

  // Properties specific to Higgs state.
  if (higgsType == 0) {
    nameSave = "f f' -> H0 f f'(Z0 Z0 fusion) (SM)";
    codeSave = 906;
    idRes    = 25;
    coup2Z   = 1.;
  }
  else if (higgsType == 1) {
    nameSave = "f f' -> h0(H1) f f' (Z0 Z0 fusion)";
    codeSave = 1006;
    idRes    = 25;
    coup2Z  = settingsPtr->parm("HiggsH1:coup2Z");
  }
  else if (higgsType == 2) {
    nameSave = "f f' -> H0(H2) f f' (Z0 Z0 fusion)";
    codeSave = 1026;
    idRes    = 35;
    coup2Z  = settingsPtr->parm("HiggsH2:coup2Z");
  }
  else if (higgsType == 3) {
    nameSave = "f f' -> A0(A3) f f' (Z0 Z0 fusion)";
    codeSave = 1046;
    idRes    = 36;
    coup2Z  = settingsPtr->parm("HiggsA3:coup2Z");
  }

  // Common fixed mass and coupling factor.
  mZS = pow2( particleDataPtr->m0(23) );
  prefac = 0.25 * mZS * pow3( 4. * M_PI / (couplingsPtr->sin2thetaW()
           * couplingsPtr->cos2thetaW()) );

  // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(idRes);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma3ff2HfftZZ::sigmaKin() {

  // Required four-vector products.
  double pp12 = 0.5 * sH;
  double pp14 = 0.5 * mH * p4cm.pNeg();
  double pp15 = 0.5 * mH * p5cm.pNeg();
  double pp24 = 0.5 * mH * p4cm.pPos();
  double pp25 = 0.5 * mH * p5cm.pPos();
  double pp45 = p4cm * p5cm;

  // Propagator factors and two possible numerators.
  double prop = pow2( (2. * pp14 + mZS) * (2. * pp25 + mZS) );
  sigma1      = prefac * pp12 * pp45  / prop;
  sigma2      = prefac * pp15 * pp24  / prop;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma3ff2HfftZZ::sigmaHat() {

  // Flavour-dependent coupling factors for two incoming flavours.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  double lf1S  = pow2( couplingsPtr->lf(id1Abs) );
  double rf1S  = pow2( couplingsPtr->rf(id1Abs) );
  double lf2S  = pow2( couplingsPtr->lf(id2Abs) );
  double rf2S  = pow2( couplingsPtr->rf(id2Abs) );
  double c1    = lf1S * lf2S + rf1S * rf2S;
  double c2    = lf1S * rf2S + rf1S * lf2S;

  // Combine couplings and kinematics factors.
  double sigma = pow3(alpEM) * (c1 * sigma1 + c2 * sigma2) * pow2(coup2Z);

  // Secondary width for H0, H1, H2 or A3.
  sigma       *= openFrac;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3ff2HfftZZ::setIdColAcol() {

  // Trivial flavours: out = in.
  setId( id1, id2, idRes, id1, id2);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9 && abs(id2) < 9 && id1*id2 > 0)
                         setColAcol( 1, 0, 2, 0, 0, 0, 1, 0, 2, 0);
  else if (abs(id1) < 9 && abs(id2) < 9)
                         setColAcol( 1, 0, 0, 2, 0, 0, 1, 0, 0, 2);
  else if (abs(id1) < 9) setColAcol( 1, 0, 0, 0, 0, 0, 1, 0, 0, 0);
  else if (abs(id2) < 9) setColAcol( 0, 0, 1, 0, 0, 0, 0, 0, 1, 0);
  else                   setColAcol( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  if ( (abs(id1) < 9 && id1 < 0) || (abs(id1) > 10 && id2 < 0) )
    swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma3ff2HfftZZ::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma3ff2HfftWW class.
// Cross section for f_1 f_2 -> H0 f_3 f_4 (W+ W- fusion of SM or BSM Higgs).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma3ff2HfftWW::initProc() {

  // Properties specific to Higgs state.
  if (higgsType == 0) {
    nameSave = "f_1 f_2 -> H0 f_3 f_4 (W+ W- fusion) (SM)";
    codeSave = 907;
    idRes    = 25;
    coup2W   = 1.;
  }
  else if (higgsType == 1) {
    nameSave = "f_1 f_2 -> h0(H1) f_3 f_4 (W+ W- fusion)";
    codeSave = 1007;
    idRes    = 25;
    coup2W   = settingsPtr->parm("HiggsH1:coup2W");
  }
  else if (higgsType == 2) {
    nameSave = "f_1 f_2 -> H0(H2) f_3 f_4 (W+ W- fusion)";
    codeSave = 1027;
    idRes    = 35;
    coup2W   = settingsPtr->parm("HiggsH2:coup2W");
  }
  else if (higgsType == 3) {
    nameSave = "f_1 f_2 -> A0(A3) f_3 f_4 (W+ W- fusion)";
    codeSave = 1047;
    idRes    = 36;
    coup2W   = settingsPtr->parm("HiggsA3:coup2W");
  }

  // Common fixed mass and coupling factor.
  mWS = pow2( particleDataPtr->m0(24) );
  prefac = mWS * pow3( 4. * M_PI / couplingsPtr->sin2thetaW() );

  // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(idRes);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma3ff2HfftWW::sigmaKin() {

  // Required four-vector products.
  double pp12  = 0.5 * sH;
  double pp14  = 0.5 * mH * p4cm.pNeg();
  double pp25  = 0.5 * mH * p5cm.pPos();
  double pp45  = p4cm * p5cm;

  // Cross section: kinematics part. Combine with couplings.
  double prop  = pow2( (2. * pp14 + mWS) * (2. * pp25 + mWS) );
  sigma0       = prefac * pp12 * pp45 * pow2(coup2W) / prop;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma3ff2HfftWW::sigmaHat() {

  // Some flavour combinations not possible.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  if ( (id1Abs%2 == id2Abs%2 && id1 * id2 > 0)
    || (id1Abs%2 != id2Abs%2 && id1 * id2 < 0) ) return 0.;

  // Basic cross section. CKM factors for final states.
  double sigma = sigma0 * pow3(alpEM) * couplingsPtr->V2CKMsum(id1Abs)
    * couplingsPtr->V2CKMsum(id2Abs);

  // Secondary width for H0, H1, H2 or A3.
  sigma       *= openFrac;

  // Spin-state extra factor 2 per incoming neutrino.
  if (id1Abs == 12 || id1Abs == 14 || id1Abs == 16) sigma *= 2.;
  if (id2Abs == 12 || id2Abs == 14 || id2Abs == 16) sigma *= 2.;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3ff2HfftWW::setIdColAcol() {

  // Pick out-flavours by relative CKM weights.
  id4 = couplingsPtr->V2CKMpick(id1);
  id5 = couplingsPtr->V2CKMpick(id2);
  setId( id1, id2, idRes, id4, id5);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9 && abs(id2) < 9 && id1*id2 > 0)
                         setColAcol( 1, 0, 2, 0, 0, 0, 1, 0, 2, 0);
  else if (abs(id1) < 9 && abs(id2) < 9)
                         setColAcol( 1, 0, 0, 2, 0, 0, 1, 0, 0, 2);
  else if (abs(id1) < 9) setColAcol( 1, 0, 0, 0, 0, 0, 1, 0, 0, 0);
  else if (abs(id2) < 9) setColAcol( 0, 0, 1, 0, 0, 0, 0, 0, 1, 0);
  else                   setColAcol( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  if ( (abs(id1) < 9 && id1 < 0) || (abs(id1) > 10 && id2 < 0) )
    swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma3ff2HfftWW::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma3gg2HQQbar class.
// Cross section for g g -> H0 Q Qbar (Q Qbar fusion of SM or BSM Higgs).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma3gg2HQQbar::initProc() {

  // Properties specific to Higgs state for the "g g -> H ttbar" process.
  // (H can be H0 SM or H1, H2, A3 from BSM).
  if (higgsType == 0 && idNew == 6) {
    nameSave = "g g -> H t tbar (SM)";
    codeSave = 908;
    idRes    = 25;
    coup2Q   = 1.;
  }
  else if (higgsType == 1 && idNew == 6) {
    nameSave = "g g -> h0(H1) t tbar";
    codeSave = 1008;
    idRes    = 25;
    coup2Q   = settingsPtr->parm("HiggsH1:coup2u");
  }
  else if (higgsType == 2 && idNew == 6) {
    nameSave = "g g -> H0(H2) t tbar";
    codeSave = 1028;
    idRes    = 35;
    coup2Q   = settingsPtr->parm("HiggsH2:coup2u");
  }
  else if (higgsType == 3 && idNew == 6) {
    nameSave = "g g -> A0(A3) t tbar";
    codeSave = 1048;
    idRes    = 36;
    coup2Q   = settingsPtr->parm("HiggsA3:coup2u");
  }

  // Properties specific to Higgs state for the "g g -> H b bbar" process.
  // (H can be H0 SM or H1, H2, A3 from BSM).
  if (higgsType == 0 && idNew == 5) {
    nameSave = "g g -> H b bbar (SM)";
    codeSave = 912;
    idRes    = 25;
    coup2Q   = 1.;
  }
  else if (higgsType == 1 && idNew == 5) {
    nameSave = "g g -> h0(H1) b bbar";
    codeSave = 1012;
    idRes    = 25;
    coup2Q   = settingsPtr->parm("HiggsH1:coup2d");
  }
  else if (higgsType == 2 && idNew == 5) {
    nameSave = "g g -> H0(H2) b bbar";
    codeSave = 1032;
    idRes    = 35;
    coup2Q   = settingsPtr->parm("HiggsH2:coup2d");
  }
  else if (higgsType == 3 && idNew == 5) {
    nameSave = "g g -> A0(A3) b bbar";
    codeSave = 1052;
    idRes    = 36;
    coup2Q   = settingsPtr->parm("HiggsA3:coup2d");
  }

  // Common mass and coupling factors.
  double mWS      = pow2(particleDataPtr->m0(24));
  prefac          = (4. * M_PI / couplingsPtr->sin2thetaW()) * pow2(4. * M_PI)
                  * 0.25 / mWS;

  // Secondary open width fraction.
  openFracTriplet = particleDataPtr->resOpenFrac(idRes, idNew, -idNew);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma3gg2HQQbar::sigmaKin() {

  // Running mass of heavy quark.
  double mQ2run = pow2( particleDataPtr->mRun(idNew, mH) );

  // Linear combination of p_Q and p_Qbar to ensure common mass.
  double mQ2  = m4 * m5;
  double epsi = 0.;
  if (m4 != m5) {
    double s45  = (p4cm + p5cm).m2Calc();
    mQ2       = 0.5 * (s4 + s5) - 0.25 * pow2(s4 - s5) / s45;
    epsi      = 0.5 * (s5 - s4) / s45;
  }

  // Set up kinematics: g(4) g(5) -> H(3) Q(1) Qbar(2) in outgoing sense.
  Vec4 pTemp[6];
  pTemp[4]   = Vec4( 0., 0., -0.5* mH, -0.5* mH);
  pTemp[5]   = Vec4( 0., 0.,  0.5* mH, -0.5* mH);
  pTemp[1]   = p4cm + epsi * (p4cm + p5cm);
  pTemp[2]   = p5cm - epsi * (p4cm + p5cm);
  pTemp[3]   = p3cm;

  // Four-product combinations.
  double z1  = pTemp[1] * pTemp[2];
  double z2  = pTemp[1] * pTemp[3];
  double z3  = pTemp[1] * pTemp[4];
  double z4  = pTemp[1] * pTemp[5];
  double z5  = pTemp[2] * pTemp[3];
  double z6  = pTemp[2] * pTemp[4];
  double z7  = pTemp[2] * pTemp[5];
  double z8  = pTemp[3] * pTemp[4];
  double z9  = pTemp[3] * pTemp[5];
  double z10 = pTemp[4] * pTemp[5];

  // Powers required as shorthand in matriz elements.
  double mQ4  = mQ2 * mQ2;
  double mQ6  = mQ2 * mQ4;
  double z1S  = z1 * z1;
  double z2S  = z2 * z2;
  double z3S  = z3 * z3;
  double z4S  = z4 * z4;
  double z5S  = z5 * z5;
  double z6S  = z6 * z6;
  double z7S  = z7 * z7;
  double z8S  = z8 * z8;
  double z9S  = z9 * z9;
  double z10S = z10 * z10;

  // Evaluate matriz elements for g + g -> Q + Qbar + H.
  // (H can be H0 SM or H1, H2, A3 from BSM).
  double fm[9][9];
  fm[1][1] = 64*mQ6+16*mQ4*s3+32*mQ4*(z1+2*z2+z4+z9+2*
  z7+z5)+8*mQ2*s3*(-z1-z4+2*z7)+16*mQ2*(z2*z9+4*z2*
  z7+z2*z5-2*z4*z7-2*z9*z7)+8*s3*z4*z7-16*z2*z9*z7;
  fm[1][2] = 16*mQ6+8*mQ4*(-2*z1+z2-2*z3-2*z4-4*z10+z9-z8+2
  *z7-4*z6+z5)+8*mQ2*(-2*z1*z2-2*z2*z4-2*z2*z10+z2*z7-2*
  z2*z6-2*z3*z7+2*z4*z7+4*z10*z7-z9*z7-z8*z7)+16*z2*z7*(z4+
  z10);
  fm[1][3] = 16*mQ6-4*mQ4*s3+8*mQ4*(-2*z1+2*z2-2*z3-4*
  z4-8*z10+z9+z8-2*z7-4*z6+2*z5)-(4*mQ2*s3)*(z1+z4+z10
  +z6)+8*mQ2*(-2*z1*z2-2*z1*z10+z1*z9+z1*z8-2*z1*z5+z2S
  -4*z2*z4-5*z2*z10+z2*z8-z2*z7-3*z2*z6+z2*z5+z3*z9+2*z3*z7
  -z3*z5+z4*z8+2*z4*z6-3*z4*z5-5*z10*z5+z9*z8+z9*z6+z9*z5+
  z8*z7-4*z6*z5+z5S)-(16*z2*z5)*(z1+z4+z10+z6);
  fm[1][4] = 16*mQ6+4*mQ4*s3+16*mQ4*(-z1+z2-z3-z4+z10-
  z9-z8+2*z7+2*z6-z5)+4*mQ2*s3*(z1+z3+z4+z10+2*z7+2*z6
  )+8*mQ2*(4*z1*z10+4*z1*z7+4*z1*z6+2*z2*z10-z2*z9-z2*z8+
  4*z2*z7+4*z2*z6-z2*z5+4*z10*z5+4*z7*z5+4*z6*z5)-(8*s3*
  z1)*(z10+z7+z6)+16*z2*z5*(z10+z7+z6);
  fm[1][5] = 8*mQ4*(-2*z1-2*z4+z10-z9)+4*mQ2*(4*z1S-2*z1*
  z2+8*z1*z3+6*z1*z10-2*z1*z9+4*z1*z8+4*z1*z7+4*z1*z6+2*z1*
  z5+z2*z10+4*z3*z4-z3*z9+2*z3*z7+3*z4*z8-2*z4*z6+2*z4*z5-4
  *z10*z7+3*z10*z5-3*z9*z6+3*z8*z7-4*z7S+4*z7*z5)+8*(z1S
  *z9-z1S*z8-z1*z2*z7+z1*z2*z6+z1*z3*z9+z1*z3*z5-z1*z4*
  z8-z1*z4*z5+z1*z10*z9+z1*z9*z7+z1*z9*z6-z1*z8*z7-z2*z3*z7
  +z2*z4*z6-z2*z10*z7-z2*z7S+z3*z7*z5-z4*z10*z5-z4*z7*z5-
  z4*z6*z5);
  fm[1][6] = 16*mQ4*(-4*z1-z4+z9-z7)+4*mQ2*s3*(-2*z1-z4-
  z7)+16*mQ2*(-2*z1S-3*z1*z2-2*z1*z4-3*z1*z9-2*z1*z7-3*
  z1*z5-2*z2*z4-2*z7*z5)-8*s3*z4*z7+8*(-z1*z2*z9-2*z1*z2
  *z5-z1*z9S-z1*z9*z5+z2S*z7-z2*z4*z5+z2*z9*z7-z2*z7*z5
  +z4*z9*z5+z4*z5S);
  fm[1][7] = 8*mQ4*(2*z3+z4+3*z10+z9+2*z8+3*z7+6*z6)+2*mQ2*
  s3*(-2*z3-z4+3*z10+3*z7+6*z6)+4*mQ2*(4*z1*z10+4*z1*
  z7+8*z1*z6+6*z2*z10+z2*z9+2*z2*z8+6*z2*z7+12*z2*z6-8*z3*
  z7+4*z4*z7+4*z4*z6+4*z10*z5+4*z9*z7+4*z9*z6-8*z8*z7+4*z7*
  z5+8*z6*z5)+4*s3*(-z1*z10-z1*z7-2*z1*z6+2*z3*z7-z4*z7-
  z4*z6)+8*z2*(z10*z5+z9*z7+z9*z6-2*z8*z7+z7*z5+2*z6*z5);
  fm[1][8] = 8*mQ4*(2*z3+z4+3*z10+2*z9+z8+3*z7+6*z6)+2*mQ2*
  s3*(-2*z3-z4+2*z10+z7+2*z6)+4*mQ2*(4*z1*z10-2*z1*z9+
  2*z1*z8+4*z1*z7+8*z1*z6+5*z2*z10+2*z2*z9+z2*z8+4*z2*z7+8*
  z2*z6-z3*z9-8*z3*z7+2*z3*z5+2*z4*z9-z4*z8+4*z4*z7+4*z4*z6
  +4*z4*z5+5*z10*z5+z9S-z9*z8+2*z9*z7+5*z9*z6+z9*z5-7*z8*
  z7+2*z8*z5+2*z7*z5+10*z6*z5)+2*s3*(-z1*z10+z3*z7-2*z4*
  z7+z4*z6)+4*(-z1*z9S+z1*z9*z8-2*z1*z9*z5-z1*z8*z5+2*z2*
  z10*z5+z2*z9*z7+z2*z9*z6-2*z2*z8*z7+3*z2*z6*z5+z3*z9*z5+
  z3*z5S+z4*z9*z5-2*z4*z8*z5+2*z4*z5S);
  fm[2][2] = 16*mQ6+16*mQ4*(-z1+z3-z4-z10+z7-z6)+16*mQ2*(
  z3*z10+z3*z7+z3*z6+z4*z7+z10*z7)-16*z3*z10*z7;
  fm[2][3] = 16*mQ6+8*mQ4*(-2*z1+z2+2*z3-4*z4-4*z10-z9+z8-2
  *z7-2*z6+z5)+8*mQ2*(-2*z1*z5+4*z3*z10-z3*z9-z3*z8-2*z3*
  z7+2*z3*z6+z3*z5-2*z4*z5-2*z10*z5-2*z6*z5)+16*z3*z5*(z10+
  z6);
  fm[2][4] = 8*mQ4*(-2*z1-2*z3+z10-z8)+4*mQ2*(4*z1S-2*z1*
  z2+8*z1*z4+6*z1*z10+4*z1*z9-2*z1*z8+4*z1*z7+4*z1*z6+2*z1*
  z5+z2*z10+4*z3*z4+3*z3*z9-2*z3*z7+2*z3*z5-z4*z8+2*z4*z6-4
  *z10*z6+3*z10*z5+3*z9*z6-3*z8*z7-4*z6S+4*z6*z5)+8*(-z1S
  *z9+z1S*z8+z1*z2*z7-z1*z2*z6-z1*z3*z9-z1*z3*z5+z1*z4
  *z8+z1*z4*z5+z1*z10*z8-z1*z9*z6+z1*z8*z7+z1*z8*z6+z2*z3*
  z7-z2*z4*z6-z2*z10*z6-z2*z6S-z3*z10*z5-z3*z7*z5-z3*z6*
  z5+z4*z6*z5);
  fm[2][5] = 16*mQ4*z10+8*mQ2*(2*z1S+2*z1*z3+2*z1*z4+2*z1
  *z10+2*z1*z7+2*z1*z6+z3*z7+z4*z6)+8*(-2*pow3(z1)-2*z1S*z3-
  2*z1S*z4-2*z1S*z10-2*z1S*z7-2*z1S*z6-2*z1*z3*z4-
  z1*z3*z10-2*z1*z3*z6-z1*z4*z10-2*z1*z4*z7-z1*z10S-z1*
  z10*z7-z1*z10*z6-2*z1*z7*z6+z3S*z7-z3*z4*z7-z3*z4*z6+z3
  *z10*z7+z3*z7S-z3*z7*z6+z4S*z6+z4*z10*z6-z4*z7*z6+z4*
  z6S);
  fm[2][6] = 8*mQ4*(-2*z1+z10-z9-2*z7)+4*mQ2*(4*z1S+2*z1*
  z2+4*z1*z3+4*z1*z4+6*z1*z10-2*z1*z9+4*z1*z8+8*z1*z6-2*z1*
  z5+4*z2*z4+3*z2*z10+2*z2*z7-3*z3*z9-2*z3*z7-4*z4S-4*z4*
  z10+3*z4*z8+2*z4*z6+z10*z5-z9*z6+3*z8*z7+4*z7*z6)+8*(z1S
  *z9-z1S*z8-z1*z2*z7+z1*z2*z6+z1*z3*z9+z1*z3*z5+z1*z4*
  z9-z1*z4*z8-z1*z4*z5+z1*z10*z9+z1*z9*z6-z1*z8*z7-z2*z3*z7
  -z2*z4*z7+z2*z4*z6-z2*z10*z7+z3*z7*z5-z4S*z5-z4*z10*z5-
  z4*z6*z5);
  fm[2][7] = 8*mQ4*(z3+2*z4+3*z10+z7+2*z6)+4*mQ2*(-4*z1*z3-
  2*z1*z4-2*z1*z10+z1*z9-z1*z8-4*z1*z7-2*z1*z6+z2*z3+2*z2*
  z4+3*z2*z10+z2*z7+2*z2*z6-6*z3*z4-6*z3*z10-2*z3*z9-2*z3*
  z7-4*z3*z6-z3*z5-6*z4S-6*z4*z10-3*z4*z9-z4*z8-4*z4*z7-2
  *z4*z6-2*z4*z5-3*z10*z9-3*z10*z8-6*z10*z7-6*z10*z6+z10*z5
  +z9*z7-2*z8*z7-2*z8*z6-6*z7*z6+z7*z5-6*z6S+2*z6*z5)+4*(
  -z1S*z9+z1S*z8-2*z1*z2*z10-3*z1*z2*z7-3*z1*z2*z6+z1*
  z3*z9-z1*z3*z5+z1*z4*z9+z1*z4*z8+z1*z4*z5+z1*z10*z9+z1*
  z10*z8-z1*z9*z6+z1*z8*z6+z2*z3*z7-3*z2*z4*z7-z2*z4*z6-3*
  z2*z10*z7-3*z2*z10*z6-3*z2*z7*z6-3*z2*z6S-2*z3*z4*z5-z3
  *z10*z5-z3*z6*z5-z4S*z5-z4*z10*z5+z4*z6*z5);
  fm[2][8] = 8*mQ4*(z3+2*z4+3*z10+z7+2*z6)+4*mQ2*(-4*z1*z3-
  2*z1*z4-2*z1*z10-z1*z9+z1*z8-4*z1*z7-2*z1*z6+z2*z3+2*z2*
  z4+z2*z10-z2*z7-2*z2*z6-6*z3*z4-6*z3*z10-2*z3*z9+z3*z8-2*
  z3*z7-4*z3*z6+z3*z5-6*z4S-6*z4*z10-2*z4*z9-4*z4*z7-2*z4
  *z6+2*z4*z5-3*z10*z9-3*z10*z8-6*z10*z7-6*z10*z6+3*z10*z5-
  z9*z6-2*z8*z7-3*z8*z6-6*z7*z6+z7*z5-6*z6S+2*z6*z5)+4*(
  z1S*z9-z1S*z8-z1*z2*z7+z1*z2*z6-3*z1*z3*z5+z1*z4*z9-
  z1*z4*z8-3*z1*z4*z5+z1*z10*z9+z1*z10*z8-2*z1*z10*z5+z1*z9
  *z6+z1*z8*z7+z1*z8*z6-z2*z4*z7+z2*z4*z6-z2*z10*z7-z2*z10*
  z6-2*z2*z7*z6-z2*z6S-3*z3*z4*z5-3*z3*z10*z5+z3*z7*z5-3*
  z3*z6*z5-3*z4S*z5-3*z4*z10*z5-z4*z6*z5);
  fm[3][3] = 64*mQ6+16*mQ4*s3+32*mQ4*(z1+z2+2*z3+z8+z6
  +2*z5)+8*mQ2*s3*(-z1+2*z3-z6)+16*mQ2*(z2*z5-2*z3*
  z8-2*z3*z6+4*z3*z5+z8*z5)+8*s3*z3*z6-16*z3*z8*z5;
  fm[3][4] = 16*mQ4*(-4*z1-z3+z8-z6)+4*mQ2*s3*(-2*z1-z3-
  z6)+16*mQ2*(-2*z1S-3*z1*z2-2*z1*z3-3*z1*z8-2*z1*z6-3*
  z1*z5-2*z2*z3-2*z6*z5)-8*s3*z3*z6+8*(-z1*z2*z8-2*z1*z2
  *z5-z1*z8S-z1*z8*z5+z2S*z6-z2*z3*z5+z2*z8*z6-z2*z6*z5
  +z3*z8*z5+z3*z5S);
  fm[3][5] = 8*mQ4*(-2*z1+z10-z8-2*z6)+4*mQ2*(4*z1S+2*z1*
  z2+4*z1*z3+4*z1*z4+6*z1*z10+4*z1*z9-2*z1*z8+8*z1*z7-2*z1*
  z5+4*z2*z3+3*z2*z10+2*z2*z6-4*z3S-4*z3*z10+3*z3*z9+2*z3
  *z7-3*z4*z8-2*z4*z6+z10*z5+3*z9*z6-z8*z7+4*z7*z6)+8*(-z1S
  *z9+z1S*z8+z1*z2*z7-z1*z2*z6-z1*z3*z9+z1*z3*z8-z1*z3
  *z5+z1*z4*z8+z1*z4*z5+z1*z10*z8-z1*z9*z6+z1*z8*z7+z2*z3*
  z7-z2*z3*z6-z2*z4*z6-z2*z10*z6-z3S*z5-z3*z10*z5-z3*z7*
  z5+z4*z6*z5);
  fm[3][6] = 16*mQ6+4*mQ4*s3+16*mQ4*(-z1-z2+2*z3+2*z4+
  z10-z9-z8-z7-z6+z5)+4*mQ2*s3*(z1+2*z3+2*z4+z10+z7+z6
  )+8*mQ2*(4*z1*z3+4*z1*z4+4*z1*z10+4*z2*z3+4*z2*z4+4*z2*
  z10-z2*z5+4*z3*z5+4*z4*z5+2*z10*z5-z9*z5-z8*z5)-(8*s3*
  z1)*(z3+z4+z10)+16*z2*z5*(z3+z4+z10);
  fm[3][7] = 8*mQ4*(3*z3+6*z4+3*z10+z9+2*z8+2*z7+z6)+2*mQ2*
  s3*(z3+2*z4+2*z10-2*z7-z6)+4*mQ2*(4*z1*z3+8*z1*z4+4*
  z1*z10+2*z1*z9-2*z1*z8+2*z2*z3+10*z2*z4+5*z2*z10+2*z2*z9+
  z2*z8+2*z2*z7+4*z2*z6-7*z3*z9+2*z3*z8-8*z3*z7+4*z3*z6+4*
  z3*z5+5*z4*z8+4*z4*z6+8*z4*z5+5*z10*z5-z9*z8-z9*z6+z9*z5+
  z8S-z8*z7+2*z8*z6+2*z8*z5)+2*s3*(-z1*z10+z3*z7-2*z3*
  z6+z4*z6)+4*(-z1*z2*z9-2*z1*z2*z8+z1*z9*z8-z1*z8S+z2S
  *z7+2*z2S*z6+3*z2*z4*z5+2*z2*z10*z5-2*z2*z9*z6+z2*z8*z7
  +z2*z8*z6-2*z3*z9*z5+z3*z8*z5+z4*z8*z5);
  fm[3][8] = 8*mQ4*(3*z3+6*z4+3*z10+2*z9+z8+2*z7+z6)+2*mQ2*
  s3*(3*z3+6*z4+3*z10-2*z7-z6)+4*mQ2*(4*z1*z3+8*z1*z4+
  4*z1*z10+4*z2*z3+8*z2*z4+4*z2*z10-8*z3*z9+4*z3*z8-8*z3*z7
  +4*z3*z6+6*z3*z5+4*z4*z8+4*z4*z6+12*z4*z5+6*z10*z5+2*z9*
  z5+z8*z5)+4*s3*(-z1*z3-2*z1*z4-z1*z10+2*z3*z7-z3*z6-z4
  *z6)+8*z5*(z2*z3+2*z2*z4+z2*z10-2*z3*z9+z3*z8+z4*z8);
  fm[4][4] = 64*mQ6+16*mQ4*s3+32*mQ4*(z1+2*z2+z3+z8+2*
  z6+z5)+8*mQ2*s3*(-z1-z3+2*z6)+16*mQ2*(z2*z8+4*z2*
  z6+z2*z5-2*z3*z6-2*z8*z6)+8*s3*z3*z6-16*z2*z8*z6;
  fm[4][5] = 16*mQ6+8*mQ4*(-2*z1+z2-2*z3-2*z4-4*z10-z9+z8-4
  *z7+2*z6+z5)+8*mQ2*(-2*z1*z2-2*z2*z3-2*z2*z10-2*z2*z7+
  z2*z6+2*z3*z6-2*z4*z6+4*z10*z6-z9*z6-z8*z6)+16*z2*z6*(z3+
  z10);
  fm[4][6] = 16*mQ6-4*mQ4*s3+8*mQ4*(-2*z1+2*z2-4*z3-2*
  z4-8*z10+z9+z8-4*z7-2*z6+2*z5)-(4*mQ2*s3)*(z1+z3+z10
  +z7)+8*mQ2*(-2*z1*z2-2*z1*z10+z1*z9+z1*z8-2*z1*z5+z2S
  -4*z2*z3-5*z2*z10+z2*z9-3*z2*z7-z2*z6+z2*z5+z3*z9+2*z3*z7
  -3*z3*z5+z4*z8+2*z4*z6-z4*z5-5*z10*z5+z9*z8+z9*z6+z8*z7+
  z8*z5-4*z7*z5+z5S)-(16*z2*z5)*(z1+z3+z10+z7);
  fm[4][7] = 8*mQ4*(-z3-2*z4-3*z10-2*z9-z8-6*z7-3*z6)+2*mQ2
  *s3*(z3+2*z4-3*z10-6*z7-3*z6)+4*mQ2*(-4*z1*z10-8*z1*
  z7-4*z1*z6-6*z2*z10-2*z2*z9-z2*z8-12*z2*z7-6*z2*z6-4*z3*
  z7-4*z3*z6+8*z4*z6-4*z10*z5+8*z9*z6-4*z8*z7-4*z8*z6-8*z7*
  z5-4*z6*z5)+4*s3*(z1*z10+2*z1*z7+z1*z6+z3*z7+z3*z6-2*
  z4*z6)+8*z2*(-z10*z5+2*z9*z6-z8*z7-z8*z6-2*z7*z5-z6*z5);
  fm[4][8] = 8*mQ4*(-z3-2*z4-3*z10-z9-2*z8-6*z7-3*z6)+2*mQ2
  *s3*(z3+2*z4-2*z10-2*z7-z6)+4*mQ2*(-4*z1*z10-2*z1*z9
  +2*z1*z8-8*z1*z7-4*z1*z6-5*z2*z10-z2*z9-2*z2*z8-8*z2*z7-4
  *z2*z6+z3*z9-2*z3*z8-4*z3*z7-4*z3*z6-4*z3*z5+z4*z8+8*z4*
  z6-2*z4*z5-5*z10*z5+z9*z8+7*z9*z6-2*z9*z5-z8S-5*z8*z7-2
  *z8*z6-z8*z5-10*z7*z5-2*z6*z5)+2*s3*(z1*z10-z3*z7+2*z3
  *z6-z4*z6)+4*(-z1*z9*z8+z1*z9*z5+z1*z8S+2*z1*z8*z5-2*z2
  *z10*z5+2*z2*z9*z6-z2*z8*z7-z2*z8*z6-3*z2*z7*z5+2*z3*z9*
  z5-z3*z8*z5-2*z3*z5S-z4*z8*z5-z4*z5S);
  fm[5][5] = 16*mQ6+16*mQ4*(-z1-z3+z4-z10-z7+z6)+16*mQ2*(
  z3*z6+z4*z10+z4*z7+z4*z6+z10*z6)-16*z4*z10*z6;
  fm[5][6] = 16*mQ6+8*mQ4*(-2*z1+z2-4*z3+2*z4-4*z10+z9-z8-2
  *z7-2*z6+z5)+8*mQ2*(-2*z1*z5-2*z3*z5+4*z4*z10-z4*z9-z4*
  z8+2*z4*z7-2*z4*z6+z4*z5-2*z10*z5-2*z7*z5)+16*z4*z5*(z10+
  z7);
  fm[5][7] = 8*mQ4*(-2*z3-z4-3*z10-2*z7-z6)+4*mQ2*(2*z1*z3+
  4*z1*z4+2*z1*z10+z1*z9-z1*z8+2*z1*z7+4*z1*z6-2*z2*z3-z2*
  z4-3*z2*z10-2*z2*z7-z2*z6+6*z3S+6*z3*z4+6*z3*z10+z3*z9+
  3*z3*z8+2*z3*z7+4*z3*z6+2*z3*z5+6*z4*z10+2*z4*z8+4*z4*z7+
  2*z4*z6+z4*z5+3*z10*z9+3*z10*z8+6*z10*z7+6*z10*z6-z10*z5+
  2*z9*z7+2*z9*z6-z8*z6+6*z7S+6*z7*z6-2*z7*z5-z6*z5)+4*(-
  z1S*z9+z1S*z8+2*z1*z2*z10+3*z1*z2*z7+3*z1*z2*z6-z1*z3
  *z9-z1*z3*z8-z1*z3*z5-z1*z4*z8+z1*z4*z5-z1*z10*z9-z1*z10*
  z8-z1*z9*z7+z1*z8*z7+z2*z3*z7+3*z2*z3*z6-z2*z4*z6+3*z2*
  z10*z7+3*z2*z10*z6+3*z2*z7S+3*z2*z7*z6+z3S*z5+2*z3*z4
  *z5+z3*z10*z5-z3*z7*z5+z4*z10*z5+z4*z7*z5);
  fm[5][8] = 8*mQ4*(-2*z3-z4-3*z10-2*z7-z6)+4*mQ2*(2*z1*z3+
  4*z1*z4+2*z1*z10-z1*z9+z1*z8+2*z1*z7+4*z1*z6-2*z2*z3-z2*
  z4-z2*z10+2*z2*z7+z2*z6+6*z3S+6*z3*z4+6*z3*z10+2*z3*z8+
  2*z3*z7+4*z3*z6-2*z3*z5+6*z4*z10-z4*z9+2*z4*z8+4*z4*z7+2*
  z4*z6-z4*z5+3*z10*z9+3*z10*z8+6*z10*z7+6*z10*z6-3*z10*z5+
  3*z9*z7+2*z9*z6+z8*z7+6*z7S+6*z7*z6-2*z7*z5-z6*z5)+4*(
  z1S*z9-z1S*z8-z1*z2*z7+z1*z2*z6+z1*z3*z9-z1*z3*z8+3*
  z1*z3*z5+3*z1*z4*z5-z1*z10*z9-z1*z10*z8+2*z1*z10*z5-z1*z9
  *z7-z1*z9*z6-z1*z8*z7-z2*z3*z7+z2*z3*z6+z2*z10*z7+z2*z10*
  z6+z2*z7S+2*z2*z7*z6+3*z3S*z5+3*z3*z4*z5+3*z3*z10*z5+
  z3*z7*z5+3*z4*z10*z5+3*z4*z7*z5-z4*z6*z5);
  fm[6][6] = 64*mQ6+16*mQ4*s3+32*mQ4*(z1+z2+2*z4+z9+z7
  +2*z5)+8*mQ2*s3*(-z1+2*z4-z7)+16*mQ2*(z2*z5-2*z4*
  z9-2*z4*z7+4*z4*z5+z9*z5)+8*s3*z4*z7-16*z4*z9*z5;
  fm[6][7] = 8*mQ4*(-6*z3-3*z4-3*z10-2*z9-z8-z7-2*z6)+2*mQ2
  *s3*(-2*z3-z4-2*z10+z7+2*z6)+4*mQ2*(-8*z1*z3-4*z1*z4
  -4*z1*z10+2*z1*z9-2*z1*z8-10*z2*z3-2*z2*z4-5*z2*z10-z2*z9
  -2*z2*z8-4*z2*z7-2*z2*z6-5*z3*z9-4*z3*z7-8*z3*z5-2*z4*z9+
  7*z4*z8-4*z4*z7+8*z4*z6-4*z4*z5-5*z10*z5-z9S+z9*z8-2*z9
  *z7+z9*z6-2*z9*z5+z8*z7-z8*z5)+2*s3*(z1*z10-z3*z7+2*z4
  *z7-z4*z6)+4*(2*z1*z2*z9+z1*z2*z8+z1*z9S-z1*z9*z8-2*
  z2S*z7-z2S*z6-3*z2*z3*z5-2*z2*z10*z5-z2*z9*z7-z2*z9*z6+
  2*z2*z8*z7-z3*z9*z5-z4*z9*z5+2*z4*z8*z5);
  fm[6][8] = 8*mQ4*(-6*z3-3*z4-3*z10-z9-2*z8-z7-2*z6)+2*mQ2
  *s3*(-6*z3-3*z4-3*z10+z7+2*z6)+4*mQ2*(-8*z1*z3-4*z1*
  z4-4*z1*z10-8*z2*z3-4*z2*z4-4*z2*z10-4*z3*z9-4*z3*z7-12*
  z3*z5-4*z4*z9+8*z4*z8-4*z4*z7+8*z4*z6-6*z4*z5-6*z10*z5-z9
  *z5-2*z8*z5)+4*s3*(2*z1*z3+z1*z4+z1*z10+z3*z7+z4*z7-2*
  z4*z6)+8*z5*(-2*z2*z3-z2*z4-z2*z10-z3*z9-z4*z9+2*z4*z8);
  fm[7][7] = 72*mQ4*z10+18*mQ2*s3*z10+8*mQ2*(z1*z10+9*
  z2*z10+7*z3*z7+2*z3*z6+2*z4*z7+7*z4*z6+z10*z5+2*z9*z7+7*
  z9*z6+7*z8*z7+2*z8*z6)+2*s3*(-z1*z10-7*z3*z7-2*z3*z6-2
  *z4*z7-7*z4*z6)+4*z2*(z10*z5+2*z9*z7+7*z9*z6+7*z8*z7+2*z8
  *z6);
  fm[7][8] = 72*mQ4*z10+2*mQ2*s3*z10+4*mQ2*(2*z1*z10+
  10*z2*z10+7*z3*z9+2*z3*z8+14*z3*z7+4*z3*z6+2*z4*z9+7*z4*
  z8+4*z4*z7+14*z4*z6+10*z10*z5+z9S+7*z9*z8+2*z9*z7+7*z9*
  z6+z8S+7*z8*z7+2*z8*z6)+2*s3*(7*z1*z10-7*z3*z7-2*z3*
  z6-2*z4*z7-7*z4*z6)+2*(-2*z1*z9S-14*z1*z9*z8-2*z1*z8S
  +2*z2*z10*z5+2*z2*z9*z7+7*z2*z9*z6+7*z2*z8*z7+2*z2*z8*z6+
  7*z3*z9*z5+2*z3*z8*z5+2*z4*z9*z5+7*z4*z8*z5);
  fm[8][8] = 72*mQ4*z10+18*mQ2*s3*z10+8*mQ2*(z1*z10+z2
  *z10+7*z3*z9+2*z3*z8+7*z3*z7+2*z3*z6+2*z4*z9+7*z4*z8+2*z4
  *z7+7*z4*z6+9*z10*z5)+2*s3*(-z1*z10-7*z3*z7-2*z3*z6-2*
  z4*z7-7*z4*z6)+4*z5*(z2*z10+7*z3*z9+2*z3*z8+2*z4*z9+7*z4*
  z8);
  double fm99 = -4*mQ4*z10-mQ2*s3*z10+4*mQ2*(-z1*z10-z2*z10+
  z3*z7+z4*z6-z10*z5+z9*z6+z8*z7)+s3*(z1*z10-z3*z7-z4*z6
  )+2*z2*(-z10*z5+z9*z6+z8*z7);
  double fm910 = -4*mQ4*z10-mQ2*s3*z10+2*mQ2*(-2*z1*z10-2*z2*
  z10+2*z3*z9+2*z3*z7+2*z4*z6-2*z10*z5+z9*z8+2*z8*z7)+s3
  *(z1*z10-z3*z7-z4*z6)+2*(-z1*z9*z8-z2*z10*z5+z2*z8*z7+z3*
  z9*z5);
  double fmxx = -4*mQ4*z10-mQ2*s3*z10+2*mQ2*(-2*z1*z10-2*z2*
  z10+2*z4*z8+2*z4*z6+2*z3*z7-2*z10*z5+z9*z8+2*z9*z6)+s3
  *(z1*z10-z3*z7-z4*z6)+2*(-z1*z9*z8-z2*z10*z5+z2*z9*z6+z4*
  z8*z5);
  fm910 = 0.5*(fmxx+fm910);
  double fm1010 = -4*mQ4*z10-mQ2*s3*z10+4*mQ2*(-z1*z10-z2*z10+
  z3*z7+z4*z6-z10*z5+z9*z3+z8*z4)+s3*(z1*z10-z3*z7-z4*z6
  )+2*z5*(-z10*z2+z9*z3+z8*z4);
  fm[7][7] -= 2. * fm99;
  fm[7][8] -= 2. * fm910;
  fm[8][8] -= 2. * fm1010;

  // Propagators.
  double ss1 = (pTemp[1] + pTemp[3]).m2Calc() - mQ2;
  double ss2 = (pTemp[1] + pTemp[4]).m2Calc() - mQ2;
  double ss3 = (pTemp[1] + pTemp[5]).m2Calc() - mQ2;
  double ss4 = (pTemp[2] + pTemp[3]).m2Calc() - mQ2;
  double ss5 = (pTemp[2] + pTemp[4]).m2Calc() - mQ2;
  double ss6 = (pTemp[2] + pTemp[5]).m2Calc() - mQ2;
  double ss7 = sH;

  // Propagator combinations.
  double dz[9];
  dz[1]      = ss1 * ss6;
  dz[2]      = ss2 * ss6;
  dz[3]      = ss2 * ss4;
  dz[4]      = ss1 * ss5;
  dz[5]      = ss3 * ss5;
  dz[6]      = ss3 * ss4;
  dz[7]      = ss7 * ss1;
  dz[8]      = ss7 * ss4;

  // Colour factors.
  double clr[9][9];
  for (int i = 1; i < 4; ++i)
  for (int j = 1; j < 4; ++j) {
    clr[i][j]     = 16. / 3.;
    clr[i][j+3]   = -2. / 3.;
    clr[i+3][j]   = -2. / 3.;
    clr[i+3][j+3] = 16. / 3.;
  }
  for (int i = 1; i < 4; ++i)
  for (int j = 1; j < 3; ++j) {
    clr[i][j+6]   = -6.;
    clr[i+3][j+6] = 6.;
    clr[j+6][i]   = -6.;
    clr[j+6][i+3] = 6.;
  }
  for (int i = 1; i < 3; ++i)
  for (int j = 1; j < 3; ++j)
    clr[i+6][j+6] = 12.;

  // Produce final result: matrix elements * colours * propagators.
  double wtSum = 0.;
  for (int i = 1; i < 9; ++i)
  for (int j = i; j < 9; ++j) {
    double fac = (j == i) ? 4. : 8.;
    wtSum += fm[i][j] * fac * clr[i][j] / (dz[i] * dz[j]);
  }
  wtSum *= -1./256.;

  // Combine factors.
  sigma  = prefac * alpEM * pow2(alpS) * mQ2run * wtSum *pow2(coup2Q);

  // Secondary width for H, Q and Qbar (latter for top only).
  // (H can be H0 SM or H1, H2, A3 from BSM).
  sigma *= openFracTriplet;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3gg2HQQbar::setIdColAcol() {

  // Pick out-flavours by relative CKM weights.
  setId( id1, id2, idRes, idNew, -idNew);

  // Colour flow topologies.
  if (rndmPtr->flat() < 0.5) setColAcol( 1, 2, 2, 3, 0, 0, 1, 0, 0, 3);
  else                    setColAcol( 1, 2, 3, 1, 0, 0, 3, 0, 0, 2);

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma3gg2HQQbar::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma3qqbar2HQQbar class.
// Cross section for q qbar -> H0 Q Qbar (Q Qbar fusion of SM Higgs).
// REDUCE output and part of the rest courtesy Z. Kunszt,
// see Z. Kunszt, Nucl. Phys. B247 (1984) 339.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma3qqbar2HQQbar::initProc() {

  // Properties specific to Higgs state for the "q qbar -> H ttbar" process.
  // (H can be H0 SM or H1, H2, A3 from BSM).

  if (higgsType == 0 && idNew == 6) {
    nameSave = "q qbar -> H t tbar (SM)";
    codeSave = 909;
    idRes    = 25;
    coup2Q   = 1.;
  }
  else if (higgsType == 1 && idNew == 6) {
    nameSave = "q qbar -> h0(H1) t tbar";
    codeSave = 1009;
    idRes    = 25;
    coup2Q   = settingsPtr->parm("HiggsH1:coup2u");
  }
  else if (higgsType == 2 && idNew == 6) {
    nameSave = "q qbar -> H0(H2) t tbar";
    codeSave = 1029;
    idRes    = 35;
    coup2Q   = settingsPtr->parm("HiggsH2:coup2u");
  }
  else if (higgsType == 3 && idNew == 6) {
    nameSave = "q qbar -> A0(A3) t tbar";
    codeSave = 1049;
    idRes    = 36;
    coup2Q   = settingsPtr->parm("HiggsA3:coup2u");
  }

 // Properties specific to Higgs state for the "q qbar -> H b bbar" process.
 // (H can be H0 SM or H1, H2, A3 from BSM).
  if (higgsType == 0 && idNew == 5) {
    nameSave = "q qbar -> H b bbar (SM)";
    codeSave = 913;
    idRes    = 25;
    coup2Q   = 1.;
  }
  else if (higgsType == 1 && idNew == 5) {
    nameSave = "q qbar -> h0(H1) b bbar";
    codeSave = 1013;
    idRes    = 25;
    coup2Q   = settingsPtr->parm("HiggsH1:coup2d");
  }
  else if (higgsType == 2 && idNew == 5) {
    nameSave = "q qbar -> H0(H2) b bbar";
    codeSave = 1033;
    idRes    = 35;
    coup2Q   = settingsPtr->parm("HiggsH2:coup2d");
  }
  else if (higgsType == 3 && idNew == 5) {
    nameSave = "q qbar -> A0(A3) b bbar";
    codeSave = 1053;
    idRes    = 36;
    coup2Q   = settingsPtr->parm("HiggsA3:coup2d");
  }

  // Common mass and coupling factors.
  double mWS      = pow2(particleDataPtr->m0(24));
  prefac          = (4. * M_PI / couplingsPtr->sin2thetaW()) * pow2(4. * M_PI)
                  * 0.25 / mWS;

  // Secondary open width fraction.
  openFracTriplet = particleDataPtr->resOpenFrac(idRes, idNew, -idNew);

}

//--------------------------------------------------------------------------

// Evaluate sigma(sHat), part independent of incoming flavour.

void Sigma3qqbar2HQQbar::sigmaKin() {

  // Running mass of heavy quark.
  double mQ2run = pow2( particleDataPtr->mRun(idNew, mH) );

  // Linear combination of p_Q and p_Qbar to ensure common mass.
  double mQ2  = m4 * m5;
  double epsi = 0.;
  if (m4 != m5) {
    double s45  = (p4cm + p5cm).m2Calc();
    mQ2       = 0.5 * (s4 + s5) - 0.25 * pow2(s4 - s5) / s45;
    epsi      = 0.5 * (s5 - s4) / s45;
  }

  // Set up kinematics: q(4) qbar(5) -> H(3) Q(1) Qbar(2) in outgoing sense.
  Vec4 pTemp[6];
  pTemp[4]   = Vec4( 0., 0., -0.5* mH, -0.5* mH);
  pTemp[5]   = Vec4( 0., 0.,  0.5* mH, -0.5* mH);
  pTemp[1]   = p4cm + epsi * (p4cm + p5cm);
  pTemp[2]   = p5cm - epsi * (p4cm + p5cm);
  pTemp[3]   = p3cm;

  // Four-product combinations.
  double z1  = pTemp[1] * pTemp[2];
  double z2  = pTemp[1] * pTemp[3];
  double z3  = pTemp[1] * pTemp[4];
  double z4  = pTemp[1] * pTemp[5];
  double z5  = pTemp[2] * pTemp[3];
  double z6  = pTemp[2] * pTemp[4];
  double z7  = pTemp[2] * pTemp[5];
  double z8  = pTemp[3] * pTemp[4];
  double z9  = pTemp[3] * pTemp[5];
  double z10 = pTemp[4] * pTemp[5];

  // Powers required as shorthand in matriz elements.
  double mQ4  = mQ2 * mQ2;

  // Evaluate matrix elements for q + qbar -> Q + Qbar + H.
  // (H can be H0 SM or H1, H2, A3 from BSM).
  double a11 = -8.*mQ4*z10-2.*mQ2*s3*z10-(8.*mQ2)*(z2*z10+z3
  *z7+z4*z6+z9*z6+z8*z7)+2.*s3*(z3*z7+z4*z6)-(4.*z2)*(z9
  *z6+z8*z7);
  double a12 = -8.*mQ4*z10+4.*mQ2*(-z2*z10-z3*z9-2.*z3*z7-z4*z8-
  2.*z4*z6-z10*z5-z9*z8-z9*z6-z8*z7)+2.*s3*(-z1*z10+z3*z7
  +z4*z6)+2.*(2.*z1*z9*z8-z2*z9*z6-z2*z8*z7-z3*z9*z5-z4*z8*
  z5);
  double a22 = -8.*mQ4*z10-2.*mQ2*s3*z10-(8.*mQ2)*(z3*z9+z3*
  z7+z4*z8+z4*z6+z10*z5)+2.*s3*(z3*z7+z4*z6)-(4.*z5)*(z3
  *z9+z4*z8);

  // Propagators and propagator combinations.
  double ss1 = (pTemp[1] + pTemp[3]).m2Calc() - mQ2;
  double ss4 = (pTemp[2] + pTemp[3]).m2Calc() - mQ2;
  double ss7 = sH;
  double dz7 = ss7 * ss1;
  double dz8 = ss7 * ss4;

  // Produce final result: matrix elements * propagators.
  a11 /=  (dz7 * dz7);
  a12 /=  (dz7 * dz8);
  a22 /=  (dz8 * dz8);
  double wtSum = -(a11 + a22 + 2.*a12) * (8./9.);

  // Combine factors.
  sigma = prefac * alpEM * pow2(alpS) * mQ2run * wtSum * pow2(coup2Q);

  // Secondary width for H, Q and Qbar (latter for top only).
  // (H can be H0 SM or H1, H2, A3 from BSM).
  sigma *= openFracTriplet;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3qqbar2HQQbar::setIdColAcol() {

  // Pick out-flavours by relative CKM weights.
  setId( id1, id2, idRes, idNew, -idNew);

  // Colour flow topologies.
  if (id1 > 0) setColAcol( 1, 0, 0, 2, 0, 0, 1, 0, 0, 2);
  else         setColAcol( 0, 1, 2, 0, 0, 0, 2, 0, 0, 1);

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma3qqbar2HQQbar::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma2qg2Hq class.
// Cross section for q g -> H q.
// (H can be H0 SM or H1, H2, A3 from BSM).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qg2Hq::initProc() {

  // Properties specific to Higgs state for the "c g -> H c" process.
  // (H can be H0 SM or H1, H2, A3 from BSM).
  if (higgsType == 0 && idNew == 4) {
    nameSave = "c g -> H c (SM)";
    codeSave = 911;
    idRes    = 25;
  }
  else if (higgsType == 1 && idNew == 4) {
    nameSave = "c g -> h0(H1) c";
    codeSave = 1011;
    idRes    = 25;
  }
  else if (higgsType == 2 && idNew == 4) {
    nameSave = "c g -> H0(H2) c";
    codeSave = 1031;
    idRes    = 35;
  }
  else if (higgsType == 3 && idNew == 4) {
    nameSave = "c g -> A0(A3) c";
    codeSave = 1051;
    idRes    = 36;
  }

 // Properties specific to Higgs state for the "b g -> H b" process.
 // (H can be H0 SM or H1, H2, A3 from BSM).
 if (higgsType == 0 && idNew == 5) {
    nameSave = "b g -> H b (SM)";
    codeSave = 911;
    idRes    = 25;
  }
  else if (higgsType == 1 && idNew == 5) {
    nameSave = "b g -> h0(H1) b";
    codeSave = 1011;
    idRes    = 25;
  }
  else if (higgsType == 2 && idNew == 5) {
    nameSave = "b g -> H0(H2) b";
    codeSave = 1031;
    idRes    = 35;
  }
  else if (higgsType == 3 && idNew == 5) {
    nameSave = "b g -> A0(A3) b";
    codeSave = 1051;
    idRes    = 36;
  }

  // Standard parameters.
  m2W       = pow2( particleDataPtr->m0(24) );
  thetaWRat = 1. / (24. * couplingsPtr->sin2thetaW());

  // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(idRes);


}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2qg2Hq::sigmaKin() {

  // Running mass provides coupling.
  double m2Run = pow2( particleDataPtr->mRun(idNew, mH) );

  // Cross section, including couplings and kinematics.
  sigma = (M_PI / sH2) * alpS * alpEM * thetaWRat * (m2Run/m2W)
    * (sH / (s4 - uH) + 2. * s4 * (s3 - uH) / pow2(s4 - uH)
      + (s4 - uH) / sH - 2. * s4 / (s4 - uH)
      + 2. * (s3 - uH)  * (s3 - s4 - sH) / ((s4 - uH) * sH) );

  // Include secondary width for H0, H1, H2 or A3. Done.
  sigma *= openFrac;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma2qg2Hq::sigmaHat() {

  // Check that specified flavour present.
  if (abs(id1) != idNew && abs(id2) != idNew) return 0.;

  // Answer.
  return sigma;

}


//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2Hq::setIdColAcol() {

  // Flavour set up for q g -> H0 q.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, idRes, idq);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21);

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2qg2Hq::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma2gg2Hglt class.
// Cross section for g g -> H g (H SM Higgs or BSM Higgs) via top loop.
// (H can be H0 SM or H1, H2, A3 from BSM).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2Hglt::initProc() {

  // Properties specific to Higgs state.
  if (higgsType == 0) {
    nameSave = "g g -> H g (SM; top loop)";
    codeSave = 914;
    idRes    = 25;
  }
  else if (higgsType == 1) {
    nameSave = "g g -> h0(H1) g (BSM; top loop)";
    codeSave = 1014;
    idRes    = 25;
  }
  else if (higgsType == 2) {
    nameSave = "g g -> H0(H2) g (BSM; top loop)";
    codeSave = 1034;
    idRes    = 35;
  }
  else if (higgsType == 3) {
    nameSave = "g g -> A0(A3) g (BSM; top loop)";
    codeSave = 1054;
    idRes    = 36;
  }

  // Normalization factor by g g -> H partial width.
  // (H can be H0 SM or H1, H2, A3 from BSM).
  double mHiggs = particleDataPtr->m0(idRes);
  widHgg = particleDataPtr->resWidthChan(idRes, mHiggs, 21, 21);

   // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(idRes);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2gg2Hglt::sigmaKin() {

  // Evaluate cross section. Secondary width for H0, H1, H2 or A3.
  sigma  = (M_PI / sH2) * (3. / 16.) * alpS * (widHgg / m3)
    * (sH2 * sH2 + tH2 * tH2 + uH2 * uH2 + pow4(s3))
    / (sH * tH * uH * s3);
  sigma *= openFrac;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2Hglt::setIdColAcol() {

  // Flavour set up for g g -> H g trivial.
  // (H can be H0 SM or H1, H2, A3 from BSM).
  setId( 21, 21, idRes, 21);

  // Colour flow topologies: random choice between two mirrors.
  if (rndmPtr->flat() < 0.5) setColAcol( 1, 2, 2, 3, 0, 0, 1, 3);
  else                    setColAcol( 1, 2, 3, 1, 0, 0, 3, 2);

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2gg2Hglt::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma2qg2Hqlt class.
// Cross section for q g -> H q (H SM or BSM Higgs) via top loop.
// (H can be H0 SM or H1, H2, A3 from BSM).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qg2Hqlt::initProc() {

  // Properties specific to Higgs state.
  if (higgsType == 0) {
    nameSave = "q g -> H q (SM; top loop)";
    codeSave = 915;
    idRes    = 25;
  }
  else if (higgsType == 1) {
    nameSave = "q g -> h0(H1) q (BSM; top loop)";
    codeSave = 1015;
    idRes    = 25;
  }
  else if (higgsType == 2) {
    nameSave = "q g -> H0(H2) q (BSM; top loop)";
    codeSave = 1035;
    idRes    = 35;
  }
  else if (higgsType == 3) {
    nameSave = "q g -> A0(A3) q (BSM; top loop)";
    codeSave = 1055;
    idRes    = 36;
  }

  // Normalization factor by g g -> H partial width.
  // (H can be H0 SM or H1, H2, A3 from BSM).
  double mHiggs = particleDataPtr->m0(idRes);
  widHgg = particleDataPtr->resWidthChan(idRes, mHiggs, 21, 21);

  // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(idRes);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat, part independent of incoming flavour).

void Sigma2qg2Hqlt::sigmaKin() {

  // Evaluate cross section. Secondary width for H0, H1, H2 or A3.
  sigma  = (M_PI / sH2) * (1. / 12.) * alpS * (widHgg / m3)
    * (sH2 + uH2) / (-tH * s3);
  sigma *= openFrac;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2Hqlt::setIdColAcol() {

  // Flavour set up for q g -> H q.
  // (H can be H0 SM or H1, H2, A3 from BSM).
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, idRes, idq);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21);

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2qg2Hqlt::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma2qqbar2Hglt class.
// Cross section for q qbar -> H g (H SM or BSM Higgs) via top loop.
// (H can be H0 SM or H1, H2, A3 from BSM).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qqbar2Hglt::initProc() {

  // Properties specific to Higgs state.
  if (higgsType == 0) {
    nameSave = "q qbar -> H g (SM; top loop)";
    codeSave = 916;
    idRes    = 25;
  }
  else if (higgsType == 1) {
    nameSave = "q qbar -> h0(H1) g (BSM; top loop)";
    codeSave = 1016;
    idRes    = 25;
  }
  else if (higgsType == 2) {
    nameSave = "q qbar -> H0(H2) g (BSM; top loop)";
    codeSave = 1036;
    idRes    = 35;
  }
  else if (higgsType == 3) {
    nameSave = "q qbar -> A0(A3) g (BSM; top loop)";
    codeSave = 1056;
    idRes    = 36;
  }

  // Normalization factor by g g -> H partial width.
  // (H can be H0 SM or H1, H2, A3 from BSM).
  double mHiggs = particleDataPtr->m0(idRes);
  widHgg = particleDataPtr->resWidthChan(idRes, mHiggs, 21, 21);

  // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(idRes);


}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2qqbar2Hglt::sigmaKin() {

  // Evaluate cross section. Secondary width for H0, H1, H2 or A3.
  sigma  = (M_PI / sH2) * (2. / 9.) * alpS * (widHgg / m3)
    * (tH2 + uH2) / (sH * s3);
  sigma *= openFrac;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2Hglt::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idRes, 21);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2qqbar2Hglt::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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


//==========================================================================

// Sigma1ffbar2Hchg class.
// Cross section for f fbar -> H+- (f is quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ffbar2Hchg::initProc() {

  // Find pointer to H+-.
  HResPtr = particleDataPtr->particleDataEntryPtr(37);

  // Store H+- mass and width for propagator.
  mRes      = HResPtr->m0();
  GammaRes  = HResPtr->mWidth();
  m2Res     = mRes*mRes;
  GamMRat   = GammaRes / mRes;

  // Couplings.
  m2W       = pow2(particleDataPtr->m0(24));
  thetaWRat = 1. / (8. * couplingsPtr->sin2thetaW());
  tan2Beta  = pow2(settingsPtr->parm("HiggsHchg:tanBeta"));

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma1ffbar2Hchg::sigmaKin() {

  // Set up Breit-Wigner. Width out only includes open channels.
  sigBW    = 4. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  widthOutPos = HResPtr->resWidthOpen( 37, mH);
  widthOutNeg = HResPtr->resWidthOpen(-37, mH);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma1ffbar2Hchg::sigmaHat() {

  // Only allow generation-diagonal states.
  int id1Abs     = abs(id1);
  int id2Abs     = abs(id2);
  int idUp       = max(id1Abs, id2Abs);
  int idDn       = min(id1Abs, id2Abs);
  if (idUp%2 != 0 || idUp - idDn != 1) return 0.;

  // Calculate mass-dependent incoming width. Total cross section.
  double m2RunUp = pow2(particleDataPtr->mRun(idUp, mH));
  double m2RunDn = pow2(particleDataPtr->mRun(idDn, mH));
  double widthIn = alpEM * thetaWRat * (mH/m2W)
     * (m2RunDn * tan2Beta + m2RunUp / tan2Beta);
  int idUpChg    = (id1Abs%2 == 0) ? id1 : id2;
  double sigma   = (idUpChg > 0) ? widthIn * sigBW * widthOutPos
                                 : widthIn * sigBW * widthOutNeg;

  // Colour factor. Answer.
  if (idUp < 9) sigma /= 3.;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2Hchg::setIdColAcol() {

  // Charge of Higgs. Fill flavours.
  int idUpChg = (abs(id1)%2 == 0) ? id1 : id2;
  int idHchg  = (idUpChg > 0) ? 37 : -37;
  setId( id1, id2, idHchg);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma1ffbar2Hchg::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma2qg2Hq class.
// Cross section for q g -> H+- q'.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qg2Hchgq::initProc() {

  // Standard parameters.
  m2W       = pow2( particleDataPtr->m0(24) );
  thetaWRat = 1. / (24. * couplingsPtr->sin2thetaW());
  tan2Beta  = pow2(settingsPtr->parm("HiggsHchg:tanBeta"));

  // Incoming flavour within same doublet. Uptype and downtype flavours.
  idOld     = (idNew%2 == 0) ? idNew - 1 : idNew + 1;
  idUp      = max(idOld, idNew);
  idDn      = min(idOld, idNew);

  // Secondary open width fraction.
  openFracPos = (idOld%2 == 0) ? particleDataPtr->resOpenFrac( 37,  idNew)
                               : particleDataPtr->resOpenFrac(-37,  idNew);
  openFracNeg = (idOld%2 == 0) ? particleDataPtr->resOpenFrac(-37, -idNew)
                               : particleDataPtr->resOpenFrac( 37, -idNew);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2qg2Hchgq::sigmaKin() {

  // Running masses provides coupling.
  double m2RunUp = pow2(particleDataPtr->mRun(idUp, mH));
  double m2RunDn = pow2(particleDataPtr->mRun(idDn, mH));

  // Cross section, including couplings and kinematics.
  sigma = (M_PI / sH2) * alpS * alpEM * thetaWRat
    * (m2RunDn * tan2Beta + m2RunUp / tan2Beta) / m2W
    * (sH / (s4 - uH) + 2. * s4 * (s3 - uH) / pow2(s4 - uH)
      + (s4 - uH) / sH - 2. * s4 / (s4 - uH)
      + 2. * (s3 - uH)  * (s3 - s4 - sH) / ((s4 - uH) * sH) );

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma2qg2Hchgq::sigmaHat() {

  // Check that specified flavour present.
  if (abs(id1) != idOld && abs(id2) != idOld) return 0.;

  // Answer.
  return (id1 == idOld || id2 == idOld) ? sigma * openFracPos
                                        : sigma * openFracNeg;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2Hchgq::setIdColAcol() {

  // Flavour set up for q g -> H+- q'.
  int idq = (id2 == 21) ? id1 : id2;
  id3     = ( (idq > 0 && idOld%2 == 0) || (idq < 0 && idOld%2 != 0) )
            ? 37 : -37;
  id4     = (idq > 0) ? idNew : -idNew;
  setId( id1, id2, id3, id4);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21);

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2qg2Hchgq::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma2ffbar2A3H12 class.
// Cross section for f fbar -> A0(H_3) h0(H_1) or A0(H_3) H0(H_2).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2A3H12::initProc() {

  // Set up whether h0(H_1) or H0(H_2).
  higgs12    = (higgsType == 1) ? 25 : 35;
  codeSave   = (higgsType == 1) ? 1081 : 1082;
  nameSave   = (higgsType == 1) ? "f fbar -> A0(H3) h0(H1)"
                                : "f fbar -> A0(H3) H0(H2)";
  coupZA3H12 = (higgsType == 1) ? settingsPtr->parm("HiggsA3:coup2H1Z")
                                : settingsPtr->parm("HiggsA3:coup2H2Z");

  // Standard parameters.
  double mZ  = particleDataPtr->m0(23);
  double GammaZ = particleDataPtr->mWidth(23);
  m2Z        = mZ * mZ;
  mGammaZ    = mZ * GammaZ;
  thetaWRat  = 1. / (4. * couplingsPtr->sin2thetaW()
             * couplingsPtr->cos2thetaW());

  // Secondary open width fraction.
  openFrac   = particleDataPtr->resOpenFrac(36, higgs12);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2ffbar2A3H12::sigmaKin() {

  // Common kinematics factora.
  sigma0 = (M_PI / sH2) * pow2(alpEM * thetaWRat * coupZA3H12)
    * (uH * tH - s3 * s4) / ( pow2(sH - m2Z) + pow2(mGammaZ) );

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma2ffbar2A3H12::sigmaHat() {

  // Couplings for incoming flavour.
  int idAbs    = abs(id1);
  double lIn   = couplingsPtr->lf(idAbs);
  double rIn   = couplingsPtr->rf(idAbs);

  // Combine to total cross section. Colour factor.
  double sigma = (pow2(lIn) + pow2(rIn)) * sigma0 * openFrac;
  if (idAbs < 9) sigma /= 3.;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2A3H12::setIdColAcol() {

  // Flavours trivial
  setId( id1, id2, 36, higgs12);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2ffbar2A3H12::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma2ffbar2HchgH12 class.
// Cross section for f fbar -> H+- h0(H_1) or H+- H0(H_2).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2HchgH12::initProc() {

  // Set up whether h0(H_1) or H0(H_2).
  higgs12    = (higgsType == 1) ? 25 : 35;
  codeSave   = (higgsType == 1) ? 1083 : 1084;
  nameSave   = (higgsType == 1) ? "f fbar' -> H+- h0(H1)"
                                : "f fbar' -> H+- H0(H2)";
  coupWHchgH12 = (higgsType == 1) ? settingsPtr->parm("HiggsHchg:coup2H1W")
                                  : settingsPtr->parm("HiggsHchg:coup2H2W");

  // Standard parameters.
  double mW  = particleDataPtr->m0(24);
  double GammaW = particleDataPtr->mWidth(24);
  m2W        = mW * mW;
  mGammaW    = mW * GammaW;
  thetaWRat  = 1. / (2. * couplingsPtr->sin2thetaW());

  // Secondary open width fraction.
  openFracPos   = particleDataPtr->resOpenFrac( 37, higgs12);
  openFracNeg   = particleDataPtr->resOpenFrac(-37, higgs12);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2ffbar2HchgH12::sigmaKin() {

  // Common kinematics factora.
  sigma0 = 0.5 * (M_PI / sH2) * pow2(alpEM * thetaWRat * coupWHchgH12)
    * (uH * tH - s3 * s4) / ( pow2(sH - m2W) + pow2(mGammaW) );

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma2ffbar2HchgH12::sigmaHat() {

  // Combine to total cross section. CKM and colour factor.
  int idUp = (abs(id1)%2 == 0) ? id1 : id2;
  double sigma = (idUp > 0) ? sigma0 * openFracPos : sigma0 * openFracNeg;
  if (abs(id1) < 9) sigma *= couplingsPtr->V2CKMid(abs(id1), abs(id2)) / 3.;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2HchgH12::setIdColAcol() {

  // Charge of Higgs. Fill flavours.
  int idUpChg = (abs(id1)%2 == 0) ? id1 : id2;
  int idHchg  = (idUpChg > 0) ? 37 : -37;
  setId( id1, id2, idHchg, higgs12);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2ffbar2HchgH12::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

// Sigma2ffbar2HposHneg class.
// Cross section for q g -> H+- q'.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2HposHneg::initProc() {

  // Standard parameters.
  double mZ = particleDataPtr->m0(23);
  double GammaZ = particleDataPtr->mWidth(23);
  m2Z       = mZ * mZ;
  mGammaZ   = mZ * GammaZ;
  thetaWRat = 1. / (4. * couplingsPtr->sin2thetaW()
            * couplingsPtr->cos2thetaW());

  // Charged Higgs coupling to gamma and Z0.
  eH        = -1.;
  lH        = -1. + 2. * couplingsPtr->sin2thetaW();

  // Secondary open width fraction.
  openFrac  = particleDataPtr->resOpenFrac(37, -37);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2ffbar2HposHneg::sigmaKin() {

  // Common kinematics factora.
  double preFac = M_PI * pow2(alpEM) * ((uH * tH - s3 * s4) / sH2);
  double propZ  = 1. / ( pow2(sH - m2Z) + pow2(mGammaZ) );

  // Separate parts for gamma*, interference and Z0.
  gamSig    = preFac * 2. * pow2(eH) / sH2;
  intSig    = preFac * 2. * eH * lH * thetaWRat * propZ * (sH - m2Z) / sH;
  resSig    = preFac * pow2(lH * thetaWRat) * propZ;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma2ffbar2HposHneg::sigmaHat() {

  // Couplings for incoming flavour.
  int idAbs    = abs(id1);
  double eIn   = couplingsPtr->ef(idAbs);
  double lIn   = couplingsPtr->lf(idAbs);
  double rIn   = couplingsPtr->rf(idAbs);

  // Combine to total cross section. Colour factor.
  double sigma = (pow2(eIn) * gamSig + eIn * (lIn + rIn) * intSig
    + (pow2(lIn) + pow2(rIn)) * resSig) * openFrac;
  if (idAbs < 9) sigma /= 3.;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2HposHneg::setIdColAcol() {

  // Flavours trivial
  setId( id1, id2, 37, -37);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma2ffbar2HposHneg::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

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

//==========================================================================

} // end namespace Pythia8
