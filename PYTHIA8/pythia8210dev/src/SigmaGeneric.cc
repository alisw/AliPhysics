// SigmaGeneric.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Johan Bijnens, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for various generic
// production processes, to be used as building blocks for some BSM processes.
// Currently represented by QCD pair production of colour triplet objects,
// with spin either 0, 1/2 or 1.

// Cross sections are only provided for fixed m3 = m4, so do some gymnastics:
// i) s34Avg picked so that beta34 same when s3, s4 -> s34Avg.
// ii) tHQ = tH - mQ^2 = -0.5 sH (1 - beta34 cos(thetaH)) for m3 = m4 = mQ,
//     but tH - uH = sH beta34 cos(thetaH) also for m3 != m4, so use
//     tH, uH selected for m3 != m4 to derive tHQ, uHQ valid for m3 = m4.

#include "Pythia8/SigmaGeneric.h"

namespace Pythia8 {

//==========================================================================

// Sigma2gg2qGqGbar class.
// Cross section for g g -> qG qGbar (generic quark of spin 0, 1/2 or 1).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2qGqGbar::initProc() {

  // Number of colours. Anomalous coupling kappa - 1 used for vector state.
  nCHV         = settingsPtr->mode("HiddenValley:Ngauge");
  kappam1      = settingsPtr->parm("HiddenValley:kappa") - 1.;
  hasKappa     = (abs(kappam1) > 1e-8);

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(idNew, -idNew);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2gg2qGqGbar::sigmaKin() {

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  double delta  = 0.25 * pow2(s3 - s4) / sH;
  double s34Avg = 0.5 * (s3 + s4) - delta;
  double tHavg  = tH - delta;
  double uHavg  = uH - delta;
  double tHQ    = -0.5 * (sH - tH + uH);
  double uHQ    = -0.5 * (sH + tH - uH);
  double tHQ2   = tHQ * tHQ;
  double uHQ2   = uHQ * uHQ;

  //  Evaluate cross section for spin 0 colour triplet.
  if (spinSave == 0) {
    sigSum = 0.5 * ( 7. / 48. + 3. * pow2(uHavg - tHavg) / (16. * sH2) )
      * ( 1. + 2. * s34Avg * tHavg / pow2(tHavg - s34Avg)
      + 2. * s34Avg * uHavg / pow2(uHavg - s34Avg)
      + 4. * pow2(s34Avg) / ((tHavg - s34Avg) * (uHavg - s34Avg)) );

    // Equal probability for two possible colour flows.
    sigTS = 0.5 * sigSum;
    sigUS = sigTS;
  }

  //  Evaluate cross section for spin 1/2 colour triplet.
  else if (spinSave == 1) {
    double tumHQ = tHQ * uHQ - s34Avg * sH;
    sigTS = ( uHQ / tHQ - 2.25 * uHQ2 / sH2 + 4.5 * s34Avg * tumHQ
      / ( sH * tHQ2) + 0.5 * s34Avg * (tHQ + s34Avg) / tHQ2
      - s34Avg*s34Avg / (sH * tHQ) ) / 6.;
    sigUS = ( tHQ / uHQ - 2.25 * tHQ2 / sH2 + 4.5 * s34Avg * tumHQ
      / ( sH * uHQ2) + 0.5 * s34Avg * (uHQ + s34Avg) / uHQ2
      - s34Avg*s34Avg / (sH * uHQ) ) / 6.;
    sigSum = sigTS + sigUS;
  }

  //  Evaluate cross section for spin 1 colour triplet.
  else {
    double tmu      = tHavg - uHavg;
    double s34Pos   = s34Avg / sH;
    double s34Pos2  = s34Pos * s34Pos;
    double s34Neg   = sH / s34Avg;
    double s34Neg2  = s34Neg * s34Neg;
    sigSum = pow2(tmu) * sH2 * (241./1536. - 1./32. * s34Pos
        + 9./16. * s34Pos2)
      + pow4(tmu) * (37./512. + 9./64. * s34Pos)
      + pow6(tmu) * (9./512. / sH2)
      + sH2 * sH2 * (133./1536. - 7./64. * s34Pos + 7./16. * s34Pos2);

    // Anomalous coupling.
    if (hasKappa)
    sigSum += pow2(tmu) * sH2 * (kappam1 * (143./384. - 7./3072 * s34Neg)
      + pow2(kappam1) * (- 1./768. * s34Neg + 185./768.)
      + pow3(kappam1) * (- 7./3072. *  s34Neg2
        - 25./3072. * s34Neg + 67./1536.)
      + pow4(kappam1) * (- 37./49152. * s34Neg2
        - 25./6144. * s34Neg + 5./1536.) )
      + pow4(tmu) * (kappam1 * 3./32.
      + pow2(kappam1) * (7./6144. * s34Neg2 - 7./768. * s34Neg + 3./128.)
      + pow3(kappam1) * (7./6144. * s34Neg2 - 7./1536. * s34Neg)
      + pow4(kappam1) * (- 1./49152. * s34Neg2 + 5./6144. * s34Neg) )
      + pow6(tmu) * pow4(kappam1) * 13./49152. / pow2(s34Avg)
      + sH2 * sH2 * ( kappam1 * 77./384.
      + pow2(kappam1) * (7./6144. * s34Neg2 + 1./96.* s34Neg + 39./256.)
      + pow3(kappam1) * (7./6144. * s34Neg2 + 13./1024. * s34Neg + 61./1536.)
      + pow4(kappam1) * (25./49152. * s34Neg2 + 5./1536. * s34Neg + 1./512.)
      );

    // Equal probability for two possible colour flows.
    sigSum /= pow2( (uHavg-s34Avg) * (tHavg-s34Avg) );
    sigTS = 0.5 * sigSum;
    sigUS = sigTS;
  }

  // Final answer, with common factors.
  sigma = (M_PI / sH2) * pow2(alpS) * sigSum * nCHV * openFracPair;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2qGqGbar::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, idNew, -idNew);

  // Two colour flow topologies.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else                 setColAcol( 1, 2, 3, 1, 3, 0, 0, 2);

}

//==========================================================================

// Sigma2qqbar2qGqGbar class.
// Cross section for q qbar -> qG qGbar (generic quark of spin 0, 1/2 or 1).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qqbar2qGqGbar::initProc() {

  // Number of colours. Coupling kappa used for vector state.
  nCHV         = settingsPtr->mode("HiddenValley:Ngauge");
  kappa        = settingsPtr->parm("HiddenValley:kappa");

   // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(idNew, -idNew);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2qqbar2qGqGbar::sigmaKin() {

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  double delta  = 0.25 * pow2(s3 - s4) / sH;
  double s34Avg = 0.5 * (s3 + s4) - delta;
  double tHavg  = tH - delta;
  double uHavg  = uH - delta;
  double tHQ    = -0.5 * (sH - tH + uH);
  double uHQ    = -0.5 * (sH + tH - uH);
  double tHQ2   = tHQ * tHQ;
  double uHQ2   = uHQ * uHQ;

  //  Evaluate cross section for spin 0 colour triplet.
  if (spinSave == 0) {
    sigSum = (1./9.) * (sH * (sH - 4. * s34Avg)
      - pow2(uHavg - tHavg)) / sH2;
  }

  //  Evaluate cross section for spin 1/2 colour triplet.
  else if (spinSave == 1) {
    sigSum = (4./9.) * ((tHQ2 + uHQ2) / sH2 + 2. * s34Avg / sH);
  }

  //  Evaluate cross section for spin 1 colour triplet.
  else {
    double tuH34 = (tHavg + uHavg) / s34Avg;
    sigSum = (1./9.) * (
      pow2(1. + kappa) * sH * s34Avg * (pow2(tuH34) - 4.)
      + (tHavg * uHavg - pow2(s34Avg)) * (8. + 2. * (1. - pow2(kappa)) * tuH34
      + pow2(kappa) * pow2(tuH34)) ) / sH2;
  }

  // Final answer, with common factors.
  sigma = (M_PI / sH2) * pow2(alpS) * sigSum * nCHV * openFracPair;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2qGqGbar::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idNew, -idNew);

  // tH defined between f and qG: must swap tHat <-> uHat if qbar q in.
  swapTU = (id1 < 0);

  // Colour flow topologies.
  if (id1 > 0) setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  else         setColAcol( 0, 2, 1, 0, 1, 0, 0, 2);

}

//==========================================================================

// Sigma2ffbar2fGfGbar class.
// Cross section for f fbar -> qG qGbar (generic quark of spin 0, 1/2 or 1)
// via gamma^*/Z^* s-channel exchange. Still under development!! ??

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2fGfGbar::initProc() {

  // Charge and number of colours. Coupling kappa used for vector state.
  if (settingsPtr->flag("HiddenValley:doKinMix"))
    eQHV2      = pow2(settingsPtr->parm("HiddenValley:kinMix"));
  else
    eQHV2      = pow2( particleDataPtr->charge(idNew) );
  nCHV         = settingsPtr->mode("HiddenValley:Ngauge");
  kappa        = settingsPtr->parm("HiddenValley:kappa");

  // Coloured or uncoloured particle.
  hasColour    = (particleDataPtr->colType(idNew) != 0);
  colFac       = (hasColour) ? 3. : 1.;

   // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(idNew, -idNew);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2fGfGbar::sigmaKin() {

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  double delta  = 0.25 * pow2(s3 - s4) / sH;
  double s34Avg = 0.5 * (s3 + s4) - delta;
  double tHavg  = tH - delta;
  double uHavg  = uH - delta;
  double tHQ    = -0.5 * (sH - tH + uH);
  double uHQ    = -0.5 * (sH + tH - uH);
  double tHQ2   = tHQ * tHQ;
  double uHQ2   = uHQ * uHQ;

  //  Evaluate cross section for spin 0 colour triplet.
  if (spinSave == 0) {
    sigSum = 0.5 * (sH * (sH - 4. * s34Avg) - pow2(uHavg - tHavg)) / sH2;
  }

  //  Evaluate cross section for spin 1/2 colour triplet.
  else if (spinSave == 1) {
    sigSum = 2. * ((tHQ2 + uHQ2) / sH2 + 2. * s34Avg / sH);
  }

  //  Evaluate cross section for spin 1 colour triplet.
  else {
    double tuH34 = (tHavg + uHavg) / s34Avg;
    sigSum = 0.5 * ( pow2(1. + kappa) * sH * s34Avg * (pow2(tuH34) - 4.)
      + (tHavg * uHavg - pow2(s34Avg)) * (8. + 2. * (1. - pow2(kappa)) * tuH34
      + pow2(kappa) * pow2(tuH34)) ) / sH2;
  }

  // Final-state charge factors.
  sigSum *= colFac * eQHV2 * (1. + alpS / M_PI);

  // Final answer, except for initial-state weight
  sigma0 = (M_PI / sH2) * pow2(alpEM) * sigSum * nCHV * openFracPair;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2fGfGbar::sigmaHat() {

  // Charge and colour factors.
  double eNow  = couplingsPtr->ef( abs(id1) );
  double sigma = sigma0 * pow2(eNow);
  if (abs(id1) < 9) sigma /= 3.;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2fGfGbar::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idNew, -idNew);

  // tH defined between f and qG: must swap tHat <-> uHat if fbar f in.
  swapTU = (id1 < 0);

  // Colour flow topologies.
  if (hasColour) {
    if (id1 > 0 && id1 < 7)       setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
    else if (id1 > -7 && id1 < 0) setColAcol( 0, 1, 1, 0, 2, 0, 0, 2);
    else                          setColAcol( 0, 0, 0, 0, 1, 0, 0, 1);
  } else {
    if (id1 > 0 && id1 < 7)       setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
    else if (id1 > -7 && id1 < 0) setColAcol( 0, 1, 1, 0, 0, 0, 0, 0);
    else                          setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  }

}

//==========================================================================

// Sigma1ffbar2Zv class.
// Cross section for f fbar -> Zv, where Zv couples both to the SM and
// to a hidden sector. Primitive coupling structure.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ffbar2Zv::initProc() {

  // Store Zv mass and width for propagator.
  idZv     = 4900023;
  mRes     = particleDataPtr->m0(idZv);
  GammaRes = particleDataPtr->mWidth(idZv);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(idZv);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat); first step when inflavours unknown.

void Sigma1ffbar2Zv::sigmaKin() {

  // Breit-Wigner, including some (guessed) spin factors.
  double sigBW    = 12. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Outgoing width: only includes channels left open.
  double widthOut = particlePtr->resWidthOpen(663, mH);

  // Temporary answer.
  sigOut = sigBW * widthOut;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat); second step when inflavours known.

double Sigma1ffbar2Zv::sigmaHat() {

  // Incoming quark or lepton; for former need two 1/3 colour factors.
  int id1Abs     = abs(id1);
  double widthIn = particlePtr->resWidthChan( mH, id1Abs, -id1Abs);
  if (id1Abs < 6) widthIn /= 9.;
  return widthIn * sigOut;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2Zv::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idZv);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 6) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma1ffbar2Zv::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Identity of mother of decaying resonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For Z' itself angular distribution as if gamma*.
  if (iResBeg == 5 && iResEnd == 5) {
    double mr     = 4. * pow2(process[6].m()) / sH;
    double cosThe = (process[3].p() - process[4].p())
      * (process[7].p() - process[6].p()) / (sH * sqrtpos(1. - mr));
    double wt     = 1. + pow2(cosThe) + mr * (1. - pow2(cosThe));
    return 0.5 * wt;
  }

  // For top decay hand over to standard routine.
  if (idMother == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // Else done.
  return 1.;

}

//==========================================================================

} // end namespace Pythia8
