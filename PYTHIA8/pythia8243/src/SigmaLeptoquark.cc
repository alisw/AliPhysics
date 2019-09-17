// SigmaLeptoquark.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// leptoquark simulation classes.

#include "Pythia8/SigmaLeptoquark.h"

namespace Pythia8 {

//==========================================================================

// Sigma1ql2LeptoQuark class.
// Cross section for q l -> LQ (leptoquark state).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ql2LeptoQuark::initProc() {

  // Store LQ mass and width for propagator.
  mRes     = particleDataPtr->m0(42);
  GammaRes = particleDataPtr->mWidth(42);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Yukawa coupling strength.
  kCoup    = settingsPtr->parm("LeptoQuark:kCoup");

  // Set pointer to particle properties and decay table.
  LQPtr    = particleDataPtr->particleDataEntryPtr(42);

  // Read out quark and lepton the LQ couples to.
  idQuark  = LQPtr->channel(0).product(0);
  idLepton = LQPtr->channel(0).product(1);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1ql2LeptoQuark::sigmaKin() {

  // Incoming width for correct quark-lepton pair.
  widthIn  = 0.25 * alpEM * kCoup * mH;

  // Set up Breit-Wigner.
  sigBW    = 4. * M_PI/ ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat) for specific incoming flavours.

double Sigma1ql2LeptoQuark::sigmaHat() {

  // Identify whether correct incoming flavours.
  int idLQ = 0;
  if      (id1 ==  idQuark && id2 ==  idLepton) idLQ =  42;
  else if (id2 ==  idQuark && id1 ==  idLepton) idLQ =  42;
  else if (id1 == -idQuark && id2 == -idLepton) idLQ = -42;
  else if (id2 == -idQuark && id1 == -idLepton) idLQ = -42;
  if (idLQ == 0) return 0.;

  // Outgoing width and total sigma. Done.
  return widthIn * sigBW * LQPtr->resWidthOpen(idLQ, mH);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ql2LeptoQuark::setIdColAcol() {

  // Flavours.
  int idq  = (abs(id1) < 9) ? id1 : id2;
  int idLQ = (idq > 0) ? 42 : -42;
  setId( id1, id2, idLQ);

  // Colour flow topology.
  if (id1 == idq) setColAcol( 1, 0, 0, 0, 1, 0);
  else            setColAcol( 0, 0, 1, 0, 1, 0);
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2qg2LeptoQuarkl class.
// Cross section for q g -> LQ l (leptoquark state).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qg2LeptoQuarkl::initProc() {

  // Store LQ mass and width for propagator.
  mRes     = particleDataPtr->m0(42);
  GammaRes = particleDataPtr->mWidth(42);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Yukawa coupling strength.
  kCoup    = settingsPtr->parm("LeptoQuark:kCoup");

  // Read out quark and lepton the LQ couples to.
  ParticleDataEntry* LQPtr = particleDataPtr->particleDataEntryPtr(42);
  idQuark  = LQPtr->channel(0).product(0);
  idLepton = LQPtr->channel(0).product(1);

   // Secondary open width fraction.
  openFracPos = LQPtr->resOpenFrac( 42);
  openFracNeg = LQPtr->resOpenFrac(-42);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2qg2LeptoQuarkl::sigmaKin() {

  //  Evaluate cross section.
  sigma0 = (M_PI / sH2) * kCoup * (alpS * alpEM / 6.) * (-tH / sH)
    * (uH2 + s3 * s3) / pow2(uH - s3);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat) for specific incoming flavours.

double Sigma2qg2LeptoQuarkl::sigmaHat() {

  // Check that correct incoming flavour.
  if (abs(id1) != idQuark && abs(id2) != idQuark) return 0.;

  // Answer, with secondary width correction.
  double sigma = sigma0;
  sigma *= (id1 == idQuark || id2 == idQuark) ? openFracPos : openFracNeg;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2LeptoQuarkl::setIdColAcol() {

  // Flavour set up for q g -> H q.
  int idq = (id2 == 21) ? id1 : id2;
  int idLQ = (idq > 0) ? 42 : -42;
  int idlp = (idq > 0) ? -idLepton : idLepton;
  setId( id1, id2, idLQ, idlp);

  // tH defined between f and LQ: must swap tHat <-> uHat if g q in.
  swapTU = (id1 == 21);

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 2, 0, 0, 0);
  else           setColAcol( 2, 1, 1, 0, 2, 0, 0, 0);
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2gg2LQLQbar class.
// Cross section for g g -> LQ LQbar (leptoquark state).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2LQLQbar::initProc() {

  // Store LQ mass and width for propagator.
  mRes     = particleDataPtr->m0(42);
  GammaRes = particleDataPtr->mWidth(42);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

   // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(42, -42);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2gg2LQLQbar::sigmaKin() {

  // Average outgoing masses and adjust kinematics accordingly.
  double delta = 0.25 * pow2(s3 - s4) / sH;
  double m2avg = 0.5 * (s3 + s4) - delta;
  double tHavg = tH - delta;
  double uHavg = uH - delta;

  //  Evaluate cross section. Secondary width for G*.
  sigma = (M_PI / sH2) * 0.5 * pow2(alpS)
    * ( 7. / 48. + 3. * pow2(uHavg - tHavg) / (16. * sH2) )
    * ( 1. + 2. * m2avg * tHavg / pow2(tHavg - m2avg)
    + 2. * m2avg * uHavg / pow2(uHavg - m2avg)
    + 4. * m2avg * m2avg / ((tHavg - m2avg) * (uHavg - m2avg)) );
  sigma *= openFrac;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2LQLQbar::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, 42, -42);

  // Colour flow topologies: random choice between two mirrors.
  if (rndmPtr->flat() < 0.5) setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else                       setColAcol( 1, 2, 3, 1, 3, 0, 0, 2);

}

//==========================================================================

// Sigma2qqbar2LQLQbar class.
// Cross section for q qbar -> LQ LQbar (leptoquark state).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qqbar2LQLQbar::initProc() {

  // Store LQ mass and width for propagator.
  mRes     = particleDataPtr->m0(42);
  GammaRes = particleDataPtr->mWidth(42);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Yukawa coupling strength.
  kCoup    = settingsPtr->parm("LeptoQuark:kCoup");

  // Read out quark and lepton the LQ couples to.
  ParticleDataEntry* LQPtr = particleDataPtr->particleDataEntryPtr(42);
  idQuark  = LQPtr->channel(0).product(0);

   // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(42, -42);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2qqbar2LQLQbar::sigmaKin() {

  // Average outgoing masses and adjust kinematics accordingly.
  double delta = 0.25 * pow2(s3 - s4) / sH;
  double m2avg = 0.5 * (s3 + s4) - delta;
  double tHavg = tH - delta;
  double uHavg = uH - delta;

  // Evaluate cross section for quark of different flavour than LQ.
  sigmaDiff = (M_PI / sH2) * (pow2(alpS) / 9.)
    * ( sH * (sH - 4. * m2avg) - pow2(uHavg - tHavg) ) / sH2;

  // Evaluate cross section for quark of same flavour as LQ.
  sigmaSame = sigmaDiff + (M_PI / sH2) * (pow2(kCoup * alpEM) / 8.)
    * (-sH * tHavg - pow2(m2avg-tHavg)) / pow2(tHavg)
    + (M_PI / sH2) * (kCoup * alpEM * alpS / 18.) * ( (m2avg - tHavg)
    * (uHavg - tHavg) + sH * (m2avg + tHavg) ) / (sH * tHavg);

  // Open fraction.
  sigmaDiff *= openFrac;
  sigmaSame *= openFrac;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2LQLQbar::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 42, -42);

  // tH defined between f and LQ: must swap tHat <-> uHat if qbar q in.
  swapTU = (id1 < 0);

  // Colour flow topologies.
  if (id1 > 0) setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  else         setColAcol( 0, 2, 1, 0, 1, 0, 0, 2);

}

//==========================================================================

} // end namespace Pythia8
