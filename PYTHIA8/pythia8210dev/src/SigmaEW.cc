// SigmaEW.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// electroweak simulation classes (including photon processes).

#include "Pythia8/SigmaEW.h"

namespace Pythia8 {

//==========================================================================

// Sigma2qg2qgamma class.
// Cross section for q g -> q gamma.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2qg2qgamma::sigmaKin() {

  // Calculate kinematics dependence.
  sigUS  = (1./3.) * (sH2 + uH2) / (-sH * uH);

  // Answer.
  sigma0 =  (M_PI/sH2) * alpS * alpEM * sigUS;

  }

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2qg2qgamma::sigmaHat() {

  // Incoming flavour gives charge factor.
  int idNow    = (id2 == 21) ? id1 : id2;
  double eNow  = couplingsPtr->ef( abs(idNow) );
  return sigma0 * pow2(eNow);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2qgamma::setIdColAcol() {

  // Construct outgoing flavours.
  id3 = (id1 == 21) ? 22 : id1;
  id4 = (id2 == 21) ? 22 : id2;
  setId( id1, id2, id3, id4);

  // Colour flow topology. Swap if first is gluon, or when antiquark.
  setColAcol( 1, 0, 2, 1, 2, 0, 0, 0);
  if (id1 == 21) swapCol1234();
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2ggamma class.
// Cross section for q qbar -> g gamma.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2qqbar2ggamma::sigmaKin() {

  // Calculate kinematics dependence.
  double sigTU = (8./9.) * (tH2 + uH2) / (tH * uH);

  // Answer.
  sigma0 = (M_PI/sH2) * alpS * alpEM * sigTU;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2qqbar2ggamma::sigmaHat() {

  // Incoming flavour gives charge factor.
  double eNow  = couplingsPtr->ef( abs(id1) );
  return sigma0 * pow2(eNow);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2ggamma::setIdColAcol() {

  // Outgoing flavours trivial.
  setId( id1, id2, 21, 22);

  // One colour flow topology. Swap if first is antiquark.
  setColAcol( 1, 0, 0, 2, 1, 2, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2gg2ggamma class.
// Cross section for g g -> g gamma.
// Proceeds through a quark box, by default using 5 massless quarks.

//--------------------------------------------------------------------------

// Initialize process, especially parton-flux object.

void Sigma2gg2ggamma::initProc() {

  // Maximum quark flavour in loop.
  int nQuarkLoop = settingsPtr->mode("PromptPhoton:nQuarkLoop");

  // Calculate charge factor from the allowed quarks in the box.
  chargeSum                       = - 1./3. + 2./3. - 1./3.;
  if (nQuarkLoop >= 4) chargeSum += 2./3.;
  if (nQuarkLoop >= 5) chargeSum -= 1./3.;
  if (nQuarkLoop >= 6) chargeSum += 2./3.;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat).

void Sigma2gg2ggamma::sigmaKin() {

  // Logarithms of Mandelstam variable ratios.
  double logST = log( -sH / tH );
  double logSU = log( -sH / uH );
  double logTU = log(  tH / uH );

  // Real and imaginary parts of separate amplitudes.
  double b0stuRe = 1. + (tH - uH) / sH * logTU
    + 0.5 * (tH2 + uH2) / sH2 * (pow2(logTU) + pow2(M_PI));
  double b0stuIm = 0.;
  double b0tsuRe = 1. + (sH - uH) / tH * logSU
    + 0.5 * (sH2 + uH2) / tH2 * pow2(logSU);
  double b0tsuIm = -M_PI * ( (sH - uH) / tH + (sH2 + uH2) / tH2 * logSU);
  double b0utsRe = 1. + (sH - tH) / uH * logST
    + 0.5 * (sH2 + tH2) / uH2 * pow2(logST);
  double b0utsIm = -M_PI * ( (sH - tH) / uH + (sH2 + tH2) / uH2 * logST);
  double b1stuRe = -1.;
  double b1stuIm = 0.;
  double b2stuRe = -1.;
  double b2stuIm = 0.;

  // Calculate kinematics dependence.
  double sigBox = pow2(b0stuRe) + pow2(b0stuIm) + pow2(b0tsuRe)
    + pow2(b0tsuIm) + pow2(b0utsRe) + pow2(b0utsIm) + 4. * pow2(b1stuRe)
    + 4. * pow2(b1stuIm) + pow2(b2stuRe) + pow2(b2stuIm);

  // Answer.
  sigma = (5. / (192. * M_PI * sH2)) * pow2(chargeSum)
    * pow3(alpS) * alpEM * sigBox;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2ggamma::setIdColAcol() {

  // Flavours and colours are trivial.
  setId( id1, id2, 21, 22);
  setColAcol( 1, 2, 2, 3, 1, 3, 0, 0);
  if (rndmPtr->flat() > 0.5) swapColAcol();

}

//==========================================================================

// Sigma2ffbar2gammagamma class.
// Cross section for q qbar -> gamma gamma.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2gammagamma::sigmaKin() {

  // Calculate kinematics dependence.
  sigTU  = 2. * (tH2 + uH2) / (tH * uH);

  // Answer contains factor 1/2 from identical photons.
  sigma0 = (M_PI/sH2) * pow2(alpEM) * 0.5 * sigTU;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2gammagamma::sigmaHat() {

  // Incoming flavour gives charge and colour factors.
  double eNow   = couplingsPtr->ef( abs(id1) );
  double colFac = (abs(id1) < 9) ? 1. / 3. : 1.;
  return  sigma0 * pow4(eNow) * colFac;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2gammagamma::setIdColAcol() {

  // Outgoing flavours trivial.
  setId( id1, id2, 22, 22);

  // No colours at all or one flow topology. Swap if first is antiquark.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2gg2gammagamma class.
// Cross section for g g -> gamma gamma.
// Proceeds through a quark box, by default using 5 massless quarks.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2gammagamma::initProc() {

  // Maximum quark flavour in loop.
  int nQuarkLoop = settingsPtr->mode("PromptPhoton:nQuarkLoop");

  // Calculate charge factor from the allowed quarks in the box.
  charge2Sum                       = 1./9. + 4./9. + 1./9.;
  if (nQuarkLoop >= 4) charge2Sum += 4./9.;
  if (nQuarkLoop >= 5) charge2Sum += 1./9.;
  if (nQuarkLoop >= 6) charge2Sum += 4./9.;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat).

void Sigma2gg2gammagamma::sigmaKin() {

  // Logarithms of Mandelstam variable ratios.
  double logST = log( -sH / tH );
  double logSU = log( -sH / uH );
  double logTU = log(  tH / uH );

  // Real and imaginary parts of separate amplitudes.
  double b0stuRe = 1. + (tH - uH) / sH * logTU
    + 0.5 * (tH2 + uH2) / sH2 * (pow2(logTU) + pow2(M_PI));
  double b0stuIm = 0.;
  double b0tsuRe = 1. + (sH - uH) / tH * logSU
    + 0.5 * (sH2 + uH2) / tH2 * pow2(logSU);
  double b0tsuIm = -M_PI * ( (sH - uH) / tH + (sH2 + uH2) / tH2 * logSU);
  double b0utsRe = 1. + (sH - tH) / uH * logST
    + 0.5 * (sH2 + tH2) / uH2 * pow2(logST);
  double b0utsIm = -M_PI * ( (sH - tH) / uH + (sH2 + tH2) / uH2 * logST);
  double b1stuRe = -1.;
  double b1stuIm = 0.;
  double b2stuRe = -1.;
  double b2stuIm = 0.;

  // Calculate kinematics dependence.
  double sigBox = pow2(b0stuRe) + pow2(b0stuIm) + pow2(b0tsuRe)
    + pow2(b0tsuIm) + pow2(b0utsRe) + pow2(b0utsIm) + 4. * pow2(b1stuRe)
    + 4. * pow2(b1stuIm) + pow2(b2stuRe) + pow2(b2stuIm);

  // Answer contains factor 1/2 from identical photons.
  sigma = (0.5 / (16. * M_PI * sH2)) * pow2(charge2Sum)
    * pow2(alpS) * pow2(alpEM) * sigBox;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2gammagamma::setIdColAcol() {

  // Flavours and colours are trivial.
  setId( id1, id2, 22, 22);
  setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

}

//==========================================================================

// Sigma2ff2fftgmZ class.
// Cross section for f f' -> f f' via t-channel gamma*/Z0 exchange
// (f is quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ff2fftgmZ::initProc() {

  // Store Z0 mass for propagator. Common coupling factor.
  gmZmode   = settingsPtr->mode("WeakZ0:gmZmode");
  mZ        = particleDataPtr->m0(23);
  mZS       = mZ*mZ;
  thetaWRat = 1. / (16. * couplingsPtr->sin2thetaW()
            * couplingsPtr->cos2thetaW());

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ff2fftgmZ::sigmaKin() {

  // Cross section part common for all incoming flavours.
  double sigma0 = (M_PI / sH2) * pow2(alpEM);

  // Kinematical functions for gamma-gamma, gamma-Z and Z-Z parts.
  sigmagmgm = sigma0 * 2. * (sH2 + uH2) / tH2;
  sigmagmZ  = sigma0 * 4. * thetaWRat * sH2 / (tH * (tH - mZS));
  sigmaZZ   = sigma0 * 2. * pow2(thetaWRat) * sH2 / pow2(tH - mZS);
  if (gmZmode == 1) {sigmagmZ = 0.; sigmaZZ = 0.;}
  if (gmZmode == 2) {sigmagmgm = 0.; sigmagmZ = 0.;}

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ff2fftgmZ::sigmaHat() {

  // Couplings for current flavour combination.
  int id1Abs = abs(id1);
  double  e1 = couplingsPtr->ef(id1Abs);
  double  v1 = couplingsPtr->vf(id1Abs);
  double  a1 = couplingsPtr->af(id1Abs);
  int id2Abs = abs(id2);
  double  e2 = couplingsPtr->ef(id2Abs);
  double  v2 = couplingsPtr->vf(id2Abs);
  double  a2 = couplingsPtr->af(id2Abs);

  // Distinguish same-sign and opposite-sign fermions.
  double epsi = (id1 * id2 > 0) ? 1. : -1.;

  // Flavour-dependent cross section.
  double sigma = sigmagmgm * pow2(e1 * e2)
    + sigmagmZ * e1 * e2 * (v1 * v2 * (1. + uH2 / sH2)
      + a1 * a2 * epsi * (1. - uH2 / sH2))
    + sigmaZZ * ((v1*v1 + a1*a1) * (v2*v2 + a2*a2) * (1. + uH2 / sH2)
      + 4. * v1 * a1 * v2 * a2 * epsi * (1. - uH2 / sH2));

  // Spin-state extra factor 2 per incoming neutrino.
  if (id1Abs == 12 || id1Abs == 14 || id1Abs == 16) sigma *= 2.;
  if (id2Abs == 12 || id2Abs == 14 || id2Abs == 16) sigma *= 2.;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ff2fftgmZ::setIdColAcol() {

  // Trivial flavours: out = in.
  setId( id1, id2, id1, id2);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9 && abs(id2) < 9 && id1*id2 > 0)
                         setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  else if (abs(id1) < 9 && abs(id2) < 9)
                         setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  else if (abs(id1) < 9) setColAcol( 1, 0, 0, 0, 1, 0, 0, 0);
  else if (abs(id2) < 9) setColAcol( 0, 0, 1, 0, 0, 0, 1, 0);
  else                   setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if ( (abs(id1) < 9 && id1 < 0) || (abs(id1) > 10 && id2 < 0) )
    swapColAcol();

}

//==========================================================================

// Sigma2ff2fftW class.
// Cross section for f_1 f_2 -> f_3 f_4 via t-channel W+- exchange
// (f is quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ff2fftW::initProc() {

  // Store W+- mass for propagator. Common coupling factor.
  mW        = particleDataPtr->m0(24);
  mWS       = mW*mW;
  thetaWRat = 1. / (4. * couplingsPtr->sin2thetaW());

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ff2fftW::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * pow2(alpEM * thetaWRat)
    * 4. * sH2 / pow2(tH - mWS);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ff2fftW::sigmaHat() {

  // Some flavour combinations not possible.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  if ( (id1Abs%2 == id2Abs%2 && id1 * id2 > 0)
    || (id1Abs%2 != id2Abs%2 && id1 * id2 < 0) ) return 0.;

  // Basic cross section.
  double sigma = sigma0;
  if (id1 * id2 < 0) sigma *= uH2 / sH2;

  // CKM factors for final states.
  sigma *= couplingsPtr->V2CKMsum(id1Abs) *  couplingsPtr->V2CKMsum(id2Abs);

  // Spin-state extra factor 2 per incoming neutrino.
  if (id1Abs == 12 || id1Abs == 14 || id1Abs == 16) sigma *= 2.;
  if (id2Abs == 12 || id2Abs == 14 || id2Abs == 16) sigma *= 2.;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ff2fftW::setIdColAcol() {

  // Pick out-flavours by relative CKM weights.
  id3 = couplingsPtr->V2CKMpick(id1);
  id4 = couplingsPtr->V2CKMpick(id2);
  setId( id1, id2, id3, id4);

  // Colour flow topologies. Swap when antiquarks.
  if      (abs(id1) < 9 && abs(id2) < 9 && id1*id2 > 0)
                         setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  else if (abs(id1) < 9 && abs(id2) < 9)
                         setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  else if (abs(id1) < 9) setColAcol( 1, 0, 0, 0, 1, 0, 0, 0);
  else if (abs(id2) < 9) setColAcol( 0, 0, 1, 0, 0, 0, 1, 0);
  else                   setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if ( (abs(id1) < 9 && id1 < 0) || (abs(id1) > 10 && id2 < 0) )
    swapColAcol();

}


//==========================================================================

// Sigma2qq2QqtW class.
// Cross section for q q' -> Q q" via t-channel W+- exchange.
// Related to Sigma2ff2ffViaW class, but with massive matrix elements.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qq2QqtW::initProc() {

  // Process name.
  nameSave                 = "q q -> Q q (t-channel W+-)";
  if (idNew == 4) nameSave = "q q -> c q (t-channel W+-)";
  if (idNew == 5) nameSave = "q q -> b q (t-channel W+-)";
  if (idNew == 6) nameSave = "q q -> t q (t-channel W+-)";
  if (idNew == 7) nameSave = "q q -> b' q (t-channel W+-)";
  if (idNew == 8) nameSave = "q q -> t' q (t-channel W+-)";

  // Store W+- mass for propagator. Common coupling factor.
  mW        = particleDataPtr->m0(24);
  mWS       = mW*mW;
  thetaWRat = 1. / (4. * couplingsPtr->sin2thetaW());

  // Secondary open width fractions, relevant for top (or heavier).
  openFracPos = particleDataPtr->resOpenFrac(idNew);
  openFracNeg = particleDataPtr->resOpenFrac(-idNew);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2qq2QqtW::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * pow2(alpEM * thetaWRat) * 4. / pow2(tH - mWS);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2qq2QqtW::sigmaHat() {

  // Some flavour combinations not possible.
  int id1Abs  = abs(id1);
  int id2Abs  = abs(id2);
  bool diff12 = (id1Abs%2 != id2Abs%2);
  if ( (!diff12 && id1 * id2 > 0)
    || ( diff12 && id1 * id2 < 0) ) return 0.;

  // Basic cross section.
  double sigma = sigma0;
  sigma *= (id1 * id2 > 0) ? sH * (sH - s3) : uH * (uH - s3);

  // Secondary width if t or tbar produced on either side.
  double openFrac1 = (id1 > 0) ? openFracPos : openFracNeg;
  double openFrac2 = (id2 > 0) ? openFracPos : openFracNeg;

  // CKM factors for final states; further impossible case.
  bool diff1N = (id1Abs%2 != idNew%2);
  bool diff2N = (id2Abs%2 != idNew%2);
  if (diff1N && diff2N)
    sigma *= ( couplingsPtr->V2CKMid(id1Abs, idNew) * openFrac1
             * couplingsPtr->V2CKMsum(id2Abs) + couplingsPtr->V2CKMsum(id1Abs)
             * couplingsPtr->V2CKMid(id2Abs, idNew) * openFrac2 );
  else if (diff1N)
    sigma *= couplingsPtr->V2CKMid(id1Abs, idNew) * openFrac1
             * couplingsPtr->V2CKMsum(id2Abs);
  else if (diff2N)
    sigma *= couplingsPtr->V2CKMsum(id1Abs)
             * couplingsPtr->V2CKMid(id2Abs, idNew) * openFrac2;
  else sigma = 0.;

  // Spin-state extra factor 2 per incoming neutrino.
  if (id1Abs == 12 || id1Abs == 14 || id1Abs == 16) sigma *= 2.;
  if (id2Abs == 12 || id2Abs == 14 || id2Abs == 16) sigma *= 2.;

  // Answer.
  return  sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qq2QqtW::setIdColAcol() {

  // For topologies like d dbar -> (t/c/u) (t/c/u)bar pick side.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  int side   = 1;
  if ( (id1Abs + idNew)%2 == 1 && (id2Abs + idNew)%2 == 1 ) {
    double prob1 = couplingsPtr->V2CKMid(id1Abs, idNew)
                   * couplingsPtr->V2CKMsum(id2Abs);
    prob1       *= (id1 > 0) ? openFracPos : openFracNeg;
    double prob2 = couplingsPtr->V2CKMid(id2Abs, idNew)
                   * couplingsPtr->V2CKMsum(id1Abs);
    prob2       *= (id2 > 0) ? openFracPos : openFracNeg;
    if (prob2 > rndmPtr->flat() * (prob1 + prob2)) side = 2;
  }
  else if ((id2Abs + idNew)%2 == 1) side = 2;

  // Pick out-flavours by relative CKM weights.
  if (side == 1) {
    // q q' -> t q" : correct order from start.
    id3 = (id1 > 0) ? idNew : -idNew;
    id4 = couplingsPtr->V2CKMpick(id2);
    setId( id1, id2, id3, id4);
  } else {
    // q q' -> q" t : stored as t q" so swap tHat <-> uHat.
    swapTU = true;
    id3 = couplingsPtr->V2CKMpick(id1);
    id4 = (id2 > 0) ? idNew : -idNew;
    setId( id1, id2, id4, id3);
  }

  // Colour flow topologies. Swap when antiquarks on side 1.
  if      (side == 1 && id1 * id2 > 0) setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  else if              (id1 * id2 > 0) setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);
  else if (side == 1)                  setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  else                                 setColAcol( 1, 0, 0, 2, 0, 2, 1, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles of W in top decay.

double Sigma2qq2QqtW::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // For top decay hand over to standard routine, else done.
  if (idNew == 6 && process[process[iResBeg].mother1()].idAbs() == 6)
       return weightTopDecay( process, iResBeg, iResEnd);
  else return 1.;

}

//==========================================================================

// Sigma1ffbar2gmZ class.
// Cross section for f fbar -> gamma*/Z0 (f is quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ffbar2gmZ::initProc() {

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode     = settingsPtr->mode("WeakZ0:gmZmode");

  // Store Z0 mass and width for propagator.
  mRes        = particleDataPtr->m0(23);
  GammaRes    = particleDataPtr->mWidth(23);
  m2Res       = mRes*mRes;
  GamMRat     = GammaRes / mRes;
  thetaWRat   = 1. / (16. * couplingsPtr->sin2thetaW()
              * couplingsPtr->cos2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(23);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1ffbar2gmZ::sigmaKin() {

  // Common coupling factors.
  double colQ = 3. * (1. + alpS / M_PI);

  // Reset quantities to sum. Declare variables in loop.
  gamSum = 0.;
  intSum = 0.;
  resSum = 0.;
  int    idAbs, onMode;
  double mf, mr, psvec, psaxi, betaf, ef2, efvf, vf2af2, colf;

  // Loop over all Z0 decay channels.
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    idAbs = abs( particlePtr->channel(i).product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf = particleDataPtr->m0(idAbs);

      // Check that above threshold. Phase space.
      if (mH > 2. * mf + MASSMARGIN) {
        mr    = pow2(mf / mH);
        betaf = sqrtpos(1. - 4. * mr);
        psvec = betaf * (1. + 2. * mr);
        psaxi = pow3(betaf);

        // Combine phase space with couplings.
        ef2    = couplingsPtr->ef2(idAbs) * psvec;
        efvf   = couplingsPtr->efvf(idAbs) * psvec;
        vf2af2 = couplingsPtr->vf2(idAbs) * psvec
               + couplingsPtr->af2(idAbs) * psaxi;
        colf   = (idAbs < 6) ? colQ : 1.;

        // Store sum of combinations. For outstate only open channels.
        onMode = particlePtr->channel(i).onMode();
        if (onMode == 1 || onMode == 2) {
          gamSum += colf * ef2;
          intSum += colf * efvf;
          resSum += colf * vf2af2;
        }

      // End loop over fermions.
      }
    }
  }

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  gamProp = 4. * M_PI * pow2(alpEM) / (3. * sH);
  intProp = gamProp * 2. * thetaWRat * sH * (sH - m2Res)
          / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  resProp = gamProp * pow2(thetaWRat * sH)
          / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp = 0.; resProp = 0.;}
  if (gmZmode == 2) {gamProp = 0.; intProp = 0.;}

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma1ffbar2gmZ::sigmaHat() {

  // Combine gamma, interference and Z0 parts.
  int idAbs = abs(id1);
  double sigma = couplingsPtr->ef2(idAbs)    * gamProp * gamSum
               + couplingsPtr->efvf(idAbs)   * intProp * intSum
               + couplingsPtr->vf2af2(idAbs) * resProp * resSum;

  // Colour factor. Answer.
  if (idAbs < 9) sigma /= 3.;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2gmZ::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 23);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for gamma*/Z0 decay angle.

double Sigma1ffbar2gmZ::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Z should sit in entry 5.
  if (iResBeg != 5 || iResEnd != 5) return 1.;

  // Couplings for in- and out-flavours.
  int idInAbs  = process[3].idAbs();
  double ei    = couplingsPtr->ef(idInAbs);
  double vi    = couplingsPtr->vf(idInAbs);
  double ai    = couplingsPtr->af(idInAbs);
  int idOutAbs = process[6].idAbs();
  double ef    = couplingsPtr->ef(idOutAbs);
  double vf    = couplingsPtr->vf(idOutAbs);
  double af    = couplingsPtr->af(idOutAbs);

  // Phase space factors. (One power of beta left out in formulae.)
  double mf    = process[6].m();
  double mr    = mf*mf / sH;
  double betaf = sqrtpos(1. - 4. * mr);

  // Coefficients of angular expression.
  double coefTran = ei*ei * gamProp * ef*ef + ei * vi * intProp * ef * vf
    + (vi*vi + ai*ai) * resProp * (vf*vf + pow2(betaf) * af*af);
  double coefLong = 4. * mr * ( ei*ei * gamProp * ef*ef
    + ei * vi * intProp * ef * vf + (vi*vi + ai*ai) * resProp * vf*vf );
  double coefAsym = betaf * ( ei * ai * intProp * ef * af
    + 4. * vi * ai * resProp * vf * af );

  // Flip asymmetry for in-fermion + out-antifermion.
  if (process[3].id() * process[6].id() < 0) coefAsym = -coefAsym;

  // Reconstruct decay angle and weight for it.
  double cosThe = (process[3].p() - process[4].p())
    * (process[7].p() - process[6].p()) / (sH * betaf);
  double wtMax = 2. * (coefTran + abs(coefAsym));
  double wt    = coefTran * (1. + pow2(cosThe))
     + coefLong * (1. - pow2(cosThe)) + 2. * coefAsym * cosThe;

  // Done.
  return (wt / wtMax);

}

//==========================================================================

// Sigma1ffbar2W class.
// Cross section for f fbar' -> W+- (f is quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ffbar2W::initProc() {

  // Store W+- mass and width for propagator.
  mRes     = particleDataPtr->m0(24);
  GammaRes = particleDataPtr->mWidth(24);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;
  thetaWRat = 1. / (12. * couplingsPtr->sin2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(24);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1ffbar2W::sigmaKin() {

  // Set up Breit-Wigner. Cross section for W+ and W- separately.
  double sigBW  = 12. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double preFac = alpEM * thetaWRat * mH;
  sigma0Pos     = preFac * sigBW * particlePtr->resWidthOpen(24, mH);
  sigma0Neg     = preFac * sigBW * particlePtr->resWidthOpen(-24, mH);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma1ffbar2W::sigmaHat() {

  // Secondary width for W+ or W-. CKM and colour factors.
  int idUp = (abs(id1)%2 == 0) ? id1 : id2;
  double sigma = (idUp > 0) ? sigma0Pos : sigma0Neg;
  if (abs(id1) < 9) sigma *= couplingsPtr->V2CKMid(abs(id1), abs(id2)) / 3.;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2W::setIdColAcol() {

  // Sign of outgoing W.
  int sign          = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 24 * sign);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for W decay angle.

double Sigma1ffbar2W::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // W should sit in entry 5.
  if (iResBeg != 5 || iResEnd != 5) return 1.;

  // Phase space factors.
  double mr1    = pow2(process[6].m()) / sH;
  double mr2    = pow2(process[7].m()) / sH;
  double betaf  = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2);

  // Sign of asymmetry.
  double eps    = (process[3].id() * process[6].id() > 0) ? 1. : -1.;

  // Reconstruct decay angle and weight for it.
  double cosThe = (process[3].p() - process[4].p())
    * (process[7].p() - process[6].p()) / (sH * betaf);
  double wtMax  = 4.;
  double wt     = pow2(1. + betaf * eps * cosThe) - pow2(mr1 - mr2);

  // Done.
  return (wt / wtMax);

}

//==========================================================================

// Sigma2ffbar2ffbarsgm class.
// Cross section f fbar -> gamma* -> f' fbar', for multiparton interactions.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2ffbarsgm::sigmaKin() {

  // Pick new flavour. Allow three leptons and five quarks.
  double colQ     = 1. + (alpS / M_PI);
  double flavWt   = 3. + colQ * 11. / 3.;
  double flavRndm = rndmPtr->flat() * flavWt;
  if (flavRndm < 3.) {
    if      (flavRndm < 1.) idNew = 11;
    else if (flavRndm < 2.) idNew = 13;
    else                    idNew = 15;
  } else {
    flavRndm = 3. * (flavRndm - 3.) / colQ;
    if      (flavRndm <  4.) idNew = 2;
    else if (flavRndm <  8.) idNew = 4;
    else if (flavRndm <  9.) idNew = 1;
    else if (flavRndm < 10.) idNew = 3;
    else                     idNew = 5;
  }
  double mNew  = particleDataPtr->m0(idNew);
  double m2New = mNew*mNew;

  // Calculate kinematics dependence. Give correct mass factors for
  // tHat, uHat defined as if massless kinematics, d(sigma)/d(Omega)
  // = beta (1 + cos^2(theta) + (1 - beta^2) sin^2(theta)).
  // Special case related to phase space form in multiparton interactions.
  double sigS = 0.;
  if (sH > 4. * m2New) {
    double beta = sqrt(1. - 4. * m2New / sH);
    sigS = beta * (2.* (tH2 + uH2) + 4. * (1. - beta * beta) * tH * uH)
      / sH2;
  }

  // Answer is proportional to number of outgoing flavours.
  sigma0 = (M_PI/sH2) * pow2(alpEM) * sigS * flavWt;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2ffbarsgm::sigmaHat() {

  // Charge and colour factors.
  double eNow  = couplingsPtr->ef( abs(id1) );
  double sigma = sigma0 * pow2(eNow);
  if (abs(id1) < 9) sigma /= 3.;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2ffbarsgm::setIdColAcol() {

  // Set outgoing flavours.
  id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, id2, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9 && idNew < 9) setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  else if (abs(id1) < 9)         setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else if (idNew < 9)            setColAcol( 0, 0, 0, 0, 1, 0, 0, 1);
  else                           setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2ffbar2ffbarsgmZ class.
// Cross section f fbar -> gamma*/Z0 -> f' fbar',
// i.e. gamma*/Z0 decay as part of the hard process.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2ffbarsgmZ::initProc() {

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode     = settingsPtr->mode("WeakZ0:gmZmode");

  // Store Z0 mass and width for propagator.
  mRes        = particleDataPtr->m0(23);
  GammaRes    = particleDataPtr->mWidth(23);
  m2Res       = mRes*mRes;
  GamMRat     = GammaRes / mRes;
  thetaWRat   = 1. / (16. * couplingsPtr->sin2thetaW()
              * couplingsPtr->cos2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(23);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2ffbarsgmZ::sigmaKin() {

  // Common coupling factor.
  colQ = 3. * (1. + alpS / M_PI);

  // Reset vectors and sums. Declare variables in loop.
  idVec.resize(0);
  gamT.resize(0);
  gamL.resize(0);
  intT.resize(0);
  intL.resize(0);
  intA.resize(0);
  resT.resize(0);
  resL.resize(0);
  resA.resize(0);
  gamSumT = 0.;
  gamSumL = 0.;
  intSumT = 0.;
  intSumL = 0.;
  intSumA = 0.;
  resSumT = 0.;
  resSumL = 0.;
  resSumA = 0.;
  int    onMode, idAbs;
  double mf, mr, betaf, ef, vf, af, colf, gamTf, gamLf, intTf, intLf,
         intAf, resTf, resLf, resAf;

  // Loop over all Z0 decay channels.
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    onMode = particlePtr->channel(i).onMode();
    idAbs = abs( particlePtr->channel(i).product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (onMode == 1 || onMode == 2) && ((idAbs > 0 && idAbs < 6)
      || ( idAbs > 10 && idAbs < 17)) ) {
      mf = particleDataPtr->m0(idAbs);

      // Check that channel above threshold. Phase space.
      if (mH > 2. * mf + MASSMARGIN) {
        mr    = pow2(mf / mH);
        betaf = sqrtpos(1. - 4. * mr);

        // Combine couplings (including colour) with phase space.
        ef    = couplingsPtr->ef(idAbs);
        vf    = couplingsPtr->vf(idAbs);
        af    = couplingsPtr->af(idAbs);
        colf  = (idAbs < 6) ? colQ : 1.;
        gamTf = colf * ef * ef * betaf;
        gamLf = gamTf * 4. * mr;
        intTf = colf * ef * vf * betaf;
        intLf = intTf * 4. * mr;
        intAf = colf * ef * af * betaf;
        resTf = colf * (vf * vf * betaf + af * af * pow3(betaf));
        resLf = colf * vf * vf * betaf * 4. * mr;
        resAf = colf * vf * af * betaf * 4.;

        // Store individual coplings and their sums.
        idVec.push_back(idAbs);
        gamT.push_back(gamTf);
        gamL.push_back(gamLf);
        intT.push_back(intTf);
        intL.push_back(intLf);
        intA.push_back(intAf);
        resT.push_back(resTf);
        resL.push_back(resLf);
        resA.push_back(resAf);
        gamSumT += gamTf;
        gamSumL += gamLf;
        intSumT += intTf;
        intSumL += intLf;
        intSumA += intAf;
        resSumT += resTf;
        resSumL += resLf;
        resSumA += resAf;

      // End loop over Z0 decay channels.
      }
    }
  }

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  gamProp = M_PI * pow2(alpEM) / sH2;
  intProp = gamProp * 2. * thetaWRat * sH * (sH - m2Res)
          / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  resProp = gamProp * pow2(thetaWRat * sH)
          / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp = 0.; resProp = 0.;}
  if (gmZmode == 2) {gamProp = 0.; intProp = 0.;}

  // Scattering angle in subsystem rest frame.
  cThe = (tH - uH) / sH;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2ffbarsgmZ::sigmaHat() {

  // Couplings for current in-flavour.
  int id1Abs = abs(id1);
  double ei  = couplingsPtr->ef(id1Abs);
  double vi  = couplingsPtr->vf(id1Abs);
  double ai  = couplingsPtr->af(id1Abs);

  // Coefficients of angular expression.
  double coefT = ei*ei * gamProp * gamSumT + ei*vi * intProp * intSumT
               + (vi*vi + ai*ai) * resProp * resSumT;
  double coefL = ei*ei * gamProp * gamSumL + ei*vi * intProp * intSumL
               + (vi*vi + ai*ai) * resProp * resSumL;
  double coefA = ei*ai * intProp * intSumA + vi*ai * resProp * resSumA;

  // Colour factor. Answer.
  double sigma = coefT * (1. + pow2(cThe)) + coefL * (1. - pow2(cThe))
               + 2. * coefA * cThe;
  if (id1Abs < 9) sigma /= 3.;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2ffbarsgmZ::setIdColAcol() {

  // Couplings for chosen in-flavour.
  int id1Abs = abs(id1);
  double ei  = couplingsPtr->ef(id1Abs);
  double vi  = couplingsPtr->vf(id1Abs);
  double ai  = couplingsPtr->af(id1Abs);

  // Contribution from each allowed out-flavour.
  sigTLA.resize(0);
  for (int i = 0; i < int(idVec.size()); ++i) {
    double coefT = ei*ei * gamProp * gamT[i] + ei*vi * intProp * intT[i]
                 + (vi*vi + ai*ai) * resProp * resT[i];
    double coefL = ei*ei * gamProp * gamL[i] + ei*vi * intProp * intL[i]
                 + (vi*vi + ai*ai) * resProp * resL[i];
    double coefA = ei*ai * intProp * intA[i] + vi*ai * resProp * resA[i];
    double sigma = coefT * (1. + pow2(cThe)) + coefL * (1. - pow2(cThe))
                 + 2. * coefA * cThe;
    sigTLA.push_back(sigma);
  }

  // Pick outgoing flavours.
  int idNew = idVec[ rndmPtr->pick(sigTLA) ];
  id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, id2, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9 && idNew < 9) setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  else if (abs(id1) < 9)         setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else if (idNew < 9)            setColAcol( 0, 0, 0, 0, 1, 0, 0, 1);
  else                           setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2ffbar2ffbarsW class.
// Cross section f_1 fbar_2 -> W+- -> f_3 fbar_4,
// i.e. W decay as part of the hard process.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2ffbarsW::initProc() {

  // Store W+- mass and width for propagator.
  mRes     = particleDataPtr->m0(24);
  GammaRes = particleDataPtr->mWidth(24);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;
  thetaWRat = 1. / (12. * couplingsPtr->sin2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(24);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2ffbarsW::sigmaKin() {

  // Set up Breit-Wigner. Common sigmaHat cross section for W+ and W-.
  double sigBW  = 12. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double preFac = alpEM * thetaWRat * mH;
  sigma0        = preFac * sigBW * particlePtr->resWidthOpen(24, mH);

  // Convert to d(sigmaHat)/d(tHat). Fixed so integrates to sigmaHat.
  sigma0       *= 3. * uH2 / (sH2 * sH);

  // Pick a decay channel.
  if (!particlePtr->preparePick(24, mH)) {
    sigma0 = 0.;
    return;
  }
  DecayChannel& channel = particlePtr->pickChannel();
  id3New = channel.product(0);
  id4New = channel.product(1);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2ffbarsW::sigmaHat() {

  // Secondary width for W+-. CKM and colour factors.
  double sigma = sigma0;
  if (abs(id1) < 9) sigma *= couplingsPtr->V2CKMid(abs(id1), abs(id2)) / 3.;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2ffbarsW::setIdColAcol() {

  // Find sign of W and its decay products. Set outgoing flavours.
  int idUp  = (abs(id1)%2 == 0) ? id1 : id2;
  id3 = (idUp > 0) ? id3New : -id3New;
  id4 = (idUp > 0) ? id4New : -id4New;
  if (id1 * id3 < 0) swap(id3, id4);
  setId( id1, id2, id3, id4);

  // Colour flow topologies. Swap when antifermions in 1 and 3.
  if (abs(id1) < 9 && abs(id3) < 9) setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  else if (abs(id1) < 9)            setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else if (abs(id3) < 9)            setColAcol( 0, 0, 0, 0, 1, 0, 0, 1);
  else                              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2ffbar2FFbarsgmZ class.
// Cross section f fbar -> gamma*/Z0 -> F Fbar.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2FFbarsgmZ::initProc() {

  // Process name.
  nameSave                 = "f fbar -> F Fbar (s-channel gamma*/Z0)";
  if (idNew == 4) nameSave = "f fbar -> c cbar (s-channel gamma*/Z0)";
  if (idNew == 5) nameSave = "f fbar -> b bbar (s-channel gamma*/Z0)";
  if (idNew == 6) nameSave = "f fbar -> t tbar (s-channel gamma*/Z0)";
  if (idNew == 7) nameSave = "f fbar -> b' b'bar (s-channel gamma*/Z0)";
  if (idNew == 8) nameSave = "f fbar -> t' t'bar (s-channel gamma*/Z0)";
  if (idNew == 15) nameSave = "f fbar -> tau+ tau- (s-channel gamma*/Z0)";
  if (idNew == 17) nameSave = "f fbar -> tau'+ tau'- (s-channel gamma*/Z0)";
  if (idNew == 18)
    nameSave   = "f fbar -> nu'_tau nu'bar_tau (s-channel gamma*/Z0)";

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode      = settingsPtr->mode("WeakZ0:gmZmode");

  // Store Z0 mass and width for propagator.
  mRes         = particleDataPtr->m0(23);
  GammaRes     = particleDataPtr->mWidth(23);
  m2Res        = mRes*mRes;
  GamMRat      = GammaRes / mRes;
  thetaWRat    = 1. / (16. * couplingsPtr->sin2thetaW()
                 * couplingsPtr->cos2thetaW());

  // Store couplings of F.
  ef           = couplingsPtr->ef(idNew);
  vf           = couplingsPtr->vf(idNew);
  af           = couplingsPtr->af(idNew);

  // Secondary open width fraction, relevant for top (or heavier).
  openFracPair = particleDataPtr->resOpenFrac(idNew, -idNew);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2FFbarsgmZ::sigmaKin() {

  // Check that above threshold.
  isPhysical     = true;
  if (mH < m3 + m4 + MASSMARGIN) {
    isPhysical   = false;
    return;
  }

  // Define average F, Fbar mass so same beta. Phase space.
  double s34Avg  = 0.5 * (s3 + s4) - 0.25 * pow2(s3 - s4) / sH;
  mr             = s34Avg / sH;
  betaf          = sqrtpos(1. - 4. * mr);

  // Final-state colour factor.
  double colF    = (idNew < 9) ? 3. * (1. + alpS / M_PI) : 1.;

  // Reconstruct decay angle so can reuse 2 -> 1 cross section.
  cosThe         = (tH - uH) / (betaf * sH);

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  gamProp = colF * M_PI * pow2(alpEM) / sH2;
  intProp = gamProp * 2. * thetaWRat * sH * (sH - m2Res)
          / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  resProp = gamProp * pow2(thetaWRat * sH)
          / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp = 0.; resProp = 0.;}
  if (gmZmode == 2) {gamProp = 0.; intProp = 0.;}

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2FFbarsgmZ::sigmaHat() {

  // Fail if below threshold.
  if (!isPhysical) return 0.;

  // Couplings for in-flavours.
  int idAbs       = abs(id1);
  double ei       = couplingsPtr->ef(idAbs);
  double vi       = couplingsPtr->vf(idAbs);
  double ai       = couplingsPtr->af(idAbs);

  // Coefficients of angular expression.
  double coefTran = ei*ei * gamProp * ef*ef + ei * vi * intProp * ef * vf
    + (vi*vi + ai*ai) * resProp * (vf*vf + pow2(betaf) * af*af);
  double coefLong = 4. * mr * ( ei*ei * gamProp * ef*ef
    + ei * vi * intProp * ef * vf + (vi*vi + ai*ai) * resProp * vf*vf );
  double coefAsym = betaf * ( ei * ai * intProp * ef * af
    + 4. * vi * ai * resProp * vf * af );

  // Combine gamma, interference and Z0 parts.
  double sigma    = coefTran * (1. + pow2(cosThe))
   + coefLong * (1. - pow2(cosThe)) + 2. * coefAsym * cosThe;

  // Top: corrections for closed decay channels.
  sigma *= openFracPair;

  // Initial-state colour factor. Answer.
  if (idAbs < 9) sigma /= 3.;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2FFbarsgmZ::setIdColAcol() {

  // Set outgoing flavours.
  id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, id2, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9 && idNew < 9) setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  else if (abs(id1) < 9)         setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else if (idNew < 9)            setColAcol( 0, 0, 0, 0, 1, 0, 0, 1);
  else                           setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles of W in top decay.

double Sigma2ffbar2FFbarsgmZ::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // For top decay hand over to standard routine, else done.
  if (idNew == 6 && process[process[iResBeg].mother1()].idAbs() == 6)
       return weightTopDecay( process, iResBeg, iResEnd);
  else return 1.;

}

//==========================================================================

// Sigma2ffbar2FfbarsW class.
// Cross section f fbar' -> W+- -> F fbar".

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2FfbarsW::initProc() {

  // Process name.
  nameSave                 = "f fbar -> F fbar (s-channel W+-)";
  if (idNew == 4) nameSave = "f fbar -> c qbar (s-channel W+-)";
  if (idNew == 5) nameSave = "f fbar -> b qbar (s-channel W+-)";
  if (idNew == 6) nameSave = "f fbar -> t qbar (s-channel W+-)";
  if (idNew == 7) nameSave = "f fbar -> b' qbar (s-channel W+-)";
  if (idNew == 8) nameSave = "f fbar -> t' qbar (s-channel W+-)";
  if (idNew == 7 && idNew2 == 6)
    nameSave = "f fbar -> b' tbar (s-channel W+-)";
  if (idNew == 8 && idNew2 == 7)
    nameSave = "f fbar -> t' b'bar (s-channel W+-)";
  if (idNew == 15 || idNew == 16)
    nameSave = "f fbar -> tau nu_taubar (s-channel W+-)";
  if (idNew == 17 || idNew == 18)
    nameSave = "f fbar -> tau'  nu'_taubar (s-channel W+-)";

  // Store W+- mass and width for propagator.
  mRes      = particleDataPtr->m0(24);
  GammaRes  = particleDataPtr->mWidth(24);
  m2Res     = mRes*mRes;
  GamMRat   = GammaRes / mRes;
  thetaWRat = 1. / (12. * couplingsPtr->sin2thetaW());

  // For t/t' want to use at least b mass.
  idPartner = idNew2;
  if ( (idNew == 6 || idNew == 8) && idNew2 == 0 ) idPartner = 5;

  // Sum of CKM weights for quarks.
  V2New     = (idNew < 9) ? couplingsPtr->V2CKMsum(idNew) : 1.;
  if (idNew2 != 0) V2New = couplingsPtr->V2CKMid(idNew, idNew2);

  // Secondary open width fractions, relevant for top or heavier.
  openFracPos = particleDataPtr->resOpenFrac( idNew, -idNew2);
  openFracNeg = particleDataPtr->resOpenFrac(-idNew,  idNew2);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2FfbarsW::sigmaKin() {

  // Check that above threshold.
  isPhysical    = true;
  if (mH < m3 + m4 + MASSMARGIN) {
    isPhysical  = false;
    return;
  }

  // Phase space factors.
  double mr1    = s3 / sH;
  double mr2    = s4 / sH;
  double betaf  = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2);

  // Reconstruct decay angle so can reuse 2 -> 1 cross section.
  double cosThe = (tH - uH) / (betaf * sH);

  // Set up Breit-Wigner and in- and out-widths.
  double sigBW  = 9. * M_PI * pow2(alpEM * thetaWRat)
                / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Initial-state colour factor.
  double colF   = (idNew < 9) ? 3. * (1. + alpS / M_PI) * V2New : 1.;

  // Angular dependence.
  double wt     = pow2(1. + betaf * cosThe) - pow2(mr1 - mr2);

  // Temporary answer.
  sigma0        = sigBW * colF * wt;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2FfbarsW::sigmaHat() {

  // Fail if below threshold.
  if (!isPhysical) return 0.;

  // CKM and colour factors.
  double sigma = sigma0;
  if (abs(id1) < 9) sigma *= couplingsPtr->V2CKMid(abs(id1), abs(id2)) / 3.;

  // Correction for secondary width in top (or heavier) decay.
  int idSame = ((abs(id1) + idNew)%2 == 0) ? id1 : id2;
  sigma *= (idSame > 0) ? openFracPos : openFracNeg;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2FfbarsW::setIdColAcol() {

  // Set outgoing flavours.
  id3 = idNew;
  id4 = (idNew2 != 0) ? idNew2 : couplingsPtr->V2CKMpick(idNew);
  if (idNew%2 == 0) {
    int idInUp = (abs(id1)%2 == 0) ? id1 : id2;
    if (idInUp > 0) id4 = -id4;
    else            id3 = -id3;
  } else {
    int idInDn = (abs(id1)%2 == 1) ? id1 : id2;
    if (idInDn > 0) id4 = -id4;
    else            id3 = -id3;
  }
  setId( id1, id2, id3, id4);

  // Swap tHat and uHat for fbar' f -> F f".
  if (id1 * id3 < 0) swapTU = true;

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9 && idNew < 9) setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  else if (abs(id1) < 9)         setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else if (idNew < 9)            setColAcol( 0, 0, 0, 0, 1, 0, 0, 1);
  else                           setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapCol12();
  if (id3 < 0) swapCol34();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles of W in top decay.

double Sigma2ffbar2FfbarsW::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // For top decay hand over to standard routine, else done.
  if (idNew == 6 && process[process[iResBeg].mother1()].idAbs() == 6)
       return weightTopDecay( process, iResBeg, iResEnd);
  else return 1.;

}

//==========================================================================

// Sigma2ffbargmZWgmZW class.
// Collects common methods for f fbar ->  gamma*/Z0/W+- gamma*/Z0/W-+.

//--------------------------------------------------------------------------

// Calculate and store internal products.

void Sigma2ffbargmZWgmZW::setupProd( Event& process, int i1, int i2,
  int i3, int i4, int i5, int i6) {

  // Store incoming and outgoing momenta,
  pRot[1] = process[i1].p();
  pRot[2] = process[i2].p();
  pRot[3] = process[i3].p();
  pRot[4] = process[i4].p();
  pRot[5] = process[i5].p();
  pRot[6] = process[i6].p();

  // Do random rotation to avoid accidental zeroes in HA expressions.
  bool smallPT = false;
  do {
    smallPT = false;
    double thetaNow = acos(2. * rndmPtr->flat() - 1.);
    double phiNow   = 2. * M_PI * rndmPtr->flat();
    for (int i = 1; i <= 6; ++i) {
      pRot[i].rot( thetaNow, phiNow);
      if (pRot[i].pT2() < 1e-4 * pRot[i].pAbs2()) smallPT = true;
    }
  } while (smallPT);

  // Calculate internal products.
  for (int i = 1; i < 6; ++i) {
    for (int j = i + 1; j <= 6; ++j) {
      hA[i][j] =
          sqrt( (pRot[i].e() - pRot[i].pz()) * (pRot[j].e() + pRot[j].pz())
        / pRot[i].pT2() ) * complex( pRot[i].px(), pRot[i].py() )
        - sqrt( (pRot[i].e() + pRot[i].pz()) * (pRot[j].e() - pRot[j].pz())
        / pRot[j].pT2() ) * complex( pRot[j].px(), pRot[j].py() );
      hC[i][j] = conj( hA[i][j] );
      if (i <= 2) {
        hA[i][j] *= complex( 0., 1.);
        hC[i][j] *= complex( 0., 1.);
      }
      hA[j][i] = - hA[i][j];
      hC[j][i] = - hC[i][j];
    }
  }

}

//--------------------------------------------------------------------------

// Evaluate the F function of Gunion and Kunszt.

complex Sigma2ffbargmZWgmZW::fGK(int j1, int j2, int j3, int j4, int j5,
  int j6) {

  return 4. * hA[j1][j3] * hC[j2][j6]
         * ( hA[j1][j5] * hC[j1][j4] + hA[j3][j5] * hC[j3][j4] );

}

//--------------------------------------------------------------------------

// Evaluate the Xi function of Gunion and Kunszt.

double Sigma2ffbargmZWgmZW::xiGK( double tHnow, double uHnow) {

  return - 4. * s3 * s4 + tHnow * (3. * tHnow + 4. * uHnow)
         + tHnow * tHnow * ( tHnow * uHnow / (s3 * s4)
           - 2. * (1. / s3 + 1./s4) * (tHnow + uHnow)
           + 2. * (s3 / s4 + s4 / s3) );

}

//--------------------------------------------------------------------------

// Evaluate the Xj function of Gunion and Kunszt.

double Sigma2ffbargmZWgmZW::xjGK( double tHnow, double uHnow) {

  return 8. * pow2(s3 + s4) - 8. * (s3 + s4) * (tHnow + uHnow)
         - 6. * tHnow * uHnow - 2. * tHnow * uHnow * ( tHnow * uHnow
           / (s3 * s4) - 2. * (1. / s3 + 1. / s4) * (tHnow + uHnow)
           + 2. * (s3 / s4 + s4 / s3) );

}

//==========================================================================

// Sigma2ffbar2gmZgmZ class.
// Cross section for f fbar -> gamma*/Z0 gamma*/Z0 (f is quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2gmZgmZ::initProc() {

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode     = settingsPtr->mode("WeakZ0:gmZmode");

  // Store Z0 mass and width for propagator.
  mRes        = particleDataPtr->m0(23);
  GammaRes    = particleDataPtr->mWidth(23);
  m2Res       = mRes*mRes;
  GamMRat     = GammaRes / mRes;
  thetaWRat   = 1. / (16. * couplingsPtr->sin2thetaW()
                * couplingsPtr->cos2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(23);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2gmZgmZ::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * pow2(alpEM) * 0.5
    * ( (tH2 + uH2 + 2. * (s3 + s4) * sH) / (tH * uH)
    - s3 * s4 * (1./tH2 + 1./uH2) );

  // Common coupling factors at the resonance masses
  double alpEM3 = couplingsPtr->alphaEM(s3);
  double alpS3  = couplingsPtr->alphaS(s3);
  double colQ3  = 3. * (1. + alpS3 / M_PI);
  double alpEM4 = couplingsPtr->alphaEM(s4);
  double alpS4  = couplingsPtr->alphaS(s4);
  double colQ4  = 3. * (1. + alpS4 / M_PI);

  // Reset quantities to sum. Declare variables in loop.
  gamSum3 = 0.;
  intSum3 = 0.;
  resSum3 = 0.;
  gamSum4 = 0.;
  intSum4 = 0.;
  resSum4 = 0.;
  int    onMode;
  double mf, mr, psvec, psaxi, betaf, ef2, efvf, vf2af2, colf;

  // Loop over all Z0 decay channels.
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    int idAbs = abs( particlePtr->channel(i).product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf     = particleDataPtr->m0(idAbs);
      onMode = particlePtr->channel(i).onMode();

      // First Z0: check that above threshold. Phase space.
      if (m3 > 2. * mf + MASSMARGIN) {
        mr    = pow2(mf / m3);
        betaf = sqrtpos(1. - 4. * mr);
        psvec = betaf * (1. + 2. * mr);
        psaxi  = pow3(betaf);

        // First Z0: combine phase space with couplings.
        ef2    = couplingsPtr->ef2(idAbs) * psvec;
        efvf   = couplingsPtr->efvf(idAbs) * psvec;
        vf2af2 = couplingsPtr->vf2(idAbs) * psvec
               + couplingsPtr->af2(idAbs) * psaxi;
        colf   = (idAbs < 6) ? colQ3 : 1.;

        // First Z0: store sum of combinations for open outstate channels.
        if (onMode == 1 || onMode == 2) {
          gamSum3 += colf * ef2;
          intSum3 += colf * efvf;
          resSum3 += colf * vf2af2;
        }
      }

      // Second Z0: check that above threshold. Phase space.
      if (m4 > 2. * mf + MASSMARGIN) {
        mr    = pow2(mf / m4);
        betaf = sqrtpos(1. - 4. * mr);
        psvec = betaf * (1. + 2. * mr);
        psaxi = pow3(betaf);

        // Second Z0: combine phase space with couplings.
        ef2    = couplingsPtr->ef2(idAbs) * psvec;
        efvf   = couplingsPtr->efvf(idAbs) * psvec;
        vf2af2 = couplingsPtr->vf2(idAbs) * psvec
               + couplingsPtr->af2(idAbs) * psaxi;
        colf   = (idAbs < 6) ? colQ4 : 1.;

        // Second Z0: store sum of combinations for open outstate channels.
        if (onMode == 1 || onMode == 2) {
          gamSum4 += colf * ef2;
          intSum4 += colf * efvf;
          resSum4 += colf * vf2af2;
        }
      }

    // End loop over fermions.
    }
  }

  // First Z0: calculate prefactors for gamma/interference/Z0 terms.
  gamProp3 = 4. * alpEM3 / (3. * M_PI * s3);
  intProp3 = gamProp3 * 2. * thetaWRat * s3 * (s3 - m2Res)
           / ( pow2(s3 - m2Res) + pow2(s3 * GamMRat) );
  resProp3 = gamProp3 * pow2(thetaWRat * s3)
           / ( pow2(s3 - m2Res) + pow2(s3 * GamMRat) );

  // First Z0: optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp3 = 0.; resProp3 = 0.;}
  if (gmZmode == 2) {gamProp3 = 0.; intProp3 = 0.;}

  // Second Z0: calculate prefactors for gamma/interference/Z0 terms.
  gamProp4 = 4. * alpEM4 / (3. * M_PI * s4);
  intProp4 = gamProp4 * 2. * thetaWRat * s4 * (s4 - m2Res)
           / ( pow2(s4 - m2Res) + pow2(s4 * GamMRat) );
  resProp4 = gamProp4 * pow2(thetaWRat * s4)
           / ( pow2(s4 - m2Res) + pow2(s4 * GamMRat) );

  // Second Z0: optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp4 = 0.; resProp4 = 0.;}
  if (gmZmode == 2) {gamProp4 = 0.; intProp4 = 0.;}

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2gmZgmZ::sigmaHat() {

  // Charge/2, left- and righthanded couplings for in-fermion.
  int idAbs = abs(id1);
  double ei = 0.5 * couplingsPtr->ef(idAbs);
  double li =       couplingsPtr->lf(idAbs);
  double ri =       couplingsPtr->rf(idAbs);

  // Combine left/right gamma, interference and Z0 parts for each Z0.
  double left3  = ei * ei * gamProp3 * gamSum3
                + ei * li * intProp3 * intSum3
                + li * li * resProp3 * resSum3;
  double right3 = ei * ei * gamProp3 * gamSum3
                + ei * ri * intProp3 * intSum3
                + ri * ri * resProp3 * resSum3;
  double left4  = ei * ei * gamProp4 * gamSum4
                + ei * li * intProp4 * intSum4
                + li * li * resProp4 * resSum4;
  double right4 = ei * ei * gamProp4 * gamSum4
                + ei * ri * intProp4 * intSum4
                + ri * ri * resProp4 * resSum4;

  // Combine left- and right-handed couplings for the two Z0's.
  double sigma = sigma0 * (left3 * left4 + right3 * right4);

  // Correct for the running-width Z0 propagators weight in PhaseSpace.
  sigma /= (runBW3 * runBW4);

  // Initial-state colour factor. Answer.
  if (idAbs < 9) sigma /= 3.;
  return  sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2gmZgmZ::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 23, 23);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate correlated decay flavours of the two gamma*/Z0.
// Unique complication, caused by gamma*/Z0 mix different left/right.

double Sigma2ffbar2gmZgmZ::weightDecayFlav( Event& process) {

  // Order so that fbar(1) f(2) -> f'(3) fbar'(4) f"(5) fbar"(6).
  i1 = (process[3].id() < 0) ? 3 : 4;
  i2 = 7 - i1;
  i3 = (process[7].id() > 0) ? 7 : 8;
  i4 = 15 - i3;
  i5 = (process[9].id() > 0) ? 9 : 10;
  i6 = 19 - i5;

  // Charge/2, left- and righthanded couplings for in- and out-fermions.
  int idAbs = process[i1].idAbs();
  double ei = 0.5 * couplingsPtr->ef(idAbs);
  double li =       couplingsPtr->lf(idAbs);
  double ri =       couplingsPtr->rf(idAbs);
  idAbs     = process[i3].idAbs();
  double e3  = 0.5 * couplingsPtr->ef(idAbs);
  double l3  =       couplingsPtr->lf(idAbs);
  double r3  =       couplingsPtr->rf(idAbs);
  idAbs      = process[i5].idAbs();
  double e4  = 0.5 * couplingsPtr->ef(idAbs);
  double l4  =       couplingsPtr->lf(idAbs);
  double r4  =       couplingsPtr->rf(idAbs);

  // Left- and righthanded couplings combined with propagators.
  c3LL = ei * ei * gamProp3 * e3 * e3
       + ei * li * intProp3 * e3 * l3
       + li * li * resProp3 * l3 * l3;
  c3LR = ei * ei * gamProp3 * e3 * e3
       + ei * li * intProp3 * e3 * r3
       + li * li * resProp3 * r3 * r3;
  c3RL = ei * ei * gamProp3 * e3 * e3
       + ei * ri * intProp3 * e3 * l3
       + ri * ri * resProp3 * l3 * l3;
  c3RR = ei * ei * gamProp3 * e3 * e3
       + ei * ri * intProp3 * e3 * r3
       + ri * ri * resProp3 * r3 * r3;
  c4LL = ei * ei * gamProp4 * e4 * e4
       + ei * li * intProp4 * e4 * l4
       + li * li * resProp4 * l4 * l4;
  c4LR = ei * ei * gamProp4 * e4 * e4
       + ei * li * intProp4 * e4 * r4
       + li * li * resProp4 * r4 * r4;
  c4RL = ei * ei * gamProp4 * e4 * e4
       + ei * ri * intProp4 * e4 * l4
       + ri * ri * resProp4 * l4 * l4;
  c4RR = ei * ei * gamProp4 * e4 * e4
       + ei * ri * intProp4 * e4 * r4
       + ri * ri * resProp4 * r4 * r4;

  // Flavour weight and maximum.
  flavWt = (c3LL + c3LR) * (c4LL + c4LR) + (c3RL + c3RR) * (c4RL + c4RR);
  double flavWtMax = (c3LL + c3LR + c3RL + c3RR) * (c4LL + c4LR + c4RL + c4RR);

  // Done.
  return flavWt / flavWtMax;

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles of the two gamma*/Z0.

double Sigma2ffbar2gmZgmZ::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Two resonance decays, but with common weight.
  if (iResBeg != 5 || iResEnd != 6) return 1.;

  // Set up four-products and internal products.
  setupProd( process, i1, i2, i3, i4, i5, i6);

  // Flip tHat and uHat if first incoming is fermion.
  double tHres   = tH;
  double uHres   = uH;
  if (process[3].id() > 0) swap( tHres, uHres);

  // Kinematics factors (norm(x) = |x|^2).
  double fGK135 = norm( fGK( 1, 2, 3, 4, 5, 6) / tHres
                      + fGK( 1, 2, 5, 6, 3, 4) / uHres );
  double fGK145 = norm( fGK( 1, 2, 4, 3, 5, 6) / tHres
                      + fGK( 1, 2, 5, 6, 4, 3) / uHres );
  double fGK136 = norm( fGK( 1, 2, 3, 4, 6, 5) / tHres
                      + fGK( 1, 2, 6, 5, 3, 4) / uHres );
  double fGK146 = norm( fGK( 1, 2, 4, 3, 6, 5) / tHres
                      + fGK( 1, 2, 6, 5, 4, 3) / uHres );
  double fGK253 = norm( fGK( 2, 1, 5, 6, 3, 4) / tHres
                      + fGK( 2, 1, 3, 4, 5, 6) / uHres );
  double fGK263 = norm( fGK( 2, 1, 6, 5, 3, 4) / tHres
                      + fGK( 2, 1, 3, 4, 6, 5) / uHres );
  double fGK254 = norm( fGK( 2, 1, 5, 6, 4, 3) / tHres
                      + fGK( 2, 1, 4, 3, 5, 6) / uHres );
  double fGK264 = norm( fGK( 2, 1, 6, 5, 4, 3) / tHres
                      + fGK( 2, 1, 4, 3, 6, 5) / uHres );

  // Weight and maximum.
  double wt     = c3LL * c4LL * fGK135 + c3LR * c4LL * fGK145
                + c3LL * c4LR * fGK136 + c3LR * c4LR * fGK146
                + c3RL * c4RL * fGK253 + c3RR * c4RL * fGK263
                + c3RL * c4RR * fGK254 + c3RR * c4RR * fGK264;
  double wtMax  = 16. * s3 * s4 * flavWt
    * ( (tHres*tHres + uHres*uHres + 2. * sH * (s3 + s4)) / (tHres * uHres)
      - s3 * s4 * (1. / (tHres*tHres) + 1. / (uHres*uHres)) );

  // Done.
  return wt / wtMax;

}

//==========================================================================

// Sigma2ffbar2ZW class.
// Cross section for f fbar' -> Z0 W+- (f is quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2ZW::initProc() {

  // Store W+- mass and width for propagator.
  mW   = particleDataPtr->m0(24);
  widW = particleDataPtr->mWidth(24);
  mWS  = mW*mW;
  mwWS = pow2(mW * widW);

  // Left-handed couplings for up/nu- and down/e-type quarks.
  lun   = (hasLeptonBeams) ? couplingsPtr->lf(12) : couplingsPtr->lf(2);
  lde   = (hasLeptonBeams) ? couplingsPtr->lf(11) : couplingsPtr->lf(1);

  // Common weak coupling factor.
  sin2thetaW = couplingsPtr->sin2thetaW();
  cos2thetaW = couplingsPtr->cos2thetaW();
  thetaWRat  = 1. / (4. * cos2thetaW);
  cotT       = sqrt(cos2thetaW / sin2thetaW);
  thetaWpt   = (9. - 8. * sin2thetaW) / 4.;
  thetaWmm   = (8. * sin2thetaW - 6.) / 4.;

  // Secondary open width fractions.
  openFracPos = particleDataPtr->resOpenFrac(23,  24);
  openFracNeg = particleDataPtr->resOpenFrac(23, -24);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2ffbar2ZW::sigmaKin() {

  // Evaluate cross section, as programmed by Merlin Kole (after tidying),
  // based on Brown, Sahdev, Mikaelian, Phys Rev. D20 (1979) 1069.
  /*
  double resBW  = 1. / (pow2(sH - mWS) + mwWS);
  double prefac = 12.0 * M_PI * pow2(alpEM) / (sH2 * 8. * sin2thetaW);
  double temp1  = tH * uH - s3 * s4;
  double temp2  = temp1 / (s3 * s4);
  double temp3  = (s3 + s4) / sH;
  double temp4  = s3 * s4 / sH;
  double partA  = temp2 * (0.25 - 0.5 * temp3
    + (pow2(s3 + s4) + 8. * s3 * s4)/(4. * sH2) )
    + (s3 + s4)/(s3 * s4) * (sH/2. - s3 - s4 + pow2(s3 - s4)/(2. * sH));
  double partB1 = lun * (0.25 * temp2 * (1. - temp3 - 4. * temp4 / tH)
    + ((s3 + s4)/(2. * s3 * s4)) * (sH - s3 - s4 + 2. * s3 * s4 / tH) );
  double partB2 = lde * ( 0.25 * temp2 * (1.- temp3 - 4. * temp4 / uH)
    + ((s3 + s4)/(2. * s3 * s4)) * (sH - s3 - s4 + 2. * s3 * s4 / uH) );
  double partE = 0.25 * temp2 + 0.5 *(s3 + s4) / temp4;
  sigma0 = prefac * (pow2(cotT) * sH2 * resBW * partA
    + 2.* sH * cotT * resBW * (sH - mWS) * (partB2 - partB1)
    + pow2(lun - lde) * partE + pow2(lde) * temp1/uH2
    + pow2(lun) * temp1/tH2 + 2. * lun * lde * sH * (s3 + s4) / (uH * tH));
  */

  // Evaluate cross section. Expression from EHLQ, with bug fix,
  // but can still give negative cross section so suspect.
  double resBW  = 1. / (pow2(sH - mWS) + mwWS);
  sigma0  = (M_PI / sH2) * 0.5 * pow2(alpEM / sin2thetaW);
  sigma0 *= sH * resBW * (thetaWpt * pT2 + thetaWmm * (s3 + s4))
    + (sH - mWS) * resBW * sH * (pT2 - s3 - s4) * (lun / tH - lde / uH)
    + thetaWRat * sH * pT2 * ( lun*lun / tH2 + lde*lde / uH2 )
    + 2. * thetaWRat * sH * (s3 + s4) * lun * lde / (tH * uH);
  // Need to protect against negative cross sections at times.
  sigma0 = max(0., sigma0);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2ZW::sigmaHat() {

  // CKM and colour factors.
  double sigma = sigma0;
  if (abs(id1) < 9) sigma *= couplingsPtr->V2CKMid(abs(id1), abs(id2)) / 3.;

  // Corrections for secondary widths in Z0 and W+- decays.
  int idUp = (abs(id1)%2 == 0) ? id1 : id2;
  sigma *= (idUp > 0) ? openFracPos : openFracNeg;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2ZW::setIdColAcol() {

  // Sign of outgoing W.
  int sign = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 23, 24 * sign);

  // tHat is defined between (f, W-) or (fbar, W+),
  // so OK for u/ubar on side 1, but must swap tHat <-> uHat if d/dbar.
  if (abs(id1)%2 == 1) swapTU = true;

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for Z0 and W+- decay angles.

double Sigma2ffbar2ZW::weightDecay( Event& process, int iResBeg, int iResEnd) {

  // Two resonance decays, but with common weight.
  if (iResBeg != 5 || iResEnd != 6) return 1.;

  // Order so that fbar(1) f(2) -> f'(3) fbar'(4) f"(5) fbar"(6)
  // with f' fbar' from W+- and f" fbar" from Z0 (note flip Z0 <-> W+-).
  int i1 = (process[3].id() < 0) ? 3 : 4;
  int i2 = 7 - i1;
  int i3 = (process[9].id() > 0) ? 9 : 10;
  int i4 = 19 - i3;
  int i5 = (process[7].id() > 0) ? 7 : 8;
  int i6 = 15 - i5;

  // Set up four-products and internal products.
  setupProd( process, i1, i2, i3, i4, i5, i6);

  // Swap tHat and uHat if incoming fermion is downtype.
  double tHres   = tH;
  double uHres   = uH;
  if (process[i2].id()%2 == 1) swap( tHres, uHres);

  //  Couplings of incoming (anti)fermions and outgoing from Z0.
  int idAbs     = process[i1].idAbs();
  double ai     = couplingsPtr->af(idAbs);
  double li1    = couplingsPtr->lf(idAbs);
  idAbs         = process[i2].idAbs();
  double li2    = couplingsPtr->lf(idAbs);
  idAbs         = process[i5].idAbs();
  double l4     = couplingsPtr->lf(idAbs);
  double r4     = couplingsPtr->rf(idAbs);

  // W propagator/interference factor.
  double Wint   = cos2thetaW * (sH - mWS) / (pow2(sH - mWS) + mwWS);

  // Combinations of couplings and kinematics (norm(x) = |x|^2).
  double aWZ    = li2 / tHres - 2. * Wint * ai;
  double bWZ    = li1 / uHres + 2. * Wint * ai;
  double fGK135 = norm( aWZ * fGK( 1, 2, 3, 4, 5, 6)
                      + bWZ * fGK( 1, 2, 5, 6, 3, 4) );
  double fGK136 = norm( aWZ * fGK( 1, 2, 3, 4, 6, 5)
                      + bWZ * fGK( 1, 2, 6, 5, 3, 4) );
  double xiT    = xiGK( tHres, uHres);
  double xiU    = xiGK( uHres, tHres);
  double xjTU   = xjGK( tHres, uHres);

  // Weight and maximum weight.
  double wt     = l4*l4 * fGK135 + r4*r4 * fGK136;
  double wtMax  = 4. * s3 * s4 * (l4*l4 + r4*r4)
                * (aWZ * aWZ * xiT + bWZ * bWZ * xiU + aWZ * bWZ * xjTU);

  // Done.
  return wt / wtMax;

}

//==========================================================================

// Sigma2ffbar2WW class.
// Cross section for f fbar -> W- W+ (f is quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2WW::initProc() {

  // Store Z0 mass and width for propagator. Common coupling factor.
  mZ           = particleDataPtr->m0(23);
  widZ         = particleDataPtr->mWidth(23);
  mZS          = mZ*mZ;
  mwZS         = pow2(mZ * widZ);
  thetaWRat    = 1. / (4. * couplingsPtr->sin2thetaW());

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(24, -24);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2ffbar2WW::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 =  (M_PI / sH2) * pow2(alpEM);

  // Z0 propagator and gamma*/Z0 interference.
  double Zprop   = sH2 / (pow2(sH - mZS) + mwZS);
  double Zint    = Zprop * (1. - mZS / sH);

  // Common coupling factors (g = gamma*, Z = Z0, f = t-channel fermion).
  cgg            = 0.5;
  cgZ            = thetaWRat * Zint;
  cZZ            = 0.5 * pow2(thetaWRat) * Zprop;
  cfg            = thetaWRat;
  cfZ            = pow2(thetaWRat) * Zint;
  cff            = pow2(thetaWRat);

  // Kinematical functions.
  double rat34   = sH * (2. * (s3 + s4) + pT2) / (s3 * s4);
  double lambdaS = pow2(sH - s3 - s4) - 4. * s3 * s4;
  double intA    = (sH - s3 - s4) * rat34 / sH;
  double intB    = 4. * (s3 + s4 - pT2);
  gSS            = (lambdaS * rat34 + 12. * sH * pT2) / sH2;
  gTT            = rat34 + 4. * sH * pT2 / tH2;
  gST            = intA + intB / tH;
  gUU            = rat34 + 4. * sH * pT2 / uH2;
  gSU            = intA + intB / uH;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2WW::sigmaHat() {

  // Flavour-specific couplings.
  int idAbs = abs(id1);
  double ei = couplingsPtr->ef(idAbs);
  double vi = couplingsPtr->vf(idAbs);
  double ai = couplingsPtr->af(idAbs);

  // Combine, with different cases for up- and down-type in-flavours.
  double sigma = sigma0;
  sigma *= (idAbs%2 == 1)
    ? (cgg * ei*ei + cgZ * ei * vi + cZZ * (vi*vi + ai*ai)) * gSS
      + (cfg * ei + cfZ * (vi + ai)) * gST + cff * gTT
    : (cgg * ei*ei + cgZ * ei * vi + cZZ * (vi*vi + ai*ai)) * gSS
      - (cfg * ei + cfZ * (vi + ai)) * gSU + cff * gUU;

  // Initial-state colour factor. Correction for secondary widths. Answer.
  if (idAbs < 9) sigma /= 3.;
  sigma *= openFracPair;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2WW::setIdColAcol() {

  // Always order W- W+, i.e. W- first.
  setId( id1, id2, -24, 24);

  // tHat is defined between (f, W-) or (fbar, W+),
  if (id1 < 0) swapTU = true;

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for W+ and W- decay angles.

double Sigma2ffbar2WW::weightDecay( Event& process, int iResBeg, int iResEnd) {

  // Two resonance decays, but with common weight.
  if (iResBeg != 5 || iResEnd != 6) return 1.;

  // Order so that fbar(1) f(2) -> f'(3) fbar'(4) f"(5) fbar"(6).
  // with f' fbar' from W- and f" fbar" from W+.
  int i1 = (process[3].id() < 0) ? 3 : 4;
  int i2 = 7 - i1;
  int i3 = (process[7].id() > 0) ? 7 : 8;
  int i4 = 15 - i3;
  int i5 = (process[9].id() > 0) ? 9 : 10;
  int i6 = 19 - i5;

  // Set up four-products and internal products.
  setupProd( process, i1, i2, i3, i4, i5, i6);

  // tHat and uHat of fbar f -> W- W+ opposite to previous convention.
  double tHres   = uH;
  double uHres   = tH;

  //  Couplings of incoming (anti)fermion.
  int idAbs     = process[i1].idAbs();
  double ai     = couplingsPtr->af(idAbs);
  double li     = couplingsPtr->lf(idAbs);
  double ri     = couplingsPtr->rf(idAbs);

  // gamma*/Z0 propagator/interference factor.
  double Zint   = mZS * (sH - mZS) / (pow2(sH - mZS) + mwZS);

  // Combinations of couplings and kinematics (norm(x) = |x|^2).
  double dWW    = (li * Zint + ai) / sH;
  double aWW    = dWW + 0.5 * (ai + 1.) / tHres;
  double bWW    = dWW + 0.5 * (ai - 1.) / uHres;
  double cWW    = ri * Zint / sH;
  double fGK135 = norm( aWW * fGK( 1, 2, 3, 4, 5, 6)
                      - bWW * fGK( 1, 2, 5, 6, 3, 4) );
  double fGK253 = norm( cWW * ( fGK( 2, 1, 5, 6, 3, 4)
                              - fGK( 2, 1, 3, 4, 5, 6) ) );
  double xiT    = xiGK( tHres, uHres);
  double xiU    = xiGK( uHres, tHres);
  double xjTU   = xjGK( tHres, uHres);

  // Weight and maximum weight.
  double wt     = fGK135 + fGK253;
  double wtMax  = 4. * s3 * s4
                * ( aWW * aWW * xiT + bWW * bWW * xiU - aWW * bWW * xjTU
                  + cWW * cWW * (xiT + xiU - xjTU) );

  // Done.
  return wt / wtMax;
}

//==========================================================================

// Sigma2ffbargmZggm class.
// Collects common methods for f fbar -> gamma*/Z0 g/gamma and permutations.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbargmZggm::initProc() {

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode     = settingsPtr->mode("WeakZ0:gmZmode");

  // Store Z0 mass and width for propagator.
  mRes        = particleDataPtr->m0(23);
  GammaRes    = particleDataPtr->mWidth(23);
  m2Res       = mRes*mRes;
  GamMRat     = GammaRes / mRes;
  thetaWRat   = 1. / (16. * couplingsPtr->sin2thetaW()
                * couplingsPtr->cos2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(23);

}

//--------------------------------------------------------------------------

// Evaluate sum of flavour couplings times phase space.

void Sigma2ffbargmZggm::flavSum() {

  // Coupling factors for Z0 subsystem.
  double alpSZ = couplingsPtr->alphaS(s3);
  double colQZ = 3. * (1. + alpSZ / M_PI);

  // Reset quantities to sum. Declare variables in loop.
  gamSum = 0.;
  intSum = 0.;
  resSum = 0.;
  int    onMode;
  double mf, mr, psvec, psaxi, betaf, ef2, efvf, vf2af2, colf;

  // Loop over all Z0 decay channels.
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    int idAbs = abs( particlePtr->channel(i).product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf = particleDataPtr->m0(idAbs);

      // Check that above threshold. Phase space.
      if (m3 > 2. * mf + MASSMARGIN) {
        mr    = pow2(mf / m3);
        betaf = sqrtpos(1. - 4. * mr);
        psvec = betaf * (1. + 2. * mr);
        psaxi = pow3(betaf);

        // Combine phase space with couplings.
        ef2    = couplingsPtr->ef2(idAbs) * psvec;
        efvf   = couplingsPtr->efvf(idAbs) * psvec;
        vf2af2 = couplingsPtr->vf2(idAbs) * psvec
               + couplingsPtr->af2(idAbs) * psaxi;
        colf   = (idAbs < 6) ? colQZ : 1.;

        // Store sum of combinations. For outstate only open channels.
        onMode = particlePtr->channel(i).onMode();
        if (onMode == 1 || onMode == 2) {
          gamSum += colf * ef2;
          intSum += colf * efvf;
          resSum += colf * vf2af2;
        }

      // End loop over fermions.
      }
    }
  }

  // Done. Return values in gamSum, intSum and resSum.

}

//--------------------------------------------------------------------------

// Calculate common parts of gamma/interference/Z0 propagator terms.

void Sigma2ffbargmZggm::propTerm() {

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  gamProp = 4. * alpEM / (3. * M_PI * s3);
  intProp = gamProp * 2. * thetaWRat * s3 * (s3 - m2Res)
          / ( pow2(s3 - m2Res) + pow2(s3 * GamMRat) );
  resProp = gamProp * pow2(thetaWRat * s3)
          / ( pow2(s3 - m2Res) + pow2(s3 * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp = 0.; resProp = 0.;}
  if (gmZmode == 2) {gamProp = 0.; intProp = 0.;}

}

//--------------------------------------------------------------------------

// Evaluate weight for gamma*/Z0 decay angle.

double Sigma2ffbargmZggm::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Z should sit in entry 5 and one more parton in entry 6.
  if (iResBeg != 5 || iResEnd != 6) return 1.;

  // In an outgoing sense fermions are labelled f(1) fbar(2) f'(3) fbar'(4)
  // where f' fbar' come from gamma*/Z0 decay.
  int i1, i2;
  int i3 = (process[7].id() > 0) ? 7 : 8;
  int i4 = 15 - i3;

  // Order so that fbar(1) f(2) -> gamma*/Z0 g/gamma.
  if (process[3].idAbs() < 20 && process[4].idAbs() < 20) {
    i1 = (process[3].id() < 0) ? 3 : 4;
    i2 = 7 - i1;

  // Order so that f(2)/fbar(1)  g/gamma -> f(1)/fbar(2) f'(3) gamma*/Z0.
  } else if (process[3].idAbs() < 20) {
    i1 = (process[3].id() < 0) ? 3 : 6;
    i2 = 9 - i1;
  } else {
    i1 = (process[4].id() < 0) ? 4 : 6;
    i2 = 10 - i1;
  }

  // Charge/2, left- and righthanded couplings for in- and out-fermion.
  int id1Abs   = process[i1].idAbs();
  double ei    = 0.5 * couplingsPtr->ef(id1Abs);
  double li    =       couplingsPtr->lf(id1Abs);
  double ri    =       couplingsPtr->rf(id1Abs);
  int id3Abs   = process[i3].idAbs();
  double ef    = 0.5 * couplingsPtr->ef(id3Abs);
  double lf    =       couplingsPtr->lf(id3Abs);
  double rf    =       couplingsPtr->rf(id3Abs);

  // Combinations of left/right for in/out, gamma*/interference/Z0.
  double clilf = ei*ei * gamProp * ef*ef + ei*li * intProp * ef*lf
               + li*li * resProp * lf*lf;
  double clirf = ei*ei * gamProp * ef*ef + ei*li * intProp * ef*rf
               + li*li * resProp * rf*rf;
  double crilf = ei*ei * gamProp * ef*ef + ei*ri * intProp * ef*lf
               + ri*ri * resProp * lf*lf;
  double crirf = ei*ei * gamProp * ef*ef + ei*ri * intProp * ef*rf
               + ri*ri * resProp * rf*rf;

  // Evaluate four-vector products.
  double p13   = process[i1].p() * process[i3].p();
  double p14   = process[i1].p() * process[i4].p();
  double p23   = process[i2].p() * process[i3].p();
  double p24   = process[i2].p() * process[i4].p();

  // Calculate weight and its maximum.
  double wt    = (clilf + crirf) * (p13*p13 + p24*p24)
               + (clirf + crilf) * (p14*p14 + p23*p23) ;
  double wtMax = (clilf + clirf + crilf + crirf)
               * (pow2(p13 + p14) + pow2(p23 + p24));

  // Done.
  return (wt / wtMax);

}

//==========================================================================

// Sigma2qqbar2gmZg class.
// Cross section for q qbar -> gamma*/Z0 g.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2qqbar2gmZg::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * (alpEM * alpS)
    * (2./9.) * (tH2 + uH2 + 2. * sH * s3) / (tH * uH);

  // Calculate flavour sums for final state.
  flavSum();

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  propTerm();

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2qqbar2gmZg::sigmaHat() {

  // Combine gamma, interference and Z0 parts.
  int idAbs    = abs(id1);
  double sigma = sigma0
               * ( couplingsPtr->ef2(idAbs)    * gamProp * gamSum
                 + couplingsPtr->efvf(idAbs)   * intProp * intSum
                 + couplingsPtr->vf2af2(idAbs) * resProp * resSum);

  // Correct for the running-width Z0 propagater weight in PhaseSpace.
  sigma       /= runBW3;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2gmZg::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 23, 21);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qg2gmZq class.
// Cross section for q g -> gamma*/Z0 q.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2qg2gmZq::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * (alpEM * alpS)
    * (1./12.) * (sH2 + uH2 + 2. * tH * s3) / (-sH * uH);

  // Calculate flavour sums for final state.
  flavSum();

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  propTerm();

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2qg2gmZq::sigmaHat() {

  // Combine gamma, interference and Z0 parts.
  int idAbs    = (id2 == 21) ? abs(id1) : abs(id2);
  double sigma = sigma0
               * ( couplingsPtr->ef2(idAbs)    * gamProp * gamSum
                 + couplingsPtr->efvf(idAbs)   * intProp * intSum
                 + couplingsPtr->vf2af2(idAbs) * resProp * resSum);

  // Correct for the running-width Z0 propagater weight in PhaseSpace.
  sigma       /= runBW3;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2gmZq::setIdColAcol() {

  // Flavour set up for q g -> gamma*/Z0 q.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, 23, idq);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21);

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2ffbar2gmZgm class.
// Cross section for f fbar -> gamma*/Z0 gamma.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2gmZgm::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * (alpEM*alpEM)
    * 0.5 * (tH2 + uH2 + 2. * sH * s3) / (tH * uH);

  // Calculate flavour sums for final state.
  flavSum();

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  propTerm();


}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2gmZgm::sigmaHat() {

  // Combine gamma, interference and Z0 parts.
  int idAbs    = abs(id1);
  double sigma = sigma0 * couplingsPtr->ef2(idAbs)
               * ( couplingsPtr->ef2(idAbs)    * gamProp * gamSum
                 + couplingsPtr->efvf(idAbs)   * intProp * intSum
                 + couplingsPtr->vf2af2(idAbs) * resProp * resSum);

  // Correct for the running-width Z0 propagater weight in PhaseSpace.
  sigma       /= runBW3;

  // Colour factor. Answer.
  if (idAbs < 9) sigma /= 3.;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2gmZgm::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 23, 22);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2fgm2gmZf class.
// Cross section for f gamma -> gamma*/Z0 f'.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2fgm2gmZf::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * (alpEM*alpEM)
    * 0.5 * (sH2 + uH2 + 2. * tH * s3) / (- sH * uH);

  // Calculate flavour sums for final state.
  flavSum();

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  propTerm();

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2fgm2gmZf::sigmaHat() {

  // Combine gamma, interference and Z0 parts.
  int idAbs    = (id2 == 22) ? abs(id1) : abs(id2);
  double sigma = sigma0 * couplingsPtr->ef2(idAbs)
               * ( couplingsPtr->ef2(idAbs)    * gamProp * gamSum
                 + couplingsPtr->efvf(idAbs)   * intProp * intSum
                 + couplingsPtr->vf2af2(idAbs) * resProp * resSum);

  // Correct for the running-width Z0 propagater weight in PhaseSpace.
  sigma         /= runBW3;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2fgm2gmZf::setIdColAcol() {

  // Flavour set up for q gamma -> gamma*/Z0 q.
  int idq = (id2 == 22) ? id1 : id2;
  setId( id1, id2, 23, idq);

  // tH defined between f and f': must swap tHat <-> uHat if q gamma in.
  swapTU = (id2 == 22);

  // Colour flow topologies. Swap when antiquarks.
  if      (abs(id1) < 9) setColAcol( 1, 0, 0, 0, 0, 0, 1, 0);
  else if (abs(id2) < 9) setColAcol( 0, 0, 1, 0, 0, 0, 1, 0);
  else                   setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2ffbarWggm class.
// Collects common methods for f fbar -> W+- g/gamma and permutations.

//--------------------------------------------------------------------------

// Evaluate weight for W+- decay angle.

double Sigma2ffbarWggm::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // W should sit in entry 5 and one more parton in entry 6.
  if (iResBeg != 5 || iResEnd != 6) return 1.;

  // In an outgoing sense fermions are labelled f(1) fbar(2) f'(3) fbar'(4)
  // where f' fbar' come from W+- decay.
  int i1, i2;
  int i3 = (process[7].id() > 0) ? 7 : 8;
  int i4 = 15 - i3;

  // Order so that fbar(1) f(2) -> W+- g/gamma.
  if (process[3].idAbs() < 20 && process[4].idAbs() < 20) {
    i1 = (process[3].id() < 0) ? 3 : 4;
    i2 = 7 - i1;

  // Order so that f(2)/fbar(1)  g/gamma -> f(1)/fbar(2) f'(3) W+-.
  } else if (process[3].idAbs() < 20) {
    i1 = (process[3].id() < 0) ? 3 : 6;
    i2 = 9 - i1;
  } else {
    i1 = (process[4].id() < 0) ? 4 : 6;
    i2 = 10 - i1;
  }

  // Evaluate four-vector products.
  double p13 = process[i1].p() * process[i3].p();
  double p14 = process[i1].p() * process[i4].p();
  double p23 = process[i2].p() * process[i3].p();
  double p24 = process[i2].p() * process[i4].p();

  // Calculate weight and its maximum.
  double wt    = pow2(p13) + pow2(p24);
  double wtMax = pow2(p13 + p14) + pow2(p23 + p24);

  // Done.
  return (wt / wtMax);

}

//==========================================================================

// Sigma2qqbar2Wg class.
// Cross section for q qbar' -> W+- g.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qqbar2Wg::initProc() {

  // Secondary open width fractions, relevant for top (or heavier).
  openFracPos = particleDataPtr->resOpenFrac(24);
  openFracNeg = particleDataPtr->resOpenFrac(-24);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2qqbar2Wg::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * (alpEM * alpS / couplingsPtr->sin2thetaW())
    * (2./9.) * (tH2 + uH2 + 2. * sH * s3) / (tH * uH);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2qqbar2Wg::sigmaHat() {

  // CKM factor. Secondary width for W+ or W-.
  double sigma = sigma0 * couplingsPtr->V2CKMid(abs(id1), abs(id2));
  int idUp     = (abs(id1)%2 == 0) ? id1 : id2;
  sigma       *= (idUp > 0) ? openFracPos : openFracNeg;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2Wg::setIdColAcol() {

  // Sign of outgoing W.
  int sign = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 24 * sign, 21);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qg2Wq class.
// Cross section for q g -> W+- q'.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qg2Wq::initProc() {

  // Secondary open width fractions, relevant for top (or heavier).
  openFracPos = particleDataPtr->resOpenFrac(24);
  openFracNeg = particleDataPtr->resOpenFrac(-24);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2qg2Wq::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * (alpEM * alpS / couplingsPtr->sin2thetaW())
    * (1./12.) * (sH2 + uH2 + 2. * tH * s3) / (-sH * uH);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2qg2Wq::sigmaHat() {

  // CKM factor. Secondary width for W+ or W-.
  int idAbs    = (id2 == 21) ? abs(id1) : abs(id2);
  double sigma = sigma0 * couplingsPtr->V2CKMsum(idAbs);
  int idUp     = (id2 == 21) ? id1 : id2;
  if (idAbs%2 == 1) idUp = -idUp;
  sigma       *= (idUp > 0) ? openFracPos : openFracNeg;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2Wq::setIdColAcol() {

  // Sign of outgoing W. Flavour of outgoing quark.
  int idq           = (id2 == 21) ? id1 : id2;
  int sign          = 1 - 2 * (abs(idq)%2);
  if (idq < 0) sign = -sign;
  id4 = couplingsPtr->V2CKMpick(idq);

  // Flavour set up for q g -> W q.
  setId( id1, id2, 24 * sign, id4);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21);

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2ffbar2Wgm class.
// Cross section for f fbar' -> W+- gamma.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2Wgm::initProc() {

  // Secondary open width fractions, relevant for top (or heavier).
  openFracPos = particleDataPtr->resOpenFrac(24);
  openFracNeg = particleDataPtr->resOpenFrac(-24);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2Wgm::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * (alpEM*alpEM / couplingsPtr->sin2thetaW())
    * 0.5 * (tH2 + uH2 + 2. * sH * s3) / (tH * uH);
}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2Wgm::sigmaHat() {

  // Extrafactor different for e nu and q qbar' instate.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  double chgUp = (id1Abs > 10) ? 0. : 2./3.;
  double sigma = sigma0 * pow2( chgUp - tH / (tH + uH) );

  // CKM and colour factors. Secondary width for W+ or W-.
  if (id1Abs < 9) sigma *= couplingsPtr->V2CKMid(id1Abs, id2Abs) / 3.;
  int idUp     = (abs(id1)%2 == 0) ? id1 : id2;
  sigma       *= (idUp > 0) ? openFracPos : openFracNeg;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2Wgm::setIdColAcol() {

  // Sign of outgoing W.
  int sign          = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 24 * sign, 22);

  // tH defined between (f,W-) or (fbar',W+).
  swapTU = (sign * id1 > 0);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2fgm2Wf class.
// Cross section for f gamma -> W+- f'.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2fgm2Wf::initProc() {

  // Secondary open width fractions, relevant for top (or heavier).
  openFracPos = particleDataPtr->resOpenFrac(24);
  openFracNeg = particleDataPtr->resOpenFrac(-24);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2fgm2Wf::sigmaKin() {

  // Cross section part common for all incoming flavours.
  sigma0 = (M_PI / sH2) * (alpEM*alpEM / couplingsPtr->sin2thetaW())
    * 0.5 * (sH2 + uH2 + 2. * tH * s3) / (pT2 * s3 - sH * uH);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2fgm2Wf::sigmaHat() {

  // Extrafactor dependent on charge of incoming fermion.
  int idAbs     = (id2 == 22) ? abs(id1) : abs(id2);
  double charge = (idAbs > 10) ? 1. : ( (idAbs%2 == 1) ? 1./3. : 2./3. );
  double sigma  = sigma0 * pow2( charge  - sH / (sH + uH) );

  // CKM factor. Secondary width for W+ or W-.
  sigma        *= couplingsPtr->V2CKMsum(idAbs);
  int idUp      = (id2 == 22) ? id1 : id2;
  if (idAbs%2 == 1) idUp = -idUp;
  sigma        *= (idUp > 0) ? openFracPos : openFracNeg;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2fgm2Wf::setIdColAcol() {

  // Sign of outgoing W. Flavour of outgoing fermion.
  int idq           = (id2 == 22) ? id1 : id2;
  int sign          = 1 - 2 * (abs(idq)%2);
  if (idq < 0) sign = -sign;
  id4 = couplingsPtr->V2CKMpick(idq);

  // Flavour set up for q gamma -> W q.
  setId( id1, id2, 24 * sign, id4);

  // tH defined between f and f': must swap tHat <-> uHat if q gamma in.
  swapTU = (id2 == 22);

  // Colour flow topologies. Swap when antiquarks.
  if      (abs(id1) < 9) setColAcol( 1, 0, 0, 0, 0, 0, 1, 0);
  else if (abs(id2) < 9) setColAcol( 0, 0, 1, 0, 0, 0, 1, 0);
  else                   setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2gmgm2ffbar class.
// Cross section for gamma gamma -> l lbar.

//--------------------------------------------------------------------------

// Initialize process wrt flavour.

void Sigma2gmgm2ffbar::initProc() {

  // Process name.
  nameSave = "gamma gamma -> f fbar";
  if (idNew ==  1) nameSave = "gamma gamma -> q qbar (uds)";
  if (idNew ==  4) nameSave = "gamma gamma -> c cbar";
  if (idNew ==  5) nameSave = "gamma gamma -> b bbar";
  if (idNew ==  6) nameSave = "gamma gamma -> t tbar";
  if (idNew == 11) nameSave = "gamma gamma -> e+ e-";
  if (idNew == 13) nameSave = "gamma gamma -> mu+ mu-";
  if (idNew == 15) nameSave = "gamma gamma -> tau+ tau-";

  // Generate massive phase space, except for u+d+s.
  idMass = 0;
  if (idNew > 3) idMass = idNew;

  // Charge and colour factor.
  ef4 = 1.;
  if (idNew == 1) ef4 = 3. * (pow4(2./3.) + 2. * pow4(1./3.));
  if (idNew == 4 || idNew == 6) ef4 = 3. * pow4(2./3.);
  if (idNew == 5) ef4 = 3. * pow4(1./3.);

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(idNew, -idNew);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence.

void Sigma2gmgm2ffbar::sigmaKin() {

  // Pick current flavour for u+d+s mix by e_q^4 weights.
  if (idNew == 1) {
    double rId = 18. * rndmPtr->flat();
    idNow = 1;
    if (rId > 1.)  idNow = 2;
    if (rId > 17.) idNow = 3;
    s34Avg = pow2(particleDataPtr->m0(idNow));
  } else {
    idNow = idNew;
    s34Avg = 0.5 * (s3 + s4) - 0.25 * pow2(s3 - s4) / sH;
  }

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  double tHQ    = -0.5 * (sH - tH + uH);
  double uHQ    = -0.5 * (sH + tH - uH);
  double tHQ2   = tHQ * tHQ;
  double uHQ2   = uHQ * uHQ;

  // Calculate kinematics dependence.
  if (sH < 4. * s34Avg) sigTU = 0.;
  else sigTU = 2. * (tHQ * uHQ - s34Avg * sH)
    * (tHQ2 + uHQ2 + 2. * s34Avg * sH) / (tHQ2 * uHQ2);

  // Answer.
  sigma = (M_PI / sH2) * pow2(alpEM) * ef4 * sigTU * openFracPair;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gmgm2ffbar::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idNow, -idNow);

  // Colour flow in singlet state.
  if (idNow < 10) setColAcol( 0, 0, 0, 0, 1, 0, 0, 1);
  else            setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//==========================================================================

} // end namespace Pythia8
