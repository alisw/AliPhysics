// SigmaDM.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// Dark Matter simulation classes.

#include "Pythia8/SigmaDM.h"

namespace Pythia8 {

//==========================================================================

// Sigma2ffbar2Zp2XX.
// Cross section for f fbar' -> Zprime -> XX. (Zprime a.k.a. DMmed(s=1).)

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ffbar2Zp2XX::initProc() {

  // Store mass and width for propagator, and couplings.
  kinMix    = settingsPtr->parm("Zp:kineticMixing");
  mRes      = particleDataPtr->m0(55);
  GammaRes  = particleDataPtr->mWidth(55);
  m2Res     = mRes*mRes;
  alpEM     = couplingsPtr->alphaEM(m2Res);
  gZp       = settingsPtr->parm("Zp:gZp");
  eps       = settingsPtr->parm("Zp:epsilon");

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(55);

  // Loop over all Zp decay channels.
  double mf, mr, psvec, psaxi, betaf;
  int decMode = settingsPtr->mode("Zp:decayMode");
  preFac = 0.0;
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    int idAbs = abs( particlePtr->channel(i).product(0));
    double vf = 0.;
    double af = 0.;

    // Turn off decay modes depending on settings.
    bool turnOff = false;
    // DM-only.
    if (idAbs != 52 && decMode == 0) turnOff = true;
    // Dijet only.
    else if (decMode == 1 && idAbs > 10) turnOff = true;
    else if (decMode > 1) {
      if (idAbs < 10 || idAbs > 20) turnOff = true;
      // Dilepton.
      if (decMode == 2 && idAbs %2 == 0 ) turnOff = true;
      // Invisible to neutrinos.
      if (decMode == 3 && idAbs %2 ) turnOff = true;
    }
    if (turnOff) {
      particlePtr->channel(i).onMode(0);
      continue;
    }

    // Set couplings: quarks.
    if (idAbs < 7) {
      if (abs(id1)%2 == 0) {
        if (kinMix) {
          vf = eps * (2./3. + couplingsPtr->vf(2));
          af = eps * couplingsPtr->af(2);
        } else {
          vf = settingsPtr->parm("Zp:vu");
          af = settingsPtr->parm("Zp:au");
        }
      } else {
        if (kinMix) {
          vf = eps * (-1./3. + couplingsPtr->vf(1));
          af = eps * couplingsPtr->af(1);
        } else {
          vf = settingsPtr->parm("Zp:vd");
          af = settingsPtr->parm("Zp:ad");
        }
      }
    }

    // Set couplings: leptons.
    if (idAbs > 10 && idAbs < 17) {
      if (abs(id1)%2 == 0) {
        if (kinMix) {
          vf = eps * couplingsPtr->vf(12);
          af = eps * couplingsPtr->af(12);
        } else {
          vf = settingsPtr->parm("Zp:vv");
          af = settingsPtr->parm("Zp:av");
        }
      } else {
        if (kinMix) {
          vf = eps * (-1. + couplingsPtr->vf(11));
          af = eps * couplingsPtr->af(11);
        } else {
          vf = settingsPtr->parm("Zp:vl");
          af = settingsPtr->parm("Zp:al");
        }
      }
    }

    // Set couplings: DM.
    if (idAbs == 52) {
      vf = settingsPtr->parm("Zp:vX");
      af = settingsPtr->parm("Zp:aX");
    }

    // Check that mass is above threshold. Calculate phase space.
    mf = particleDataPtr->m0(idAbs);
    if (mRes > 2. * mf + MASSMARGIN) {
      mr    = pow2(mf / mRes);
      betaf = sqrtpos(1. - 4. * mr);
      psvec = betaf * (1. + 2. * mr);
      psaxi = pow3(betaf);
      // For kinetic mixing, coupling to SM fermions is through EM coupling.
      double fac = (kinMix && idAbs != 52) ? 4.0 * M_PI * alpEM : pow2(gZp);
      if (idAbs < 10) fac *= 3.0;
      preFac += fac * (vf * vf * psvec + af * af * psaxi ) ;

    }
  }

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1ffbar2Zp2XX::sigmaKin() {

  sigma0 = mRes / ( pow2(sH - m2Res) + pow2(mRes * GammaRes) );

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma1ffbar2Zp2XX::sigmaHat() {

  // Check for allowed flavour combinations.
  if (id1 + id2 != 0 || abs(id1) > 6 ) return 0.;

  // Set couplings of incoming quarks.
  double vf, af;
  if (abs(id1)%2 == 0) {
    if (kinMix) {
      vf = eps * couplingsPtr->vf(2);
      af = eps * couplingsPtr->af(2);
    } else {
      vf = settingsPtr->parm("Zp:vu");
      af = settingsPtr->parm("Zp:au");
    }
  } else {
    if (kinMix) {
      vf = eps * couplingsPtr->vf(1);
      af = eps * couplingsPtr->af(1);
    } else {
      vf = settingsPtr->parm("Zp:vd");
      af = settingsPtr->parm("Zp:ad");
    }
  }

  // Combine for cross section.
  double coup = pow2(gZp);
  if (kinMix) coup = 4.0 * M_PI * alpEM;
  double vf2af2 = coup * (vf * vf + af * af);
  double sigma = preFac * sigma0 * vf2af2;

  // Colour factor.
  if (abs(id1) < 7) sigma /= 3;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2Zp2XX::setIdColAcol() {

  setId(id1, id2, 55);
  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2Zpg2XXj.
// Cross section for q qbar -> Zprime -> XX + jet. (Zprime a.k.a. DMmed(s=1).)

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qqbar2Zpg2XXj::initProc() {

  // Store mass and width for propagator, and couplings.
  kinMix    = settingsPtr->flag("Zp:kineticMixing");
  mRes      = particleDataPtr->m0(55);
  GammaRes  = particleDataPtr->mWidth(55);
  m2Res     = mRes*mRes;
  alpEM     = couplingsPtr->alphaEM(m2Res);
  gZp       = settingsPtr->parm("Zp:gZp");
  eps       = settingsPtr->parm("Zp:epsilon");

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(55);

  // Turn off all decay modes except DM (dark photon).
  preFac = 0.0;
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    int idAbs = abs( particlePtr->channel(i).product(0));
    if (idAbs < 20) particlePtr->channel(i).onMode(0);
  }
  preFac = particleDataPtr->resOpenFrac(52, -52);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2qqbar2Zpg2XXj::sigmaKin() {

  double propZp = s3 / ( pow2(s3 - m2Res) + pow2(mRes * GammaRes) );
  double alpD = pow2(gZp) / 4.0 / M_PI;
  if (kinMix) alpD = alpEM;
  sigma0 = (M_PI / sH2) * (alpD * alpS) * propZp
    * (2./9.) * (tH2 + uH2 + 2. * sH * s3) / (tH * uH);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma2qqbar2Zpg2XXj::sigmaHat() {

  // Check for allowed flavour combinations
  if (id1 + id2 != 0 || abs(id1) > 6 ) return 0.;

  // Set couplings of incoming quarks.
  double vf, af;
  if (abs(id1)%2 == 0) {
    if (kinMix) {
      vf = eps * couplingsPtr->vf(2);
      af = eps * couplingsPtr->af(2);
    } else {
      vf = settingsPtr->parm("Zp:vu");
      af = settingsPtr->parm("Zp:au");
    }
  } else {
    if (kinMix) {
      vf = eps * couplingsPtr->vf(1);
      af = eps * couplingsPtr->af(1);
    } else {
      vf = settingsPtr->parm("Zp:vd");
      af = settingsPtr->parm("Zp:ad");
    }
  }

  // Combine for cross section.
  double vf2af2 = vf * vf + af * af;
  double sigma = sigma0 * vf2af2 * preFac;

    // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2Zpg2XXj::setIdColAcol() {

  setId(id1, id2, 55, 21);

  // Colour flow topologies.
  if (id1 > 0) setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  else         setColAcol( 0, 2, 1, 0, 0, 0, 1, 2);

}

//==========================================================================

// Sigma2ffbar2Zp2H.
// Cross section for f fbar' -> Zprime H.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2ZpH::initProc() {

  // Store mass and width for propagator, and couplings.
  kinMix    = settingsPtr->flag("Zp:kineticMixing");
  mRes      = particleDataPtr->m0(55);
  GammaRes  = particleDataPtr->mWidth(55);
  m2Res     = mRes*mRes;
  coupZpH   = settingsPtr->parm("Zp:coupH");
  gZp       = settingsPtr->parm("Zp:gZp");
  eps       = settingsPtr->parm("Zp:epsilon");
  if (kinMix) coupZpH = eps;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(55);

  // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(55, 25);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2ffbar2ZpH::sigmaKin() {

  double denom = (pow2(sH - m2Res) + pow2(mRes * GammaRes));
  sigma0 = (M_PI / sH2) * 8. * pow2(gZp * coupZpH)
    * (tH * uH - s3 * s4 + 2. * sH * s4);
  sigma0 /= denom;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma2ffbar2ZpH::sigmaHat() {

  // Check for allowed flavour combinations.
  if (id1 + id2 != 0 ) return 0.;

  // Coupling a_f^2 + v_f^2 to s-channel Zp and colour factor.
  double vf = 0., af = 0.;
  if (abs(id1)%2 == 0) {
    if (kinMix) {
      vf = eps * couplingsPtr->vf(2);
      af = eps * couplingsPtr->af(2);
    } else {
      vf = settingsPtr->parm("Zp:vu");
      af = settingsPtr->parm("Zp:au");
    }
  } else {
    if (kinMix) {
      vf = eps * couplingsPtr->vf(1);
      af = eps * couplingsPtr->af(1);
    } else {
      vf = settingsPtr->parm("Zp:vd");
      af = settingsPtr->parm("Zp:ad");
    }
  }

  // Combine for cross section, inbcluding colour factor.
  double sigma = sigma0 * (vf * vf + af * af);
  if (abs(id1) < 9) sigma /= 3.;

  // Secondary width for Zp and H.
  sigma       *= openFrac;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2ZpH::setIdColAcol() {

  setId(id1, id2, 55, 25);
  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2ffbar2S2XX.
// Cross section for f fbar' -> S -> XX.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1gg2S2XX::initProc() {

  // Store mass and width for propagator.
  mRes      = particleDataPtr->m0(54);
  GammaRes  = particleDataPtr->mWidth(54);
  m2Res     = mRes*mRes;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(54);

  // Turn off all decay modes except DM.
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    int idAbs = abs( particlePtr->channel(i).product(0));
    if (idAbs != 52) particlePtr->channel(i).onMode(0);
  }

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1gg2S2XX::sigmaKin() {

  double propS = sH / ( pow2(sH - m2Res) + pow2(mRes * GammaRes) );
  sigma0        = 8. * M_PI * propS;

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma1gg2S2XX::sigmaHat() {

  // Check for allowed flavour combinations
  if (id1 != id2 || abs(id1) != 21 ) return 0.;

  // Incoming gg -> S width, colour-corrected.
  double widthIn  = particlePtr->resWidthChan( mRes, 21, 21) / 64.;

  // Width out only includes open channels.
  double widthOut = particlePtr->resWidthChan( mRes, 52, -52);

  // Done.
  double sigma = widthIn * sigma0 * widthOut;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1gg2S2XX::setIdColAcol() {

  setId(id1, id2, 54);
  setColAcol( 1, 2, 2, 1, 0, 0);

}

//==========================================================================

// Sigma2gg2Sg2XXj.
// Cross section for g g -> S g -> XX + jet.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2Sg2XXj::initProc() {

  // Store mass and width for propagator.
  mRes      = particleDataPtr->m0(54);
  GammaRes  = particleDataPtr->mWidth(54);
  m2Res     = mRes*mRes;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(54);

  // Turn off all decay modes except DM.
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    int idAbs = abs( particlePtr->channel(i).product(0));
    if (idAbs != 52) particlePtr->channel(i).onMode(0);
  }

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2gg2Sg2XXj::sigmaKin() {

  double wid = particlePtr->resWidthChan(m3, 21, 21);
  sigma0  = (M_PI / sH2) * (3. / 16.) * alpS * (wid / m3)
    * (sH2 * sH2 + tH2 * tH2 + uH2 * uH2 + pow2(sH2))
    / (sH * tH * uH * sH);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma2gg2Sg2XXj::sigmaHat() {

  return sigma0 * particlePtr->resWidthChan(mRes, 52, -52);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2Sg2XXj::setIdColAcol() {

  setId(id1, id2, 54, 21);

  if( rndmPtr->flat() < 0.5)
    setColAcol( 1, 2, 3, 1, 0, 0, 3, 2);
  else
    setColAcol( 1, 2, 2, 3, 0, 0, 1, 3);

}

//==========================================================================

void Sigma2qqbar2DY::initProc() {

  // Type of model, notably outgoing DM particle pair.
  type  = settingsPtr->mode("DM:DYtype");
  nplet = settingsPtr->mode("DM:Nplet");
  if (type == 1) {
    nameSave = "q qbar -> Sl(DM) Sl(DM)*";
    id3 = 56;
    id4 = -56;
  } else if (type == 2) {
    nameSave = "q qbar -> X+ X-";
    id3 = 57;
    id4 = -57;
  } else if (type == 3) {
    nameSave = "q qbar -> X++ X--";
    id3 = 59;
    id4 = -59;
  } else if (type == 4) {
    nameSave = "q qbar' -> X2 X+ + c.c.";
    id3 = 57;
    id4 = 58;
    isUD = true;
  }

  // Set particle masses based on couplings (next-to-minimal DM).
  M1     = settingsPtr->parm("DM:M1");
  M2     = settingsPtr->parm("DM:M2");
  Lambda = settingsPtr->parm("DM:Lambda");

  // Mixing parameters.
  double vev = 174.0;
  double mixing = vev / Lambda;
  if (type > 1) mixing *= sqrt(2) * vev;
  if (type > 2) mixing *= pow2(vev) / pow2(Lambda) / sqrt(12);
  double term1 = sqrt(pow2(M2 - M1) + pow2(mixing));
  double sin2th = 0.5 * (1 - abs(M2 - M1) / term1);

  // Switch to find n-plet neutral partner.
  if (type >= 2) {

    // Z couplings for singly and doubly charged partners.
    coupW11 = sqrt(sin2th); // chi+ chi1
    coupW12 = sqrt(1.0 - sin2th); // chi+ chi2
    coupW2  = 1.0; // chi++ chi+
    if (nplet == 3) {
      coupW11 *= sqrt(3);
      coupW12 *= sqrt(3);
      coupW2  *= sqrt(3);
    }
    if (type == 4 && coupW12 < coupW11) {
      id4 = 52;
    }
  }

  // Set propagator mass.
  if (!isUD) {
    mRes      = particleDataPtr->m0(23);
    GammaRes  = particleDataPtr->mWidth(23);
    m2Res     = mRes*mRes;
  } else {
    mRes      = particleDataPtr->m0(24);
    GammaRes  = particleDataPtr->mWidth(24);
    m2Res     = mRes*mRes;
  }
  xW = couplingsPtr->sin2thetaW();

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(id3, id4);

}
//--------------------------------------------------------------------------

void Sigma2qqbar2DY::sigmaKin() {

  // Z/W propagator.
  double sV = sH - m2Res;
  double d  = pow2(sV) + pow2(mRes * GammaRes);
  propRes   = complex( sV / d, mRes * GammaRes / d);

  sigma0 =  M_PI/(4 * sH2) * openFracPair * pow2(alpEM);

}

//--------------------------------------------------------------------------

double Sigma2qqbar2DY::sigmaHat() {

  // In-pair must be opposite-sign.
  if (id1 * id2 > 0) return 0.0;

  // Common factor for LR and RL contributions.
  double facTU = uH*tH-s3*s4;
  int    id1A  = abs(id1);
  double eQ    = (id1A%2 == 0) ? 2./3. : -1./3. ;
  double eSl   = -1. ;
  double LqZ   = couplingsPtr->lf(abs(id1));
  double RqZ   = couplingsPtr->rf(abs(id1));

  double sumColS = 0., sumColT = 0., sumInterference = 0.;
  double LL = 0., RR = 0.;

  // Couplings for singly charged partners
  if (nplet == 1) { LL = 1.0 - 2.0 * xW; RR = -2.0 * xW; }
  if (nplet == 2 || nplet == 3) { LL = 2.0 - 2.0 * xW; RR = -2.0 * xW; }

  // Coupling for doubly charged partner.
  if (type == 3) { LL = 4.0 - 2.0 * xW; RR = -2.0 * xW; }

  // s-channel Z/photon and interference.
  if (abs(id1) == abs(id2) && abs(id3) == abs(id4)) {
    double CoupZ;
    CoupZ =  couplingsPtr->rf(11);

    // Scalar lepton production.
    if (type == 1) {
      // s-channel Z.
      sumColS += sigma0 * facTU / 16.0 / pow2(xW) / pow2(1.0-xW)
        * norm(propRes) * CoupZ  * ( pow2(LqZ) + pow2(RqZ) );
      // gamma: factor 2 since contributes to both ha != hb helicities.
      sumColS += (abs(CoupZ) > 0.0) ? 2. * pow2(eQ) * pow2(eSl) * sigma0
        * facTU / pow2(sH) : 0.0;
      // Z/gamma interference.
      sumInterference += eQ * eSl * sigma0 * facTU / 2.0 / xW / (1.-xW)
        * sqrt(norm(propRes)) / sH * CoupZ * (LqZ + RqZ);
    }

    // Production of (singly and doubly) charged partners.
    if (type > 1 && type < 4) {
      double fac = (uH - s3) * (uH - s4) + (tH - s3) * (tH - s4)
        + 2 * m3 * m4 * sH;
      sumColS += sigma0 * fac * norm(propRes) * (pow2(LL) + pow2(RR))
        * ( pow2(LqZ) + pow2(RqZ) );
      sumColS += (abs(CoupZ) > 0.0) ? 2. * pow2(eQ) * pow2(eSl) * sigma0 * fac
        / pow2(sH) : 0.0;
      sumInterference += eQ * eSl * sigma0 * fac / 2.0 / xW / (1.-xW)
        * sqrt(norm(propRes)) / sH * CoupZ * (LqZ + RqZ);
    }

  // Production of neutral and singly-charged partners.
  } else if (type == 4) {
    if (!isUD || id1 * id2 > 0) return 0.0;
    if (abs(id1)%2 + abs(id2)%2 != 1) return 0.0 ;

    // s-channel W contribution
    double coupW = coupW11 > coupW12 ? coupW11 : coupW12;
    double fac = pow2(coupW) * norm(propRes) / 2.0;
    sumColS = (uH - s3) * (uH - s4) + (tH - s3) * (tH - s4) + 2 * m3 * m4 * sH;
    sumColS *= sigma0 * fac / xW ;
  }

  // t-channel for lepton colliders only (TODO)

  // Cross section.
  double sigma = sumColS + sumColT + sumInterference;

  return sigma;

}

//--------------------------------------------------------------------------

void Sigma2qqbar2DY::setIdColAcol() {

  // Normal or charge conjugate process.
  int up = abs(id1)%2 == 0 ? id1 : id2;
  if (up < 0 && abs(id3) == 57 && id4 == 58)
    setId(id1, id2, -57, 58);
  else
    setId(id1, id2, id3, id4);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

} // end namespace Pythia8
