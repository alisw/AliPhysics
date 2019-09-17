// SigmaExtraDim.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Stefan Ask for the *LED* routines.
// Function definitions (not found in the header) for the
// extra-dimensional simulation classes.

#include "Pythia8/SigmaExtraDim.h"

namespace Pythia8 {

//==========================================================================

// ampLedS (amplitude) method for LED graviton tree level exchange.
// Based on Eq. (8) in JHEP 1105 (2011) 092, arXiv:1101.4919.

complex ampLedS(double x, double n, double L, double M) {

  complex cS(0., 0.);
  if (n <= 0) return cS;

  // Constants.
  double exp1 = n - 2;
  double exp2 = n + 2;
  double rC = sqrt(pow(M_PI,n)) * pow(L,exp1)
            / (GammaReal(n/2.) * pow(M,exp2));

  // Base functions, F1 and F2.
  complex I(0., 1.);
  if (x < 0) {
    double sqrX = sqrt(-x);
    if (int(n) % 2 == 0) {
      cS = -log(fabs(1 - 1/x));
    } else {
      cS = (2.*atan(sqrX) - M_PI)/sqrX;
    }
  } else if ((x > 0) && (x < 1)) {
    double sqrX = sqrt(x);
    if (int(n) % 2 == 0) {
      cS = -log(fabs(1 - 1/x)) - M_PI*I;
    } else {
      double rat = (sqrX + 1)/(sqrX - 1);
      cS = log(fabs(rat))/sqrX - M_PI*I/sqrX;
    }
  } else if (x > 1){
    double sqrX = sqrt(x);
    if (int(n) % 2 == 0) {
      cS = -log(fabs(1 - 1/x));
    } else {
      double rat = (sqrX + 1)/(sqrX - 1);
      cS = log(fabs(rat))/sqrX;
    }
  }

  // Recursive part.
  int nL;
  int nD;
  if (int(n) % 2 == 0) {
    nL = int(n/2.);
    nD = 2;
  } else {
    nL = int((n + 1)/2.);
    nD = 1;
  }
  for (int i=1; i<nL; ++i) {
    cS = x*cS - 2./nD;
    nD += 2;
  }

  return rC*cS;
}

//--------------------------------------------------------------------------

// Common method, "Mandelstam polynomial", for LED dijet processes.

double funLedG(double x, double y) {
  double ret = pow(x,4) + 10. * pow(x,3) * y + 42. * pow2(x) * pow2(y)
             + 64. * x * pow(y,3) + 32. * pow(y,4);
  return ret;
}

//==========================================================================

// Sigma1gg2GravitonStar class.
// Cross section for g g -> G* (excited graviton state).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1gg2GravitonStar::initProc() {

  // Store G* mass and width for propagator.
  idGstar  = 5100039;
  mRes     = particleDataPtr->m0(idGstar);
  GammaRes = particleDataPtr->mWidth(idGstar);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // SMinBulk = off/on, use universal coupling (kappaMG)
  // or individual (Gxx) between graviton and SM particles.
  eDsmbulk   = settingsPtr->flag("ExtraDimensionsG*:SMinBulk");
  eDvlvl = false;
  if (eDsmbulk) eDvlvl = settingsPtr->flag("ExtraDimensionsG*:VLVL");
  kappaMG    = settingsPtr->parm("ExtraDimensionsG*:kappaMG");
  for (int i = 0; i < 27; ++i) eDcoupling[i] = 0.;
  double tmPcoup = settingsPtr->parm("ExtraDimensionsG*:Gqq");
  for (int i = 1; i <= 4; ++i)  eDcoupling[i] = tmPcoup;
  eDcoupling[5] = settingsPtr->parm("ExtraDimensionsG*:Gbb");
  eDcoupling[6] = settingsPtr->parm("ExtraDimensionsG*:Gtt");
  tmPcoup = settingsPtr->parm("ExtraDimensionsG*:Gll");
  for (int i = 11; i <= 16; ++i) eDcoupling[i] = tmPcoup;
  eDcoupling[21] = settingsPtr->parm("ExtraDimensionsG*:Ggg");
  eDcoupling[22] = settingsPtr->parm("ExtraDimensionsG*:Ggmgm");
  eDcoupling[23] = settingsPtr->parm("ExtraDimensionsG*:GZZ");
  eDcoupling[24] = settingsPtr->parm("ExtraDimensionsG*:GWW");
  eDcoupling[25] = settingsPtr->parm("ExtraDimensionsG*:Ghh");

  // Set pointer to particle properties and decay table.
  gStarPtr = particleDataPtr->particleDataEntryPtr(idGstar);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1gg2GravitonStar::sigmaKin() {

  // Incoming width for gluons.
  double widthIn  = mH / (160. * M_PI);

  // RS graviton coupling
  if (eDsmbulk) widthIn *= 2. * pow2(eDcoupling[21] * mH);
  else          widthIn *= pow2(kappaMG * mH / mRes);

  // Set up Breit-Wigner. Width out only includes open channels.
  double sigBW    = 5. * M_PI/ ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double widthOut = gStarPtr->resWidthOpen(idGstar, mH);

  // Modify cross section in wings of peak. Done.
  sigma           = widthIn * sigBW * widthOut;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1gg2GravitonStar::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, idGstar);

  // Colour flow topology.
  setColAcol( 1, 2, 2, 1, 0, 0);

}

//--------------------------------------------------------------------------

// Evaluate weight for G* decay angle.
// SA: Angle dist. for decay G* -> W/Z/h, based on
// Phys.Rev. D65 (2002) 075008, [arXiv:hep-ph/0103308v3]

double Sigma1gg2GravitonStar::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // G* should sit in entry 5.
  if (iResBeg != 5 || iResEnd != 5) return 1.;

  // Phase space factors. Reconstruct decay angle.
  double mr1    = pow2(process[6].m()) / sH;
  double mr2    = pow2(process[7].m()) / sH;
  double betaf  = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2);
  double cosThe = (process[3].p() - process[4].p())
    * (process[7].p() - process[6].p()) / (sH * betaf);

  // Default is isotropic decay.
  double wt     = 1.;

  // Angular weight for g + g -> G* -> f + fbar.
  if (process[6].idAbs() < 19) {
    wt = 1. - pow4(cosThe);

  // Angular weight for g + g -> G* -> g + g or gamma + gamma.
  } else if (process[6].id() == 21 || process[6].id() == 22) {
    wt = (1. + 6. * pow2(cosThe) + pow4(cosThe)) / 8.;

  // Angular weight for g + g -> G* -> Z + Z or W + W.
  } else if (process[6].id() == 23 || process[6].id() == 24) {
    double beta2 = pow2(betaf);
    double cost2 = pow2(cosThe);
    double cost4 = pow2(cost2);
    wt = pow2(beta2 - 2.)*(1. - 2.*cost2 + cost4);
    // Longitudinal W/Z only.
    if(eDvlvl) {
      wt /= 4.;
    // Transverse W/Z contributions as well.
    } else {
      double beta4 = pow2(beta2);
      double beta8 = pow2(beta4);
      wt += 2.*pow2(beta4 - 1.)*beta4*cost4;
      wt += 2.*pow2(beta2 - 1.)*(1. - 2.*beta4*cost2 + beta8*cost4);
      wt += 2.*(1. + 6.*beta4*cost2 + beta8*cost4);
      wt += 8.*(1. - beta2)*(1. - cost4);
      wt /= 18.;
    }

  // Angular weight for g + g -> G* -> h + h
  } else if (process[6].id() == 25) {
    double beta2 = pow2(betaf);
    double cost2 = pow2(cosThe);
    double cost4 = pow2(cost2);
    wt = pow2(beta2 - 2.)*(1. - 2.*cost2 + cost4);
    wt /= 4.;
  }

  // Done.
  return wt;

}

//==========================================================================

// Sigma1ffbar2GravitonStar class.
// Cross section for f fbar -> G* (excited graviton state).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ffbar2GravitonStar::initProc() {

  // Store G* mass and width for propagator.
  idGstar  = 5100039;
  mRes     = particleDataPtr->m0(idGstar);
  GammaRes = particleDataPtr->mWidth(idGstar);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // SMinBulk = off/on, use universal coupling (kappaMG)
  // or individual (Gxx) between graviton and SM particles.
  eDsmbulk   = settingsPtr->flag("ExtraDimensionsG*:SMinBulk");
  eDvlvl = false;
  if (eDsmbulk) eDvlvl = settingsPtr->flag("ExtraDimensionsG*:VLVL");
  kappaMG    = settingsPtr->parm("ExtraDimensionsG*:kappaMG");
  for (int i = 0; i < 27; ++i) eDcoupling[i] = 0.;
  double tmPcoup = settingsPtr->parm("ExtraDimensionsG*:Gqq");
  for (int i = 1; i <= 4; ++i)  eDcoupling[i] = tmPcoup;
  eDcoupling[5] = settingsPtr->parm("ExtraDimensionsG*:Gbb");
  eDcoupling[6] = settingsPtr->parm("ExtraDimensionsG*:Gtt");
  tmPcoup = settingsPtr->parm("ExtraDimensionsG*:Gll");
  for (int i = 11; i <= 16; ++i) eDcoupling[i] = tmPcoup;
  eDcoupling[21] = settingsPtr->parm("ExtraDimensionsG*:Ggg");
  eDcoupling[22] = settingsPtr->parm("ExtraDimensionsG*:Ggmgm");
  eDcoupling[23] = settingsPtr->parm("ExtraDimensionsG*:GZZ");
  eDcoupling[24] = settingsPtr->parm("ExtraDimensionsG*:GWW");
  eDcoupling[25] = settingsPtr->parm("ExtraDimensionsG*:Ghh");

  // Set pointer to particle properties and decay table.
  gStarPtr = particleDataPtr->particleDataEntryPtr(idGstar);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1ffbar2GravitonStar::sigmaKin() {

  // Incoming width for fermions, disregarding colour factor.
  double widthIn  = mH / (80. * M_PI);

  // Set up Breit-Wigner. Width out only includes open channels.
  double sigBW    = 5. * M_PI/ ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double widthOut = gStarPtr->resWidthOpen(idGstar, mH);

  // Modify cross section in wings of peak. Done.
  sigma0          = widthIn * sigBW * widthOut;
}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part dependent of incoming flavour.

double Sigma1ffbar2GravitonStar::sigmaHat() {

  double sigma = sigma0;

  // RS graviton coupling
  if (eDsmbulk) sigma *= 2. * pow2(eDcoupling[min( abs(id1), 26)] * mH);
  else          sigma *= pow2(kappaMG * mH / mRes);

  // If initial quarks, 1/N_C
  if (abs(id1) < 9) sigma /= 3.;

  return sigma;
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2GravitonStar::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idGstar);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for G* decay angle.
// SA: Angle dist. for decay G* -> W/Z/h, based on
// Phys.Rev. D65 (2002) 075008, [arXiv:hep-ph/0103308v3]

double Sigma1ffbar2GravitonStar::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // G* should sit in entry 5.
  if (iResBeg != 5 || iResEnd != 5) return 1.;

  // Phase space factors. Reconstruct decay angle.
  double mr1    = pow2(process[6].m()) / sH;
  double mr2    = pow2(process[7].m()) / sH;
  double betaf  = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2);
  double cosThe = (process[3].p() - process[4].p())
    * (process[7].p() - process[6].p()) / (sH * betaf);

  // Default is isotropic decay.
  double wt     = 1.;

  // Angular weight for f + fbar -> G* -> f + fbar.
  if (process[6].idAbs() < 19) {
    wt = (1. - 3. * pow2(cosThe) + 4. * pow4(cosThe)) / 2.;

  // Angular weight for f + fbar -> G* -> g + g or gamma + gamma.
  } else if (process[6].id() == 21 || process[6].id() == 22) {
    wt = 1. - pow4(cosThe);

  // Angular weight for f + fbar -> G* -> Z + Z or W + W.
  }  else if (process[6].id() == 23 || process[6].id() == 24) {
    double beta2 = pow2(betaf);
    double cost2 = pow2(cosThe);
    double cost4 = pow2(cost2);
    wt = pow2(beta2 - 2.)*cost2*(1. - cost2);
    // Longitudinal W/Z only.
    if (eDvlvl) {
      wt /= 4.;
    // Transverse W/Z contributions as well.
    } else {
      wt += pow2(beta2 - 1.)*cost2*(1. - cost2);
      wt += 2.*(1. - cost4);
      wt += (1. - beta2)*(1. - 3.*cost2 + 4.*cost4);
      wt /= 8.;
    }

  // Angular weight for f + fbar -> G* -> h + h
  } else if (process[6].id() == 25) {
    double beta2 = pow2(betaf);
    double cost2 = pow2(cosThe);
    wt = pow2(beta2 - 2.)*cost2*(1. - cost2);
    wt /= 4.;
  }

  // Done.
  return wt;

}

//==========================================================================

// Sigma1qqbar2KKgluonStar class.
// Cross section for q qbar -> g^*/KK-gluon^* (excited KK-gluon state).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1qqbar2KKgluonStar::initProc() {

  // Store kk-gluon* mass and width for propagator.
  idKKgluon = 5100021;
  mRes      = particleDataPtr->m0(idKKgluon);
  GammaRes  = particleDataPtr->mWidth(idKKgluon);
  m2Res     = mRes*mRes;
  GamMRat   = GammaRes / mRes;

  // KK-gluon gv/ga couplings and interference.
  for (int i = 0; i < 10; ++i) { eDgv[i] = 0.; eDga[i] = 0.; }
  double tmPgL = settingsPtr->parm("ExtraDimensionsG*:KKgqL");
  double tmPgR = settingsPtr->parm("ExtraDimensionsG*:KKgqR");
  for (int i = 1; i <= 4; ++i) {
    eDgv[i] = 0.5 * (tmPgL + tmPgR);
    eDga[i] = 0.5 * (tmPgL - tmPgR);
  }
  tmPgL = settingsPtr->parm("ExtraDimensionsG*:KKgbL");
  tmPgR = settingsPtr->parm("ExtraDimensionsG*:KKgbR");
  eDgv[5] = 0.5 * (tmPgL + tmPgR); eDga[5] = 0.5 * (tmPgL - tmPgR);
  tmPgL = settingsPtr->parm("ExtraDimensionsG*:KKgtL");
  tmPgR = settingsPtr->parm("ExtraDimensionsG*:KKgtR");
  eDgv[6] = 0.5 * (tmPgL + tmPgR); eDga[6] = 0.5 * (tmPgL - tmPgR);
  interfMode    = settingsPtr->mode("ExtraDimensionsG*:KKintMode");

  // Set pointer to particle properties and decay table.
  gStarPtr = particleDataPtr->particleDataEntryPtr(idKKgluon);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1qqbar2KKgluonStar::sigmaKin() {

  // Incoming width for fermions.
  double widthIn  = alpS * mH * 4 / 27;
  double widthOut = alpS * mH / 6;

  // Loop over all decay channels.
  sumSM  = 0.;
  sumInt = 0.;
  sumKK  = 0.;

  for (int i = 0; i < gStarPtr->sizeChannels(); ++i) {
    int idAbs = abs( gStarPtr->channel(i).product(0) );

    // Only contributions quarks.
    if ( idAbs > 0 && idAbs <= 6 ) {
      double mf = particleDataPtr->m0(idAbs);

      // Check that above threshold. Phase space.
      if (mH > 2. * mf + MASSMARGIN) {
        double mr    = pow2(mf / mH);
        double beta  = sqrtpos(1. - 4. * mr);

        // Store sum of combinations. For outstate only open channels.
        int onMode = gStarPtr->channel(i).onMode();
        if (onMode == 1 || onMode == 2) {
          sumSM  += beta * (1. + 2. * mr);
          sumInt += beta * eDgv[min(idAbs, 9)] * (1. + 2. * mr);
          sumKK  += beta * (pow2(eDgv[min(idAbs, 9)]) * (1. + 2.*mr)
                          + pow2(eDga[min(idAbs, 9)]) * (1. - 4.*mr));
        }
      }
    }
  }

  // Set up Breit-Wigner. Width out only includes open channels.
  sigSM  = widthIn * 12. * M_PI *  widthOut / sH2;
  sigInt = 2. * sigSM * sH * (sH - m2Res)
         / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  sigKK  = sigSM * sH2 / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Optionally only keep g* or gKK term.
  if (interfMode == 1) {sigInt = 0.; sigKK = 0.;}
  if (interfMode == 2) {sigSM  = 0.; sigInt = 0.;}

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part dependent of incoming flavour.

double Sigma1qqbar2KKgluonStar::sigmaHat() {

  // RS graviton coupling.
  double sigma = sigSM * sumSM
               + eDgv[min(abs(id1), 9)] * sigInt * sumInt
               + ( pow2(eDgv[min(abs(id1), 9)])
                 + pow2(eDga[min(abs(id1), 9)]) ) * sigKK * sumKK;

  return sigma;
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1qqbar2KKgluonStar::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idKKgluon);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 1, 2);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for KK-gluon* decay angle (based on ffbar2gmZ).

double Sigma1qqbar2KKgluonStar::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // g* should sit in entry 5.
  if (iResBeg != 5 || iResEnd != 5) return 1.;

  // Couplings for in- and out-flavours (alpS already included).
  int idInAbs  = process[3].idAbs();
  double vi    = eDgv[min(idInAbs, 9)];
  double ai    = eDga[min(idInAbs, 9)];
  int idOutAbs = process[6].idAbs();
  double vf    = eDgv[min(idOutAbs, 9)];
  double af    = eDga[min(idOutAbs, 9)];

  // Phase space factors. (One power of beta left out in formulae.)
  double mf    = process[6].m();
  double mr    = mf*mf / sH;
  double betaf = sqrtpos(1. - 4. * mr);

  // Coefficients of angular expression.
  double coefTran = sigSM + vi * sigInt * vf
    + (vi*vi + ai*ai) * sigKK * (vf*vf + pow2(betaf) * af*af);
  double coefLong = 4. * mr * ( sigSM + vi * sigInt * vf
                              + (vi*vi + ai*ai) * sigKK * vf*vf );
  double coefAsym = betaf * ( ai * sigInt * af
    + 4. * vi * ai * sigKK * vf * af );

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

// Sigma2gg2GravitonStarg class.
// Cross section for g g -> G* g (excited graviton state).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2GravitonStarg::initProc() {

  // Store G* mass and width for propagator.
  idGstar  = 5100039;
  mRes     = particleDataPtr->m0(idGstar);
  GammaRes = particleDataPtr->mWidth(idGstar);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Overall coupling strength kappa * m_G*.
  kappaMG  = settingsPtr->parm("ExtraDimensionsG*:kappaMG");

   // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(idGstar);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2gg2GravitonStarg::sigmaKin() {

  //  Evaluate cross section. Secondary width for G*.
  sigma = (3. * pow2(kappaMG) * alpS) / (32. * sH * m2Res)
    * ( pow2(tH2 + tH * uH + uH2) / (sH2 * tH * uH)
    + 2. * (tH2 / uH + uH2 / tH) / sH + 3. * (tH / uH + uH / tH)
    + 2. * (sH / uH + sH/tH) + sH2 / (tH * uH) );
  sigma *= openFrac;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2GravitonStarg::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, idGstar, 21);

  // Colour flow topologies: random choice between two mirrors.
  if (rndmPtr->flat() < 0.5) setColAcol( 1, 2, 2, 3, 0, 0, 1, 3);
  else                    setColAcol( 1, 2, 3, 1, 0, 0, 3, 2);

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles: currently G* assumed isotropic.

double Sigma2gg2GravitonStarg::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // No equations for G* decay so assume isotropic.
  return 1.;

}

//==========================================================================

// Sigma2qg2GravitonStarq class.
// Cross section for q g -> G* q (excited graviton state).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qg2GravitonStarq::initProc() {

  // Store G* mass and width for propagator.
  idGstar  = 5100039;
  mRes     = particleDataPtr->m0(idGstar);
  GammaRes = particleDataPtr->mWidth(idGstar);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Overall coupling strength kappa * m_G*.
  kappaMG  = settingsPtr->parm("ExtraDimensionsG*:kappaMG");

   // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(idGstar);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2qg2GravitonStarq::sigmaKin() {

  //  Evaluate cross section. Secondary width for G*.
  sigma = -(pow2(kappaMG) * alpS) / (192. * sH * m2Res )
    * ( 4. * (sH2 + uH2) / (tH * sH) + 9. * (sH + uH) / sH + sH / uH
    + uH2 / sH2 + 3. * tH * (4. + sH / uH + uH / sH) / sH
    + 4. * tH2 * (1. / uH + 1. / sH) / sH + 2. * tH2 * tH / (uH * sH2) );
  sigma *= openFrac;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2GravitonStarq::setIdColAcol() {

  // Flavour set up for q g -> H q.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, idGstar, idq);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21);

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles: currently G* assumed isotropic.

double Sigma2qg2GravitonStarq::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // No equations for G* decay so assume isotropic.
  return 1.;

}

//==========================================================================

// Sigma2qqbar2GravitonStarg class.
// Cross section for q qbar -> G* g (excited graviton state).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qqbar2GravitonStarg::initProc() {

  // Store G* mass and width for propagator.
  idGstar  = 5100039;
  mRes     = particleDataPtr->m0(idGstar);
  GammaRes = particleDataPtr->mWidth(idGstar);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Overall coupling strength kappa * m_G*.
  kappaMG  = settingsPtr->parm("ExtraDimensionsG*:kappaMG");

   // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(idGstar);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma2qqbar2GravitonStarg::sigmaKin() {

  // Evaluate cross section. Secondary width for G*.
  sigma = (pow2(kappaMG) * alpS) / (72. * sH * m2Res)
    * ( 4. * (tH2 + uH2) / sH2 + 9. * (tH + uH) / sH
    + (tH2 / uH + uH2 / tH) / sH + 3. * (4. + tH / uH + uH/ tH)
    + 4. * (sH / uH + sH / tH) + 2. * sH2 / (tH * uH) );
  sigma *= openFrac;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2GravitonStarg::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idGstar, 21);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles: currently G* assumed isotropic.

double Sigma2qqbar2GravitonStarg::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // No equations for G* decay so assume isotropic.
  return 1.;

}

//==========================================================================

// NOAM: Sigma2ffbar2TEVffbar class.
// Cross section for, f fbar -> gammaKK/ZKK -> F Fbar.
// Process provided by N. Hod et al. and is described in arXiv:XXXX.YYYY

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2ffbar2TEVffbar::initProc() {

  // Process name.
  if (idNew == 1) nameSave = "f fbar -> d dbar (s-channel gamma_KK/Z_KK)";
  if (idNew == 2) nameSave = "f fbar -> u ubar (s-channel gamma_KK/Z_KK)";
  if (idNew == 3) nameSave = "f fbar -> s sbar (s-channel gamma_KK/Z_KK)";
  if (idNew == 4) nameSave = "f fbar -> c cbar (s-channel gamma_KK/Z_KK)";
  if (idNew == 5) nameSave = "f fbar -> b bbar (s-channel gamma_KK/Z_KK)";
  if (idNew == 6) nameSave = "f fbar -> t tbar (s-channel gamma_KK/Z_KK)";
  if (idNew == 11) nameSave = "f fbar -> e+ e- (s-channel gamma_KK/Z_KK)";
  if (idNew == 12) nameSave = "f fbar -> nue nuebar (s-channel gamma_KK/Z_KK)";
  if (idNew == 13) nameSave = "f fbar -> mu+ mu- (s-channel gamma_KK/Z_KK)";
  if (idNew == 14) nameSave
    = "f fbar -> numu numubar (s-channel gamma_KK/Z_KK)";
  if (idNew == 15) nameSave = "f fbar -> tau+ tau- (s-channel gamma_KK/Z_KK)";
  if (idNew == 16) nameSave
    = "f fbar -> nutau nutaubar (s-channel gamma_KK/Z_KK)";

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode = settingsPtr->mode("ExtraDimensionsTEV:gmZmode");

  // Pick number of KK excitations
  nexcitationmax  = settingsPtr->mode("ExtraDimensionsTEV:nMax");

  // Initialize the widths of the KK propogators.
  // partial width of the KK photon
  wgmKKFactor = 0.;
  // total width of the KK photon
  wgmKKn      = 0.;
  // will be proportional to "wZ0" + ttbar addition
  wZKKn       = 0.;

  // Store Z0 mass and width for propagator.
  wZ0 = particleDataPtr->mWidth(23);
  mRes  = particleDataPtr->m0(23);
  m2Res = mRes*mRes;

  // Store the top mass for the ttbar width calculations
  mTop  = particleDataPtr->m0(6);
  m2Top = mTop*mTop;

  // Store the KK mass parameter, equivalent to the mass of the first KK
  // excitation: particleDataPtr->m0(5000023);
  mStar = (double)settingsPtr->parm("ExtraDimensionsTEV:mStar");

  // Get alphaEM - relevant for the calculation of the widths
  alphaemfixed = settingsPtr->parm("StandardModel:alphaEM0");

  // initialize imaginari number
  mI = complex(0.,1.);

  // Sum all partial widths of the KK photon except for the ttbar channel
  // which is handeled afterwards seperately
  if (gmZmode>=0 && gmZmode<=5) {
    for (int i=1 ; i<17 ; i++) {
      if (i==7) { i=11; }
      // skip the ttbar decay and add its contribution later
      if (i==6) { continue; }
      if (i<9) {
        wgmKKFactor += ( (alphaemfixed / 6.) * 4.
                    * couplingsPtr->ef(i) * couplingsPtr->ef(i) * 3. );
      }
      else {
        wgmKKFactor += (alphaemfixed / 6.) * 4.
                    * couplingsPtr->ef(i) * couplingsPtr->ef(i);
      }
    }
  }

  // Get the helicity-couplings of the Z0 to all the fermions except top
  gMinusF  = ( couplingsPtr->t3f(idNew) - couplingsPtr->ef(idNew)
           * couplingsPtr->sin2thetaW() )
           / sqrt( couplingsPtr->sin2thetaW()*couplingsPtr->cos2thetaW() );
  gPlusF   = -1. * couplingsPtr->ef(idNew) * couplingsPtr->sin2thetaW()
           / sqrt( couplingsPtr->sin2thetaW() * couplingsPtr->cos2thetaW() );
  // Get the helicity-couplings of the Z0 to the top quark
  gMinusTop  = ( couplingsPtr->t3f(6) - couplingsPtr->ef(6)
             * couplingsPtr->sin2thetaW() )
             / sqrt( couplingsPtr->sin2thetaW()*couplingsPtr->cos2thetaW() );

  gPlusTop   = -1. * couplingsPtr->ef(6) * couplingsPtr->sin2thetaW()
             / sqrt( couplingsPtr->sin2thetaW() * couplingsPtr->cos2thetaW() );
  // calculate the constant factor of the unique ttbar decay width
  ttbarwFactorA = pow2(gMinusTop) + pow2(gPlusTop);
  ttbarwFactorB = 6.*gMinusTop*gPlusTop - pow2(gMinusTop) - pow2(gPlusTop);

  // Secondary open width fraction, relevant for top (or heavier).
  openFracPair = 1.;
  if ((idNew >=6 && idNew <=8) || idNew == 17 || idNew == 18)
    openFracPair = particleDataPtr->resOpenFrac(idNew, -idNew);

}

//--------------------------------------------------------------------------

// For improving the phase-space sampling (there can be 2 resonances)

int Sigma2ffbar2TEVffbar::resonanceB() {

  return 23;

}

//--------------------------------------------------------------------------

// For improving the phase-space sampling (there can be 2 resonances)

int Sigma2ffbar2TEVffbar::resonanceA() {

  if (gmZmode>=3) {
    phaseSpacemHatMin  = settingsPtr->parm("PhaseSpace:mHatMin");
    phaseSpacemHatMax  = settingsPtr->parm("PhaseSpace:mHatMax");
    double mResFirstKKMode = sqrt(pow2(particleDataPtr->m0(23)) + pow2(mStar));
    if (mResFirstKKMode/2. <= phaseSpacemHatMax
        || 3*mResFirstKKMode/2. >= phaseSpacemHatMin) { return 5000023; }
    else { return 23; }
  // no KK terms at all
  } else { return 23; }

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2ffbar2TEVffbar::sigmaKin() {

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

  // Reconstruct decay angle so can reuse 2 -> 1 cross section.
  cosThe         = (tH - uH) / (betaf * sH);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2ffbar2TEVffbar::sigmaHat() {

  // Fail if below threshold.
  if (!isPhysical) return 0.;

  // Couplings for in/out-flavours.
  int idAbs = abs(id1);

  // The couplings of the Z0 to the fermions for in/out flavors
  gMinusf  = ( couplingsPtr->t3f(idAbs) - couplingsPtr->ef(idAbs)
               * couplingsPtr->sin2thetaW() )
           / sqrt( couplingsPtr->sin2thetaW()*couplingsPtr->cos2thetaW() );
  gPlusf   = -1. * couplingsPtr->ef(idAbs)*couplingsPtr->sin2thetaW()
           / sqrt( couplingsPtr->sin2thetaW()*couplingsPtr->cos2thetaW() );

  // Initialize the some values
  helicityME2 = 0.;
  coefAngular = 0.;
  gf=0.;
  gF=0.;
  gammaProp = complex(0.,0.);
  resProp   = complex(0.,0.);
  gmPropKK  = complex(0.,0.);
  ZPropKK   = complex(0.,0.);
  totalProp = complex(0.,0.);

  // Sum all initial and final helicity states this corresponds to an
  // unpolarized beams and unmeasured polarization final-state
  for (double helicityf=-0.5 ; helicityf<=0.5 ; helicityf++) {
    for (double helicityF=-0.5 ; helicityF<=0.5 ; helicityF++) {
          // the couplings for the initial-final helicity configuration
      gF = (helicityF == +0.5) ? gMinusF : gPlusF;
      gf = (helicityf == +0.5) ? gMinusf : gPlusf;
      // 0=SM gmZ,  1=SM gm,  2=SM Z,  3=SM+KK gmZ,  4=KK gm,  5=KK Z
      switch(gmZmode) {
        // SM photon and Z0 only
        case 0:
          gammaProp = couplingsPtr->ef(idAbs)*couplingsPtr->ef(idNew)/sH;
          resProp   = gf*gF/( sH - m2Res + mI*sH*(wZ0/mRes) );
          break;
        // SM photon only
        case 1:
          gammaProp = couplingsPtr->ef(idAbs)*couplingsPtr->ef(idNew)/sH;
          break;
        // SM Z0 only
        case 2:
          resProp   = gf*gF/( sH - m2Res + mI*sH*(wZ0/mRes) );
          break;
        // KK photon and Z
        case 3:
          gammaProp = couplingsPtr->ef(idAbs)*couplingsPtr->ef(idNew)/sH;
          resProp   = gf*gF/( sH - m2Res + mI*sH*(wZ0/mRes) );
          ZPropKK   = complex(0.,0.);
          gmPropKK  = complex(0.,0.);
                  // Sum all KK excitations contributions
          for(int nexcitation = 1; nexcitation <= nexcitationmax;
            nexcitation++) {
            mZKKn   = sqrt(m2Res + pow2(mStar * nexcitation));
            m2ZKKn  = m2Res + pow2(mStar * nexcitation);
            mgmKKn  = mStar * nexcitation;
            m2gmKKn = (mStar*nexcitation)*(mStar*nexcitation);
            // calculate the width of the n'th excitation of the KK Z
            // (proportional to the Z0 width + ttbar partial width)
            ttbarwZKKn = 2.*(alphaemfixed*3./6.)*mZKKn
                        * sqrt(1.-4.*m2Top/m2ZKKn)
                        * (ttbarwFactorA+(m2Top/m2ZKKn)*ttbarwFactorB);
            wZKKn       = 2.*wZ0*mZKKn/mRes+ttbarwZKKn;
            // calculate the width of the n'th excitation of the
            // KK photon
            ttbarwgmKKn = 2.*(alphaemfixed*3./6.)*mgmKKn
                        * sqrt(1.-4.*m2Top/m2gmKKn)
                        * 2.*pow2(couplingsPtr->ef(6))*(1.+2.*(m2Top/m2gmKKn));
            wgmKKn       = wgmKKFactor*mgmKKn+ttbarwgmKKn;
            // the propogators
            gmPropKK += (2.*couplingsPtr->ef(idAbs)*couplingsPtr->ef(idNew))
                      / (sH-m2gmKKn+mI*sH*wgmKKn/mgmKKn);
            ZPropKK  += (2.*gf*gF)/(sH-m2ZKKn+mI*sH*wZKKn/mZKKn );
          }
          break;
        // SM photon and Z0 with KK photon only
        case 4:
          gammaProp = couplingsPtr->ef(idAbs)*couplingsPtr->ef(idNew)/sH;
          resProp   = gf*gF/( sH - m2Res + mI*sH*(wZ0/mRes) );
          gmPropKK  = complex(0.,0.);
          for (int nexcitation = 1; nexcitation <= nexcitationmax;
            nexcitation++ ) {
            mgmKKn  = mStar * nexcitation;
            m2gmKKn = (mStar*nexcitation)*(mStar*nexcitation);

            ttbarwgmKKn = 2.*(alphaemfixed*3./6.)*mgmKKn
                        * sqrt(1.-4.*m2Top/m2gmKKn)
                        * 2.*pow2(couplingsPtr->ef(6))
                        * (1.+2.*(m2Top/m2gmKKn));
            wgmKKn         = wgmKKFactor*mgmKKn+ttbarwgmKKn;
            gmPropKK += (2.*couplingsPtr->ef(idAbs)*couplingsPtr->ef(idNew))
                      / (sH-m2gmKKn+mI*sH*wgmKKn/mgmKKn);
          }
          break;
        // SM photon and Z0 with KK Z only
        case 5:
          gammaProp = couplingsPtr->ef(idAbs)*couplingsPtr->ef(idNew)/sH;
          resProp   = gf*gF/( sH - m2Res + mI*sH*(wZ0/mRes) );
          ZPropKK   = complex(0.,0.);
          for (int nexcitation = 1; nexcitation <= nexcitationmax;
            nexcitation++ ) {
            mZKKn   = sqrt(m2Res + pow2(mStar * nexcitation));
            m2ZKKn  = m2Res + pow2(mStar * nexcitation);

            ttbarwZKKn = 2.*(alphaemfixed*3./6.)*mZKKn
                          * sqrt(1.-4.*m2Top/m2ZKKn)
                          * (ttbarwFactorA+(m2Top/m2ZKKn)*ttbarwFactorB);
            wZKKn       = 2.*wZ0*mZKKn/mRes+ttbarwZKKn;
            ZPropKK    += (2.*gf*gF)/(sH-m2ZKKn+mI*sH*wZKKn/mZKKn );
          }
          break;
        default: break;
      // end run over initial and final helicity states
      }

          // sum all contributing amplitudes
      totalProp = gammaProp + resProp + ZPropKK + gmPropKK;

          // angular distribution for the helicity configuration
      coefAngular = 1. + 4. * helicityF * helicityf * cosThe;

          // the squared helicity matrix element
      helicityME2 += real(totalProp*conj(totalProp))*pow2(coefAngular);
    }
  }

  // calculate the coefficient of the squared helicity matrix element.
  coefTot = (2./sH) * 2*M_PI * pow2(alpEM)/(4.*sH) * pow2(sH)/4.;

  // the full squared helicity matrix element.
  double sigma = helicityME2 * coefTot;

  // Top: corrections for closed decay channels.
  sigma *= openFracPair;

  // Initial-state colour factor. Answer.
  if (idAbs < 9) sigma /= 3.;

  // Final-state colour factor. Answer.
  if (idNew < 9) sigma *= 3.*(1.+alpS/M_PI);

  return sigma;
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2TEVffbar::setIdColAcol() {

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

double Sigma2ffbar2TEVffbar::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // For top decay hand over to standard routine, else done.
  if (idNew == 6 && process[process[iResBeg].mother1()].idAbs() == 6)
       return weightTopDecay( process, iResBeg, iResEnd);
  else return 1.;

}

//==========================================================================

// Sigma2gg2LEDUnparticleg class.
// Cross section for g g -> U/G g (real graviton emission in
// large extra dimensions or unparticle emission).

//--------------------------------------------------------------------------

void Sigma2gg2LEDUnparticleg::initProc() {

  // Init model parameters.
  eDidG    = 5000039;
  if (eDgraviton) {
    eDspin     = (settingsPtr->flag("ExtraDimensionsLED:GravScalar")) ? 0 : 2;
    eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
    eDdU       = 0.5 * eDnGrav + 1;
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsLED:MD");
    eDlambda   = 1;
    eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
    eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");
    eDcf       = settingsPtr->parm("ExtraDimensionsLED:c");
  } else {
    eDspin     = settingsPtr->mode("ExtraDimensionsUnpart:spinU");
    eDdU       = settingsPtr->parm("ExtraDimensionsUnpart:dU");
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsUnpart:LambdaU");
    eDlambda   = settingsPtr->parm("ExtraDimensionsUnpart:lambda");
    eDcutoff   = settingsPtr->mode("ExtraDimensionsUnpart:CutOffMode");
  }

  // The A(dU) or S'(n) value.
  double tmpAdU = 0;
  if (eDgraviton) {
    tmpAdU  = 2 * M_PI * sqrt( pow(M_PI, double(eDnGrav)) )
            / GammaReal(0.5 * eDnGrav);
    if (eDspin == 0) {  // Scalar graviton
      tmpAdU *= sqrt( pow(2., double(eDnGrav)) );
      eDcf   *= eDcf;
    }
  } else {
    tmpAdU = 16 * pow2(M_PI) * sqrt(M_PI) / pow(2. * M_PI, 2. * eDdU)
      * GammaReal(eDdU + 0.5) / (GammaReal(eDdU - 1.) * GammaReal(2. * eDdU));
  }

  // Cross section related constants
  // and ME dependent powers of lambda / LambdaU.
  double tmpExp   = eDdU - 2;
  double tmpLS    = pow2(eDLambdaU);
  eDconstantTerm = tmpAdU / (2 * 16 * pow2(M_PI) * tmpLS * pow(tmpLS,tmpExp));
  if (eDgraviton) {
    eDconstantTerm /= tmpLS;
  } else if (eDspin == 0) {
    eDconstantTerm *= pow2(eDlambda) / tmpLS;
  } else {
    eDconstantTerm = 0;
    infoPtr->errorMsg("Error in Sigma2gg2LEDUnparticleg::initProc: "
                      "Incorrect spin value (turn process off)!");
  }

}

//--------------------------------------------------------------------------

void Sigma2gg2LEDUnparticleg::sigmaKin() {

  // Set graviton mass.
  mG        = m3;
  mGS       = mG*mG;

  // Set mandelstam variables and ME expressions.
  if (eDgraviton) {

    double A0 = 1/sH;
    if (eDspin == 0) {  // Scalar graviton
      double tmpTerm1 = uH + tH;
      double tmpTerm2 = uH + sH;
      double tmpTerm3 = tH + sH;
      double T0 = pow(tmpTerm1,4) + pow(tmpTerm2,4) + pow(tmpTerm3,4)
                + 12. * sH * tH * uH * mGS;
      eDsigma0 = eDcf * A0 * T0 / (sH2 * tH * uH);
    } else {
      double xH = tH/sH;
      double yH = mGS/sH;
      double xHS = pow2(xH);
      double yHS = pow2(yH);
      double xHC = pow(xH,3);
      double yHC = pow(yH,3);
      double xHQ = pow(xH,4);
      double yHQ = pow(yH,4);

      double T0 = 1/(xH*(yH-1-xH));
      double T1 = 1 + 2*xH + 3*xHS + 2*xHC + xHQ;
      double T2 = -2*yH*(1 + xHC);
      double T3 = 3*yHS*(1 + xHS);
      double T4 = -2*yHC*(1 + xH);
      double T5 = yHQ;

      eDsigma0 = A0 * T0 *( T1 + T2 + T3 + T4 + T5 );
    }

  } else if (eDspin == 0) {

    double A0  = 1/pow2(sH);
    double sHQ = pow(sH,4);
    double tHQ = pow(tH,4);
    double uHQ = pow(uH,4);

    eDsigma0 = A0 * (pow(mGS,4) + sHQ + tHQ + uHQ) / (sH * tH * uH);

  }

  // Mass measure, (m^2)^(d-2).
  double tmpExp = eDdU - 2;
  eDsigma0 *= pow(mGS, tmpExp);

  // Constants.
  eDsigma0 *= eDconstantTerm;

}

//--------------------------------------------------------------------------

double Sigma2gg2LEDUnparticleg::sigmaHat() {

  // Mass spectrum weighting.
  double sigma = eDsigma0 / runBW3;

  // SM couplings...
  if (eDgraviton) {
    sigma *= 16 * M_PI * alpS * 3 / 16;
  } else if (eDspin == 0) {
    sigma *= 6 * M_PI * alpS;
  }

  // Truncate sH region or use form factor.
  // Form factor uses either pythia8 renormScale2
  // or E_jet in cms.
  if (eDcutoff == 1) {
    if (sH > pow2(eDLambdaU) ) { sigma *= pow(eDLambdaU,4)/pow2(sH); }
  } else if ( (eDgraviton && (eDspin == 2))
           && ((eDcutoff == 2) || (eDcutoff == 3)) ) {
    double tmPmu = sqrt(Q2RenSave);
    if (eDcutoff == 3) tmPmu = (sH + s4 - s3) / (2 * mH);
    double tmPformfact = tmPmu / (eDtff * eDLambdaU);
    double tmPexp = double(eDnGrav) + 2;
    sigma *= 1 / (1 + pow(tmPformfact, tmPexp));
  }

  return sigma;
}

//--------------------------------------------------------------------------

void Sigma2gg2LEDUnparticleg::setIdColAcol() {

 // Flavours trivial.
  setId( 21, 21, eDidG, 21);

  // Colour flow topologies: random choice between two mirrors.
  if (rndmPtr->flat() < 0.5) setColAcol( 1, 2, 2, 3, 0, 0, 1, 3);
  else                    setColAcol( 1, 2, 3, 1, 0, 0, 3, 2);

}

//==========================================================================

// Sigma2qg2LEDUnparticleq class.
// Cross section for q g -> U/G q (real graviton emission in
// large extra dimensions or unparticle emission).

//--------------------------------------------------------------------------

void Sigma2qg2LEDUnparticleq::initProc() {

  // Init model parameters.
  eDidG    = 5000039;
  if (eDgraviton) {
    eDspin     = (settingsPtr->flag("ExtraDimensionsLED:GravScalar")) ? 0 : 2;
    eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
    eDdU       = 0.5 * eDnGrav + 1;
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsLED:MD");
    eDlambda   = 1;
    eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
    eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");
    eDgf       = settingsPtr->parm("ExtraDimensionsLED:g");
    eDcf       = settingsPtr->parm("ExtraDimensionsLED:c");
  } else {
    eDspin     = settingsPtr->mode("ExtraDimensionsUnpart:spinU");
    eDdU       = settingsPtr->parm("ExtraDimensionsUnpart:dU");
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsUnpart:LambdaU");
    eDlambda   = settingsPtr->parm("ExtraDimensionsUnpart:lambda");
    eDcutoff   = settingsPtr->mode("ExtraDimensionsUnpart:CutOffMode");
  }

  // The A(dU) or S'(n) value.
  double tmpAdU = 0;
  if (eDgraviton) {
    tmpAdU  = 2 * M_PI * sqrt( pow(M_PI, double(eDnGrav)) )
            / GammaReal(0.5 * eDnGrav);
    // Scalar graviton
    if (eDspin == 0) {
      tmpAdU *= 2. * sqrt( pow(2., double(eDnGrav)) );
      eDcf   *= 4. * eDcf / pow2(eDLambdaU);
      double tmpExp = 2. * double(eDnGrav) / (double(eDnGrav) + 2.);
      eDgf   *= eDgf / pow(2. * M_PI, tmpExp);
    }
  } else {
    tmpAdU = 16 * pow2(M_PI) * sqrt(M_PI) / pow(2. * M_PI, 2. * eDdU)
      * GammaReal(eDdU + 0.5) / (GammaReal(eDdU - 1.) * GammaReal(2. * eDdU));
  }

  // Cross section related constants
  // and ME dependent powers of lambda / LambdaU.
  double tmpExp   = eDdU - 2;
  double tmpLS    = pow2(eDLambdaU);
  eDconstantTerm = tmpAdU / (2 * 16 * pow2(M_PI) * tmpLS * pow(tmpLS,tmpExp));
  if (eDgraviton && (eDspin == 2)) {
    eDconstantTerm /= tmpLS;
  } else if (eDspin == 1) {
    eDconstantTerm *= pow2(eDlambda);
  } else if (eDspin == 0) {
    eDconstantTerm *= pow2(eDlambda);
  } else {
    eDconstantTerm = 0;
    infoPtr->errorMsg("Error in Sigma2qg2LEDUnparticleq::initProc: "
                      "Incorrect spin value (turn process off)!");
  }


}

//--------------------------------------------------------------------------

void Sigma2qg2LEDUnparticleq::sigmaKin() {

  // Set graviton mass.
  mG        = m3;
  mGS       = mG*mG;

  // Set mandelstam variables and ME expressions.
  if (eDgraviton) {

    double A0 = 1/sH;
    // Scalar graviton
    if (eDspin == 0) {
      A0 /= sH;
      double T0 = -(uH2 + pow2(mGS)) / (sH * tH);
      double T1 = -(tH2 + sH2)/ uH;
      eDsigma0 = A0 * (eDgf * T0 + eDcf * T1);
    } else {
      double xH = tH/sH;
      double yH = mGS/sH;
      double x2H = xH/(yH - 1 - xH);
      double y2H = yH/(yH - 1 - xH);
      double x2HS = pow2(x2H);
      double y2HS = pow2(y2H);
      double x2HC = pow(x2H,3);
      double y2HC = pow(y2H,3);

      double T0 = -(yH - 1 - xH);
      double T20 = 1/(x2H*(y2H-1-x2H));
      double T21 = -4*x2H*(1 + x2H)*(1 + 2*x2H + 2*x2HS);
      double T22 = y2H*(1 + 6*x2H + 18*x2HS + 16*x2HC);
      double T23 = -6*y2HS*x2H*(1+2*x2H);
      double T24 = y2HC*(1 + 4*x2H);

      eDsigma0 = A0 * T0 * T20 * ( T21 + T22 + T23 + T24 );
    }

  } else if (eDspin == 1) {

    double A0  = 1/pow2(sH);
    double tmpTerm1 = tH - mGS;
    double tmpTerm2 = sH - mGS;

    eDsigma0 = A0 * (pow2(tmpTerm1) + pow2(tmpTerm2)) / (sH*tH);

  } else if (eDspin == 0) {

    double A0  = 1/pow2(sH);
    // Sign correction by Tom
    eDsigma0 = A0 * (pow2(tH) + pow2(mGS)) / (sH*uH);

  }

  // Mass measure, (m^2)^(d-2).
  double tmpExp = eDdU - 2;
  eDsigma0 *= pow(mGS, tmpExp);

  // Constants.
  eDsigma0 *= eDconstantTerm;

}

//--------------------------------------------------------------------------

double Sigma2qg2LEDUnparticleq::sigmaHat() {

  // Mass spactrum weighting.
  double sigma = eDsigma0 /runBW3;

  // SM couplings...
  if (eDgraviton) {
    sigma *= 16 * M_PI * alpS / 96;
  } else if (eDspin == 1) {
    sigma *= - 4 * M_PI * alpS / 3;
  } else if (eDspin == 0) {
    sigma *= - 2 * M_PI * alpS / 3;
  }

  // Truncate sH region or use form factor.
  // Form factor uses either pythia8 renormScale2
  // or E_jet in cms.
  if (eDcutoff == 1) {
    if (sH > pow2(eDLambdaU) ) { sigma *= pow(eDLambdaU,4)/pow2(sH); }
  } else if ( (eDgraviton && (eDspin == 2))
           && ((eDcutoff == 2) || (eDcutoff == 3)) ) {
    double tmPmu = sqrt(Q2RenSave);
    if (eDcutoff == 3) tmPmu = (sH + s4 - s3) / (2 * mH);
    double tmPformfact = tmPmu / (eDtff * eDLambdaU);
    double tmPexp = double(eDnGrav) + 2;
    sigma *= 1 / (1 + pow(tmPformfact, tmPexp));
  }

  return sigma;
}

//--------------------------------------------------------------------------

void Sigma2qg2LEDUnparticleq::setIdColAcol() {

  // Flavour set up for q g -> G* q.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, eDidG, idq);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21);

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2LEDUnparticleg class.
// Cross section for q qbar -> U/G g (real graviton emission in
// large extra dimensions or unparticle emission).

//--------------------------------------------------------------------------

void Sigma2qqbar2LEDUnparticleg::initProc() {

  // Init model parameters.
  eDidG    = 5000039;
  if (eDgraviton) {
    eDspin     = (settingsPtr->flag("ExtraDimensionsLED:GravScalar")) ? 0 : 2;
    eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
    eDdU       = 0.5 * eDnGrav + 1;
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsLED:MD");
    eDlambda   = 1;
    eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
    eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");
    eDgf       = settingsPtr->parm("ExtraDimensionsLED:g");
    eDcf       = settingsPtr->parm("ExtraDimensionsLED:c");
  } else {
    eDspin     = settingsPtr->mode("ExtraDimensionsUnpart:spinU");
    eDdU       = settingsPtr->parm("ExtraDimensionsUnpart:dU");
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsUnpart:LambdaU");
    eDlambda   = settingsPtr->parm("ExtraDimensionsUnpart:lambda");
    eDcutoff   = settingsPtr->mode("ExtraDimensionsUnpart:CutOffMode");
  }

  // The A(dU) or S'(n) value.
  double tmpAdU = 0;
  if (eDgraviton) {
    tmpAdU  = 2 * M_PI * sqrt( pow(M_PI, double(eDnGrav)) )
            / GammaReal(0.5 * eDnGrav);
    // Scalar graviton
    if (eDspin == 0) {
      tmpAdU *= 2. * sqrt( pow(2., double(eDnGrav)) );
      eDcf   *= 4. * eDcf / pow2(eDLambdaU);
      double tmpExp = 2. * double(eDnGrav) / (double(eDnGrav) + 2.);
      eDgf   *= eDgf / pow(2. * M_PI, tmpExp);
    }
  } else {
    tmpAdU = 16 * pow2(M_PI) * sqrt(M_PI) / pow(2. * M_PI, 2. * eDdU)
      * GammaReal(eDdU + 0.5) / (GammaReal(eDdU - 1.) * GammaReal(2. * eDdU));
  }

  // Cross section related constants
  // and ME dependent powers of lambda / LambdaU.
  double tmpExp   = eDdU - 2;
  double tmpLS    = pow2(eDLambdaU);
  eDconstantTerm = tmpAdU / (2 * 16 * pow2(M_PI) * tmpLS * pow(tmpLS,tmpExp));
  if (eDgraviton && (eDspin == 2)) {
    eDconstantTerm /= tmpLS;
  } else if (eDspin == 1) {
    eDconstantTerm *= pow2(eDlambda);
  } else if (eDspin == 0) {
    eDconstantTerm *= pow2(eDlambda);
  } else {
    eDconstantTerm = 0;
    infoPtr->errorMsg("Error in Sigma2qqbar2LEDUnparticleg::initProc: "
                      "Incorrect spin value (turn process off)!");
  }

}

//--------------------------------------------------------------------------

void Sigma2qqbar2LEDUnparticleg::sigmaKin() {

  // Set graviton mass.
  mG        = m3;
  mGS       = mG*mG;

  // Set mandelstam variables and ME expressions.
  if (eDgraviton) {

    double A0 = 1/sH;
    // Scalar graviton
    if (eDspin == 0) {
      A0 /= sH;
      double tmpTerm1 = uH + tH;
      double T0 = (2. * mGS * sH + pow2(tmpTerm1)) / (uH * tH);
      double T1 = (tH2 + uH2) / sH;
      eDsigma0 = A0 * (eDgf * T0 + eDcf * T1);
    } else {
      double xH = tH/sH;
      double yH = mGS/sH;
      double xHS = pow2(xH);
      double yHS = pow2(yH);
      double xHC = pow(xH,3);
      double yHC = pow(yH,3);

      double T0 = 1/(xH*(yH-1-xH));
      double T1 = -4*xH*(1 + xH)*(1 + 2*xH + 2*xHS);
      double T2 = yH*(1 + 6*xH + 18*xHS + 16*xHC);
      double T3 = -6*yHS*xH*(1+2*xH);
      double T4 = yHC*(1 + 4*xH);

      eDsigma0 = A0 * T0 *( T1 + T2 + T3 + T4 );
    }

  } else if (eDspin == 1) {

    double A0  = 1/pow2(sH);
    double tmpTerm1 = tH - mGS;
    double tmpTerm2 = uH - mGS;

    eDsigma0 = A0 * (pow2(tmpTerm1) + pow2(tmpTerm2)) / (tH * uH);

  } else if (eDspin == 0) {

    double A0  = 1/pow2(sH);

    eDsigma0 = A0 * (pow2(sH) - pow2(mGS)) / (tH * uH);

  }

  // Mass measure, (m^2)^(d-2).
  double tmpExp = eDdU - 2;
  eDsigma0 *= pow(mGS, tmpExp);

  // Constants.
  eDsigma0 *= eDconstantTerm;

}

//--------------------------------------------------------------------------

double Sigma2qqbar2LEDUnparticleg::sigmaHat() {

  // Mass spactrum weighting.
  double sigma = eDsigma0 /runBW3;

  // SM couplings...
  if (eDgraviton) {
    sigma *= 16 * M_PI * alpS / 36;
  } else if (eDspin == 1) {
    sigma *= 4 * M_PI * 8 * alpS / 9;
  } else if (eDspin == 0) {
    sigma *= 4 * M_PI * 4 * alpS / 9;
  }

  // Truncate sH region or use form factor.
  // Form factor uses either pythia8 renormScale2
  // or E_jet in cms.
  if (eDcutoff == 1) {
    if (sH > pow2(eDLambdaU) ) { sigma *= pow(eDLambdaU,4)/pow2(sH); }
  } else if ( (eDgraviton && (eDspin == 2))
           && ((eDcutoff == 2) || (eDcutoff == 3)) ) {
    double tmPmu = sqrt(Q2RenSave);
    if (eDcutoff == 3) tmPmu = (sH + s4 - s3) / (2 * mH);
    double tmPformfact = tmPmu / (eDtff * eDLambdaU);
    double tmPexp = double(eDnGrav) + 2;
    sigma *= 1 / (1 + pow(tmPformfact, tmPexp));
  }

  return sigma;
}

//--------------------------------------------------------------------------

void Sigma2qqbar2LEDUnparticleg::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, eDidG, 21);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2ffbar2LEDUnparticleZ class.
// Cross section for f fbar -> U/G Z (real LED graviton or unparticle
// emission).

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// FIXRATIO:
// Ratio between the two possible coupling constants of the spin-2 ME.
// A value different from one give rise to an IR divergence which makes
// the event generation very slow, so this values is fixed to 1 until
// investigated further.
const double Sigma2ffbar2LEDUnparticleZ::FIXRATIO = 1.;

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDUnparticleZ::initProc() {

  // Init model parameters.
  eDidG        = 5000039;
  if (eDgraviton) {
    eDspin     = 2;
    eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
    eDdU       = 0.5 * eDnGrav + 1;
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsLED:MD");
    eDlambda   = 1;
    eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
    eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");
  } else {
    eDspin     = settingsPtr->mode("ExtraDimensionsUnpart:spinU");
    eDdU       = settingsPtr->parm("ExtraDimensionsUnpart:dU");
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsUnpart:LambdaU");
    eDlambda   = settingsPtr->parm("ExtraDimensionsUnpart:lambda");
    eDratio    = FIXRATIO;
    //         = settingsPtr->parm("ExtraDimensionsUnpart:ratio");
    eDcutoff   = settingsPtr->mode("ExtraDimensionsUnpart:CutOffMode");
  }

  // Store Z0 mass and width for propagator.
  mZ        = particleDataPtr->m0(23);
  widZ      = particleDataPtr->mWidth(23);
  mZS       = mZ*mZ;
  mwZS      = pow2(mZ * widZ);

  // Init spin-2 parameters
  if ( eDspin != 2 ){
    eDgraviton = false;
    eDlambdaPrime = 0;
  } else if (eDgraviton) {
    eDlambda = 1;
    eDratio = 1;
    eDlambdaPrime = eDlambda;
  } else {
    eDlambdaPrime = eDratio * eDlambda;
  }

  // The A(dU) or S'(n) value
  double tmpAdU = 16 * pow2(M_PI) * sqrt(M_PI) / pow(2. * M_PI, 2. * eDdU)
    * GammaReal(eDdU + 0.5) / (GammaReal(eDdU - 1.) * GammaReal(2. * eDdU));

  if (eDgraviton) {
    tmpAdU  = 2 * M_PI * sqrt( pow(M_PI, double(eDnGrav)) )
            / GammaReal(0.5 * eDnGrav);
  }

  // Standard 2 to 2 cross section related constants
  double tmpTerm1 = 1/(2 * 16 * pow2(M_PI));
  double tmpLS    = pow2(eDLambdaU);

  // Spin dependent constants from ME.
  double tmpTerm2 = 0;
  if ( eDspin == 0 ) {
    tmpTerm2 = 2 * pow2(eDlambda);
  } else if (eDspin == 1) {
    tmpTerm2 = 4 * pow2(eDlambda);
  } else if (eDspin == 2) {
    tmpTerm2 = pow2(eDlambda)/(4 * 3 * tmpLS);
  }

  // Unparticle phase space related
  double tmpExp2 = eDdU - 2;
  double tmpTerm3 = tmpAdU / (tmpLS * pow(tmpLS, tmpExp2));

  // All in total
  eDconstantTerm = tmpTerm1 * tmpTerm2 * tmpTerm3;

  // Secondary open width fraction.
  openFrac = particleDataPtr->resOpenFrac(23);

}

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDUnparticleZ::sigmaKin() {

  // Set graviton mass and some powers of mandelstam variables
  mU        = m3;
  mUS       = mU*mU;

  sHS = pow2(sH);
  tHS = pow2(tH);
  uHS = pow2(uH);
  tHC = pow(tH,3);
  uHC = pow(uH,3);
  tHQ = pow(tH,4);
  uHQ = pow(uH,4);
  tHuH = tH+uH;

  // Evaluate (m**2, t, u) part of differential cross section.
  // Extra 1/sHS comes from standard 2 to 2 cross section
  // phase space factors.

  if ( eDspin == 0 ) {

    double A0 = 1/sHS;
    double T1 = - sH/tH - sH/uH;
    double T2 = - (1 - mZS/tH)*(1 - mUS/tH);
    double T3 = - (1 - mZS/uH)*(1 - mUS/uH);
    double T4 = 2*(1 - mUS/tH)*(1 - mUS/uH);

    eDsigma0 = A0 * ( T1 + T2 + T3 + T4);

  } else if ( eDspin == 1 ) {

    double A0 = 1/sHS;
    double T1 = 0.5 * (tH/uH + uH/tH);
    double T2 =  pow2(mZS + mUS)/(tH * uH);
    double T3 = - 0.5 * mUS * (mZS/tHS + mZS/uHS) ;
    double T4 = - (mZS+mUS)*(1/tH + 1/uH);

    eDsigma0 = A0 * ( T1 + T2 + T3 + T4 );

  } else if ( eDspin == 2 ) {

    double A0   = 1 / ( sHS * uHS * tHS * pow2(sH-mZS) );
    double F0 = 2*tHS*uHS*( 16*pow(mZS,3) +  mUS*(7*tHS + 12*tH*uH + 7*uHS)
              - 3*(3*tHC + 11*tHS*uH + 11*tH*uHS + 3*uHC)
              + 6*pow(mZS,2)*(7*mUS - 2*tHuH) + mZS*(14*pow(mUS,2)
              - 15*tHS - 44*tH*uH - 15*uHS + 2*mUS*tHuH) );
    double F2 = 2*tHS*uHS*tHuH*( -8*pow(mZS,2)*tHuH
              + 4*mZS*(tHS + 3*tH*uH + uHS)
              + 3*(tHC + 5*tHS*uH + 5*tH*uHS + uHC) );
    double F4 = -2*tHS*uHS*pow(tHuH,3)*(tHS + uHS - mZS*tHuH);

    double G0 = 4*tH*uH*( 6*pow(mZS,3)*(mUS - tH - uH)*tHuH
              + pow(mZS,2)*( 9*tHC + 7*tHS*uH + 7*tH*uHS + 9*uHC
              + 15*pow2(mUS)*tHuH - 2*mUS*(12*tHS + 19*tH*uH + 12*uHS) )
              + tH*uH*( 6*pow(mUS,3) - 9*pow(mUS,2)*tHuH - mUS*(tHS
              + 12*tH*uH + uHS) + 6*(tHC + 6*tHS*uH + 6*tH*uHS + uHC) )
              + mZS*(-3*tHQ + 25*tHC*uH + 58*tHS*uHS + 25*tH*uHC
              - 3*uHQ + 6*pow(mUS,3)*tHuH
              - pow(mUS,2)*(15*tHS + 2*tH*uH + 15*uHS) + 2*mUS*(6*tHC
              - 11*tHS*uH - 11*tH*uHS + 6*uHC)) );
    double G2 = -4*tHS*uHS*tHuH*( -10*pow2(mZS)*tHuH + 2*mZS*(3*tHS
              + 7*tH*uH + 3*uHS) + 3*(tHC + 5*tHS*uH + 5*tH*uHS + uHC) );
    double G4 = -2*F4;

    double H0 = 24*pow(mZS,3)*tH*uH*pow2(-mUS + tHuH)
              - 6*pow(mZS,2)*tH*uH*( -9*pow(mUS,3) + 24*pow(mUS,2)*tHuH
              - mUS*(21*tHS + 38*tH*uH + 21*uHS)
              + 2*(3*tHC + 5*tHS*uH + 5*tH*uHS + 3*uHC) )
              - mZS*( 3*pow(mUS,4)*(tHS - 12*tH*uH + uHS)
              - 2*tH*uH*pow2(tHuH)*(6*tHS - 29*tH*uH + 6*uHS)
              - 6*pow(mUS,3)*(tHC - 16*tHS*uH - 16*tH*uHS + uHC)
              + 54*mUS*tH*uH*(tHC + tHS*uH + tH*uHS + uHC)
              + pow2(mUS)*(3*tHQ - 102*tHC*uH - 166*tHS*uHS
              - 102*tH*uHC + 3*uHQ) )
              + tH*uH*( 6*pow(mUS,5) - 18*pow(mUS,4)*tHuH
              - 12*pow(mUS,2)*pow(tHuH,3)
              + 3*pow(mUS,3)*(7*tHS + 12*tH*uH + 7*uHS)
              - 18*tH*uH*(tHC + 5*tHS*uH + 5*tH*uHS + uHC)
              + mUS*(3*tHQ + 32*tHC*uH + 78*tHS*uHS + 32*tH*uHC + 3*uHQ) );
    double H2 = 2*tHS*uHS*pow2(tHuH)*( -12*pow2(mZS) + 8*mZS*tHuH
              + 3*(tHS + 4*tH*uH + uHS) );
    double H4 = F4;

    eDsigma0 = A0*( F0 + 1/mUS*F2 + 1/pow2(mUS)*F4
             + eDratio*(G0 + 1/mUS*G2 + 1/pow2(mUS)*G4)
             + pow2(eDratio)*(H0 + 1/mUS*H2 + 1/pow2(mUS)*H4) );

  } else {

    eDsigma0 = 0;

  }

}

//--------------------------------------------------------------------------

double Sigma2ffbar2LEDUnparticleZ::sigmaHat() {

  // Electroweak couplings.
  int idAbs    = abs(id1);
  // Note: 1/2 * (g_L^2 + g_R^2) = (g_v^2 + g_a^2)
  double facEWS  = 4 * M_PI * alpEM
                   / (couplingsPtr->sin2thetaW() * couplingsPtr->cos2thetaW())
                   * ( 0.25 * 0.25 * couplingsPtr->vf2af2(idAbs) );

  // Mass Spectrum, (m^2)^(d-2)
  double tmpExp = eDdU - 2;
  double facSpect = pow(mUS, tmpExp);

  // Total cross section
  double sigma = eDconstantTerm * facEWS * facSpect * eDsigma0;

  // Secondary width for Z0
  sigma       *= openFrac;

  // If f fbar are quarks (1/N_c)
  if (idAbs < 9) sigma /= 3.;

  // Related to mass spactrum weighting.
  sigma /= runBW3;

  // Truncate sH region or use form factor.
  // Form factor uses either pythia8 renormScale2
  // or E_jet in cms.
  if (eDcutoff == 1) {
    if (sH > pow2(eDLambdaU) ) { sigma *= pow(eDLambdaU,4)/pow2(sH); }
  } else if (eDgraviton && ((eDcutoff == 2) || (eDcutoff == 3))) {
    double tmPmu = sqrt(Q2RenSave);
    if (eDcutoff == 3) tmPmu = (sH + s4 - s3) / (2 * mH);
    double tmPformfact = tmPmu / (eDtff * eDLambdaU);
    double tmPexp = double(eDnGrav) + 2;
    sigma *= 1 / (1 + pow(tmPformfact, tmPexp));
  }

  return sigma;

}

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDUnparticleZ::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, eDidG, 23);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2ffbar2LEDUnparticlegamma class.
// Cross section for f fbar -> U/G gamma (real LED graviton or unparticle
// emission).

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// FIXRATIO:
// Ratio between the two possible coupling constants of the spin-2 ME.
// A value different from one give rise to an IR divergence which makes
// the event generation very slow, so this values is fixed to 1 until
// investigated further.
const double Sigma2ffbar2LEDUnparticlegamma::FIXRATIO = 1.;

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDUnparticlegamma::initProc() {

  // WARNING: Keep in mind that this class uses the photon limit
  //          of the Z+G/U ME code. This might give rise to some
  //          confusing things, e.g. mZ = particleDataPtr->m0(22);

  // Init model parameters.
  eDidG        = 5000039;
  if (eDgraviton) {
    eDspin     = 2;
    eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
    eDdU       = 0.5 * eDnGrav + 1;
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsLED:MD");
    eDlambda   = 1;
    eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
    eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");
  } else {
    eDspin     = settingsPtr->mode("ExtraDimensionsUnpart:spinU");
    eDdU       = settingsPtr->parm("ExtraDimensionsUnpart:dU");
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsUnpart:LambdaU");
    eDlambda   = settingsPtr->parm("ExtraDimensionsUnpart:lambda");
    eDratio    = FIXRATIO;
    //         = settingsPtr->parm("ExtraDimensionsUnpart:ratio");
    eDcutoff   = settingsPtr->mode("ExtraDimensionsUnpart:CutOffMode");
  }

  // Store Z0 mass.
  mZ        = particleDataPtr->m0(22);
  mZS       = mZ*mZ;

  // Init spin-2 parameters
  if ( eDspin != 2 ){
    eDgraviton = false;
    eDlambdaPrime = 0;
  } else if (eDgraviton) {
    eDlambda = 1;
    eDratio = 1;
    eDlambdaPrime = eDlambda;
  } else {
    eDlambdaPrime = eDratio * eDlambda;
  }

  // The A(dU) or S'(n) value
  double tmpAdU = 16 * pow2(M_PI) * sqrt(M_PI) / pow(2. * M_PI, 2. * eDdU)
    * GammaReal(eDdU + 0.5) / (GammaReal(eDdU - 1.) * GammaReal(2. * eDdU));

  if (eDgraviton) {
    tmpAdU  = 2 * M_PI * sqrt( pow(M_PI, double(eDnGrav)) )
            / GammaReal(0.5 * eDnGrav);
  }

  // Standard 2 to 2 cross section related constants
  double tmpTerm1 = 1/(2 * 16 * pow2(M_PI));
  double tmpLS    = pow2(eDLambdaU);

  // Spin dependent constants from ME.
  double tmpTerm2 = 0;
  if ( eDspin == 0 ) {
    tmpTerm2 = 2 * pow2(eDlambda);
  } else if (eDspin == 1) {
    tmpTerm2 = 4 * pow2(eDlambda);
  } else if (eDspin == 2) {
    tmpTerm2 = pow2(eDlambda)/(4 * 3 * tmpLS);
  }

  // Unparticle phase space related
  double tmpExp2 = eDdU - 2;
  double tmpTerm3 = tmpAdU / (tmpLS * pow(tmpLS, tmpExp2));

  // All in total
  eDconstantTerm = tmpTerm1 * tmpTerm2 * tmpTerm3;

}

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDUnparticlegamma::sigmaKin() {

  // Set graviton mass and some powers of mandelstam variables
  mU        = m3;
  mUS       = mU*mU;

  sHS = pow2(sH);
  tHS = pow2(tH);
  uHS = pow2(uH);
  tHC = pow(tH,3);
  uHC = pow(uH,3);
  tHQ = pow(tH,4);
  uHQ = pow(uH,4);
  tHuH = tH+uH;

  // Evaluate (m**2, t, u) part of differential cross section.
  // Extra 1/sHS comes from standard 2 to 2 cross section
  // phase space factors.

  if ( eDspin == 0 ) {

    double A0 = 1/sHS;
    double T1 = - sH/tH - sH/uH;
    double T2 = - (1 - mZS/tH)*(1 - mUS/tH);
    double T3 = - (1 - mZS/uH)*(1 - mUS/uH);
    double T4 = 2*(1 - mUS/tH)*(1 - mUS/uH);

    eDsigma0 = A0 * ( T1 + T2 + T3 + T4);

  } else if ( eDspin == 1 ) {

    double A0 = 1/sHS;
    double T1 = 0.5 * (tH/uH + uH/tH);
    double T2 =  pow2(mZS + mUS)/(tH * uH);
    double T3 = - 0.5 * mUS * (mZS/tHS + mZS/uHS) ;
    double T4 = - (mZS+mUS)*(1/tH + 1/uH);

    eDsigma0 = A0 * ( T1 + T2 + T3 + T4 );

  } else if ( eDspin == 2 ) {

    double A0 = 1 / ( sHS * uHS * tHS * pow2(sH-mZS) );
    double F0 = 2*tHS*uHS*( 16*pow(mZS,3) +  mUS*(7*tHS + 12*tH*uH + 7*uHS)
              - 3*(3*tHC + 11*tHS*uH + 11*tH*uHS + 3*uHC)
              + 6*pow(mZS,2)*(7*mUS - 2*tHuH) + mZS*(14*pow(mUS,2)
              - 15*tHS - 44*tH*uH - 15*uHS + 2*mUS*tHuH) );
    double F2 = 2*tHS*uHS*tHuH*( -8*pow(mZS,2)*tHuH
              + 4*mZS*(tHS + 3*tH*uH + uHS)
              + 3*(tHC + 5*tHS*uH + 5*tH*uHS + uHC) );
    double F4 = -2*tHS*uHS*pow(tHuH,3)*(tHS + uHS - mZS*tHuH);

    double G0 = 4*tH*uH*( 6*pow(mZS,3)*(mUS - tH - uH)*tHuH
              + pow(mZS,2)*( 9*tHC + 7*tHS*uH + 7*tH*uHS + 9*uHC
              + 15*pow2(mUS)*tHuH - 2*mUS*(12*tHS + 19*tH*uH + 12*uHS) )
              + tH*uH*( 6*pow(mUS,3) - 9*pow(mUS,2)*tHuH
              - mUS*(tHS + 12*tH*uH + uHS)
              + 6*(tHC + 6*tHS*uH + 6*tH*uHS + uHC) )
              + mZS*(-3*tHQ + 25*tHC*uH + 58*tHS*uHS + 25*tH*uHC
              - 3*uHQ + 6*pow(mUS,3)*tHuH
              - pow(mUS,2)*(15*tHS + 2*tH*uH + 15*uHS)
              + 2*mUS*(6*tHC - 11*tHS*uH - 11*tH*uHS + 6*uHC)) );
    double G2 = -4*tHS*uHS*tHuH*( -10*pow2(mZS)*tHuH
              + 2*mZS*(3*tHS + 7*tH*uH + 3*uHS)
              + 3*(tHC + 5*tHS*uH + 5*tH*uHS + uHC) );
    double G4 = -2*F4;

    double H0 = 24*pow(mZS,3)*tH*uH*pow2(-mUS + tHuH)
              - 6*pow(mZS,2)*tH*uH*( -9*pow(mUS,3) + 24*pow(mUS,2)*tHuH
              - mUS*(21*tHS + 38*tH*uH + 21*uHS)
              + 2*(3*tHC + 5*tHS*uH + 5*tH*uHS + 3*uHC) )
              - mZS*( 3*pow(mUS,4)*(tHS - 12*tH*uH + uHS)
              - 2*tH*uH*pow2(tHuH)*(6*tHS - 29*tH*uH + 6*uHS)
              - 6*pow(mUS,3)*(tHC - 16*tHS*uH - 16*tH*uHS + uHC)
              + 54*mUS*tH*uH*(tHC + tHS*uH + tH*uHS + uHC)
              + pow2(mUS)*(3*tHQ - 102*tHC*uH - 166*tHS*uHS
              - 102*tH*uHC + 3*uHQ) )
              + tH*uH*( 6*pow(mUS,5) - 18*pow(mUS,4)*tHuH
              - 12*pow(mUS,2)*pow(tHuH,3)
              + 3*pow(mUS,3)*(7*tHS + 12*tH*uH + 7*uHS)
              - 18*tH*uH*(tHC + 5*tHS*uH + 5*tH*uHS + uHC)
              + mUS*(3*tHQ + 32*tHC*uH + 78*tHS*uHS + 32*tH*uHC + 3*uHQ) );
    double H2 = 2*tHS*uHS*pow2(tHuH)*( -12*pow2(mZS) + 8*mZS*tHuH
              + 3*(tHS + 4*tH*uH + uHS) );
    double H4 = F4;

    eDsigma0 = A0*( F0 + 1/mUS*F2 + 1/pow2(mUS)*F4
             + eDratio*(G0 + 1/mUS*G2 + 1/pow2(mUS)*G4)
             + pow2(eDratio)*(H0 + 1/mUS*H2 + 1/pow2(mUS)*H4) );

  } else {

    eDsigma0 = 0;

  }

}

//--------------------------------------------------------------------------

double Sigma2ffbar2LEDUnparticlegamma::sigmaHat() {

  // Electroweak couplings..
  int idAbs    = abs(id1);
  double facEWS = 4 * M_PI * alpEM * couplingsPtr->ef2(idAbs);

  // Mass Spectrum, (m^2)^(d-2)
  double tmpExp = eDdU - 2;
  double facSpect = pow(mUS, tmpExp);

  // Total cross section
  double sigma = eDconstantTerm * facEWS * facSpect * eDsigma0;

  // If f fbar are quarks
  if (idAbs < 9) sigma /= 3.;

  // Related to mass spactrum weighting.
  sigma /= runBW3;

  // Truncate sH region or use form factor.
  // Form factor uses either pythia8 renormScale2
  // or E_jet in cms.
  if (eDcutoff == 1) {
    if (sH > pow2(eDLambdaU) ) { sigma *= pow(eDLambdaU,4)/pow2(sH); }
  } else if (eDgraviton && ((eDcutoff == 2) || (eDcutoff == 3))) {
    double tmPmu = sqrt(Q2RenSave);
    if (eDcutoff == 3) tmPmu = (sH + s4 - s3) / (2 * mH);
    double tmPformfact = tmPmu / (eDtff * eDLambdaU);
    double tmPexp = double(eDnGrav) + 2;
    sigma *= 1 / (1 + pow(tmPformfact, tmPexp));
  }

  return sigma;

}

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDUnparticlegamma::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, eDidG, 22);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2ffbar2LEDgammagamma class.
// Cross section for f fbar -> (LED G*/U*) -> gamma gamma
// (virtual graviton/unparticle exchange).

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDgammagamma::initProc() {

  // Init model parameters.
  if (eDgraviton) {
    eDspin     = 2;
    eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
    eDdU       = 2;
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsLED:LambdaT");
    eDlambda   = 1;
    eDnegInt   = settingsPtr->mode("ExtraDimensionsLED:NegInt");
    eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
    eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");
  } else {
    eDspin     = settingsPtr->mode("ExtraDimensionsUnpart:spinU");
    eDdU       = settingsPtr->parm("ExtraDimensionsUnpart:dU");
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsUnpart:LambdaU");
    eDlambda   = settingsPtr->parm("ExtraDimensionsUnpart:lambda");
    eDnegInt   = 0;
  }

  // Model dependent constants.
  if (eDgraviton) {
    eDlambda2chi = 4*M_PI;
    if (eDnegInt == 1) eDlambda2chi *= -1.;
  } else {
    double tmPAdU = 16 * pow2(M_PI) * sqrt(M_PI) / pow(2. * M_PI, 2. * eDdU)
      * GammaReal(eDdU + 0.5) / (GammaReal(eDdU - 1.) * GammaReal(2. * eDdU));
    double tmPdUpi = eDdU * M_PI;
    eDlambda2chi = pow2(eDlambda) * tmPAdU / (2 * sin(tmPdUpi));
  }

  // Model parameter check (if not applicable, sigma = 0).
  // Note: SM contribution still generated.
  if ( !(eDspin==0 || eDspin==2) ) {
    eDlambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2ffbar2LEDgammagamma::initProc: "
                      "Incorrect spin value (turn process off)!");
  } else if ( !eDgraviton && (eDdU >= 2)) {
    eDlambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2ffbar2LEDgammagamma::initProc: "
                      "This process requires dU < 2 (turn process off)!");
  }

}

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDgammagamma::sigmaKin() {

  // Mandelstam variables.
  double sHS = pow2(sH);
  double sHQ = pow(sH, 4);
  double tHS = pow2(tH);
  double uHS = pow2(uH);

  // Form factor.
  double tmPeffLambdaU = eDLambdaU;
  if (eDgraviton && ((eDcutoff == 2) || (eDcutoff == 3))) {
    double tmPffterm = sqrt(Q2RenSave) / (eDtff * eDLambdaU);
    double tmPexp = double(eDnGrav) + 2;
    double tmPformfact = 1 + pow(tmPffterm, tmPexp);
    tmPeffLambdaU *= pow(tmPformfact,0.25);
  }

  // ME from spin-0 and spin-2 unparticles
  // including extra 1/sHS from 2-to-2 phase space.
  if (eDspin == 0) {
    double tmPsLambda2 = sH / pow2(tmPeffLambdaU);
    double tmPexp = 2 * eDdU - 1;
    eDterm1 = pow(tmPsLambda2,tmPexp);
    eDterm1 /= sHS;
  } else {
    eDterm1 = (uH / tH + tH / uH);
    eDterm1 /= sHS;
    double tmPsLambda2 = sH / pow2(tmPeffLambdaU);
    double tmPexp = eDdU;
    eDterm2 = pow(tmPsLambda2,tmPexp) * (uHS + tHS) / sHS;
    eDterm2 /= sHS;
    tmPexp = 2 * eDdU;
    eDterm3 = pow(tmPsLambda2,tmPexp) * tH * uH * (uHS + tHS) / sHQ;
    eDterm3 /= sHS;
  }

}

//--------------------------------------------------------------------------

double Sigma2ffbar2LEDgammagamma::sigmaHat() {

  // Incoming fermion flavor.
  int idAbs      = abs(id1);

  // Couplings and constants.
  // Note: ME already contain 1/2 for identical
  //       particles in the final state.
  double sigma = 0;
  if (eDspin == 0) {
    sigma = pow2(eDlambda2chi) * eDterm1 / 8;
  } else {
    double tmPe2Q2 = 4 * M_PI * alpEM * couplingsPtr->ef2(idAbs);
    double tmPdUpi = eDdU * M_PI;
    sigma = pow2(tmPe2Q2) * eDterm1
          - tmPe2Q2 * eDlambda2chi * cos(tmPdUpi) * eDterm2
          + pow2(eDlambda2chi) * eDterm3 / 4;
  }

  // dsigma/dt, 2-to-2 phase space factors.
  sigma /= 16 * M_PI;

  // If f fbar are quarks.
  if (idAbs < 9) sigma /= 3.;

  return sigma;
}

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDgammagamma::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 22, 22);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2gg2LEDgammagamma class.
// Cross section for g g -> (LED G*/U*) -> gamma gamma
// (virtual graviton/unparticle exchange).

//--------------------------------------------------------------------------

void Sigma2gg2LEDgammagamma::initProc() {

  // Init model parameters.
  if (eDgraviton) {
    eDspin     = 2;
    eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
    eDdU       = 2;
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsLED:LambdaT");
    eDlambda   = 1;
    eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
    eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");
  } else {
    eDspin     = settingsPtr->mode("ExtraDimensionsUnpart:spinU");
    eDdU       = settingsPtr->parm("ExtraDimensionsUnpart:dU");
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsUnpart:LambdaU");
    eDlambda   = settingsPtr->parm("ExtraDimensionsUnpart:lambda");
  }

  // Model dependent constants.
  if (eDgraviton) {
    eDlambda2chi = 4 * M_PI;

  } else {
    double tmPAdU = 16 * pow2(M_PI) * sqrt(M_PI) / pow(2. * M_PI, 2. * eDdU)
      * GammaReal(eDdU + 0.5) / (GammaReal(eDdU - 1.) * GammaReal(2. * eDdU));
    double tmPdUpi = eDdU * M_PI;
    eDlambda2chi = pow2(eDlambda) * tmPAdU / (2 * sin(tmPdUpi));
  }

  // Model parameter check (if not applicable, sigma = 0).
  if ( !(eDspin==0 || eDspin==2) ) {
    eDlambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2gg2LEDgammagamma::initProc: "
                      "Incorrect spin value (turn process off)!");
  } else if ( !eDgraviton && (eDdU >= 2)) {
    eDlambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2gg2LEDgammagamma::initProc: "
                      "This process requires dU < 2 (turn process off)!");
  }

}

//--------------------------------------------------------------------------

void Sigma2gg2LEDgammagamma::sigmaKin() {

  // Mandelstam variables.
  double sHS = pow2(sH);
  double sHQ = pow(sH, 4);
  double tHQ = pow(tH, 4);
  double uHQ = pow(uH, 4);

  // Form factor.
  double tmPeffLambdaU = eDLambdaU;
  if (eDgraviton && ((eDcutoff == 2) || (eDcutoff == 3))) {
    double tmPffterm = sqrt(Q2RenSave) / (eDtff * eDLambdaU);
    double tmPexp = double(eDnGrav) + 2;
    double tmPformfact = 1 + pow(tmPffterm, tmPexp);
    tmPeffLambdaU *= pow(tmPformfact,0.25);
  }

  // ME from spin-0 and spin-2 unparticles.
  if (eDspin == 0) {
    double tmPsLambda2 = sH / pow2(tmPeffLambdaU);
    double tmPexp = 2 * eDdU;
    eDsigma0 = pow(tmPsLambda2,tmPexp);
  } else {
    double tmPsLambda2 = sH / pow2(tmPeffLambdaU);
    double tmPexp = 2 * eDdU;
    eDsigma0 = pow(tmPsLambda2,tmPexp) * (uHQ + tHQ) / sHQ;
  }

  // extra 1/sHS from 2-to-2 phase space.
  eDsigma0 /= sHS;

}

//--------------------------------------------------------------------------

double Sigma2gg2LEDgammagamma::sigmaHat() {

  // Couplings and constants.
  // Note: ME already contain 1/2 for identical
  //       particles in the final state.
  double sigma = eDsigma0;
  if (eDspin == 0) {
    sigma *= pow2(eDlambda2chi) / 256;
  } else {
    sigma *= pow2(eDlambda2chi) / 32;
  }

  // dsigma/dt, 2-to-2 phase space factors.
  sigma /= 16 * M_PI;

  return sigma;
}

//--------------------------------------------------------------------------

void Sigma2gg2LEDgammagamma::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, 22, 22);

  // Colour flow topologies.
  setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

}

//==========================================================================

// Sigma2ffbar2LEDllbar class.
// Cross section for f fbar -> (LED G*/U*) -> l lbar
// (virtual graviton/unparticle exchange).
// Does not include t-channel contributions relevant for e^+e^- to e^+e^-

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDllbar::initProc() {

  // Init model parameters.
  if (eDgraviton) {
    eDspin     = 2;
    eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
    eDdU       = 2;
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsLED:LambdaT");
    eDlambda   = 1;
    eDnegInt   = settingsPtr->mode("ExtraDimensionsLED:NegInt");
    eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
    eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");
  } else {
    eDspin     = settingsPtr->mode("ExtraDimensionsUnpart:spinU");
    eDdU       = settingsPtr->parm("ExtraDimensionsUnpart:dU");
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsUnpart:LambdaU");
    eDlambda   = settingsPtr->parm("ExtraDimensionsUnpart:lambda");
    eDnxx      = settingsPtr->mode("ExtraDimensionsUnpart:gXX");
    eDnxy      = settingsPtr->mode("ExtraDimensionsUnpart:gXY");
    eDnegInt   = 0;
  }

  eDmZ  = particleDataPtr->m0(23);
  eDmZS = eDmZ * eDmZ;
  eDGZ  = particleDataPtr->mWidth(23);
  eDGZS = eDGZ * eDGZ;

  // Model dependent constants.
  if (eDgraviton) {
    eDlambda2chi = 4*M_PI;
    if (eDnegInt == 1) eDlambda2chi *= -1.;
  } else {
    double tmPAdU = 16 * pow2(M_PI) * sqrt(M_PI) / pow(2. * M_PI, 2. * eDdU)
      * GammaReal(eDdU + 0.5) / (GammaReal(eDdU - 1.) * GammaReal(2. * eDdU));
    double tmPdUpi = eDdU * M_PI;
    eDlambda2chi = pow2(eDlambda) * tmPAdU / (2 * sin(tmPdUpi));
  }

  // Model parameter check (if not applicable, sigma = 0).
  // Note: SM contribution still generated.
  if ( !(eDspin==1 || eDspin==2) ) {
    eDlambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2ffbar2LEDllbar::initProc: "
                      "Incorrect spin value (turn process off)!");
  } else if ( !eDgraviton && (eDdU >= 2)) {
    eDlambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2ffbar2LEDllbar::initProc: "
                      "This process requires dU < 2 (turn process off)!");
  }

}

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDllbar::sigmaKin() {

  // Mandelstam variables.
  double tHS = pow2(tH);
  double uHS = pow2(uH);
  double tHC = pow(tH,3);
  double uHC = pow(uH,3);
  double tHQ = pow(tH,4);
  double uHQ = pow(uH,4);

  // Form factor.
  double tmPeffLambdaU = eDLambdaU;
  if (eDgraviton && ((eDcutoff == 2) || (eDcutoff == 3))) {
    double tmPffterm = sqrt(Q2RenSave) / (eDtff * eDLambdaU);
    double tmPexp = double(eDnGrav) + 2;
    double tmPformfact = 1 + pow(tmPffterm, tmPexp);
    tmPeffLambdaU *= pow(tmPformfact,0.25);
  }

  // ME from spin-1 and spin-2 unparticles
  eDdenomPropZ = pow2(sH - eDmZS) + eDmZS * eDGZS;
  eDrePropZ = (sH - eDmZS) / eDdenomPropZ;
  eDimPropZ = -eDmZ * eDGZ / eDdenomPropZ;
  eDrePropGamma = 1 / sH;
  if (eDspin == 1) {
    double tmPsLambda2 = sH / pow2(tmPeffLambdaU);
    double tmPexp = eDdU - 2;
    eDabsMeU  = eDlambda2chi * pow(tmPsLambda2,tmPexp)
              / pow2(tmPeffLambdaU);
  } else {
    double tmPsLambda2 = sH / pow2(tmPeffLambdaU);
    double tmPexp = eDdU - 2;
    double tmPA = -eDlambda2chi * pow(tmPsLambda2,tmPexp)
                 / (8 * pow(tmPeffLambdaU,4));
    eDabsAS = pow2(tmPA);
    eDreA   = tmPA * cos(M_PI * eDdU);
    eDreABW = tmPA * ((sH - eDmZS) * cos(M_PI * eDdU) + eDmZ * eDGZ
            * sin(M_PI * eDdU)) / eDdenomPropZ;
    eDpoly1 = tHQ + uHQ - 6*tHC*uH - 6*tH*uHC + 18*tHS*uHS;
    double tmPdiffUT = uH - tH;
    eDpoly2 = pow(tmPdiffUT,3);
    eDpoly3 = tHC - 3*tHS*uH - 3*tH*uHS + uHC;
  }

}

//--------------------------------------------------------------------------

double Sigma2ffbar2LEDllbar::sigmaHat() {

  // Incoming fermion flavor.
  int idAbs      = abs(id1);

  // Couplings and constants.
  // Qq = couplingsPtr->ef(idAbs), quark, i.e. id > 0.
  // Ql = couplingsPtr->ef(11), electron.
  double tmPe2QfQl = 4 * M_PI * alpEM * couplingsPtr->ef(idAbs)
                      * couplingsPtr->ef(11);
  double tmPgvq = 0.25 * couplingsPtr->vf(idAbs);
  double tmPgaq = 0.25 * couplingsPtr->af(idAbs);
  double tmPgLq = tmPgvq  + tmPgaq;
  double tmPgRq = tmPgvq  - tmPgaq;
  double tmPgvl = 0.25 * couplingsPtr->vf(11);
  double tmPgal = 0.25 * couplingsPtr->af(11);
  double tmPgLl = tmPgvl  + tmPgal;
  double tmPgRl = tmPgvl  - tmPgal;
  double tmPe2s2c2 = 4 * M_PI * alpEM
    / (couplingsPtr->sin2thetaW() * couplingsPtr->cos2thetaW());

  // LL, RR, LR, RL  couplings.
  vector<double> tmPcoupZ;
  tmPcoupZ.push_back(tmPe2s2c2 * tmPgLq * tmPgLl);
  tmPcoupZ.push_back(tmPe2s2c2 * tmPgRq * tmPgRl);
  tmPcoupZ.push_back(tmPe2s2c2 * tmPgRq * tmPgLl);
  tmPcoupZ.push_back(tmPe2s2c2 * tmPgLq * tmPgRl);
  vector<double> tmPcoupU;
  if (eDnxx == 1) {
    // LL
    tmPcoupU.push_back(-1);
    // RR
    tmPcoupU.push_back(-1);
  } else if (eDnxx == 2) {
    // LL
    tmPcoupU.push_back(0);
    // RR
    tmPcoupU.push_back(0);
  } else {
    // LL
    tmPcoupU.push_back(1);
    // RR
    tmPcoupU.push_back(1);
  }
  if (eDnxy == 1) {
    // RL
    tmPcoupU.push_back(-1);
    // LR
    tmPcoupU.push_back(-1);
  } else if (eDnxy == 2) {
    // RL
    tmPcoupU.push_back(0);
    // LR
    tmPcoupU.push_back(0);
  } else {
    // RL
    tmPcoupU.push_back(1);
    // LR
    tmPcoupU.push_back(1);
  }

  // Matrix elements
  double tmPMES = 0;
  if (eDspin == 1) {

    for (unsigned int i = 0; i<tmPcoupZ.size(); ++i) {
      double tmPMS = pow2(tmPcoupU[i] * eDabsMeU)
        + pow2(tmPe2QfQl * eDrePropGamma)
        + pow2(tmPcoupZ[i]) / eDdenomPropZ
        + 2 * cos(M_PI * eDdU) * tmPcoupU[i] * eDabsMeU
            * tmPe2QfQl * eDrePropGamma
        + 2 * cos(M_PI * eDdU) * tmPcoupU[i] * eDabsMeU
            * tmPcoupZ[i] * eDrePropZ
        + 2 * tmPe2QfQl * eDrePropGamma
            * tmPcoupZ[i] * eDrePropZ
        - 2 * sin(M_PI * eDdU) * tmPcoupU[i] * eDabsMeU
            * tmPcoupZ[i] * eDimPropZ;

      if (i<2) { tmPMES += 4 * pow2(uH) * tmPMS; }
      else if (i<4) { tmPMES += 4 * pow2(tH) * tmPMS; }
    }

  } else {

    for (unsigned int i = 0; i<tmPcoupZ.size(); ++i) {
      double tmPMS = pow2(tmPe2QfQl * eDrePropGamma)
        + pow2(tmPcoupZ[i]) / eDdenomPropZ
        + 2 * tmPe2QfQl * eDrePropGamma  * tmPcoupZ[i] * eDrePropZ;

      if (i<2) { tmPMES += 4 * pow2(uH) * tmPMS; }
      else if (i<4) { tmPMES += 4 * pow2(tH) * tmPMS; }
    }
    tmPMES += 8 * eDabsAS * eDpoly1;
    tmPMES += 16 * tmPe2QfQl * eDrePropGamma * eDreA * eDpoly2;
    tmPMES += 16 * tmPe2s2c2 * eDreABW * (tmPgaq * tmPgal * eDpoly3
                                          + tmPgvq * tmPgvl * eDpoly2);

  }

  // dsigma/dt, 2-to-2 phase space factors.
  double sigma = 0.25 * tmPMES;  // 0.25, is the spin average
  sigma /= 16 * M_PI * pow2(sH);

  // If f fbar are quarks.
  if (idAbs < 9) sigma /= 3.;

  // sigma(ffbar->llbar) = 3 * sigma(ffbar->eebar)
  sigma *= 3.;

  return sigma;
}

//--------------------------------------------------------------------------

void Sigma2ffbar2LEDllbar::setIdColAcol() {

  double tmPrand = rndmPtr->flat();
  // Flavours trivial.
  if (tmPrand < 0.33333333) {      setId( id1, id2, 11, -11); }
  else if (tmPrand < 0.66666667) { setId( id1, id2, 13, -13); }
  else {                            setId( id1, id2, 15, -15); }

  // tH defined between f and f': must swap tHat <-> uHat if id1 is fbar.
  swapTU = (id2 > 0);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2gg2LEDllbar class.
// Cross section for g g -> (LED G*/U*) -> l lbar
// (virtual graviton/unparticle exchange).

//--------------------------------------------------------------------------

void Sigma2gg2LEDllbar::initProc() {

  // Init model parameters.
  if (eDgraviton) {
    eDspin     = 2;
    eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
    eDdU       = 2;
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsLED:LambdaT");
    eDlambda   = 1;
    eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
    eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");
  } else {
    eDspin     = settingsPtr->mode("ExtraDimensionsUnpart:spinU");
    eDdU       = settingsPtr->parm("ExtraDimensionsUnpart:dU");
    eDLambdaU  = settingsPtr->parm("ExtraDimensionsUnpart:LambdaU");
    eDlambda   = settingsPtr->parm("ExtraDimensionsUnpart:lambda");
  }

  // Model dependent constants.
  if (eDgraviton) {
    eDlambda2chi = 4 * M_PI;

  } else {
    double tmPAdU = 16 * pow2(M_PI) * sqrt(M_PI) / pow(2. * M_PI, 2. * eDdU)
      * GammaReal(eDdU + 0.5) / (GammaReal(eDdU - 1.) * GammaReal(2. * eDdU));
    double tmPdUpi = eDdU * M_PI;
    eDlambda2chi = pow2(eDlambda) * tmPAdU / (2 * sin(tmPdUpi));
  }

  // Model parameter check (if not applicable, sigma = 0).
  if ( !(eDspin==2) ) {
    eDlambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2gg2LEDllbar::initProc: "
                      "Incorrect spin value (turn process off)!");
  } else if ( !eDgraviton && (eDdU >= 2)) {
    eDlambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2gg2LEDllbar::initProc: "
                      "This process requires dU < 2 (turn process off)!");
  }

}

//--------------------------------------------------------------------------

void Sigma2gg2LEDllbar::sigmaKin() {

  // Form factor.
  double tmPeffLambdaU = eDLambdaU;
  if (eDgraviton && ((eDcutoff == 2) || (eDcutoff == 3))) {
    double tmPffterm = sqrt(Q2RenSave) / (eDtff * eDLambdaU);
    double tmPexp = double(eDnGrav) + 2;
    double tmPformfact = 1 + pow(tmPffterm, tmPexp);
    tmPeffLambdaU *= pow(tmPformfact,0.25);
  }

  // ME from spin-2 unparticle.
  double tmPsLambda2 = sH / pow2(tmPeffLambdaU);
  double tmPexp = eDdU - 2;
  double tmPA = -eDlambda2chi * pow(tmPsLambda2,tmPexp)
               / (8 * pow(tmPeffLambdaU,4));
  eDsigma0 = 4 * pow2(tmPA) * uH * tH * (pow2(uH) + pow2(tH));

  // extra 1/sHS from 2-to-2 phase space.
  eDsigma0 /= 16 * M_PI * pow2(sH);

  // sigma(ffbar->llbar) = 3 * sigma(ffbar->eebar)
  eDsigma0 *= 3.;

}

//--------------------------------------------------------------------------

void Sigma2gg2LEDllbar::setIdColAcol() {

  double tmPrand = rndmPtr->flat();
  // Flavours trivial.
  if (tmPrand < 0.33333333) {      setId( 21, 21, 11, -11); }
  else if (tmPrand < 0.66666667) { setId( 21, 21, 13, -13); }
  else {                            setId( 21, 21, 15, -15); }

  // Colour flow topologies.
  setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

}

//==========================================================================

// Sigma2gg2LEDgg class.
// Cross section for g g -> (LED G*) -> g g.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2LEDgg::initProc() {

  // Read model parameters.
  eDopMode   = settingsPtr->mode("ExtraDimensionsLED:opMode");
  eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
  eDMD       = settingsPtr->parm("ExtraDimensionsLED:MD");
  eDLambdaT  = settingsPtr->parm("ExtraDimensionsLED:LambdaT");
  eDnegInt   = settingsPtr->mode("ExtraDimensionsLED:NegInt");
  eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
  eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence.

void Sigma2gg2LEDgg::sigmaKin() {

  // Get S(x) values for G amplitude.
  complex sS(0., 0.);
  complex sT(0., 0.);
  complex sU(0., 0.);
  if (eDopMode == 0) {
    sS = ampLedS( sH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sT = ampLedS( tH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sU = ampLedS( uH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
  } else {
    // Form factor.
    double effLambda = eDLambdaT;
    if ((eDcutoff == 2) || (eDcutoff == 3)) {
      double ffterm = sqrt(Q2RenSave) / (eDtff * eDLambdaT);
      double exp    = double(eDnGrav) + 2.;
      double formfa = 1. + pow(ffterm, exp);
      effLambda *= pow(formfa,0.25);
    }
    sS = 4.*M_PI/pow(effLambda,4);
    sT = 4.*M_PI/pow(effLambda,4);
    sU = 4.*M_PI/pow(effLambda,4);
    if (eDnegInt == 1) {
      sS *= -1.;
      sT *= -1.;
      sU *= -1.;
    }
  }

  // Calculate kinematics dependence.
  double sH3 = sH*sH2;
  double tH3 = tH*tH2;
  double uH3 = uH*uH2;

  sigTS  = (128. * pow2(M_PI) * pow2(alpS)) * (9./4.)
         * (tH2 / sH2 + 2. * tH / sH + 3. + 2. * sH / tH + sH2 / tH2)
         + 24.*M_PI*alpS*( (sH3/tH + tH2 + 3.*(sH*tH + sH2))*sS.real()
                         + (tH3/sH + sH2 + 3.*(tH*sH + tH2))*sT.real())
         + pow2(uH2)*( 4.*real(sS*conj(sS)) + sS.real()*sT.real()
                     + sS.imag()*sT.imag() + 4.*real(sT*conj(sT)));


  sigUS  = (128. * pow2(M_PI) * pow2(alpS)) * (9./4.)
         * (uH2 / sH2 + 2. * uH / sH + 3. + 2. * sH / uH + sH2 / uH2)
         + 24.*M_PI*alpS*( (sH3/uH + uH2 + 3.*(sH*uH + sH2))*sS.real()
                         + (uH3/sH + sH2 + 3.*(uH*sH + uH2))*sU.real())
         + pow2(tH2)*( 4.*real(sS*conj(sS)) + sS.real()*sU.real()
                     + sS.imag()*sU.imag() + 4.*real(sU*conj(sU)));

  sigTU  = (128. * pow2(M_PI) * pow2(alpS)) * (9./4.)
         * (tH2 / uH2 + 2. * tH / uH + 3. + 2. * uH / tH + uH2 / tH2)
         + 24.*M_PI*alpS*( (tH3/uH + uH2 + 3.*(tH*uH + tH2))*sT.real()
                         + (uH3/tH + tH2 + 3.*(uH*tH + uH2))*sU.real())
         + pow2(sH2)*( 4.*real(sT*conj(sT)) + sT.real()*sU.real()
                     + sT.imag()*sU.imag() + 4.*real(sU*conj(sU)));

  sigSum = sigTS + sigUS + sigTU;

  // Answer contains factor 1/2 from identical gluons.
  sigma  = 0.5 * sigSum / (128. * M_PI * sH2);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2LEDgg::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, 21, 21);

  // Three colour flow topologies, each with two orientations.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 4, 4, 3);
  else if (sigRand < sigTS + sigUS)
                       setColAcol( 1, 2, 3, 1, 3, 4, 4, 2);
  else                 setColAcol( 1, 2, 3, 4, 1, 4, 3, 2);
  if (rndmPtr->flat() > 0.5) swapColAcol();

}

//==========================================================================

// Sigma2gg2LEDqqbar class.
// Cross section for g g -> (LED G*) -> q qbar.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2LEDqqbar::initProc() {

  // Read number of quarks to be considered in massless approximation
  // as well as model parameters.
  nQuarkNew  = settingsPtr->mode("ExtraDimensionsLED:nQuarkNew");
  eDopMode   = settingsPtr->mode("ExtraDimensionsLED:opMode");
  eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
  eDMD       = settingsPtr->parm("ExtraDimensionsLED:MD");
  eDLambdaT  = settingsPtr->parm("ExtraDimensionsLED:LambdaT");
  eDnegInt   = settingsPtr->mode("ExtraDimensionsLED:NegInt");
  eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
  eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence.

void Sigma2gg2LEDqqbar::sigmaKin() {

  // Get S(x) values for G amplitude.
  complex sS(0., 0.);
  complex sT(0., 0.);
  complex sU(0., 0.);
  if (eDopMode == 0) {
    sS = ampLedS( sH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sT = ampLedS( tH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sU = ampLedS( uH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
  } else {
    // Form factor.
    double effLambda = eDLambdaT;
    if ((eDcutoff == 2) || (eDcutoff == 3)) {
      double ffterm = sqrt(Q2RenSave) / (eDtff * eDLambdaT);
      double exp    = double(eDnGrav) + 2.;
      double formfa = 1. + pow(ffterm, exp);
      effLambda *= pow(formfa,0.25);
    }
    sS = 4.*M_PI/pow(effLambda,4);
    sT = 4.*M_PI/pow(effLambda,4);
    sU = 4.*M_PI/pow(effLambda,4);
    if (eDnegInt == 1) {
      sS *= -1.;
      sT *= -1.;
      sU *= -1.;
    }
  }

  // Pick new flavour.
  idNew = 1 + int( nQuarkNew * rndmPtr->flat() );
  mNew  = particleDataPtr->m0(idNew);
  m2New = mNew*mNew;

  // Calculate kinematics dependence.
  sigTS = 0.;
  sigUS = 0.;
  if (sH > 4. * m2New) {
    double tH3 = tH*tH2;
    double uH3 = uH*uH2;
    sigTS = (16. * pow2(M_PI) * pow2(alpS))
          * ((1./6.) * uH / tH - (3./8.) * uH2 / sH2)
          - 0.5 * M_PI * alpS * uH2 * sS.real()
          + (3./16.) * uH3 * tH * real(sS*conj(sS));
    sigUS = (16. * pow2(M_PI) * pow2(alpS))
          * ((1./6.) * tH / uH - (3./8.) * tH2 / sH2)
          - 0.5 * M_PI * alpS * tH2 * sS.real()
          + (3./16.) * tH3 * uH * real(sS*conj(sS));
  }
  sigSum = sigTS + sigUS;

  // Answer is proportional to number of outgoing flavours.
  sigma  = nQuarkNew * sigSum / (16. * M_PI * sH2);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2LEDqqbar::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idNew, -idNew);

  // Two colour flow topologies.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else                 setColAcol( 1, 2, 3, 1, 3, 0, 0, 2);

}

//==========================================================================

// Sigma2qg2LEDqg class.
// Cross section for q g -> (LED G*) -> q g.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qg2LEDqg::initProc() {

  // Read model parameters.
  eDopMode   = settingsPtr->mode("ExtraDimensionsLED:opMode");
  eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
  eDMD       = settingsPtr->parm("ExtraDimensionsLED:MD");
  eDLambdaT  = settingsPtr->parm("ExtraDimensionsLED:LambdaT");
  eDnegInt   = settingsPtr->mode("ExtraDimensionsLED:NegInt");
  eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
  eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence.

void Sigma2qg2LEDqg::sigmaKin() {

  // Get S(x) values for G amplitude.
  complex sS(0., 0.);
  complex sT(0., 0.);
  complex sU(0., 0.);
  if (eDopMode == 0) {
    sS = ampLedS( sH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sT = ampLedS( tH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sU = ampLedS( uH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
  } else {
    // Form factor.
    double effLambda = eDLambdaT;
    if ((eDcutoff == 2) || (eDcutoff == 3)) {
      double ffterm = sqrt(Q2RenSave) / (eDtff * eDLambdaT);
      double exp    = double(eDnGrav) + 2.;
      double formfa = 1. + pow(ffterm, exp);
      effLambda *= pow(formfa,0.25);
    }
    sS = 4.*M_PI/pow(effLambda,4);
    sT = 4.*M_PI/pow(effLambda,4);
    sU = 4.*M_PI/pow(effLambda,4);
    if (eDnegInt == 1) {
      sS *= -1.;
      sT *= -1.;
      sU *= -1.;
    }
  }

  // Calculate kinematics dependence.
  double sH3 = sH*sH2;
  double uH3 = uH*uH2;
  sigTS  = (16. * pow2(M_PI) * pow2(alpS))
         * (uH2 / tH2 - (4./9.) * uH / sH)
         + (4./3.) * M_PI * alpS * uH2 * sT.real()
         - 0.5 * uH3 * sH * real(sT*conj(sT));
  sigTU  = (16. * pow2(M_PI) * pow2(alpS))
         * (sH2 / tH2 - (4./9.) * sH / uH)
         + (4./3.) * M_PI * alpS * sH2 * sT.real()
         - 0.5 * sH3 * uH * real(sT*conj(sT));
  sigSum = sigTS + sigTU;

  // Answer.
  sigma  = sigSum / (16. * M_PI * sH2);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2LEDqg::setIdColAcol() {

  // Outgoing = incoming flavours.
  setId( id1, id2, id1, id2);

  // Two colour flow topologies. Swap if first is gluon, or when antiquark.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 2, 1, 3, 0, 2, 3);
  else                 setColAcol( 1, 0, 2, 3, 2, 0, 1, 3);
  if (id1 == 21) swapCol1234();
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qq2LEDqq class.
// Cross section for q q(bar)' -> (LED G*) -> q q(bar)'

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qq2LEDqq::initProc() {

  // Read model parameters.
  eDopMode   = settingsPtr->mode("ExtraDimensionsLED:opMode");
  eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
  eDMD       = settingsPtr->parm("ExtraDimensionsLED:MD");
  eDLambdaT  = settingsPtr->parm("ExtraDimensionsLED:LambdaT");
  eDnegInt   = settingsPtr->mode("ExtraDimensionsLED:NegInt");
  eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
  eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour.

void Sigma2qq2LEDqq::sigmaKin() {

  // Get S(x) values for G amplitude.
  complex sS(0., 0.);
  complex sT(0., 0.);
  complex sU(0., 0.);
  if (eDopMode == 0) {
    sS = ampLedS( sH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sT = ampLedS( tH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sU = ampLedS( uH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
  } else {
    // Form factor.
    double effLambda = eDLambdaT;
    if ((eDcutoff == 2) || (eDcutoff == 3)) {
      double ffterm = sqrt(Q2RenSave) / (eDtff * eDLambdaT);
      double exp    = double(eDnGrav) + 2.;
      double formfa = 1. + pow(ffterm, exp);
      effLambda *= pow(formfa,0.25);
    }
    sS = 4.*M_PI/pow(effLambda,4);
    sT = 4.*M_PI/pow(effLambda,4);
    sU = 4.*M_PI/pow(effLambda,4);
    if (eDnegInt == 1) {
      sS *= -1.;
      sT *= -1.;
      sU *= -1.;
    }
  }

  // Calculate kinematics dependence for different terms.
  sigT   = (4./9.) * (sH2 + uH2) / tH2;
  sigU   = (4./9.) * (sH2 + tH2) / uH2;
  sigTU  = - (8./27.) * sH2 / (tH * uH);
  sigST  = - (8./27.) * uH2 / (sH * tH);
  // Graviton terms.
  sigGrT1 = funLedG(tH, uH) * real(sT*conj(sT)) / 8.;
  sigGrT2 = funLedG(tH, sH) * real(sT*conj(sT)) / 8.;
  sigGrU  = funLedG(uH, tH) * real(sU*conj(sU)) / 8.;
  sigGrTU = (8./9.) * M_PI * alpS * sH2
          * ((4.*uH + tH)*sT.real()/uH + (4.*tH + uH)*sU.real()/tH)
          + (sT.real()*sU.real() + sT.imag()*sU.imag())
          * (4.*tH + uH)*(4.*uH + tH) * sH2 / 48.;
  sigGrST = (8./9.) * M_PI * alpS * uH2
          * ((4.*tH + sH)*sS.real()/tH + (4.*sH + tH)*sT.real()/sH)
          + (sS.real()*sT.real() + sS.imag()*sT.imag())
          * (4.*sH + tH)*(4.*tH + sH) * uH2 / 48.;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence.

double Sigma2qq2LEDqq::sigmaHat() {

  // Combine cross section terms; factor 1/2 when identical quarks.
  if (id2 ==  id1) {
    sigSum  = (16. * pow2(M_PI) * pow2(alpS)) * (sigT + sigU + sigTU)
            + sigGrT1 + sigGrU + sigGrTU;
    sigSum *= 0.5;
  } else if (id2 == -id1) {
    sigSum = (16. * pow2(M_PI) * pow2(alpS)) * (sigT + sigST)
           + sigGrT2 + sigGrST;
  } else {
    sigSum = 16. * pow2(M_PI) * pow2(alpS) * sigT + sigGrT1;
  }

  // Answer.
  return sigSum / (16. * M_PI * sH2);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qq2LEDqq::setIdColAcol() {

  // Outgoing = incoming flavours.
  setId( id1, id2, id1, id2);

  // Colour flow topologies. Swap when antiquarks.
  double sigTtot = sigT + sigGrT2;
  double sigUtot = sigU + sigGrU;
  if (id1 * id2 > 0)  setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);
  else                setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  if (id2 == id1 && (sigTtot + sigUtot) * rndmPtr->flat() > sigTtot)
                      setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2LEDgg class.
// Cross section for q qbar -> (LED G*) -> g g.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qqbar2LEDgg::initProc() {

  // Read model parameters.
  eDopMode   = settingsPtr->mode("ExtraDimensionsLED:opMode");
  eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
  eDMD       = settingsPtr->parm("ExtraDimensionsLED:MD");
  eDLambdaT  = settingsPtr->parm("ExtraDimensionsLED:LambdaT");
  eDnegInt   = settingsPtr->mode("ExtraDimensionsLED:NegInt");
  eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
  eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence.

void Sigma2qqbar2LEDgg::sigmaKin() {

  // Get S(x) values for G amplitude.
  complex sS(0., 0.);
  complex sT(0., 0.);
  complex sU(0., 0.);
  if (eDopMode == 0) {
    sS = ampLedS( sH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sT = ampLedS( tH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sU = ampLedS( uH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
  } else {
    // Form factor.
    double effLambda = eDLambdaT;
    if ((eDcutoff == 2) || (eDcutoff == 3)) {
      double ffterm = sqrt(Q2RenSave) / (eDtff * eDLambdaT);
      double exp    = double(eDnGrav) + 2.;
      double formfa = 1. + pow(ffterm, exp);
      effLambda *= pow(formfa,0.25);
    }
    sS = 4.*M_PI/pow(effLambda,4);
    sT = 4.*M_PI/pow(effLambda,4);
    sU = 4.*M_PI/pow(effLambda,4);
    if (eDnegInt == 1) {
      sS *= -1.;
      sT *= -1.;
      sU *= -1.;
    }
  }

  // Calculate kinematics dependence.
  double tH3 = tH*tH2;
  double uH3 = uH*uH2;
  sigTS  = (16. * pow2(M_PI) * pow2(alpS))
         * ((1./6.) * uH / tH - (3./8.) * uH2 / sH2)
         - 0.5 * M_PI * alpS * uH2 * sS.real()
         + (3./16.) * uH3 * tH * real(sS*conj(sS));
  sigUS  = (16. * pow2(M_PI) * pow2(alpS))
         * ((1./6.) * tH / uH - (3./8.) * tH2 / sH2)
         - 0.5 * M_PI * alpS * tH2 * sS.real()
         + (3./16.) * tH3 * uH * real(sS*conj(sS));

  sigSum = sigTS + sigUS;

  // Answer contains factor 1/2 from identical gluons.
  sigma  = (64./9.) * 0.5 * sigSum / (16. * M_PI * sH2);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2LEDgg::setIdColAcol() {

  // Outgoing flavours trivial.
  setId( id1, id2, 21, 21);

  // Two colour flow topologies. Swap if first is antiquark.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 0, 2, 1, 3, 3, 2);
  else                 setColAcol( 1, 0, 0, 2, 3, 2, 1, 3);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2LEDqqbarNew class.
// Cross section q qbar -> (LED G*) -> q' qbar'.

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qqbar2LEDqqbarNew::initProc() {

  // Read number of quarks to be considered in massless approximation
  // as well as model parameters.
  nQuarkNew  = settingsPtr->mode("ExtraDimensionsLED:nQuarkNew");
  eDopMode   = settingsPtr->mode("ExtraDimensionsLED:opMode");
  eDnGrav    = settingsPtr->mode("ExtraDimensionsLED:n");
  eDMD       = settingsPtr->parm("ExtraDimensionsLED:MD");
  eDLambdaT  = settingsPtr->parm("ExtraDimensionsLED:LambdaT");
  eDcutoff   = settingsPtr->mode("ExtraDimensionsLED:CutOffMode");
  eDtff      = settingsPtr->parm("ExtraDimensionsLED:t");

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence.

void Sigma2qqbar2LEDqqbarNew::sigmaKin() {

  // Get S(x) values for G amplitude.
  complex sS(0., 0.);
  complex sT(0., 0.);
  complex sU(0., 0.);
  if (eDopMode == 0) {
    sS = ampLedS( sH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sT = ampLedS( tH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
    sU = ampLedS( uH/pow2(eDLambdaT), eDnGrav, eDLambdaT, eDMD);
  } else {
    // Form factor.
    double effLambda = eDLambdaT;
    if ((eDcutoff == 2) || (eDcutoff == 3)) {
      double ffterm = sqrt(Q2RenSave) / (eDtff * eDLambdaT);
      double exp    = double(eDnGrav) + 2.;
      double formfa = 1. + pow(ffterm, exp);
      effLambda *= pow(formfa,0.25);
    }
    sS = 4.*M_PI/pow(effLambda,4);
    sT = 4.*M_PI/pow(effLambda,4);
    sU = 4.*M_PI/pow(effLambda,4);
  }

  // Pick new flavour.
  idNew = 1 + int( nQuarkNew * rndmPtr->flat() );
  mNew  = particleDataPtr->m0(idNew);
  m2New = mNew*mNew;

  // Calculate kinematics dependence.
  sigS                      = 0.;
  if (sH > 4. * m2New) {
    sigS = (16. * pow2(M_PI) * pow2(alpS))
         * (4./9.) * (tH2 + uH2) / sH2
         + funLedG(sH, tH) * real(sS*conj(sS)) / 8.;
  }
  // Answer is proportional to number of outgoing flavours.
  sigma = nQuarkNew * sigS / (16. * M_PI * sH2);

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2LEDqqbarNew::setIdColAcol() {

  // Set outgoing flavours ones.
  id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, id2, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

} // end namespace Pythia8
