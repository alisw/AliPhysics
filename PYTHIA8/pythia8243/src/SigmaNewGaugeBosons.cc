// SigmaNewGaugeBosons.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// leptoquark simulation classes.

#include "Pythia8/SigmaNewGaugeBosons.h"

namespace Pythia8 {


//==========================================================================

// Sigma1ffbarZprimeWprime class.
// Collects common methods for f fbar -> Z'/W' -> WW/WZ -> 4 fermions.
// Copied from SigmaEW for gauge-boson-pair production.

//--------------------------------------------------------------------------

// Calculate and store internal products.

void Sigma1ffbarZprimeWprime::setupProd( Event& process, int i1, int i2,
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

complex Sigma1ffbarZprimeWprime::fGK(int j1, int j2, int j3, int j4,
  int j5, int j6) {

  return 4. * hA[j1][j3] * hC[j2][j6]
         * ( hA[j1][j5] * hC[j1][j4] + hA[j3][j5] * hC[j3][j4] );

}

//--------------------------------------------------------------------------

// Evaluate the Xi function of Gunion and Kunszt.

double Sigma1ffbarZprimeWprime::xiGK( double tHnow, double uHnow,
  double s3now, double s4now) {

  return - 4. * s3now * s4now + tHnow * (3. * tHnow + 4. * uHnow)
    + tHnow * tHnow * ( tHnow * uHnow / (s3now * s4now)
    - 2. * (1. / s3now + 1./s4now) * (tHnow + uHnow)
    + 2. * (s3now / s4now + s4now / s3now) );

}

//--------------------------------------------------------------------------

// Evaluate the Xj function of Gunion and Kunszt.

double Sigma1ffbarZprimeWprime::xjGK( double tHnow, double uHnow,
  double s3now, double s4now) {

  return 8. * pow2(s3now + s4now) - 8. * (s3now + s4now) * (tHnow + uHnow)
    - 6. * tHnow * uHnow - 2. * tHnow * uHnow * ( tHnow * uHnow
    / (s3now * s4now) - 2. * (1. / s3now + 1. / s4now) * (tHnow + uHnow)
    + 2. * (s3now / s4now + s4now / s3now) );

}

//==========================================================================

// Sigma1ffbar2gmZZprime class.
// Cross section for f fbar -> gamma*/Z0/Z'0 (f is quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ffbar2gmZZprime::initProc() {

  // Allow to pick only parts of full gamma*/Z0/Z'0 expression.
  gmZmode     = settingsPtr->mode("Zprime:gmZmode");

  // Store Z'0 mass and width for propagator.
  mRes        = particleDataPtr->m0(32);
  GammaRes    = particleDataPtr->mWidth(32);
  m2Res       = mRes*mRes;
  GamMRat     = GammaRes / mRes;
  sin2tW      = couplingsPtr->sin2thetaW();
  cos2tW      = 1. - sin2tW;
  thetaWRat   = 1. / (16. * sin2tW * cos2tW);

  // Properties of Z0 resonance also needed.
  mZ          = particleDataPtr->m0(23);
  GammaZ      = particleDataPtr->mWidth(23);
  m2Z         = mZ*mZ;
  GamMRatZ    = GammaZ / mZ;

  // Ensure that arrays initially are empty.
  for (int i = 0; i < 20; ++i) afZp[i] = 0.;
  for (int i = 0; i < 20; ++i) vfZp[i] = 0.;

  // Store first-generation axial and vector couplings.
  afZp[1]     = settingsPtr->parm("Zprime:ad");
  afZp[2]     = settingsPtr->parm("Zprime:au");
  afZp[11]    = settingsPtr->parm("Zprime:ae");
  afZp[12]    = settingsPtr->parm("Zprime:anue");
  vfZp[1]     = settingsPtr->parm("Zprime:vd");
  vfZp[2]     = settingsPtr->parm("Zprime:vu");
  vfZp[11]    = settingsPtr->parm("Zprime:ve");
  vfZp[12]    = settingsPtr->parm("Zprime:vnue");

  // Determine if the 4th generation should be included
  bool coupZp2gen4    = settingsPtr->flag("Zprime:coup2gen4");

  maxZpGen = (coupZp2gen4) ? 8 : 6;

  // Second and third (and possibly 4th) generation could be carbon copy
  // of this...
  if (settingsPtr->flag("Zprime:universality")) {
    for (int i = 3; i <= maxZpGen; ++i) {
      afZp[i]    = afZp[i-2];
      vfZp[i]    = vfZp[i-2];
      afZp[i+10] = afZp[i+8];
      vfZp[i+10] = vfZp[i+8];
    }

  // ... or could have different couplings.
  } else {
    afZp[3]   = settingsPtr->parm("Zprime:as");
    afZp[4]   = settingsPtr->parm("Zprime:ac");
    afZp[5]   = settingsPtr->parm("Zprime:ab");
    afZp[6]   = settingsPtr->parm("Zprime:at");
    afZp[13]  = settingsPtr->parm("Zprime:amu");
    afZp[14]  = settingsPtr->parm("Zprime:anumu");
    afZp[15]  = settingsPtr->parm("Zprime:atau");
    afZp[16]  = settingsPtr->parm("Zprime:anutau");
    vfZp[3]   = settingsPtr->parm("Zprime:vs");
    vfZp[4]   = settingsPtr->parm("Zprime:vc");
    vfZp[5]   = settingsPtr->parm("Zprime:vb");
    vfZp[6]   = settingsPtr->parm("Zprime:vt");
    vfZp[13]  = settingsPtr->parm("Zprime:vmu");
    vfZp[14]  = settingsPtr->parm("Zprime:vnumu");
    vfZp[15]  = settingsPtr->parm("Zprime:vtau");
    vfZp[16]  = settingsPtr->parm("Zprime:vnutau");
    if( coupZp2gen4 ) {
      afZp[7]   = settingsPtr->parm("Zprime:abPrime");
      afZp[8]   = settingsPtr->parm("Zprime:atPrime");
      vfZp[7]   = settingsPtr->parm("Zprime:vbPrime");
      vfZp[8]   = settingsPtr->parm("Zprime:vtPrime");
      afZp[17]  = settingsPtr->parm("Zprime:atauPrime");
      afZp[18]  = settingsPtr->parm("Zprime:anutauPrime");
      vfZp[17]  = settingsPtr->parm("Zprime:vtauPrime");
      vfZp[18]  = settingsPtr->parm("Zprime:vnutauPrime");
    }
  }

  // Coupling for Z' -> W+ W- and decay angular admixture.
  coupZpWW    = settingsPtr->parm("Zprime:coup2WW");
  anglesZpWW  = settingsPtr->parm("Zprime:anglesWW");

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(32);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1ffbar2gmZZprime::sigmaKin() {

  // Common coupling factors.
  double colQ = 3. * (1. + alpS / M_PI);

  // Reset quantities to sum. Declare variables in loop.
  gamSum   = 0.;
  gamZSum  = 0.;
  ZSum     = 0.;
  gamZpSum = 0.;
  ZZpSum   = 0.;
  ZpSum    = 0.;
  int    idAbs, onMode;
  double mf, mr, ps, kinFacA, kinFacV, ef, af, vf, apf, vpf,
         ef2, efvf, vaf2, efvpf, vafvapf, vapf2, colf;

  // Loop over all open Z'0 decay channels.
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    onMode = particlePtr->channel(i).onMode();
    if (onMode != 1 && onMode != 2) continue;
    idAbs = abs( particlePtr->channel(i).product(0) );

    // Contributions from three/four fermion generations,
    // and optionally also from excited fermions.
    if ( (idAbs >  0 && idAbs <= maxZpGen)
      || (idAbs > 10 && idAbs <= maxZpGen+10)
      || (idAbs > 4000000 && idAbs <= 4000006)
      || (idAbs > 4000010 && idAbs <= 4000016) ) {
      int idAbs4 = (idAbs < 4000000) ? idAbs : idAbs - 4000000;
      mf = particleDataPtr->m0(idAbs);

      // Check that above threshold.
      if (mH > 2. * mf + MASSMARGIN) {
        mr        = pow2(mf / mH);
        ps        = sqrtpos(1. - 4. * mr);

        // Couplings of gamma^*/Z^0/Z'^0  to final flavour
        ef        = couplingsPtr->ef(idAbs4);
        af        = couplingsPtr->af(idAbs4);
        vf        = couplingsPtr->vf(idAbs4);
        apf       = afZp[idAbs4];
        vpf       = vfZp[idAbs4];

        // Combine couplings with kinematical factors.
        kinFacA   = pow3(ps);
        kinFacV   = ps * (1. + 2. * mr);
        ef2       = ef * ef * kinFacV;
        efvf      = ef * vf * kinFacV;
        vaf2      = vf * vf * kinFacV + af * af * kinFacA;
        efvpf     = ef * vpf * kinFacV;
        vafvapf   = vf * vpf * kinFacV + af * apf * kinFacA;
        vapf2     = vpf * vpf * kinFacV + apf * apf * kinFacA;

        // Colour factor. Additionally secondary width for heavy particles.
        colf      = (idAbs4 < 9) ? colQ : 1.;
        if ( (idAbs > 5 && idAbs < 9) || (idAbs > 17 && idAbs < 19)
          || idAbs > 4000000)
          colf *= particleDataPtr->resOpenFrac(idAbs, -idAbs);

        // Store sum of combinations.
        gamSum   += colf * ef2;
        gamZSum  += colf * efvf;
        ZSum     += colf * vaf2;
        gamZpSum += colf * efvpf;
        ZZpSum   += colf * vafvapf;
        ZpSum    += colf * vapf2;
      }

    // Optional contribution from W+ W-.
    } else if (idAbs == 24) {
      mf = particleDataPtr->m0(idAbs);
      if (mH > 2. * mf + MASSMARGIN) {
        mr        = pow2(mf / mH);
        ps        = sqrtpos(1. - 4. * mr);
        ZpSum    += pow2(coupZpWW * cos2tW) * pow3(ps)
                  * (1. + 20. * mr + 12. * mr*mr)
                  * particleDataPtr->resOpenFrac(24, -24);
      }
    }
  }

  // Calculate prefactors for gamma/Z0/Z'0 cross section terms.
  double propZ  = sH / ( pow2(sH - m2Z) + pow2(sH * GamMRatZ) );
  double propZp = sH / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  gamNorm   = 4. * M_PI * pow2(alpEM) / (3. * sH);
  gamZNorm  = gamNorm * 2. * thetaWRat * (sH - m2Z) * propZ;
  ZNorm     = gamNorm * pow2(thetaWRat) * sH * propZ;
  gamZpNorm = gamNorm * 2. * thetaWRat * (sH - m2Res) * propZp;
  ZZpNorm   = gamNorm * 2. * pow2(thetaWRat) * ((sH - m2Z) * (sH - m2Res)
              + sH * GamMRatZ * sH * GamMRat) * propZ * propZp;
  ZpNorm    = gamNorm * pow2(thetaWRat) * sH * propZp;

  // Optionally only keep some of gamma*, Z0 and Z' terms.
  if (gmZmode == 1) {gamZNorm = 0; ZNorm = 0.; gamZpNorm = 0.;
    ZZpNorm = 0.; ZpNorm = 0.;}
  if (gmZmode == 2) {gamNorm = 0.; gamZNorm = 0.; gamZpNorm = 0.;
    ZZpNorm = 0.; ZpNorm = 0.;}
  if (gmZmode == 3) {gamNorm = 0.; gamZNorm = 0.; ZNorm = 0.;
    gamZpNorm = 0.; ZZpNorm = 0.;}
  if (gmZmode == 4) {gamZpNorm = 0.; ZZpNorm = 0.; ZpNorm = 0.;}
  if (gmZmode == 5) {gamZNorm = 0.; ZNorm = 0.; ZZpNorm = 0.;}
  if (gmZmode == 6) {gamNorm = 0.; gamZNorm = 0.; gamZpNorm = 0.;}

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma1ffbar2gmZZprime::sigmaHat() {

  // Couplings to an incoming flavour.
  int idAbs      = abs(id1);
  double ei      = couplingsPtr->ef(idAbs);
  double ai      = couplingsPtr->af(idAbs);
  double vi      = couplingsPtr->vf(idAbs);
  double api     = afZp[idAbs];
  double vpi     = vfZp[idAbs];
  double ei2     = ei * ei;
  double eivi    = ei * vi;
  double vai2    = vi * vi + ai * ai;
  double eivpi   = ei * vpi;
  double vaivapi = vi * vpi + ai * api;;
  double vapi2   = vpi * vpi + api * api;

  // Combine gamma, interference and Z0 parts.
  double sigma = ei2 * gamNorm * gamSum + eivi * gamZNorm * gamZSum
               + vai2 * ZNorm * ZSum + eivpi * gamZpNorm * gamZpSum
               + vaivapi * ZZpNorm * ZZpSum + vapi2 * ZpNorm * ZpSum;

  // Colour factor. Answer.
  if (idAbs < 9) sigma /= 3.;
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2gmZZprime::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 32);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for gamma*/Z0/Z'0 decay angle.

double Sigma1ffbar2gmZZprime::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Default values, in- and out-flavours in process.
  double wt    = 1.;
  double wtMax = 1.;
  int idInAbs  = process[3].idAbs();
  int idOutAbs = process[6].idAbs();

  // Angular weight for outgoing fermion pair.
  if (iResBeg == 5 && iResEnd == 5 && (idOutAbs <= maxZpGen
    || (idOutAbs > 10 && idOutAbs <= maxZpGen+10) || idOutAbs > 4000000) ) {

    // Couplings for in- and out-flavours.
    double ei  = couplingsPtr->ef(idInAbs);
    double vi  = couplingsPtr->vf(idInAbs);
    double ai  = couplingsPtr->af(idInAbs);
    double vpi = vfZp[idInAbs];
    double api = afZp[idInAbs];
    int idOutAbs4 = (idOutAbs < 4000000) ? idOutAbs : idOutAbs - 4000000;
    double ef  = couplingsPtr->ef(idOutAbs4);
    double vf  = couplingsPtr->vf(idOutAbs4);
    double af  = couplingsPtr->af(idOutAbs4);
    double vpf = vfZp[idOutAbs4];
    double apf = afZp[idOutAbs4];

    // Phase space factors. (One power of beta left out in formulae.)
    double mr1 = pow2(process[6].m()) / sH;
    double mr2 = pow2(process[7].m()) / sH;
    double ps  = sqrtpos(pow2(1. - mr1 - mr2) - 4. * mr1 * mr2);
    double mrAvg = 0.5 * (mr1 + mr2) - 0.25 * pow2(mr1 - mr2);

    // Coefficients of angular expression.
    double coefTran = ei*ei * gamNorm * ef*ef + ei * vi * gamZNorm * ef * vf
      + (vi*vi + ai*ai) * ZNorm * (vf*vf + ps*ps * af*af)
      + ei * vpi * gamZpNorm * ef * vpf
      + (vi * vpi + ai * api) * ZZpNorm * (vf * vpf + ps*ps * af * apf)
      + (vpi*vpi + api*api) * ZpNorm * (vpf*vpf + ps*ps * apf*apf);
    double coefLong = 4. * mrAvg * ( ei*ei * gamNorm * ef*ef
      + ei * vi * gamZNorm * ef * vf + (vi*vi + ai*ai) * ZNorm * vf*vf
      + ei * vpi * gamZpNorm * ef * vpf
      + (vi * vpi + ai * api) * ZZpNorm * vf * vpf
      + (vpi*vpi + api*api) * ZpNorm * vpf*vpf );
    double coefAsym = ps * ( ei * ai * gamZNorm * ef * af
      + 4. * vi * ai * ZNorm * vf * af + ei * api * gamZpNorm * ef * apf
      + (vi * api + vpi * ai) * ZZpNorm * (vf * apf + vpf * af)
      + 4. * vpi * api * ZpNorm * vpf * apf );

    // Flip asymmetry for in-fermion + out-antifermion.
    if (process[3].id() * process[6].id() < 0) coefAsym = -coefAsym;

    // Reconstruct decay angle and weight for it.
    double cosThe = (process[3].p() - process[4].p())
      * (process[7].p() - process[6].p()) / (sH * ps);
    wt    = coefTran * (1. + pow2(cosThe))
       + coefLong * (1. - pow2(cosThe)) + 2. * coefAsym * cosThe;
    wtMax = 2. * (coefTran + abs(coefAsym));
  }

  // Angular weight for Z' -> W+ W-.
  else if (iResBeg == 5 && iResEnd == 5 && idOutAbs == 24) {
    double mr1 = pow2(process[6].m()) / sH;
    double mr2 = pow2(process[7].m()) / sH;
    double ps  = sqrtpos(pow2(1. - mr1 -mr2) - 4. * mr1 * mr2);
    double cCos2 = - (1./16.) * ps*ps * (1. - 2. * mr1 - 2. * mr2
      + mr1*mr1 + mr2*mr2 + 10. * mr1 * mr2);
    double cFlat = -cCos2 + 0.5 * (mr1 + mr2)
      * (1. - 2. * mr1 - 2. * mr2 + pow2(mr1 - mr2));

    // Reconstruct decay angle and weight for it.
    double cosThe = (process[3].p() - process[4].p())
      * (process[7].p() - process[6].p()) / (sH * ps);
    wt    = cFlat + cCos2 * cosThe*cosThe;
    wtMax = cFlat + max(0., cCos2);
  }

  // Angular weight for f + fbar -> Z' -> W+ + W- -> 4 fermions.
  else if (iResBeg == 6 && iResEnd == 7 && idOutAbs == 24) {

    // Order so that fbar(1) f(2) -> f'(3) fbar'(4) f"(5) fbar"(6).
    // with f' fbar' from W- and f" fbar" from W+.
    int i1 = (process[3].id() < 0) ? 3 : 4;
    int i2 = 7 - i1;
    int i3 = (process[8].id() > 0) ? 8 : 9;
    int i4 = 17 - i3;
    int i5 = (process[10].id() > 0) ? 10 : 11;
    int i6 = 21 - i5;
    if (process[6].id() > 0) {swap(i3, i5); swap(i4, i6);}

    // Decay distribution like in f fbar -> Z^* -> W+ W-.
    if (rndmPtr->flat() > anglesZpWW) {

      // Set up four-products and internal products.
      setupProd( process, i1, i2, i3, i4, i5, i6);

      // tHat and uHat of fbar f -> W- W+, and their squared masses.
      int iNeg     = (process[6].id() < 0) ? 6 : 7;
      int iPos     = 13 - iNeg;
      double tHres = (process[i1].p() - process[iNeg].p()).m2Calc();
      double uHres = (process[i1].p() - process[iPos].p()).m2Calc();
      double s3now = process[iNeg].m2();
      double s4now = process[iPos].m2();

      // Kinematics combinations (norm(x) = |x|^2).
      double fGK135 = norm(fGK( 1, 2, 3, 4, 5, 6) - fGK( 1, 2, 5, 6, 3, 4) );
      double fGK253 = norm(fGK( 2, 1, 5, 6, 3, 4) - fGK( 2, 1, 3, 4, 5, 6) );
      double xiT    = xiGK( tHres, uHres, s3now, s4now);
      double xiU    = xiGK( uHres, tHres, s3now, s4now);
      double xjTU   = xjGK( tHres, uHres, s3now, s4now);

      //  Couplings of incoming (anti)fermion. Combine with kinematics.
      int idAbs     = process[i1].idAbs();
      double li     = 0.5 * (vfZp[idAbs] + afZp[idAbs]);
      double ri     = 0.5 * (vfZp[idAbs] - afZp[idAbs]);
      wt            = li*li * fGK135 + ri*ri * fGK253;
      wtMax         = 4. * s3now * s4now * (li*li + ri*ri)
                    * (xiT + xiU - xjTU);

    // Decay distribution like in f fbar -> h^0 -> W+ W-.
    } else {
      double p35  = 2. * process[i3].p() * process[i5].p();
      double p46  = 2. * process[i4].p() * process[i6].p();
      wt          = 16. * p35 * p46;
      wtMax       = sH2;
    }
  }

  // Angular weight in top decay by standard routine.
  else if (process[process[iResBeg].mother1()].idAbs() == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // Angular weight for fourth generation or excited fermions not implemented.

  // Done.
  return (wt / wtMax);

}

//==========================================================================

// Sigma1ffbar2Wprime class.
// Cross section for f fbar' -> W'+- (f is quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ffbar2Wprime::initProc() {

  // Store W+- mass and width for propagator.
  mRes     = particleDataPtr->m0(34);
  GammaRes = particleDataPtr->mWidth(34);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;
  thetaWRat = 1. / (12. * couplingsPtr->sin2thetaW());

  // Axial and vector couplings of fermions.
  aqWp      = settingsPtr->parm("Wprime:aq");
  vqWp      = settingsPtr->parm("Wprime:vq");
  alWp      = settingsPtr->parm("Wprime:al");
  vlWp      = settingsPtr->parm("Wprime:vl");

  // Coupling for W' -> W Z and decay angular admixture.
  coupWpWZ    = settingsPtr->parm("Wprime:coup2WZ");
  anglesWpWZ  = settingsPtr->parm("Wprime:anglesWZ");

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(34);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1ffbar2Wprime::sigmaKin() {

  // Set up Breit-Wigner. Cross section for W+ and W- separately.
  double sigBW  = 12. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double preFac = alpEM * thetaWRat * mH;
  sigma0Pos     = preFac * sigBW * particlePtr->resWidthOpen(34, mH);
  sigma0Neg     = preFac * sigBW * particlePtr->resWidthOpen(-34, mH);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma1ffbar2Wprime::sigmaHat() {

  // Secondary width for W+ or W-. CKM and colour factors.
  int idUp = (abs(id1)%2 == 0) ? id1 : id2;
  double sigma = (idUp > 0) ? sigma0Pos : sigma0Neg;
  if (abs(id1) < 7) sigma *= couplingsPtr->V2CKMid(abs(id1), abs(id2)) / 3.;

  // Couplings.
  if (abs(id1) < 7) sigma *= 0.5 * (aqWp * aqWp + vqWp * vqWp);
  else              sigma *= 0.5 * (alWp * alWp + vlWp * vlWp);

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2Wprime::setIdColAcol() {

  // Sign of outgoing W.
  int sign          = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 34 * sign);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for W decay angle.

double Sigma1ffbar2Wprime::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Default values, in- and out-flavours in process.
  double wt    = 1.;
  double wtMax = 1.;
  int idInAbs  = process[3].idAbs();
  int idOutAbs = process[6].idAbs();

  // Angular weight for outgoing fermion pair.
  if (iResBeg == 5 && iResEnd == 5 &&
    (idOutAbs < 7 || ( idOutAbs > 10 && idOutAbs < 17)) ) {

    // Couplings for in- and out-flavours.
    double ai  = (idInAbs < 9) ? aqWp : alWp;
    double vi  = (idInAbs < 9) ? vqWp : vlWp;
    double af  = (idOutAbs < 9) ? aqWp : alWp;
    double vf  = (idOutAbs < 9) ? vqWp : vlWp;

    // Asymmetry expression.
    double coefAsym = 8. * vi * ai * vf * af
      / ((vi*vi + ai*ai) * (vf*vf + af*af));

    // Flip asymmetry for in-fermion + out-antifermion.
    if (process[3].id() * process[6].id() < 0) coefAsym = -coefAsym;

    // Phase space factors.
    double mr1 = pow2(process[6].m()) / sH;
    double mr2 = pow2(process[7].m()) / sH;
    double ps  = sqrtpos(pow2(1. - mr1 - mr2) - 4. * mr1 * mr2);

    // Reconstruct decay angle and weight for it.
    double cosThe = (process[3].p() - process[4].p())
      * (process[7].p() - process[6].p()) / (sH * ps);
    wt    = 1. + coefAsym * cosThe + cosThe * cosThe;
    wtMax = 2. + abs(coefAsym);
  }

  // Angular weight for W' -> W Z.
  else if (iResBeg == 5 && iResEnd == 5 && idOutAbs == 24) {
    double mr1 = pow2(process[6].m()) / sH;
    double mr2 = pow2(process[7].m()) / sH;
    double ps  = sqrtpos(pow2(1. - mr1 - mr2) - 4. * mr1 * mr2);
    double cCos2 = - (1./16.) * ps*ps * (1. - 2. * mr1 - 2. * mr2
      + mr1*mr1 + mr2*mr2 + 10. * mr1 * mr2);
    double cFlat = -cCos2 + 0.5 * (mr1 + mr2)
      * (1. - 2. * mr1 - 2. * mr2 + pow2(mr1 - mr2));

    // Reconstruct decay angle and weight for it.
    double cosThe = (process[3].p() - process[4].p())
      * (process[7].p() - process[6].p()) / (sH * ps);
    wt    = cFlat + cCos2 * cosThe*cosThe;
    wtMax = cFlat + max(0., cCos2);
  }

  // Angular weight for f + fbar -> W' -> W + Z -> 4 fermions.
  else if (iResBeg == 6 && iResEnd == 7
    && (idOutAbs == 24 || idOutAbs == 23)) {

    // Order so that fbar(1) f(2) -> f'(3) fbar'(4) f"(5) fbar"(6).
    // with f' fbar' from W and f" fbar" from Z.
    int i1 = (process[3].id() < 0) ? 3 : 4;
    int i2 = 7 - i1;
    int i3 = (process[8].id() > 0) ? 8 : 9;
    int i4 = 17 - i3;
    int i5 = (process[10].id() > 0) ? 10 : 11;
    int i6 = 21 - i5;
    if (process[6].id() == 23) {swap(i3, i5); swap(i4, i6);}

    // Decay distribution like in f fbar -> Z^* -> W+ W-.
    if (rndmPtr->flat() > anglesWpWZ) {

      // Set up four-products and internal products.
      setupProd( process, i1, i2, i3, i4, i5, i6);

      // tHat and uHat of fbar f -> W Z, and their squared masses.
      int iW       = (process[6].id() == 23) ? 7 : 6;
      int iZ       = 13 - iW;
      double tHres = (process[i1].p() - process[iW].p()).m2Calc();
      double uHres = (process[i1].p() - process[iZ].p()).m2Calc();
      double s3now = process[iW].m2();
      double s4now = process[iZ].m2();

      // Kinematics combinations (norm(x) = |x|^2).
      double fGK135 = norm(fGK( 1, 2, 3, 4, 5, 6) - fGK( 1, 2, 5, 6, 3, 4) );
      double fGK136 = norm(fGK( 1, 2, 3, 4, 6, 5) - fGK( 1, 2, 6, 5, 3, 4) );
      double xiT    = xiGK( tHres, uHres, s3now, s4now);
      double xiU    = xiGK( uHres, tHres, s3now, s4now);
      double xjTU   = xjGK( tHres, uHres, s3now, s4now);

      //  Couplings of outgoing fermion from Z. Combine with kinematics.
      int idAbs     = process[i5].idAbs();
      double lfZ    = couplingsPtr->lf(idAbs);
      double rfZ    = couplingsPtr->rf(idAbs);
      wt            = lfZ*lfZ * fGK135 + rfZ*rfZ * fGK136;
      wtMax         = 4. * s3now * s4now * (lfZ*lfZ + rfZ*rfZ)
                    * (xiT + xiU - xjTU);

    // Decay distribution like in f fbar -> H^+- -> W+- Z0.
    } else {
      double p35  = 2. * process[i3].p() * process[i5].p();
      double p46  = 2. * process[i4].p() * process[i6].p();
      wt          = 16. * p35 * p46;
      wtMax       = sH2;
    }
  }

  // Angular weight in top decay by standard routine.
  else if (process[process[iResBeg].mother1()].idAbs() == 6)
    return weightTopDecay( process, iResBeg, iResEnd);

  // Done.
  return (wt / wtMax);

}


//==========================================================================

// Sigma1ffbar2Rhorizontal class.
// Cross section for f fbar' -> R^0 (f is a quark or lepton).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma1ffbar2Rhorizontal::initProc() {

  // Store R^0 mass and width for propagator.
  mRes     = particleDataPtr->m0(41);
  GammaRes = particleDataPtr->mWidth(41);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;
  thetaWRat = 1. / (12. * couplingsPtr->sin2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(41);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour.

void Sigma1ffbar2Rhorizontal::sigmaKin() {

  // Set up Breit-Wigner. Cross section for W+ and W- separately.
  double sigBW  = 12. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double preFac = alpEM * thetaWRat * mH;
  sigma0Pos     = preFac * sigBW * particlePtr->resWidthOpen(41, mH);
  sigma0Neg     = preFac * sigBW * particlePtr->resWidthOpen(-41, mH);

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence.

double Sigma1ffbar2Rhorizontal::sigmaHat() {

  // Check for allowed flavour combinations, one generation apart.
  if (id1 * id2 > 0 || abs(id1 + id2) != 2) return 0.;

  // Find whether R0 or R0bar. Colour factors.
  double sigma = (id1 + id2 > 0) ? sigma0Pos : sigma0Neg;
  if (abs(id1) < 7) sigma /= 3.;

  // Answer.
  return sigma;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2Rhorizontal::setIdColAcol() {

  // Outgoing R0 or R0bar.
  id3 = (id1 +id2 > 0) ? 41 : -41;
  setId( id1, id2, id3);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

} // end namespace Pythia8
