// SigmaLeftRightSym.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// left-right-symmetry simulation classes. 

#include "SigmaLeftRightSym.h"

namespace Pythia8 {

//==========================================================================

// Sigma1ffbar2ZRight class.
// Cross section for f fbar -> Z_R^0 (righthanded gauge boson). 

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma1ffbar2ZRight::initProc() {

  // Store Z_R mass and width for propagator. 
  idZR     = 9900023;
  mRes     = particleDataPtr->m0(idZR);
  GammaRes = particleDataPtr->mWidth(idZR);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;
  sin2tW   = couplingsPtr->sin2thetaW();

  // Set pointer to particle properties and decay table.
  ZRPtr    = particleDataPtr->particleDataEntryPtr(idZR);

} 

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour. 

void Sigma1ffbar2ZRight::sigmaKin() { 

  // Set up Breit-Wigner. Width out only includes open channels. 
  double sigBW    = 12. * M_PI/ ( pow2(sH - m2Res) + pow2(sH * GamMRat) );    
  double widthOut = ZRPtr->resWidthOpen(idZR, mH);

  // Prefactor for incoming widths. Combine. Done. 
  double preFac   = alpEM * mH / ( 48. * sin2tW * (1. - sin2tW) 
                  * (1. - 2. * sin2tW) ); 
  sigma0          = preFac * sigBW * widthOut;    

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence. 

double Sigma1ffbar2ZRight::sigmaHat() { 

  // Vector and axial couplings of incoming fermion pair.
  int    idAbs = abs(id1); 
  double af = 0.;
  double vf = 0.;
  if (idAbs < 9 && idAbs%2 == 1) {
    af      = -1. + 2. * sin2tW; 
    vf      = -1. + 4. * sin2tW / 3.;
  } else if (idAbs < 9) {   
    af      = 1. - 2. * sin2tW; 
    vf      = 1. - 8. * sin2tW / 3.;
  } else if (idAbs < 19 && idAbs%2 == 1) {   
    af      = -1. + 2. * sin2tW; 
    vf      = -1. + 4. * sin2tW;
  }

  // Colour factor. Answer.
  double sigma = (vf*vf + af*af) * sigma0;
  if (idAbs < 9) sigma /= 3.;
  return sigma;    

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2ZRight::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idZR);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for Z_R decay angle.
  
double Sigma1ffbar2ZRight::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);

  // Z_R should sit in entry 5.
  if (iResBeg != 5 || iResEnd != 5) return 1.;

  // Couplings for in- and out-flavours.
  double ai, vi, af, vf;
  int idInAbs   = process[3].idAbs();
  if (idInAbs < 9 && idInAbs%2 == 1) {
    ai          = -1. + 2. * sin2tW; 
    vi          = -1. + 4. * sin2tW / 3.;
  } else if (idInAbs < 9) {   
    ai          = 1. - 2. * sin2tW; 
    vi          = 1. - 8. * sin2tW / 3.;
  } else {   
    ai          = -1. + 2. * sin2tW; 
    vi          = -1. + 4. * sin2tW;
  }
  int idOutAbs = process[6].idAbs();
  if (idOutAbs < 9 && idOutAbs%2 == 1) {
    af          = -1. + 2. * sin2tW; 
    vf          = -1. + 4. * sin2tW / 3.;
  } else if (idOutAbs < 9) {   
    af          = 1. - 2. * sin2tW; 
    vf          = 1. - 8. * sin2tW / 3.;
  } else {   
    af          = -1. + 2. * sin2tW; 
    vf          = -1. + 4. * sin2tW;
  }

  // Phase space factors. Reconstruct decay angle.
  double mr1    = pow2(process[6].m()) / sH;
  double mr2    = pow2(process[7].m()) / sH;
  double betaf  = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2); 
  double cosThe = (process[3].p() - process[4].p()) 
    * (process[7].p() - process[6].p()) / (sH * betaf);

  // Angular weight and its maximum.
  double wt1    = (vi*vi + ai*ai) * (vf*vf + af*af * betaf*betaf);
  double wt2    = (1. - betaf*betaf) * (vi*vi + ai*ai) * vf*vf;
  double wt3    = betaf * 4. * vi * ai * vf * af;  
  if (process[3].id() * process[6].id() < 0) wt3 = -wt3;
  double wt     = wt1 * (1. + cosThe*cosThe) + wt2 * (1. - cosThe*cosThe)
                + 2. * wt3 * cosThe;
  double wtMax  = 2. * (wt1 + abs(wt3));

  // Done.
  return wt / wtMax;

}

//==========================================================================

// Sigma1ffbar2WRight class.
// Cross section for f fbar' -> W_R^+- (righthanded gauge boson). 

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma1ffbar2WRight::initProc() {

  // Store W_R^+- mass and width for propagator. 
  idWR     = 9900024;
  mRes     = particleDataPtr->m0(idWR);
  GammaRes = particleDataPtr->mWidth(idWR);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;
  thetaWRat = 1. / (12. * couplingsPtr->sin2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(idWR);

} 

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour. 

void Sigma1ffbar2WRight::sigmaKin() { 

  // Common coupling factors.
  double colQ   = 3. * (1. + alpS / M_PI);

  // Reset quantities to sum. Declare variables inside loop.
  double widOutPos = 0.; 
  double widOutNeg = 0.; 
  int    id1Now, id2Now, id1Abs, id2Abs, id1Neg, id2Neg, onMode;
  double widNow, widSecPos, widSecNeg, mf1, mf2, mr1, mr2, kinFac;

  // Loop over all W_R^+- decay channels. 
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    id1Now      = particlePtr->channel(i).product(0);
    id2Now      = particlePtr->channel(i).product(1);
    id1Abs      = abs(id1Now);
    id2Abs      = abs(id2Now);

    // Check that above threshold. Phase space.
    mf1 = particleDataPtr->m0(id1Abs);
    mf2 = particleDataPtr->m0(id2Abs);
    if (mH > mf1 + mf2 + MASSMARGIN) {
      mr1    = pow2(mf1 / mH);
      mr2    = pow2(mf2 / mH);
      kinFac = (1. - 0.5 * (mr1 + mr2) - 0.5 * pow2(mr1 - mr2))
             * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

      // Combine kinematics with colour factor and CKM couplings.
      widNow = kinFac;
      if (id1Abs < 9) widNow *= colQ * couplingsPtr->V2CKMid(id1Abs, id2Abs);
 
      // Secondary width from top and righthanded neutrino decay.
      id1Neg    = (id1Abs < 19) ? -id1Now : id1Abs; 
      id2Neg    = (id2Abs < 19) ? -id2Now : id2Abs; 
      widSecPos = particleDataPtr->resOpenFrac(id1Now, id2Now); 
      widSecNeg = particleDataPtr->resOpenFrac(id1Neg, id2Neg); 

      // Add weight for channels on for all, W_R^+ and W_R^-, respectively.
      onMode = particlePtr->channel(i).onMode();
      if (onMode == 1 || onMode == 2) widOutPos += widNow * widSecPos;
      if (onMode == 1 || onMode == 3) widOutNeg += widNow * widSecNeg;

    // End loop over fermions.
    }
  }

  // Set up Breit-Wigner. Cross section for W_R^+ and W_R^- separately.
  double sigBW = 12. * M_PI * pow2(alpEM * thetaWRat) * sH
               / ( pow2(sH - m2Res) + pow2(sH * GamMRat) ); 
  sigma0Pos = sigBW * widOutPos;    
  sigma0Neg = sigBW * widOutNeg;    

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence. 

double Sigma1ffbar2WRight::sigmaHat() {

  // Secondary width for W_R^+ or W_R^-. CKM and colour factors.
  int idUp = (abs(id1)%2 == 0) ? id1 : id2;
  double sigma = (idUp > 0) ? sigma0Pos : sigma0Neg;
  if (abs(id1) < 9) sigma *= couplingsPtr->V2CKMid(abs(id1), abs(id2)) / 3.;

  // Answer.
  return sigma;    

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ffbar2WRight::setIdColAcol() {

  // Sign of outgoing W_R.
  int sign          = (abs(id1)%2 == 0) ? 1 : -1;
  if (id1 < 0) sign = -sign;
  setId( id1, id2, idWR * sign);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for W_R decay angle.
  
double Sigma1ffbar2WRight::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);

  // W_R should sit in entry 5.
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

// Sigma1ll2Hchgchg class.
// Cross section for l l -> H_L^++-- or H_R^++-- (doubly charged Higgs).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma1ll2Hchgchg::initProc() {

  // Set process properties: H_L^++-- or H_R^++--.
  if (leftRight == 1) {
    idHLR    = 9900041;
    codeSave = 3121;
    nameSave = "l l -> H_L^++--";
  } else {
    idHLR    = 9900042;
    codeSave = 3141;
    nameSave = "l l -> H_R^++--";
  }

  // Read in Yukawa matrix for couplings to a lepton pair.
  yukawa[1][1]  = settingsPtr->parm("LeftRightSymmmetry:coupHee");
  yukawa[2][1]  = settingsPtr->parm("LeftRightSymmmetry:coupHmue");
  yukawa[2][2]  = settingsPtr->parm("LeftRightSymmmetry:coupHmumu");
  yukawa[3][1]  = settingsPtr->parm("LeftRightSymmmetry:coupHtaue");
  yukawa[3][2]  = settingsPtr->parm("LeftRightSymmmetry:coupHtaumu");
  yukawa[3][3]  = settingsPtr->parm("LeftRightSymmmetry:coupHtautau");

  // Store H_L/R mass and width for propagator. 
  mRes     = particleDataPtr->m0(idHLR);
  GammaRes = particleDataPtr->mWidth(idHLR);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Set pointer to particle properties and decay table.
  particlePtr = particleDataPtr->particleDataEntryPtr(idHLR);

} 

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence. 

double Sigma1ll2Hchgchg::sigmaHat() {

  // Initial state must consist of two identical-sign leptons.
  if (id1 * id2 < 0) return 0.;
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);  
  if (id1Abs != 11 && id1Abs != 13 && id1Abs != 15) return 0.;
  if (id2Abs != 11 && id2Abs != 13 && id2Abs != 15) return 0.;

  // Set up Breit-Wigner, inwidth and outwidth.
  double sigBW  = 8. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) ); 
  double widIn  = pow2(yukawa[(id1Abs-9)/2][(id2Abs-9)/2]) 
                * mH / (8. * M_PI);
  int idSgn     = (id1 < 0) ? idHLR : -idHLR;
  double widOut = particlePtr->resWidthOpen( idSgn, mH);

  // Answer.
  return widIn * sigBW * widOut;    

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma1ll2Hchgchg::setIdColAcol() {

  // Sign of outgoing H_L/R.
  int idSgn     = (id1 < 0) ? idHLR : -idHLR;
  setId( id1, id2, idSgn);

  // No colours whatsoever.
  setColAcol( 0, 0, 0, 0, 0, 0);

}

//--------------------------------------------------------------------------

// Evaluate weight for H_L/R sequential decay angles.
  
double Sigma1ll2Hchgchg::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);
 
  // Else isotropic decay.
  return 1.;

}

//==========================================================================

// Sigma2lgm2Hchgchgl class.
// Cross section for l gamma -> H_L^++-- l or H_R^++-- l 
// (doubly charged Higgs).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2lgm2Hchgchgl::initProc() {

  // Set process properties: H_L^++-- or H_R^++-- and e/mu/tau.
  idHLR        = (leftRight == 1) ? 9900041 : 9900042;
  codeSave     = (leftRight == 1) ? 3122 : 3142;
  if (idLep == 13) codeSave += 2;
  if (idLep == 15) codeSave += 4;
  if      (codeSave == 3122) nameSave = "l^+- gamma -> H_L^++-- e^-+";
  else if (codeSave == 3123) nameSave = "l^+- gamma -> H_L^++-- mu^-+";
  else if (codeSave == 3124) nameSave = "l^+- gamma -> H_L^++-- tau^-+";
  else if (codeSave == 3142) nameSave = "l^+- gamma -> H_R^++-- e^-+";
  else if (codeSave == 3143) nameSave = "l^+- gamma -> H_R^++-- mu^-+";
  else                       nameSave = "l^+- gamma -> H_R^++-- tau^-+";

  // Read in relevantYukawa matrix for couplings to a lepton pair.
  if (idLep == 11) {
    yukawa[1]  = settingsPtr->parm("LeftRightSymmmetry:coupHee");
    yukawa[2]  = settingsPtr->parm("LeftRightSymmmetry:coupHmue");
    yukawa[3]  = settingsPtr->parm("LeftRightSymmmetry:coupHtaue");
  } else if (idLep == 13) { 
    yukawa[1]  = settingsPtr->parm("LeftRightSymmmetry:coupHmue");
    yukawa[2]  = settingsPtr->parm("LeftRightSymmmetry:coupHmumu");
    yukawa[3]  = settingsPtr->parm("LeftRightSymmmetry:coupHtaumu");
  } else {
    yukawa[1]  = settingsPtr->parm("LeftRightSymmmetry:coupHtaue");
    yukawa[2]  = settingsPtr->parm("LeftRightSymmmetry:coupHtaumu");
    yukawa[3]  = settingsPtr->parm("LeftRightSymmmetry:coupHtautau");
  }

  // Secondary open width fractions.
  openFracPos  = particleDataPtr->resOpenFrac( idHLR);
  openFracNeg  = particleDataPtr->resOpenFrac(-idHLR);

} 

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence. 

double Sigma2lgm2Hchgchgl::sigmaHat() {

  // Initial state must consist of a lepton and a photon.
  int idIn     = (id2 == 22) ? id1 : id2;  
  int idInAbs  = abs(idIn);
  if (idInAbs != 11 && idInAbs != 13 && idInAbs != 15) return 0.;

  // Incoming squared lepton mass.
  double s1    = pow2( particleDataPtr->m0(idInAbs) );
  
  // Kinematical expressions.
  double smm1  = 8. * (sH + tH - s3) * (sH + tH - 2. * s3 - s1 - s4) 
               / pow2(uH - s3);
  double smm2  = 2. * ( (2. * s3 - 3. * s1) * s4 + (s1 - 2. * s4) * tH
               - (tH - s4) * sH ) / pow2(tH - s4);
  double smm3  = 2. * ( (2. * s3 - 3. * s4 + tH) * s1 
               - (2. * s1 - s4 + tH) * sH ) / pow2(sH - s1);
  double smm12 = 4. * ( (2. * s1 - s4 - 2. * s3 + tH) * sH 
               + (tH - 3. * s3 - 3. * s4) * tH + (2. * s3 - 2. * s1 
               + 3. * s4) * s3 ) / ( (uH - s3) * (tH - s4) );
  double smm13 = -4. * ( (tH + s1 - 2. * s4) * tH - (s3 + 3. * s1 - 2. * s4)
               * s3 + (s3 + 3. * s1 + tH) * sH - pow2(tH - s3 + sH) )
               / ( (uH - s3) * (sH - s1) );
  double smm23 = -4. * ( (s1 - s4 + s3) * tH - s3*s3 + s3 * (s1 + s4)
               - 3. * s1 * s4 - (s1 - s4 - s3 + tH) * sH)
               / ( (sH - s1) * (tH - s4) );
  double sigma = alpEM * pow2(sH / (sH - s1) ) * (smm1 + smm2 + smm3
               + smm12 + smm13 + smm23) / (4. * sH2);

  // Lepton Yukawa and secondary widths.
  sigma       *= pow2(yukawa[(idInAbs-9)/2]);
  sigma       *= (idIn < 0) ? openFracPos : openFracNeg; 

  // Answer.
  return sigma;    

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2lgm2Hchgchgl::setIdColAcol() {

  // Sign of outgoing H_L/R.
  int idIn     = (id2 == 22) ? id1 : id2;
  int idSgn    = (idIn < 0) ? idHLR : -idHLR;
  int idOut    = (idIn < 0) ? idLep : -idLep;
  setId( id1, id2, idSgn, idOut);

  // tHat is defined between incoming lepton and outgoing Higgs.
  if (id1 == 22) swapTU = true;

  // No colours whatsoever.
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//--------------------------------------------------------------------------

// Evaluate weight for H_L/R sequential decay angles.
  
double Sigma2lgm2Hchgchgl::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);
 
  // Else isotropic decay.
  return 1.;

}

//==========================================================================

// Sigma3ff2HchgchgfftWW class.
// Cross section for  f_1 f_2 -> H_(L/R)^++-- f_3 f_4 (W+- W+- fusion).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma3ff2HchgchgfftWW::initProc() {

  // Set process properties: H_L^++-- or H_R^++--.
  if (leftRight == 1) {
    idHLR      = 9900041;
    codeSave   = 3125;
    nameSave   = "f_1 f_2 -> H_L^++-- f_3 f_4 (W+- W+- fusion)";
  } else {
    idHLR      = 9900042;
    codeSave   = 3145;
    nameSave   = "f_1 f_2 -> H_R^++-- f_3 f_4 (W+- W+- fusion)";
  }

  // Common fixed mass and coupling factor.
  double mW    = particleDataPtr->m0(24); 
  double mWR   = particleDataPtr->m0(9900024);
  mWS          = (leftRight == 1) ? pow2(mW) : pow2(mWR);
  double gL    = settingsPtr->parm("LeftRightSymmmetry:gL");
  double gR    = settingsPtr->parm("LeftRightSymmmetry:gR");
  double vL    = settingsPtr->parm("LeftRightSymmmetry:vL"); 
  prefac       = (leftRight == 1) ? pow2(pow4(gL) * vL) 
                                  : 2. * pow2(pow3(gR) * mWR); 
  // Secondary open width fractions.
  openFracPos  = particleDataPtr->resOpenFrac( idHLR);
  openFracNeg  = particleDataPtr->resOpenFrac(-idHLR);

} 

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), part independent of incoming flavour. 

void Sigma3ff2HchgchgfftWW::sigmaKin() { 

  // Required four-vector products.
  double pp12  = 0.5 * sH;
  double pp14  = 0.5 * mH * p4cm.pNeg();
  double pp15  = 0.5 * mH * p5cm.pNeg();
  double pp24  = 0.5 * mH * p4cm.pPos();
  double pp25  = 0.5 * mH * p5cm.pPos();
  double pp45  = p4cm * p5cm;

  // Cross section: kinematics part. Combine with couplings.
  double propT = 1. / ( (2. * pp14 + mWS) * (2. * pp25 + mWS) );  
  double propU = 1. / ( (2. * pp24 + mWS) * (2. * pp15 + mWS) );
  sigma0TU     = prefac * pp12 * pp45 * pow2(propT + propU); 
  sigma0T      = prefac * pp12 * pp45 * 2. * pow2(propT); 

}

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence. 

double Sigma3ff2HchgchgfftWW::sigmaHat() { 

  // Do not allow creation of righthanded neutrinos for H_R.
  int id1Abs   = abs(id1);
  int id2Abs   = abs(id2);
  if ( leftRight == 2 && (id1Abs > 10 || id2Abs > 10) ) return 0.; 

  // Many flavour combinations not possible because of charge.
  int chg1     = (( id1Abs%2 == 0 && id1 > 0) 
               || (id1Abs%2 == 1 && id1 < 0) ) ? 1 : -1;
  int chg2     = (( id2Abs%2 == 0 && id2 > 0) 
               || (id2Abs%2 == 1 && id2 < 0) ) ? 1 : -1;
  if (abs(chg1 + chg2) != 2) return 0.;

  // Basic cross section. CKM factors for final states.
  double sigma = (id2 == id1 && id1Abs > 10) ? sigma0TU : sigma0T;
  sigma       *= couplingsPtr->V2CKMsum(id1Abs) 
               * couplingsPtr->V2CKMsum(id2Abs);

  // Secondary width for H0.
  sigma       *= (chg1 + chg2 == 2) ? openFracPos : openFracNeg;

  // Spin-state extra factor 2 per incoming neutrino.
  if (id1Abs == 12 || id1Abs == 14 || id1Abs == 16) sigma *= 2.; 
  if (id2Abs == 12 || id2Abs == 14 || id2Abs == 16) sigma *= 2.; 

  // Answer.
  return sigma;
 
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3ff2HchgchgfftWW::setIdColAcol() {

  // Pick out-flavours by relative CKM weights.
  int id1Abs   = abs(id1);
  int id2Abs   = abs(id2);
  id4          = couplingsPtr->V2CKMpick(id1);
  id5          = couplingsPtr->V2CKMpick(id2);

  // Find charge of Higgs .
  id3 = (( id1Abs%2 == 0 && id1 > 0) || (id1Abs%2 == 1 && id1 < 0) ) 
      ? idHLR : -idHLR;
  setId( id1, id2, id3, id4, id5);

  // Colour flow topologies. Swap when antiquarks.
  if (id1Abs < 9 && id2Abs < 9 && id1*id2 > 0) 
                       setColAcol( 1, 0, 2, 0, 0, 0, 1, 0, 2, 0);
  else if (id1Abs < 9 && id2Abs < 9)
                       setColAcol( 1, 0, 0, 2, 0, 0, 1, 0, 0, 2); 
  else if (id1Abs < 9) setColAcol( 1, 0, 0, 0, 0, 0, 1, 0, 0, 0); 
  else if (id2Abs < 9) setColAcol( 0, 0, 1, 0, 0, 0, 0, 0, 1, 0); 
  else                 setColAcol( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  if ( (id1Abs < 9 && id1 < 0) || (id1Abs > 10 && id2 < 0) ) 
    swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles.

double Sigma3ff2HchgchgfftWW::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);

  // Else done.
  return 1.; 

}

//==========================================================================

// Sigma2ffbar2HchgchgHchgchg class.
// Cross section for f fbar -> H_(L/R)^++ H_(L/R)^-- (doubly charged Higgs).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2ffbar2HchgchgHchgchg::initProc() {

  // Set process properties: H_L^++ H_L^-- or H_R^++ H_R^--.
  if (leftRight == 1) {
    idHLR      = 9900041;
    codeSave   = 3126;
    nameSave   = "f fbar -> H_L^++ H_L^--";
  } else {
    idHLR      = 9900042;
    codeSave   = 3146;
    nameSave   = "f fbar -> H_R^++ H_R^--";
  }

  // Read in Yukawa matrix for couplings to a lepton pair.
  yukawa[1][1] = settingsPtr->parm("LeftRightSymmmetry:coupHee");
  yukawa[2][1] = settingsPtr->parm("LeftRightSymmmetry:coupHmue");
  yukawa[2][2] = settingsPtr->parm("LeftRightSymmmetry:coupHmumu");
  yukawa[3][1] = settingsPtr->parm("LeftRightSymmmetry:coupHtaue");
  yukawa[3][2] = settingsPtr->parm("LeftRightSymmmetry:coupHtaumu");
  yukawa[3][3] = settingsPtr->parm("LeftRightSymmmetry:coupHtautau");

  // Electroweak parameters.
  mRes         = particleDataPtr->m0(23);
  GammaRes     = particleDataPtr->mWidth(23);
  m2Res        = mRes*mRes;
  GamMRat      = GammaRes / mRes;
  sin2tW       = couplingsPtr->sin2thetaW();
  preFac       = (1. - 2. * sin2tW) / ( 8. * sin2tW * (1. - sin2tW) );

  // Open fraction from secondary widths.
  openFrac     = particleDataPtr->resOpenFrac( idHLR, -idHLR);

} 

//--------------------------------------------------------------------------

// Evaluate sigmaHat(sHat), including incoming flavour dependence. 

double Sigma2ffbar2HchgchgHchgchg::sigmaHat() {

  // Electroweak couplings to gamma^*/Z^0.
  int    idAbs   = abs(id1);
  double ei      = couplingsPtr->ef(idAbs);
  double vi      = couplingsPtr->vf(idAbs); 
  double ai      = couplingsPtr->af(idAbs); 

  // Part via gamma^*/Z^0 propagator. No Z^0 coupling to H_R.
  double resProp = 1. / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double sigma   = 8. * pow2(alpEM) * ei*ei / sH2;
  if (leftRight == 1) sigma += 8. * pow2(alpEM) 
    * (2. * ei * vi * preFac * (sH - m2Res) * resProp / sH
    + (vi * vi + ai * ai) * pow2(preFac) * resProp);      

  // Part via t-channel lepton + interference; sum over possibilities. 
  if (idAbs == 11 || idAbs == 13 || idAbs == 15) {
    double yuk2Sum;
    if (idAbs == 11) yuk2Sum 
      = pow2(yukawa[1][1]) + pow2(yukawa[2][1]) + pow2(yukawa[3][1]);
    else if (idAbs == 13) yuk2Sum 
      = pow2(yukawa[2][1]) + pow2(yukawa[2][2]) + pow2(yukawa[3][2]);
    else yuk2Sum 
      = pow2(yukawa[3][1]) + pow2(yukawa[3][2]) + pow2(yukawa[3][3]);
    yuk2Sum /= 4. * M_PI;
    sigma   += 8. * alpEM * ei * yuk2Sum / (sH * tH) 
      + 4. * pow2(yuk2Sum) / tH2;
    if (leftRight == 1) sigma += 8. * alpEM * (vi + ai) * yuk2Sum 
      * preFac * (sH - m2Res) * resProp / tH;   
  }  

  // Common kinematical factor. Colour factor.
  sigma *= M_PI * (tH * uH - s3 * s4) / sH2;
  if (idAbs < 9) sigma /= 3.;  

  // Answer.
  return sigma;    

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2ffbar2HchgchgHchgchg::setIdColAcol() {

  // Outgoing flavours trivial.
  setId( id1, id2, idHLR, -idHLR);

  // tHat is defined between incoming fermion and outgoing H--.
  if (id1 > 0) swapTU = true;

  // No colours at all or one flow topology. Swap if first is antiquark.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for H_L/R sequential decay angles.
  
double Sigma2ffbar2HchgchgHchgchg::weightDecay( Event& process, 
  int iResBeg, int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);
 
  // Else isotropic decay.
  return 1.;

}

//==========================================================================

} // end namespace Pythia8
