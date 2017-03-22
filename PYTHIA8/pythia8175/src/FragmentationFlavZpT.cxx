// FragmentationFlavZpT.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// StringFlav, StringZ and StringPT classes.

#include "FragmentationFlavZpT.h"

namespace Pythia8 {

//==========================================================================

// The StringFlav class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Offset for different meson multiplet id values.
const int StringFlav::mesonMultipletCode[6] 
  = { 1, 3, 10003, 10001, 20003, 5};

// Clebsch-Gordan coefficients for baryon octet and decuplet are
// fixed once and for all, so only weighted sum needs to be edited.
// Order: ud0 + u, ud0 + s, uu1 + u, uu1 + d, ud1 + u, ud1 + s.
const double StringFlav::baryonCGOct[6] 
  = { 0.75, 0.5, 0., 0.1667, 0.0833, 0.1667};
const double StringFlav::baryonCGDec[6] 
  = { 0.,  0.,  1., 0.3333, 0.6667, 0.3333};

//--------------------------------------------------------------------------

// Initialize data members of the flavour generation.

void StringFlav::init(Settings& settings, Rndm* rndmPtrIn) {

  // Save pointer.
  rndmPtr         = rndmPtrIn;

  // Basic parameters for generation of new flavour.
  probQQtoQ       = settings.parm("StringFlav:probQQtoQ");
  probStoUD       = settings.parm("StringFlav:probStoUD");
  probSQtoQQ      = settings.parm("StringFlav:probSQtoQQ");
  probQQ1toQQ0    = settings.parm("StringFlav:probQQ1toQQ0");

  // Parameters derived from above.
  probQandQQ      = 1. + probQQtoQ;
  probQandS       = 2. + probStoUD;
  probQandSinQQ   = 2. + probSQtoQQ * probStoUD;
  probQQ1corr     = 3. * probQQ1toQQ0;
  probQQ1corrInv  = 1. / probQQ1corr;
  probQQ1norm     = probQQ1corr / (1. + probQQ1corr);

  // Parameters for normal meson production.
  for (int i = 0; i < 4; ++i) mesonRate[i][0] = 1.;
  mesonRate[0][1] = settings.parm("StringFlav:mesonUDvector");
  mesonRate[1][1] = settings.parm("StringFlav:mesonSvector");
  mesonRate[2][1] = settings.parm("StringFlav:mesonCvector");
  mesonRate[3][1] = settings.parm("StringFlav:mesonBvector");

  // Parameters for L=1 excited-meson production.
  mesonRate[0][2] = settings.parm("StringFlav:mesonUDL1S0J1");
  mesonRate[1][2] = settings.parm("StringFlav:mesonSL1S0J1");
  mesonRate[2][2] = settings.parm("StringFlav:mesonCL1S0J1");
  mesonRate[3][2] = settings.parm("StringFlav:mesonBL1S0J1");
  mesonRate[0][3] = settings.parm("StringFlav:mesonUDL1S1J0");
  mesonRate[1][3] = settings.parm("StringFlav:mesonSL1S1J0");
  mesonRate[2][3] = settings.parm("StringFlav:mesonCL1S1J0");
  mesonRate[3][3] = settings.parm("StringFlav:mesonBL1S1J0");
  mesonRate[0][4] = settings.parm("StringFlav:mesonUDL1S1J1");
  mesonRate[1][4] = settings.parm("StringFlav:mesonSL1S1J1");
  mesonRate[2][4] = settings.parm("StringFlav:mesonCL1S1J1");
  mesonRate[3][4] = settings.parm("StringFlav:mesonBL1S1J1");
  mesonRate[0][5] = settings.parm("StringFlav:mesonUDL1S1J2");
  mesonRate[1][5] = settings.parm("StringFlav:mesonSL1S1J2");
  mesonRate[2][5] = settings.parm("StringFlav:mesonCL1S1J2");
  mesonRate[3][5] = settings.parm("StringFlav:mesonBL1S1J2");

  // Store sum over multiplets for Monte Carlo generation.
  for (int i = 0; i < 4; ++i) mesonRateSum[i] 
    = mesonRate[i][0] + mesonRate[i][1] + mesonRate[i][2] 
    + mesonRate[i][3] + mesonRate[i][4] + mesonRate[i][5];

  // Parameters for uubar - ddbar - ssbar meson mixing.
  for (int spin = 0; spin < 6; ++spin) { 
    double theta;
    if      (spin == 0) theta = settings.parm("StringFlav:thetaPS");
    else if (spin == 1) theta = settings.parm("StringFlav:thetaV");
    else if (spin == 2) theta = settings.parm("StringFlav:thetaL1S0J1");
    else if (spin == 3) theta = settings.parm("StringFlav:thetaL1S1J0");
    else if (spin == 4) theta = settings.parm("StringFlav:thetaL1S1J1");
    else                theta = settings.parm("StringFlav:thetaL1S1J2");
    double alpha = (spin == 0) ? 90. - (theta + 54.7) : theta + 54.7;
    alpha *= M_PI / 180.;
    // Fill in (flavour, spin)-dependent probability of producing
    // the lightest or the lightest two mesons of the nonet. 
    mesonMix1[0][spin] = 0.5;
    mesonMix2[0][spin] = 0.5 * (1. + pow2(sin(alpha)));
    mesonMix1[1][spin] = 0.;
    mesonMix2[1][spin] = pow2(cos(alpha));
  }

  // Additional suppression of eta and etaPrime.
  etaSup      = settings.parm("StringFlav:etaSup");
  etaPrimeSup = settings.parm("StringFlav:etaPrimeSup");

  // Sum of baryon octet and decuplet weights.
  decupletSup = settings.parm("StringFlav:decupletSup");
  for (int i = 0; i < 6; ++i) baryonCGSum[i]
    = baryonCGOct[i] + decupletSup * baryonCGDec[i];

  // Maximum SU(6) weight for ud0, ud1, uu1 types.
  baryonCGMax[0] = max( baryonCGSum[0], baryonCGSum[1]);
  baryonCGMax[1] = baryonCGMax[0]; 
  baryonCGMax[2] = max( baryonCGSum[2], baryonCGSum[3]);
  baryonCGMax[3] = baryonCGMax[2]; 
  baryonCGMax[4] = max( baryonCGSum[4], baryonCGSum[5]);
  baryonCGMax[5] = baryonCGMax[4]; 

  // Popcorn baryon parameters.
  popcornRate    = settings.parm("StringFlav:popcornRate");
  popcornSpair   = settings.parm("StringFlav:popcornSpair"); 
  popcornSmeson  = settings.parm("StringFlav:popcornSmeson");
  
  // Suppression of leading (= first-rank) baryons.
  suppressLeadingB = settings.flag("StringFlav:suppressLeadingB");
  lightLeadingBSup = settings.parm("StringFlav:lightLeadingBSup");
  heavyLeadingBSup = settings.parm("StringFlav:heavyLeadingBSup");

  // Begin calculation of derived parameters for baryon production.

  // Enumerate distinguishable diquark types (in diquark first is popcorn q).
  enum Diquark {ud0, ud1, uu1, us0, su0, us1, su1, ss1};

  // Maximum SU(6) weight by diquark type. 
  double barCGMax[8];
  barCGMax[ud0] = baryonCGMax[0]; 
  barCGMax[ud1] = baryonCGMax[4]; 
  barCGMax[uu1] = baryonCGMax[2]; 
  barCGMax[us0] = baryonCGMax[0]; 
  barCGMax[su0] = baryonCGMax[0]; 
  barCGMax[us1] = baryonCGMax[4]; 
  barCGMax[su1] = baryonCGMax[4]; 
  barCGMax[ss1] = baryonCGMax[2]; 

  // Diquark SU(6) survival = Sum_quark (quark tunnel weight) * SU(6).
  double dMB[8];
  dMB[ud0] = 2. * baryonCGSum[0] + probStoUD * baryonCGSum[1];
  dMB[ud1] = 2. * baryonCGSum[4] + probStoUD * baryonCGSum[5];
  dMB[uu1] = baryonCGSum[2] + (1. + probStoUD) * baryonCGSum[3];
  dMB[us0] = (1. + probStoUD) * baryonCGSum[0] + baryonCGSum[1];
  dMB[su0] = dMB[us0];
  dMB[us1] = (1. + probStoUD) * baryonCGSum[4] + baryonCGSum[5];  
  dMB[su1] = dMB[us1];
  dMB[ss1] = probStoUD * baryonCGSum[2] + 2. * baryonCGSum[3];
  for (int i = 1; i < 8; ++i) dMB[i] = dMB[i] / dMB[0];

  // Tunneling factors for diquark production; only half a pair = sqrt.
  double probStoUDroot    = sqrt(probStoUD);
  double probSQtoQQroot   = sqrt(probSQtoQQ);
  double probQQ1toQQ0root = sqrt(probQQ1toQQ0);
  double qBB[8];
  qBB[ud1] = probQQ1toQQ0root;
  qBB[uu1] = probQQ1toQQ0root;
  qBB[us0] = probSQtoQQroot;
  qBB[su0] = probStoUDroot * probSQtoQQroot;
  qBB[us1] = probQQ1toQQ0root * qBB[us0];
  qBB[su1] = probQQ1toQQ0root * qBB[su0];
  qBB[ss1] = probStoUDroot * pow2(probSQtoQQroot) * probQQ1toQQ0root;

  // spin * (vertex factor) * (half-tunneling factor above).
  double qBM[8];
  qBM[ud1] = 3. * qBB[ud1];
  qBM[uu1] = 6. * qBB[uu1];
  qBM[us0] = probStoUD * qBB[us0];
  qBM[su0] = qBB[su0]; 
  qBM[us1] = probStoUD * 3. * qBB[us1];
  qBM[su1] = 3. * qBB[su1];
  qBM[ss1] = probStoUD * 6. * qBB[ss1];

  // Combine above two into total diquark weight for q -> B Bbar.
  for (int i = 1; i < 8; ++i) qBB[i] = qBB[i] * qBM[i];

  // Suppression from having strange popcorn meson.
  qBM[us0] *= popcornSmeson;
  qBM[us1] *= popcornSmeson;
  qBM[ss1] *= popcornSmeson;

  // Suppression for a heavy quark of a diquark to fit into a baryon
  // on the other side of popcorn meson: (0) s/u for q -> B M;
  // (1) s/u for rank 0 diquark su -> M B; (2) ditto for s -> c/b.
  double uNorm = 1. + qBM[ud1] + qBM[uu1] + qBM[us0] + qBM[us1];
  scbBM[0] = (2. * (qBM[su0] + qBM[su1]) + qBM[ss1]) / uNorm;  
  scbBM[1] = scbBM[0] * popcornSpair * qBM[su0] / qBM[us0];
  scbBM[2] = (1. + qBM[ud1]) * (2. + qBM[us0]) / uNorm; 

  // Include maximum of Clebsch-Gordan coefficients.
  for (int i = 1; i < 8; ++i) dMB[i] *= qBM[i];
  for (int i = 1; i < 8; ++i) qBM[i] *= barCGMax[i] / barCGMax[0];
  for (int i = 1; i < 8; ++i) qBB[i] *= barCGMax[i] / barCGMax[0];

  // Popcorn fraction for normal diquark production.
  double qNorm = uNorm * popcornRate / 3.;
  double sNorm = scbBM[0] * popcornSpair;
  popFrac = qNorm * (1. + qBM[ud1] + qBM[uu1] + qBM[us0] + qBM[us1] 
    + sNorm * (qBM[su0] + qBM[su1] + 0.5 * qBM[ss1])) / (1. +  qBB[ud1] 
    + qBB[uu1] + 2. * (qBB[us0] + qBB[us1]) + 0.5 * qBB[ss1]);  

  // Popcorn fraction for rank 0 diquarks, depending on number of s quarks.
  popS[0] = qNorm * qBM[ud1] / qBB[ud1];
  popS[1] = qNorm * 0.5 * (qBM[us1] / qBB[us1] 
    + sNorm * qBM[su1] / qBB[su1]);
  popS[2] = qNorm * sNorm * qBM[ss1] / qBB[ss1]; 

  // Recombine diquark weights to flavour and spin ratios. Second index: 
  // 0 = s/u popcorn quark ratio.
  // 1, 2 = s/u ratio for vertex quark if popcorn quark is u/d or s.
  // 3 = q/q' vertex quark ratio if popcorn quark is light and = q.
  // 4, 5, 6 = (spin 1)/(spin 0) ratio for su, us and ud.  

  // Case 0: q -> B B. 
  dWT[0][0] = (2. * (qBB[su0] + qBB[su1]) + qBB[ss1])
    / (1. + qBB[ud1] + qBB[uu1] + qBB[us0] + qBB[us1]);
  dWT[0][1] = 2. * (qBB[us0] + qBB[us1]) / (1. + qBB[ud1] + qBB[uu1]);
  dWT[0][2] = qBB[ss1] / (qBB[su0] + qBB[su1]);
  dWT[0][3] = qBB[uu1] / (1. + qBB[ud1] + qBB[uu1]);
  dWT[0][4] = qBB[su1] / qBB[su0];
  dWT[0][5] = qBB[us1] / qBB[us0];
  dWT[0][6] = qBB[ud1];

  // Case 1: q -> B M B.
  dWT[1][0] = (2. * (qBM[su0] + qBM[su1]) + qBM[ss1])
    / (1. + qBM[ud1] + qBM[uu1] + qBM[us0] + qBM[us1]);
  dWT[1][1] = 2. * (qBM[us0] + qBM[us1]) / (1. + qBM[ud1] + qBM[uu1]);
  dWT[1][2] = qBM[ss1] / (qBM[su0] + qBM[su1]);
  dWT[1][3] = qBM[uu1] / (1. + qBM[ud1] + qBM[uu1]);
  dWT[1][4] = qBM[su1] / qBM[su0];
  dWT[1][5] = qBM[us1] / qBM[us0];
  dWT[1][6] = qBM[ud1];
  
  // Case 2: qq -> M B; diquark inside chain.
  dWT[2][0] = (2. * (dMB[su0] + dMB[su1]) + dMB[ss1])
    / (1. + dMB[ud1] + dMB[uu1] + dMB[us0] + dMB[us1]);
  dWT[2][1] = 2. * (dMB[us0] + dMB[us1]) / (1. + dMB[ud1] + dMB[uu1]);
  dWT[2][2] = dMB[ss1] / (dMB[su0] + dMB[su1]);
  dWT[2][3] = dMB[uu1] / (1. + dMB[ud1] + dMB[uu1]);
  dWT[2][4] = dMB[su1] / dMB[su0];
  dWT[2][5] = dMB[us1] / dMB[us0];
  dWT[2][6] = dMB[ud1];

}

//--------------------------------------------------------------------------

// Pick a new flavour (including diquarks) given an incoming one.

FlavContainer StringFlav::pick(FlavContainer& flavOld) {

  // Initial values for new flavour.
  FlavContainer flavNew;
  flavNew.rank = flavOld.rank + 1;

  // For original diquark assign popcorn quark and whether popcorn meson.
  int idOld = abs(flavOld.id);
  if (flavOld.rank == 0 && idOld > 1000) assignPopQ(flavOld);

  // Diquark exists, to be forced into baryon now.
  bool doOldBaryon    = (idOld > 1000 && flavOld.nPop == 0);
  // Diquark exists, but do meson now.
  bool doPopcornMeson = flavOld.nPop > 0;
  // Newly created diquark gives baryon now, antibaryon later.
  bool doNewBaryon    = false;

  // Choose whether to generate a new meson or a new baryon.
  if (!doOldBaryon && !doPopcornMeson && probQandQQ * rndmPtr->flat() > 1.) {
    doNewBaryon = true;
    if ((1. + popFrac) * rndmPtr->flat() > 1.) flavNew.nPop = 1;
  }

  // Optional suppression of first-rank baryon.
  if (flavOld.rank == 0 && doNewBaryon && suppressLeadingB) {
    double leadingBSup = (idOld < 4) ? lightLeadingBSup : heavyLeadingBSup; 
    if (rndmPtr->flat() > leadingBSup) {
      doNewBaryon = false;
      flavNew.nPop = 0;
    }
  }

  // Single quark for new meson or for baryon where diquark already exists.
  if (!doPopcornMeson && !doNewBaryon) {
    flavNew.id = pickLightQ();
    if ( (flavOld.id > 0 && flavOld.id < 9) || flavOld.id < -1000 ) 
      flavNew.id = -flavNew.id; 

    // Done for simple-quark case. 
    return flavNew;
  }

  // Case: 0 = q -> B B, 1 = q -> B M B, 2 = qq -> M B.
  int iCase = flavNew.nPop;
  if (flavOld.nPop == 1) iCase = 2; 

  // Flavour of popcorn quark (= q shared between B and Bbar).
  if (doNewBaryon) { 
    double sPopWT = dWT[iCase][0];
    if (iCase == 1) sPopWT *= scbBM[0] * popcornSpair;
    double rndmFlav = (2. + sPopWT) * rndmPtr->flat();
    flavNew.idPop = 1;
    if (rndmFlav > 1.) flavNew.idPop = 2;
    if (rndmFlav > 2.) flavNew.idPop = 3;
  } else flavNew.idPop = flavOld.idPop;
  
  // Flavour of vertex quark.
  double sVtxWT = dWT[iCase][1];
  if (flavNew.idPop >= 3) sVtxWT = dWT[iCase][2]; 
  if (flavNew.idPop > 3) sVtxWT *= 0.5 * (1. + 1./dWT[iCase][4]);
  double rndmFlav = (2. + sVtxWT) * rndmPtr->flat();
  flavNew.idVtx = 1;
  if (rndmFlav > 1.) flavNew.idVtx = 2;
  if (rndmFlav > 2.) flavNew.idVtx = 3;

  // Special case for light flavours, possibly identical. 
  if (flavNew.idPop < 3 && flavNew.idVtx < 3) {  
    flavNew.idVtx = flavNew.idPop;
    if (rndmPtr->flat() > dWT[iCase][3]) flavNew.idVtx = 3 - flavNew.idPop;
  }

  // Pick 2 * spin + 1.
  int spin = 3;
  if (flavNew.idVtx != flavNew.idPop) {
    double spinWT = dWT[iCase][6];
    if (flavNew.idVtx == 3) spinWT = dWT[iCase][5];
    if (flavNew.idPop >= 3) spinWT = dWT[iCase][4];
    if ((1. + spinWT) * rndmPtr->flat() < 1.) spin = 1;
  }

  // Form outgoing diquark. Done. 
  flavNew.id = 1000 * max(flavNew.idVtx, flavNew.idPop) 
    + 100 * min(flavNew.idVtx, flavNew.idPop) + spin;
  if ( (flavOld.id < 0 && flavOld.id > -9) || flavOld.id > 1000 ) 
    flavNew.id = -flavNew.id; 
  return flavNew;

}

//--------------------------------------------------------------------------

// Combine two flavours (including diquarks) to produce a hadron.
// The weighting of the combination may fail, giving output 0.

int StringFlav::combine(FlavContainer& flav1, FlavContainer& flav2) {

  // Recognize largest and smallest flavour.
  int id1Abs = abs(flav1.id); 
  int id2Abs = abs(flav2.id);
  int idMax = max(id1Abs, id2Abs);
  int idMin = min(id1Abs, id2Abs);
 
  // Construct a meson.
  if (idMax < 9 || idMin > 1000) {

    // Popcorn meson: use only vertex quarks. Fail if none.
    if (idMin > 1000) {
      id1Abs = flav1.idVtx;
      id2Abs = flav2.idVtx;
      idMax = max(id1Abs, id2Abs);
      idMin = min(id1Abs, id2Abs);
      if (idMin == 0) return 0;
    }

    // Pick spin state and preliminary code.
    int flav = (idMax < 3) ? 0 : idMax - 2;
    double rndmSpin = mesonRateSum[flav] * rndmPtr->flat();
    int spin = -1;
    do rndmSpin -= mesonRate[flav][++spin];
    while (rndmSpin > 0.);
    int idMeson = 100 * idMax + 10 * idMin + mesonMultipletCode[spin];

    // For nondiagonal mesons distinguish particle/antiparticle.
    if (idMax != idMin) {
      int sign = (idMax%2 == 0) ? 1 : -1;
      if ( (idMax == id1Abs && flav1.id < 0) 
        || (idMax == id2Abs && flav2.id < 0) ) sign = -sign;
      idMeson *= sign;  

    // For light diagonal mesons include uubar - ddbar - ssbar mixing.
    } else if (flav < 2) {
      double rMix = rndmPtr->flat();
      if      (rMix < mesonMix1[flav][spin]) idMeson = 110;
      else if (rMix < mesonMix2[flav][spin]) idMeson = 220;
      else                                   idMeson = 330;
      idMeson += mesonMultipletCode[spin];

      // Additional suppression of eta and eta' may give failure.
      if (idMeson == 221 && etaSup < rndmPtr->flat()) return 0;
      if (idMeson == 331 && etaPrimeSup < rndmPtr->flat()) return 0;
    }

    // Finished for mesons.
    return idMeson;
  }

  // SU(6) factors for baryon production may give failure.
  int idQQ1 = idMax / 1000;
  int idQQ2 = (idMax / 100) % 10;
  int spinQQ = idMax % 10;
  int spinFlav = spinQQ - 1;
  if (spinFlav == 2 && idQQ1 != idQQ2) spinFlav = 4;
  if (idMin != idQQ1 && idMin != idQQ2) spinFlav++;   
  if (baryonCGSum[spinFlav] < rndmPtr->flat() * baryonCGMax[spinFlav]) 
    return 0;

  // Order quarks to form baryon. Pick spin.
  int idOrd1 = max( idMin, max( idQQ1, idQQ2) ); 
  int idOrd3 = min( idMin, min( idQQ1, idQQ2) ); 
  int idOrd2 = idMin + idQQ1 + idQQ2 - idOrd1 - idOrd3;
  int spinBar = (baryonCGSum[spinFlav] * rndmPtr->flat() 
    < baryonCGOct[spinFlav]) ? 2 : 4;
  
  // Distinguish Lambda- and Sigma-like. 
  bool LambdaLike = false;
  if (spinBar == 2 && idOrd1 > idOrd2 && idOrd2 > idOrd3) {
    LambdaLike = (spinQQ == 1);
    if (idOrd1 != idMin && spinQQ == 1) LambdaLike = (rndmPtr->flat() < 0.25); 
    else if (idOrd1 != idMin)           LambdaLike = (rndmPtr->flat() < 0.75);
  }

  // Form baryon code and return with sign.  
  int idBaryon = (LambdaLike) 
    ? 1000 * idOrd1 + 100 * idOrd3 + 10 * idOrd2 + spinBar
    : 1000 * idOrd1 + 100 * idOrd2 + 10 * idOrd3 + spinBar;
   return (flav1.id > 0) ? idBaryon : -idBaryon;

}

//--------------------------------------------------------------------------

// Assign popcorn quark inside an original (= rank 0) diquark.

void StringFlav::assignPopQ(FlavContainer& flav) {

  // Safety check that intended to do something.
  int idAbs = abs(flav.id);
  if (flav.rank > 0 || idAbs < 1000) return;

  // Make choice of popcorn quark.
  int id1 = (idAbs/1000)%10;
  int id2 = (idAbs/100)%10;
  double pop2WT = 1.;
       if (id1 == 3) pop2WT = scbBM[1];
  else if (id1 >  3) pop2WT = scbBM[2];  
       if (id2 == 3) pop2WT /= scbBM[1];
  else if (id2 >  3) pop2WT /= scbBM[2];
  // Agrees with Patrik code, but opposite to intention??
  flav.idPop = ((1. + pop2WT) * rndmPtr->flat() > 1.) ? id2 : id1; 
  flav.idVtx = id1 + id2 - flav.idPop;

  // Also determine if to produce popcorn meson.
  flav.nPop = 0;
  double popWT = popS[0];
  if (id1 == 3) popWT = popS[1];
  if (id2 == 3) popWT = popS[2];
  if (idAbs%10 == 1) popWT *= sqrt(probQQ1toQQ0);
  if ((1. + popWT) * rndmPtr->flat() > 1.) flav.nPop = 1;
  
}

//--------------------------------------------------------------------------

// Combine two quarks to produce a diquark. 
// Normally according to production composition, but nonvanishing idHad
// means diquark from known hadron content, so use SU(6) wave fucntion.

int StringFlav::makeDiquark(int id1, int id2, int idHad) {

  // Initial values.
  int idMin = min( abs(id1), abs(id2));
  int idMax = max( abs(id1), abs(id2));
  int spin = 1;

  // Select spin of diquark formed from two valence quarks in proton.
  // (More hadron cases??)
  if (abs(idHad) == 2212 || abs(idHad) == 2112) {
    if (idMin == 1 && idMax == 2 && rndmPtr->flat() < 0.75) spin = 0;  

  // Else select spin of diquark according to production composition.
  } else {
    if (idMin != idMax && rndmPtr->flat() > probQQ1norm) spin = 0; 
  }   

  // Combined diquark code.
  int idNewAbs = 1000 * idMax + 100 * idMin + 2 * spin + 1; 
  return (id1 > 0) ? idNewAbs : -idNewAbs; 

}
 
//==========================================================================

// The StringZ class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// When a or c are close to special cases, default to these.
const double StringZ::CFROMUNITY = 0.01;
const double StringZ::AFROMZERO  = 0.02;
const double StringZ::AFROMC     = 0.01;

// Do not take exponent of too large or small number.
const double StringZ::EXPMAX     = 50.; 

//--------------------------------------------------------------------------

// Initialize data members of the string z selection.

void StringZ::init(Settings& settings, ParticleData& particleData, 
  Rndm* rndmPtrIn) {

  // Save pointer.
  rndmPtr       = rndmPtrIn;

  // c and b quark masses.
  mc2           = pow2( particleData.m0(4)); 
  mb2           = pow2( particleData.m0(5)); 

  // Paramaters of Lund/Bowler symmetric fragmentation function.
  aLund         = settings.parm("StringZ:aLund");
  bLund         = settings.parm("StringZ:bLund");
  aExtraDiquark = settings.parm("StringZ:aExtraDiquark");
  rFactC        = settings.parm("StringZ:rFactC");
  rFactB        = settings.parm("StringZ:rFactB");
  rFactH        = settings.parm("StringZ:rFactH");

  // Flags and parameters of nonstandard Lund fragmentation functions.
  useNonStandC  = settings.flag("StringZ:useNonstandardC");
  useNonStandB  = settings.flag("StringZ:useNonstandardB");
  useNonStandH  = settings.flag("StringZ:useNonstandardH");
  aNonC         = settings.parm("StringZ:aNonstandardC");
  aNonB         = settings.parm("StringZ:aNonstandardB");
  aNonH         = settings.parm("StringZ:aNonstandardH");
  bNonC         = settings.parm("StringZ:bNonstandardC");
  bNonB         = settings.parm("StringZ:bNonstandardB");
  bNonH         = settings.parm("StringZ:bNonstandardH");

  // Flags and parameters of Peterson/SLAC fragmentation function.
  usePetersonC  = settings.flag("StringZ:usePetersonC");
  usePetersonB  = settings.flag("StringZ:usePetersonB");
  usePetersonH  = settings.flag("StringZ:usePetersonH");
  epsilonC      = settings.parm("StringZ:epsilonC");
  epsilonB      = settings.parm("StringZ:epsilonB");
  epsilonH      = settings.parm("StringZ:epsilonH");

  // Parameters for joining procedure.
  stopM         = settings.parm("StringFragmentation:stopMass");
  stopNF        = settings.parm("StringFragmentation:stopNewFlav");
  stopS         = settings.parm("StringFragmentation:stopSmear");

}

//--------------------------------------------------------------------------

// Generate the fraction z that the next hadron will take, 
// using either Lund/Bowler or, for heavy, Peterson/SLAC functions.
// Note: for a heavy new coloured particle we assume pT negligible.

double StringZ::zFrag( int idOld, int idNew, double mT2) {

  // Find if old or new flavours correspond to diquarks.
  int idOldAbs = abs(idOld);
  int idNewAbs = abs(idNew);
  bool isOldDiquark = (idOldAbs > 1000 && idOldAbs < 10000);
  bool isNewDiquark = (idNewAbs > 1000 && idNewAbs < 10000);

  // Find heaviest quark in fragmenting parton/diquark.
  int idFrag = idOldAbs;
  if (isOldDiquark) idFrag = max( idOldAbs / 1000, (idOldAbs / 100) % 10);
  
  // Use Peterson where explicitly requested for heavy flavours.
  if (idFrag == 4 && usePetersonC) return zPeterson( epsilonC);
  if (idFrag == 5 && usePetersonB) return zPeterson( epsilonB);
  if (idFrag >  5 && usePetersonH) {
    double epsilon = epsilonH * mb2 / mT2; 
    return zPeterson( epsilon);
  }

  // Nonstandard a and b values implemented for heavy flavours.
  double aNow = aLund;
  double bNow = bLund;
  if (idFrag == 4 && useNonStandC) {
    aNow = aNonC;
    bNow = bNonC;
  } else if (idFrag == 5 && useNonStandB) {
    aNow = aNonB;
    bNow = bNonB;
  } else if (idFrag >  5 && useNonStandH) {
    aNow = aNonH;
    bNow = bNonH;
  }

  // Shape parameters of Lund symmetric fragmentation function.
  double aShape = aNow;
  if (isOldDiquark) aShape += aExtraDiquark;
  double bShape = bNow * mT2;
  double cShape = 1.;
  if (isOldDiquark) cShape -= aExtraDiquark;
  if (isNewDiquark) cShape += aExtraDiquark;
  if (idFrag == 4) cShape += rFactC * bNow * mc2;
  if (idFrag == 5) cShape += rFactB * bNow * mb2;
  if (idFrag >  5) cShape += rFactH * bNow * mT2;
  return zLund( aShape, bShape, cShape);

}

//--------------------------------------------------------------------------

// Generate a random z according to the Lund/Bowler symmetric
// fragmentation function f(z) = (1 -z)^a * exp(-b/z) / z^c.
// Normalized so that f(z_max) = 1  it can also be written as
// f(z) = exp( a * ln( (1 - z) / (1 - z_max) ) + b * (1/z_max - 1/z) 
//           + c * ln(z_max/z) ).  

double StringZ::zLund( double a, double b, double c) {

  // Special cases for c = 1, a = 0 and a = c. 
  bool cIsUnity = (abs( c - 1.) < CFROMUNITY);
  bool aIsZero = (a < AFROMZERO);
  bool aIsC = (abs(a - c) < AFROMC);

  // Determine position of maximum.
  double zMax;
  if (aIsZero) zMax = (c > b) ? b / c : 1.; 
  else if (aIsC) zMax = b / (b + c);
  else { zMax = 0.5 * (b + c - sqrt( pow2(b - c) + 4. * a * b)) / (c - a);
         if (zMax > 0.9999 && b > 100.) zMax = min(zMax, 1. - a / b); }   
        
  // Subdivide z range if distribution very peaked near either endpoint.
  bool peakedNearZero = (zMax < 0.1);
  bool peakedNearUnity = (zMax > 0.85 && b > 1.);

  // Find integral of trial function everywhere bigger than f. 
  // (Dummy start values.)
  double fIntLow = 1.; 
  double fIntHigh = 1.; 
  double fInt = 2.; 
  double zDiv = 0.5; 
  double zDivC = 0.5;
  // When z_max is small use that f(z)
  //   < 1     for z < z_div = 2.75 * z_max,
  //   < (z_div/z)^c for z > z_div (=> logarithm for c = 1, else power).   
  if (peakedNearZero) {
    zDiv = 2.75 * zMax;
    fIntLow = zDiv; 
    if (cIsUnity) fIntHigh = -zDiv * log(zDiv);
    else { zDivC = pow( zDiv, 1. - c);
           fIntHigh = zDiv * (1. - 1./zDivC) / (c - 1.);} 
    fInt = fIntLow + fIntHigh;
  // When z_max large use that f(z)
  //   < exp( b * (z - z_div) ) for z < z_div with z_div messy expression,
  //   < 1   for z > z_div.
  // To simplify expressions the integral is extended to z =  -infinity.
  } else if (peakedNearUnity) {
    double rcb = sqrt(4. + pow2(c / b));
    zDiv = rcb - 1./zMax - (c / b) * log( zMax * 0.5 * (rcb + c / b) );  
    if (!aIsZero) zDiv += (a/b) * log(1. - zMax);
    zDiv = min( zMax, max(0., zDiv));
    fIntLow = 1. / b;
    fIntHigh = 1. - zDiv; 
    fInt = fIntLow + fIntHigh;
  }

  // Choice of z, preweighted for peaks at low or high z. (Dummy start values.)
  double z = 0.5;
  double fPrel = 1.; 
  double fVal = 1.; 
  do { 
    // Choice of z flat good enough for distribution peaked in the middle;
    // if not this z can be reused as a random number in general.
    z = rndmPtr->flat();
    fPrel = 1.;
    // When z_max small use flat below z_div and 1/z^c above z_div.
    if (peakedNearZero) {
      if (fInt * rndmPtr->flat() < fIntLow) z = zDiv * z;
      else if (cIsUnity) {z = pow( zDiv, z); fPrel = zDiv / z;}
      else { z = pow( zDivC + (1. - zDivC) * z, 1. / (1. - c) );
             fPrel = pow( zDiv / z, c); }
    // When z_max large use exp( b * (z -z_div) ) below z_div 
    // and flat above it.
    } else if (peakedNearUnity) {
      if (fInt * rndmPtr->flat() < fIntLow) { 
        z = zDiv + log(z) / b;
        fPrel = exp( b * (z - zDiv) ); 
      } else z = zDiv + (1. - zDiv) * z; 
    }  

    // Evaluate actual f(z) (if in physical range) and correct.
    if (z > 0 && z < 1) {
      double fExp = b * (1. / zMax - 1. / z)+ c * log(zMax / z);
      if (!aIsZero) fExp += a * log( (1. - z) / (1. - zMax) );
      fVal = exp( max( -EXPMAX, min( EXPMAX, fExp) ) ) ;
    } else fVal = 0.;
  } while (fVal < rndmPtr->flat() * fPrel);

  // Done.
  return z;

}

//--------------------------------------------------------------------------

// Generate a random z according to the Peterson/SLAC formula
// f(z) = 1 / ( z * (1 - 1/z - epsilon/(1-z))^2 )
//      = z * (1-z)^2 / ((1-z)^2 + epsilon * z)^2.

double StringZ::zPeterson( double epsilon) {

  double z, fVal;

  // For large epsilon pick z flat and reject, 
  // knowing that 4 * epsilon * f(z) < 1 everywhere.
  if (epsilon > 0.01) { 
    do { 
      z = rndmPtr->flat();
      fVal = 4. * epsilon * z * pow2(1. - z) 
        / pow2( pow2(1. - z) + epsilon * z);
    } while (fVal < rndmPtr->flat());
    return z; 
  } 
  
  // Else split range, using that 4 * epsilon * f(z) 
  //   < 4 * epsilon / (1 - z)^2 for 0 < z < 1 - 2 * sqrt(epsilon)
  //   < 1                       for 1 - 2 * sqrt(epsilon) < z < 1
  double epsRoot = sqrt(epsilon);
  double epsComb = 0.5 / epsRoot - 1.;
  double fIntLow = 4. * epsilon * epsComb;
  double fInt = fIntLow + 2. * epsRoot;
  do { 
    if (rndmPtr->flat() * fInt < fIntLow) {
      z = 1. - 1. / (1. + rndmPtr->flat() * epsComb);
      fVal = z * pow2( pow2(1. - z) / (pow2(1. - z) + epsilon * z) );
    } else {
      z = 1. - 2. * epsRoot * rndmPtr->flat();
      fVal = 4. * epsilon * z * pow2(1. - z) 
        / pow2( pow2(1. - z) + epsilon * z);
    }
  } while (fVal < rndmPtr->flat());
  return z; 

} 
 
//==========================================================================

// The StringPT class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// To avoid division by zero one must have sigma > 0.
const double StringPT::SIGMAMIN     = 0.2;

//--------------------------------------------------------------------------

// Initialize data members of the string pT selection.

void StringPT::init(Settings& settings,  ParticleData& , Rndm* rndmPtrIn) {

  // Save pointer.
  rndmPtr        = rndmPtrIn;

  // Parameters of the pT width and enhancement.
  double sigma     = settings.parm("StringPT:sigma");
  sigmaQ           = sigma / sqrt(2.);
  enhancedFraction = settings.parm("StringPT:enhancedFraction");
  enhancedWidth    = settings.parm("StringPT:enhancedWidth");

  // Parameter for pT suppression in MiniStringFragmentation.
  sigma2Had        = 2. * pow2( max( SIGMAMIN, sigma) );
  
}

//--------------------------------------------------------------------------

// Generate Gaussian pT such that <p_x^2> = <p_x^2> = sigma^2 = width^2/2,
// but with small fraction multiplied up to a broader spectrum.

pair<double, double> StringPT::pxy() {

  double sigma = sigmaQ;
  if (rndmPtr->flat() < enhancedFraction) sigma *= enhancedWidth;
  pair<double, double> gauss2 = rndmPtr->gauss2();
  return pair<double, double>(sigma * gauss2.first, sigma * gauss2.second);

}
  
//==========================================================================

} // end namespace Pythia8
