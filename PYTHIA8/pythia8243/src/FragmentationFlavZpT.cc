// FragmentationFlavZpT.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// StringFlav, StringZ and StringPT classes.

#include "Pythia8/FragmentationFlavZpT.h"

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

void StringFlav::init(Settings& settings, ParticleData* particleDataPtrIn,
  Rndm* rndmPtrIn, Info* infoPtrIn) {

  // Save pointers.
  rndmPtr         = rndmPtrIn;
  particleDataPtr = particleDataPtrIn;
  infoPtr         = infoPtrIn;

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

  // Spin parameters for combining two quarks to a diquark.
  vector<double> pQQ1tmp = settings.pvec("StringFlav:probQQ1toQQ0join");
  for (int i = 0; i < 4; ++i)
    probQQ1join[i] = 3. * pQQ1tmp[i] / (1. + 3. * pQQ1tmp[i]);

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
    // Fill in rates for multiplication.
    mesMixRate1[0][spin] = mesonMix1[0][spin];
    mesMixRate2[0][spin] = mesonMix2[0][spin] - mesonMix1[0][spin];
    mesMixRate3[0][spin] = 1.0 - mesMixRate1[0][spin] - mesMixRate2[0][spin];
    mesMixRate1[1][spin] = mesonMix1[1][spin];
    mesMixRate2[1][spin] = mesonMix2[1][spin] - mesonMix1[1][spin];
    mesMixRate3[1][spin] = 1.0 - mesMixRate1[1][spin] - mesMixRate2[1][spin];
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

  // Use Gaussian model but with mT2 suppression?
  mT2suppression   = settings.flag("StringPT:mT2suppression");
  sigmaHad         = (sqrt(2.0)*settings.parm("StringPT:sigma"));
  widthPreStrange  = settings.parm("StringPT:widthPreStrange");
  widthPreDiquark  = settings.parm("StringPT:widthPreDiquark");
  useWidthPre      = (widthPreStrange > 1.0) || (widthPreDiquark > 1.0);

  // Enhanded-rate prefactor for MPIs and/or nearby string pieces.
  closePacking     = settings.flag("StringPT:closePacking");
  exponentMPI      = settings.parm("StringPT:expMPI");
  exponentNSP      = settings.parm("StringPT:expNSP");

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

  // Use thermal model?
  thermalModel = settings.flag("StringPT:thermalModel");
  if (thermalModel || mT2suppression) {

    // Temperature parameters for thermal model.
    temperature      = settings.parm("StringPT:temperature");
    tempPreFactor    = settings.parm("StringPT:tempPreFactor");

    // Hadron multiplets in thermal model.
    mesonNonetL1     = settings.flag("StringFlav:mesonNonetL1");
    nNewQuark        = settings.mode("StringFlav:nQuark");

    // Fill list of possible hadrons that are allowed to be produced.
    // Also include a list of "emergency" hadrons that are needed to get
    // rid of all possible endpoint (di)quarks.
    vector<int> hadIDsProd, hadIDsHvyC, hadIDsHvyB;

    // Baryon octet and decuplet.
    int baryonLight[18] = { 2112, 2212, 3112, 3122, 3212, 3222, 3312, 3322,
                            1114, 2114, 2214, 2224, 3114, 3214, 3224, 3314,
                            3324, 3334 };
    int baryonHvyC[22]  = { 4112, 4122, 4132, 4212, 4222, 4232, 4312, 4322,
                            4332, 4412, 4422, 4432,
                            4114, 4214, 4224, 4314, 4324, 4334, 4414, 4424,
                            4434, 4444 };
    int baryonHvyB[35]  = { 5112, 5122, 5132, 5142, 5212, 5222, 5232, 5242,
                            5312, 5322, 5332, 5342, 5412, 5422, 5432, 5442,
                            5512, 5522, 5532, 5542,
                            5114, 5214, 5224, 5314, 5324, 5334, 5414, 5424,
                            5434, 5444, 5514, 5524, 5534, 5544, 5554 };
    for (int i = 0; i < 18; i++) hadIDsProd.push_back( baryonLight[i] );
    // Check how many heavy baryons to include.
    if (nNewQuark > 4) {
      for (int i = 0; i < 35; i++) hadIDsProd.push_back( baryonHvyB[i] );
    } else {
      // Only include lightest combinations.
      int bBar[9] = { 5112, 5122, 5132, 5212, 5222, 5232, 5312, 5322, 5332 };
      for (int i = 0; i < 9; i++) {
        hadIDsHvyB.push_back(  bBar[i] );
        hadIDsHvyB.push_back( -bBar[i] );
      }
    }
    if (nNewQuark > 3) {
      for (int i = 0; i < 22; i++) hadIDsProd.push_back( baryonHvyC[i] );
    } else {
      // Only include lightest combinations.
      int cBar[9] = { 4112, 4122, 4132, 4212, 4222, 4232, 4312, 4322, 4332 };
      for (int i = 0; i < 9; i++) {
        hadIDsHvyC.push_back(  cBar[i] );
        hadIDsHvyC.push_back( -cBar[i] );
      }
    }
    // Antibaryons.
    int sizeNow = int(hadIDsProd.size());
    for (int i = 0; i < sizeNow; i++) hadIDsProd.push_back( -hadIDsProd[i] );

    // Mesons nonets. Take pseudoscalar PDG codes as basis.
    int mesonPSLight[9] = { 311, 321, 211, -311, -321, -211, 111, 221, 331 };
    int mesonPSHvyC[7]  = { 411, 421, 431, -411, -421, -431, 441 };
    int mesonPSHvyB[9]  = { 511, 521, 531, 541, -511, -521, -531, -541, 551 };
    vector<int> mesonPS;
    for (int i = 0; i < 9; i++) mesonPS.push_back( mesonPSLight[i] );
    // Check how many heavy mesons to include. If not included in ordinary
    // production, fill minimal list with "emergency" hadrons
    if (nNewQuark > 4) {
      for (int i = 0; i < 9; i++) mesonPS.push_back( mesonPSHvyB[i] );
    } else {
      // Include all possible combinations, only pseudoscalar as they
      // are the lightest ones.
      int bMes[10] = { 511, 521, 531, 541, -511, -521, -531, -541, 551 };
      for (int i = 0; i < 10; i++) hadIDsHvyB.push_back( bMes[i] );
    }
    if (nNewQuark > 3) {
      for (int i = 0; i < 7; i++) mesonPS.push_back( mesonPSHvyC[i] );
    } else {
      // Include all possible combinations, only pseudoscalar as they
      // are the lightest ones.
      int cMes[8] = { 411, 421, 431, -411, -421, -431, 441 };
      for (int i = 0; i < 8; i++) hadIDsHvyC.push_back( cMes[i] );
    }
    int nMeson = int(mesonPS.size());
    // Pseudoscalar nonet J=0, S=0, L=0.
    for (int i = 0; i < nMeson; i++)
      hadIDsProd.push_back( mesonPS[i] );
    // Vector nonet J=1, S=1, L=0.
    for (int i = 0; i < nMeson; i++)
      hadIDsProd.push_back( mesonPS[i] + (mesonPS[i] > 0 ? 2 : -2) );
    // Include L=1 nonets?
    if (mesonNonetL1) {
      // Pseudovector nonet J=1, S=0, L=1.
      for (int i = 0; i < nMeson; i++)
        hadIDsProd.push_back( mesonPS[i] + (mesonPS[i] > 0 ? 10002 : -10002) );
      // Scalar nonet J=0, S=1, L=1.
      for (int i = 0; i < nMeson; i++)
        hadIDsProd.push_back( mesonPS[i] + (mesonPS[i] > 0 ? 10000 : -10000) );
      // Pseudovector nonet J=1, S=1, L=1.
      for (int i = 0; i < nMeson; i++)
        hadIDsProd.push_back( mesonPS[i] + (mesonPS[i] > 0 ? 20002 : -20002) );
      // Tensor nonet J=2, S=1, L=1.
      for (int i = 0; i < nMeson; i++)
        hadIDsProd.push_back( mesonPS[i] + (mesonPS[i] > 0 ? 4 : -4) );
    }

    // Fill list of all hadrons ids (ordinary and "emergency").
    vector<int> hadIDsAll;
    for (int i = 0; i < int(hadIDsProd.size()); i++)
      hadIDsAll.push_back( hadIDsProd[i] );
    for (int i = 0; i < int(hadIDsHvyC.size()); i++)
      hadIDsAll.push_back( hadIDsHvyC[i] );
    for (int i = 0; i < int(hadIDsHvyB.size()); i++)
      hadIDsAll.push_back( hadIDsHvyB[i] );

    // Fill map with IDs of hadron constituents for all hadrons.
    for (int i = 0; i < int(hadIDsAll.size()); i++) {
      int id    = hadIDsAll[i];
      int idAbs = abs(id);
      vector< pair<int,int> > quarkCombis;
      // Baryon can be split into q + qq in several different ways.
      if (particleDataPtr->isBaryon(id)) {
        bool isOctet   = ( (idAbs % 10) == 2 );
        int  q3        = (idAbs/10)   % 10;
        int  q2        = (idAbs/100)  % 10;
        int  q1        = (idAbs/1000) % 10;
        bool threeFlav = q1 != q2 && q1 != q3 && q2 != q3;
        // Baryon octet J=1/2.
        if (isOctet) {
          if (threeFlav) {
            // Add (q2+q3)_0/1 + q1.
            // if (q2 < q3) (q2+q3)_0 and if (q2 > q3) (q2+q3)_1.
            int j = (q2 < q3) ? 1 : 3;
            int qn[2]  = { min( q3, q2), max( q3, q2) };
            addQuarkDiquark(quarkCombis, q1,
              1000 * qn[1] + 100 * qn[0] + j, id);
            // Add other combinations. Can be both, J=0 or J=1.
            for (j = 1; j < 4; j += 2) {
              // (q1+q3)j + q2
              addQuarkDiquark(quarkCombis, q2, 1000 * q1 + 100 * q3 + j, id);
              // (q1+q2)j + q3
              addQuarkDiquark(quarkCombis, q3, 1000 * q1 + 100 * q2 + j, id);
            }
          } else {
            // Quarks with the same flavour form J=1,
            // all other combinations can be both, J=0 or J=1.
            for (int j = 1; j < 4; j += 2) {
              // (q1+q2)1 + q3
              if ( j == 3 || q1 != q2 )
                addQuarkDiquark(quarkCombis, q3, 1000 * q1 + 100 * q2 + j, id);
              // (q1+q3)1 + q2
              if ( j == 3 || q1 != q3 )
                addQuarkDiquark(quarkCombis, q2, 1000 * q1 + 100 * q3 + j, id);
              // (q2+q3)1 + q1
              if ( j == 3 || q2 != q3 )
                addQuarkDiquark(quarkCombis, q1, 1000 * q2 + 100 * q3 + j, id);
            }
          }
        }
        // Baryon decuplet J=3/2.
        else {
          // All quark pairs form diquarks with J=1.
          // (q1+q2)1 + q3
          addQuarkDiquark(quarkCombis, q3, 1000 * q1 + 100 * q2 + 3, id);
          // (q1+q3)1 + q2
          addQuarkDiquark(quarkCombis, q2, 1000 * q1 + 100 * q3 + 3, id);
          // (q2+q3)1 + q1
          addQuarkDiquark(quarkCombis, q1, 1000 * q2 + 100 * q3 + 3, id);
        }
      // Mesons usually have a trivial subdivision into quark + antiquark.
      // Mixing of diagonal mesons is taken into account later.
      } else {
        int q1        = (idAbs/100) % 10;
        bool uptype1  = (q1 % 2 == 0);
        int q2        = (idAbs/10)  % 10;
        bool uptype2  = (q2 % 2 == 0);
        int quark     = q1;
        int antiQuark = q2;
        // id > 0: downtype+uptype: up = quark, down = antiquark (default)
        // id > 0: same type -> larger id decides
        if ( uptype2 && !uptype1 ) swap( quark, antiQuark);
        if ( (q1 > q2 && !uptype1 && !uptype2)
          || (q2 > q1 &&  uptype2 &&   uptype1) ) swap( quark, antiQuark);
        if (id < 0) swap( quark, antiQuark);
        quarkCombis.push_back( make_pair( quark, -antiQuark) );
      }
      hadronConstIDs[id] = quarkCombis;
    }

    // Copy into smaller versions (one for ordinary production, two for
    // "emergency")
    map< int, vector< pair<int,int> > > hadConstIDsC, hadConstIDsB,
                                        hadConstIDsProd;
    for (int i=0; i<int(hadIDsAll.size()); i++) {
      int id = hadIDsAll[i];
      if (find(hadIDsProd.begin(), hadIDsProd.end(), id) != hadIDsProd.end())
        hadConstIDsProd[id] = hadronConstIDs[id];
      if (find(hadIDsHvyC.begin(), hadIDsHvyC.end(), id) != hadIDsHvyC.end())
        hadConstIDsC[id] = hadronConstIDs[id];
      if (find(hadIDsHvyB.begin(), hadIDsHvyB.end(), id) != hadIDsHvyB.end())
        hadConstIDsB[id] = hadronConstIDs[id];
    }
    map< int, map< int, vector< pair<int,int> > > > hadConstIDsHvy;
    hadConstIDsHvy[4] = hadConstIDsC;
    hadConstIDsHvy[5] = hadConstIDsB;

    // List with all possible initial (di)quarks we could get.
    int inIDs[26]    = { 1, 2, 3, 4, 5, 1103, 2203, 3303, 2101, 2103, 3101,
                         3103, 3201, 3203, 4101, 4103, 4201, 4203, 4301,
                         4303, 5101, 5103, 5201, 5203, 5301, 5303 };
    int inIDsHvyC[2] = { 4403, -4403 };
    int inIDsHvyB[6] = { 5503, -5503, 5401, -5401, 5403, -5403 };
    vector<int> incomingIDs;
    for (int i = 0; i < 26; i++) {
      incomingIDs.push_back( inIDs[i]);
      incomingIDs.push_back(-inIDs[i]);
    }
    // If we include heavy quark hadrons we include the following diquarks in
    // addition.
    if (nNewQuark > 3) {
      for (int i = 0; i < 2; i++) incomingIDs.push_back(inIDsHvyC[i]);
      if (nNewQuark > 4) {
        for (int i = 0; i < 6; i++) incomingIDs.push_back( inIDsHvyB[i]);
      }
    }
    int nIncome = int(incomingIDs.size());

    // Loop over list with all possible initial (di)quarks.
    // Fill map possibleHadrons with
    // key = initial (di)quark id, value = list of possible hadron ids
    //                                     + nr in hadronConstIDs.
    for (int iIDin = 0; iIDin < nIncome; iIDin++) {
      int idIn    = incomingIDs[iIDin];
      int idInAbs = abs(idIn);
      map< int, vector< pair<int,int> > > hadConstIDsNow = hadConstIDsProd;
      // For heavy quarks add "emergency" list, if needed.
      for (int iHvy = nNewQuark+1; iHvy <= 5; iHvy++) {
        if (particleDataPtr->nQuarksInCode(idInAbs, iHvy) > 0)
          for (map< int, vector< pair<int,int> > >::iterator
               it = hadConstIDsHvy[iHvy].begin();
               it != hadConstIDsHvy[iHvy].end(); ++it)
            hadConstIDsNow[it->first] = it->second;
      }
      // Fill list: first parameter of pair is hadron ID, second is nr of
      // hadron constituents in the list.
      vector< pair<int,int> > possibleHadronIDs;
      // Loop through list with hadrons and their (di)quark content,
      // check if possible to produce given the choice of initial (di)quark.
      for (map< int, vector< pair<int,int> > >::iterator
           it = hadConstIDsNow.begin(); it != hadConstIDsNow.end(); ++it) {
        vector< pair<int,int> > constituentIDs = it->second;
        int nConst   = int(constituentIDs.size());
        int hadronID = it->first;
        // Loop over constituent IDs.
        for (int iConst = 0; iConst < nConst; iConst++) {
          int ID1 = constituentIDs[iConst].first;
          int ID2 = constituentIDs[iConst].second;
          if ( (ID1 == idIn) || (ID2 == idIn) ) {
            possibleHadronIDs.push_back( make_pair(hadronID,iConst) );
            // To include uubar-ddbar-ssbar mixing include all diagonal mesons.
            if ( (idInAbs < 4) && (ID1 == -ID2) ) {
              if (idInAbs == 1) {
                possibleHadronIDs.push_back( make_pair(hadronID+110,iConst) );
                possibleHadronIDs.push_back( make_pair(hadronID+220,iConst) );
              } else if (idInAbs == 2) {
                possibleHadronIDs.push_back( make_pair(hadronID-110,iConst) );
                possibleHadronIDs.push_back( make_pair(hadronID+110,iConst) );
              } else if (idInAbs == 3) {
                possibleHadronIDs.push_back( make_pair(hadronID-220,iConst) );
                possibleHadronIDs.push_back( make_pair(hadronID-110,iConst) );
              }
            }
          }
        }
      }
      if (int(possibleHadronIDs.size()) < 1)
        infoPtr->errorMsg("Error in StringFlav::init: no possible "
          "hadrons found");
      possibleHadrons[idIn] = possibleHadronIDs;
    }

    // Calculate baryon octet and decuplet weighting factors
    // based on Clebsch-Gordan coefficients and spin counting.
    // Parameters: qDi1 qDi2 q3 spin.
    // Zero for flavour=0 and same flavour diquarks with J=0.
    for (int q1 = 0; q1 < 6; q1++) {
      for (int q2 = 0; q2 < 6; q2++) {
        baryonOctWeight[q1][q1][q2][0] = 0.0; // qq0 + r
        baryonDecWeight[q1][q1][q2][0] = 0.0; // qq0 + r
        for (int spin = 0; spin < 1; spin++) {
          baryonOctWeight[ 0][q1][q2][spin] = 0.0;
          baryonOctWeight[q1][ 0][q2][spin] = 0.0;
          baryonOctWeight[q1][q2][ 0][spin] = 0.0;
          baryonDecWeight[ 0][q1][q2][spin] = 0.0;
          baryonDecWeight[q1][ 0][q2][spin] = 0.0;
          baryonDecWeight[q1][q2][ 0][spin] = 0.0;
        }
      }
    }
    // Clebsch-Gordon for the rest.
    for (int q1 = 1; q1 < 6; q1++) {
      baryonOctWeight[q1][q1][q1][1] = 0.0; // qq1 + q
      baryonDecWeight[q1][q1][q1][1] = 1.0;
      for (int q2 = 1; q2 < 6; q2++) if (q1!=q2) {
        baryonOctWeight[q1][q1][q2][1] = 0.1667; // qq1 + r
        baryonDecWeight[q1][q1][q2][1] = 0.3333;
        baryonOctWeight[q1][q2][q1][0] = 0.75;   // qr0 + q
        baryonDecWeight[q1][q2][q1][0] = 0.0;
        baryonOctWeight[q2][q1][q1][0] = 0.75;   // rq0 + q
        baryonDecWeight[q2][q1][q1][0] = 0.0;
        baryonOctWeight[q1][q2][q1][1] = 0.0833; // qr1 + q
        baryonDecWeight[q1][q2][q1][1] = 0.6667;
        baryonOctWeight[q2][q1][q1][1] = 0.0833; // rq1 + q
        baryonDecWeight[q2][q1][q1][1] = 0.6667;
        for (int q3 = 0; q3 < 6; q3++) if ((q1 != q3) && (q2 != q3)) {
          baryonOctWeight[q1][q2][q3][0] = 0.5;    // qr0 + s
          baryonDecWeight[q1][q2][q3][0] = 0.0;
          baryonOctWeight[q1][q2][q3][1] = 0.1667; // qr1 + s
          baryonDecWeight[q1][q2][q3][1] = 0.3333;
        }
      }
    }
    // Spin 1 diquarks get extra factor of 3. And all factors
    // get relative baryon-to-meson ratio.
    double BtoMratio = settings.parm("StringFlav:BtoMratio");
    for (int q1 = 0; q1 < 6; q1++) {
      for (int q2 = 0; q2 < 6; q2++) {
        for (int q3 = 0; q3 < 6; q3++) {
          for (int spin = 0; spin < 2; spin++) {
            baryonOctWeight[q1][q2][q3][spin] *= BtoMratio;
            baryonDecWeight[q1][q2][q3][spin] *= BtoMratio;
            if (spin == 1) {
              baryonOctWeight[q1][q2][q3][1] *= 3.0;
              baryonDecWeight[q1][q2][q3][1] *= 3.0;
            }
          }
        }
      }
    }

    // Go through the list of possible hadrons and calculate the prefactor
    // that will multiply the rate.
    double strSup = settings.parm("StringFlav:StrangeSuppression");
    for (int iIDin = 0; iIDin < nIncome; iIDin++) {
      int idIn      = incomingIDs[iIDin];
      int idInAbs   = abs(idIn);
      vector< pair<int,int> > possibleHadronsNow = possibleHadrons[idIn];
      vector<double> prefactors;
      for (int iHad = 0; iHad < int(possibleHadronsNow.size()); iHad++) {
        double prefacNow = 1.0;
        // Get hadron and constituents.
        int hadronID    = possibleHadronsNow[iHad].first;
        int hadronIDabs = abs(hadronID);
        int iConst      = possibleHadronsNow[iHad].second;
        int ID1         = hadronConstIDs[hadronID][iConst].first;
        int ID2         = hadronConstIDs[hadronID][iConst].second;
        // Extra suppression factor for s/c/b quarks.
        double nHeavy   = 0.0;
        for (int i = 3; i <= 5; i++) {
          nHeavy += particleDataPtr->nQuarksInCode( ID1, i);
          nHeavy += particleDataPtr->nQuarksInCode( ID2, i);
        }
        prefacNow      *= pow(strSup, nHeavy);
        if (particleDataPtr->isMeson(hadronID)) {
          // Extra factor according to last digit for spin counting.
          prefacNow *= (abs(hadronID) % 10);
          // Include correct uubar-ddbar-ssbar mixing factor;
          if ( (idInAbs < 4) && (ID1 == -ID2) ) {
            int flav = ( (idInAbs < 3) ? 0 : 1 );
            // Get spin used as counter for the different multiplets
            int spin = getMesonSpinCounter(hadronID);
            double mesonMix[3] = { mesMixRate1[flav][spin],
                                   mesMixRate2[flav][spin],
                                   mesMixRate3[flav][spin] };
            prefacNow *= mesonMix[abs(ID1)-1];
          }
        } else {
          // Check if baryon is octet or decuplet.
          bool isOct = ((hadronIDabs % 10) == 2);
          // Make sure ID2 is diquark.
          if (abs(ID2) < abs(ID1)) swap(ID1,ID2);
          // Extract quark flavours and spin from diquark.
          int Q1      = ( (abs(ID2)/1000) % 10 );
          int Q2      = ( (abs(ID2)/100)  % 10 );
          int diqSpin = ( ((abs(ID2) % 10) == 1) ? 0 : 1 );
          // Single quark.
          int Q3      = abs(ID1);
          // Find Clebsch-Gordan: q1 in DQ | q2 in DQ | q3 | S of DQ
          if (isOct) prefacNow *= baryonOctWeight[Q1][Q2][Q3][diqSpin];
          else       prefacNow *= baryonDecWeight[Q1][Q2][Q3][diqSpin];
          // Special cases for Lamda (312) and Sigma (321) or the like.
          if ( isOct && (Q1!=Q2) && (Q1!=Q3) && (Q2!=Q3) ) {
            // Extract the two lightest quarks from hadron.
            int Qhad1   = ( (hadronIDabs/10)  % 10 );
            int Qhad2   = ( (hadronIDabs/100) % 10 );
            int QhadMin = min(Qhad1,Qhad2);
            int QhadMax = max(Qhad1,Qhad2);
            // Extract the two quarks from the diquark.
            int QdiqMin = min(Q1,Q2);
            int QdiqMax = max(Q1,Q2);
            // Don't do anything if (12) or (21) is diquark.
            if ( !((QdiqMin == QhadMin) && (QdiqMax == QhadMax)) ) {
              // Sigma (321)
              if (Qhad2 > Qhad1) prefacNow *= ( (diqSpin == 0) ? 0.75 : 0.25 );
              // Lamda (312)
              else               prefacNow *= ( (diqSpin == 0) ? 0.25 : 0.27 );
            }
          }
        }
        // Save prefactor.
        prefactors.push_back(prefacNow);
      }
      possibleRatePrefacs[idIn] = prefactors;
    }

    // Now the same again for joining the last two (di)quarks into hadron.
    for (int iIDin1 = 0; iIDin1 < nIncome; iIDin1++) {
      int idIn1     = incomingIDs[iIDin1];
      int idIn1Abs  = abs(idIn1);
      // Loop over possible partners, start with next quark.
      for (int iIDin2 = iIDin1+1; iIDin2 < nIncome; iIDin2++) {
        int idIn2      = incomingIDs[iIDin2];
        int idIn2Abs   = abs(idIn2);
        int idInNow[2] = { min(idIn1,idIn2), max(idIn1,idIn2) };
        pair<int,int> inPair = pair<int,int>(idInNow[0], idInNow[1]);
        // Skip all combinations with two diquarks.
        if ( (idIn1Abs > 1000) && (idIn2Abs > 1000) ) continue;
        // Skip all combinations with two quarks or two antiquarks.
        if ( ( ((idIn1 > 0) && (idIn2 > 0)) || ((idIn1 < 0) && (idIn2 < 0)) )
             && (idIn1Abs < 10) && (idIn2Abs < 10) ) continue;
        // Skip all combinations with quark-antidiquark and
        // antiquark-diquark. (1 = diquark, 2 = quark not possible).
        if ( ((idIn2 >  1000) && (idIn1Abs < 10) && (idIn1 < 0)) ||
             ((idIn2 < -1000) && (idIn1Abs < 10) && (idIn1 > 0)) ) continue;
        // If we are not including heavy quarks skip combinations
        // of heavy quark - diquark with heavy quark.
        if ((idIn1Abs < 10) && (idIn2Abs > 1000)) {
          vector< pair<int,int> > hvyCombs;
          if (nNewQuark < 5) {
            hvyCombs.push_back(make_pair(4,4));
            if (nNewQuark < 4) {
              hvyCombs.push_back(make_pair(5,4));
              hvyCombs.push_back(make_pair(4,5));
              hvyCombs.push_back(make_pair(5,5));
            }
          }
          bool skip = false;
          for (int iComb = 0; iComb < int(hvyCombs.size()); iComb++) {
            int idNow[2] = { hvyCombs[iComb].first, hvyCombs[iComb].second };
            if ( (particleDataPtr->nQuarksInCode(idIn2Abs,idNow[0]) > 0) &&
                 (idIn1Abs == idNow[1]) ) skip = true;
          }
          if (skip) continue;
        }
        // Now decide which list of possible hadrons to use.
        // As we might have to use the special list for heavy quarks we
        // use the maximum of the absolute ids in case of two quarks and
        // check the maximum flavour in case of quark - diquark pair.
        int idUse;
        if ( (idIn1Abs < 10) && (idIn2Abs < 10) ) { // quark - quark
          idUse = ( (idIn1Abs > idIn2Abs) ? idIn1 : idIn2 );
        } else { // quark - diquark
          // Check if diquark contains a heavier flavour then the quark.
          bool useDiquark = false;
          for (int plus = 1; plus < 5; plus++)
            if (particleDataPtr->nQuarksInCode(idIn2Abs, idIn1Abs + plus) > 0)
              useDiquark = true;
          idUse = ( useDiquark ? idIn2 : idIn1 );
        }
        vector<double> possibleRatePrefacsNow = possibleRatePrefacs[idUse];
        vector< pair<int,int> > possibleHadronsNow = possibleHadrons[idUse];
        // New list to fill.
        vector< pair<int,int> > possibleHadronsNew;
        vector<double> possibleRatePrefacsNew;
        // Now loop over possible hadrons and check if other (di)quark
        // in constituents matches idIn2.
        for (int iHad = 0; iHad < int(possibleHadronsNow.size()); iHad++) {
          // Get constituents.
          int hadronID = possibleHadronsNow[iHad].first;
          int iConst   = possibleHadronsNow[iHad].second;
          int ID1      = hadronConstIDs[hadronID][iConst].first;
          int ID2      = hadronConstIDs[hadronID][iConst].second;
          if ( ((ID1 == idIn1) && (ID2 == idIn2)) ||
               ((ID1 == idIn2) && (ID2 == idIn1)) ) {
            // Can take this combination.
            possibleHadronsNew.push_back(possibleHadronsNow[iHad]);
            possibleRatePrefacsNew.push_back(possibleRatePrefacsNow[iHad]);
          }
        }
        if (int(possibleHadronsNew.size()) < 1)
          infoPtr->errorMsg("Error in StringFlav::init: no possible "
            "hadrons found for last two");
        // Save.
        possibleRatePrefacsLast[inPair] = possibleRatePrefacsNew;
        possibleHadronsLast[inPair]     = possibleHadronsNew;
      }
    }
  }

  // Initialize winning parameters.
  hadronIDwin   = 0;
  idNewWin      = 0;
  hadronMassWin = -1.0;

}

//--------------------------------------------------------------------------

// Pick a new flavour (including diquarks) given an incoming one for
// Gaussian pTq^2 distribution.

FlavContainer StringFlav::pickGauss(FlavContainer& flavOld) {

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

// Pick a hadron, based on generated pT value and initial (di)quark.
// Check all possible hadrons and calculate their relative suppression
// based on exp(-mThadron/T), possibly multiplied by spin counting, meson
// mixing or baryon weighting factors.
// First return value is hadron ID, second new (di)quark ID.

FlavContainer StringFlav::pickThermal(FlavContainer& flavOld,
  double pT, double nNSP) {

  // Initial values for new flavour.
  FlavContainer flavNew;
  flavNew.rank = flavOld.rank + 1;

  int idIn        = flavOld.id;
  int idInAbs     = abs(idIn);
  double temprNow = temperature;
  // Temperature increase to work against asymmetry. Apply for
  // s/c/b and diquarks.
  if (idInAbs > 2) temprNow *= tempPreFactor;
  // Enhanded-rate prefactor for MPIs and/or nearby string pieces.
  if (closePacking) {
    temprNow     *= pow(max(1.0,double(infoPtr->nMPI())), exponentMPI);
    temprNow     *= pow(max(1.0,nNSP), exponentNSP);
  }
  // Get Gaussian width in case of mT2 suppression.
  double sigmaNow = sigmaHad;

  // Prefactor for strange quarks and diquarks.
  if (useWidthPre) {
    if (abs(idIn) > 10) sigmaNow *= widthPreDiquark;
    sigmaNow     *= pow(widthPreStrange,
                    particleDataPtr->nQuarksInCode(idIn,3) );
  }

  // Enhanded-rate prefactor for MPIs and/or nearby string pieces.
  if (closePacking) {
    sigmaNow     *= pow(max(1.0,double(infoPtr->nMPI())), exponentMPI);
    sigmaNow     *= pow(max(1.0,nNSP), exponentNSP);
  }

  // Get the list of allowed hadrons and constituents for that
  // initial (di)quark. First parameter of pair is hadron ID, second
  // is nr of hadron constituents in the list.
  vector<double> possibleRatePrefacsNow      = possibleRatePrefacs[idIn];
  vector< pair<int,int> > possibleHadronsNow = possibleHadrons[idIn];
  int nPossHads = int(possibleHadronsNow.size());
  if (nPossHads < 1) {
    infoPtr->errorMsg("Error in StringFlav::pickThermal: no possible "
      "hadrons found");
    return 0;
  }

  // Vector with hadron masses. Is -1.0 if m0 is use for calculating
  // the suppression rate and mSel if mSel is used.
  vector<double> possibleHadronMasses;

  // Calculate rates/suppression factors for given pT.
  vector<double> rates;
  double rateSum = 0.0;
  for (int iHad = 0; iHad < nPossHads; iHad++) {
    int hadronID = possibleHadronsNow[iHad].first;
    // Pick mass and calculate suppression factor.
    double mass  = particleDataPtr->mSel(hadronID);
    possibleHadronMasses.push_back(mass);
    double rate  = exp( -sqrt(pow2(pT)+pow2(mass))/temprNow );
    // mT2 suppression with Gaussian pT?
    if (mT2suppression) rate = exp( -(pow2(pT)+pow2(mass))/pow2(sigmaNow) );
    // Multiply rate with prefactor.
    rate *= possibleRatePrefacsNow[iHad];
    // Save rate and add to sum
    rates.push_back(rate);
    rateSum += rate;
  }
  // Normalize rates
  for (int iHad = 0; iHad < nPossHads; iHad++) rates[iHad] /= rateSum;

  // Get accumulated rates
  vector<double> accumRates;
  for (int iHad = 0; iHad < nPossHads; iHad++) accumRates.push_back(0);
  for (int iHad1 = 0; iHad1 < nPossHads; iHad1++)
    for (int iHad2 = 0; iHad2 <= iHad1; iHad2++)
      accumRates[iHad1] += rates[iHad2];

  // Random number to decide which hadron to pick
  double rand       = rndmPtr->flat();
  int hadronID      = 0;
  int iConst        = 0;
  double hadronMass = -1.0;
  for (int iHad = 0; iHad < nPossHads; iHad++) {
    if (rand <= accumRates[iHad]) {
      hadronID   = possibleHadronsNow[iHad].first;
      iConst     = possibleHadronsNow[iHad].second;
      hadronMass = possibleHadronMasses[iHad];
      break;
    }
  }

  // Get flavour of (di)quark to use next time.
  int idNext = 0;
  vector< pair<int,int> > constituentIDs = hadronConstIDs[hadronID];
  // Mesons
  if (particleDataPtr->isMeson(hadronID)) {
    int ID1 = constituentIDs[0].first;
    int ID2 = constituentIDs[0].second;
    // Special case for diagonal meson, flavour remains
    if (ID1 == -ID2) idNext = idIn;
    else idNext = (idIn == ID1 ? -ID2 : -ID1);
  }
  // Baryons
  else {
    int ID1 = constituentIDs[iConst].first;
    int ID2 = constituentIDs[iConst].second;
    if (ID1 == idIn) idNext = -ID2;
    if (ID2 == idIn) idNext = -ID1;
  }

  // Save new flavour and hadron.
  flavNew.id    = -idNext;  // id used to build hadron
  hadronIDwin   = hadronID;
  idNewWin      = idNext;   // id used in next step
  hadronMassWin = hadronMass;

  // Done.
  return flavNew;

}

//--------------------------------------------------------------------------

// Combine two flavours (including diquarks) to produce a hadron.
// The weighting of the combination may fail, giving output 0.

int StringFlav::combine(FlavContainer& flav1, FlavContainer& flav2) {

  // Recognize largest and smallest flavour.
  int id1Abs = abs(flav1.id);
  int id2Abs = abs(flav2.id);
  int idMax  = max(id1Abs, id2Abs);
  int idMin  = min(id1Abs, id2Abs);

  // Construct a meson.
  if (idMax < 9 || idMin > 1000) {

    // Popcorn meson: use only vertex quarks. Fail if none.
    if (idMin > 1000) {
      id1Abs = flav1.idVtx;
      id2Abs = flav2.idVtx;
      idMax  = max(id1Abs, id2Abs);
      idMin  = min(id1Abs, id2Abs);
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

// Combine two flavours (including diquarks) to produce a hadron. Function
// called in case of combining the two remaining flavours into last hadron.

int StringFlav::combineLastThermal(FlavContainer& flav1, FlavContainer& flav2,
  double pT, double nNSP) {

  // Decide randomly on whether to treat flav1 or flav2 as incoming.
  int idIn[2]    = { flav1.id, flav2.id };
  if (rndmPtr->flat() < 0.5) swap(idIn[0], idIn[1]);
  int idInNow[2] = { min(idIn[0],idIn[1]), max(idIn[0],idIn[1]) };

  int idInAbs     = abs(idIn[0]);
  double temprNow = temperature;
  // Temperature increase to work against asymmetry. Apply for
  // s/c/b and diquarks.
  if (idInAbs > 2) temprNow *= tempPreFactor;

  // Enhanded-rate prefactor for MPIs and/or nearby string pieces.
  if (closePacking) {
    temprNow     *= pow(max(1.0,double(infoPtr->nMPI())), exponentMPI);
    temprNow     *= pow(max(1.0,nNSP), exponentNSP);
  }

  // Get Gaussian width in case of mT2 suppression.
  double sigmaNow = sigmaHad;
  // Prefactor for strange quarks and diquarks.
  if (useWidthPre) {
    if (abs(idInAbs) > 10) sigmaNow *= widthPreDiquark;
    sigmaNow     *= pow(widthPreStrange,
                    particleDataPtr->nQuarksInCode(idInAbs,3) );
  }

  // Enhanded-rate prefactor for MPIs and/or nearby string pieces.
  if (closePacking) {
    sigmaNow     *= pow(max(1.0,double(infoPtr->nMPI())), exponentMPI);
    sigmaNow     *= pow(max(1.0,nNSP), exponentNSP);
  }

  // Get the list of allowed hadrons and constituents for that combination
  // of (di)quarks. First parameter of pair is hadron ID, second
  // is nr of hadron constituents in the list.
  pair<int,int> inPr = pair<int,int>(idInNow[0], idInNow[1]);
  vector<double> possibleRatePrefacsNow      = possibleRatePrefacsLast[inPr];
  vector< pair<int,int> > possibleHadronsNow = possibleHadronsLast[inPr];
  int nPossHads = int(possibleHadronsNow.size());
  if (nPossHads < 1) {
    infoPtr->errorMsg("Error in StringFlav::combineLastThermal: no possible "
      "hadrons found for last two");
    return 0;
  }

  // Vector with hadron masses. Is -1.0 if m0 is use for calculating
  // the suppression rate and mSel if mSel is used.
  vector<double> possibleHadronMasses;

  // Calculate rates/suppression factors for given pT.
  vector<double> rates;
  double rateSum = 0.0;
  for (int iHad = 0; iHad < nPossHads; iHad++) {
    int hadronID = possibleHadronsNow[iHad].first;
    // Pick mass and calculate suppression factor.
    double mass = particleDataPtr->mSel(hadronID);
    possibleHadronMasses.push_back(mass);
    double rate = exp( -sqrt(pow2(pT)+pow2(mass))/temprNow );
    // mT2 suppression with Gaussian pT?
    if (mT2suppression) rate = exp( -(pow2(pT)+pow2(mass))/pow2(sigmaNow) );
    // Multiply rate with prefactor.
    rate *= possibleRatePrefacsNow[iHad];
    // Save rate and add to sum
    rates.push_back(rate);
    rateSum += rate;
  }
  // Normalize rates
  for (int iHad = 0; iHad < nPossHads; iHad++) rates[iHad] /= rateSum;

  // Get accumulated rates
  vector<double> accumRates;
  for (int iHad = 0; iHad < nPossHads; iHad++) accumRates.push_back(0);
  for (int iHad1 = 0; iHad1 < nPossHads; iHad1++)
    for (int iHad2 = 0; iHad2 <= iHad1; iHad2++)
      accumRates[iHad1] += rates[iHad2];

  // Random number to decide which hadron to pick
  double rand       = rndmPtr->flat();
  int hadronID      = 0;
  double hadronMass = -1.0;
  for (int iHad = 0; iHad < nPossHads; iHad++) {
    if (rand <= accumRates[iHad]) {
      hadronID   = possibleHadronsNow[iHad].first;
      hadronMass = possibleHadronMasses[iHad];
      break;
    }
  }

  // Save hadron.
  hadronIDwin   = hadronID;
  hadronMassWin = hadronMass;

  // Done.
  return hadronIDwin;
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
  if      (id1 == 3) pop2WT = scbBM[1];
  else if (id1 >  3) pop2WT = scbBM[2];
  if      (id2 == 3) pop2WT /= scbBM[1];
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
// means diquark from known hadron content, so use SU(6) wave function.

int StringFlav::makeDiquark(int id1, int id2, int idHad) {

  // Initial values.
  int idMin = min( abs(id1), abs(id2));
  int idMax = max( abs(id1), abs(id2));
  int spin = 1;

  // Select spin of diquark formed from two valence quarks in proton.
  // (More hadron cases??)
  if (abs(idHad) == 2212 || abs(idHad) == 2112) {
    if (idMin == 1 && idMax == 2 && rndmPtr->flat() < 0.75) spin = 0;

  // Else select spin of diquark according to assumed spin-1 suppression.
  } else if (idMin != idMax) {
    if (rndmPtr->flat() > probQQ1join[min(idMax,5) - 2]) spin = 0;
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
  Rndm* rndmPtrIn, Info* infoPtrIn) {

  // Save pointers.
  rndmPtr       = rndmPtrIn;
  infoPtr       = infoPtrIn;

  // c and b quark masses.
  mc2           = pow2( particleData.m0(4));
  mb2           = pow2( particleData.m0(5));

  // Paramaters of Lund/Bowler symmetric fragmentation function.
  aLund         = settings.parm("StringZ:aLund");
  bLund         = settings.parm("StringZ:bLund");
  aExtraSQuark  = settings.parm("StringZ:aExtraSQuark");
  aExtraDiquark = settings.parm("StringZ:aExtraDiquark");
  rFactC        = settings.parm("StringZ:rFactC");
  rFactB        = settings.parm("StringZ:rFactB");
  rFactH        = settings.parm("StringZ:rFactH");

  // Alternative parameterisation of Lund FF using average z(rho) instead of b.
  if (settings.flag("StringZ:deriveBLund")) {
    if (!deriveBLund(settings, particleData)) {
      infoPtr->errorMsg("Error in StringZ::init: Derivation of b parameter "
        " failed. Reverting to default.");
      settings.resetParm("StringZ:bLund");
    }
  }

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

// Alternative parameterisation of the Lund function. Derive the bLund
// parameter given the average z for fixed a and mT2.

bool StringZ::deriveBLund(Settings& settings, ParticleData& particleData) {

  // Set up using reference mT2 = mRho^2 + 2*sigmaPT^2
  double mRef   = particleData.m0(113);
  double mT2ref = pow2(mRef) + 2.*pow2(settings.parm("stringPT:sigma"));
  double avgZ   = settings.parm("StringZ:avgZLund");
  double a      = settings.parm("StringZ:aLund");

  // Set up average FF encapsulator. Args: a, b (dummy), c, mT2
  LundFFAvg lundFFAvg;
  vector<double> args(4);
  args[0] = a;
  args[1] = 1.;
  args[2] = 1.;
  args[3] = mT2ref;
  double bNow = 0.;

  // Check if derived b fell inside the nominal range for bLund
  bool check = lundFFAvg.brent(bNow, avgZ, 1, 0.01, 20.0, args, 1.e-6, 1000);
  if (check) {
    settings.parm("StringZ:bLund", bNow, false);

    // Print out derived value for b (and mT2ref), noting if outside range.
    cout << fixed << setprecision(2) << "\n <z(rho)> = " << setw(5)
         << avgZ << " for aLund = "<< a <<" & mT2ref = " << setw(5) << mT2ref
         << " GeV^2 gave bLund = " << setw(5) << bNow << " GeV^-2:";
    if ( bNow == settings.parm("StringZ:bLund") ) cout <<" accepted" << endl;
    else {
      // If outside range, tell user but force anyway so fits can see behaviour
      cout << " accepted (forced)" << endl;
      settings.parm("StringZ:bLund", bNow, true);
    }

    // No further calls needed since b parameter updated in settings database.
    settings.flag("StringZ:deriveBLund", false);
  }
  return check;
}

//--------------------------------------------------------------------------

// Generate the fraction z that the next hadron will take,
// using either Lund/Bowler or, for heavy, Peterson/SLAC functions.
// Note: for a heavy new coloured particle we assume pT negligible.

double StringZ::zFrag( int idOld, int idNew, double mT2) {

  // Find if old or new flavours correspond to diquarks.
  int idOldAbs = abs(idOld);
  int idNewAbs = abs(idNew);
  bool isOldSQuark = (idOldAbs == 3);
  bool isNewSQuark = (idNewAbs == 3);
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
  if (isOldSQuark)  aShape += aExtraSQuark;
  if (isOldDiquark) aShape += aExtraDiquark;
  double bShape = bNow * mT2;
  double cShape = 1.;
  if (isOldSQuark)  cShape -= aExtraSQuark;
  if (isNewSQuark)  cShape += aExtraSQuark;
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

void StringPT::init(Settings& settings,  ParticleData* particleDataPtrIn,
  Rndm* rndmPtrIn, Info* infoPtrIn) {

  // Save pointer.
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  infoPtr         = infoPtrIn;

  // Parameters of the pT width and enhancement.
  double sigma     = settings.parm("StringPT:sigma");
  sigmaQ           = sigma / sqrt(2.);
  enhancedFraction = settings.parm("StringPT:enhancedFraction");
  enhancedWidth    = settings.parm("StringPT:enhancedWidth");
  widthPreStrange  = settings.parm("StringPT:widthPreStrange");
  widthPreDiquark  = settings.parm("StringPT:widthPreDiquark");
  useWidthPre      = (widthPreStrange > 1.0) || (widthPreDiquark > 1.0);

  // Temperature for thermal model.
  thermalModel     = settings.flag("StringPT:thermalModel");
  temperature      = settings.parm("StringPT:temperature");
  tempPreFactor    = settings.parm("StringPT:tempPreFactor");

  // Upper estimate of thermal spectrum: fraction at x = pT_quark/T < 1.
  fracSmallX       = 0.6 / (0.6 + (1.2/0.9) * exp(-0.9));

  // Enhanded-width prefactor for MPIs and/or nearby string pieces.
  closePacking     = settings.flag("StringPT:closePacking");
  exponentMPI      = settings.parm("StringPT:expMPI");
  exponentNSP      = settings.parm("StringPT:expNSP");

  // Parameter for pT suppression in MiniStringFragmentation.
  sigma2Had        = 2. * pow2( max( SIGMAMIN, sigma) );

}

//--------------------------------------------------------------------------

// Generate quark pT according to fitting functions, such that
// hadron pT is generated according to exp(-pT/T) d^2pT.

pair<double, double> StringPT::pxyThermal(int idIn, double nNSP) {

  double temprNow = temperature;
  // Temperature increase to work against asymmetry. Apply for
  // s/c/b and diquarks.
  if (abs(idIn) > 2) temprNow *= tempPreFactor;

  // Enhanded-width prefactor for MPIs and/or nearby string pieces.
  if (closePacking) {
    temprNow *= pow(max(1.0,double(infoPtr->nMPI())), exponentMPI);
    temprNow *= pow(max(1.0,nNSP), exponentNSP);
  }

  // Pick x = pT_quark/T according to K_{1/4}(x)/x^{1/4} * x dx.
  double xrand, approx, wanted;
  do {
    xrand = (rndmPtr->flat() < fracSmallX) ? rndmPtr->flat()
          : 1. - log(rndmPtr->flat()) / 0.9;
    approx = (xrand < 1.) ? 0.6 : 1.2 * exp(-0.9 * xrand);
    wanted = BesselK14(xrand) * pow( xrand, 0.75);
  } while (rndmPtr->flat() * approx > wanted);

  // Find pT_quark. Random number to decide on angle.
  double pTquark = xrand * temprNow;
  double phi     = 2.0 * M_PI * rndmPtr->flat();

  // Done.
  return pair<double, double>( pTquark * cos(phi), pTquark * sin(phi) );

}

//--------------------------------------------------------------------------

// Generate Gaussian pT such that <p_x^2> = <p_x^2> = sigma^2 = width^2/2,
// but with small fraction multiplied up to a broader spectrum.

pair<double, double> StringPT::pxyGauss(int idIn, double nNSP) {

  // Normal (classical) width selection.
  double sigma = sigmaQ;
  if (rndmPtr->flat() < enhancedFraction) sigma *= enhancedWidth;

  // Prefactor for strange quarks and diquarks.
  if (useWidthPre) {
    if (abs(idIn) > 10) sigma *= widthPreDiquark;
    sigma *= pow(widthPreStrange, particleDataPtr->nQuarksInCode(idIn, 3) );
  }

  // Enhanded-width prefactor for MPIs and/or nearby string pieces.
  if (closePacking) {
    sigma *= pow(max(1.0,double(infoPtr->nMPI())), exponentMPI);
    sigma *= pow(max(1.0,nNSP), exponentNSP);
  }

  // Generate (p_x, p_y) pair.
  pair<double, double> gauss2 = rndmPtr->gauss2();
  return pair<double, double>(sigma * gauss2.first, sigma * gauss2.second);

}

//--------------------------------------------------------------------------

// Evaluate Bessel function K_{1/4}(x).
// Use power series for x < 2.5 and asymptotic expansion for x > 2.5.
// Number of terms picked to have accuracy better than 1 per mille.
// Based on M. Abramowitz and I.A. Stegun, eqs. 9.6.2, 9.6.10, 9.7.2.

double StringPT::BesselK14(double x) {

  // Power series expansion of K_{1/4} : k = 0 term.
  if (x < 2.5) {
    double xRat  = 0.25 * x * x;
    double prodP = pow( 0.5 * x, -0.25) / 1.2254167024;
    double prodN = pow( 0.5 * x,  0.25) / 0.9064024771;
    double sum   = prodP - prodN;

    // Power series expansion of K_{1/4} : m > 0 terms.
    for (int k = 1; k < 6; ++k) {
      prodP *= xRat / (k * (k - 0.25));
      prodN *= xRat / (k * (k + 0.25));
      sum   += prodP - prodN;
    }
    sum *= M_PI * sqrt(0.5);
    return sum;

  // Asymptotic expansion of K_{1/4}.
  } else {
    double asym  = sqrt(M_PI * 0.5 / x) * exp(-x);
    double term1 = -         0.75 / ( 8. * x);
    double term2 = -term1 *  8.75 / (16. * x);
    double term3 = -term2 * 24.75 / (24. * x);
    double term4 = -term3 * 48.75 / (32. * x);
    asym *= 1. + term1 + term2 + term3 + term4;
    return asym;
  }
}

//==========================================================================

} // end namespace Pythia8
