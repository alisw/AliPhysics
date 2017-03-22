// SigmaQCD.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// QCD simulation classes. 

#include "SigmaQCD.h"

namespace Pythia8 {

//==========================================================================

// Sigma0AB2AB class.
// Cross section for elastic scattering A B -> A B.

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma0AB2AB::setIdColAcol() {

  // Flavours and colours are trivial. 
  setId( idA, idB, idA, idB);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
}

//==========================================================================

// Sigma0AB2XB class.
// Cross section for single diffractive scattering A B -> X B.

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma0AB2XB::setIdColAcol() {

  // Flavours and colours are trivial. 
  int idX          = 10* (abs(idA) / 10) + 9900000; 
  if (idA < 0) idX = -idX;
  setId( idA, idB, idX, idB);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//==========================================================================

// Sigma0AB2AX class.
// Cross section for single diffractive scattering A B -> A X.

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma0AB2AX::setIdColAcol() {

  // Flavours and colours are trivial. 
  int idX          = 10* (abs(idB) / 10) + 9900000; 
  if (idB < 0) idX = -idX;
  setId( idA, idB, idA, idX);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//==========================================================================

// Sigma0AB2XX class.
// Cross section for double diffractive scattering A B -> X X.

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma0AB2XX::setIdColAcol() {

  // Flavours and colours are trivial. 
  int          idX1 = 10* (abs(idA) / 10) + 9900000; 
  if (idA < 0) idX1 = -idX1;
  int          idX2 = 10* (abs(idB) / 10) + 9900000; 
  if (idB < 0) idX2 = -idX2;
  setId( idA, idB, idX1, idX2);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//==========================================================================

// Sigma0AB2AXB class.
// Cross section for central scattering A B -> A X B.

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma0AB2AXB::setIdColAcol() {
  
  // Central diffractive state represented by rho_diffr0. Colours trivial.
  int idX = 9900110; 
  setId( idA, idB, idA, idB,idX);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  
}

//==========================================================================

// Sigma2gg2gg class.
// Cross section for g g -> g g.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence.

void Sigma2gg2gg::sigmaKin() {

  // Calculate kinematics dependence.
  sigTS  = (9./4.) * (tH2 / sH2 + 2. * tH / sH + 3. + 2. * sH / tH 
           + sH2 / tH2);
  sigUS  = (9./4.) * (uH2 / sH2 + 2. * uH / sH + 3. + 2. * sH / uH 
           + sH2 / uH2);
  sigTU  = (9./4.) * (tH2 / uH2 + 2. * tH / uH + 3. + 2. * uH / tH 
           + uH2 / tH2);
  sigSum = sigTS + sigUS + sigTU;

  // Answer contains factor 1/2 from identical gluons.
  sigma  = (M_PI / sH2) * pow2(alpS) * 0.5 * sigSum;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2gg::setIdColAcol() {

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

// Sigma2gg2qqbar class.
// Cross section for g g -> q qbar (q = u, d, s, i.e. almost massless).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2gg2qqbar::initProc() {

  // Read number of quarks to be considered in massless approximation.
  nQuarkNew       = settingsPtr->mode("HardQCD:nQuarkNew");

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence. 

void Sigma2gg2qqbar::sigmaKin() { 

  // Pick new flavour.
  idNew = 1 + int( nQuarkNew * rndmPtr->flat() ); 
  mNew  = particleDataPtr->m0(idNew);
  m2New = mNew*mNew;
  
  // Calculate kinematics dependence.
  sigTS = 0.;
  sigUS = 0.;
  if (sH > 4. * m2New) {
    sigTS = (1./6.) * uH / tH - (3./8.) * uH2 / sH2;
    sigUS = (1./6.) * tH / uH - (3./8.) * tH2 / sH2; 
  }
  sigSum = sigTS + sigUS;

  // Answer is proportional to number of outgoing flavours.
  sigma  = (M_PI / sH2) * pow2(alpS) * nQuarkNew * sigSum;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2qqbar::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idNew, -idNew);

  // Two colour flow topologies.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else                 setColAcol( 1, 2, 3, 1, 3, 0, 0, 2); 

}

//==========================================================================

// Sigma2qg2qg class.
// Cross section for q g -> q g.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence. 

void Sigma2qg2qg::sigmaKin() { 

  // Calculate kinematics dependence.
  sigTS  = uH2 / tH2 - (4./9.) * uH / sH;
  sigTU  = sH2 / tH2 - (4./9.) * sH / uH;
  sigSum = sigTS + sigTU;

  // Answer.
  sigma  = (M_PI / sH2) * pow2(alpS) * sigSum;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2qg::setIdColAcol() {

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

// Sigma2qq2qq class.
// Cross section for q qbar' -> q qbar' or q q' -> q q' 
// (qbar qbar' -> qbar qbar'), q' may be same as q.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qq2qq::sigmaKin() { 

  // Calculate kinematics dependence for different terms.
  sigT   = (4./9.) * (sH2 + uH2) / tH2;
  sigU   = (4./9.) * (sH2 + tH2) / uH2;
  sigTU  = - (8./27.) * sH2 / (tH * uH);
  sigST  = - (8./27.) * uH2 / (sH * tH);

}

//--------------------------------------------------------------------------


// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qq2qq::sigmaHat() {  

  // Combine cross section terms; factor 1/2 when identical quarks.
  if      (id2 ==  id1) sigSum = 0.5 * (sigT + sigU + sigTU);
  else if (id2 == -id1) sigSum = sigT + sigST;
  else                      sigSum = sigT;

  // Answer.
  return (M_PI/sH2) * pow2(alpS) * sigSum;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qq2qq::setIdColAcol() {

  // Outgoing = incoming flavours.
  setId( id1, id2, id1, id2);

  // Colour flow topologies. Swap when antiquarks.
  if (id1 * id2 > 0)  setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);
  else                setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  if (id2 == id1 && (sigT + sigU) * rndmPtr->flat() > sigT)
                      setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2gg class.
// Cross section for q qbar -> g g.

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence. 

void Sigma2qqbar2gg::sigmaKin() { 

  // Calculate kinematics dependence.
  sigTS  = (32./27.) * uH / tH - (8./3.) * uH2 / sH2;
  sigUS  = (32./27.) * tH / uH - (8./3.) * tH2 / sH2;
  sigSum = sigTS + sigUS;

  // Answer contains factor 1/2 from identical gluons.
  sigma  = (M_PI / sH2) * pow2(alpS) * 0.5 * sigSum;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2gg::setIdColAcol() {

  // Outgoing flavours trivial.
  setId( id1, id2, 21, 21);

  // Two colour flow topologies. Swap if first is antiquark.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 0, 2, 1, 3, 3, 2);
  else                 setColAcol( 1, 0, 0, 2, 3, 2, 1, 3); 
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2qqbarNew class.
// Cross section q qbar -> q' qbar'.

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qqbar2qqbarNew::initProc() {

  // Read number of quarks to be considered in massless approximation.
  nQuarkNew       = settingsPtr->mode("HardQCD:nQuarkNew");

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence. 

void Sigma2qqbar2qqbarNew::sigmaKin() { 

  // Pick new flavour.
  idNew = 1 + int( nQuarkNew * rndmPtr->flat() ); 
  mNew  = particleDataPtr->m0(idNew);
  m2New = mNew*mNew;

  // Calculate kinematics dependence.
  sigS                      = 0.;
  if (sH > 4. * m2New) sigS = (4./9.) * (tH2 + uH2) / sH2; 

  // Answer is proportional to number of outgoing flavours.
  sigma = (M_PI / sH2) * pow2(alpS) * nQuarkNew * sigS;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2qqbarNew::setIdColAcol() {

  // Set outgoing flavours ones.
  id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, id2, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2gg2QQbar class.
// Cross section g g -> Q Qbar (Q = c, b or t).
// Only provided for fixed m3 = m4 so do some gymnastics:
// i) s34Avg picked so that beta34 same when s3, s4 -> s34Avg.
// ii) tHQ = tH - mQ^2 = -0.5 sH (1 - beta34 cos(thetaH)) for m3 = m4 = mQ,
//     but tH - uH = sH beta34 cos(thetaH) also for m3 != m4, so use
//     tH, uH selected for m3 != m4 to derive tHQ, uHQ valid for m3 = m4.   

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2gg2QQbar::initProc() {

  // Process name.
  nameSave                 = "g g -> Q Qbar";
  if (idNew == 4) nameSave = "g g -> c cbar";
  if (idNew == 5) nameSave = "g g -> b bbar";
  if (idNew == 6) nameSave = "g g -> t tbar";
  if (idNew == 7) nameSave = "g g -> b' b'bar";
  if (idNew == 8) nameSave = "g g -> t' t'bar";

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(idNew, -idNew);

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence. 

void Sigma2gg2QQbar::sigmaKin() { 

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  double s34Avg = 0.5 * (s3 + s4) - 0.25 * pow2(s3 - s4) / sH; 
  double tHQ    = -0.5 * (sH - tH + uH);
  double uHQ    = -0.5 * (sH + tH - uH); 
  double tHQ2   = tHQ * tHQ;
  double uHQ2   = uHQ * uHQ;

  // Calculate kinematics dependence.
  double tumHQ = tHQ * uHQ - s34Avg * sH;
  sigTS = ( uHQ / tHQ - 2.25 * uHQ2 / sH2 + 4.5 * s34Avg * tumHQ 
    / ( sH * tHQ2) + 0.5 * s34Avg * (tHQ + s34Avg) / tHQ2 
    - s34Avg*s34Avg / (sH * tHQ) ) / 6.;
  sigUS = ( tHQ / uHQ - 2.25 * tHQ2 / sH2 + 4.5 * s34Avg * tumHQ 
    / ( sH * uHQ2) + 0.5 * s34Avg * (uHQ + s34Avg) / uHQ2 
    - s34Avg*s34Avg / (sH * uHQ) ) / 6.;
  sigSum = sigTS + sigUS;

  // Answer.
  sigma = (M_PI / sH2) * pow2(alpS) * sigSum * openFracPair;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2QQbar::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idNew, -idNew);

  // Two colour flow topologies.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else                 setColAcol( 1, 2, 3, 1, 3, 0, 0, 2); 

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles of W in top decay.

double Sigma2gg2QQbar::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // For top decay hand over to standard routine, else done.
  if (idNew == 6 && process[process[iResBeg].mother1()].idAbs() == 6) 
       return weightTopDecay( process, iResBeg, iResEnd);
  else return 1.; 

}

//==========================================================================

// Sigma2qqbar2QQbar class.
// Cross section q qbar -> Q Qbar (Q = c, b or t).
// Only provided for fixed m3 = m4 so do some gymnastics:
// i) s34Avg picked so that beta34 same when s3, s4 -> s34Avg.
// ii) tHQ = tH - mQ^2 = -0.5 sH (1 - beta34 cos(thetaH)) for m3 = m4 = mQ,
//     but tH - uH = sH beta34 cos(thetaH) also for m3 != m4, so use
//     tH, uH selected for m3 != m4 to derive tHQ, uHQ valid for m3 = m4.   

//--------------------------------------------------------------------------

// Initialize process, especially parton-flux object. 
  
void Sigma2qqbar2QQbar::initProc() {

  // Process name.
  nameSave                 = "q qbar -> Q Qbar";
  if (idNew == 4) nameSave = "q qbar -> c cbar";
  if (idNew == 5) nameSave = "q qbar -> b bbar";
  if (idNew == 6) nameSave = "q qbar -> t tbar";
  if (idNew == 7) nameSave = "q qbar -> b' b'bar";
  if (idNew == 8) nameSave = "q qbar -> t' t'bar";

  // Secondary open width fraction.
  openFracPair = particleDataPtr->resOpenFrac(idNew, -idNew);

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat) - no incoming flavour dependence. 

void Sigma2qqbar2QQbar::sigmaKin() { 

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  double s34Avg = 0.5 * (s3 + s4) - 0.25 * pow2(s3 - s4) / sH; 
  double tHQ    = -0.5 * (sH - tH + uH);
  double uHQ    = -0.5 * (sH + tH - uH); 
  double tHQ2   = tHQ * tHQ;
  double uHQ2   = uHQ * uHQ;

  // Calculate kinematics dependence.
  double sigS = (4./9.) * ((tHQ2 + uHQ2) / sH2 + 2. * s34Avg / sH); 

  // Answer.
  sigma = (M_PI / sH2) * pow2(alpS) * sigS * openFracPair;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2QQbar::setIdColAcol() {

  // Set outgoing flavours.
  id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, id2, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}

//--------------------------------------------------------------------------

// Evaluate weight for decay angles of W in top decay.

double Sigma2qqbar2QQbar::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // For top decay hand over to standard routine, else done.
  if (idNew == 6 && process[process[iResBeg].mother1()].idAbs() == 6) 
       return weightTopDecay( process, iResBeg, iResEnd);
  else return 1.; 

}


//==========================================================================

// Sigma3gg2ggg class.
// Cross section for g g -> g g g.

//--------------------------------------------------------------------------

// Evaluate |M|^2 - no incoming flavour dependence.

void Sigma3gg2ggg::sigmaKin() {

  // Calculate all four-vector products.
  Vec4 p1cm( 0., 0.,  0.5 * mH, 0.5 * mH);
  Vec4 p2cm( 0., 0., -0.5 * mH, 0.5 * mH);
  pp[1][2] = p1cm * p2cm;
  pp[1][3] = p1cm * p3cm;
  pp[1][4] = p1cm * p4cm;
  pp[1][5] = p1cm * p5cm;
  pp[2][3] = p2cm * p3cm;
  pp[2][4] = p2cm * p4cm;
  pp[2][5] = p2cm * p5cm;
  pp[3][4] = p3cm * p4cm;
  pp[3][5] = p3cm * p5cm;
  pp[4][5] = p4cm * p5cm;
  for (int i = 1; i < 5; ++i)
    for (int j = i + 1; j < 6; ++j) pp[j][i] = pp[i][j];       
  
  // Cross section, in three main sections.
  double num1 = cycle(1,2,3,4,5) + cycle(1,2,3,5,4) + cycle(1,2,4,3,5) 
              + cycle(1,2,4,5,3) + cycle(1,2,5,3,4) + cycle(1,2,5,4,3)
              + cycle(1,3,2,4,5) + cycle(1,3,2,5,4) + cycle(1,3,4,2,5)
              + cycle(1,3,5,2,4) + cycle(1,4,2,3,5) + cycle(1,4,3,2,5);
  double num2 = pow4(pp[1][2]) + pow4(pp[1][3]) + pow4(pp[1][4]) 
              + pow4(pp[1][5]) + pow4(pp[2][3]) + pow4(pp[2][4])
              + pow4(pp[2][5]) + pow4(pp[3][4]) + pow4(pp[3][5])
              + pow4(pp[4][5]);
  double den  = pp[1][2] * pp[1][3] * pp[1][4] * pp[1][5] * pp[2][3]
              * pp[2][4] * pp[2][5] * pp[3][4] * pp[3][5] * pp[4][5];

  // Answer has a factor 6 due to identical gluons
  // This is cancelled by phase space factor (1 / 6)
  sigma = pow3(4. * M_PI * alpS) * (27./16.) * num1 * num2 / den;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3gg2ggg::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, 21, 21, 21);

  // Three colour flow topologies, each with two orientations.
  /*
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 4, 4, 3);
  else if (sigRand < sigTS + sigUS) 
                       setColAcol( 1, 2, 3, 1, 3, 4, 4, 2);
  else                 setColAcol( 1, 2, 3, 4, 1, 4, 3, 2); 
  if (rndmPtr->flat() > 0.5) swapColAcol();
  */

  // Temporary solution. 
  setColAcol( 1, 2, 2, 3, 1, 4, 4, 5, 5, 3);  
}


//==========================================================================

// Sigma3qqbar2ggg class.
// Cross section for q qbar -> g g g.

//--------------------------------------------------------------------------

// Evaluate |M|^2 - no incoming flavour dependence.
void Sigma3qqbar2ggg::sigmaKin() {

  // Setup four-vectors
  pCM[0] = Vec4( 0., 0.,  0.5 * mH, 0.5 * mH);
  pCM[1] = Vec4( 0., 0., -0.5 * mH, 0.5 * mH);
  pCM[2] = p3cm;
  pCM[3] = p4cm;
  pCM[4] = p5cm;

  // Calculate |M|^2
  // Answer has a factor 6 due to identical gluons, 
  // which is cancelled by phase space factor (1 / 6)
  sigma = m2Calc();

}

//--------------------------------------------------------------------------

// |M|^2

inline double Sigma3qqbar2ggg::m2Calc() {

  // Calculate four-products
  double sHnow  = (pCM[0] + pCM[1]).m2Calc();
  double sHhalf = sH / 2.;

  // qbar (p+) + q(p-) -> g(k1) g(k2) g(k3)
  // a_i = (p+ . ki), i = 1, 2, 3
  // b_i = (p- . ki), i = 1, 2, 3
  a[0] = pCM[0] * pCM[2];
  a[1] = pCM[0] * pCM[3];
  a[2] = pCM[0] * pCM[4];
  b[0] = pCM[1] * pCM[2];
  b[1] = pCM[1] * pCM[3];
  b[2] = pCM[1] * pCM[4];

  pp[0][1] = pCM[2] * pCM[3];
  pp[1][2] = pCM[3] * pCM[4];
  pp[2][0] = pCM[4] * pCM[2];

  // ab[i][j] = a_i * b_j + a_j * b_i
  ab[0][1] = a[0] * b[1] + a[1] * b[0];
  ab[1][2] = a[1] * b[2] + a[2] * b[1];
  ab[2][0] = a[2] * b[0] + a[0] * b[2];

  // Cross section
  double num1 = a[0] * b[0] * (a[0] * a[0] + b[0] * b[0]) +
                a[1] * b[1] * (a[1] * a[1] + b[1] * b[1]) +
                a[2] * b[2] * (a[2] * a[2] + b[2] * b[2]);
  double den1 = a[0] * a[1] * a[2] * b[0] * b[1] * b[2];
  double num2 = - ( ab[0][1] / pp[0][1] )
                - ( ab[1][2] / pp[1][2] )
                - ( ab[2][0] / pp[2][0] );
  double num3 = a[2] * b[2] * ab[0][1] / (pp[1][2] * pp[2][0] ) +
                a[0] * b[0] * ab[1][2] / (pp[2][0] * pp[0][1] ) +
                a[1] * b[1] * ab[2][0] / (pp[0][1] * pp[1][2] );

  // Final answer
  return  pow3(4. * M_PI * alpS) * (8. / 324.) * (num1 / den1) *
          ( sHhalf + 9. * (sHhalf + num2) + (2. * 81. / sHnow) * num3 );

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3qqbar2ggg::setIdColAcol(){

  // Flavours are trivial.
  setId( id1, id2, 21, 21, 21);

  // Temporary solution. 
  setColAcol( 1, 0, 0, 2, 1, 3, 3, 4, 4, 2);  
  if (id1 < 0) swapColAcol();
}

//--------------------------------------------------------------------------

// Map a final state configuration

inline void Sigma3qqbar2ggg::mapFinal() {
  switch (config) {
  case 0: pCM[2] = p3cm; pCM[3] = p4cm; pCM[4] = p5cm; break;
  case 1: pCM[2] = p3cm; pCM[3] = p5cm; pCM[4] = p4cm; break;
  case 2: pCM[2] = p4cm; pCM[3] = p3cm; pCM[4] = p5cm; break;
  case 3: pCM[2] = p4cm; pCM[3] = p5cm; pCM[4] = p3cm; break;
  case 4: pCM[2] = p5cm; pCM[3] = p3cm; pCM[4] = p4cm; break;
  case 5: pCM[2] = p5cm; pCM[3] = p4cm; pCM[4] = p3cm; break;
  }
}

//==========================================================================

// Sigma3qg2qgg class.
// Cross section for q g -> q g g.
// Crossed relation from q qbar -> g g g:
//   qbar(p+) q(p-) -> g(k1) g(k2) g(k3)

//--------------------------------------------------------------------------

// Evaluate |M|^2 - no incoming flavour dependence
// Note: two different contributions from gq and qg incoming

void Sigma3qg2qgg::sigmaKin() {

  // Pick a final state configuration
  pickFinal();

  // gq and qg incoming
  for (int i = 0; i < 2; i++) {

    // Map incoming four-vectors to p+, p-, k1, k2, k3
    pCM[0] = Vec4( 0., 0.,  0.5 * mH, 0.5 * mH);
    pCM[1] = Vec4( 0., 0., -0.5 * mH, 0.5 * mH);
    mapFinal();

    // Crossing
    swap(pCM[i], pCM[2]);

    // |M|^2
    // XXX - Extra factor of (3) from selecting a final state
    // configuration (already a factor of 2 in the original
    // answer due to two identical final state gluons)???
    // Extra factor of (3 / 8) as average over incoming gluon
    sigma[i] = 3. * (3. / 8.) * m2Calc();

  } // for (int i = 0; i < 2; i++)

}

//--------------------------------------------------------------------------

// Evaluate |M|^2 - incoming flavour dependence
// Pick from two configurations calculated previously

double Sigma3qg2qgg::sigmaHat() {
  // gq or qg incoming
  return (id1 == 21) ? sigma[0] : sigma[1];
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3qg2qgg::setIdColAcol(){
  // Outgoing flavours; only need to know where the quark is
  int qIdx    = config / 2;
  int idTmp[3] = { 21, 21, 21 };
  idTmp[qIdx]  = (id1 == 21) ? id2 : id1;
  setId( id1, id2, idTmp[0], idTmp[1], idTmp[2]);

  // Temporary solution
  if      (qIdx == 0) setColAcol(1, 0, 2, 1, 4, 0, 3, 4, 2, 3);
  else if (qIdx == 1) setColAcol(1, 0, 2, 1, 3, 4, 4, 0, 2, 3);
  else                setColAcol(1, 0, 2, 1, 3, 4, 2, 3, 4, 0);
  // gq or qg incoming
  if (id1 == 21) {
    swap( colSave[1], colSave[2]);
    swap(acolSave[1], acolSave[2]);
  }
  // qbar rather than q incoming
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//==========================================================================

// Sigma3gg2qqbarg class.
// Cross section for g g -> q qbar g
// Crossed relation from q qbar -> g g g

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma3gg2qqbarg::initProc() {

  // Read number of quarks to be considered in massless approximation.
  nQuarkNew       = settingsPtr->mode("HardQCD:nQuarkNew");

}

//--------------------------------------------------------------------------
 
// Evaluate |M|^2 - no incoming flavour dependence.

void Sigma3gg2qqbarg::sigmaKin() {

  // Incoming four-vectors
  pCM[0] = Vec4( 0., 0.,  0.5 * mH, 0.5 * mH);
  pCM[1] = Vec4( 0., 0., -0.5 * mH, 0.5 * mH);

  // Pick and map a final state configuration
  pickFinal();
  mapFinal();

  // Crossing
  swap(pCM[0], pCM[2]);
  swap(pCM[1], pCM[3]);

  // |M|^2
  // Extra factor of (6.) from picking a final state configuration
  // Extra factor of nQuarkNew
  // Extra factor of (3. / 8.) ^ 2 as averaging over two incoming gluons
  sigma = 6. * nQuarkNew * (3. / 8.) * (3. / 8.) * m2Calc();

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3gg2qqbarg::setIdColAcol(){

  // Pick new flavour
  int idNew = 1 + int( nQuarkNew * rndmPtr->flat() ); 

  // Outgoing flavours; easiest just to map by hand
  switch (config) {
  case 0: id3 =  idNew; id4 = -idNew; id5 =  21;    break;
  case 1: id3 =  idNew; id4 =  21;    id5 = -idNew; break;
  case 2: id3 = -idNew; id4 =  idNew; id5 =  21;    break;
  case 3: id3 =  21;    id4 =  idNew; id5 = -idNew; break;
  case 4: id3 = -idNew; id4 =  21;    id5 =  idNew; break;
  case 5: id3 =  21;    id4 = -idNew; id5 =  idNew; break;
  }
  setId(id1, id2, id3, id4, id5);

  // Temporary solution
  switch (config) {
  case 0: setColAcol( 1, 2, 2, 3, 4, 0, 0, 3, 1, 4 ); break;
  case 1: setColAcol( 1, 2, 2, 3, 4, 0, 1, 4, 0, 3 ); break;
  case 2: setColAcol( 1, 2, 2, 3, 0, 3, 4, 0, 1, 4 ); break;
  case 3: setColAcol( 1, 2, 2, 3, 1, 4, 4, 0, 0, 3 ); break;
  case 4: setColAcol( 1, 2, 2, 3, 0, 3, 1, 4, 4, 0 ); break;
  case 5: setColAcol( 1, 2, 2, 3, 1, 4, 0, 3, 4, 0 ); break;
  }
  
}

//==========================================================================

// Sigma3qq2qqgDiff class.
// Cross section for q q' -> q q' g, q != q'

//--------------------------------------------------------------------------

// Evaluate |M|^2 - no incoming flavour dependence

void Sigma3qq2qqgDiff::sigmaKin() {

  // q1(p+) q2(p-) -> q1(q+) q2(q-) g(k)

  // Incoming four-vectors
  pCM[0] = Vec4( 0., 0.,  0.5 * mH, 0.5 * mH);
  pCM[1] = Vec4( 0., 0., -0.5 * mH, 0.5 * mH);
  // Pick and map a final state configuration
  pickFinal();
  mapFinal();

  // |M|^2
  // Extra factor of (6.) from picking a final state configuration
  sigma = 6. * m2Calc();
}

//--------------------------------------------------------------------------

// |M|^2

inline double Sigma3qq2qqgDiff::m2Calc() {

  // Four-products
  s  = (pCM[0] + pCM[1]).m2Calc();
  t  = (pCM[0] - pCM[2]).m2Calc();
  u  = (pCM[0] - pCM[3]).m2Calc();
  up = (pCM[1] - pCM[2]).m2Calc();
  sp = (pCM[2] + pCM[3]).m2Calc();
  tp = (pCM[1] - pCM[3]).m2Calc();

  // |M|^2
  double num1 = (s * s + sp * sp + u * u + up * up) / (t * tp);
  double den1 = (pCM[0] * pCM[4]) * (pCM[1] * pCM[4]) *
                (pCM[2] * pCM[4]) * (pCM[3] * pCM[4]);
  double num2 = (u + up) * (s * sp + t * tp - u * up) +
                u * (s * t + sp * tp) + up * (s * tp + sp * t);
  double num3 = (s + sp) * (s * sp - t * tp - u * up) +
                2. * t * tp * (u + up) + 2. * u * up * (t + tp);
    
  // (N^2 - 1)^2 / 4N^3 = 16. / 27.
  // (N^2 - 1)   / 4N^3 =  2. / 27.
  return (1. / 8.) * pow3(4. * M_PI * alpS) * num1 / den1 *
         ( (16. / 27.) * num2 - (2. / 27.) * num3 );

}

//--------------------------------------------------------------------------

// Evaluate |M|^2 - incoming flavour dependence

double Sigma3qq2qqgDiff::sigmaHat() {
  // Different incoming flavours only
  if (abs(id1) == abs(id2)) return 0.;
  return sigma;
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3qq2qqgDiff::setIdColAcol(){
  
  // Outgoing flavours; easiest just to map by hand
  switch (config) {
  case 0: id3 = id1; id4 = id2; id5 = 21;  break;
  case 1: id3 = id1; id4 = 21;  id5 = id2; break;
  case 2: id3 = id2; id4 = id1; id5 = 21;  break;
  case 3: id3 = 21;  id4 = id1; id5 = id2; break;
  case 4: id3 = id2; id4 = 21;  id5 = id1; break;
  case 5: id3 = 21;  id4 = id2; id5 = id1; break;
  }
  setId(id1, id2, id3, id4, id5);

  // Temporary solution; id1 and id2 can be q/qbar independently
  int cols[5][2];
  if (id1 > 0) {
    cols[0][0] = 1; cols[0][1] = 0;
    cols[2][0] = 1; cols[2][1] = 0;
  } else {
    cols[0][0] = 0; cols[0][1] = 1;
    cols[2][0] = 0; cols[2][1] = 1;
  }
  if (id2 > 0) {
    cols[1][0] = 2; cols[1][1] = 0;
    cols[3][0] = 3; cols[3][1] = 0;
    cols[4][0] = 2; cols[4][1] = 3;
  } else {
    cols[1][0] = 0; cols[1][1] = 2;
    cols[3][0] = 0; cols[3][1] = 3;
    cols[4][0] = 3; cols[4][1] = 2;
  }
  // Map correct final state configuration
  int i3 = 0, i4 = 0, i5 = 0;
  switch (config) {
  case 0: i3 = 2; i4 = 3; i5 = 4; break;
  case 1: i3 = 2; i4 = 4; i5 = 3; break;
  case 2: i3 = 3; i4 = 2; i5 = 4; break;
  case 3: i3 = 4; i4 = 2; i5 = 3; break;
  case 4: i3 = 3; i4 = 4; i5 = 2; break;
  case 5: i3 = 4; i4 = 3; i5 = 2; break;
  }
  // Put colours in place
  setColAcol(cols[0][0],  cols[0][1],  cols[1][0],  cols[1][1],
             cols[i3][0], cols[i3][1], cols[i4][0], cols[i4][1],
             cols[i5][0], cols[i5][1]);

}

//--------------------------------------------------------------------------

// Map a final state configuration

inline void Sigma3qq2qqgDiff::mapFinal() {
  switch (config) {
  case 0: pCM[2] = p3cm; pCM[3] = p4cm; pCM[4] = p5cm; break;
  case 1: pCM[2] = p3cm; pCM[3] = p5cm; pCM[4] = p4cm; break;
  case 2: pCM[2] = p4cm; pCM[3] = p3cm; pCM[4] = p5cm; break;
  case 3: pCM[2] = p4cm; pCM[3] = p5cm; pCM[4] = p3cm; break;
  case 4: pCM[2] = p5cm; pCM[3] = p3cm; pCM[4] = p4cm; break;
  case 5: pCM[2] = p5cm; pCM[3] = p4cm; pCM[4] = p3cm; break;
  }
}


//==========================================================================

// Sigma3qqbar2qqbargDiff
// Cross section for q qbar -> q' qbar' g
// Crossed relation from q q' -> q q' g, q != q'

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma3qqbar2qqbargDiff::initProc() {

  // Read number of quarks to be considered in massless approximation.
  nQuarkNew       = settingsPtr->mode("HardQCD:nQuarkNew");

}

//--------------------------------------------------------------------------

// Evaluate |M|^2 - no incoming flavour dependence.

void Sigma3qqbar2qqbargDiff::sigmaKin() {
  // Overall 6 possibilities for final state ordering
  // To keep symmetry between final states, always map to:
  //  1) q1(p+)    qbar1(p-)  -> qbar2(q+)  q2(q-)     g(k)
  //  2) qbar1(p+) q1(p-)     -> q2(q+)     qbar2(q-)  g(k)
  // Crossing p- and q+ gives:
  //  1) q1(p+)    q2(-q+)    -> q1(-p-)    q2(q-)     g(k)
  //  2) qbar1(p+) qbar2(-q+) -> qbar1(-p-) qbar2(q-)  g(k)

  // Incoming four-vectors
  pCM[0] = Vec4( 0., 0.,  0.5 * mH, 0.5 * mH);
  pCM[1] = Vec4( 0., 0., -0.5 * mH, 0.5 * mH);
  // Pick and map a final state configuration
  pickFinal();
  mapFinal();

  // Crossing
  swap(pCM[1], pCM[2]);
  pCM[1] = -pCM[1];
  pCM[2] = -pCM[2];

  // |M|^2
  // Extra factor of (6.) from picking a final state configuration
  // Extra factor of (nQuarkNew - 1) from new q/qbar pairs
  // XXX - Extra factor of (2.) from second possible crossing???
  sigma = 6. * (nQuarkNew - 1) * 2. * m2Calc();
  
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3qqbar2qqbargDiff::setIdColAcol(){

  // Pick new q qbar flavour with incoming flavour disallowed
  int idNew = 1 + int( (nQuarkNew - 1) * rndmPtr->flat() ); 
  if (idNew >= abs(id1)) ++idNew;
  // For qbar q incoming, q+ is always mapped to q2
  // For q qbar incoming, q+ is always mapped to qbar2
  if (id1 > 0) idNew = -idNew;

  // Outgoing flavours; easiest just to map by hand
  switch (config) {
  case 0: id3 =  idNew; id4 = -idNew; id5 =  21;    break;
  case 1: id3 =  idNew; id4 =  21;    id5 = -idNew; break;
  case 2: id3 = -idNew; id4 =  idNew; id5 =  21;    break;
  case 3: id3 =  21;    id4 =  idNew; id5 = -idNew; break;
  case 4: id3 = -idNew; id4 =  21;    id5 =  idNew; break;
  case 5: id3 =  21;    id4 = -idNew; id5 =  idNew; break;
  }
  setId(id1, id2, id3, id4, id5);

  // Temporary solution; start with q qbar -> qbar q g
  int cols[5][2];
  cols[0][0] = 1; cols[0][1] = 0;
  cols[1][0] = 0; cols[1][1] = 2;
  cols[2][0] = 0; cols[2][1] = 3;
  cols[3][0] = 1; cols[3][1] = 0;
  cols[4][0] = 3; cols[4][1] = 2;
  // Map into correct place
  int i3 = 0, i4 = 0, i5 = 0;
  switch (config) {
  case 0: i3 = 2; i4 = 3; i5 = 4; break;
  case 1: i3 = 2; i4 = 4; i5 = 3; break;
  case 2: i3 = 3; i4 = 2; i5 = 4; break;
  case 3: i3 = 4; i4 = 2; i5 = 3; break;
  case 4: i3 = 3; i4 = 4; i5 = 2; break;
  case 5: i3 = 4; i4 = 3; i5 = 2; break;
  }
  setColAcol(cols[0][0], cols[0][1], cols[1][0], cols[1][1],
             cols[i3][0], cols[i3][1], cols[i4][0], cols[i4][1],
             cols[i5][0], cols[i5][1]);
  // Swap for qbar q incoming
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma3qg2qqqbarDiff class.
// Cross section for q g -> q q' qbar'
// Crossed relation from q q' -> q q' g, q != q'
//   q1(p+) q2(p-) -> q1(q+) q2(q-) g(k)

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma3qg2qqqbarDiff::initProc() {

  // Read number of quarks to be considered in massless approximation.
  nQuarkNew       = settingsPtr->mode("HardQCD:nQuarkNew");

}

//--------------------------------------------------------------------------

// Evaluate |M|^2 - no incoming flavour dependence

void Sigma3qg2qqqbarDiff::sigmaKin() {

  // Pick a final state configuration
  pickFinal();

  // gq or qg incoming
  for (int i = 0; i < 2; i++) {

    // Map incoming four-vectors
    pCM[0] = Vec4( 0., 0.,  0.5 * mH, 0.5 * mH);
    pCM[1] = Vec4( 0., 0., -0.5 * mH, 0.5 * mH);
    mapFinal();

    // Crossing (note extra -ve sign in total sigma)
    swap(pCM[i], pCM[4]);
    pCM[i] = -pCM[i];
    pCM[4] = -pCM[4];

    // |M|^2
    // Extra factor of (6) from picking a final state configuration
    // Extra factor of (3 / 8) as averaging over incoming gluon
    // Extra factor of (nQuarkNew - 1) due to new q/qbar pair
    sigma[i] = -6. * (3. / 8.) * (nQuarkNew - 1) * m2Calc();

  }
  
}

//--------------------------------------------------------------------------

// Evaluate |M|^2 - incoming flavour dependence

double Sigma3qg2qqqbarDiff::sigmaHat() {
  // gq or qg incoming
  return (id1 == 21) ? sigma[0] : sigma[1];
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3qg2qqqbarDiff::setIdColAcol(){
  // Pick new q qbar flavour with incoming flavour disallowed
  int sigmaIdx = (id1 == 21) ? 0 : 1;
  int idIn     = (id1 == 21) ? id2 : id1;
  int idNew    = 1 + int( (nQuarkNew - 1) * rndmPtr->flat() );
  if (idNew >= abs(idIn)) ++idNew;

  // qbar instead of q incoming means swap outgoing q/qbar pair
  int id3Tmp = idIn, id4Tmp = idNew, id5Tmp = -idNew;
  if (idIn < 0)  swap(id4Tmp, id5Tmp);
  // If g q incoming rather than q g, idIn and idNew 
  // should be exchanged (see sigmaKin)
  if (sigmaIdx == 0) swap(id3Tmp, id4Tmp);
  // Outgoing flavours; now just map as if q g incoming
  switch (config) {
  case 0: id3 = id3Tmp; id4 = id4Tmp; id5 = id5Tmp; break;
  case 1: id3 = id3Tmp; id4 = id5Tmp; id5 = id4Tmp; break;
  case 2: id3 = id4Tmp; id4 = id3Tmp; id5 = id5Tmp; break;
  case 3: id3 = id5Tmp; id4 = id3Tmp; id5 = id4Tmp; break;
  case 4: id3 = id4Tmp; id4 = id5Tmp; id5 = id3Tmp; break;
  case 5: id3 = id5Tmp; id4 = id4Tmp; id5 = id3Tmp; break;
  }
  setId(id1, id2, id3, id4, id5);

  // Temporary solution; start with either
  // g q1    -> q1    q2 qbar2
  // g qbar1 -> qbar1 qbar2 q2
  int cols[5][2];
  cols[0][0] = 1; cols[0][1] = 2;
  if (idIn > 0) {
    cols[1][0] = 3; cols[1][1] = 0;
    cols[2][0] = 1; cols[2][1] = 0;
    cols[3][0] = 3; cols[3][1] = 0;
    cols[4][0] = 0; cols[4][1] = 2;
  } else {
    cols[1][0] = 0; cols[1][1] = 3;
    cols[2][0] = 0; cols[2][1] = 2;
    cols[3][0] = 0; cols[3][1] = 3;
    cols[4][0] = 1; cols[4][1] = 0;
  }
  // Swap incoming if q/qbar g instead
  if (id2 == 21) {
    swap(cols[0][0], cols[1][0]);
    swap(cols[0][1], cols[1][1]);
  }
  // Map final state
  int i3 = 0, i4 = 0, i5 = 0;
  if (sigmaIdx == 0) {
    switch (config) {
    case 0: i3 = 3; i4 = 2; i5 = 4; break;
    case 1: i3 = 3; i4 = 4; i5 = 2; break;
    case 2: i3 = 2; i4 = 3; i5 = 4; break;
    case 3: i3 = 4; i4 = 3; i5 = 2; break;
    case 4: i3 = 2; i4 = 4; i5 = 3; break;
    case 5: i3 = 4; i4 = 2; i5 = 3; break;
    }
  } else {
    switch (config) {
    case 0: i3 = 2; i4 = 3; i5 = 4; break;
    case 1: i3 = 2; i4 = 4; i5 = 3; break;
    case 2: i3 = 3; i4 = 2; i5 = 4; break;
    case 3: i3 = 4; i4 = 2; i5 = 3; break;
    case 4: i3 = 3; i4 = 4; i5 = 2; break;
    case 5: i3 = 4; i4 = 3; i5 = 2; break;
    }
  }
  setColAcol(cols[0][0],  cols[0][1],  cols[1][0],  cols[1][1],
             cols[i3][0], cols[i3][1], cols[i4][0], cols[i4][1],
             cols[i5][0], cols[i5][1]);
}

//==========================================================================

// Sigma3qq2qqgSame class.
// Cross section for q q' -> q q' g, q == q'.

//--------------------------------------------------------------------------

// Evaluate |M|^2 - no incoming flavour dependence

void Sigma3qq2qqgSame::sigmaKin() {
  // q1(p+) q2(p-) -> q1(q+) q2(q-) g(k)

  // Incoming four-vectors
  pCM[0] = Vec4( 0., 0.,  0.5 * mH, 0.5 * mH);
  pCM[1] = Vec4( 0., 0., -0.5 * mH, 0.5 * mH);
  // Pick/map a final state configuration
  pickFinal();
  mapFinal();

  // |M|^2
  // Extra factor (3) from picking final state configuration
  // (original answer already has a factor 2 from identical
  // quarks in the final state)
  sigma = 3. * m2Calc();

}

//--------------------------------------------------------------------------

// |M|^2

inline double Sigma3qq2qqgSame::m2Calc() {

  // Four-products
  s  = (pCM[0] + pCM[1]).m2Calc();
  t  = (pCM[0] - pCM[2]).m2Calc();
  u  = (pCM[0] - pCM[3]).m2Calc();
  sp = (pCM[2] + pCM[3]).m2Calc();
  tp = (pCM[1] - pCM[3]).m2Calc();
  up = (pCM[1] - pCM[2]).m2Calc();

  // |M|^2
  ssp  = s * sp;
  ttp  = t * tp;
  uup  = u * up;
  s_sp = s + sp;
  t_tp = t + tp;
  u_up = u + up;

  double den1 = (pCM[0] * pCM[4]) * (pCM[1] * pCM[4]) *
                (pCM[2] * pCM[4]) * (pCM[3] * pCM[4]);

  double fac1 = s * (t * u + tp * up) + sp * (t * up + tp * u);
  double fac2 = ssp - ttp - uup;
  double fac3 = 2. * (ttp * u_up + uup * t_tp);

  double num1 = u_up * (ssp + ttp - uup) + fac1;
  double num2 = s_sp * fac2 + fac3;
  double num3 = (s * s + sp * sp + u * u + up * up) / (t * tp);

  double num4 = t_tp * (ssp - ttp + uup) + fac1;
  double num5 = (s * s + sp * sp + t * t + tp * tp) / (u * up);

  double num6 = s_sp * fac2 - fac3 - 2. * fac1;
  double num7 = (s * s + sp * sp) * fac2;
  double den7 = (ttp * uup);
  
  // C1 = (N^2 - 1)^2 / 4N^3 = 16. / 27.
  // C2 = (N^2 - 1)   / 4N^3 =  2. / 27.
  // C3 = (N^4 - 1)   / 8N^4 = 10. / 81.
  // C4 = (N^2 - 1)^2 / 8N^4 =  8. / 81.
  return (1. / 8.) * pow3(4. * M_PI * alpS) *
         ( ( (16. / 27.) * num1 - (2. / 27.) * num2 ) * num3 +
         ( (16. / 27.) * num4 - (2. / 27.) * num2 ) * num5 +
         ( (10. / 81.) * num2 + (8. / 81.) * num6 ) *
         ( num7 / den7 ) ) / den1;

}

//--------------------------------------------------------------------------

// Evaluate |M|^2 - incoming flavour dependence

double Sigma3qq2qqgSame::sigmaHat() {
  // q q / qbar qbar incoming states only
  if (id1 != id2) return 0.;
  return sigma;
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3qq2qqgSame::setIdColAcol(){

  // Need to know where the gluon was mapped (pCM[4])
  int gIdx = 0;
  switch (config) {
  case 3: case 5: gIdx = 0; break;
  case 1: case 4: gIdx = 1; break;
  case 0: case 2: gIdx = 2; break;
  }

  // Outgoing flavours
  int idTmp[3] = { id1, id1, id1 };
  idTmp[gIdx]  = 21;
  setId(id1, id2, idTmp[0], idTmp[1], idTmp[2]);

  // Temporary solution; start with q q -> q q g
  setColAcol(1, 0, 2, 0, 1, 0, 3, 0, 2, 3);
  // Map gluon
  swap( colSave[5],  colSave[gIdx + 3]);
  swap(acolSave[5], acolSave[gIdx + 3]);
  // Swap if qbar qbar incoming
  if (id1 < 0) swapColAcol();
  
}

//--------------------------------------------------------------------------

// Map a final state configuration
inline void Sigma3qq2qqgSame::mapFinal() {
  switch (config) {
  case 0: pCM[2] = p3cm; pCM[3] = p4cm; pCM[4] = p5cm; break;
  case 1: pCM[2] = p3cm; pCM[3] = p5cm; pCM[4] = p4cm; break;
  case 2: pCM[2] = p4cm; pCM[3] = p3cm; pCM[4] = p5cm; break;
  case 3: pCM[2] = p4cm; pCM[3] = p5cm; pCM[4] = p3cm; break;
  case 4: pCM[2] = p5cm; pCM[3] = p3cm; pCM[4] = p4cm; break;
  case 5: pCM[2] = p5cm; pCM[3] = p4cm; pCM[4] = p3cm; break;
  }
}

//==========================================================================

// Sigma3qqbar2qqbargSame class.
// Cross section for q qbar -> q qbar g
// Crossed relation from q(bar) q(bar) -> q(bar) q(bar) g:
//   q1(p+) q2(p-) -> q1(q+) q2(q-) g(k)
//--------------------------------------------------------------------------

// Evaluate |M|^2 - no incoming flavour dependence

void Sigma3qqbar2qqbargSame::sigmaKin() {

  // Incoming four-vectors
  pCM[0] = Vec4( 0., 0.,  0.5 * mH, 0.5 * mH);
  pCM[1] = Vec4( 0., 0., -0.5 * mH, 0.5 * mH);

  // Pick and map a final state configuration
  pickFinal();
  mapFinal();

  // Crossing
  swap(pCM[1], pCM[3]);
  pCM[1] = -pCM[1];
  pCM[3] = -pCM[3];

  // |M|^2
  // Extra factor of (6) from picking a final state configuration
  sigma = 6. * m2Calc();

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3qqbar2qqbargSame::setIdColAcol(){
  // Outgoing flavours; easiest to map by hand
  switch (config) {
  case 0: id3 = id1; id4 = id2; id5 = 21;  break;
  case 1: id3 = id1; id4 = 21;  id5 = id2; break;
  case 2: id3 = id2; id4 = id1; id5 = 21;  break;
  case 3: id3 = 21;  id4 = id1; id5 = id2; break;
  case 4: id3 = id2; id4 = 21;  id5 = id1; break;
  case 5: id3 = 21;  id4 = id2; id5 = id1; break;
  }
  setId(id1, id2, id3, id4, id5);

  // Temporary solution; start with q qbar -> q qbar g
  int cols[5][2];
  cols[0][0] = 1; cols[0][1] = 0;
  cols[1][0] = 0; cols[1][1] = 2;
  cols[2][0] = 1; cols[2][1] = 0;
  cols[3][0] = 0; cols[3][1] = 3;
  cols[4][0] = 3; cols[4][1] = 2;
  // Map final state
  int i3 = 0, i4 = 0, i5 = 0;
  switch (config) {
  case 0: i3 = 2; i4 = 3; i5 = 4; break;
  case 1: i3 = 2; i4 = 4; i5 = 3; break;
  case 2: i3 = 3; i4 = 2; i5 = 4; break;
  case 3: i3 = 4; i4 = 2; i5 = 3; break;
  case 4: i3 = 3; i4 = 4; i5 = 2; break;
  case 5: i3 = 4; i4 = 3; i5 = 2; break;
  }
  setColAcol(cols[0][0],  cols[0][1],  cols[1][0],  cols[1][1],
             cols[i3][0], cols[i3][1], cols[i4][0], cols[i4][1],
             cols[i5][0], cols[i5][1]);
  // Swap for qbar q incoming
  if (id1 < 0) swapColAcol();
}

//==========================================================================

// Sigma3qg2qqqbarSame class.
// Cross section for q g -> q q qbar.
// Crossed relation from q(bar) q(bar) -> q(bar) q(bar) g:
//   q1(p+) q1(p-) -> q1(q+) q1(q-) g(k)

//--------------------------------------------------------------------------

// Evaluate |M|^2 - no incoming flavour dependence

void Sigma3qg2qqqbarSame::sigmaKin() {

  // Pick a final state configuration
  pickFinal();

  // gq and qg incoming
  for (int i = 0; i < 2; i++) {

    // Map incoming four-vectors
    pCM[0] = Vec4( 0., 0.,  0.5 * mH, 0.5 * mH);
    pCM[1] = Vec4( 0., 0., -0.5 * mH, 0.5 * mH);
    mapFinal();

    // Crossing (note extra -ve sign in total sigma)
    swap(pCM[i], pCM[4]);
    pCM[i] = -pCM[i];
    pCM[4] = -pCM[4];

    // |M|^2
    // XXX - Extra factor of (3) from picking a final state configuration???
    // Extra factor of (3 / 8) as averaging over incoming gluon
    sigma[i] = -3. * (3. / 8.) * m2Calc();

  }

}

//--------------------------------------------------------------------------

// Evaluate |M|^2 - incoming flavour dependence

double Sigma3qg2qqqbarSame::sigmaHat() {
  // gq or qg incoming
  return (id1 == 21) ? sigma[0] : sigma[1];
}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma3qg2qqqbarSame::setIdColAcol(){

  // Pick outgoing flavour configuration
  int idIn = (id1 == 21) ? id2 : id1;

  // Outgoing flavours; easiest just to map by hand
  switch (config) {
  case 0: id3 =  idIn; id4 =  idIn;   id5 = -idIn; break;
  case 1: id3 =  idIn; id4 = -idIn;   id5 =  idIn; break;
  case 2: id3 =  idIn; id4 =  idIn;   id5 = -idIn; break;
  case 3: id3 = -idIn; id4 =  idIn;   id5 =  idIn; break;
  case 4: id3 =  idIn; id4 = -idIn;   id5 =  idIn; break;
  case 5: id3 = -idIn; id4 =  idIn;   id5 =  idIn; break;
  }
  setId(id1, id2, id3, id4, id5);

  // Temporary solution; start with either
  // g q1    -> q1    q2 qbar2
  // g qbar1 -> qbar1 qbar2 q2
  int cols[5][2];
  cols[0][0] = 1; cols[0][1] = 2;
  if (idIn > 0) {
    cols[1][0] = 3; cols[1][1] = 0;
    cols[2][0] = 1; cols[2][1] = 0;
    cols[3][0] = 3; cols[3][1] = 0;
    cols[4][0] = 0; cols[4][1] = 2;
  } else {
    cols[1][0] = 0; cols[1][1] = 3;
    cols[2][0] = 0; cols[2][1] = 2;
    cols[3][0] = 0; cols[3][1] = 3;
    cols[4][0] = 1; cols[4][1] = 0;
  }
  // Swap incoming if q/qbar g instead
  if (id2 == 21) {
    swap(cols[0][0], cols[1][0]);
    swap(cols[0][1], cols[1][1]);
  }
  // Map final state
  int i3 = 0, i4 = 0, i5 = 0;
  switch (config) {
  case 0: i3 = 2; i4 = 3; i5 = 4; break;
  case 1: i3 = 2; i4 = 4; i5 = 3; break;
  case 2: i3 = 3; i4 = 2; i5 = 4; break;
  case 3: i3 = 4; i4 = 2; i5 = 3; break;
  case 4: i3 = 3; i4 = 4; i5 = 2; break;
  case 5: i3 = 4; i4 = 3; i5 = 2; break;
  }
  setColAcol(cols[0][0],  cols[0][1],  cols[1][0],  cols[1][1],
             cols[i3][0], cols[i3][1], cols[i4][0], cols[i4][1],
             cols[i5][0], cols[i5][1]);
}

//==========================================================================

} // end namespace Pythia8
