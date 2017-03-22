// SigmaOnia.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// charmonia/bottomonia simulation classes. 

#include "SigmaOnia.h"

namespace Pythia8 {

//==========================================================================

// Sigma2gg2QQbar3S11g class.
// Cross section g g -> QQbar[3S1(1)] g (Q = c or b).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2gg2QQbar3S11g::initProc() {

  // Produced state. Process name. Onium matrix element.
  idHad = (idNew == 4) ? 443 : 553;
  nameSave = (idNew == 4) ? "g g -> ccbar[3S1(1)] g" 
    : "g g -> bbbar[3S1(1)] g";
  oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:OJpsi3S11")
    : settingsPtr->parm("Bottomonium:OUpsilon3S11");

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence. 

void Sigma2gg2QQbar3S11g::sigmaKin() { 

  // Calculate kinematics dependence.
  double stH = sH + tH;
  double tuH = tH + uH;
  double usH = uH + sH;
  double sig = (10. * M_PI / 81.) * m3 * ( pow2(sH * tuH) 
    + pow2(tH * usH) + pow2(uH * stH) ) / pow2( stH * tuH * usH );

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2QQbar3S11g::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 21);

  // Two orientations of colour flow..
  setColAcol( 1, 2, 2, 3, 0, 0, 1, 3);
  if (rndmPtr->flat() > 0.5) swapColAcol();

}

//==========================================================================

// Sigma2gg2QQbar3PJ1g class.
// Cross section g g -> QQbar[3PJ(1)] g (Q = c or b, J = 0, 1 or 2).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2gg2QQbar3PJ1g::initProc() {

  // Produced state. Process name. Onium matrix element.
  idHad = 0;
  nameSave = "illegal process";
  if (jSave == 0) {
    idHad = (idNew == 4) ? 10441 : 10551;
    nameSave = (idNew == 4) ? "g g -> ccbar[3P0(1)] g" 
      : "g g -> bbbar[3P0(1)] g";
  } else if (jSave == 1) {
    idHad = (idNew == 4) ? 20443 : 20553;
    nameSave = (idNew == 4) ? "g g -> ccbar[3P1(1)] g" 
      : "g g -> bbbar[3P1(1)] g";
  } else if (jSave == 2) {
    idHad = (idNew == 4) ? 445 : 555;
    nameSave = (idNew == 4) ? "g g -> ccbar[3P2(1)] g" 
      : "g g -> bbbar[3P2(1)] g";
  } 
  oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:Ochic03P01")
    : settingsPtr->parm("Bottomonium:Ochib03P01");

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence. 

void Sigma2gg2QQbar3PJ1g::sigmaKin() { 

  // Useful derived kinematics quantities.
  double pRat  = (sH * uH + uH * tH + tH * sH)/ sH2;
  double qRat  = tH * uH / sH2;
  double rRat  = s3 / sH;
  double pRat2 = pRat * pRat;
  double pRat3 = pRat2 * pRat;
  double pRat4 = pRat3 * pRat;
  double qRat2 = qRat * qRat;
  double qRat3 = qRat2 * qRat;
  double qRat4 = qRat3 * qRat;
  double rRat2 = rRat * rRat;
  double rRat3 = rRat2 * rRat;
  double rRat4 = rRat3 * rRat;

  // Calculate kinematics dependence.
  double sig = 0.;
  if (jSave == 0) {
    sig = (8. * M_PI / (9. * m3 * sH)) 
      * ( 9. * rRat2 * pRat4 * (rRat4 - 2. * rRat2 * pRat + pRat2)
      - 6. * rRat * pRat3 * qRat * (2. * rRat4 - 5. * rRat2 * pRat 
      + pRat2) - pRat2 * qRat2 * (rRat4 + 2. * rRat2 * pRat - pRat2)
      + 2. * rRat * pRat * qRat3 * (rRat2 - pRat) + 6. * rRat2 * qRat4)
      / (qRat * pow4(qRat - rRat * pRat));
  } else if (jSave == 1) {
    sig =  (8. * M_PI / (3.* m3 * sH)) * pRat2 
      * (rRat * pRat2 * (rRat2 - 4. * pRat)
      + 2. * qRat * (-rRat4 + 5. * rRat2 * pRat + pRat2)
      - 15. * rRat * qRat2) / pow4(qRat - rRat * pRat);
  } else if (jSave == 2) {
    sig = (8. * M_PI / (9. * m3 * sH)) 
      * (12. * rRat2 * pRat4 * (rRat4 - 2. * rRat2 * pRat + pRat2)
      - 3. * rRat * pRat3 * qRat * (8. * rRat4 - rRat2 * pRat + 4. * pRat2)
      + 2. * pRat2 * qRat2 * (-7. * rRat4 + 43. * rRat2 * pRat + pRat2)
      + rRat * pRat * qRat3 * (16. * rRat2 - 61. * pRat)
      + 12. * rRat2 * qRat4) / (qRat * pow4(qRat-rRat * pRat));
  }

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2QQbar3PJ1g::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 21);

  // Two orientations of colour flow.
  setColAcol( 1, 2, 2, 3, 0, 0, 1, 3);
  if (rndmPtr->flat() > 0.5) swapColAcol();

}

//==========================================================================

// Sigma2qg2QQbar3PJ1q class.
// Cross section q g -> QQbar[3PJ(1)] q (Q = c or b, J = 0, 1 or 2).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qg2QQbar3PJ1q::initProc() {

  // Produced state. Process name. Onium matrix element.
  idHad = 0;
  nameSave = "illegal process";
  if (jSave == 0) {
    idHad = (idNew == 4) ? 10441 : 10551;
    nameSave = (idNew == 4) ? "q g -> ccbar[3P0(1)] q" 
      : "q g -> bbbar[3P0(1)] q";
  } else if (jSave == 1) {
    idHad = (idNew == 4) ? 20443 : 20553;
    nameSave = (idNew == 4) ? "q g -> ccbar[3P1(1)] q" 
      : "q g -> bbbar[3P1(1)] q";
  } else if (jSave == 2) {
   idHad = (idNew == 4) ? 445 : 555;
    nameSave = (idNew == 4) ? "q g -> ccbar[3P2(1)] q" 
      : "q g -> bbbar[3P2(1)] q";
  } 
  oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:Ochic03P01")
    : settingsPtr->parm("Bottomonium:Ochib03P01");

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence. 

void Sigma2qg2QQbar3PJ1q::sigmaKin() { 

  // Calculate kinematics dependence.
  double usH = uH + sH;
  double sig = 0.;
  if (jSave == 0) {
    sig = - (16. * M_PI / 81.) * pow2(tH - 3. * s3) * (sH2 + uH2)
      / (m3 * tH * pow4(usH));
  } else if (jSave == 1) {
    sig = - (32. * M_PI / 27.) * (4. * s3 * sH * uH + tH * (sH2 + uH2))
      / (m3 * pow4(usH));
  } else if (jSave == 2) {
    sig = - (32. *M_PI / 81.) * ( (6. * s3*s3 + tH2) * pow2(usH)
      - 2. * sH * uH * (tH2 + 6. * s3 * usH)) / (m3 * tH * pow4(usH));
  }

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2QQbar3PJ1q::setIdColAcol() {

  // Flavours are trivial.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, idHad, idq);

  // tH defined between q_in and q_out: must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21); 

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2QQbar3PJ1g class.
// Cross section q qbar -> QQbar[3PJ(1)] g (Q = c or b, J = 0, 1 or 2).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qqbar2QQbar3PJ1g::initProc() {

  // Produced state. Process name. Onium matrix element.
  idHad = 0;
  nameSave = "illegal process";
  if (jSave == 0) {
    idHad = (idNew == 4) ? 10441 : 10551;
    nameSave = (idNew == 4) ? "q qbar -> ccbar[3P0(1)] g" 
      : "q qbar -> bbbar[3P0(1)] g";
  } else if (jSave == 1) {
    idHad = (idNew == 4) ? 20443 : 20553;
    nameSave = (idNew == 4) ? "q qbar -> ccbar[3P1(1)] g" 
      : "q qbar -> bbbar[3P1(1)] g";
  } else if (jSave == 2) {
   idHad = (idNew == 4) ? 445 : 555;
    nameSave = (idNew == 4) ? "q qbar -> ccbar[3P2(1)] g" 
      : "q qbar -> bbbar[3P2(1)] g";
  } 
  oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:Ochic03P01")
    : settingsPtr->parm("Bottomonium:Ochib03P01");

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence. 

void Sigma2qqbar2QQbar3PJ1g::sigmaKin() { 

  // Calculate kinematics dependence.
  double tuH = tH + uH;
  double sig = 0.;
  if (jSave == 0) {
    sig =(128. * M_PI / 243.) * pow2(sH - 3. * s3) * (tH2 + uH2)
      / (m3 * sH * pow4(tuH));
  } else if (jSave == 1) {
    sig = (256. * M_PI / 81.) * (4. * s3 * tH * uH + sH * (tH2 + uH2))
      / (m3 * pow4(tuH));
  } else if (jSave == 2) {
    sig = (256. * M_PI / 243.) * ( (6. * s3*s3 + sH2) * pow2(tuH)
      - 2. * tH * uH * (sH2 + 6. * s3 * tuH) )/ (m3 * sH * pow4(tuH));
  }

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2QQbar3PJ1g::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 21);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2gg2QQbarX8g class.
// Cross section g g -> QQbar[X(8)] g (Q = c or b, X = 3S1, 1S0 or 3PJ).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2gg2QQbarX8g::initProc() {

  // Produced state. Process name. Onium matrix element.
  idHad = 0;
  nameSave = "illegal process";
  if (stateSave == 0) {
    idHad = (idNew == 4) ? 9900443 : 9900553;
    nameSave = (idNew == 4) ? "g g -> ccbar[3S1(8)] g" 
      : "g g -> bbbar[3S1(8)] g";
    oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:OJpsi3S18")
      : settingsPtr->parm("Bottomonium:OUpsilon3S18");
  } else if (stateSave == 1) {
    idHad = (idNew == 4) ? 9900441 : 9900551;
    nameSave = (idNew == 4) ? "g g -> ccbar[1S0(8)] g" 
      : "g g -> bbbar[1S0(8)] g";
    oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:OJpsi1S08")
      : settingsPtr->parm("Bottomonium:OUpsilon1S08");
  } else if (stateSave == 2) {
    idHad = (idNew == 4) ? 9910441 : 9910551;
    nameSave = (idNew == 4) ? "g g -> ccbar[3PJ(8)] g" 
      : "g g -> bbbar[3PJ(8)] g";
    oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:OJpsi3P08")
      : settingsPtr->parm("Bottomonium:OUpsilon3P08");
  } 

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence. 

void Sigma2gg2QQbarX8g::sigmaKin() { 

  // Calculate kinematics dependence.
  double stH = sH + tH;
  double tuH = tH + uH;
  double usH = uH + sH;
  double sig = 0.;
  if (stateSave == 0) {
    sig = (M_PI / 72.) * m3 * ( 27. * (pow2(stH) + pow2(tuH)
      + pow2(usH)) / (s3*s3) - 16. ) * ( pow2(sH * tuH) 
      + pow2(tH * usH) + pow2(uH * stH) ) / pow2( stH * tuH * usH );
  } else if (stateSave == 1) {
    sig = (5. * M_PI / 16.) * m3 * ( pow2(uH / (tuH * usH))  
      + pow2(sH / (stH * usH)) + pow2(tH / (stH * tuH)) ) * ( 12. 
      + (pow4(stH) + pow4(tuH) + pow4(usH)) / (s3 * sH * tH * uH) ); 
  } else if (stateSave == 2) {
    double sH3 = sH2 * sH;
    double sH4 = sH3 * sH;
    double sH5 = sH4 * sH;
    double sH6 = sH5 * sH;
    double sH7 = sH6 * sH;
    double sH8 = sH7 * sH;
    double tH3 = tH2 * tH;
    double tH4 = tH3 * tH;
    double tH5 = tH4 * tH;
    double tH6 = tH5 * tH;
    double tH7 = tH6 * tH;
    double tH8 = tH7 * tH;
    double ssttH = sH * sH + sH * tH + tH * tH;
    sig = 5. * M_PI * (3. * sH * tH * stH * pow4(ssttH)
      - s3 * pow2(ssttH) * (7. * sH6 + 36. * sH5 * tH + 45. * sH4 * tH2
        + 28. * sH3 * tH3 + 45. * sH2 * tH4 + 36. * sH * tH5 + 7. * tH6)
      + pow2(s3) * stH * (35. *sH8 + 169. * sH7 * tH + 299. * sH6 * tH2
        + 401. * sH5 * tH3 + 418. * sH4 * tH4 + 401. * sH3 * tH5 
        + 299. * sH2 * tH6 + 169. * sH * tH7 + 35. * tH8)
      - pow3(s3) * (84. *sH8+432. *sH7*tH+905. *sH6*tH2 
        + 1287. * sH5 * tH3 + 1436. * sH4 * tH4 +1287. * sH3 * tH5 
        + 905. * sH2 * tH6 + 432. * sH * tH7 + 84. * tH8)
      + pow4(s3) * stH * (126. * sH6 + 451. * sH5 * tH +677. * sH4 * tH2
        + 836. * sH3 * tH3 + 677. * sH2 * tH4 + 451. * sH * tH5 
        + 126. * tH6)
      - pow5(s3) * 3. * (42. * sH6 + 171. * sH5 * tH + 304. * sH4 * tH2
        + 362. * sH3 * tH3 + 304. * sH2 * tH4 + 171. * sH * tH5 
        + 42. * tH6)
      + pow3(s3 * s3) * 2. * stH * (42. * sH4 + 106. * sH3 * tH 
        + 119. * sH2 * tH2 + 106. * sH * tH3 + 42. * tH4)
      - pow4(s3) * pow3(s3) * (35. * sH4 + 99. * sH3 * tH 
        + 120. * sH2 * tH2 + 99.  * sH * tH3 + 35.  * tH4)
      + pow4(s3 * s3) * 7. * stH * ssttH)
      / (sH * tH * uH * s3 * m3 * pow3(stH * tuH * usH));
  } 

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2QQbarX8g::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 21);

  // Split total contribution into different colour flows just like in
  // g g -> g g (with kinematics recalculated for massless partons).
  double sHr    = - (tH + uH);
  double sH2r   = sHr * sHr;
  double sigTS  = tH2/sH2r + 2.*tH/sHr + 3. + 2.*sHr/tH + sH2r/tH2;
  double sigUS  = uH2/sH2r + 2.*uH/sHr + 3. + 2.*sHr/uH + sH2r/uH2;
  double sigTU  = tH2/uH2 + 2.*tH/uH + 3. + 2.*uH/tH + uH2/tH2;
  double sigSum = sigTS + sigUS + sigTU;

  // Three colour flow topologies, each with two orientations.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 4, 4, 3);
  else if (sigRand < sigTS + sigUS) 
                       setColAcol( 1, 2, 3, 1, 3, 4, 4, 2);
  else                 setColAcol( 1, 2, 3, 4, 1, 4, 3, 2); 
  if (rndmPtr->flat() > 0.5) swapColAcol();


}

//==========================================================================

// Sigma2qg2QQbarX8q class.
// Cross section q g -> QQbar[X(8)] q (Q = c or b, X = 3S1, 1S0 or 3PJ).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qg2QQbarX8q::initProc() {

  // Produced state. Process name. Onium matrix element.
  idHad = 0;
  nameSave = "illegal process";
  if (stateSave == 0) {
    idHad = (idNew == 4) ? 9900443 : 9900553;
    nameSave = (idNew == 4) ? "q g -> ccbar[3S1(8)] q" 
      : "q g -> bbbar[3S1(8)] q";
    oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:OJpsi3S18")
      : settingsPtr->parm("Bottomonium:OUpsilon3S18");
  } else if (stateSave == 1) {
    idHad = (idNew == 4) ? 9900441 : 9900551;
    nameSave = (idNew == 4) ? "q g -> ccbar[1S0(8)] q" 
      : "q g -> bbbar[1S0(8)] q";
    oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:OJpsi1S08")
      : settingsPtr->parm("Bottomonium:OUpsilon1S08");
  } else if (stateSave == 2) {
    idHad = (idNew == 4) ? 9910441 : 9910551;
    nameSave = (idNew == 4) ? "q g -> ccbar[3PJ(8)] q" 
      : "q g -> bbbar[3PJ(8)] q";
    oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:OJpsi3P08")
      : settingsPtr->parm("Bottomonium:OUpsilon3P08");
  } 

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence. 

void Sigma2qg2QQbarX8q::sigmaKin() { 

  // Calculate kinematics dependence.
  double stH  = sH + tH;
  double tuH  = tH + uH;
  double usH  = uH + sH;
  double stH2 = stH * stH;
  double tuH2 = tuH * tuH;
  double usH2 = usH * usH;
  double sig  = 0.;
  if (stateSave == 0) {
    sig = - (M_PI / 27.)* (4. * (sH2 + uH2) - sH * uH) * (stH2 +tuH2)
      / (s3 * m3 * sH * uH * usH2);
  } else if (stateSave == 1) {
    sig = - (5. * M_PI / 18.) * (sH2 + uH2) / (m3 * tH * usH2);
  } else if (stateSave == 2) {
    sig = - (10. * M_PI / 9.) * ( (7. * usH + 8. * tH) * (sH2 + uH2)
      + 4. * tH * (2. * pow2(s3) - stH2 - tuH2) ) 
      / (s3 * m3 * tH * usH2 * usH);
  } 

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2QQbarX8q::setIdColAcol() {

  // Flavours are trivial.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, idHad, idq);

  // tH defined between q_in and q_out: must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21); 

  // Split total contribution into different colour flows just like in
  // q g -> q g (with kinematics recalculated for massless partons).
  double sHr    = - (tH + uH);
  double sH2r   = sHr * sHr;
  double sigTS  = uH2/tH2 - (4./9.) * uH/sHr;
  double sigTU  = sH2r/tH2 - (4./9.) * sHr/uH;
  double sigSum = sigTS + sigTU;

  // Two colour flow topologies. Swap if first is gluon, or when antiquark.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 2, 1, 2, 3, 3, 0);
  else                 setColAcol( 1, 0, 2, 3, 1, 3, 2, 0); 
  if (id1 == 21) swapCol12();
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2QQbarX8g class.
// Cross section q qbar -> QQbar[X(8)] g (Q = c or b, X = 3S1, 1S0 or 3PJ).

//--------------------------------------------------------------------------

// Initialize process. 
  
void Sigma2qqbar2QQbarX8g::initProc() {

  // Produced state. Process name. Onium matrix element.
  idHad = 0;
  nameSave = "illegal process";
  if (stateSave == 0) {
    idHad = (idNew == 4) ? 9900443 : 9900553;
    nameSave = (idNew == 4) ? "q qbar -> ccbar[3S1(8)] g" 
      : "q qbar -> bbbar[3S1(8)] g";
    oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:OJpsi3S18")
      : settingsPtr->parm("Bottomonium:OUpsilon3S18");
  } else if (stateSave == 1) {
    idHad = (idNew == 4) ? 9900441 : 9900551;
    nameSave = (idNew == 4) ? "q qbar -> ccbar[1S0(8)] g" 
      : "q qbar -> bbbar[1S0(8)] g";
    oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:OJpsi1S08")
      : settingsPtr->parm("Bottomonium:OUpsilon1S08");
  } else if (stateSave == 2) {
    idHad = (idNew == 4) ? 9910441 : 9910551;
    nameSave = (idNew == 4) ? "q qbar -> ccbar[3PJ(8)] g" 
      : "q qbar -> bbbar[3PJ(8)] g";
    oniumME = (idNew == 4) ? settingsPtr->parm("Charmonium:OJpsi3P08")
      : settingsPtr->parm("Bottomonium:OUpsilon3P08");
  } 

} 

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence. 

void Sigma2qqbar2QQbarX8g::sigmaKin() { 

  // Calculate kinematics dependence.
  double stH  = sH + tH;
  double tuH  = tH + uH;
  double usH  = uH + sH;
  double stH2 = stH * stH;
  double tuH2 = tuH * tuH;
  double usH2 = usH * usH;
  double sig  = 0.;
  if (stateSave == 0) {
    sig = (8. * M_PI / 81.) * (4. * (tH2 + uH2) - tH * uH) 
      * (stH2 + usH2) / (s3 * m3 * tH * uH * tuH2);
  } else if (stateSave == 1) {
    sig = (20. * M_PI / 27.) * (tH2 + uH2) / (m3 * sH * tuH2);
  } else if (stateSave == 2) {
    sig = (80. * M_PI / 27.) * ( (7. * tuH + 8. * sH) * (tH2 + uH2)
      + 4. * sH * (2. * pow2(s3) - stH2 -usH2) ) 
      / (s3 * m3 * sH * tuH2 * tuH);
  } 

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;  

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2QQbarX8g::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 21);

  // Split total contribution into different colour flows just like in
  // q qbar -> g g (with kinematics recalculated for massless partons).
  double sHr    = - (tH + uH);
  double sH2r   = sHr * sHr;
  double sigTS  = (4. / 9.) * uH / tH - uH2 / sH2r;
  double sigUS  = (4. / 9.) * tH / uH - tH2 / sH2r;
  double sigSum = sigTS + sigUS;

  // Two colour flow topologies. Swap if first is antiquark.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 0, 2, 1, 3, 3, 2);
  else                 setColAcol( 1, 0, 0, 2, 3, 2, 1, 3); 
  if (id1 < 0) swapColAcol();

}

//==========================================================================

} // end namespace Pythia8

