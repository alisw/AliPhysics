// SimpleWeakShowerMEs.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// SimpleWeakShowerMEs class.

#include "Pythia8/SimpleWeakShowerMEs.h"

namespace Pythia8 {

//==========================================================================

// The SimpleWeakShowerMEs class.
// The 2 -> 2 MEs contain the correct kinematics, but for some of
// the 2 -> 3 MEs some couplings have been switched off.
// This class is used for ME corrections in the weak shower and for merging
// whenever weak reclusterings are allowed.

//--------------------------------------------------------------------------

// Calculate the 2->2 ME qg -> qg, omitting an overall factor of g_s^4 / 9.

double SimpleWeakShowerMEs::getMEqg2qg(double sH, double tH, double uH) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;
  return (sH2 + uH2) * (9. / tH2 - 4. / (sH * uH));

}

//--------------------------------------------------------------------------

// Calculate the 2->2 ME qq' -> qq', omitting an overall factor of g_s^4 / 9s.

double SimpleWeakShowerMEs::getMEqq2qq(double sH, double tH, double uH,
  bool sameID) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;
  if (sameID) return 2. * ( (sH2 + uH2) / tH2 + (sH2 + tH2) / uH2
    - 2. * sH2 / (3. * tH * uH) );
  else return 4. * (sH2 + uH2) / tH2;

}

//--------------------------------------------------------------------------

// Calculate the 2->2 ME gg -> gg, omitting an overall factor of g_s^4 / 9.

double SimpleWeakShowerMEs::getMEgg2gg(double sH, double tH, double uH) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;
  return (81. / 8.) * ( (tH2 + uH2) / sH2 + (sH2 + uH2) / tH2
    + (sH2 + tH2) / uH2 + 3. );

}

//--------------------------------------------------------------------------

// Calculate the 2->2 ME gg -> qqbar, omitting an overall factor of g_s^4 / 9.

double SimpleWeakShowerMEs::getMEgg2qqbar(double sH, double tH, double uH) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;
  return (tH2 + uH2) * ( 3. / (2. * tH * uH) - 27. / (8. * sH2) );

}

//--------------------------------------------------------------------------

// Calculate the 2->2 ME qqbar -> gg, omitting an overall factor of g_s^4 / 9.

double SimpleWeakShowerMEs::getMEqqbar2gg(double sH, double tH, double uH) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;
  return (tH2 + uH2) * ( 16. / (3. * tH * uH) - 12. / sH2 );

}

//--------------------------------------------------------------------------

// Calculate the 2->2 ME qqbar -> qqbar, omitting an overall factor of
// g_s^4 / 9.

double SimpleWeakShowerMEs::getMEqqbar2qqbar(double sH, double tH, double uH,
  bool sameID) {

  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;
  if (sameID) return 4. * (tH2 + uH2) / sH2 - (8./3.) * uH2 / (sH * tH)
    + 4. * (sH2 + uH2) / tH2;
  else return 4. * (tH2 + uH2) / sH2;

}

//--------------------------------------------------------------------------

// Calculate the 2->3 ME ug -> ugZ, omitting an overall factor of
// \frac{(g_s^4 * 4\pi \alpha_{em} * (lU^2 + rU^2))}
//      {9 * \cos^2(\theta_W) \sin^2(\theta_W)}.
// p1 = incoming quark, p2 = incoming gluon, p3 = outgoing gluon,
// p4 = outgoing Z, p5 = outgoing quark.

double SimpleWeakShowerMEs::getMEqg2qgZ(Vec4 p1, Vec4 p2, Vec4 p3, Vec4 p4,
  Vec4 p5) {

  double p12  = p1*p2;
  double p13  = p1*p3;
  double p14  = p1*p4;
  double p23  = p2*p3;
  double p24  = p2*p4;
  double p44  = p4*p4;
  double p12P = (p1 + p2).m2Calc();
  double p13M = (p1 - p3).m2Calc();
  double p14M = (p1 - p4).m2Calc();
  double p23M = (p2 - p3).m2Calc();
  double p25M = (p2 - p5).m2Calc();
  double p35P = (p3 + p5).m2Calc();
  double p45P = (p4 + p5).m2Calc();
  double p44S = p44*p44;

  double dia1 = (-4*p12*(p12* (p44 - 2*p24) -p23 * (p44 - 2*p24) +2*p13*p24))
    / pow2(p12P * p45P);

  double dia2 = -(((p12 - p13 - p23) * (p12 * (p44 - 2*p14)
    - p13 * (p44 - 2*p14) + 2*p14*p23)) / (p12P * p13M * pow2(p45P)));

  double dia3 = (8*(p44 - 2*p12) * p12 * (p12 - p23 - p24))
    / (pow2(p12P) * p35P * p45P);

  double dia4 = -(8*pow3(p12) + 2*pow2(p12) * (p44 - 6*p13 - 6*p14 - 6*p23
    - 4*p24) + 4*pow2(p13) * (p44 - 2*p14 - p24) - 2*p14*p23
    * (-p44 + 2*p14 + 2*p23 + 4*p24) + 2*p13 * (-4*pow2(p14)
    + 2*p14 * (p44 - 3*p23 - 3*p24) + 2*p23 * (p44 - p24) + p44*p24)
    + p12* (p44S + 4*pow2(p13) + 4*pow2(p14) + 4*pow2(p23) - 4*p14
    * (p44 - 4*p23 - 3*p24) - 4*p13* (p44 - 4*p14 - 4*p23 - 3*p24)
    - 2*p44*p24 + 4*p23*p24)) / (4. * p12P * p13M * p25M * p45P);

  double dia5 = (-9*(2*pow3(p12) + (p13 + p23) * (-(p14* p23)
    + p13 * (p44 - 2*p14 - p24)) + pow2(p12)* (2*p44 - 4*p13 - 3*p14
    - 4*p23 - 3*p24) + p12 * (2*pow2(p13) + p23 * (-p44 + 4*p14
    + 2*p23 + 3*p24) + p13 * (-2*p44 + 5*p14 + 4*p23 + 4*p24))))
    / (p12P * p23M* pow2(p45P));

  double dia6 = (-4*(2*pow3(p12) + pow2(p12)* (7*p44 - 4*p13 - 6*p14
    - 2*p23 - 4*p24) + p12 * (2*pow2(p13) + 4*pow2(p14) - 2*p44*p23
    - 3*p44*p24 + 2*p23* p24 + 2*pow2(p24) + p14 * (-8*p44 + 6*p23 + 2*p24)
    + p13* (-7*p44 + 10*p14 + 6*p23 + 6*p24)) + 2*(pow2(p13) * (p44
    - 2*p14 - p24) + p14* (p14* (p44 - 2*p23) + p23 * (p44 - p24) + p44*p24)
    + p13 * (-2*pow2(p14) + p23 * (p44 - 2*p24) + p14* (2*p44 - 3*p23
    - 2*p24) + (p44 - p24)* p24)))) / (p12P * p14M * p35P * p45P);

  double dia7 = ((p12 - p13 - p14) * (p44 - 2*p23) * (p12 - p13 - p14
    - p23 - p24)) / (p12P * p14M * p25M * p45P);

  double dia8 = (9*(2*pow2((p1*p2))* (p44 - 10*p14 - 16*p23 + 4*p24)
    + p12* (p44S + 20*pow2(p14) - 6*p44*p23 + 16*pow2(p23)
    + p13* (-18*p44 + 28*p14 + 32*p23 - 8*p24) + 2*p44*p24
    + 16*p23* p24 - 8*pow2(p24) + p14* (-8*p44 + 52*p23 + 4*p24))
    + 2*(4*pow2(p13)* (p44 - p14 - 2*p23) - p14*p23* (-3*p44 + 10*p14
    + 8*p23 + 8*p24) + p13* (-4*pow2(p14) - 8*pow2(p23) + 4*p23* (p44 - p24)
    + 2*p14* (2*p44 - 8*p23 + p24) + p24* (p44 + 4*p24)))))
    / (8. * p12P * p23M * p14M * p45P);

  double dia9 = (-4*p13* (2*pow2(p12) + p12* (p44 - 4*p13 - 2*p14
    - 4*p23 - 2*p24) + 2*(p13 + p23)* (p13 + p14 + p23 + p24)))
    / (pow2(p13M) * pow2(p45P));

  double dia10 = (4*pow3(p12) - 4*pow2(p12) * (3*p13 + 2*p14 + 2*p23 + p24)
    + p12* (p44S + 8*pow2(p13) + 4*pow2(p14) + 4*pow2(p23)
    - 4*p14* (p44 - 3*p23 - 3*p24) - 2*p44*p24 + 4*p23*p24
    + p13* (-6*p44 + 20*p14 + 4*p23 + 12*p24)) + 2*(p14*p23* (p44
    - 2*p14 - 2*p23 - 4*p24) + 2*pow2(p13)* (p44 - 2*p14 - 2*p24)
    + p13* (-4*pow2(p14) + 2*p14* (p44 - 2*p23 - 3*p24)
    + 2*p23* (p44 - p24) + p44*p24))) / (4. * p12P * p13M * p35P* p45P);

  double dia11 = (4*p13* (p44 + 2*p13)* (p44 + 2*p12 - 2*p14 - 2*p24))
   / (pow2(p13M)* p25M * p45P);

  double dia12 = (-9*(-2*pow3(p12) + 2*pow3(p13) - 2*p14* pow2(p23)
    + pow2(p12)* (p44 + 6*p13 - 2*p14 + 4*p23 + 2*p24) + p13*p23
    * (-p44 - 2*p14 + 2*p23 + 6*p24) + pow2(p13)* (p44 + 4*p23 + 6*p24)
    - p12* (6*pow2(p13) + p23* (p44 - 4*p14 + 2*p23 + 2*p24)
    + p13* (-2*p14 + 8*(p23 + p24))))) / (2.* p13M * p23M * pow2(p45P));

  double dia13 = -((p12 - p13 - p14)* (p44 - 2*p23)* (p44 - 2*p24))
    / (2. * p13M * p14M * p35P * p45P);

  double dia14 = (2*(-8*pow2(p13)* (p44 - p14) - 2*pow2(p12)
    * (p44 - 2*p14 + 2*p23 - 2*p24) + 2*p14* (-p44S - 2*pow2(p23)
    + p23 * (p44 + 2*p14 - 2*p24) + 2*p44*p24) + p12* (p44S
    - 4*pow2(p14) + 4*pow2(p23) + 2*p13* (3*p44 - 6*p14 - 2*p23 - 2*p24)
    - 4*pow2(p24)) + p13* (-2*p44S + 8*pow2(p14) + 2*(p44 + 2*p23)* p24
    + 4*pow2(p24) + p14* (-8*p44 + 8*p23 + 4*p24))))
    / (p13M * p14M * p25M * p45P);

  double dia15 = (-9*(8*pow3(p12) + 2*pow2(p12)* (p44 - 8*p13 - 2*p14
    - 4*p23 - 8*p24) + p12* (3*p44S + 8*pow2(p13) - 4*pow2(p14)
    - 6*p44*p23 - 10*p44*p24 + 24*p23* p24 + 8*pow2(p24) + 4*p14
    * (p23 + 3*p24) + 2*p13 * (3*p44 - 2*p14 + 8*p23 + 12*p24))
    + 2*(p14*p23 * (p44 + 2*p14 - 8*p24) + 4*pow2(p13)* (p14
    + p23 - p24) + p13* (4*pow2(p14) - 4*pow2(p23)
    + 4*p23* (p44 - 4*p24) - 10*p14* p24 + (3*p44 - 4*p24)* p24))))
    / (8. * p13M * p23M * p14M * p45P);

  double dia16 = (4*p12*p23* (p44 + 2*p12 - 2*p14 - 2*p24))
    / (pow2(p12P) * pow2(p35P));

  double dia17 = ((p12 - p13 - p14)* (p44 - 2*p24)
    * (p12 - p13 - p14 - p23 - p24)) / (p12P * p13M * p25M * p35P);

  double dia18 = (9*(-20*pow3(p12) + 4*pow2(p12)* (4*p44 + p13 + 2*p14
    + p23) + 2*p14*p23* (p44 + 2*p14 - 4*p23 - 2*p24) + 8*pow2(p13)
    * (p44 - p14 - p24) + p13* (-8*pow2(p14) + 8*p44*p23
    + 4*p14* (2*p44 - 2*p23 - 5*p24) + 6*p44*p24 - 4*pow2(p24))
    + p12* (p44S - 4*pow2(p14) - 14*p44*p23 + 8*pow2(p23)
    - 8*p14* (p44 - p23 - 2*p24) - 12*p44*p24 + 12*p23* p24
    + 4*pow2(p24) + p13 * (-22*p44 + 20*p14 - 8*p23 + 24*p24))))
    / (8. * p12P * p23M * p35P * p45P);

  double dia19 = (-8*p13* (p44 + 2*p12 - 2*p14 - 2*p24)*(p12 - p14 - p24))
    / (p12P * p14M * pow2(p35P));

  double dia20 = (-8*pow3(p12) - p12* (4*pow2(p13) + 4*pow2(p14)
    - 4*p14* (p44 - 3*p23 - 4*p24) - 4*p13* (p44 - 4*p14 - 3*p23 - 4*p24)
    + (p44 - 2*p24)* (p44 - 2*p23 - 2*p24)) - 2*pow2(p12)
    * (p44 - 6*p13 - 6*p14 - 4*p23 - 6*p24) + 2*p14*p23
    * (-p44 + 2*p14 + 2*p24) + pow2(p13) * (-4*p44 + 8*p14 + 4*p24)
    + p13* (8*pow2(p14) - 4*p14* (p44 - 3*p23 - 3*p24) - 2*(p44 - 2*p24)
    * (2*p23 + p24))) / (4. * p12P * p14M * p25M* p35P);

  double dia21 = (-9*(14*pow3(p12) + pow2(p12)* (7*p44 - 10*p13 - 30*p14
    - 10*p23 - 22*p24) + 2*p12 * (2*pow2(p13) + 12*pow2(p14) - p44*p23
    - 3*p44*p24 + 2*p23* p24 + 4*pow2(p24) + p13* (-4*p44 + 10*p14
    + 6*p23 + 5*p24) + p14* (-5*p44 + 7*p23 + 16*p24)) + 2*(2*pow2(p13)
    * (p44 - 2*p14 - p24) + p13* (-4*pow2(p14) + 2*p23* (p44 - 2*p24)
    + 2*p14* (p44 - 3*p23 - 2*p24) + p44*p24) + p14* (-4*pow2(p14)
    + 2*p14* (p44 - p23 - 3*p24) + p23* (p44 - 2*p24)
    + 2*(p44 - p24)* p24)))) / (4. * p12P * p23M * p14M* p35P);

  double dia22 = (-8*p13*p23* (-p12 + p23 + p24))
    / (pow2(p13M)* pow2(p25M));

  double dia23 = (-9*(2*pow3(p12) - pow2(p12)* (p44 + 10*p13 + 6*p14)
    + p12* (p44S + 8*pow2(p13) + 4*pow2(p14) + 3*p44*p23 - 2*pow2(p23)
    - p44*p24 - 4*p23 * p24 - 2*pow2(p24) + 2*p13* (6*p14 + 7*p23 + 2*p24)
    + p14* (-2*p44 + 4*p23 + 6*p24)) - 2*(4*pow3(p13) + 2*pow2(p13)
    * (p44 + 2*p14 + p23 - p24) + p14*p23* (2*p14 - p23 + p24)
    + p13 * (p44S + 4*pow2(p14) + (-2*p44 + p23)* p24 - pow2(p24)
    + p14* (-2*p44 + 4*p23 + 2*p24))))) / (4. * p13M * p23M * p25M * p45P);

  double dia24 = (-2*pow2(p12)* (p44 + p13 + p23 - p24) + p12* (2*pow2(p13)
    - p13* (p44 - 6*p14) + (p44 + 2*p14 + 2*p23 - 2*p24)* (p23 + p24))
    + 2*(pow2(p13)* (p44 - 2*p14 - p24) - p14*p23 * (p23 + p24)
    + p13* (-2*pow2(p14) + p14 * (p44 - p23 - 2*p24) + p24* (p23 + p24))))
    / (2. * p13M * p14M * p25M* p35P);

  double dia25 = (8*p12* (p44 + 2*p12 - 2*p23 - 2*p24)*(p12 - p23 - p24))
    / (p13M * p14M * pow2(p25M));

  double dia26 = (9*(4*pow3(p12) + 2*pow2(p13)* (p44 - 2*p14 - 3*p24)
    - 2*pow2(p12)* (p13 + 6*p23 + 2*p24) + p12* (p44S + 6*pow2(p13)
    - 2*p44*p14 + 4*pow2(p14) - 2*p44*p23 + 8*pow2(p23) + 3*p13
    * (p44 + 2*p14 - 2*p23 - 2*p24) - 2*p44*p24 + 8*p23*p24)
    + p14* (-p44S + 2*p14* (p44 - 2*p23 - 2*p24) + 4*p23* p24
    + 4*pow2(p24)) - p13 * (p44S + 4*pow2(p14) - 4*pow2(p23)
    + 2*p23 * (p44 - 6*p24) + 2*p44*p24 - 8*pow2(p24)
    + p14* (-4*p44 + 2*p23 + 8*p24)))) / (4. * p13M * p23M * p14M * p25M);

  double dia27 = (-9*(3*pow3(p12) - pow3(p13) + p14 * pow2(p23)
    + pow2(p12)* (p44 - 7*p13 - 3*p14 - 6*p23 - 2*p24) - p13*p23
    * (2*p44 + p23 - 2*p24) - pow2(p13)* (p14 + 2*p23 - 2*p24)
    + p12 * (5*pow2(p13) + p13* (p44 + 4*p14 + 8*p23)
    + p23 * (p44 + 2*p14 + 3*p23 + 2*p24)))) / (pow2(p23M) * pow2(p45P));

  double dia28 = (9*(2*pow3(p12) - pow2(p12)* (p44 + 4*p13 + 10*p14
    + 8*p23 - 2*p24) + 2*pow2(p13)* (p44 - 2*p14 - p24) + p14* (p44S
    + 2*p44*p23 + 2*pow2(p23) - 2*p14* (p44 + 6*p23) - 4*p44*p24)
    + p12* (-p44S + 2*pow2(p13) + 8*pow2(p14) + p13 * (-5*p44 + 14*p14)
    - 3*p44*p23 - 2*pow2(p23) + 2*p14* (p44 + 10*p23) + 6*p44*p24
    + 2*p23* p24 - 4*pow2(p24)) + p13* (-4*pow2(p14) + pow2(p44 - 2*p24)
    + 2*p23* (p44 + p24) + p14* (-6*p23 + 4*p24))))
    / (4. * p23M * p14M * p35P * p45P);

  double dia29 = (-9*(6*pow3(p12) - 2*pow2(p13)* (p44 - 2*p14 + 3*p24)
    - pow2(p12)* (p44 + 12*p13 + 2*p14 + 8*p23 + 10*p24) + p14* (p44S
    - 2*p14* (p44 - 6*p23) - 6*p44*p23 - 2*pow2(p23) - 4*p44*p24)
    + p13* (p44S + 4*pow2(p14) - 8*pow2(p23) - 2*p44* p24 - 4*pow2(p24)
    - 2*p14* (2*p44 - 5*p23 + 4*p24) - 2*p23* (2*p44 + 5*p24))
    + p12* (6* pow2(p13) - 4*pow2(p14) + 3*p44*p23 + 2*pow2(p23)
    + 6*p23* p24 + 4*pow2(p24) + p14 * (6*p44 - 4*p23 + 4*p24)
    + p13* (7*p44 - 2*p14 + 16*p23 + 16*p24))))
    / (4. * p23M * p14M * p25M * p45P);

  double dia30 = (-9*(4*pow3(p12) + 2*pow2(p12)* (p44 - 4*p13 + 2*p14
    - 6*p24) + p12* (p44S + 4*pow2(p13) - 8*pow2(p14) - 6*p44*p23
    + 12*pow2(p23) - 6*p44 *p24 + 4*p23* p24 + 8*pow2(p24) + 2*p14
    * (p44 - 4*p23 + 2*p24) + 8*p13* (p44 - p14 + 2*p23 + 2*p24))
    - 2*(-2*p14* (p44 + 2*p14 - 3*p23)* p23 + pow2(p13)* (p44 - 2*p14
    + 2*p24) + p13 * (-2*pow2(p14) + 8*pow2(p23) - 3*p23* (p44 - 2*p24)
    + p24* (-p44 + 4*p24) + p14* (p44 + 6*p24)))))
    / (2. * pow2(p23M) * p14M * p45P);

  double dia31 = (-2*(p13* (p44 - 2*p14) + p14* (p44 + 2*p12 - 2*p14
    - 2*p23 - 2*p24))* (p44 + 2*p12 - 2*p14 - 2*p24))
    / (pow2(p14M) * pow2(p35P));

  double dia32 = -((p44 - 2*p14)* (p12* (p44 - 2*p14) - p13* (p44 - 2*p14)
    + 2*p14*p23)) / (2. * pow2(p14M) * p25M * p35P);

  double dia33 = (-9*(4*pow3(p14) + 2*pow2(p12)* (p44 + 2*p14)
    - 4*pow2(p14)* (p44 + p23 - 3*p24) - 2*p44*p13* p24
    + p14* (p44S + 2*p44*p23 + (-6*p44 + 4*p13)* p24 + 8*pow2(p24))
    + p12 * (2*p13* (p44 - 2*p14) + p44*(3*p44 + 2*p23 - 2*p24)
    - 2*p14* (3*p44 + 2*p23 + 6*p24)))) / (4. * p23M * pow2(p14M) * p35P);

  double dia34 = (-4*(p12 - p23 - p24)* (p44*p12-2*p14*p24))
    / (pow2(p14M)* pow2(p25M));

  double dia35 = (-9*(2*pow2(p12)* (p44 + 2*p14) - 2*p13 * (p44 - 2*p14)
    * (p44 - 2*p14 + 2*p23 + p24) + 2*p14* (4*pow2(p23) + 8*p23* p24
    + p24* (-p44 + 2*p14 + 4*p24)) + p12* (2*p13* (p44 - 2*p14)
    + p44*(p44 - 2*p23 - 2*p24) - 2*p14* (p44 + 6*p23 + 6*p24))))
    / (4. * p23M * pow2(p14M) * p25M);

  double dia36 = (-9*(4*pow2(p12)* (p44 + 2*p14) - p13 * (p44 - 2*p14)
    * (p44 - 2*p14 + 8*p23 + 4*p24) + p14* (-p44S - 4*pow2(p14)
    + 16*pow2(p23) - 4*p23 * (p44 - 4*p24) - 8*p44*p24 + 16*pow2(p24)
    + 4*p14* (p44 + 2*p23 + 4*p24)) + p12* (4*p13 * (p44 - 2*p14)
    - 4*pow2(p14) + p44*(3*p44 + 4*p23 - 4*p24) - 4*p14* (p44 + 6*p23
    + 6*p24)))) / (4. * pow2(p23M) * pow2(p14M));

  return dia1 + dia2 + dia3 + dia4 + dia5 + dia6 + dia7 + dia8 + dia9
   + dia10 + dia11 + dia12 + dia13 + dia14 + dia15 + dia16 + dia17 + dia18
   + dia19 + dia20 + dia21 + dia22 + dia23 + dia24 + dia25 + dia26 + dia27
   + dia28 + dia29 + dia30 + dia31 + dia32 + dia33 + dia34 + dia35 + dia36;

}

//--------------------------------------------------------------------------

// Calculate the 2->3 ME ud -> udZ, with the coupling between Z and d
// set to zero, and omitting an overall factor of
// \frac{(g_s^4 * 4\pi \alpha_{em} * (lU^2 + rU^2))}
//      {9 * \cos^2(\theta_W) \sin^2(\theta_W)}..
// p1 = incoming u, p2 = incoming d, p3 = outgoing Z,
// p4 = outgoing d, p5 = outgoing u.

double SimpleWeakShowerMEs::getMEqq2qqZ(Vec4 p1, Vec4 p2, Vec4 p3, Vec4 p4,
  Vec4 p5) {

  double p12  = p1*p2;
  double p13  = p1*p3;
  double p14  = p1*p4;
  double p23  = p2*p3;
  double p24  = p2*p4;
  double p33  = p3*p3;
  double p13M = (p1 - p3).m2Calc();
  double p24M = (p2 - p4).m2Calc();
  double p35P = (p3 + p5).m2Calc();
  double p33S = p33*p33;

  double dia1 = (-4*(2*pow3(p12) + pow2(p12) * (p33 - 2*p13 - 4*p14
    - 2*p23 - 4*p24) + p14* (2*p14* p23 - (p33 - 2*p23)* p24)
    + p12 * (2*pow2(p14) + 2*p24* (p13 + p23 + p24)
    + p14* (p33 + 2*p13 + 4*p24)))) / pow2(p24M* p35P);

  double dia2 = (-2*(4*pow3(p12) + 4*pow2(p12)
    * (p33 - 2*p14 - 3*p23) + p12 * (p33S - 4*pow2(p13) + 4*pow2(p14)
    - 6*p33*p23 + 8*pow2(p23) + 4*p13 * (p23 - p24) - 4*p33*p24 + 4*p23
    * p24 + 4*pow2(p24) + 4*p14* (p33 + 4*p23 + 4*p24)) + 2*(-2*pow2(p14)
    * p23 + p13* (p33 + 2*p13 - 2*p24) * p24 + p14* (-4*pow2(p23) + p23
    * (p33 - 6*p13 - 6*p24) + 2*(p33 - p13 - 2*p24)* p24))))
    / (p13M * pow2(p24M) * p35P);

  double dia3 = (-2*(2*pow2(p12)* (p33 + 2*p13)
    - 2*p33*p14* (p23 + p24) + 4*pow2(p13)* (2*p23 + p24) + p12
    * (-4*pow2(p13) + p33*(p33 + 2*p14 - 2*p23) - 4*p13* (p14 + 3*p23
    + 2*p24)) + p13* (8*pow2(p23) + 2*p24 * (-p33 + 2*p14 + 2*p24) + p23
    * (-4*p33 + 4*p14 + 8*p24)))) / pow2(p13M * p24M);

  return dia1 + dia2 + dia3;

}

//==========================================================================

} // end namespace Pythia8
