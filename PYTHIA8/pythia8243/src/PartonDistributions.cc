// PartonDistributions.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the PDF, LHAPDF,
// GRV94L, CTEQ5L,  MSTWpdf, CTEQ6pdf, GRVpiL, PomFix, PomH1FitAB,
// PomH1Jets, Lepton, NNPDF and CJKL classes.

#include "Pythia8/PartonDistributions.h"

namespace Pythia8 {

//==========================================================================

// Base class for parton distribution functions.

//--------------------------------------------------------------------------

// Resolve valence content for assumed meson. Possibly modified later.

void PDF::setValenceContent() {

  // Subdivide meson by flavour content.
  if (idBeamAbs < 100 || idBeamAbs > 1000) return;
  int idTmp1 = idBeamAbs/100;
  int idTmp2 = (idBeamAbs/10)%10;

  // Find which is quark and which antiquark.
  if (idTmp1%2 == 0) {
    idVal1 =  idTmp1;
    idVal2 = -idTmp2;
  } else {
    idVal1 =  idTmp2;
    idVal2 = -idTmp1;
  }
  if (idBeam < 0) {
    idVal1 = -idVal1;
    idVal2 = -idVal2;
  }

  // Special case for Pomeron, to start off.
  if (idBeamAbs == 990) {
    idVal1 =  1;
    idVal2 = -1;
  }

  // Photon not fixed until at Process-/PartonLevel.
  if (idBeamAbs == 22) {
    idVal1 =  10;
    idVal2 = -10;
  }

}

//--------------------------------------------------------------------------

// Standard parton densities.

double PDF::xf(int id, double x, double Q2) {

  // Need to update if flavour, x or Q2 changed.
  // Use idSav = 9 to indicate that ALL flavours are up-to-date.
  // Assume that flavour and antiflavour always updated simultaneously.
  if ( (abs(idSav) != abs(id) && idSav != 9) || x != xSav || Q2 != Q2Sav)
    {idSav = id; xfUpdate(id, x, Q2); xSav = x; Q2Sav = Q2;}

  // Baryon beams: only p and pbar for now.
  if (idBeamAbs == 2212) {
    int idNow = (idBeam > 0) ? id : -id;
    int idAbs = abs(id);
    if (idNow ==  0 || idAbs == 21) return max(0., xg);
    if (idNow ==  1) return max(0., xd);
    if (idNow == -1) return max(0., xdbar);
    if (idNow ==  2) return max(0., xu);
    if (idNow == -2) return max(0., xubar);
    if (idNow ==  3) return max(0., xs);
    if (idNow == -3) return max(0., xsbar);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Baryon beams: n and nbar by isospin conjugation of p and pbar.
  } else if (idBeamAbs == 2112) {
    int idNow = (idBeam > 0) ? id : -id;
    int idAbs = abs(id);
    if (idNow ==  0 || idAbs == 21) return max(0., xg);
    if (idNow ==  1) return max(0., xu);
    if (idNow == -1) return max(0., xubar);
    if (idNow ==  2) return max(0., xd);
    if (idNow == -2) return max(0., xdbar);
    if (idNow ==  3) return max(0., xs);
    if (idNow == -3) return max(0., xsbar);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Nondiagonal meson beams: only pi+ and pi- for now.
  // Some LHAPDF sets are stored with u d as valence, so use dbar = u.
  } else if (idBeamAbs == 211) {
    int idNow = (idBeam > 0) ? id : -id;
    int idAbs = abs(id);
    if (idNow ==  0 || idAbs == 21) return max(0., xg);
    if (idNow ==  1) return max(0., xubar );
    if (idNow == -1) return max(0., xu );
    if (idNow ==  2) return max(0., xu);
    if (idNow == -2) return max(0., xubar);
    if (idNow ==  3) return max(0., xs);
    if (idNow == -3) return max(0., xsbar);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Diagonal meson beams: only pi0, Pomeron for now.
  } else if (idBeam == 111 || idBeam == 990) {
    int idAbs = abs(id);
    if (id ==  0 || idAbs == 21) return max(0., xg);
    if (id == idVal1 || id == idVal2) return max(0., xu);
    if (idAbs <=  2) return max(0., xubar);
    if (idAbs ==  3) return max(0., xs);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Photon beam.
  } else if (idBeam == 22) {
    int idAbs = abs(id);
    if (id ==  0 || idAbs == 21) return max(0., xg);
    if (id ==  1)    return max(0., xd);
    if (id == -1)    return max(0., xdbar);
    if (id ==  2)    return max(0., xu);
    if (id == -2)    return max(0., xubar);
    if (id ==  3)    return max(0., xs);
    if (id == -3)    return max(0., xsbar);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Photon beam inside lepton beam.
  } else if ( ( idBeamAbs == 11 || idBeamAbs == 13 || idBeamAbs == 15 )
    && hasGammaInLepton ) {
    int idAbs = abs(id);
    if (idAbs ==  0 || idAbs == 21) return max(0., xg);
    if (idAbs ==  1) return max(0., xd);
    if (idAbs ==  2) return max(0., xu);
    if (idAbs ==  3) return max(0., xs);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Nuclear PDFs
  } else if ( idBeamAbs > 100000000 ) {
    int idAbs = abs(id);
    if (idAbs ==  0 || idAbs == 21) return max(0., xg);
    if (id ==  1)    return max(0., xd);
    if (id == -1)    return max(0., xdbar);
    if (id ==  2)    return max(0., xu);
    if (id == -2)    return max(0., xubar);
    if (id ==  3)    return max(0., xs);
    if (id == -3)    return max(0., xsbar);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0;

  // Lepton beam.
  } else {
    if (id == idBeam ) return max(0., xlepton);
    if (abs(id) == 22) return max(0., xgamma);
    return 0.;
  }

}

//--------------------------------------------------------------------------

// Only valence part of parton densities.

double PDF::xfVal(int id, double x, double Q2) {

  // Need to update if flavour, x or Q2 changed.
  // Use idSav = 9 to indicate that ALL flavours are up-to-date.
  // Assume that flavour and antiflavour always updated simultaneously.
  if ( (abs(idSav) != abs(id) && idSav != 9) || x != xSav || Q2 != Q2Sav)
    {idSav = id; xfUpdate(id, x, Q2); xSav = x; Q2Sav = Q2;}

  // Baryon and nondiagonal meson beams: only p, pbar, n, nbar, pi+, pi-.
  if (idBeamAbs == 2212) {
    int idNow = (idBeam > 0) ? id : -id;
    if (idNow == 1) return max(0., xdVal);
    if (idNow == 2) return max(0., xuVal);
    return 0.;
  } else if (idBeamAbs == 2112) {
    int idNow = (idBeam > 0) ? id : -id;
    if (idNow == 1) return max(0., xuVal);
    if (idNow == 2) return max(0., xdVal);
    return 0.;
  } else if (idBeamAbs == 211) {
    int idNow = (idBeam > 0) ? id : -id;
    if (idNow == 2 || idNow == -1) return max(0., xuVal);
    return 0.;

  // Diagonal meson beams: only pi0, Pomeron for now.
  } else if (idBeam == 111 || idBeam == 990) {
    if (id == idVal1 || id == idVal2) return max(0., xuVal);
    return 0.;

  // Photon beam.
  } else if (idBeam == 22) {
    int idAbs = abs(id);
    if (id == idVal1 || id == idVal2) {
      if (idAbs == 1) return max(0., xdVal);
      if (idAbs == 2) return max(0., xuVal);
      if (idAbs == 3) return max(0., xsVal);
      if (idAbs == 4) return max(0., xcVal);
      if (idAbs == 5) return max(0., xbVal);
    }
    return 0.;

  // Lepton beam.
  } else {
    if (id == idBeam) return max(0., xlepton);
    return 0.;
  }

}

//--------------------------------------------------------------------------

// Only sea part of parton densities.

double PDF::xfSea(int id, double x, double Q2) {

  // Need to update if flavour, x or Q2 changed.
  // Use idSav = 9 to indicate that ALL flavours are up-to-date.
  // Assume that flavour and antiflavour always updated simultaneously.
  if ( (abs(idSav) != abs(id) && idSav != 9) || x != xSav || Q2 != Q2Sav)
    {idSav = id; xfUpdate(id, x, Q2); xSav = x; Q2Sav = Q2;}

  // Hadron beams.
  if (idBeamAbs > 100) {
    int idNow = (idBeam > 0) ? id : -id;
    int idAbs = abs(id);
    if (idNow == 0 || idAbs == 21) return max(0., xg);
    if (idBeamAbs == 2212) {
      if (idNow ==  1) return max(0., xdSea);
      if (idNow == -1) return max(0., xdbar);
      if (idNow ==  2) return max(0., xuSea);
      if (idNow == -2) return max(0., xubar);
    } else if (idBeamAbs == 2112) {
      if (idNow ==  1) return max(0., xuSea);
      if (idNow == -1) return max(0., xubar);
      if (idNow ==  2) return max(0., xdSea);
      if (idNow == -2) return max(0., xdbar);
    } else {
      if (idAbs <=  2) return max(0., xuSea);
    }
    if (idNow ==  3) return max(0., xs);
    if (idNow == -3) return max(0., xsbar);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Photon beam.
  } else if (idBeamAbs == 22) {
    int idAbs = abs(id);
    if ( id == 0 || idAbs == 21 ) return max(0., xg);
    if ( idAbs == 22 ) return max(0., xgamma);

    // If a valence parton return only the sea part.
    // Otherwise return the total PDF.
    if ( id == idVal1 || id == idVal2 ) {
      if (idAbs ==  1) return max(0., xdSea);
      if (idAbs ==  2) return max(0., xuSea);
      if (idAbs ==  3) return max(0., xsSea);
      if (idAbs ==  4) return max(0., xcSea);
      if (idAbs ==  5) return max(0., xbSea);
    } else {
      if (idAbs ==  1) return max(0., xd);
      if (idAbs ==  2) return max(0., xu);
      if (idAbs ==  3) return max(0., xs);
      if (idAbs ==  4) return max(0., xc);
      if (idAbs ==  5) return max(0., xb);
    }
    return 0.;

  // Lepton beam.
  } else {
    if (abs(id) == 22) return max(0., xgamma);
    return 0.;
  }

}

//==========================================================================

// Gives the GRV 94 L (leading order) parton distribution function set
// in parametrized form. Authors: M. Glueck, E. Reya and A. Vogt.
// Ref: M. Glueck, E. Reya and A. Vogt, Z.Phys. C67 (1995) 433.

void GRV94L::xfUpdate(int , double x, double Q2) {

  // Common expressions. Constrain Q2 for which parametrization is valid.
  double mu2  = 0.23;
  double lam2 = 0.2322 * 0.2322;
  double s    = (Q2 > mu2) ? log( log(Q2/lam2) / log(mu2/lam2) ) : 0.;
  double ds   = sqrt(s);
  double s2   = s * s;
  double s3   = s2 * s;

  // uv :
  double nu  =  2.284 + 0.802 * s + 0.055 * s2;
  double aku =  0.590 - 0.024 * s;
  double bku =  0.131 + 0.063 * s;
  double au  = -0.449 - 0.138 * s - 0.076 * s2;
  double bu  =  0.213 + 2.669 * s - 0.728 * s2;
  double cu  =  8.854 - 9.135 * s + 1.979 * s2;
  double du  =  2.997 + 0.753 * s - 0.076 * s2;
  double uv  = grvv (x, nu, aku, bku, au, bu, cu, du);

  // dv :
  double nd  =  0.371 + 0.083 * s + 0.039 * s2;
  double akd =  0.376;
  double bkd =  0.486 + 0.062 * s;
  double ad  = -0.509 + 3.310 * s - 1.248 * s2;
  double bd  =  12.41 - 10.52 * s + 2.267 * s2;
  double cd  =  6.373 - 6.208 * s + 1.418 * s2;
  double dd  =  3.691 + 0.799 * s - 0.071 * s2;
  double dv  = grvv (x, nd, akd, bkd, ad, bd, cd, dd);

  // udb :
  double alx =  1.451;
  double bex =  0.271;
  double akx =  0.410 - 0.232 * s;
  double bkx =  0.534 - 0.457 * s;
  double agx =  0.890 - 0.140 * s;
  double bgx = -0.981;
  double cx  =  0.320 + 0.683 * s;
  double dx  =  4.752 + 1.164 * s + 0.286 * s2;
  double ex  =  4.119 + 1.713 * s;
  double esx =  0.682 + 2.978 * s;
  double udb = grvw (x, s, alx, bex, akx, bkx, agx, bgx, cx,
    dx, ex, esx);

  // del :
  double ne  =  0.082 + 0.014 * s + 0.008 * s2;
  double ake =  0.409 - 0.005 * s;
  double bke =  0.799 + 0.071 * s;
  double ae  = -38.07 + 36.13 * s - 0.656 * s2;
  double be  =  90.31 - 74.15 * s + 7.645 * s2;
  double ce  =  0.;
  double de  =  7.486 + 1.217 * s - 0.159 * s2;
  double del = grvv (x, ne, ake, bke, ae, be, ce, de);

  // sb :
  double sts =  0.;
  double als =  0.914;
  double bes =  0.577;
  double aks =  1.798 - 0.596 * s;
  double as  = -5.548 + 3.669 * ds - 0.616 * s;
  double bs  =  18.92 - 16.73 * ds + 5.168 * s;
  double dst =  6.379 - 0.350 * s  + 0.142 * s2;
  double est =  3.981 + 1.638 * s;
  double ess =  6.402;
  double sb  = grvs (x, s, sts, als, bes, aks, as, bs, dst, est, ess);

  // cb :
  double stc =  0.888;
  double alc =  1.01;
  double bec =  0.37;
  double akc =  0.;
  double ac  =  0.;
  double bc  =  4.24  - 0.804 * s;
  double dct =  3.46  - 1.076 * s;
  double ect =  4.61  + 1.49  * s;
  double esc =  2.555 + 1.961 * s;
  double chm = grvs (x, s, stc, alc, bec, akc, ac, bc, dct, ect, esc);

  // bb :
  double stb =  1.351;
  double alb =  1.00;
  double beb =  0.51;
  double akb =  0.;
  double ab  =  0.;
  double bb  =  1.848;
  double dbt =  2.929 + 1.396 * s;
  double ebt =  4.71  + 1.514 * s;
  double esb =  4.02  + 1.239 * s;
  double bot = grvs (x, s, stb, alb, beb, akb, ab, bb, dbt, ebt, esb);

  // gl :
  double alg =  0.524;
  double beg =  1.088;
  double akg =  1.742 - 0.930 * s;
  double bkg =                     - 0.399 * s2;
  double ag  =  7.486 - 2.185 * s;
  double bg  =  16.69 - 22.74 * s  + 5.779 * s2;
  double cg  = -25.59 + 29.71 * s  - 7.296 * s2;
  double dg  =  2.792 + 2.215 * s  + 0.422 * s2 - 0.104 * s3;
  double eg  =  0.807 + 2.005 * s;
  double esg =  3.841 + 0.316 * s;
  double gl  = grvw (x, s, alg, beg, akg, bkg, ag, bg, cg,
    dg, eg, esg);

  // Update values
  xg    = gl;
  xu    = uv + 0.5*(udb - del);
  xd    = dv + 0.5*(udb + del);
  xubar = 0.5*(udb - del);
  xdbar = 0.5*(udb + del);
  xs    = sb;
  xsbar = sb;
  xc    = chm;
  xb    = bot;

  // Subdivision of valence and sea.
  xuVal = uv;
  xuSea = xubar;
  xdVal = dv;
  xdSea = xdbar;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//--------------------------------------------------------------------------

double GRV94L::grvv (double x, double n, double ak, double bk, double a,
   double b, double c, double d) {

  double dx = sqrt(x);
  return n * pow(x, ak) * (1. + a * pow(x, bk) + x * (b + c * dx)) *
    pow(1. - x, d);

}

//--------------------------------------------------------------------------

double GRV94L::grvw (double x, double s, double al, double be, double ak,
  double bk, double a, double b, double c, double d, double e, double es) {

  double lx = log(1./x);
  return (pow(x, ak) * (a + x * (b + x * c)) * pow(lx, bk) + pow(s, al)
    * exp(-e + sqrt(es * pow(s, be) * lx))) * pow(1. - x, d);

}

//--------------------------------------------------------------------------

double GRV94L::grvs (double x, double s, double sth, double al, double be,
  double ak, double ag, double b, double d, double e, double es) {

  if(s <= sth) {
    return 0.;
  } else {
    double dx = sqrt(x);
    double lx = log(1./x);
    return pow(s - sth, al) / pow(lx, ak) * (1. + ag * dx + b * x) *
      pow(1. - x, d) * exp(-e + sqrt(es * pow(s, be) * lx));
  }

}

//==========================================================================

// Gives the CTEQ 5 L (leading order) parton distribution function set
// in parametrized form. Parametrization by J. Pumplin.
// Ref: CTEQ Collaboration, H.L. Lai et al., Eur.Phys.J. C12 (2000) 375.

// The range of (x, Q) covered by this parametrization of the QCD
// evolved parton distributions is 1E-6 < x < 1, 1.1 GeV < Q < 10 TeV.
// In the current implementation, densities are frozen at borders.

void CTEQ5L::xfUpdate(int , double x, double Q2) {

  // Constrain x and Q2 to range for which parametrization is valid.
  double Q = sqrt( max( 1., min( 1e8, Q2) ) );
  x = max( 1e-6, min( 1.-1e-10, x) );

  // Derived kinematical quantities.
  double y = - log(x);
  double u = log( x / 0.00001);
  double x1 = 1. - x;
  double x1L = log(1. - x);
  double sumUbarDbar = 0.;

  // Parameters of parametrizations.
  const double Qmin[8] = { 0., 0., 0., 0., 0., 0., 1.3, 4.5};
  const double alpha[8] = { 0.2987216, 0.3407552, 0.4491863, 0.2457668,
    0.5293999, 0.3713141, 0.03712017, 0.004952010 };
  const double ut1[8] = { 4.971265, 2.612618, -0.4656819, 3.862583,
    0.1895615, 3.753257, 4.400772, 5.562568 };
  const double ut2[8] = { -1.105128, -1.258304e5, -274.2390, -1.265969,
    -3.069097, -1.113085, -1.356116, -1.801317 };
  const double am[8][9][3] = {
    // d.
    { {  0.5292616E+01, -0.2751910E+01, -0.2488990E+01 },
      {  0.9714424E+00,  0.1011827E-01, -0.1023660E-01 },
      { -0.1651006E+02,  0.7959721E+01,  0.8810563E+01 },
      { -0.1643394E+02,  0.5892854E+01,  0.9348874E+01 },
      {  0.3067422E+02,  0.4235796E+01, -0.5112136E+00 },
      {  0.2352526E+02, -0.5305168E+01, -0.1169174E+02 },
      { -0.1095451E+02,  0.3006577E+01,  0.5638136E+01 },
      { -0.1172251E+02, -0.2183624E+01,  0.4955794E+01 },
      {  0.1662533E-01,  0.7622870E-02, -0.4895887E-03 } },
    // u.
    { {  0.9905300E+00, -0.4502235E+00,  0.1624441E+00 },
      {  0.8867534E+00,  0.1630829E-01, -0.4049085E-01 },
      {  0.8547974E+00,  0.3336301E+00,  0.1371388E+00 },
      {  0.2941113E+00, -0.1527905E+01,  0.2331879E+00 },
      {  0.3384235E+02,  0.3715315E+01,  0.8276930E+00 },
      {  0.6230115E+01,  0.3134639E+01, -0.1729099E+01 },
      { -0.1186928E+01, -0.3282460E+00,  0.1052020E+00 },
      { -0.8545702E+01, -0.6247947E+01,  0.3692561E+01 },
      {  0.1724598E-01,  0.7120465E-02,  0.4003646E-04 } },
    // g.
    { {  0.1193572E+03, -0.3886845E+01, -0.1133965E+01 },
      { -0.9421449E+02,  0.3995885E+01,  0.1607363E+01 },
      {  0.4206383E+01,  0.2485954E+00,  0.2497468E+00 },
      {  0.1210557E+03, -0.3015765E+01, -0.1423651E+01 },
      { -0.1013897E+03, -0.7113478E+00,  0.2621865E+00 },
      { -0.1312404E+01, -0.9297691E+00, -0.1562531E+00 },
      {  0.1627137E+01,  0.4954111E+00, -0.6387009E+00 },
      {  0.1537698E+00, -0.2487878E+00,  0.8305947E+00 },
      {  0.2496448E-01,  0.2457823E-02,  0.8234276E-03 } },
    // ubar + dbar.
    { {  0.2647441E+02,  0.1059277E+02, -0.9176654E+00 },
      {  0.1990636E+01,  0.8558918E-01,  0.4248667E-01 },
      { -0.1476095E+02, -0.3276255E+02,  0.1558110E+01 },
      { -0.2966889E+01, -0.3649037E+02,  0.1195914E+01 },
      { -0.1000519E+03, -0.2464635E+01,  0.1964849E+00 },
      {  0.3718331E+02,  0.4700389E+02, -0.2772142E+01 },
      { -0.1872722E+02, -0.2291189E+02,  0.1089052E+01 },
      { -0.1628146E+02, -0.1823993E+02,  0.2537369E+01 },
      { -0.1156300E+01, -0.1280495E+00,  0.5153245E-01 } },
    // dbar/ubar.
    { { -0.6556775E+00,  0.2490190E+00,  0.3966485E-01 },
      {  0.1305102E+01, -0.1188925E+00, -0.4600870E-02 },
      { -0.2371436E+01,  0.3566814E+00, -0.2834683E+00 },
      { -0.6152826E+01,  0.8339877E+00, -0.7233230E+00 },
      { -0.8346558E+01,  0.2892168E+01,  0.2137099E+00 },
      {  0.1279530E+02,  0.1021114E+00,  0.5787439E+00 },
      {  0.5858816E+00, -0.1940375E+01, -0.4029269E+00 },
      { -0.2795725E+02, -0.5263392E+00,  0.1290229E+01 },
      {  0.0000000E+00,  0.0000000E+00,  0.0000000E+00 } },
    // sbar.
    { {  0.1580931E+01, -0.2273826E+01, -0.1822245E+01 },
      {  0.2702644E+01,  0.6763243E+00,  0.7231586E-02 },
      { -0.1857924E+02,  0.3907500E+01,  0.5850109E+01 },
      { -0.3044793E+02,  0.2639332E+01,  0.5566644E+01 },
      { -0.4258011E+01, -0.5429244E+01,  0.4418946E+00 },
      {  0.3465259E+02, -0.5532604E+01, -0.4904153E+01 },
      { -0.1658858E+02,  0.2923275E+01,  0.2266286E+01 },
      { -0.1149263E+02,  0.2877475E+01, -0.7999105E+00 },
      {  0.0000000E+00,  0.0000000E+00,  0.0000000E+00 } },
    // cbar.
    { { -0.8293661E+00, -0.3982375E+01, -0.6494283E-01 },
      {  0.2754618E+01,  0.8338636E+00, -0.6885160E-01 },
      { -0.1657987E+02,  0.1439143E+02, -0.6887240E+00 },
      { -0.2800703E+02,  0.1535966E+02, -0.7377693E+00 },
      { -0.6460216E+01, -0.4783019E+01,  0.4913297E+00 },
      {  0.3141830E+02, -0.3178031E+02,  0.7136013E+01 },
      { -0.1802509E+02,  0.1862163E+02, -0.4632843E+01 },
      { -0.1240412E+02,  0.2565386E+02, -0.1066570E+02 },
      {  0.0000000E+00,  0.0000000E+00,  0.0000000E+00 } },
    // bbar.
    { { -0.6031237E+01,  0.1992727E+01, -0.1076331E+01 },
      {  0.2933912E+01,  0.5839674E+00,  0.7509435E-01 },
      { -0.8284919E+01,  0.1488593E+01, -0.8251678E+00 },
      { -0.1925986E+02,  0.2805753E+01, -0.3015446E+01 },
      { -0.9480483E+01, -0.9767837E+00, -0.1165544E+01 },
      {  0.2193195E+02, -0.1788518E+02,  0.9460908E+01 },
      { -0.1327377E+02,  0.1201754E+02, -0.6277844E+01 },
      {  0.0000000E+00,  0.0000000E+00,  0.0000000E+00 },
      {  0.0000000E+00,  0.0000000E+00,  0.0000000E+00 } } };

  // Loop over 8 different parametrizations. Check if inside allowed region.
  for (int i = 0; i < 8; ++i) {
    double answer = 0.;
    if (Q > max(Qmin[i], alpha[i])) {

      // Evaluate answer.
      double tmp = log(Q / alpha[i]);
      double sb = log(tmp);
      double sb1 = sb - 1.2;
      double sb2 = sb1*sb1;
      double af[9];
      for (int j = 0; j < 9; ++j)
        af[j] = am[i][j][0] + sb1 * am[i][j][1] + sb2 * am[i][j][2];
      double part1 = af[1] * pow( y, 1. + 0.01 * af[4]) * (1. + af[8] * u);
      double part2 = af[0] * x1 + af[3] * x;
      double part3 = x * x1 * (af[5] + af[6] * x1 + af[7] * x * x1);
      double part4 = (ut2[i] < -100.) ? ut1[i] * x1L + af[2] * x1L
                   : ut1[i] * x1L + af[2] * log(x1 + exp(ut2[i]));
      answer       = x * exp( part1 + part2 + part3 + part4);
      answer      *= 1. - Qmin[i] / Q;
    }

    // Store results.
    if (i == 0) xd = x * answer;
    else if (i == 1) xu = x * answer;
    else if (i == 2) xg = x * answer;
    else if (i == 3) sumUbarDbar = x * answer;
    else if (i == 4) { xubar = sumUbarDbar / (1. + answer);
      xdbar = sumUbarDbar * answer / (1. + answer); }
    else if (i == 5) {xs = x * answer; xsbar = xs;}
    else if (i == 6) xc = x * answer;
    else if (i == 7) xb = x * answer;
  }

  // Subdivision of valence and sea.
  xuVal = xu - xubar;
  xuSea = xubar;
  xdVal = xd - xdbar;
  xdSea = xdbar;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//==========================================================================

// The MSTWpdf class.
// MSTW 2008 PDF's, specifically the LO one.
// Original C++ version by Jeppe Andersen.
// Modified by Graeme Watt <watt(at)hep.ucl.ac.uk>.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of parton flavours, x and Q2 grid points,
// bins below c and b thresholds.
const int    MSTWpdf::np     = 12;
const int    MSTWpdf::nx     = 64;
const int    MSTWpdf::nq     = 48;
const int    MSTWpdf::nqc0   = 4;
const int    MSTWpdf::nqb0   = 14;

// Range of (x, Q2) grid.
const double MSTWpdf::xmin   = 1e-6;
const double MSTWpdf::xmax   = 1.0;
const double MSTWpdf::qsqmin = 1.0;
const double MSTWpdf::qsqmax = 1e9;

// Array of x values.
const double MSTWpdf::xxInit[65] = {0., 1e-6, 2e-6, 4e-6, 6e-6, 8e-6,
  1e-5, 2e-5, 4e-5, 6e-5, 8e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4,
  1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2, 1.4e-2, 2e-2, 3e-2, 4e-2, 6e-2,
  8e-2, 0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.30,
  0.325, 0.35, 0.375, 0.40, 0.425, 0.45, 0.475, 0.50, 0.525, 0.55,
  0.575, 0.60, 0.625, 0.65, 0.675, 0.70, 0.725, 0.75, 0.775, 0.80,
  0.825, 0.85, 0.875, 0.90, 0.925, 0.95, 0.975, 1.0 };

// Array of Q values.
const double MSTWpdf::qqInit[49] = {0., 1.0, 1.25, 1.5, 0., 0., 2.5, 3.2,
  4.0, 5.0, 6.4, 8.0, 10., 12., 0., 0., 26.0, 40.0, 64.0, 1e2, 1.6e2,
  2.4e2, 4e2, 6.4e2, 1e3, 1.8e3, 3.2e3, 5.6e3, 1e4, 1.8e4, 3.2e4, 5.6e4,
  1e5, 1.8e5, 3.2e5, 5.6e5, 1e6, 1.8e6, 3.2e6, 5.6e6, 1e7, 1.8e7, 3.2e7,
  5.6e7, 1e8, 1.8e8, 3.2e8, 5.6e8, 1e9 };

//--------------------------------------------------------------------------

// Initialize PDF: select data file and open stream.

void MSTWpdf::init(int iFitIn, string xmlPath, Info* infoPtr) {

  // Choice of fit among possibilities.
  iFit = iFitIn;

  // Select which data file to read for current fit.
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";
  string fileName = "  ";
  if (iFit == 1) fileName = "mrstlostar.00.dat";
  if (iFit == 2) fileName = "mrstlostarstar.00.dat";
  if (iFit == 3) fileName = "mstw2008lo.00.dat";
  if (iFit == 4) fileName = "mstw2008nlo.00.dat";

  // Open data file.
  ifstream data_file( (xmlPath + fileName).c_str() );
  if (!data_file.good()) {
    printErr("Error in MSTWpdf::init: did not find data file ", infoPtr);
    isSet = false;
    return;
  }

  // Initialization with a stream.
  init(data_file, infoPtr);
  data_file.close();

}

//--------------------------------------------------------------------------

// Initialize PDF: read in data grid from stream and set up interpolation.

void MSTWpdf::init(istream& data_file, Info* infoPtr) {

  // Check that data stream is available.
  if (!data_file.good()) {
    printErr("Error in MSTWpdf::init: cannot read from stream", infoPtr);
    isSet = false;
    return;
  }

  // Counters and temporary variables.
  int i,n,m,k,l,j;
  double dtemp;

  // Variables used for initialising c_ij array:
  double f[np+1][nx+1][nq+1];
  double f1[np+1][nx+1][nq+1]; // derivative w.r.t. x
  double f2[np+1][nx+1][nq+1]; // derivative w.r.t. q
  double f12[np+1][nx+1][nq+1];// cross derivative
  double f21[np+1][nx+1][nq+1];// cross derivative
  int wt[16][16]={{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                  {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
                  {-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
                  {2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
                  {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
                  {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
                  {0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
                  {0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
                  {-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
                  {0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
                  {9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
                  {-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
                  {2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
                  {0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
                  {-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
                  {4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}};
  double xxd,d1d2,cl[16],x[16],d1,d2,y[5],y1[5],y2[5],y12[5];
  double mc2,mb2,eps=1e-6; // q^2 grid points at mc2+eps, mb2+eps

  // Read distance, tolerance, heavy quark masses
  // and alphaS values from file.
  char comma;
  int nExtraFlavours;
  data_file.ignore(256,'\n');
  data_file.ignore(256,'\n');
  data_file.ignore(256,'='); data_file >> distance >> tolerance;
  data_file.ignore(256,'='); data_file >> mCharm;
  data_file.ignore(256,'='); data_file >> mBottom;
  data_file.ignore(256,'='); data_file >> alphaSQ0;
  data_file.ignore(256,'='); data_file >> alphaSMZ;
  data_file.ignore(256,'='); data_file >> alphaSorder >> comma >> alphaSnfmax;
  data_file.ignore(256,'='); data_file >> nExtraFlavours;
  data_file.ignore(256,'\n');
  data_file.ignore(256,'\n');
  data_file.ignore(256,'\n');

  // Use c and b quark masses for outlay of qq array.
  for (int iqq = 0; iqq < 49; ++iqq) qq[iqq] = qqInit[iqq];
  mc2=mCharm*mCharm;
  mb2=mBottom*mBottom;
  qq[4]=mc2;
  qq[5]=mc2+eps;
  qq[14]=mb2;
  qq[15]=mb2+eps;

  // Check that the heavy quark masses are sensible.
  if (mc2 < qq[3] || mc2 > qq[6]) {
    printErr("Error in MSTWpdf::init: invalid mCharm", infoPtr);
    isSet = false;
    return;
  }
  if (mb2 < qq[13] || mb2 > qq[16]) {
    printErr("Error in MSTWpdf::init: invalid mBottom", infoPtr);
    isSet = false;
    return;
  }

  // The nExtraFlavours variable is provided to aid compatibility
  // with future grids where, for example, a photon distribution
  // might be provided (cf. the MRST2004QED PDFs).
  if (nExtraFlavours < 0 || nExtraFlavours > 1) {
    printErr("Error in MSTWpdf::init: invalid nExtraFlavours", infoPtr);
    isSet = false;
    return;
  }

  // Now read in the grids from the grid file.
  for (n=1;n<=nx-1;n++)
    for (m=1;m<=nq;m++) {
      for (i=1;i<=9;i++)
        data_file >> f[i][n][m];
      if (alphaSorder==2) { // only at NNLO
        data_file >> f[10][n][m]; // = chm-cbar
        data_file >> f[11][n][m]; // = bot-bbar
      }
      else {
        f[10][n][m] = 0.; // = chm-cbar
        f[11][n][m] = 0.; // = bot-bbar
      }
      if (nExtraFlavours>0)
        data_file >> f[12][n][m];   // = photon
      else
        f[12][n][m] = 0.; // photon
      if (data_file.eof()) {
        printErr("Error in MSTWpdf::init: could not read data stream",
          infoPtr);
        isSet = false;
        return;
      }
    }

  // Check that ALL the file contents have been read in.
  data_file >> dtemp;
  if (!data_file.eof()) {
    printErr("Error in MSTWpdf::init: could not read data stream", infoPtr);
    isSet = false;
    return;
  }

  // PDFs are identically zero at x = 1.
  for (i=1;i<=np;i++)
    for (m=1;m<=nq;m++)
      f[i][nx][m]=0.0;

  // Set up the new array in log10(x) and log10(qsq).
  for (i=1;i<=nx;i++)
    xx[i]=log10(xxInit[i]);
  for (m=1;m<=nq;m++)
    qq[m]=log10(qq[m]);

  // Now calculate the derivatives used for bicubic interpolation.
  for (i=1;i<=np;i++) {

    // Start by calculating the first x derivatives
    // along the first x value:
    for (m=1;m<=nq;m++) {
      f1[i][1][m]=polderivative1(xx[1],xx[2],xx[3],f[i][1][m],f[i][2][m],
        f[i][3][m]);
      // Then along the rest (up to the last):
      for (k=2;k<nx;k++)
        f1[i][k][m]=polderivative2(xx[k-1],xx[k],xx[k+1],f[i][k-1][m],
          f[i][k][m],f[i][k+1][m]);
      // Then for the last column:
      f1[i][nx][m]=polderivative3(xx[nx-2],xx[nx-1],xx[nx],f[i][nx-2][m],
        f[i][nx-1][m],f[i][nx][m]);
    }

    // Then calculate the qq derivatives.  At NNLO there are
    // discontinuities in the PDFs at mc2 and mb2, so calculate
    // the derivatives at q^2 = mc2, mc2+eps, mb2, mb2+eps in
    // the same way as at the endpoints qsqmin and qsqmax.
    for (m=1;m<=nq;m++) {
      if (m==1 || m==nqc0+1 || m==nqb0+1) {
        for (k=1;k<=nx;k++)
          f2[i][k][m]=polderivative1(qq[m],qq[m+1],qq[m+2],
                                     f[i][k][m],f[i][k][m+1],f[i][k][m+2]);
      }
      else if (m==nq || m==nqc0 || m==nqb0) {
        for (k=1;k<=nx;k++)
          f2[i][k][m]=polderivative3(qq[m-2],qq[m-1],qq[m],
                                     f[i][k][m-2],f[i][k][m-1],f[i][k][m]);
      }
      else {
        // The rest:
        for (k=1;k<=nx;k++)
          f2[i][k][m]=polderivative2(qq[m-1],qq[m],qq[m+1],
                                    f[i][k][m-1],f[i][k][m],f[i][k][m+1]);
      }
    }

    // Now, calculate the cross derivatives.
    // Calculate these as the average between (d/dx)(d/dy) and (d/dy)(d/dx).

    // First calculate (d/dx)(d/dy).
    // Start by calculating the first x derivatives
    // along the first x value:
    for (m=1;m<=nq;m++)
      f12[i][1][m]=polderivative1(xx[1],xx[2],xx[3],f2[i][1][m],
        f2[i][2][m],f2[i][3][m]);
    // Then along the rest (up to the last):
    for (k=2;k<nx;k++) {
      for (m=1;m<=nq;m++)
        f12[i][k][m]=polderivative2(xx[k-1],xx[k],xx[k+1],f2[i][k-1][m],
          f2[i][k][m],f2[i][k+1][m]);
    }
    // Then for the last column:
    for (m=1;m<=nq;m++)
      f12[i][nx][m]=polderivative3(xx[nx-2],xx[nx-1],xx[nx],
        f2[i][nx-2][m],f2[i][nx-1][m],f2[i][nx][m]);

    // Now calculate (d/dy)(d/dx).
    for (m=1;m<=nq;m++) {
      if (m==1 || m==nqc0+1 || m==nqb0+1) {
        for (k=1;k<=nx;k++)
          f21[i][k][m]=polderivative1(qq[m],qq[m+1],qq[m+2],
                                      f1[i][k][m],f1[i][k][m+1],f1[i][k][m+2]);
      }
      else if (m==nq || m==nqc0 || m==nqb0) {
        for (k=1;k<=nx;k++)
          f21[i][k][m]=polderivative3(qq[m-2],qq[m-1],qq[m],
                                      f1[i][k][m-2],f1[i][k][m-1],f1[i][k][m]);
      }
      else {
        // The rest:
        for (k=1;k<=nx;k++)
          f21[i][k][m]=polderivative2(qq[m-1],qq[m],qq[m+1],
                                     f1[i][k][m-1],f1[i][k][m],f1[i][k][m+1]);
      }
    }

    // Now take the average of (d/dx)(d/dy) and (d/dy)(d/dx).
    for (k=1;k<=nx;k++) {
      for (m=1;m<=nq;m++) {
        f12[i][k][m] = 0.5*(f12[i][k][m]+f21[i][k][m]);
      }
    }

    // Now calculate the coefficients c_ij.
    for (n=1;n<=nx-1;n++) {
      for (m=1;m<=nq-1;m++) {
        d1=xx[n+1]-xx[n];
        d2=qq[m+1]-qq[m];
        d1d2=d1*d2;

        y[1]=f[i][n][m];
        y[2]=f[i][n+1][m];
        y[3]=f[i][n+1][m+1];
        y[4]=f[i][n][m+1];

        y1[1]=f1[i][n][m];
        y1[2]=f1[i][n+1][m];
        y1[3]=f1[i][n+1][m+1];
        y1[4]=f1[i][n][m+1];

        y2[1]=f2[i][n][m];
        y2[2]=f2[i][n+1][m];
        y2[3]=f2[i][n+1][m+1];
        y2[4]=f2[i][n][m+1];

        y12[1]=f12[i][n][m];
        y12[2]=f12[i][n+1][m];
        y12[3]=f12[i][n+1][m+1];
        y12[4]=f12[i][n][m+1];

        for (k=1;k<=4;k++) {
          x[k-1]=y[k];
          x[k+3]=y1[k]*d1;
          x[k+7]=y2[k]*d2;
          x[k+11]=y12[k]*d1d2;
        }

        for (l=0;l<=15;l++) {
          xxd=0.0;
          for (k=0;k<=15;k++) xxd+= wt[l][k]*x[k];
          cl[l]=xxd;
        }

        l=0;
        for (k=1;k<=4;k++)
          for (j=1;j<=4;j++) c[i][n][m][k][j]=cl[l++];
      } //m
    } //n
  } // i

}

//--------------------------------------------------------------------------

// Update PDF values.

void MSTWpdf::xfUpdate(int , double x, double Q2) {

  // Update using MSTW routine.
  double q    = sqrtpos(Q2);
  // Quarks:
  double dn   = parton(1,x,q);
  double up   = parton(2,x,q);
  double str  = parton(3,x,q);
  double chm  = parton(4,x,q);
  double bot  = parton(5,x,q);
  // Valence quarks:
  double dnv  = parton(7,x,q);
  double upv  = parton(8,x,q);
  double sv   = parton(9,x,q);
  double cv   = parton(10,x,q);
  double bv   = parton(11,x,q);
  // Antiquarks = quarks - valence quarks:
  double dsea = dn - dnv;
  double usea = up - upv;
  double sbar = str - sv;
  double cbar = chm - cv;
  double bbar = bot - bv;
  // Gluon:
  double glu  = parton(0,x,q);
  // Photon (= zero unless considering QED contributions):
  double phot = parton(13,x,q);

  // Transfer to Pythia notation.
  xg     = glu;
  xu     = up;
  xd     = dn;
  xubar  = usea;
  xdbar  = dsea;
  xs     = str;
  xsbar  = sbar;
  xc     = 0.5 * (chm + cbar);
  xb     = 0.5 * (bot + bbar);
  xgamma = phot;

  // Subdivision of valence and sea.
  xuVal  = upv;
  xuSea  = xubar;
  xdVal  = dnv;
  xdSea  = xdbar;

  // idSav = 9 to indicate that all flavours reset.
  idSav  = 9;

}

//--------------------------------------------------------------------------

// Returns the PDF value for parton of flavour 'f' at x,q.

double MSTWpdf::parton(int f,double x,double q) {

  double qsq;
  int ip;
  int interpolate(1);
  double parton_pdf=0,parton_pdf1=0,anom;
  double xxx,qqq;

  qsq=q*q;

  // If mc2 < qsq < mc2+eps, then qsq = mc2+eps.
  if (qsq>pow(10.,qq[nqc0]) && qsq<pow(10.,qq[nqc0+1])) {
    qsq = pow(10.,qq[nqc0+1]);
  }

  // If mb2 < qsq < mb2+eps, then qsq = mb2+eps.
  if (qsq>pow(10.,qq[nqb0]) && qsq<pow(10.,qq[nqb0+1])) {
    qsq = pow(10.,qq[nqb0+1]);
  }

  if (x<xmin) {
    interpolate=0;
    if (x<=0.) return 0.;
  }
  else if (x>xmax) return 0.;

  if (qsq<qsqmin) {
    interpolate=-1;
    if (q<=0.) return 0.;
  }
  else if (qsq>qsqmax) {
    interpolate=0;
  }

  if (f==0) ip=1;
  else if (f>=1 && f<=5) ip=f+1;
  else if (f<=-1 && f>=-5) ip=-f+1;
  else if (f>=7 && f<=11) ip=f;
  else if (f==13) ip=12;
  else if (abs(f)==6 || f==12) return 0.;
  else return 0.;

  // Interpolation in log10(x), log10(qsq):
  xxx=log10(x);
  qqq=log10(qsq);

  if (interpolate==1) { // do usual interpolation
    parton_pdf=parton_interpolate(ip,xxx,qqq);
    if (f<=-1 && f>=-5) // antiquark = quark - valence
      parton_pdf -= parton_interpolate(ip+5,xxx,qqq);
  }
  else if (interpolate==-1) { // extrapolate to low Q^2

    if (x<xmin) { // extrapolate to low x
      parton_pdf = parton_extrapolate(ip,xxx,log10(qsqmin));
      parton_pdf1 = parton_extrapolate(ip,xxx,log10(1.01*qsqmin));
      if (f<=-1 && f>=-5) { // antiquark = quark - valence
        parton_pdf -= parton_extrapolate(ip+5,xxx,log10(qsqmin));
        parton_pdf1 -= parton_extrapolate(ip+5,xxx,log10(1.01*qsqmin));
      }
    }
    else { // do usual interpolation
      parton_pdf = parton_interpolate(ip,xxx,log10(qsqmin));
      parton_pdf1 = parton_interpolate(ip,xxx,log10(1.01*qsqmin));
      if (f<=-1 && f>=-5) { // antiquark = quark - valence
        parton_pdf -= parton_interpolate(ip+5,xxx,log10(qsqmin));
        parton_pdf1 -= parton_interpolate(ip+5,xxx,log10(1.01*qsqmin));
      }
    }
    // Calculate the anomalous dimension, dlog(xf)/dlog(qsq),
    // evaluated at qsqmin.  Then extrapolate the PDFs to low
    // qsq < qsqmin by interpolating the anomalous dimenion between
    // the value at qsqmin and a value of 1 for qsq << qsqmin.
    // If value of PDF at qsqmin is very small, just set
    // anomalous dimension to 1 to prevent rounding errors.
    if (fabs(parton_pdf) >= 1.e-5)
      anom = max(-2.5, (parton_pdf1-parton_pdf)/parton_pdf/0.01);
    else anom = 1.;
    parton_pdf = parton_pdf*pow(qsq/qsqmin,anom*qsq/qsqmin+1.-qsq/qsqmin);

  }
  else { // extrapolate outside PDF grid to low x or high Q^2
    parton_pdf = parton_extrapolate(ip,xxx,qqq);
    if (f<=-1 && f>=-5) // antiquark = quark - valence
      parton_pdf -= parton_extrapolate(ip+5,xxx,qqq);
  }

  return parton_pdf;
}

//--------------------------------------------------------------------------

// Interpolate PDF value inside data grid.

double MSTWpdf::parton_interpolate(int ip, double xxx, double qqq) {

  double g, t, u;
  int    n, m, l;

  n=locate(xx,nx,xxx); // 0: below xmin, nx: above xmax
  m=locate(qq,nq,qqq); // 0: below qsqmin, nq: above qsqmax

  t=(xxx-xx[n])/(xx[n+1]-xx[n]);
  u=(qqq-qq[m])/(qq[m+1]-qq[m]);

  // Assume PDF proportional to (1-x)^p as x -> 1.
  if (n==nx-1) {
    double g0=((c[ip][n][m][1][4]*u+c[ip][n][m][1][3])*u
    +c[ip][n][m][1][2])*u+c[ip][n][m][1][1]; // value at xx[n]
    double g1=((c[ip][n-1][m][1][4]*u+c[ip][n-1][m][1][3])*u
           +c[ip][n-1][m][1][2])*u+c[ip][n-1][m][1][1]; // value at xx[n-1]
    double p = 1.0;
    if (g0>0.0&&g1>0.0) p = log(g1/g0)/log((xx[n+1]-xx[n-1])/(xx[n+1]-xx[n]));
    if (p<=1.0) p=1.0;
    g=g0*pow((xx[n+1]-xxx)/(xx[n+1]-xx[n]),p);
  }

  // Usual interpolation.
  else {
    g=0.0;
    for (l=4;l>=1;l--) {
      g=t*g+((c[ip][n][m][l][4]*u+c[ip][n][m][l][3])*u
         +c[ip][n][m][l][2])*u+c[ip][n][m][l][1];
    }
  }

  return g;
}

//--------------------------------------------------------------------------

// Extrapolate PDF value outside data grid.


double MSTWpdf::parton_extrapolate(int ip, double xxx, double qqq) {

  double parton_pdf=0.;
  int n,m;

  n=locate(xx,nx,xxx); // 0: below xmin, nx: above xmax
  m=locate(qq,nq,qqq); // 0: below qsqmin, nq: above qsqmax

  if (n==0&&(m>0&&m<nq)) { // if extrapolation in small x only

    double f0,f1;
    f0=parton_interpolate(ip,xx[1],qqq);
    f1=parton_interpolate(ip,xx[2],qqq);
    if ( f0>1e-3 && f1>1e-3 ) { // if values are positive, keep them so
      f0=log(f0);
      f1=log(f1);
      parton_pdf=exp(f0+(f1-f0)/(xx[2]-xx[1])*(xxx-xx[1]));
    } else // otherwise just extrapolate in the value
      parton_pdf=f0+(f1-f0)/(xx[2]-xx[1])*(xxx-xx[1]);

  } if (n>0&&m==nq) { // if extrapolation into large q only

    double f0,f1;
    f0=parton_interpolate(ip,xxx,qq[nq]);
    f1=parton_interpolate(ip,xxx,qq[nq-1]);
    if ( f0>1e-3 && f1>1e-3 ) { // if values are positive, keep them so
      f0=log(f0);
      f1=log(f1);
      parton_pdf=exp(f0+(f0-f1)/(qq[nq]-qq[nq-1])*(qqq-qq[nq]));
    } else // otherwise just extrapolate in the value
      parton_pdf=f0+(f0-f1)/(qq[nq]-qq[nq-1])*(qqq-qq[nq]);

  } if (n==0&&m==nq) { // if extrapolation into large q AND small x

    double f0,f1;
    f0=parton_extrapolate(ip,xx[1],qqq);
    f1=parton_extrapolate(ip,xx[2],qqq);
    if ( f0>1e-3 && f1>1e-3 ) { // if values are positive, keep them so
      f0=log(f0);
      f1=log(f1);
      parton_pdf=exp(f0+(f1-f0)/(xx[2]-xx[1])*(xxx-xx[1]));
    } else // otherwise just extrapolate in the value
      parton_pdf=f0+(f1-f0)/(xx[2]-xx[1])*(xxx-xx[1]);

  }

  return parton_pdf;
}

//--------------------------------------------------------------------------

// Returns an integer j such that x lies inbetween xloc[j] and xloc[j+1].
// unit offset of increasing ordered array xloc assumed.
// n is the length of the array (xloc[n] highest element).

int MSTWpdf::locate(double xloc[],int n,double x) {
  int ju,jm,jl(0),j;
  ju=n+1;

  while (ju-jl>1) {
    jm=(ju+jl)/2; // compute a mid point.
    if ( x>= xloc[jm])
      jl=jm;
    else ju=jm;
  }
  if (x==xloc[1]) j=1;
  else if (x==xloc[n]) j=n-1;
  else j=jl;

  return j;
}

//--------------------------------------------------------------------------

// Returns the estimate of the derivative at x1 obtained by a polynomial
// interpolation using the three points (x_i,y_i).

double MSTWpdf::polderivative1(double x1, double x2, double x3, double y1,
  double y2, double y3) {

  return (x3*x3*(y1-y2)+2.0*x1*(x3*(-y1+y2)+x2*(y1-y3))+x2*x2*(-y1+y3)
          +x1*x1*(-y2+y3))/((x1-x2)*(x1-x3)*(x2-x3));

}

//--------------------------------------------------------------------------

// Returns the estimate of the derivative at x2 obtained by a polynomial
// interpolation using the three points (x_i,y_i).

double MSTWpdf::polderivative2(double x1, double x2, double x3, double y1,
  double y2, double y3) {

  return (x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*(y2-y3))+x2*x2*(y1-y3)
          +x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3));

}

//--------------------------------------------------------------------------

// Returns the estimate of the derivative at x3 obtained by a polynomial
// interpolation using the three points (x_i,y_i).

double MSTWpdf::polderivative3(double x1, double x2, double x3, double y1,
  double y2, double y3) {

  return (x3*x3*(-y1+y2)+2.0*x2*x3*(y1-y3)+x1*x1*(y2-y3)+x2*x2*(-y1+y3)
          +2.0*x1*x3*(-y2+y3))/((x1-x2)*(x1-x3)*(x2-x3));

}

//==========================================================================

// The CTEQ6pdf class.
// Code for handling CTEQ6L, CTEQ6L1, CTEQ66.00, CT09MC1, CT09MC2, CT09MCS.
// Also handles ACTW Pomeron sets B, D and SG (alpha = 1.14) and D (= 1.19).

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Stay away from xMin, xMax, Qmin, Qmax limits.
const double CTEQ6pdf::EPSILON = 1e-6;

// Assumed approximate power of small-x behaviour for interpolation.
const double CTEQ6pdf::XPOWER = 0.3;

//--------------------------------------------------------------------------

// Initialize PDF: select data file and open stream.

void CTEQ6pdf::init(int iFitIn, string xmlPath, Info* infoPtr) {

  // Choice of fit among possibilities.
  iFit = iFitIn;

  // Select which data file to read for current fit.
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";
  string fileName = "  ";
  if (iFit == 1) fileName = "cteq6l.tbl";
  if (iFit == 2) fileName = "cteq6l1.tbl";
  if (iFit == 3) fileName = "ctq66.00.pds";
  if (iFit == 4) fileName = "ct09mc1.pds";
  if (iFit == 5) fileName = "ct09mc2.pds";
  if (iFit == 6) fileName = "ct09mcs.pds";
  // Ditto for current Pomeron fit.
  if (iFit == 11) fileName = "pomactwb14.pds";
  if (iFit == 12) fileName = "pomactwd14.pds";
  if (iFit == 13) fileName = "pomactwsg14.pds";
  if (iFit == 14) fileName = "pomactwd19.pds";
  bool isPdsGrid = (iFit > 2);

  // Open data file.
  ifstream pdfgrid( (xmlPath + fileName).c_str() );
  if (!pdfgrid.good()) {
    printErr("Error in CTEQ6pdf::init: did not find data file", infoPtr);
    isSet = false;
    return;
  }

  // Initialization with a stream.
  init( pdfgrid, isPdsGrid, infoPtr);
  pdfgrid.close();

}

//--------------------------------------------------------------------------

// Initialize PDF: read in data grid from stream and set up interpolation.

void CTEQ6pdf::init(istream& pdfgrid, bool isPdsGrid, Info* infoPtr) {

  // Check that data stream is available.
  if (!pdfgrid.good()) {
    printErr("Error in CTEQ6pdf::init: cannot read from stream", infoPtr);
    isSet = false;
    return;
  }

  // Read in common information.
  int    iDum;
  double orderTmp, nQTmp, qTmp, rDum;
  string line;
  getline( pdfgrid, line);
  getline( pdfgrid, line);
  getline( pdfgrid, line);
  istringstream is1(line);
  is1 >> orderTmp >> nQTmp >> lambda >> mQ[1] >> mQ[2] >> mQ[3]
     >> mQ[4] >> mQ[5] >> mQ[6];
  order  = int(orderTmp + 0.5);
  nQuark = int(nQTmp + 0.5);
  getline( pdfgrid, line);

  // Read in information for the .pds grid format.
  if (isPdsGrid) {
    getline( pdfgrid, line);
    istringstream is2(line);
    is2 >> iDum >> iDum >> iDum >> nfMx >> mxVal >> iDum;
    if (mxVal > 4) mxVal = 3;
    getline( pdfgrid, line);
    getline( pdfgrid, line);
    istringstream is3(line);
    is3 >> nX >> nT >> iDum >> nG >> iDum;
    for (int i = 0; i < nG + 2; ++i) getline( pdfgrid, line);
    getline( pdfgrid, line);
    istringstream is4(line);
    is4 >> qIni >> qMax;
    for (int iT = 0; iT <= nT; ++iT) {
      getline( pdfgrid, line);
      istringstream is5(line);
      is5 >> qTmp;
      tv[iT] = log( log( qTmp/lambda));
    }
    getline( pdfgrid, line);
    getline( pdfgrid, line);
    istringstream is6(line);
    is6 >> xMin >> rDum;
    int nPackX = 6;
    xv[0] = 0.;
    for (int iXrng = 0; iXrng < int( (nX + nPackX - 1) / nPackX); ++iXrng) {
      getline( pdfgrid, line);
      istringstream is7(line);
      for (int iX = nPackX * iXrng + 1; iX <= nPackX * (iXrng + 1); ++iX)
      if (iX <= nX) is7 >> xv[iX];
    }
  }

  // Read in information for the .tbl grid format.
  else {
    mxVal = 2;
    getline( pdfgrid, line);
    istringstream is2(line);
    is2 >> nX >> nT >> nfMx;
    getline( pdfgrid, line);
    getline( pdfgrid, line);
    istringstream is3(line);
    is3 >> qIni >> qMax;
    int    nPackT = 6;
    for (int iTrng = 0; iTrng < int( (nT + nPackT) / nPackT); ++iTrng) {
      getline( pdfgrid, line);
      istringstream is4(line);
      for (int iT = nPackT * iTrng; iT < nPackT * (iTrng + 1); ++iT)
      if (iT <= nT) {
        is4 >> qTmp;
        tv[iT] = log( log( qTmp / lambda) );
      }
    }
    getline( pdfgrid, line);
    getline( pdfgrid, line);
    istringstream is5(line);
    is5 >> xMin;
    int nPackX = 6;
    for (int iXrng = 0; iXrng < int( (nX + nPackX) / nPackX); ++iXrng) {
      getline( pdfgrid, line);
      istringstream is6(line);
      for (int iX = nPackX * iXrng; iX < nPackX * (iXrng + 1); ++iX)
      if (iX <= nX) is6 >> xv[iX];
    }
  }

  // Read in the grid proper.
  getline( pdfgrid, line);
  int nBlk  = (nX + 1) * (nT + 1);
  int nPts  = nBlk * (nfMx + 1 + mxVal);
  int nPack = (isPdsGrid) ? 6 : 5;
  for (int iRng = 0; iRng < int( (nPts + nPack - 1) / nPack); ++iRng) {
    getline( pdfgrid, line);
    istringstream is8(line);
    for (int i = nPack * iRng + 1; i <= nPack * (iRng + 1); ++i)
      if (i <= nPts) is8 >> upd[i];
  }

  // Initialize x grid mapped to x^0.3.
  xvpow[0] = 0.;
  for (int iX = 1; iX <= nX; ++iX)  xvpow[iX] = pow(xv[iX], XPOWER);

  // Set x and Q borders with some margin.
  xMinEps = xMin * (1. + EPSILON);
  xMaxEps = 1. - EPSILON;
  qMinEps = qIni * (1. + EPSILON);
  qMaxEps = qMax * (1. - EPSILON);

  // Initialize (x, Q) values of previous call.
  xLast = 0.;
  qLast = 0.;

}

//--------------------------------------------------------------------------

// Update PDF values.

void CTEQ6pdf::xfUpdate(int , double x, double Q2) {

  // Update using CTEQ6 routine, within allowed (x, q) range.
  double xEps = doExtraPol ? x : max( xMinEps, x);
  double qEps = max( qMinEps, min( qMaxEps, sqrtpos(Q2) ) );

  // Gluon:
  double glu  = xEps * parton6( 0, xEps, qEps);
  // Sea quarks (note wrong order u, d). ACTW has no s, c, b.
  double bot  = (iFit > 10) ? 0. : xEps * parton6( 5, xEps, qEps);
  double chm  = (iFit > 10) ? 0. : xEps * parton6( 4, xEps, qEps);
  double str  = xEps * parton6( 3, xEps, qEps);
  double usea = xEps * parton6(-1, xEps, qEps);
  double dsea = xEps * parton6(-2, xEps, qEps);
  // Valence quarks:
  double upv  = xEps * parton6( 1, xEps, qEps) - usea;
  double dnv  = xEps * parton6( 2, xEps, qEps) - dsea;

  // Check that rescaling *only* occurs for ACTW Pomeron PDF sets.
  if (iFit < 10) rescale = 1.;

  // Transfer to Pythia notation.
  xg     = rescale * glu;
  xu     = rescale * (upv + usea);
  xd     = rescale * (dnv + dsea);
  xubar  = rescale * usea;
  xdbar  = rescale * dsea;
  xs     = rescale * str;
  xsbar  = rescale * str;
  xc     = rescale * chm;
  xb     = rescale * bot;
  xgamma = 0.;

  // Subdivision of valence and sea.
  xuVal  = rescale * upv;
  xuSea  = rescale * usea;
  xdVal  = rescale * dnv;
  xdSea  = rescale * dsea;

  // idSav = 9 to indicate that all flavours reset.
  idSav  = 9;

}

//--------------------------------------------------------------------------

// Returns the PDF value for parton of flavour iParton at x, q.

double CTEQ6pdf::parton6(int iParton, double x, double q) {

  // Put zero for large x. Parton table and interpolation variables.
  if (x > xMaxEps) return 0.;
  int    iP = (iParton > mxVal) ? -iParton : iParton;
  double ss = pow( x, XPOWER);
  double tt = log( log(q / lambda) );

  // Find location in grid.Skip if same as in latest call.
  if (x != xLast || q != qLast) {

    // Binary search in x grid.
    iGridX  = 0;
    iGridLX = -1;
    int ju  = nX + 1;
    int jm  = 0;
    while (ju - iGridLX > 1 && jm >= 0) {
      jm = (ju + iGridLX) / 2;
      if (x >= xv[jm]) iGridLX = jm;
      else ju = jm;
    }

    // Separate acceptable from unacceptable grid points.
    if      (iGridLX <= -1)     return 0.;
    else if (iGridLX == 0)      iGridX = 0;
    else if (iGridLX <= nX - 2) iGridX = iGridLX - 1;
    else if (iGridLX == nX - 1) iGridX = iGridLX - 2;
    else                        return 0.;

    // Expressions for interpolation in x Grid.
    if (iGridLX > 1 && iGridLX < nX - 1) {
      double svec1 = xvpow[iGridX];
      double svec2 = xvpow[iGridX+1];
      double svec3 = xvpow[iGridX+2];
      double svec4 = xvpow[iGridX+3];
      double s12   = svec1 - svec2;
      double s13   = svec1 - svec3;
      xConst[8]    = svec2 - svec3;
      double s24   = svec2 - svec4;
      double s34   = svec3 - svec4;
      xConst[6]    = ss - svec2;
      xConst[7]    = ss - svec3;
      xConst[0]    = s13 / xConst[8];
      xConst[1]    = s12 / xConst[8];
      xConst[2]    = s34 / xConst[8];
      xConst[3]    = s24 / xConst[8];
      double s1213 = s12 + s13;
      double s2434 = s24 + s34;
      double sdet  = s12 * s34 - s1213 * s2434;
      double tmp   = xConst[6] * xConst[7] / sdet;
      xConst[4]    = (s34 * xConst[6] - s2434 * xConst[7]) * tmp / s12;
      xConst[5]    = (s1213 * xConst[6] - s12 * xConst[7]) * tmp / s34;
    }

    // Expression for extrapolation in x Grid.
    dlx = (iGridLX == 0 && doExtraPol)
        ? log(x / xv[1]) / log(xv[2] / xv[1]) : 1.;

    // Binary search in Q grid.
    iGridQ  = 0;
    iGridLQ = -1;
    ju      = nT + 1;
    jm      = 0;
    while (ju - iGridLQ > 1 && jm >= 0) {
      jm = (ju + iGridLQ) / 2;
      if (tt >= tv[jm]) iGridLQ = jm;
      else ju = jm;
    }
    if      (iGridLQ == 0)      iGridQ = 0;
    else if (iGridLQ <= nT - 2) iGridQ = iGridLQ - 1;
    else                        iGridQ = nT - 3;

    // Expressions for interpolation in Q Grid.
    if (iGridLQ > 0 && iGridLQ < nT - 1) {
      double tvec1 = tv[iGridQ];
      double tvec2 = tv[iGridQ+1];
      double tvec3 = tv[iGridQ+2];
      double tvec4 = tv[iGridQ+3];
      double t12   = tvec1 - tvec2;
      double t13   = tvec1 - tvec3;
      tConst[8]    =   tvec2 - tvec3;
      double t24   = tvec2 - tvec4;
      double t34   = tvec3 - tvec4;
      tConst[6]    = tt - tvec2;
      tConst[7]    = tt - tvec3;
      double tmp1  = t12 + t13;
      double tmp2  = t24 + t34;
      double tdet  = t12 * t34 - tmp1 * tmp2;
      tConst[0]    = t13 / tConst[8];
      tConst[1]    = t12 / tConst[8];
      tConst[2]    = t34 / tConst[8];
      tConst[3]    = t24 / tConst[8];
      tConst[4]    = (t34 * tConst[6] - tmp2 * tConst[7]) / t12
                     * tConst[6] * tConst[7] / tdet;
      tConst[5]    = (tmp1 * tConst[6] - t12 * tConst[7]) / t34
                     * tConst[6] * tConst[7] / tdet;
    }

    // Save x and q values so do not have to redo same again.
    xLast = x;
    qLast = q;
  }

  // Jump to here if x and q are the same as for the last call.
  int jtmp = ( (iP + nfMx) * (nT + 1) + (iGridQ - 1) ) * (nX + 1) + iGridX + 1;

  // Interpolate in x space for four different q values.
  // Also option for extrapolation to small x values.
  for(int it = 1; it <= 4; ++it) {
    int j1 = jtmp + it * (nX + 1);
    if (iGridLX == 0 && doExtraPol) {
      fVec[it] = upd[j1+1] * pow( upd[j1+2] / upd[j1+1], dlx );
    } else if (iGridX == 0) {
      double fij[5];
      fij[1] = 0.;
      fij[2] = upd[j1+1] * pow2(xv[1]);
      fij[3] = upd[j1+2] * pow2(xv[2]);
      fij[4] = upd[j1+3] * pow2(xv[3]);
      double fX = polint4F( &xvpow[0], &fij[1], ss);
      fVec[it] = (x > 0.) ? fX / pow2(x) : 0.;
    } else if (iGridLX==nX-1) {
      fVec[it] = polint4F( &xvpow[nX-3], &upd[j1], ss);
    } else {
      double sf2 = upd[j1+1];
      double sf3 = upd[j1+2];
      double g1  =  sf2 * xConst[0] - sf3 * xConst[1];
      double g4  = -sf2 * xConst[2] + sf3 * xConst[3];
      fVec[it]   = (xConst[4] * (upd[j1] - g1) + xConst[5] * (upd[j1+3] - g4)
                 + sf2 * xConst[7] - sf3 * xConst[6]) / xConst[8];
    }
  }

  // Interpolate in q space for x-interpolated values found above.
  double ff;
  if( iGridLQ <= 0 ) {
    ff = polint4F( &tv[0], &fVec[1], tt);
  } else if (iGridLQ >= nT - 1) {
    ff=polint4F( &tv[nT-3], &fVec[1], tt);
  } else {
    double tf2 = fVec[2];
    double tf3 = fVec[3];
    double g1  =  tf2 * tConst[0] - tf3 * tConst[1];
    double g4  = -tf2 * tConst[2] + tf3 * tConst[3];
    ff         = (tConst[4] * (fVec[1] - g1) + tConst[5] * (fVec[4] - g4)
               + tf2 * tConst[7] - tf3 * tConst[6]) / tConst[8];
  }

  // Done.
  return ff;
}

//--------------------------------------------------------------------------

// The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
// but assuming N=4, and ignoring the error estimation.
// Suggested by Z. Sullivan.

double CTEQ6pdf::polint4F(double xa[],double ya[],double x) {

  double y, h1, h2, h3, h4, w, den, d1, c1, d2, c2, d3, c3, cd1, cc1,
         cd2, cc2, dd1, dc1;

  h1  = xa[0] - x;
  h2  = xa[1] - x;
  h3  = xa[2] - x;
  h4  = xa[3] - x;

  w   = ya[1] - ya[0];
  den = w / (h1 - h2);
  d1  = h2 * den;
  c1  = h1 * den;

  w   = ya[2] - ya[1];
  den = w / (h2 - h3);
  d2  = h3 * den;
  c2  = h2 * den;

  w   = ya[3] - ya[2];
  den = w / (h3 - h4);
  d3  = h4 * den;
  c3  = h3 * den;

  w   = c2 - d1;
  den = w / (h1 - h3);
  cd1 = h3 * den;
  cc1 = h1 * den;

  w   = c3 - d2;
  den = w / (h2 - h4);
  cd2 = h4 * den;
  cc2 = h2 * den;

  w   = cc2 - cd1;
  den = w / (h1 - h4);
  dd1 = h4 * den;
  dc1 = h1 * den;

  if      (h3 + h4 < 0.) y = ya[3] + d3 + cd2 + dd1;
  else if (h2 + h3 < 0.) y = ya[2] + d2 + cd1 + dc1;
  else if (h1 + h2 < 0.) y = ya[1] + c2 + cd1 + dc1;
  else                   y = ya[0] + c1 + cc1 + dc1;

  return y;

}

//==========================================================================

// SA Unresolved proton: equivalent photon spectrum from
// V.M. Budnev, I.F. Ginzburg, G.V. Meledin and V.G. Serbo,
// Phys. Rept. 15 (1974/1975) 181.

// Constants:
const double ProtonPoint::ALPHAEM = 0.00729735;
const double ProtonPoint::Q2MAX   = 2.0;
const double ProtonPoint::Q20     = 0.71;
const double ProtonPoint::A       = 7.16;
const double ProtonPoint::B       = -3.96;
const double ProtonPoint::C       = 0.028;

//--------------------------------------------------------------------------

// Gives a generic Q2-independent equivalent photon spectrum.

void ProtonPoint::xfUpdate(int , double x, double /*Q2*/ ) {

  // Photon spectrum
  double tmpQ2Min = 0.88 * pow2(x) / (1. - x);
  double phiMax = phiFunc(x, Q2MAX / Q20);
  double phiMin = phiFunc(x, tmpQ2Min / Q20);

  double fgm = 0.;
  if (phiMax < phiMin) {
    printErr("Error in ProtonPoint::xfUpdate: phiMax - phiMin < 0!", infoPtr);
  } else {
    // Corresponds to: x*f(x)
    fgm = (ALPHAEM / M_PI) * (1 - x) * (phiMax - phiMin);
  }

  // Update values
  xg     = 0.;
  xu     = 0.;
  xd     = 0.;
  xubar  = 0.;
  xdbar  = 0.;
  xs     = 0.;
  xsbar  = 0.;
  xc     = 0.;
  xb     = 0.;
  xgamma = fgm;

  // Subdivision of valence and sea.
  xuVal = 0.;
  xuSea = 0;
  xdVal = 0.;
  xdSea = 0;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//--------------------------------------------------------------------------

// Function related to Q2 integration.

double ProtonPoint::phiFunc(double x, double Q) {

  double tmpV = 1. + Q;
  double tmpSum1 = 0;
  double tmpSum2 = 0;
  for (int k=1; k<4; ++k) {
    tmpSum1 += 1. / (k * pow(tmpV, k));
    tmpSum2 += pow(B, k) / (k * pow(tmpV, k));
  }

  double tmpY = pow2(x) / (1 - x);
  double funVal = (1 + A * tmpY) * (-1.*log(tmpV / Q) + tmpSum1)
                + (1 - B) * tmpY / (4 * Q * pow(tmpV, 3))
                + C * (1 + tmpY/4.)* (log((tmpV - B)/tmpV) + tmpSum2);

  return funVal;

}

//==========================================================================

// Gives the GRV 1992 pi+ (leading order) parton distribution function set
// in parametrized form. Authors: Glueck, Reya and Vogt.
// Ref: M. Glueck, E. Reya and A. Vogt, Z. Phys. C53 (1992) 651.
// Allowed variable range: 0.25 GeV^2 < Q^2 < 10^8 GeV^2 and 10^-5 < x < 1.

void GRVpiL::xfUpdate(int , double x, double Q2) {

  // Common expressions. Constrain Q2 for which parametrization is valid.
  double mu2  = 0.25;
  double lam2 = 0.232 * 0.232;
  double s    = (Q2 > mu2) ? log( log(Q2/lam2) / log(mu2/lam2) ) : 0.;
  double s2   = s * s;
  double x1   = 1. - x;
  double xL   = -log(x);
  double xS   = sqrt(x);

  // uv, dbarv.
  double uv = (0.519 + 0.180 * s - 0.011 * s2) * pow(x, 0.499 - 0.027 * s)
    * (1. + (0.381 - 0.419 * s) * xS) * pow(x1, 0.367 + 0.563 * s);

  // g.
  double gl = ( pow(x, 0.482 + 0.341 * sqrt(s))
    * ( (0.678 + 0.877 * s - 0.175 * s2) + (0.338 - 1.597 * s) * xS
    + (-0.233 * s + 0.406 * s2) * x) + pow(s, 0.599)
    * exp(-(0.618 + 2.070 * s) + sqrt(3.676 * pow(s, 1.263) * xL) ) )
    * pow(x1, 0.390 + 1.053 * s);

  // sea: u, d, s.
  double ub = pow(s, 0.55) * (1. - 0.748 * xS + (0.313 + 0.935 * s) * x)
    * pow(x1, 3.359) * exp(-(4.433 + 1.301 * s) + sqrt((9.30 - 0.887 * s)
    * pow(s, 0.56) * xL) ) / pow(xL, 2.538 - 0.763 * s);

  // c.
  double chm = (s < 0.888) ? 0. : pow(s - 0.888, 1.02) * (1. + 1.008 * x)
    * pow(x1, 1.208 + 0.771 * s) * exp(-(4.40 + 1.493 * s)
    + sqrt( (2.032 + 1.901 * s) * pow(s, 0.39) * xL) );

  // b.
  double bot = (s < 1.351) ? 0. : pow(s - 1.351, 1.03)
    * pow(x1, 0.697 + 0.855 * s) * exp(-(4.51 + 1.490 * s)
    + sqrt( (3.056 + 1.694 * s) * pow(s, 0.39) * xL) );

  // Update values.
  xg    = rescale * gl;
  xu    = rescale * (uv + ub);
  xd    = rescale * ub;
  xubar = rescale * ub;
  xdbar = rescale * (uv + ub);
  xs    = rescale * ub;
  xsbar = rescale * ub;
  xc    = rescale * chm;
  xb    = rescale * bot;

  // Subdivision of valence and sea.
  xuVal = rescale * uv;
  xuSea = rescale * ub;
  xdVal = rescale * uv;
  xdSea = rescale * ub;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//==========================================================================

// Pomeron PDF: simple Q2-independent parametrizations N x^a (1 - x)^b.

//--------------------------------------------------------------------------

// Calculate normalization factors once and for all.

void PomFix::init() {

  normGluon = GammaReal(PomGluonA + PomGluonB + 2.)
            / (GammaReal(PomGluonA + 1.) * GammaReal(PomGluonB + 1.));
  normQuark = GammaReal(PomQuarkA + PomQuarkB + 2.)
            / (GammaReal(PomQuarkA + 1.) * GammaReal(PomQuarkB + 1.));

}

//--------------------------------------------------------------------------

// Gives a generic Q2-independent Pomeron PDF.

void PomFix::xfUpdate(int , double x, double) {

  // Gluon and quark distributions.
  double gl = normGluon * pow(x, PomGluonA) * pow( (1. - x), PomGluonB);
  double qu = normQuark * pow(x, PomQuarkA) * pow( (1. - x), PomQuarkB);

  // Update values
  xg    = (1. - PomQuarkFrac) * gl;
  xu    = (PomQuarkFrac / (4. + 2. * PomStrangeSupp) ) * qu;
  xd    = xu;
  xubar = xu;
  xdbar = xu;
  xs    = PomStrangeSupp * xu;
  xsbar = xs;
  xc    = 0.;
  xb    = 0.;

  // Subdivision of valence and sea.
  xuVal = 0.;
  xuSea = xu;
  xdVal = 0.;
  xdSea = xd;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//==========================================================================

// Pomeron PDF: the H1 2006 Fit A and Fit B Q2-dependent parametrizations.

//--------------------------------------------------------------------------

// Initialize PDF: select data file and open stream.

void PomH1FitAB::init( int iFit, string xmlPath, Info* infoPtr) {

  // Open files from which grids should be read in.
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";
  string         dataFile = "pomH1FitBlo.data";
  if (iFit == 1) dataFile = "pomH1FitA.data";
  if (iFit == 2) dataFile = "pomH1FitB.data";
  ifstream is( (xmlPath + dataFile).c_str() );
  if (!is.good()) {
    printErr("Error in PomH1FitAB::init: did not find data file", infoPtr);
    isSet = false;
    return;
  }

  // Initialization with a stream.
  init( is, infoPtr );
  is.close();

}

//--------------------------------------------------------------------------

// Initialize PDF: read in data grid from stream and set up interpolation.

void PomH1FitAB::init( istream& is, Info* infoPtr) {

  // Check that data stream is available.
  if (!is.good()) {
    printErr("Error in PomH1FitAB::init: cannot read from stream", infoPtr);
    isSet = false;
    return;
  }

  // Lower and upper bounds. Bin widths for logarithmic spacing.
  nx    = 100;
  xlow  = 0.001;
  xupp  = 0.99;
  dx    = log(xupp / xlow) / (nx - 1.);
  nQ2   = 30;
  Q2low = 1.0;
  Q2upp = 30000.;
  dQ2   = log(Q2upp / Q2low) / (nQ2 - 1.);

  // Read in quark data grid.
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < nQ2; ++j)
      is >> quarkGrid[i][j];

  // Read in gluon data grid.
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < nQ2; ++j)
      is >> gluonGrid[i][j];

  // Check for errors during read-in of file.
  if (!is) {
    printErr("Error in PomH1FitAB::init: could not read data stream", infoPtr);
    isSet = false;
    return;
  }

  // Done.
  isSet = true;
  return;

}

//--------------------------------------------------------------------------

void PomH1FitAB::xfUpdate(int , double x, double Q2) {

  // Retrict input to validity range.
  double xt  = min( xupp, max( xlow, x) );
  double Q2t = min( Q2upp, max( Q2low, Q2) );

  // Lower grid point and distance above it.
  double dlx  = log( xt / xlow) / dx;
  int i       = min( nx - 2,  int(dlx) );
  dlx        -= i;
  double dlQ2 = log( Q2t / Q2low) / dQ2;
  int j       = min( nQ2 - 2, int(dlQ2) );
  dlQ2       -= j;

  // Extrapolate to small x values for quark and gluon PDF.
  double qu, gl;
  if (x < xlow && doExtraPol) {
    dlx = log( x / xlow) / dx;
    qu = (1. - dlQ2) * quarkGrid[0][j]
       * pow( quarkGrid[1][j] / quarkGrid[0][j], dlx)
       +        dlQ2 * quarkGrid[0][j + 1]
       * pow( quarkGrid[1][j + 1] / quarkGrid[0][j + 1], dlx);
    gl = (1. - dlQ2) * gluonGrid[0][j]
       * pow( gluonGrid[1][j] / gluonGrid[0][j], dlx)
       +        dlQ2 * gluonGrid[0][j + 1]
       * pow( gluonGrid[1][j + 1] / gluonGrid[0][j + 1], dlx);

  } else {
    // Interpolate to derive quark PDF.
    qu = (1. - dlx) * (1. - dlQ2) * quarkGrid[i][j]
       +       dlx  * (1. - dlQ2) * quarkGrid[i + 1][j]
       + (1. - dlx) *       dlQ2  * quarkGrid[i][j + 1]
       +       dlx  *       dlQ2  * quarkGrid[i + 1][j + 1];

    // Interpolate to derive gluon PDF.
    gl = (1. - dlx) * (1. - dlQ2) * gluonGrid[i][j]
       +       dlx  * (1. - dlQ2) * gluonGrid[i + 1][j]
       + (1. - dlx) *       dlQ2  * gluonGrid[i][j + 1]
       +       dlx  *       dlQ2  * gluonGrid[i + 1][j + 1];
  }

  // Update values.
  xg    = rescale * gl;
  xu    = rescale * qu;
  xd    = xu;
  xubar = xu;
  xdbar = xu;
  xs    = xu;
  xsbar = xu;
  xc    = 0.;
  xb    = 0.;

  // Subdivision of valence and sea.
  xuVal = 0.;
  xuSea = xu;
  xdVal = 0.;
  xdSea = xu;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//==========================================================================

// Pomeron PDF: the H1 2007 Jets Q2-dependent parametrization.

//--------------------------------------------------------------------------

// Initialize PDF: select data file and open stream.

void PomH1Jets::init( int , string xmlPath, Info* infoPtr) {

  // Open files from which grids should be read in.
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";
  ifstream is( (xmlPath + "pomH1Jets.data").c_str() );
  if (!is.good()) {
    printErr("Error in PomH1Jets::init: did not find data file", infoPtr);
    isSet = false;
    return;
  }

  // Initialization with a stream.
  init( is, infoPtr);
  is.close();

}

//--------------------------------------------------------------------------

// Initialize PDF: read in data grid from stream and set up interpolation.

void PomH1Jets::init( istream& is, Info* infoPtr) {

  // Check that data stream is available.
  if (!is.good()) {
    printErr("Error in PomH1Jets::init: cannot read from stream", infoPtr);
    isSet = false;
    return;
  }

  // Read in x and Q grids. Do interpolation logarithmically in Q2.
  for (int i = 0; i < 100; ++i) {
    is >> setw(13) >> xGrid[i];
  }
  for (int j = 0; j < 88; ++j) {
    is >> setw(13) >> Q2Grid[j];
    Q2Grid[j] = log( Q2Grid[j] );
  }

  // Read in  gluon data grid.
  for (int j = 0; j < 88; ++j) {
    for (int i = 0; i < 100; ++i) {
      is >> setw(13) >> gluonGrid[i][j];
    }
  }

  // Read in singlet data grid.
  for (int j = 0; j < 88; ++j) {
    for (int i = 0; i < 100; ++i) {
      is >> setw(13) >> singletGrid[i][j];
    }
  }

  // Read in charm data grid.
  for (int j = 0; j < 88; ++j) {
    for (int i = 0; i < 100; ++i) {
      is >> setw(13) >> charmGrid[i][j];
    }
  }

  // Check for errors during read-in of files.
  if (!is) {
    printErr("Error in PomH1Jets::init: could not read data file", infoPtr);
    isSet = false;
    return;
  }

  // Done.
  isSet = true;

}

//--------------------------------------------------------------------------

void PomH1Jets::xfUpdate(int , double x, double Q2) {

  // Find position in x array.
  double xLog = log(x);
  int    i    = 0;
  double dx   = 0.;
  if (xLog <= xGrid[0]);
  else if (xLog >= xGrid[99]) {
    i  = 98;
    dx = 1.;
  } else {
    while (xLog > xGrid[i]) ++i;
    --i;
    dx = (xLog - xGrid[i]) / (xGrid[i + 1] - xGrid[i]);
  }

  // Find position in y array.
  double Q2Log = log(Q2);
  int    j     = 0;
  double dQ2   = 0.;
  if (Q2Log <= Q2Grid[0]);
  else if (Q2Log >= Q2Grid[87]) {
    j   = 86;
    dQ2 = 1.;
  } else {
    while (Q2Log > Q2Grid[j]) ++j;
    --j;
    dQ2 = (Q2Log - Q2Grid[j]) / (Q2Grid[j + 1] - Q2Grid[j]);
  }

  // Extrapolate to small x values for gluon, singlet and charm PDF.
  double gl, sn, ch;
  if (xLog < xGrid[0] && doExtraPol) {
    double dlx = (xLog - xGrid[0]) / (xGrid[1] - xGrid[0]) ;
    gl = (1. - dQ2) * gluonGrid[0][j]
       * pow( gluonGrid[1][j] / gluonGrid[0][j], dlx)
       +        dQ2 * gluonGrid[0][j + 1]
       * pow( gluonGrid[1][j + 1] / gluonGrid[0][j + 1], dlx);
    sn = (1. - dQ2) * singletGrid[0][j]
       * pow( singletGrid[1][j] / singletGrid[0][j], dlx)
       +        dQ2 * singletGrid[0][j + 1]
       * pow( singletGrid[1][j + 1] / singletGrid[0][j + 1], dlx);
    ch = (1. - dQ2) * charmGrid[0][j]
       * pow( charmGrid[1][j] / charmGrid[0][j], dlx)
       +        dQ2 * charmGrid[0][j + 1]
       * pow( charmGrid[1][j + 1] / charmGrid[0][j + 1], dlx);

  } else {
    // Interpolate to derive gluon PDF.
    gl = (1. - dx) * (1. - dQ2) * gluonGrid[i][j]
       +       dx  * (1. - dQ2) * gluonGrid[i + 1][j]
       + (1. - dx) *       dQ2  * gluonGrid[i][j + 1]
       +       dx  *       dQ2  * gluonGrid[i + 1][j + 1];

    // Interpolate to derive singlet PDF. (Sum of u, d, s, ubar, dbar, sbar.)
    sn = (1. - dx) * (1. - dQ2) * singletGrid[i][j]
       +       dx  * (1. - dQ2) * singletGrid[i + 1][j]
       + (1. - dx) *       dQ2  * singletGrid[i][j + 1]
       +       dx  *       dQ2  * singletGrid[i + 1][j + 1];

    // Interpolate to derive charm PDF. (Charge-square times c and cbar.)
    ch = (1. - dx) * (1. - dQ2) * charmGrid[i][j]
       +       dx  * (1. - dQ2) * charmGrid[i + 1][j]
       + (1. - dx) *       dQ2  * charmGrid[i][j + 1]
       +       dx  *       dQ2  * charmGrid[i + 1][j + 1];
  }

  // Update values.
  xg    = rescale * gl;
  xu    = rescale * sn / 6.;
  xd    = xu;
  xubar = xu;
  xdbar = xu;
  xs    = xu;
  xsbar = xu;
  xc    = rescale * ch * 9./8.;
  xb    = 0.;

  // Subdivision of valence and sea.
  xuVal = 0.;
  xuSea = xu;
  xdVal = 0.;
  xdSea = xd;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//==========================================================================

void PomHISASD::xfUpdate(int, double x, double Q2) {

  // Check that pomeron momentum fraction is available.
  if ( xPomNow < 0.0 || xPomNow > 1.0 || !pPDFPtr )
    printErr("Error in PomHISASD::xfUpdate: no xPom available.", infoPtr);

  double xx = xPomNow * x;
  double fac = newfac * pow(1.0 - x, hixpow) / log(1.0 / xx);
  if ( fac == 0.0 ) fac = 1.0;

  xd = xdbar = fac * pPDFPtr->xfSea(1, xx, Q2);
  xu = xubar = fac * pPDFPtr->xfSea(2, xx, Q2);
  xs = xsbar = fac * pPDFPtr->xfSea(3, xx, Q2);
  xc =         fac * pPDFPtr->xfSea(4, xx, Q2);
  xb =         fac * pPDFPtr->xfSea(5, xx, Q2);
  xg =         fac * pPDFPtr->xfSea(21, xx, Q2);
  xlepton = xgamma = 0.0;

  // Subdivision of valence and sea.
  xuVal = 0.;
  xuSea = xu;
  xdVal = 0.;
  xdSea = xd;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//==========================================================================

// Gives electron (or muon, or tau) parton distribution.

// Constants: alphaEM(0), m_e, m_mu, m_tau.
const double Lepton::ALPHAEM = 0.00729735;
const double Lepton::ME      = 0.0005109989;
const double Lepton::MMU     = 0.10566;
const double Lepton::MTAU    = 1.77699;

void Lepton::xfUpdate(int id, double x, double Q2) {

  // Squared mass of lepton species: electron, muon, tau.
  if (!isInit) {
    double             mLep = ME;
    if (abs(id) == 13) mLep = MMU;
    if (abs(id) == 15) mLep = MTAU;
    m2Lep  = pow2( mLep );
    isInit = true;
  }

  // Electron inside electron, see R. Kleiss et al., in Z physics at
  // LEP 1, CERN 89-08, p. 34
  double xLog = log(max(1e-10,x));
  double xMinusLog = log( max(1e-10, 1. - x) );
  double Q2Log = log( max(3., Q2/m2Lep) );
  double beta = (ALPHAEM / M_PI) * (Q2Log - 1.);
  double delta = 1. + (ALPHAEM / M_PI) * (1.5 * Q2Log + 1.289868)
    + pow2(ALPHAEM / M_PI) * (-2.164868 * Q2Log*Q2Log
    + 9.840808 * Q2Log - 10.130464);
  double fPrel =  beta * pow(1. - x, beta - 1.) * sqrtpos( delta )
     - 0.5 * beta * (1. + x) + 0.125 * beta*beta * ( (1. + x)
     * (-4. * xMinusLog + 3. * xLog) - 4. * xLog / (1. - x) - 5. - x);

  // Zero distribution for very large x and rescale it for intermediate.
  if (x > 1. - 1e-10) fPrel = 0.;
  else if (x > 1. - 1e-7) fPrel *= pow(1000.,beta) / (pow(1000.,beta) - 1.);
  xlepton = x * fPrel;

  // Photons with restricted virtuality.
  double sCM = infoPtr->s();
  double m2s = 4 * m2Lep / sCM;
  double Q2minGamma = 2. * m2Lep * pow2(x)
    / ( 1. - x - m2s + sqrt(1. - m2s) * sqrt( pow2(1. - x) - m2s ) );
  xgamma = (0.5 * ALPHAEM / M_PI) * (1. + pow2(1. - x))
         * log( Q2maxGamma / Q2minGamma );

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//==========================================================================

// The NNPDF class.
// Code for handling NNPDF2.3 QCD+QED LO
// Code provided by Juan Rojo and Stefano Carrazza.

//--------------------------------------------------------------------------

// Freeze PDFs below XMINGRID
const double NNPDF::fXMINGRID = 1e-9;

//--------------------------------------------------------------------------

// Initialize PDF: select data file and open stream.

void NNPDF::init(int iFitIn, string xmlPath, Info* infoPtr) {

  // Choice of fit among possibilities.
  iFit = iFitIn;

  // Select which data file to read for current fit.
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";
  string fileName = "  ";
  // NNPDF2.3 LO QCD+QED, for two values of alphas
  if (iFit == 1) fileName = "NNPDF23_lo_as_0130_qed_mem0.grid";
  if (iFit == 2) fileName = "NNPDF23_lo_as_0119_qed_mem0.grid";
  // NNPDF2.3 NLO QCD+QED
  if (iFit == 3) fileName = "NNPDF23_nlo_as_0119_qed_mc_mem0.grid";
  // NNPDF2.4 NLO QCD+QED
  if (iFit == 4) fileName = "NNPDF23_nnlo_as_0119_qed_mc_mem0.grid";

  // Open data file.
  fstream f;
  f.open( (xmlPath + fileName).c_str(),ios::in);
  if (f.fail()) {
    printErr("Error in NNPDF::init: did not find data file ", infoPtr);
    isSet = false;
    return;
  }

  // Initialization with a stream.
  init(f, infoPtr);
  f.close();

}

//--------------------------------------------------------------------------

// Initialize PDF: read in data grid from stream and set up interpolation.

void NNPDF::init(istream& f, Info* infoPtr) {

  // Check that data stream is available.
  if (!f.good()) {
    printErr("Error in NNPDF::init: cannot read from stream", infoPtr);
    isSet = false;
    return;
  }

  // Reading grid: removing header.
  string tmp;
  for (;;) {
    getline(f,tmp);
    if (tmp.find("NNPDF20intqed") != string::npos) {
      getline(f,tmp);
      break;
    }
  }

  // Get nx and x grid.
  f >> fNX;
  fXGrid = new double[fNX];
  for (int ix = 0; ix < fNX; ix++) f >> fXGrid[ix];
  fLogXGrid = new double[fNX];
  for (int ix = 0; ix < fNX; ix++) fLogXGrid[ix] = log(fXGrid[ix]);

  // Get nQ2 and Q2 grid (ignorming first value).
  f >> fNQ2;
  f >> tmp;
  fQ2Grid = new double[fNQ2];
  for (int iq = 0; iq < fNQ2; iq++) f >> fQ2Grid[iq];
  fLogQ2Grid = new double[fNQ2];
  for (int iq = 0; iq < fNQ2; iq++) fLogQ2Grid[iq] = log(fQ2Grid[iq]);

  // Prepare grid array.
  fPDFGrid = new double**[fNFL];
  for (int i = 0; i < fNFL; i++) {
    fPDFGrid[i] = new double*[fNX];
    for (int j = 0; j < fNX; j++) {
      fPDFGrid[i][j] = new double[fNQ2];
      for (int z = 0; z < fNQ2; z++) fPDFGrid[i][j][z] = 0.0;
    }
  }

  // Check values of number of grid entries.
  if (fNX<= 0 || fNX>100 || fNQ2<=0 || fNQ2>50) {
    cout << "Error in NNPDF::init, Invalid grid values" << endl
         << "fNX = " << fNX << endl << "fNQ2 = " << fNQ2 << endl
         << "fNFL = " <<fNFL << endl;
    isSet = false;
    return;
  }

  // Ignore replica number. Read PDF grid points.
  f >> tmp;
  for (int ix = 0; ix < fNX; ix++)
    for (int iq = 0; iq < fNQ2; iq++)
      for (int fl = 0; fl < fNFL; fl++)
        f >> fPDFGrid[fl][ix][iq];

  // Other vectors.
  fRes = new double[fNFL];

}

//--------------------------------------------------------------------------

void NNPDF::xfUpdate(int , double x, double Q2) {

  // Update using NNPDF routine, within allowed (x, q) range.
  xfxevolve(x,Q2);

  // Then transfer to Pythia8 notation.
  xg     = fRes[6];
  xu     = fRes[8];
  xd     = fRes[7];
  xubar  = fRes[4];
  xdbar  = fRes[5];
  xs     = fRes[9];
  xsbar  = fRes[3];
  xc     = fRes[10];
  xb     = fRes[11];
  xgamma = fRes[13];

  // Subdivision of valence and sea.
  xuVal  = xu - xubar;
  xuSea  = xubar;
  xdVal  = xd - xdbar;
  xdSea  = xdbar;

  // idSav = 9 to indicate that all flavours reset.
  idSav  = 9;

}

//--------------------------------------------------------------------------

void NNPDF::xfxevolve(double x, double Q2) {

  // Freeze outside x-Q2 grid.
  if (x < fXMINGRID || x > fXGrid[fNX-1]) {
    if (x < fXMINGRID)  x = fXMINGRID;
    if (x > fXGrid[fNX-1]) x = fXGrid[fNX-1];
  }
  if (Q2 < fQ2Grid[0] || Q2 > fQ2Grid[fNQ2-1]) {
    if (Q2 < fQ2Grid[0]) Q2 = fQ2Grid[0];
    if (Q2 > fQ2Grid[fNQ2-1]) Q2 = fQ2Grid[fNQ2-1];
  }

  // Find nearest points in the x-Q2 grid.
  int minx = 0;
  int maxx = fNX;
  while (maxx-minx > 1) {
    int midx = (minx+maxx)/2;
    if (x < fXGrid[midx]) maxx = midx;
    else minx = midx;
  }
  int ix = minx;
  int minq = 0;
  int maxq = fNQ2;
  while (maxq-minq > 1) {
    int midq = (minq+maxq)/2;
    if (Q2 < fQ2Grid[midq]) maxq = midq;
    else minq = midq;
  }
  int iq2 = minq;

  // Assign grid for interpolation. M,N -> order of polyN interpolation.
  int    ix1a[fM], ix2a[fN];
  double x1a[fM], x2a[fN];
  double ya[fM][fN];

  for (int i = 0; i < fM; i++) {
    if (ix+1 >= fM/2 && ix+1 <= (fNX-fM/2)) ix1a[i] = ix+1 - fM/2 + i;
    if (ix+1 < fM/2) ix1a[i] = i;
    if (ix+1 > (fNX-fM/2)) ix1a[i] = (fNX-fM) + i;
    // Check grids.
    if (ix1a[i] < 0 || ix1a[i] >= fNX) {
      cout << "Error in grids! i, ixia[i] = " << i << "\t" << ix1a[i] << endl;
      return;
    }
  }

  for (int j = 0; j < fN; j++) {
    if (iq2+1 >= fN/2 && iq2+1 <= (fNQ2-fN/2)) ix2a[j] = iq2+1 - fN/2 + j;
    if (iq2+1 < fN/2) ix2a[j] = j;
    if (iq2+1 > (fNQ2-fN/2)) ix2a[j] = (fNQ2-fN) + j;
    // Check grids.
    if (ix2a[j] < 0 || ix2a[j] >= fNQ2) {
      cout << "Error in grids! j, ix2a[j] = " << j << "\t" << ix2a[j] << endl;
      return;
    }
  }

  const double xch = 1e-1;
  double x1;
  if (x < xch) x1 = log(x);
  else x1 = x;
  double x2 = log(Q2);

  for (int ipdf = 0; ipdf < fNFL; ipdf++) {
    fRes[ipdf] = 0.0;
    for (int i = 0; i < fM; i++) {
      if (x < xch) x1a[i] = fLogXGrid[ix1a[i]];
      else         x1a[i] = fXGrid[ix1a[i]];

      for (int j = 0; j < fN; j++) {
        x2a[j] = fLogQ2Grid[ix2a[j]];
        ya[i][j] = fPDFGrid[ipdf][ix1a[i]][ix2a[j]];
      }
    }

    // 2D polynomial interpolation.
    double y = 0, dy = 0;
    polin2(x1a,x2a,ya,x1,x2,y,dy);
    fRes[ipdf] = y;
  }

}

//--------------------------------------------------------------------------

// 1D polynomial interpolation.

void NNPDF::polint(double xa[], double yal[], int n, double x,
  double& y, double& dy) {

  int ns = 0;
  double dif = abs(x-xa[0]);
  double c[fM > fN ? fM : fN];
  double d[fM > fN ? fM : fN];

  for (int i = 0; i < n; i++) {
    double dift = abs(x-xa[i]);
    if (dift < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = yal[i];
    d[i] = yal[i];
  }
  y = yal[ns];
  ns--;
  for (int m = 1; m < n; m++) {
    for (int i = 0; i < n-m; i++) {
      double ho = xa[i]-x;
      double hp = xa[i+m]-x;
      double w = c[i+1]-d[i];
      double den = ho-hp;
      if (den == 0) {
        cout << "NNPDF::polint, failure" << endl;
        return;
      }
      den = w/den;
      d[i] = hp*den;
      c[i] = ho*den;
    }
    if (2*(ns+1) < n-m) dy = c[ns+1];
    else {
      dy = d[ns];
      ns--;
    }
    y+=dy;
  }
}

//--------------------------------------------------------------------------

// 2D polynomial interpolation.

void NNPDF::polin2(double x1al[], double x2al[], double yal[][fN],
  double x1, double x2, double& y, double& dy) {

  double yntmp[fN];
  double ymtmp[fM];

  for (int j = 0; j < fM; j++) {
    for (int k = 0; k < fN; k++) yntmp[k] = yal[j][k];
    polint(x2al,yntmp,fN,x2,ymtmp[j],dy);
  }
  polint(x1al,ymtmp,fM,x1,y,dy);

}

//==========================================================================

// LHAPDF plugin interface.

//--------------------------------------------------------------------------

// Constructor.

LHAPDF::LHAPDF(int idIn, string pSet, Info* infoPtrIn) :
  pdfPtr(0), infoPtr(infoPtrIn), lib(0) {
  isSet = false;
  if (!infoPtr) return;

  // Determine the plugin library name.
  if (pSet.size() < 8) {
    printErr("Error in LHAPDF::LHAPDF: invalid pSet " + pSet, infoPtr);
    return;
  }
  libName = pSet.substr(0, 7);
  if (libName != "LHAPDF5" && libName != "LHAPDF6") {
    printErr("Error in LHAPDF::LHAPDF: invalid pSet " + pSet, infoPtr);
    return;
  }
  libName = "libpythia8lhapdf" + libName.substr(6) + ".so";

  // Load the plugin library.
  const char* error(0);
  map<string, pair<void*, int> >::iterator plugin =
    infoPtr->plugins.find(libName);
  if (plugin == infoPtr->plugins.end()) {
    lib   = dlopen(libName.c_str(), RTLD_LAZY);
    error = dlerror();
  }
  if (error) {
    printErr("Error in LHAPDF::init: " + string(error), infoPtr);
    return;
  }
  if (plugin == infoPtr->plugins.end())
    infoPtr->plugins[libName] = pair<void*, int>(lib, 1);
  else {
    lib = plugin->second.first;
    ++plugin->second.second;
  }
  dlerror();

  // Determine the PDF set and member.
  string   set = pSet.substr(8);
  int      mem = 0;
  size_t   pos = set.find_last_of("/");
  if (pos != string::npos) {
    istringstream memStream(set.substr(pos + 1));
    memStream >> mem;
  }
  set = set.substr(0, pos);

  // Load the PDF.
  NewLHAPDF* newLHAPDF = (NewLHAPDF*)symbol("newLHAPDF");
  if (!newLHAPDF) return;
  pdfPtr = newLHAPDF(idIn, set, mem, infoPtr);
  isSet = true;

}

//--------------------------------------------------------------------------

// Destructor.

LHAPDF::~LHAPDF() {
  if (!infoPtr) return;
  if (!isSet)   return;

  // Delete the PDF.
  DeleteLHAPDF* deleteLHAPDF = (DeleteLHAPDF*)symbol("deleteLHAPDF");
  if (deleteLHAPDF) deleteLHAPDF(pdfPtr);

  // Close the plugin library if not needed by other instances.
  map<string, pair<void*, int> >::iterator plugin =
    infoPtr->plugins.find(libName);
  if (plugin == infoPtr->plugins.end()) return;
  --plugin->second.second;
  if (plugin->second.first && plugin->second.second == 0) {
    dlclose(plugin->second.first);
    dlerror();
    infoPtr->plugins.erase(plugin);
  }

}

//--------------------------------------------------------------------------

// Access a plugin library symbol.

LHAPDF::Symbol LHAPDF::symbol(string symName) {
  Symbol sym(0);
  const char* error(0);
  if (!infoPtr) return sym;

  // Load the symbol.
  sym = (Symbol)dlsym(lib, symName.c_str());
  error = dlerror();
  if (error) printErr("Error in LHAPDF::symbol: " + string(error), infoPtr);
  dlerror();
  return sym;

}

//==========================================================================

// Gives the CJKL leading order parton distribution function set
// in parametrized form for the real photons. Authors: F.Cornet, P.Jankowski,
// M.Krawczyk and A.Lorca, Phys. Rev. D68: 014010, 2003.
// Valid for 10^(-5) < x < 1 and 1 < Q^2 < 2*10^5 GeV^2.
// Below Q^2 = 1 a logarithmic approximation in Q^2 is used.

// Constants related to the fit.
const double CJKL::ALPHAEM = 0.007297353080;
const double CJKL::Q02     = 0.25;
const double CJKL::Q2MIN   = 0.05;
const double CJKL::Q2REF   = 1.0;
const double CJKL::LAMBDA  = 0.221;
const double CJKL::MC      = 1.3;
const double CJKL::MB      = 4.3;

//--------------------------------------------------------------------------

void CJKL::xfUpdate(int , double x, double Q2) {

  // Parameters:
  double lambda2 = pow2(LAMBDA);

  // When below reference scale calculate first with the reference scale and
  // later scale with log(Q^2).
  double Q2Save = Q2;
  bool belowRef = (Q2 < Q2REF);
  if ( belowRef) Q2 = Q2REF;

  // Evolution variable.
  double s = log( log(Q2/lambda2)/log(Q02/lambda2) );
  double plLog = 9.0/(4.0*M_PI)*log(Q2/lambda2);

  // Point-like contributions.
  double plGluon   = pointlikeG(x,s);
  double plUp      = pointlikeU(x,s);
  double plDown    = pointlikeD(x,s);
  double plStrange = plDown;

  // Hadron-like contributions.
  double hlGluon   = hadronlikeG(x,s);
  double hlVal     = hadronlikeVal(x,s);
  double hlSea     = hadronlikeSea(x,s);

  // Heavy quarks. Undo the ACOT_X rescaling for DIS kinematics.
  double xMaxC     = 1 - 6.76/(6.76 + Q2);
  double xMaxB     = 1 - 73.96/(73.96 + Q2);
  double plCharm   = pointlikeC(x*xMaxC,s,Q2)*xMaxC;
  double plBottom  = pointlikeB(x*xMaxB,s,Q2)*xMaxB;
  double hlCharm   = hadronlikeC(x*xMaxC,s,Q2)*xMaxC;
  double hlBottom  = hadronlikeB(x*xMaxB,s,Q2)*xMaxB;

  // Sum different contributions together.
  xg     = ALPHAEM*( plLog*plGluon   + hlGluon );
  xu     = ALPHAEM*( plLog*plUp      + 0.5*hlVal + hlSea );
  xd     = ALPHAEM*( plLog*plDown    + 0.5*hlVal + hlSea );
  xubar  = xu;
  xdbar  = xd;
  xs     = ALPHAEM*( plLog*plStrange + hlSea );
  xsbar  = xs;
  xc     = ALPHAEM*( plLog*plCharm   + hlCharm );
  xb     = ALPHAEM*( plLog*plBottom  + hlBottom );
  xgamma = 0;

  // Subdivision of valence and sea.
  xuVal = ALPHAEM*( plLog*plUp + 0.5*hlVal );
  xuSea = ALPHAEM*( hlSea );
  xdVal = ALPHAEM*( plLog*plDown + 0.5*hlVal );
  xdSea = ALPHAEM*( hlSea );
  xsVal = ALPHAEM*( plLog*plStrange );
  xsSea = ALPHAEM*( hlSea );
  xcVal = ALPHAEM*( plLog*plCharm );
  xcSea = ALPHAEM*( hlCharm );
  xbVal = ALPHAEM*( plLog*plBottom );
  xbSea = ALPHAEM*( hlBottom );

  // When below valid Q^2 values approximate scale evolution with log(Q^2).
  // Approximation derived by integrating xf over x and calculating the
  // derivative at Q2REF.
  if ( belowRef) {
    double logApprox = max( log(Q2Save/Q2MIN) / log(Q2REF/Q2MIN ), 0.);

    // Scale the PDFs according to log(Q^2) approx.
    xg    *= logApprox;
    xd    *= logApprox;
    xu    *= logApprox;
    xubar *= logApprox;
    xdbar *= logApprox;
    xs    *= logApprox;
    xsbar *= logApprox;
    xc    *= logApprox;
    xb    *= logApprox;
    xuVal *= logApprox;
    xuSea *= logApprox;
    xdVal *= logApprox;
    xdSea *= logApprox;
    xsVal *= logApprox;
    xsSea *= logApprox;
    xcVal *= logApprox;
    xcSea *= logApprox;
    xbVal *= logApprox;
    xbSea *= logApprox;
  }

  // idSav = 9 to indicate that all flavours reset.
  idSav  = 9;

}

//--------------------------------------------------------------------------

// Returns the x-dependence decoupled from the logarithmic scale
// dependence to approximate the PDFs from below for ISR.
// Currently flat in x (no second argument), could be improved.

double CJKL::gammaPDFxDependence(int id, double) {
  if      (abs(id) == 1) return 0.013 * ALPHAEM;
  else if (abs(id) == 2) return 0.026 * ALPHAEM;
  else if (abs(id) == 3) return 0.010 * ALPHAEM;
  else if (abs(id) == 4) return 0.020 * ALPHAEM;
  else if (abs(id) == 5) return 0.010 * ALPHAEM;
  else                   return 0;
}

//--------------------------------------------------------------------------

// Returns the reference scale for the logarithmic scale dependence to
// approximate the PDFs in ISR. Mass squared for heavy quarks and 0.2
// for others.

double CJKL::gammaPDFRefScale(int id) {
  if      (abs(id) == 4) return pow2(MC);
  else if (abs(id) == 5) return pow2(MB);
  else                   return 0.20;
}

//--------------------------------------------------------------------------

// Set valence content of the photon beam using parametrized Q2-dependence.

int CJKL::sampleGammaValFlavor(double Q2) {

  // Freeze the scale below the initial scale.
  if(Q2 < Q02) Q2 = Q02;

  // Calculate the x-integrated valence part of hadron-like contribution.
  double lambda2 = pow2(LAMBDA);
  double s  = log( log(Q2/lambda2)/log(Q02/lambda2) );
  double a  =  1.0898  + 0.38087 * s;
  double b  =  0.42654 - 1.2128  * s;
  double c  = -1.6576  + 1.7075  * s;
  double d  =  0.96155 + 1.8441  * s;
  double aa =  0.78391 - 0.06872 * s;
  double a1 = tgamma(1+aa)*tgamma(1+d)/tgamma(2+aa+d);
  double b1 = tgamma(1.5+aa)*tgamma(1+d)/tgamma(2.5+aa+d);
  double c1 = tgamma(2+aa)*tgamma(1+d)/tgamma(3+aa+d);
  double xfValHad = ALPHAEM*a*(a1 + b*b1 + c*c1);

  // Set the reference scales and charges.
  double mq2[5] = { Q02, Q02, Q02, pow2(MC), pow2(MB) };
  double eq2[5] = { 1.0/9.0, 4.0/9.0, 1.0/9.0, 4.0/9.0, 1.0/9.0 };

  // For u- and d-quarks valence contribution from hadron-like part.
  double qEvo[5] = { xfValHad/2, xfValHad/2, 0, 0, 0 };
  double qEvoTot = 0;

  // Normalization of the point-like part.
  double plNorm = 0.000936;

  // Logarithmic Q^2 evolution of gamma -> qqbar splitting for each flavor.
  for(int i = 0;i < 5;++i) {
    qEvo[i] += plNorm*eq2[i]*max(0.0,log(Q2/mq2[i]));
    qEvoTot += qEvo[i];
  }

  // Sample the valence flavor.
  double qEvoRand = qEvoTot*rndmPtr->flat();
  for(int i = 0; i < 5; ++i) {
    qEvoRand -= qEvo[i];
    if(qEvoRand <= 0.0) {
      idVal1 = i+1;
      idVal2 = -idVal1;
      break;
    }
  }

  return idVal1;
}

//--------------------------------------------------------------------------

// Sum of integrated PDFs \int dx x f(x,Q^2) at given scale Q^2.
// Integrals parametrized as a0 + a1*log(Q^2/Q0^2).

double CJKL::xfIntegratedTotal(double Q2) {

  // Freeze the scale below the initial scale.
  if(Q2 < Q02) Q2 = Q02;

  // Set the reference scales and relative contributions.
  // Gluons and u/d quarks has some non-perturbative contribution, others
  // only radiative contributions. Derived by fitting by eye to
  // a0 + a1*log(Q^2/Q0^2).
  double fq0[6] = { 0.0018, 0.0006, 0.0006, 0., 0., 0. };
  double mq2[6] = { Q02, Q02, Q02, Q02, pow2(MC), pow2(MB) };
  double eq2[6] = { 3.0/9.0, 1.0/9.0, 4.0/9.0, 1.0/9.0, 4.0/9.0, 1.0/9.0 };
  double a1 = 0.000981;

  // Logarithmic Q^2 evolution for each flavor. quarks two times, gluon
  // coefficents scaled appropriately.
  double xIntegrated = 0;
  for(int i = 0;i < 6;++i) {
    xIntegrated += fq0[i] + 2*a1*eq2[i]*max(0.0,log(Q2/mq2[i]));
  }

  return xIntegrated;

}

//--------------------------------------------------------------------------

// Returns the point-like part of the gluon.

double CJKL::pointlikeG(double x, double s) {

  // Exponents.
  double alpha1 = -0.43865;
  double alpha2 =  2.7174;
  double beta   =  0.36752;

  // Scale dependent parameters.
  double a  =  0.086893 - 0.34992  * s;
  double b  =  0.010556 + 0.049525 * s;
  double c  = -0.099005 + 0.34830  * s;
  double d  =  1.0648   + 0.143421 * s;
  double e  =  3.6717   + 2.5071   * s;
  double f  =  2.1944   + 1.9358   * s;
  double aa =  0.23679  - 0.11849  * s;
  double bb = -0.19994  + 0.028124 * s;

  // Point-like gluon parametrization.
  return max(0.0,( pow(s,alpha1)*pow(x,aa)*( a + b*sqrt(x) + c*pow(x,bb) )
    + pow(s,alpha2)*exp( -e + sqrt( f*pow(s,beta)*log(1.0/x) ) ) )
    * pow(1-x,d) );
}

//--------------------------------------------------------------------------

// Returns the point-like part of the u-quark.

double CJKL::pointlikeU(double x, double s) {

  // Exponents.
  double alpha1 = -1.0711;
  double alpha2 =  3.1320;
  double beta   =  0.69243;

  // Scale dependent parameters.
  double a  = -0.058266  + 0.20506  * s;
  double b  =  0.0097377 - 0.10617  * s;
  double c  = -0.0068345 + 0.15211  * s;
  double d  =  0.22297   + 0.013567 * s;
  double e  =  6.4289    + 2.2802   * s;
  double f  =  1.7302    + 0.76997  * s;
  double aa =  0.87940   - 0.110241 * s;
  double bb =  2.6878    - 0.040252 * s;

  // Point-like u-quark parametrization.
  return max(0.0, ( pow(s,alpha1)*pow(x,aa)*( a + b*sqrt(x) + c*pow(x,bb) )
    + pow(s,alpha2)*exp( -e + sqrt( f*pow(s,beta)*log(1.0/x) ) ) )
    * pow(1-x,d) );
}

//--------------------------------------------------------------------------

// Returns the point-like part of the d-quark.

double CJKL::pointlikeD(double x, double s) {

  // Exponents.
  double alpha1 = -1.1357;
  double alpha2 =  3.1187;
  double beta   =  0.66290;

  // Scale dependent parameters.
  double a  =  0.098814  - 0.067300  * s;
  double b  = -0.092892  + 0.049949  * s;
  double c  = -0.0066140 + 0.020427  * s;
  double d  = -0.31385   - 0.0037558 * s;
  double e  =  6.4671    + 2.2834    * s;
  double f  =  1.6996    + 0.84262   * s;
  double aa =  11.777    + 0.034760  * s;
  double bb = -11.124    - 0.20135   * s;

  // Regulate the x->1 divergence of (1-x)^d in the parameterization.
  if(x > 0.995) x = 0.995;

  // Point-like d-quark parametrization.
  return max( 0.0, ( pow(s,alpha1)*pow(x,aa)*( a + b*sqrt(x) + c*pow(x,bb) )
    + pow(s,alpha2)*exp( -e + sqrt( f*pow(s,beta)*log(1.0/x) ) ) )
    * pow(1-x,d) );
}

//--------------------------------------------------------------------------

// Returns the point-like part of the c-quark.

double CJKL::pointlikeC(double x, double s, double Q2) {

  // Scaled variable for c quarks with m = 1.3 GeV.
  double y = x + 1 - Q2/(Q2 + 6.76);

  // Kinematic boundary.
  if (y >= 1.0) return 0;

  // Declaration of parameters.
  double alpha1, alpha2, beta, a, b, c, d, e, f, aa, bb;

  // Parameters for Q^2 <= 10 GeV^2.
  if (Q2 <= 10) {

    // Exponents.
    alpha1 = 2.9808;
    alpha2 = 28.682;
    beta   = 2.4863;

    // Scale dependent parameters.
    a  = -0.18826   + 0.13565  * s;
    b  =  0.18508   - 0.11764  * s;
    c  = -0.0014153 - 0.011510 * s;
    d  = -0.48961   + 0.18810  * s;
    e  =  0.20911   - 2.8544   * s + 14.256 *s*s;
    f  =  2.7644    + 0.93717  * s;
    aa = -7.6307    + 5.6807   * s;
    bb =  394.58    - 541.82   * s + 200.82 *s*s;

  // Parameters for Q^2 > 10 GeV^2.
  } else {

    // Exponents.
    alpha1 = -1.8095;
    alpha2 =  7.9399;
    beta   =  0.041563;

    // Scale dependent parameters.
    a  = -0.54831  + 0.33412  * s;
    b  =  0.19484  + 0.041562 * s;
    c  = -0.39046  + 0.37194  * s;
    d  =  0.12717  + 0.059280 * s;
    e  =  8.7191   + 3.0194   * s;
    f  =  4.2616   + 0.73993  * s;
    aa = -0.30307  + 0.29430  * s;
    bb =  7.2383   - 1.5995   * s;
  }

  // Point-like c-quark parametrization.
  return max( 0.0, ( pow(s,alpha1)*pow(y,aa)*( a + b*sqrt(y) + c*pow(y,bb) )
    + pow(s,alpha2)*exp( -e + sqrt( f*pow(s,beta)*log(1.0/x) ) ) )
    * pow(1-y,d) );
}

//--------------------------------------------------------------------------

// Returns the point-like part of the b-quark.

double CJKL::pointlikeB(double x, double s, double Q2) {

  //Scaled variable for b quarks with m = 4.3 GeV.
  double y = x + 1 - Q2/(Q2 + 73.96);

  // Kinematic boundary.
  if (y >= 1.0) return 0;

  // Declaration of parameters.
  double alpha1, alpha2, beta, a, b, c, d, e, f, aa, bb;

  // Parameters for Q^2 <= 100 GeV^2.
  if (Q2 <= 100) {

    // Exponents.
    alpha1 =  2.2849;
    alpha2 =  6.0408;
    beta   = -0.11577;

    // Scale dependent parameters.
    a  = -0.26971   + 0.17942   * s;
    b  =  0.27033   - 0.18358   * s + 0.0061059 *s*s;
    c  =  0.0022862 - 0.0016837 * s;
    d  =  0.30807   - 0.10490   * s;
    e  =  14.812    - 1.2977    * s;
    f  =  1.7148    + 2.3532    * s + 0.053734  *sqrt(s);
    aa =  3.8140    - 1.0514    * s;
    bb =  2.2292    + 20.194    * s;

  // Parameters for Q^2 > 100 GeV^2.
  } else {

    // Exponents.
    alpha1 = -5.0607;
    alpha2 =  16.590;
    beta   =  0.87190;

    // Scale dependent parameters.
    a  = -0.72790  + 0.36549  * s;
    b  = -0.62903  + 0.56817  * s;
    c  = -2.4467   + 1.6783   * s;
    d  =  0.56575  - 0.19120  * s;
    e  =  1.4687   + 9.6071   * s;
    f  =  1.1706   + 0.99674  * s;
    aa = -0.084651 - 0.083206 * s;
    bb =  9.6036   - 3.4864   * s;
  }

  // Point-like b-quark parametrization.
  return max( 0.0, ( pow(s,alpha1)*pow(y,aa)*( a + b*sqrt(y) + c*pow(y,bb) )
    + pow(s,alpha2)*exp( -e + sqrt( f*pow(s,beta)*log(1.0/x) ) ) )
    * pow(1-y,d) );
}

//--------------------------------------------------------------------------

// Returns the hadron-like part of the gluon pdf.

double CJKL::hadronlikeG(double x, double s) {

  // Exponents.
  double alpha = 0.59945;
  double beta  = 1.1285;

  // Scale dependent parameters.
  double a  = -0.19898 + 0.57414 * s;
  double b  =  1.9942  - 1.8306  * s;
  double c  = -1.9848  + 1.4136  * s;
  double d  =  0.21294 + 2.7450  * s;
  double e  =  1.2287  + 2.4447  * s;
  double f  =  4.9230  + 0.18526 * s;
  double aa = -0.34948 + 0.47058 * s;

  // Hadron-like gluon parametrization.
  return max( 0.0, pow(1-x,d)*( pow(x,aa)*( a + b*sqrt(x) + c*x )
    + pow(s,alpha)*exp( -e + sqrt( f*pow(s,beta)*log(1.0/x) ) ) ) );
}

//--------------------------------------------------------------------------

// Returns the hadron-like part of the valence quarks.

double CJKL::hadronlikeVal(double x, double s) {

  // Scale dependent parameters.
  double a  =  1.0898  + 0.38087  * s;
  double b  =  0.42654 - 1.2128   * s;
  double c  = -1.6576  + 1.7075   * s;
  double d  =  0.96155 + 1.8441   * s;
  double aa =  0.78391 - 0.068720 * s;

  // Hadron-like valence quarks parametrization.
  return max( 0.0, pow(1-x,d)*pow(x,aa)*a*( 1 + b*sqrt(x) + c*x ) );
}

//--------------------------------------------------------------------------

// Returns the hadron-like part of the sea quarks.

double CJKL::hadronlikeSea(double x, double s) {

  // Exponents.
  double alpha = 0.71660;
  double beta  = 1.0497;

  // Scale dependent parameters.
  double a  =  0.60478 + 0.036160 * s;
  double b  =  4.2106  - 0.85835  * s;
  double d  =  4.1494  + 0.34866  * s;
  double e  =  4.5179  + 1.9219   * s;
  double f  =  5.2812  - 0.15200  * s;
  double aa =  0.72289 - 0.21562  * s;

  // Pre-calculate the logarithm.
  double logx = log(1.0/x);

  // Hadron-like sea quark parametrization.
  return max( 0.0, pow(1-x,d)*pow(s,alpha)*( 1 + a*sqrt(x) + b*x )
    * exp( -e + sqrt( f*pow(s,beta)*logx ) )*pow(logx,-aa) );
}

//--------------------------------------------------------------------------

// Returns the hadron-like part of the c-quarks.

double CJKL::hadronlikeC(double x, double s, double Q2) {

  //Scaled variable for c quarks with m = 1.3 GeV.
  double y = x + 1 - Q2/(Q2 + 6.76);

  // Kinematic boundary.
  if (y >= 1.0) return 0;

  // Pre-calculate the logarithm.
  double logx = log(1.0/x);

  // Declaration of parameters.
  double alpha, beta, a, b, d, e, f, aa;

  // Parameters for Q^2 <= 10 GeV^2.
  if (Q2 <= 10) {

    // Exponents.
    alpha = 5.6729;
    beta  = 1.4575;

    // Scale dependent parameters.
    a  = -2586.4 + 1910.1  * s;
    b  =  2695.0 - 1688.2  * s;
    d  =  1.5146 + 3.1028  * s;
    e  = -3.9185 + 11.738  * s;
    f  =  3.6126 - 1.0291  * s;
    aa =  1.6248 - 0.70433 * s;

  // Parameters for Q^2 > 10 GeV^2.
  } else {

    // Exponents.
    alpha = -1.6470;
    beta  =  0.72738;

    // Scale dependent parameters.
    a  = -2.0561  + 0.75576 * s;
    b  =  2.1266  + 0.66383 * s;
    d  =  3.0301  - 1.7499  * s + 1.6466  *s*s;
    e  =  4.1282  + 1.6929  * s - 0.26292 *s*s;
    f  =  0.89599 + 1.2761  * s - 0.15061 *s*s;
    aa = -0.78809 + 0.90278 * s;

  }

  // Hadron-like c-quark parametrization. Note typo in the CJKL paper.
  return max( 0.0, pow(1-y,d)*pow(s,alpha)*( 1 + a*sqrt(y) + b*y )
    * exp( -e + f*sqrt( pow(s,beta)*logx ) )*pow(logx,-aa) );
}

//--------------------------------------------------------------------------

// Returns the hadron-like part of the b-quarks.

double CJKL::hadronlikeB(double x, double s, double Q2) {

  // Scaled variable for b quarks with m = 4.3 GeV.
  double y = x + 1 - Q2/(Q2 + 73.96);

  // Kinematic boundary.
  if (y >= 1.0) return 0;

  // Pre-calculate the logarithm.
  double logx = log(1.0/x);

  // Declaration of parameters.
  double alpha, beta, a, b, d, e, f, aa;

  // Parameters for Q^2 <= 100 GeV^2.
  if (Q2 <= 100) {

    // Exponents.
    alpha = -10.210;
    beta  = -2.2296;

    // Scale dependent parameters.
    a  = -99.613  + 171.25   * s;
    b  =  492.61  - 420.45   * s;
    d  =  3.3917  + 0.084256 * s;
    e  =  5.6829  - 0.23571  * s;
    f  = -2.0137  + 4.6955   * s;
    aa =  0.82278 + 0.081818 * s;

  // Parameters for Q^2 > 100 GeV^2.
  } else {

    // Exponents.
    alpha = 2.4198;
    beta  = 0.40703;

    // Scale dependent parameters.
    a  = -2.1109  + 1.2711  * s;
    b  =  9.0196  - 3.6082  * s;
    d  =  3.6455  - 4.1353  * s + 2.3615  *s*s;
    e  =  4.6196  + 2.4212  * s;
    f  =  0.66454 + 1.1109  * s;
    aa = -0.98933 + 0.42366 * s + 0.15817 *s*s;
  }

  // Hadron-like b-quark parametrization. Note typo in the CJKL paper.
  return max( 0.0, pow(1-y,d)*pow(s,alpha)*( 1 + a*sqrt(y) + b*y )
    * exp( -e + f*sqrt( pow(s,beta)*logx ) )*pow(logx,-aa) );
}

//==========================================================================

// The LHAGrid1 class.
// Codes to read files in the LHAPDF6 lhagrid1 format,
// assuming that the same x grid is used for all Q subgrids.
// Results are not identical with LHAPDF6, owing to different interpolation.

//--------------------------------------------------------------------------

// Initialize PDF: select data file and open stream.

void LHAGrid1::init(string pdfWord, string xmlPath, Info* infoPtr) {

  // Identify whether file number or name.
  if (pdfWord.length() > 9 && toLower(pdfWord).substr(0,9) == "lhagrid1:")
    pdfWord = pdfWord.substr(9, pdfWord.length() - 9);
  istringstream pdfStream(pdfWord);
  int pdfSet = 0;
  pdfStream >> pdfSet;

  // Input is file name.
  string dataFile = "";
  if ( xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";
  if (pdfWord[0] == '/') dataFile = pdfWord;
  else if (pdfSet == 0) dataFile = xmlPath + pdfWord;

  // Input is fit number. Current selection for NNPDF3.1 and modified NNLO.
  else if (pdfSet == 17) dataFile = xmlPath + "NNPDF31_lo_as_0130_0000.dat";
  else if (pdfSet == 18) dataFile = xmlPath + "NNPDF31_lo_as_0118_0000.dat";
  else if (pdfSet == 19) dataFile = xmlPath
    + "NNPDF31_nlo_as_0118_luxqed_0000.dat";
  else if (pdfSet == 20) dataFile = xmlPath
    + "NNPDF31_nnlo_as_0118_luxqed_0000.dat";
  else if (pdfSet == 21) dataFile = xmlPath
    + "NNPDF31sx_nlonllx_as_0118_LHCb_luxqed_0000.dat";
  else if (pdfSet == 22) dataFile = xmlPath
    + "NNPDF31sx_nnlonllx_as_0118_LHCb_luxqed_0000.dat";

  // Pomeron PDFs, currently the GKG18 sets.
  else if (pdfSet == 112) dataFile = xmlPath + "GKG18_DPDF_FitA_LO_0000.dat";
  else if (pdfSet == 113) dataFile = xmlPath + "GKG18_DPDF_FitB_LO_0000.dat";
  else if (pdfSet == 114) dataFile = xmlPath + "GKG18_DPDF_FitA_NLO_0000.dat";
  else if (pdfSet == 115) dataFile = xmlPath + "GKG18_DPDF_FitB_NLO_0000.dat";

  // Open files from which grids should be read in.
  ifstream is( dataFile.c_str() );
  if (!is.good()) {
    printErr("Error in LHAGrid1::init: did not find data file", infoPtr);
    isSet = false;
    return;
  }

  // Initialization with a stream.
  init( is, infoPtr);
  is.close();

}

//--------------------------------------------------------------------------

// Initialize PDF: read in data grid from stream and set up interpolation.

void LHAGrid1::init(istream& is, Info* infoPtr) {

  // Check that data stream is available.
  if (!is.good()) {
    printErr("Error in LHAGrid1::init: cannot read from stream", infoPtr);
    isSet = false;
    return;
  }

  // Some local variables.
  string line;
  vector<string> idlines, pdflines;
  int nqNow, idNow, idNowMap;
  double xNow, qNow, pdfNow;

  // Skip lines of header, until ---. Probe for next subgrid in Q space.
  nqSub = 0;
  do getline( is, line);
  while (line.find("---") == string::npos);
  if (!is.good()) {
    printErr("Error in LHAGrid1::init: could not read data file", infoPtr);
    isSet = false;
    return;
  }
  while (getline( is, line)) {
    ++nqSub;

    // Read in x grid; save for first, check it matches for later ones.
    istringstream isx(line);
    if (nqSub == 1) {
      while (isx >> xNow) {
        xGrid.push_back( xNow);
        lnxGrid.push_back( log(xNow));
      }
      nx   = xGrid.size();
      xMin = xGrid.front();
      xMax = xGrid.back();
    } else {
      int ixc = -1;
      while (isx >> xNow)
      if ( abs(log(xNow) - lnxGrid[++ixc]) > 1e-5) {
        printErr("Error in LHAGrid1::init: mismatched subgrid x spacing",
          infoPtr);
        isSet = false;
        return;
      }
    }

    // Read in Q grid; append as needed. Check that subgrids match.
    getline( is, line);
    istringstream isq(line);
    nqNow = 0;
    while (isq >> qNow) {
      ++nqNow;
      qGrid.push_back( qNow);
      lnqGrid.push_back( log(qNow));
    }
    if (nqSub > 1) {
      if (abs(qGrid[nq] / qGrid[nq-1] - 1.) > 1e-5) {
        printErr("Error in LHAGrid1::init: mismatched subgrid Q borders",
          infoPtr);
        isSet = false;
        return;
      }
      qGrid[nq-1] = 0.5 * (qGrid[nq-1] + qGrid[nq]);
      qGrid[nq]   = qGrid[nq-1];
    }
    nq   = qGrid.size();
    qMin = qGrid.front();
    qMax = qGrid.back();
    nqSum.push_back(nq);
    qDiv.push_back(qMax);

    // Read in and store flavour mapping and pdf data. Separator line.
    getline( is, line);
    idlines.push_back( line);
    for (int ixq = 0; ixq < nx * nqNow; ++ixq) {
      getline( is, line);
      pdflines.push_back( line);
    }
    getline( is, line);
  }

  // Create array big enough to hold (flavour, x, Q) grid.
  pdfGrid = new double**[12];
  for (int iid = 0; iid < 12; ++iid) {
    pdfGrid[iid] = new double*[nx];
    for (int ix = 0; ix < nx; ++ix) {
      pdfGrid[iid][ix] = new double[nq];
      for (int iq = 0; iq < nq; ++iq) pdfGrid[iid][ix][iq] = 0.;
    }
  }

  // Second pass through the Q subranges.
  int iln = -1;
  for (int iqSub = 0; iqSub < nqSub; ++iqSub) {
    vector<int> idGridMap;

    // Study flavour grid and decide flavour mapping.
    istringstream isid( idlines[iqSub] );
    while (isid >> idNow) {
      idNowMap = -1;
      if (idNow == 21 || idNow == 0) idNowMap = 0;
      if (idNow > 0 && idNow < 6) idNowMap = idNow;
      if (idNow < 0 && idNow > -6) idNowMap = 5 - idNow;
      if (idNow == 22) idNowMap = 11;
      idGridMap.push_back( idNowMap);
    }
    int nid = idGridMap.size();

    // Read in data grid, line by line.
    int iq0 = (iqSub == 0) ? 0 : nqSum[iqSub - 1];
    for (int ix = 0; ix < nx; ++ix)
    for (int iq = iq0; iq < nqSum[iqSub]; ++iq) {
      istringstream ispdf( pdflines[++iln] );
      for (int iid = 0; iid < nid; ++iid) {
        ispdf >> pdfNow;
        if (idGridMap[iid] >= 0) pdfGrid[idGridMap[iid]][ix][iq] = pdfNow;
      }
    }
  }

  // For extrapolation to small x: create array for b values of x^b shape.
  pdfSlope = new double*[12];
  for (int iid = 0; iid < 12; ++iid) {
    pdfSlope[iid] = new double[nq];
    for (int iq = 0; iq < nq; ++iq) { pdfSlope[iid][iq] =
      ( min( pdfGrid[iid][0][iq], pdfGrid[iid][1][iq]) > 1e-5)
      ? ( log(pdfGrid[iid][1][iq]) - log(pdfGrid[iid][0][iq]) )
      / (lnxGrid[1] - lnxGrid[0]) : 0.;
    }
  }

}

//--------------------------------------------------------------------------

void LHAGrid1::xfUpdate(int , double x, double Q2) {

  // No PDF values if not properly set up.
  if (!isSet) {
    xg = xu = xd = xubar = xdbar = xs = xsbar = xc = xb = xgamma
    = xuVal = xuSea = xdVal = xdSea = 0.;
    return;
  }

  // Update within allowed (x, q) range.
  xfxevolve( x, Q2);

  // Then transfer to Pythia8 notation.
  xg     = pdfVal[0];
  xu     = pdfVal[2];
  xd     = pdfVal[1];
  xubar  = pdfVal[7];
  xdbar  = pdfVal[6];
  xs     = pdfVal[3];
  xsbar  = pdfVal[8];
  xc     = 0.5 * (pdfVal[4] + pdfVal[9]);
  xb     = 0.5 * (pdfVal[5] + pdfVal[10]);
  xgamma = pdfVal[11];

  // Subdivision of valence and sea.
  xuVal  = xu - xubar;
  xuSea  = xubar;
  xdVal  = xd - xdbar;
  xdSea  = xdbar;

  // idSav = 9 to indicate that all flavours reset.
  idSav  = 9;

}

//--------------------------------------------------------------------------

void LHAGrid1::xfxevolve(double x, double Q2) {

  // Find if (x, Q) inside our outside grid.
  double q = sqrt(Q2);
  int inx  = (x <= xMin) ? -1 : ((x >= xMax) ? 1 : 0);
  int inq  = (q <= qMin) ? -1 : ((q >= qMax) ? 1 : 0);

  // Set up default for x interpolation.
  int    minx  = 0;
  int    maxx  = nx - 1;
  int    m3x   = 0;
  double wx[4] = {1., 1., 1., 1.};

  // Find grid value on either side of x.
  if (inx == 0) {
    int midx;
    while (maxx - minx > 1) {
      midx = (minx + maxx) / 2;
      if (x < xGrid[midx]) maxx = midx;
      else                 minx = midx;
    }

    // Weights for cubic interpolation in ln(x).
    double lnx = log(x);
    if      (minx == 0)      m3x = 0;
    else if (maxx == nx - 1) m3x = nx - 4;
    else                     m3x = minx - 1;
    for (int i3 = 0; i3 < 4; ++i3)
    for (int j = 0; j < 4; ++j) if (j != i3)
      wx[i3] *= (lnx - lnxGrid[m3x+j]) / (lnxGrid[m3x+i3] - lnxGrid[m3x+j]);
  }

  // Find q subgrid and set up default for q interpolation.
  int    iqDiv = 0;
  for (int iqSub = 1; iqSub < nqSub; ++iqSub)
    if (q > qDiv[iqSub - 1]) iqDiv = iqSub;
  int    minS  = (iqDiv == 0) ? 0 : nqSum[iqDiv - 1];
  int    maxS  = nqSum[iqDiv] - 1;
  int    minq  = minS;
  int    maxq  = maxS;
  int    n3q   = 4;
  int    m3q   = 0.;
  double wq[4] = {1., 1., 1., 1.};

  // Find grid value on either side of q.
  if (inq == 0) {
    int midq;
    while (maxq - minq > 1) {
      midq = (minq + maxq) / 2;
      if (q < qGrid[midq]) maxq = midq;
      else                 minq = midq;
    }

    // Weights for linear or cubic interpolation in ln(q).
    double lnq = log(q);
    if (maxS - minS < 3) {
      n3q = 2;
      m3q = minq;
      wq[1] = (lnq - lnqGrid[minq]) / (lnqGrid[maxq] - lnqGrid[minq]);
      wq[0] = 1. - wq[1];
    } else {
      if      (minq == minS) m3q = minS;
      else if (maxq == maxS) m3q = maxS - 3;
      else                   m3q = minq - 1;
      for (int i3 = 0; i3 < 4; ++i3)
      for (int j = 0; j < 4; ++j) if (j != i3)
        wq[i3] *= (lnq - lnqGrid[m3q+j]) / (lnqGrid[m3q+i3] - lnqGrid[m3q+j]);
    }

  // Freeze at border of q range.
  } else {
    n3q = 1;
    if (inq == 1) m3q = nq - 1;
  }

  // Interpolate between grid elements, normally bicubic, or simpler in ln(q).
  for (int iid = 0; iid < 12; ++iid) pdfVal[iid] = 0.;
  if (inx == 0) {
    for (int iid = 0; iid < 12; ++iid)
    for (int i3x = 0; i3x < 4; ++i3x)
    for (int i3q = 0; i3q < n3q; ++i3q)
      pdfVal[iid] += wx[i3x] * wq[i3q] * pdfGrid[iid][m3x+i3x][m3q+i3q];

  // Special: extrapolate to small x. (Let vanish at large x, so no such code.)
  } else if (inx == -1) {
    for (int iid = 0; iid < 12; ++iid)
    for (int i3q = 0; i3q < n3q; ++i3q)
      pdfVal[iid] += wq[i3q] * pdfGrid[iid][0][m3q+i3q]
        * (doExtraPol ? pow( x / xMin, pdfSlope[iid][m3q+i3q]) : 1.);
  }

}

//==========================================================================

// Convolution with photon flux from leptons and photon PDFs.
// Contains a pointer to a photon PDF set and samples the
// convolution integral event-by-event basis.
// Includes also a overestimate for the PDF set in order to set up
// the phase-space sampling correctly.

// Constants related to the fit.
const double Lepton2gamma::ALPHAEM = 0.007297353080;
const double Lepton2gamma::Q2MIN   = 1.;

//--------------------------------------------------------------------------

// Update PDFs and sample a value for x_gamma.

void Lepton2gamma::xfUpdate(int , double x, double Q2) {

  // Find the maximum x value at given Q2max and sqrt(s).
  double sCM = infoPtr->s();
  double xGamMax = ( 2. - 2. * Q2max / sCM - 8. * m2lepton / sCM )
    / ( 1. + sqrt( (1. + 4. * m2lepton / Q2max) * (1. - 4. * m2lepton/sCM) ) );

  // If outside allowed x values set PDFs to zero.
  if ( x > xGamMax ) {
    xg     = 0.;
    xd     = 0.;
    xu     = 0.;
    xs     = 0.;
    xc     = 0.;
    xb     = 0.;
    xubar  = 0.;
    xdbar  = 0.;
    xsbar  = 0.;
    xGm    = 1.;
    return;
  }

  // Pre-calculate some logs.
  double log2x    = pow2( log( Q2max / (m2lepton * pow2(x)) ) );
  double log2xMax = pow2( log( Q2max / (m2lepton * pow2(xGamMax)) ) );

  // Sample x_gamma.
  if ( sampleXgamma) {
    xGm = sqrt( (Q2max / m2lepton)
      * exp( -sqrt( log2x + rndmPtr->flat() * (log2xMax - log2x) ) ) );
  }

  // Evaluate the PDFs at x/x_gamma.
  double xInGamma = x/xGm;
  double xgGm = gammaPDFPtr->xf(21, xInGamma, Q2);
  double xdGm = gammaPDFPtr->xf(1 , xInGamma, Q2);
  double xuGm = gammaPDFPtr->xf(2 , xInGamma, Q2);
  double xsGm = gammaPDFPtr->xf(3 , xInGamma, Q2);
  double xcGm = gammaPDFPtr->xf(4 , xInGamma, Q2);
  double xbGm = gammaPDFPtr->xf(5 , xInGamma, Q2);

  // Calculate the Q^2_min for sampled x_gamma.
  double m2s   = 4. * m2lepton / sCM;
  double Q2min = 2. * m2lepton * pow2(xGm)
    / ( 1. - xGm - m2s + sqrt(1. - m2s) * sqrt( pow2(1. - xGm) - m2s ) );

  // Correct with weight.
  double alphaLog = (ALPHAEM / (2. * M_PI)) * (1. + pow2(1. - xGm) )
    * 0.25 * (log2x - log2xMax) * log(Q2max / Q2min)
    / log( Q2max / ( m2lepton * pow2(xGm) ) );

  // Calculate the PDF value.
  xg = alphaLog * xgGm;
  xd = alphaLog * xdGm;
  xu = alphaLog * xuGm;
  xs = alphaLog * xsGm;
  xc = alphaLog * xcGm;
  xb = alphaLog * xbGm;
  xubar  = xu;
  xdbar  = xd;
  xsbar  = xs;

  // Photon inside electron not currently implemented (Use point-like lepton).
  xgamma = 0;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//--------------------------------------------------------------------------

// Approximate the maximum of convoluted PDF to correctly set up the
// sampling of the phase space.

double Lepton2gamma::xfMax(int id, double x, double Q2) {

  // Find the maximum x value at given Q2max and sqrt(s).
  double sCM = infoPtr->s();
  double xGamMax = ( 2. - 2. * Q2max / sCM - 8. * m2lepton / sCM )
    / ( 1. + sqrt( (1. + 4. * m2lepton / Q2max) * (1. - 4. * m2lepton/sCM) ) );

  // Set PDFs to zero outside allowed x values.
  if ( x > xGamMax ) return 0;

  // Pre-calculate some logs.
  double log2x    = pow2( log( Q2max / (m2lepton * pow2(x)) ) );
  double log2xMax = pow2( log( Q2max / (m2lepton * pow2(xGamMax)) ) );

  // Find approximate x-behaviour for each flavour. Optimized for CJKL.
  double xApprox = 0.;
  int idAbs = abs(id);
  if      (idAbs == 21 || idAbs == 0) xApprox = 2.35;
  else if (idAbs == 1) xApprox = (pow(x, 0.2) + pow(1. - x, -0.15)) * 0.8;
  else if (idAbs == 2) xApprox = (pow(x, 1.0) + pow(1. - x, -0.4))  * 0.4;
  else if (idAbs == 3) xApprox = (pow(x, 0.2) + pow(1. - x, -0.5))  * 0.5;
  else if (idAbs == 4) xApprox = (pow(x, 1.0) + pow(1. - x, -0.4))  * 0.7;
  else if (idAbs == 5) xApprox = (pow(x, 0.2) + pow(1.  -x, -0.5))  * 0.5;
  else xApprox = 0.;

  // Direct photons in usual lepton PDFs.
  if ( idAbs == 22 ) return 0;

  // Return the approximation.
  return (ALPHAEM / (2. * M_PI)) * (log2x - log2xMax) * 0.5
    * gammaPDFPtr->xf(id, x, Q2) / xApprox;
}

//--------------------------------------------------------------------------

// Return PDF without sampling x_gamma values to compute cross section with
// rescaled sHat. Not very elegant but no need to modify the xfUpdate call.

double Lepton2gamma::xfSame(int id, double x, double Q2) {
  sampleXgamma = false;
  xfUpdate(id, x, Q2);
  double xfNow = xf(id, x, Q2);
  sampleXgamma = true;
  return xfNow;
}

//==========================================================================

// Approximated photon flux that used for process sampling with external flux.

const double EPAexternal::ALPHAEM = 0.007297353080;

// Initilaize kinematics and find the normalization.

void EPAexternal::init() {

  // Collision kinematics.
  double sCM = pow2(infoPtr->eCM());
  double m2s = 4. * m2 / sCM;

  // Photon kinematics.
  xMin  = pow2(settingsPtr->parm("Photon:Wmin")) / sCM;
  xMax  = 1.0;

  // Select which overestimate is used for sampling.
  approxMode = settingsPtr->mode("PDF:lepton2gammaApprox");

  // Approximation suited for lepton beams.
  if (approxMode == 1) {

    // Derive kinematics.
    Q2min = 2. * m2 * pow2(xMin) / ( 1. - xMin - m2s
          + sqrt(1. - m2s) * sqrt( pow2(1. - xMin) - m2s) );
    Q2max = settingsPtr->parm("Photon:Q2max");
    xMax  = 2. * ( 1. - Q2max / sCM - m2s )
          / ( 1. + sqrt( (1. + 4. * m2 / Q2max) * (1. - m2s) ) );
    bool sampleQ2 = settingsPtr->flag("Photon:sampleQ2");

    // Initial values for normalization.
    double ratio, ratioMax = 0.0;
    norm = 1.0;

    // Scan through x and Q2 grid to find normalization.
    // Mainly required for flux from heavy ions with large charge.
    for (int i = 0; i < 10; ++i) {
      double xi = xMin + (xMax - xMin)*i/(10.);

      // If sampling for Q2 also, scan through the Q2 grid as well.
      if (sampleQ2) {
        for (int j = 0; j < 10; ++j) {
          double Q2j = Q2min * exp( log(Q2max/Q2min)*j/(10. - 1.0));
          ratio = xfFlux(22,xi,Q2j) / xfApprox(22,xi,Q2j);
          if (ratio > ratioMax) ratioMax = ratio;
        }

      // If not, scanning x-grid suffice.
      } else {
        ratio = xfFlux(22,xi) / xf(22,xi,1.);
        if (ratio > ratioMax) ratioMax = ratio;
      }
    }

    // Store the found normalization.
    norm = ratioMax;

  // Sampling optimized for heavy-ions with flux proportional to modified
  // bessel functions. Divided into regions with x^pow and exp(-A*x).
  } else if (approxMode == 2) {

    // Find the parameters for the overestimate and derive further variables.
    double mBeam = settingsPtr->parm("PDF:gammaFluxApprox2bMin");
    double bMin  = settingsPtr->parm("PDF:gammaFluxApprox2mBeam");
    xPow         = settingsPtr->parm("PDF:gammaFluxApprox2xPow");
    xCut         = settingsPtr->parm("PDF:gammaFluxApprox2xCut");
    bmhbarc      = bMin * mBeam / HBARC;

    // Normalizations for the two regions from the flux.
    norm1 = xMin < xCut ? pow(xMin, -1. + xPow) * xfFlux(22,xMin) : 0.0;
    norm2 = xMin < xCut ? exp( 2. * bmhbarc * xCut) * xfFlux(22,xCut) / xCut
                        : exp( 2. * bmhbarc * xMin) * xfFlux(22,xMin) / xMin;

    // Integrals of the two regions for cross section approximation.
    integral1 = xMin < xCut ? norm1 / (1. - xPow)
              * ( pow(xCut, 1. - xPow) - pow(xMin, 1. - xPow) ) : 0.;
    integral2 = xMin < xCut ? norm2 * 0.5 / bmhbarc
              * ( exp(-2. * bmhbarc * xCut) - exp(-2. * bmhbarc) )
              : norm2 * 0.5 / bmhbarc
              * ( exp(-2. * bmhbarc * xMin) - exp(-2. * bmhbarc) );
  }

}

//--------------------------------------------------------------------------

// Approximate the differential photon flux with alphaEM/PI/x/Q2.
// Derived from EPA for leptons but provides leading (small-x)
// behaviour for hadrons as well.

void EPAexternal::xfUpdate(int , double x, double Q2) {

  // Calculate (Q2-integrated) approximation for xfGamma.
  double alphaLog = (approxMode == 1) ?
    norm * ALPHAEM / M_PI * log (Q2max/Q2min) : 1.;

  // Integrated in Q2, to be used for direct process sampling.
  if (approxMode == 1) {
    xgamma = alphaLog;
  } else if (approxMode == 2) {
    if (x < xCut) xgamma = norm1 * pow(x, 1. - xPow);
    else          xgamma = norm2 * x * exp(-2. * bmhbarc * x);
  }

  // Approximate the convolution with photon PDFs.
  if (gammaPDFPtr != 0) {

    // To preserve x/xGamma < 1.
    xHadr            = x;
    double alphaLogX = 0.;

    // Integrals for the overestimates.
    if (approxMode == 1) {
      alphaLogX = alphaLog * log (xMax / xHadr);
    } else if (approxMode == 2) {
      double integral1tmp = xHadr < xCut ? norm1 / (1. - xPow)
        * ( pow(xCut, 1. - xPow) - pow(xHadr, 1. - xPow) ) : 0.;
      double xMinTmp = xHadr < xCut ? xCut : xHadr;
      double integral2tmp = norm2 * 0.5 / bmhbarc
        * ( exp(-2. * bmhbarc * xMinTmp) - exp(-2. * bmhbarc) );
      alphaLogX = integral1tmp + integral2tmp;
    }

    // Multiply the approximated flux with PDFs.
    xg = alphaLogX * gammaPDFPtr->xf(21, x, Q2);
    xd = alphaLogX * gammaPDFPtr->xf( 1, x, Q2);
    xu = alphaLogX * gammaPDFPtr->xf( 2, x, Q2);
    xs = alphaLogX * gammaPDFPtr->xf( 3, x, Q2);
    xc = alphaLogX * gammaPDFPtr->xf( 4, x, Q2);
    xb = alphaLogX * gammaPDFPtr->xf( 5, x, Q2);
    xdbar = xd;
    xubar = xu;
    xsbar = xs;
  }

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//--------------------------------------------------------------------------

// The approximated photon flux x*f^{gamma}(x,Q2).

double EPAexternal::xfApprox(int , double x, double Q2) {

  // Differential in Q2 for leptons.
  if (approxMode == 1) {
    return norm * ALPHAEM / M_PI / Q2;

  // Piece-wise approximation for heavy ions.
  } else if (approxMode == 2) {
    if (x < xCut) return norm1 * pow(x, 1. - xPow);
    else          return norm2 * x * exp(-2. * bmhbarc * x);
  }

  return 0.;
}

//--------------------------------------------------------------------------

// Accurate flux, provided externally.

double EPAexternal::xfFlux(int id, double x, double Q2) {

  // The external flux, check that pointer exists.
  if ( gammaFluxPtr != 0 ) return gammaFluxPtr->xf(id, x, Q2);
  else return 0.;
}

//--------------------------------------------------------------------------

// Photon PDFs used for the convolution with the flux.

double EPAexternal::xfGamma(int id, double x, double Q2) {

  // Return xf from the photon PDF.
  if ( gammaPDFPtr != 0 ) return gammaPDFPtr->xf(id, x, Q2);
  else return 0.;
}

//--------------------------------------------------------------------------

// Sample the x_gamma value according to given photon flux approximation.

double EPAexternal::sampleXgamma(double xMinIn) {

  // Sample with lepton-type flux.
  double xMinSample = (xMinIn < 0.) ? xMin : xMinIn;
  if (approxMode == 1) {
    return xMinSample * pow(xMax / xMinSample, rndmPtr->flat());

  // Sample with photon flux for nuclei.
  } else if (approxMode == 2) {

    // Calculate the integrals of over estimates.
    double integral1tmp = xMinSample < xCut ? norm1 / (1. - xPow)
      * ( pow(xCut, 1. - xPow) - pow(xMinSample, 1. - xPow) ) : 0.;
    double integral2tmp = norm2 * 0.5 / bmhbarc
      * ( exp(-2. * bmhbarc * xMinSample) - exp(-2. * bmhbarc) );
    double integral1Frac = integral1tmp / (integral1tmp + integral2tmp);

    // Select the sampling region.
    int samplingRegion = 1;
    if ( xMinSample > xCut || integral1Frac < rndmPtr->flat() )
      samplingRegion = 2;

    // Sample x.
    double xGm = (samplingRegion == 1)
      ? pow( pow(xMinSample,1. - xPow) + rndmPtr->flat()
      * ( pow(xCut, 1. - xPow) - pow(xMinSample, 1. - xPow) ), 1./(1. - xPow))
      : -0.5 / bmhbarc * log( exp(-2. * bmhbarc * xMinSample) - rndmPtr->flat()
      * ( exp(-2. * bmhbarc * xMinSample) -  exp(-2. * bmhbarc) ) );
    return xGm;
  }

  // Return zero for undefined cases.
  return 0.;

}

//--------------------------------------------------------------------------

// Return integrated over-estimate for photon flux to approximate soft
// cross sections.

double EPAexternal::intFluxApprox() {

  // Check the used approximation and return the integral.
  if ( approxMode == 1 )
    return ALPHAEM / M_PI * norm * log (xMax/xMin) * log(Q2max/Q2min);
  else if (approxMode == 2) return integral1 + integral2;
  else return 0.;

}

//==========================================================================

// Inherited class for nuclear PDFs. Needs a proton PDF as a baseline.

void nPDF::initNPDF(PDF* protonPDFPtrIn) {

  // Derive mass number and number of protons.
  a = (idBeam/10) % 1000;
  z = (idBeam/10000) % 1000;

  // Normalized number of protons and neutrons in nuclei.
  resetMode();

  // Initialize proton PDF pointer.
  protonPDFPtr = protonPDFPtrIn;

  // No modifications yet.
  ruv = 1.;
  rdv = 1.;
  ru  = 1.;
  rd  = 1.;
  rs  = 1.;
  rc  = 1.;
  rb  = 1.;
  rg  = 1.;
}

//--------------------------------------------------------------------------

// Updates the nPDF using provided proton PDF and nuclear modification.

void nPDF::xfUpdate(int id, double x, double Q2) {

  if (protonPDFPtr == 0) {
    printErr("Error in nPDF: No free proton PDF pointer set.");
    return;
  }

  // Update the proton PDFs and nuclear modifications.
  this->rUpdate(id, x, Q2);

  // u(bar) and d(bar) pdfs for proton.
  double xfd  = protonPDFPtr->xf( 1, x, Q2);
  double xfu  = protonPDFPtr->xf( 2, x, Q2);
  double xfdb = protonPDFPtr->xf(-1, x, Q2);
  double xfub = protonPDFPtr->xf(-2, x, Q2);

  // Neutron nPDFs using isospin symmetry.
  xd     = za * (rdv * (xfd - xfdb) + rd * xfdb)
         + na * (ruv * (xfu - xfub) + ru * xfub);
  xu     = za * (ruv * (xfu - xfub) + ru * xfub)
         + na * (rdv * (xfd - xfdb) + rd * xfdb);
  xdbar  = za * xfdb * rd + na * xfub * ru;
  xubar  = za * xfub * ru + na * xfdb * rd;
  xs     = rs * protonPDFPtr->xf( 3, x, Q2);
  xsbar  = rs * protonPDFPtr->xf(-3, x, Q2);
  xc     = rc * protonPDFPtr->xf( 4, x, Q2);
  xb     = rb * protonPDFPtr->xf( 5, x, Q2);
  xg     = rg * protonPDFPtr->xf(21, x, Q2);
  xgamma = 0.;

  // idSav = 9 to indicate that all flavours reset.
  idSav  = 9;

}

//==========================================================================

// Nuclear modifications of the PDFs from EPS09 fit, either LO or NLO.
// Ref: K.J. Eskola, H. Paukkunen and C.A. Salgado, JHEP 0904 (2009) 065.
// Grids files of different nuclei can be found from
// https://www.jyu.fi/science/en/physics/research/highenergy/urhic/npdfs/eps09

// Constants related to the fit.
const double EPS09::Q2MIN = 1.69;
const double EPS09::Q2MAX = 1000000.;
const double EPS09::XMIN  = 0.000001;
const double EPS09::XMAX  = 1.;
const double EPS09::XCUT  = 0.1;
const int EPS09::XSTEPS   = 50;
const int EPS09::Q2STEPS  = 50;

//--------------------------------------------------------------------------

// Initialize EPS09 nPDFs with given order (1=LO, 2=NLO) and error set.

void EPS09::init(int iOrderIn, int iSetIn, string xmlPath) {

  // Save the order and error set number.
  iOrder = iOrderIn;
  iSet   = iSetIn;

  // Select which data file to read for current fit.
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";
  stringstream fileSS;

  if (iOrder == 1) fileSS << xmlPath << "EPS09LOR_" << getA();
  if (iOrder == 2) fileSS << xmlPath << "EPS09NLOR_" << getA();
  string gridFile = fileSS.str();

  // Open grid file.
  ifstream fileStream( gridFile.c_str() );
  if (!fileStream.good()) {
    printErr("Error in EPS09::init: did not find grid file " + gridFile,
             infoPtr);
    isSet = false;
    return;
  }

  // Dump additional grid information here.
  double dummy;

  // Read in the interpolation grid.
  for (int i = 0;i < 31; ++i) {
    for (int j = 0;j < 51; ++j) {
      fileStream >> dummy;
      for (int k = 0;k < 51; ++k) {
        for (int l = 0;l < 8; ++l) fileStream >> grid[i][j][k][l];
      }
    }
  }
  fileStream.close();

}

//--------------------------------------------------------------------------

// Interpolation from the grid.

void EPS09::rUpdate(int , double x, double Q2) {

  // Freeze the x and Q2 values if outside the grid.
  if( x  < XMIN )  x  = XMIN;
  if( x  > XMAX )  x  = XMAX;
  if( Q2 < Q2MIN ) Q2 = Q2MIN;
  if( Q2 > Q2MAX ) Q2 = Q2MAX;

  // Calculate the position in log(log Q^2) grid:
  double dQ2 = Q2STEPS * log( log(Q2) / log(Q2MIN) )
    / log( log(Q2MAX) / log(Q2MIN) );
  int iQ2 = int(dQ2);

  // Set the Q2 index to interval [1,...,49].
  if      ( iQ2 < 1 )           iQ2 = 1;
  else if ( iQ2 > Q2STEPS - 1 ) iQ2 = Q2STEPS - 1;

  // Calculate the three nearest points in log(log Q^2) grid.
  double Q2Near[3];
  Q2Near[0] = iQ2 - 1;
  Q2Near[1] = iQ2 + 0;
  Q2Near[2] = iQ2 + 1;

  // Interpolate the grid values.
  for ( int iFlavour = 0; iFlavour < 8; ++iFlavour) {

    // Calculate the position in log(x) or x grid.
    int ix;
    int nxlog = XSTEPS/2;
    int nxlin = XSTEPS - nxlog;
    if ( x <= XCUT ) ix = int( nxlog * log(x / XMIN) / log( XCUT / XMIN ) );
    else ix = int( ( x - XCUT ) * nxlin / ( XMAX - XCUT ) + nxlog );

    // Set the x-index to interval [1,...,48].
    if ( ix < 1 ) ix = 1;

    // Do not use the last grid points for interpolation.
    if ( iFlavour == 0 || iFlavour == 1 || iFlavour == 7)
      if ( ix >= XSTEPS - 4 ) ix = XSTEPS - 4;
    if ( iFlavour > 1 && iFlavour < 7 )
      if ( ix >= XSTEPS - 7 ) ix = XSTEPS - 7;

    // Calculate the four nearest points in log-x or lin-x grid.
    double xNear[4];
    for(int i = 0;i < 4;i++) {
      if ( ix - 1 + i < nxlog ) {
        xNear[i] = XMIN * exp( ( double( ix - 1 + i ) / nxlog )
          * log( XCUT / XMIN ) );
      } else {
        xNear[i] = ( double( ix - 1 + i - nxlog) / nxlin )
          * ( XMAX - XCUT ) + XCUT;
      }
    }

    // Grid points used for interpolation.
    double xGrid[4];
    double Q2Grid[3];

    // Read in the relevant values from table and interpolate in x.
    for ( int j = 0; j < 3; ++j) {
      xGrid[0]  = grid[iSet - 1][iQ2 - 1 + j][ix - 1][iFlavour];
      xGrid[1]  = grid[iSet - 1][iQ2 - 1 + j][ix][iFlavour];
      xGrid[2]  = grid[iSet - 1][iQ2 - 1 + j][ix + 1][iFlavour];
      xGrid[3]  = grid[iSet - 1][iQ2 - 1 + j][ix + 2][iFlavour];
      Q2Grid[j] = polInt(xGrid, xNear, 4, x);
    }

    // Interpolate in Q2.
    double result = polInt(Q2Grid, Q2Near, 3, dQ2);

    // Save the values.
    if (iFlavour == 0) ruv = max(result, 0.);
    if (iFlavour == 1) rdv = max(result, 0.);
    if (iFlavour == 2) ru  = max(result, 0.);
    if (iFlavour == 3) rd  = max(result, 0.);
    if (iFlavour == 4) rs  = max(result, 0.);
    if (iFlavour == 5) rc  = max(result, 0.);
    if (iFlavour == 6) rb  = max(result, 0.);
    if (iFlavour == 7) rg  = max(result, 0.);

  }

}

//--------------------------------------------------------------------------

// Polynomial interpolation with Newton's divided difference method.

double EPS09::polInt(double* fi, double* xi, int n, double x) {

  for(int i = 1;i < n;i++) {
    for(int j = n-1;j > i - 1;j--) {
      fi[j] = (fi[j] - fi[j-1])/(xi[j] - xi[j-i]);
    }
  }
  double f = fi[n-1];
  for(int i = n-2;i > -1;i--) {
    f = (x - xi[i])*f + fi[i];
  }

  return f;

}

//==========================================================================

// Nuclear modifications of the PDFs from EPPS16 NLO fit.
// Ref: K.J. Eskola, P. Paakkinen, H. Paukkunen and C.A. Salgado,
// Eur.Phys.J. C77 (2017) no.3, 163 [arXiv:1612.05741]
// Grids files for different nuclei can be found from
// https://www.jyu.fi/science/en/physics/research/highenergy/urhic/npdfs/
// epps16-nuclear-pdfs.

// Constants related to the fit.
const double EPPS16::Q2MIN = 1.69;
const double EPPS16::Q2MAX = 100000000.;
const double EPPS16::XMIN  = 0.0000001;
const double EPPS16::XMAX  = 1.;
const int EPPS16::XSTEPS   = 80;
const int EPPS16::Q2STEPS  = 30;
const int EPPS16::NINTQ2   = 4;
const int EPPS16::NINTX    = 4;
const int EPPS16::NSETS    = 41;

//--------------------------------------------------------------------------

// Initialize EPPS16 nPDFs with given order (1=LO, 2=NLO) and error set.

void EPPS16::init(int iSetIn, string xmlPath) {

  // Save the error set number and derive useful values.
  iSet           = iSetIn;
  logQ2min       = log(Q2MIN);
  loglogQ2maxmin = log( log(Q2MAX)/logQ2min );
  logX2min       = log(XMIN) - 2. * (1. - XMIN);

  // Select which data file to read for current fit.
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";
  stringstream fileSS;
  fileSS << xmlPath << "EPPS16NLOR_" << getA();
  string gridFile = fileSS.str();

  // Open grid file.
  ifstream fileStream( gridFile.c_str() );
  if (!fileStream.good()) {
    printErr("Error in EPPS16::init: did not find grid file " + gridFile,
             infoPtr);
    isSet = false;
    return;
  }

  // Dump additional grid information here.
  double dummy;

  // Read in the interpolation grid.
  for (int i = 0;i < NSETS; ++i) {
    for (int j = 0;j < Q2STEPS+1; ++j) {
      fileStream >> dummy;
      for (int k = 0;k < XSTEPS; ++k) {
        for (int l = 0;l < 8; ++l) fileStream >> grid[i][j][k][l];
      }
    }
  }
  fileStream.close();

}

//--------------------------------------------------------------------------

// Interpolation from the grid.

void EPPS16::rUpdate(int , double x, double Q2) {

  // Freeze the x and Q2 values if outside the grid.
  if( x  < XMIN )  x  = XMIN;
  if( x  > XMAX )  x  = XMAX;
  if( Q2 < Q2MIN ) Q2 = Q2MIN;
  if( Q2 > Q2MAX ) Q2 = Q2MAX;

  // Do not use the points at mass threshold for interpolation.
  int cThreshold = 0;
  int bThreshold = 0;

  // Calculate the position in log(log Q^2) grid.
  double dQ2 = Q2STEPS * log( log(Q2) / logQ2min ) / loglogQ2maxmin;
  int    iQ2 = int(dQ2);

  // Set the Q2 index to interval [1,...,28].
  if      ( iQ2 < 1 )           iQ2 = 1;
  else if ( iQ2 > Q2STEPS - 3 ) iQ2 = Q2STEPS - 2;

  // Calculate the position in x grid.
  double dx = XSTEPS * ( 1. - (log(x) - 2. * (1. - x) ) / logX2min );
  int    ix = int(dx);

  // Set the x-index interval.
  if ( ix < 1 ) ix = 1;

  // Interpolate the grid values.
  for ( int iFlavour = 0; iFlavour < 8; ++iFlavour) {

    // Do not use the last grid points for interpolation.
    if ( (iFlavour > 1) && (iFlavour < 7) ) {
      if ( ix > XSTEPS - 6 ) ix = XSTEPS - 6;
    } else {
      if ( ix > XSTEPS - 4 ) ix = XSTEPS - 4;
    }

    // Calculate the four nearest points in x grid.
    double xNear[4];
    for(int i = 0;i < 4;i++) xNear[i] = ix - 1 + i;

    // Reject point Q=1.3 GeV from interpoilation for charm.
    if ( (iFlavour == 5) && (iQ2 == 1) ) {
      cThreshold = iQ2;
      iQ2        = 2;
    }

    // Reject points Q<4.75 GeV from interpoilation for bottom.
    if ( (iFlavour == 6) && (iQ2 < 17) && (iQ2 > 1) ) {
      bThreshold = iQ2;
      iQ2        = 17;
    }

    // Calculate the three nearest points in log(log Q^2) grid.
    double Q2Near[4];
    for(int i = 0;i < 4;i++) Q2Near[i] = iQ2 - 1 + i;

    // Grid points used for interpolation.
    double xGrid[4];
    double Q2Grid[4];

    // Read in the relevant values from table and interpolate in x.
    for ( int j = 0; j < 4; ++j) {
      xGrid[0]  = grid[iSet - 1][iQ2 - 1 + j][ix - 1][iFlavour];
      xGrid[1]  = grid[iSet - 1][iQ2 - 1 + j][ix][iFlavour];
      xGrid[2]  = grid[iSet - 1][iQ2 - 1 + j][ix + 1][iFlavour];
      xGrid[3]  = grid[iSet - 1][iQ2 - 1 + j][ix + 2][iFlavour];
      Q2Grid[j] = polInt(xGrid, xNear, NINTX, dx);
    }

    // Interpolate in Q2.
    double result = polInt(Q2Grid, Q2Near, NINTQ2, dQ2);

    // Save the values, for b non-zero only above the mass threshold.
    if (iFlavour == 0) ruv = result;
    if (iFlavour == 1) rdv = result;
    if (iFlavour == 2) ru  = result;
    if (iFlavour == 3) rd  = result;
    if (iFlavour == 4) rs  = result;
    if (iFlavour == 5) rc  = result;
    if (iFlavour == 6) rb  = ( sqrt(Q2) < 4.75 ) ? 0. : result;
    if (iFlavour == 7) rg  = result;

    // Revert back to original interpolation points.
    if (cThreshold > 0) {
      iQ2        = cThreshold;
      cThreshold = 0;
    } else if (bThreshold > 0) {
      iQ2        = bThreshold;
      bThreshold = 0;
    }

  }

}

//--------------------------------------------------------------------------

// Polynomial interpolation with Newton's divided difference method.

double EPPS16::polInt(double* fi, double* xi, int n, double x) {

  for(int i = 1;i < n;i++) {
    for(int j = n-1;j > i - 1;j--) {
      fi[j] = (fi[j] - fi[j-1])/(xi[j] - xi[j-i]);
    }
  }
  double f = fi[n-1];
  for(int i = n-2;i > -1;i--) {
    f = (x - xi[i])*f + fi[i];
  }

  return f;

}

//==========================================================================

} // end namespace Pythia8
