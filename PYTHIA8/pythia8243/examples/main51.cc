// main51.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test of LHAPDF interface and whether PDF's behave sensibly.
// October 2017: updated to test external LHAPDF6 vs internal LHAGrid1.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// Integration to check momentum sum rule.

double integrate(PDF* nowPDF, double Q2) {

  // Number of points, x ranges and initial values.
  int    nLin  = 980;
  int    nLog  = 1000;
  double xLin  = 0.02;
  double xLog  = 1e-8;
  double dxLin = (1. - xLin) / nLin;
  double dxLog = log(xLin / xLog) / nLog;
  double sum   = 0.;
  double x, sumNow;

  // Integration at large x in linear steps.
  for (int iLin = 0; iLin < nLin; ++iLin) {
    x      = xLin + (iLin + 0.5) * dxLin;
    sumNow = nowPDF->xf( 21, x, Q2) + nowPDF->xf( 22, x, Q2);
    for (int i = 1; i < 6; ++i)
      sumNow += nowPDF->xf( i, x, Q2) + nowPDF->xf( -i, x, Q2);
    sum   += dxLin * sumNow;
  }

  // Integration at small x in logarithmic steps.
  for (int iLog = 0; iLog < nLog; ++iLog) {
    x      = xLog * pow( xLin / xLog, (iLog + 0.5) / nLog );
    sumNow = nowPDF->xf( 21, x, Q2) + nowPDF->xf( 22, x, Q2);
    for (int i = 1; i < 6; ++i)
      sumNow += nowPDF->xf( i, x, Q2) + nowPDF->xf( -i, x, Q2);
    sum   += dxLog * x * sumNow;
  }

  // Done.
  return sum;

}

//==========================================================================

int main() {

  // Info member for possible error printouts etc.
  Info info;

  // Pointers to external LHAPDF6 and internal LHAGrid1 PDF packages,
  // for the same  NNPDF3.1 QCD+QED NNLOPDF set, the central member.
  PDF* extPDF = new LHAPDF( 2212, "LHAPDF6:NNPDF31_nnlo_as_0118_luxqed",
    &info);
  PDF* intPDF = new LHAGrid1( 2212, "20", "../share/Pythia8/xmldoc/", &info);

  // Alternative: compare two Pomeron PDF's. Boost second by factor 2.
  //PDF* extPDF = new PomFix( 990, -0.2, 2.5, 0., 3., 0.4, 0.5);
  //PDF* intPDF = new PomH1Jets( 990, 2.);
  //PDF* extPDF = new PomH1FitAB( 990, 2);
  //PDF* intPDF = new PomH1FitAB( 990, 3);

  // Allow extrapolation of PDF's beyond x and Q2 boundaries, at own risk.
  // Default behaviour is to freeze PDF's at boundaries.
  intPDF->setExtrapolate(true);
  extPDF->setExtrapolate(true);

  // Histogram F(x, Q2) = (9/4) x*g(x, Q2) + sum_{i = q, qbar} x*f_i(x, Q2)
  // for range 10^{-8} < x < 1 logarithmic in x and for Q2 = 4 and 100.
  Hist extF4("F( x, Q2 = 4) external", 80 , 1e-8, 1., true);
  Hist intF4("F( x, Q2 = 4) internal", 80 , 1e-8, 1., true);
  Hist ratF4("F( x, Q2 = 4) internal/external", 80 , 1e-8, 1., true);
  Hist extF100("F( x, Q2 = 100) external", 80 , 1e-8, 1., true);
  Hist intF100("F( x, Q2 = 100) internal", 80 , 1e-8, 1., true);
  Hist ratF100("F( x, Q2 = 100) internal/external", 80 , 1e-8, 1., true);

  // Loop over the two Q2 values.
  for (int iQ = 0; iQ < 2; ++iQ) {
    double Q2 = (iQ == 0) ? 4. : 100;

    // Loop over x values, in a logarithmic scale.
    for (int iX = 0; iX < 80; ++iX) {
      double xLog = -(0.1 * iX + 0.05);
      double x = pow( 10., xLog);

      // Evaluate external summed PDF, with colour factor 9/4 for gluons.
      double extSum = 2.25 * extPDF->xf( 21, x, Q2);
      for (int i = 1; i < 6; ++i)
        extSum += extPDF->xf( i, x, Q2) + extPDF->xf( -i, x, Q2);
      if (iQ == 0) extF4.fill ( x, extSum );
      else       extF100.fill ( x, extSum );

      // Evaluate internal summed PDF, with colour factor 9/4 for gluons.
      double intSum = 2.25 * intPDF->xf( 21, x, Q2);
      for (int i = 1; i < 6; ++i)
        intSum += intPDF->xf( i, x, Q2) + intPDF->xf( -i, x, Q2);
      if (iQ == 0) intF4.fill ( x, intSum );
      else       intF100.fill ( x, intSum );

    // End loops over x and Q2 values.
    }
  }

  // Show F(x, Q2) and their ratio internal/external.
  ratF4 = intF4 / extF4;
  ratF100 = intF100 / extF100;
  cout << extF4 << intF4 << ratF4 << extF100 << intF100 << ratF100;

  // Histogram momentum sum as a function of Q2.
  Hist extXSum("momentum sum(Q2) - 1 external", 100, 1e-2, 1e8, true);
  Hist intXSum("momentum sum(Q2) - 1 internal", 100, 1e-2, 1e8, true);
  Hist difXSum("momentum sum(Q2) internal - external", 100, 1e-2, 1e8, true);

  // Loop over Q2 values.
  for (int iQ = 0; iQ < 100; ++iQ) {
    double log10Q2 = -2.0 + 0.1 * iQ + 0.05;
    double Q2 = pow( 10., log10Q2);

    // Evaluate external and internal momentum sums.
    double extSum = integrate( extPDF, Q2);
    extXSum.fill( Q2, extSum - 1.);
    double intSum = integrate( intPDF, Q2);
    intXSum.fill( Q2, intSum - 1.);
  }

  // Show momentum sum as a function of Q2.
  difXSum = intXSum - extXSum;
  cout << extXSum << intXSum << difXSum;

  // Write Python code that can generate a PDF file with the distributions.
  // Note: curve and histogram style deliberately mixed for clarity.
  HistPlot hpl("main51plot");
  hpl.frame( "out51plot", "Summed PDF distribution at $Q^2 = 4$", "$x$",
    "$(9/4)x g(x, Q^2) + \\sum_{q} (xq(x, Q^2) + x\\overline{q}(x, Q^2)$");
  hpl.add( extF4, "-", "LHAPDF6");
  hpl.add( intF4, "h,red", "internal");
  hpl.plot();
  hpl.frame( "", "Summed PDF distribution at $Q^2 = 100$", "$x$",
    "$(9/4)x g(x, Q^2) + \\sum_{q} (xq(x, Q^2) + x\\overline{q}(x, Q^2)$");
  hpl.add( extF100, "-", "LHAPDF6");
  hpl.add( intF100, "h,red", "internal");
  hpl.plot();
  hpl.frame( "", "Momentum sum as a function of $Q^2$", "$Q^2$",
    "$( \\int_0^1 \\sum_{q,\\overline{q},g, \\gamma} xf_i(x, Q^2)"
    "\\, \\mathrm{d}x ) - 1$");
  hpl.add( extXSum, "-", "LHAPDF6");
  hpl.add( intXSum, "h,red", "internal");
  hpl.plot();

  // Done.
  delete extPDF;
  delete intPDF;
  return 0;
}
