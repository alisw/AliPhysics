// PythiaStdlib.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Gamma function.

#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Convert string to lowercase for case-insensitive comparisons.
// By default remove any initial and trailing blanks or escape characters.

string toLower(const string& name, bool trim) {

  // Copy string without initial and trailing blanks or escape characters.
  string temp = name;
  if (trim) {
    if (name.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return "";
    int firstChar = name.find_first_not_of(" \n\t\v\b\r\f\a");
    int lastChar  = name.find_last_not_of(" \n\t\v\b\r\f\a");
    temp          = name.substr( firstChar, lastChar + 1 - firstChar);
  }

  // Convert to lowercase letter by letter.
  for (int i = 0; i < int(temp.length()); ++i) temp[i] = tolower(temp[i]);
  return temp;

}

//--------------------------------------------------------------------------

// The Gamma function for real arguments, using the Lanczos approximation.
// Code based on http://en.wikipedia.org/wiki/Lanczos_approximation

double GammaCoef[9] = {
     0.99999999999980993,     676.5203681218851,   -1259.1392167224028,
      771.32342877765313,   -176.61502916214059,    12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

double GammaReal(double x) {

  // Reflection formula (recursive!) for x < 0.5.
  if (x < 0.5) return M_PI / (sin(M_PI * x) * GammaReal(1 - x));

  // Iterate through terms.
  double z = x - 1.;
  double gamma = GammaCoef[0];
  for (int i = 1; i < 9; ++i) gamma += GammaCoef[i] / (z + i);

  // Answer.
  double t = z + 7.5;
  gamma *= sqrt(2. * M_PI) * pow(t, z + 0.5) * exp(-t);
  return gamma;

}

//--------------------------------------------------------------------------

// Polynomial approximation for modified Bessel function of the first kind.
// Based on Abramowitz & Stegun, Handbook of mathematical functions (1964).

double besselI0(double x){

  // Parametrize in terms of t.
  double result = 0.;
  double t  = x / 3.75;
  double t2 = pow2(t);

  // Only positive values relevant.
  if      ( t < 0.) ;
  else if ( t < 1.) {
    result = 1.0 + 3.5156229 * t2 + 3.0899424 * pow2(t2)
           + 1.2067492 * pow3(t2) + 0.2659732 * pow4(t2)
           + 0.0360768 * pow5(t2) + 0.0045813 * pow6(t2);
  } else {
    double u = 1. / t;
    result = exp(x) / sqrt(x) * ( 0.39894228 + 0.01328592 * u
           + 0.00225319 * pow2(u) - 0.00157565 * pow3(u)
           + 0.00916281 * pow4(u) - 0.02057706 * pow5(u)
           + 0.02635537 * pow6(u) - 0.01647633 * pow7(u)
           + 0.00392377 * pow8(u) );
  }

  return result;
}

//--------------------------------------------------------------------------

// Polynomial approximation for modified Bessel function of the first kind.
// Based on Abramowitz & Stegun, Handbook of mathematical functions (1964).

double besselI1(double x){

  // Parametrize in terms of t.
  double result = 0.;
  double t = x / 3.75;
  double t2 = pow2(t);

  // Only positive values relevant.
  if      ( t < 0.) ;
  else if ( t < 1.) {
    result = x * ( 0.5 + 0.87890594 * t2 + 0.51498869 * pow2(t2)
           + 0.15084934 * pow3(t2) + 0.02658733 * pow4(t2)
           + 0.00301532 * pow5(t2) + 0.00032411 * pow6(t2) );
  } else {
    double u = 1. / t;
    result = exp(x) / sqrt(x) * ( 0.39894228 - 0.03988024 * u
           - 0.00368018 * pow2(u) + 0.00163801 * pow3(u)
           - 0.01031555 * pow4(u) + 0.02282967 * pow5(u)
           - 0.02895312 * pow6(u) + 0.01787654 * pow7(u)
           - 0.00420059 * pow8(u) );
  }

  return result;
}

//--------------------------------------------------------------------------

// Polynomial approximation for modified Bessel function of a second kind.
// Based on Abramowitz & Stegun, Handbook of mathematical functions (1964).

double besselK0(double x){

  double result = 0.;

  // Polynomial approximation valid ony for x > 0.
  if      ( x < 0.) ;
  else if ( x < 2.) {
    double x2 = pow2(0.5 * x);
    result = -log(0.5 * x) * besselI0(x) - 0.57721566
           + 0.42278420 * x2 + 0.23069756 * pow2(x2)
           + 0.03488590 * pow3(x2) + 0.00262698 * pow4(x2)
           + 0.00010750 * pow5(x2) + 0.00000740 * pow6(x2);
  } else {
    double z = 2. / x;
    result = exp(-x) / sqrt(x) * ( 1.25331414 - 0.07832358 * z
           + 0.02189568 * pow2(z) - 0.01062446 * pow3(z)
           + 0.00587872 * pow4(z) - 0.00251540 * pow5(z)
           + 0.00053208 * pow6(z) );
  }

  return result;
}

//--------------------------------------------------------------------------

// Polynomial approximation for modified Bessel function of a second kind.
// Based on Abramowitz & Stegun, Handbook of mathematical functions (1964).

double besselK1(double x){

  double result = 0.;

  // Polynomial approximation valid ony for x > 0.
  if      ( x < 0.) ;
  else if ( x < 2.) {
    double x2 = pow2(0.5 * x);
    result = log(0.5 * x) * besselI1(x) + 1./x * ( 1. + 0.15443144 * x2
           - 0.67278579 * pow2(x2) - 0.18156897 * pow3(x2)
           - 0.01919402 * pow4(x2) - 0.00110404 * pow5(x2)
           - 0.00004686 * pow6(x2) );
  } else {
    double z = 2. / x;
    result = exp(-x) / sqrt(x) * ( 1.25331414 + 0.23498619 * z
           - 0.03655620 * pow2(z) + 0.01504268 * pow3(z)
           - 0.00780353 * pow4(z) + 0.00325614 * pow5(z)
           - 0.00068245 * pow6(z) );
  }

  return result;
}

//==========================================================================

// FunctionEncapsulator class.
// This class serves to encapsulate a function of an arbitrary number of
// arguments, all given as doubles. The intended use is for numerical
// integration (a Gaussian quadrature routine is implemented in
// the base class), root finding, etc, without having to use function
// pointers (so will probably become obsolete when moving to C++11.)

//--------------------------------------------------------------------------

// Definition of the function to be encapsulated; base class returns 0.

double FunctionEncapsulator::f(vector<double>) {
  return 0.0;
}

//--------------------------------------------------------------------------

// Integrate the encapsulated function, f, over argument number iArg,
// from xLo, to xHi, using Gaussian quadrature, with tolerance tol.
// Return false if precision target could not be reached.
// Adapted from the CERNLIB DGAUSS routine by K.S. Kolbig.

bool FunctionEncapsulator::integrateGauss(double& result, int iArg, double xLo,
  double xHi, vector<double> args, double tol) {

  // Initialize.
  result = 0.0;

  // Sanity check.
  if (iArg >= int(args.size())) return false;

  // Boundary check: return zero if xLo >= xHi.
  if (xLo >= xHi) return true;

  // 8-point unweighted.
  static double x8[4]={  0.96028985649753623, 0.79666647741362674,
    0.52553240991632899, 0.18343464249564980};
  static double w8[4]={  0.10122853629037626, 0.22238103445337447,
    0.31370664587788729, 0.36268378337836198};
  // 16-point unweighted.
  static double x16[8]={ 0.98940093499164993, 0.94457502307323258,
    0.86563120238783174, 0.75540440835500303, 0.61787624440264375,
    0.45801677765722739, 0.28160355077925891, 0.09501250983763744};
  static double w16[8]={  0.027152459411754095, 0.062253523938647893,
    0.095158511682492785, 0.12462897125553387, 0.14959598881657673,
    0.16915651939500254,  0.18260341504492359, 0.18945061045506850};

  // Set up integration region.
  double c = 0.001/abs(xHi-xLo);
  double zLo = xLo;
  double zHi = xHi;

  bool nextbin = true;
  while ( nextbin ) {

    double zMid = 0.5*(zHi+zLo); // midpoint
    double zDel = 0.5*(zHi-zLo); // midpoint, relative to zLo

    // Calculate 8-point and 16-point quadratures.
    double s8=0.0;
    for (int i=0;i<4;i++) {
      double dz = zDel * x8[i];
      args[iArg] = zMid+dz;
      double f1 = f(args);
      args[iArg] = zMid-dz;
      double f2 = f(args);
      s8 += w8[i]*(f1 + f2);
    }
    s8 *= zDel;
    double s16=0.0;
    for (int i=0;i<8;i++) {
      double dz = zDel * x16[i];
      args[iArg] = zMid+dz;
      double f1 = f(args);
      args[iArg] = zMid-dz;
      double f2 = f(args);
      s16 += w16[i]*(f1 + f2);
    }
    s16 *= zDel;

    // Precision in this bin OK, add to cumulative and go to next.
    if (abs(s16-s8) < tol*(1+abs(s16))) {
      nextbin=true;
      result += s16;
      // Next bin: LO = end of current, HI = end of integration region.
      zLo=zHi;
      zHi=xHi;
      if ( zLo == zHi ) nextbin = false;

    // Precision in this bin not OK, subdivide.
    } else {
      if (1.0 + c*abs(zDel) == 1.0) {
        // Cannot subdivide further at double precision. Fail.
        cout << "\n FunctionEncapsulator::integrateGauss(): cannot "
             << "reach desired tolerance at double precision." << endl;
        result = 0.0 ;
        return false;
      }
      zHi = zMid;
      nextbin = true;
    }
  }

  return true;
}

//--------------------------------------------------------------------------

// Solve f(args) = targetValue for argument iArg, on interval from xLo to xHi,
// using Brent's method, with tolerance tol and a maximum number of iterations
// maxIter. Return false if precision target could not be reached.

bool FunctionEncapsulator::brent(double& solution, double targetValue,
  int iArg, double xLo, double xHi, vector<double> argsIn, double tol,
  int maxIter) {

  // Initialize.
  solution = 0.0;

  // Sanity and range checks.
  if (iArg >= int(argsIn.size())) return false;
  if (xLo > xHi) return false;

  vector<double> args(argsIn);
  // Evaluate function - targetValue at lower boundary.
  args[iArg] = xLo;
  double f1 = f(args) - targetValue;
  if (abs(f1) < tol) {
    solution = xLo;
    return true;
  }
  // Evaluate function - targetValue at upper boundary.
  args[iArg] = xHi;
  double f2 = f(args) - targetValue;
  if (abs(f2) < tol) {
    solution = xHi;
    return true;
  }

  // Check if root is bracketed.
  if ( f1 * f2 > 0.0) return false;

  // Start searching for root.
  double x1 = xLo;
  double x2 = xHi;
  double x3 = 0.5 * (xLo + xHi);

  int iter=0;
  while(++iter < maxIter) {
    // Now check at x = x3.
    args[iArg] = x3;
    double f3 = f(args) - targetValue;
    // Check if tolerance on f has been reached.
    if (abs(f3) < tol) {
      solution = x3;
      return true;
    }
    // Is root bracketed in lower or upper half?
    if (f1 * f3 < 0.0) xHi = x3;
    else xLo = x3;
    // Check if tolerance on x has been reached.
    if ((xHi - xLo) < tol * (abs(xHi) < 1.0 ? xHi : 1.0)) {
      solution = 0.5 * (xLo + xHi);
      return true;
    }

    // Work out next step to take in x.
    double den = (f2 - f1) * (f3 - f1) * (f2 - f3);
    double num = x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 * (f2 - f3)
               + f1 * x2 * (f3 - f1);
    double dx = xHi - xLo;
    if (den != 0.0) dx = f3 * num / den;

    // First attempt, using gradient
    double x = x3 + dx;

    // If this was too far, just step to the middle
    if ((xHi - x) * (x - xLo) < 0.0) {
      dx = 0.5 * (xHi - xLo);
      x = xLo + dx;
    }
    if (x < x3) {
        x2 = x3;
        f2 = f3;
    }
    else {
        x1 = x3;
        f1 = f3;
    }
    x3 = x;
  }

  // Maximum number of iterations exceeded.
  return false;

}

//==========================================================================

} // end namespace Pythia8
