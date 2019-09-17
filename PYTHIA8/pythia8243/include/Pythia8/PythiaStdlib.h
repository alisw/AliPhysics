// PythiaStdlib.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for functionality pulled in from Stdlib,
// plus a few useful utilities (small powers; positive square root,
// convert strings to lowercase, Gamma function).

#ifndef Pythia8_PythiaStdlib_H
#define Pythia8_PythiaStdlib_H

// Stdlib header files for mathematics.
#include <cmath>
#include <cstdlib>
#include <algorithm>

// Stdlib header files for strings and containers.
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <set>
#include <list>

// Stdlib header file for dynamic library loading.
#define dlsym __
#include <dlfcn.h>
#undef dlsym

// Redefine dlsym to suppress compiler warnings.
extern "C" void *(*dlsym(void *handle, const char *symbol))();

// Stdlib header file for input and output.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

// Define pi if not yet done.
#ifndef M_PI
#define M_PI 3.1415926535897932385
#endif

// By this declaration you do not need to use std:: qualifier everywhere.
//using namespace std;

// Alternatively you can specify exactly which std:: methods will be used.
// Now made default so std does not spill outside namespace Pythia8.
namespace Pythia8 {

// Generic utilities and mathematical functions.
using std::swap;
using std::max;
using std::min;
using std::abs;
using std::sort;

// Strings and containers.
using std::pair;
using std::make_pair;
using std::string;
using std::vector;
using std::map;
using std::multimap;
using std::deque;
using std::set;
using std::multiset;
using std::list;

// Input/output streams.
using std::cin;
using std::cout;
using std::cerr;
using std::istream;
using std::ostream;
using std::fstream;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::ios;

// Input/output formatting.
using std::endl;
using std::fixed;
using std::scientific;
using std::left;
using std::right;
using std::setw;
using std::setprecision;

} // end namespace Pythia8

namespace Pythia8 {

// Define conversion hbar * c = 0.2 GeV * fm = 1.
#ifndef HBARC
#define HBARC 0.19732698
#endif

// Define conversion between fm and mm, in both directions.
#ifndef FM2MM
#define FM2MM 1e-12
#endif
#ifndef MM2FM
#define MM2FM 1e12
#endif

// Powers of small integers - for balance speed/code clarity.
inline double pow2(const double& x) {return x*x;}
inline double pow3(const double& x) {return x*x*x;}
inline double pow4(const double& x) {return x*x*x*x;}
inline double pow5(const double& x) {return x*x*x*x*x;}
inline double pow6(const double& x) {return x*x*x*x*x*x;}
inline double pow7(const double& x) {return x*x*x*x*x*x*x;}
inline double pow8(const double& x) {return x*x*x*x*x*x*x*x;}

// Avoid problem with negative square root argument (from roundoff).
inline double sqrtpos(const double& x) {return sqrt( max( 0., x));}

// Convert a string to lowercase for case-insensitive comparisons.
// By default remove any initial and trailing blanks or escape characters.
string toLower(const string& name, bool trim = true);

// Variant of above, with in-place replacement.
inline void toLowerRep(string& name, bool trim = true) {
  name = toLower( name, trim);}

// The Gamma function for real argument.
double GammaReal(double x);

// Modified Bessel functions of the first and second kinds.
double besselI0(double x);
double besselI1(double x);
double besselK0(double x);
double besselK1(double x);

// Base class to encapsulate a (double) function of an arbitrary number
// of (double) arguments (to avoid using function pointers).
class FunctionEncapsulator  {
public:
  FunctionEncapsulator() {};
  virtual ~FunctionEncapsulator() { };
  virtual double f(vector<double> args);
  // Integrate over function argument iArg, using Gaussian quadrature.
  bool integrateGauss(double& result, int iArg, double xLo, double xHi,
    vector<double> args, double tol=1e-6);
  // Solve f(args) = target for argument iArg, using Brent's method
  bool brent(double& solution, double target, int iArg, double xLo, double xHi,
    vector<double> argsIn, double tol=1e-6, int maxIter=10000);
};

} // end namespace Pythia8

#endif // Pythia8_PythiaStdlib_H
