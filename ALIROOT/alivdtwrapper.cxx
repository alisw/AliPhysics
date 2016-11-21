// A set of mathematical wrapper functions on top of vdt (fast math)
// This file can be compiled into a preloadable library that can be used 
// to overwrite calls to libm

// The wrapper functions add some basic checks to be floating point exception safe
// (since vdt has restricted range) and dispatches to libm in corner-cases

// History:
// Nov 2016  - first version sandro.wenzel@cern.ch

//*************** STUFF THAT COMES FROM LIBM AND WHICH IS NEEDED BY VDT **********
// the following ensures that no symbols from cmath + math conflict with
// the symbols defined further below
#define _GLIBCXX_CMATH 1
#define _MATH_H 1

// these are (the only) symbols (from math.h) needed by vdt
# define M_PI		3.14159265358979323846	/* pi */
# define M_PI_2		1.57079632679489661923	/* pi/2 */
# define M_PI_4		0.78539816339744830962	/* pi/4 */
// dirty trick to implement fabs (needed by some vdt functions)
namespace std {
  template <typename T> T fabs(T x){
    return (x<T(0))? -x : x;
  }
}
using std::fabs;
//**********************************************************************************

#include <climits>
#include <dlfcn.h>

#include "vdt/exp.h"    // includes the vdt exp function
#include "vdt/sin.h"    // includes the vdt exp function
#include "vdt/cos.h"    // includes the vdt exp function
#include "vdt/atan2.h"  // includes the vdt exp function
#include "vdt/log.h"    // includes the vdt exp function
#include "vdt/sincos.h" // includes the vdt exp function

// function pointers, pointing to original libm symbols
typedef double (*fptr)(double);
static void *libm = dlopen("libm.so.6", RTLD_NOW);
static fptr sinfunc = (fptr) dlsym(libm, "sin");
static fptr cosfunc = (fptr) dlsym(libm, "cos");
static fptr logfunc = (fptr) dlsym(libm, "log");
static fptr sincosfunc = (fptr) dlsym(libm, "sincos");
static fptr atan2func = (fptr) dlsym(libm, "atan2");
static fptr expfunc = (fptr) dlsym(libm, "exp");

extern "C" double exp(double);
double exp(double x){
  const double lm((1.*INT_MAX - 0.5)/vdt::details::LOG2E);
  if (fabs(x) < lm){
    return vdt::fast_exp(x);
  }
  else {
    return expfunc(2.*x);
  }
}

extern "C" double log(double);
double log(double x){
  return vdt::fast_log(x);
}

static const double CosRangeConstant = 1./1.2732395447351628;
extern "C" double sin(double);
double sin(double x){
  // vdt relies on a range of inputs that can be cast to int
  const double lm(INT_MAX*CosRangeConstant);
  if(fabs(x) < lm){
    return vdt::fast_sin(x);
  }
  else {
    return sinfunc(x); 
  }
}

extern "C" double cos(double);
double cos(double x){
  const double lm(INT_MAX*CosRangeConstant);
  if(fabs(x) < lm){
    return vdt::fast_cos(x);
  }
  else {
    return cosfunc(x);
  }
}

extern "C" void sincos(double, double *, double *);
void sincos(double x, double *s, double *c){
  vdt::fast_sincos(x, *s, *c);
}

extern "C" void sincosf(float x, float *s, float *c);
void sincosf(float x, float *s, float *c) {
  vdt::fast_sincosf(x, *s, *c);
}

extern "C" double atan2(double, double);
double atan2(double x, double y) {
  return vdt::fast_atan2(x,y);
}


