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
#include "vdt/sin.h"    // includes the vdt sin function
#include "vdt/cos.h"    // includes the vdt cos function
#include "vdt/atan2.h"  // includes the vdt atan2 function
#include "vdt/log.h"    // includes the vdt log function
#include "vdt/sincos.h" // includes the vdt sincos function

//#define TEXTLOGGER 1
#ifdef TEXTLOGGER 
#include <iostream>
#endif

const double kDiff=1E-6;
#define printLog(name, arg, value1, value2) \
  if(fabs(value2-value1) > kDiff){ std::cerr << name << " arg " << arg << " libm " << value1 << " vdt " << value2 << " diff " << value1-value2 << "\n"; } 
#define printLog2Args(name, arg1, arg2, value1, value2)			\
  if(fabs(value2-value1) > kDiff){ std::cerr << name << " arg " << arg1 << " , " << arg2 << " libm " << value1 << " vdt " << value2 << " diff " << value1-value2 << "\n"; } 

// function pointers, pointing to original libm symbols
static void *libm = dlopen("libm.so.6", RTLD_NOW);

typedef double (*fptr)(double);
static fptr sinfunc = (fptr) dlsym(libm, "sin");
static fptr cosfunc = (fptr) dlsym(libm, "cos");
static fptr logfunc = (fptr) dlsym(libm, "log");
static fptr expfunc = (fptr) dlsym(libm, "exp");

typedef float (*fptrf)(float);
static fptrf sinfuncf = (fptrf) dlsym(libm, "sinf");
static fptrf cosfuncf = (fptrf) dlsym(libm, "cosf");
static fptrf logfuncf = (fptrf) dlsym(libm, "logf");
static fptrf expfuncf = (fptrf) dlsym(libm, "expf");

typedef void (*sincosfptr)(double, double*, double*);
static sincosfptr sincosfunc = (sincosfptr) dlsym(libm, "sincos");

typedef void (*sincosffptr)(float, float*, float*);
static sincosffptr sincosffunc = (sincosffptr) dlsym(libm, "sincosf");

typedef double (*atan2fptr)(double, double);
static atan2fptr atan2func = (atan2fptr) dlsym(libm, "atan2");

extern "C" double exp(double);
double exp(double x){
  const double lm((1.*INT_MAX - 0.5)/vdt::details::LOG2E);
  if (fabs(x) < lm){
#ifdef TEXTLOGGER
    printLog("#VDTCHECK-EXP", x, expfunc(x), vdt::fast_exp(x));
#endif
    return vdt::fast_exp(x);
  }
  else {
    return expfunc(x);
  }
}

extern "C" double log(double);
double log(double x){
#ifdef TEXTLOGGER
    printLog("#VDTCHECK-LOG", x, logfunc(x), vdt::fast_log(x));
#endif
  return vdt::fast_log(x);
}

static const double CosRangeConstant = 1./1.2732395447351628;
extern "C" double sin(double);
double sin(double x){
  // vdt relies on a range of inputs that can be cast to int
  const double lm(INT_MAX*CosRangeConstant);
  if(fabs(x) < lm){
#ifdef TEXTLOGGER
    printLog("#VDTCHECK-SIN", x, sinfunc(x), vdt::fast_sin(x));
#endif
    return vdt::fast_sin(x);
  }
  else {
    return sinfunc(x); 
  }
}

extern "C" float sinf(float);
float sinf(float x){
  // vdt relies on a range of inputs that can be cast to int
  const double lm(INT_MAX*CosRangeConstant);
  if(fabs(x) < static_cast<float>(lm)){
#ifdef TEXTLOGGER
    printLog("#VDTCHECK-SINF", x, sinfuncf(x), vdt::fast_sinf(x));
#endif
    return vdt::fast_sinf(x);
  }
  else {
    return sinfuncf(x);
  }
}


extern "C" double cos(double);
double cos(double x){
  const double lm(INT_MAX*CosRangeConstant);
  if(fabs(x) < lm){
#ifdef TEXTLOGGER
    printLog("#VDTCHECK-COS", x, cosfunc(x), vdt::fast_cos(x));
#endif
    return vdt::fast_cos(x);
  }
  else {
    return cosfunc(x);
  }
}

extern "C" float cosf(float);
float cosf(float x){
  const double lm(INT_MAX*CosRangeConstant);
  if(fabs(x) < static_cast<float>(lm)){
#ifdef TEXTLOGGER
    printLog("#VDTCHECK-COSF", x, cosfuncf(x), vdt::fast_cosf(x));
#endif
    return vdt::fast_cosf(x);
  }
  else {
    return cosfuncf(x);
  }
}

extern "C" void sincos(double, double *, double *);
void sincos(double x, double *s, double *c){
  vdt::fast_sincos(x, *s, *c);
#ifdef TEXTLOGGER
  double libms, libmc;
  sincosfunc(x, &libms, &libmc);
  printLog("#VDTCHECK-SINCOS-SINPART", x, libms, *s);
  printLog("#VDTCHECK-SINCOS-COSPART", x, libmc, *c);
#endif

}

extern "C" void sincosf(float x, float *s, float *c);
void sincosf(float x, float *s, float *c) {
  vdt::fast_sincosf(x, *s, *c);
#ifdef TEXTLOGGER
  float libms, libmc;
  sincosffunc(x, &libms, &libmc);
  printLog("#VDTCHECK-SINCOSF-SINPART", x, libms, *s);
  printLog("#VDTCHECK-SINCOSF-COSPART", x, libmc, *c);
#endif
}

extern "C" double atan2(double, double);
double atan2(double x, double y) {
#ifdef TEXTLOGGER
  printLog2Args("#VDTCHECK-ATAN2", x, y, atan2func(x,y), vdt::fast_atan2(x,y));
#endif
  return vdt::fast_atan2(x,y);
}


