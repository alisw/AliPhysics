// -*- C++ -*-
// CLASSDOC OFF
// $Id$
// ---------------------------------------------------------------------------
// CLASSDOC ON
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This file contains definitions of some usefull utilities and macros.
//
#ifndef _CLHEP_H_
#define _CLHEP_H_

#include <stdlib.h>
#include <limits.h>
#include <math.h>

#if defined(CLHEP_TARGET_H)
#include CLHEP_TARGET_H
#else
#include "CLHEP/config/CLHEP-default.h"
#endif

// CLASSDOC OFF
// **** You should probably not touch anything below this line: ****

typedef double HepDouble;
typedef int    HepInt;
typedef float  HepFloat;

#ifdef HEP_HAVE_BOOL
typedef bool HepBoolean;
#else
typedef int HepBoolean;
#ifndef false
const HepBoolean hep_false = 0;
#define false hep_false
#endif
#ifndef true
const HepBoolean hep_true  = 1;
#define true hep_true
#endif
#endif /* HEP_HAVE_BOOL */

#ifdef HEP_SHORT_NAMES
typedef HepBoolean Boolean;
#endif

#ifndef M_PI_2
#define M_PI_2	1.57079632679489661923
#endif

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#ifndef M_2PI
#define M_2PI   6.28318530717958647692
#endif

#ifndef CLHEP_MAX_MIN_DEFINED
#define CLHEP_MAX_MIN_DEFINED
template <class T>
inline const T& min(const T& a, const T& b) {
  // Break this into two lines to avoid an incorrect warning with
  // Cfront-based compilers.
  const T& retval = b < a ? b : a;
  return retval;
}

template <class T>
inline const T& max(const T& a, const T& b) {
  // Break this into two lines to avoid an incorrect warning with
  // Cfront-based compilers.
  const T& retval = a < b ? b : a;
  return retval;
}
#endif

#ifndef CLHEP_SQR_ABS_DEFINED
#define CLHEP_SQR_ABS_DEFINED
template <class T>
inline T sqr(const T& x) {
  return x*x;
}

template <class T>
inline T abs(const T& a) {
  return a < 0 ? -a : a;
}
#endif

#ifdef HEP_DEBUG_INLINE
#define HEP_NO_INLINE_IN_DECLARATION
#endif

#ifdef HEP_NO_INLINE_IN_DECLARATION
#define HEP_NO_INLINE_IN_TEMPLATE_DECLARATION
#endif

// Default to generate random matrix
//
#ifndef HEP_USE_RANDOM
#define HEP_USE_RANDOM
#endif

// GNU g++ compiler can optimize when returning an object.
// However g++ on HP cannot deal with this.
//
#undef HEP_GNU_OPTIMIZED_RETURN

// CLASSDOC ON
#endif /* _CLHEP_H_ */
