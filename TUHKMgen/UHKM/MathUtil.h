#ifndef NAMathUtil_h
#define NAMathUtil_h 1

#include <iostream>
#include <math.h>
#include <Rtypes.h>
//
//##############################################################################
//
//            Nikolai Amelin (C) JINR/Dubna 1999
//
//##############################################################################
//

const Double_t GeV   = 1.;
const Double_t fermi = 1.;
const Double_t hbarc = 0.197*GeV*fermi;
const Double_t N_PI  = 3.14159265359;

const Double_t N_INFINITY = 9.0E99;
const Double_t N_SMALL = 1.E-10;
   

template <class T> inline void SwapObj(T* a, T* b)
{
   T tmp= *a;
   *a = *b;
   *b = tmp;
}

template <class T> inline void Swap(T& a, T& b)
{
   T tmp = a;
   a = b;
   b = tmp;
}

template <class T> inline T Min(T a, T b)
{
   return (a < b) ? a : b;
}

template <class T> inline T Max(T a, T b)
{
   return (a > b) ? a : b;
}

template <class T> inline T Abs(T a)
{
   return (a > 0) ? a : -a;
}

template <class T> inline T Sign(T A, T B)
{
   return (B > 0) ? Abs(A) : -Abs(A);
}
template <class T> inline T min(T a, T b) 
  {
  return (a < b)?a:b;
  }

template <class T> inline T max(T a, T b) 
  {
  return (a > b)?a:b;
  }
/*
inline Double_t Rand(void)
{
   return ((Double_t)(rand() + 1))/(Double_t)(RAND_MAX + 2);//Visual C++
//   return ((Double_t)(-rand() + 1))/(Double_t)(RAND_MAX + 2);// Linux
}


inline Double_t sqr(Double_t Value) { return Value*Value;}
*/
#endif
