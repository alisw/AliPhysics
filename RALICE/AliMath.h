#ifndef ALIMATH_H
#define ALIMATH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
 
class AliMath : public TObject
{
 public:
  AliMath();                              // Default constructor
  ~AliMath();                             // Destructor
  Double_t Gamma(Double_t z);             // Standard gamma function Gamma(z)
  Double_t Gamma(Double_t a,Double_t x);  // Incomplete gamma function P(a,x)
  Double_t LnGamma(Double_t z);           // Compute ln[Gamma(z)]
  Double_t Erf(Double_t x);               // Error function erf(x)
  Double_t Erfc(Double_t x);              // Complementary error function erfc(x)
  Double_t Prob(Double_t chi2,Int_t ndf); // Compute Chi-squared probability
  Double_t BesselI(Int_t n,Double_t x);   // Compute integer order modified Bessel function I_n(x)
  Double_t BesselK(Int_t n,Double_t x);   // Compute integer order modified Bessel function K_n(x)
 
 protected:
  Double_t GamSer(Double_t a,Double_t x); // Compute P(a,x) via serial representation
  Double_t GamCf(Double_t a,Double_t x);  // Compute P(a,x) via continued fractions
  Double_t BesselI0(Double_t x);          // Compute modified Bessel function I_0(x)
  Double_t BesselK0(Double_t x);          // Compute modified Bessel function K_0(x)
  Double_t BesselI1(Double_t x);          // Compute modified Bessel function I_1(x)
  Double_t BesselK1(Double_t x);          // Compute modified Bessel function K_1(x)
 
 ClassDef(AliMath,1) // Various mathematical tools for physics analysis.
 
};
#endif
