#ifndef ALIMATH_H
#define ALIMATH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include <math.h>
 
#include "TObject.h"
 
class AliMath : public TObject
{
 public:
  AliMath();                                                 // Default constructor
  virtual ~AliMath();                                        // Destructor
  AliMath(const AliMath& m);                                 // Copy constructor
  Double_t Gamma(Double_t z) const;                          // Standard gamma function Gamma(z)
  Double_t Gamma(Double_t a,Double_t x) const;               // Incomplete gamma function P(a,x)
  Double_t LnGamma(Double_t z) const;                        // Compute ln[Gamma(z)]
  Double_t Erf(Double_t x) const;                            // Error function erf(x)
  Double_t Erfc(Double_t x) const;                           // Complementary error function erfc(x)
  Double_t Prob(Double_t chi2,Int_t ndf,Int_t mode=1) const; // Compute Chi-squared probability
  Double_t BesselI(Int_t n,Double_t x) const;                // Compute integer order mod. Bessel function I_n(x)
  Double_t BesselK(Int_t n,Double_t x) const;                // Compute integer order mod. Bessel function K_n(x)
 
 protected:
  Double_t GamSer(Double_t a,Double_t x) const; // Compute P(a,x) via serial representation
  Double_t GamCf(Double_t a,Double_t x) const;  // Compute P(a,x) via continued fractions
  Double_t BesselI0(Double_t x) const;          // Compute modified Bessel function I_0(x)
  Double_t BesselK0(Double_t x) const;          // Compute modified Bessel function K_0(x)
  Double_t BesselI1(Double_t x) const;          // Compute modified Bessel function I_1(x)
  Double_t BesselK1(Double_t x) const;          // Compute modified Bessel function K_1(x)
 
 ClassDef(AliMath,3) // Various mathematical tools for physics analysis.
 
};
#endif
