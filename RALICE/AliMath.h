#ifndef ALIMATH_H
#define ALIMATH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
// Class AliMath
// Various mathematical tools which may be very convenient while
// performing physics analysis.
//
// Example : Probability of a Chi-squared value
// =========
//
// AliMath M;
// Float_t chi2=20  ;          // The chi-squared value
// Int_t ndf=12;               // The number of degrees of freedom
// Float_t p=M.Prob(chi2,ndf); // The probability that at least a Chi-squared
//                             // value of chi2 will be observed, even for a
//                             // correct model
//
//--- NvE 14-nov-1998 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////
 
#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
 
class AliMath : public TObject
{
 public:
  AliMath();                            // Default constructor
  ~AliMath();                           // Destructor
  Float_t Gamma(Float_t z);             // Standard gamma function Gamma(z)
  Float_t Gamma(Float_t a,Float_t x);   // Incomplete gamma function P(a,x)
  Float_t LnGamma(Float_t z);           // Compute ln[Gamma(z)]
  Float_t Erf(Float_t x);               // Error function erf(x)
  Float_t Erfc(Float_t x);              // Complementary error function erfc(x)
  Float_t Prob(Float_t chi2,Int_t ndf); // Compute Chi-squared probability
 
 protected:
  Float_t GamSer(Float_t a,Float_t x);  // Compute P(a,x) via serial representation
  Float_t GamCf(Float_t a,Float_t x);   // Compute P(a,x) via continued fractions
 
 ClassDef(AliMath,1) // Class definition to enable ROOT I/O
 
};
#endif
