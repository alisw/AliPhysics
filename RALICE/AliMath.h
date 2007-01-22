#ifndef ALIMATH_H
#define ALIMATH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include <math.h>
 
#include "TObject.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
 
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
  TF1 Chi2Dist(Int_t ndf) const;                             // Provide the Chi-squared distribution function
  TF1 StudentDist(Double_t ndf) const;                       // Provide the Student's T distribution function
  TF1 FratioDist(Int_t ndf1,Int_t ndf2) const;               // Provide the Fratio distribution function
  TF1 BinomialDist(Int_t n,Double_t p) const;                // Provide the Binomial distribution function
  TF1 NegBinomialDist(Int_t k,Double_t p) const;             // Provide the Negative Binomial distribution function
  TF1 PoissonDist(Double_t mu) const;                        // Provide the Poisson distribution function
  Double_t Chi2Pvalue(Double_t chi2,Int_t ndf,Int_t sides=0,Int_t sigma=0,Int_t mode=1) const; // Chi-squared P-value
  Double_t StudentPvalue(Double_t t,Double_t ndf,Int_t sides=0,Int_t sigma=0) const; // Student's P-value
  Double_t FratioPvalue(Double_t f,Int_t ndf1,Int_t ndf2,Int_t sides=0,Int_t sigma=0) const; // F ratio P-value
  Double_t BinomialPvalue(Int_t k,Int_t n,Double_t p,Int_t sides=0,Int_t sigma=0,Int_t mode=0) const; // Bin. P-value
  Double_t PoissonPvalue(Int_t k,Double_t mu,Int_t sides=0,Int_t sigma=0) const; // Poisson P-value
  Double_t NegBinomialPvalue(Int_t n,Int_t k,Double_t p,Int_t sides=0,Int_t sigma=0,Int_t mode=0) const; // NegBin. P-value
  Double_t Nfac(Int_t n,Int_t mode=0) const;    // Compute n!
  Double_t LnNfac(Int_t n,Int_t mode=2) const;  // Compute ln(n!) 
  Double_t LogNfac(Int_t n,Int_t mode=2) const; // Compute log_10(n!) 

 protected:
  Double_t GamSer(Double_t a,Double_t x) const; // Compute P(a,x) via serial representation
  Double_t GamCf(Double_t a,Double_t x) const;  // Compute P(a,x) via continued fractions
  Double_t BesselI0(Double_t x) const;          // Compute modified Bessel function I_0(x)
  Double_t BesselK0(Double_t x) const;          // Compute modified Bessel function K_0(x)
  Double_t BesselI1(Double_t x) const;          // Compute modified Bessel function I_1(x)
  Double_t BesselK1(Double_t x) const;          // Compute modified Bessel function K_1(x)
 
 ClassDef(AliMath,5) // Various mathematical tools for physics analysis.
 
};
#endif
