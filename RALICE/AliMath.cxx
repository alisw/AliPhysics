/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.2  1999/09/29 09:24:28  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////
// Class AliMath
// Various mathematical tools which may be very convenient while
// performing physics analysis.
//
// Example : Probability of a Chi-squared value
// =========
//
// AliMath M;
// Float_t chi2=20;            // The chi-squared value
// Int_t ndf=12;               // The number of degrees of freedom
// Float_t p=M.Prob(chi2,ndf); // The probability that at least a Chi-squared
//                             // value of chi2 will be observed, even for a
//                             // correct model
//
//--- Author: Nick van Eijndhoven 14-nov-1998 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliMath.h"
 
ClassImp(AliMath) // Class implementation to enable ROOT I/O
 
AliMath::AliMath()
{
// Default constructor
}
///////////////////////////////////////////////////////////////////////////
AliMath::~AliMath()
{
// Destructor
}
///////////////////////////////////////////////////////////////////////////
Float_t AliMath::Gamma(Float_t z)
{
// Computation of gamma(z) for all z>0.
//
// The algorithm is based on the article by C.Lanczos [1] as denoted in
// Numerical Recipes 2nd ed. on p. 207 (W.H.Press et al.).
//
// [1] C.Lanczos, SIAM Journal of Numerical Analysis B1 (1964), 86.
//
//--- Nve 14-nov-1998 UU-SAP Utrecht
 
 if (z<=0.)
 {
  cout << "*Gamma(z)* Wrong argument z = " << z << endl;
  return 0;
 }
 
 Float_t v=LnGamma(z);
 return exp(v);
}
///////////////////////////////////////////////////////////////////////////
Float_t AliMath::Gamma(Float_t a,Float_t x)
{
// Computation of the incomplete gamma function P(a,x)
//
// The algorithm is based on the formulas and code as denoted in
// Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
//
//--- Nve 14-nov-1998 UU-SAP Utrecht
 
 if (a<=0.)
 {
  cout << "*Gamma(a,x)* Invalid argument a = " << a << endl;
  return 0;
 }
 
 if (x<=0.)
 {
  if (x<0) cout << "*Gamma(a,x)* Invalid argument x = " << x << endl;
  return 0;
 }
 
 if (x<(a+1.))
 {
  return GamSer(a,x);
 }
 else
 {
  return GamCf(a,x);
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliMath::LnGamma(Float_t z)
{
// Computation of ln[gamma(z)] for all z>0.
//
// The algorithm is based on the article by C.Lanczos [1] as denoted in
// Numerical Recipes 2nd ed. on p. 207 (W.H.Press et al.).
//
// [1] C.Lanczos, SIAM Journal of Numerical Analysis B1 (1964), 86.
//
// The accuracy of the result is better than 2e-10.
//
//--- Nve 14-nov-1998 UU-SAP Utrecht
 
 if (z<=0.)
 {
  cout << "*LnGamma(z)* Wrong argument z = " << z << endl;
  return 0;
 }
 
 // Coefficients for the series expansion
 Double_t c[7];
 c[0]=  2.5066282746310005;
 c[1]= 76.18009172947146;
 c[2]=-86.50532032941677;
 c[3]= 24.01409824083091;
 c[4]= -1.231739572450155;
 c[5]=  0.1208650973866179e-2;
 c[6]= -0.5395239384953e-5;
 
 Double_t x=z;
 Double_t y=x;
 Double_t tmp=x+5.5;
 tmp=(x+0.5)*log(tmp)-tmp;
 Double_t ser=1.000000000190015;
 for (Int_t i=1; i<7; i++)
 {
  y+=1.;
  ser+=c[i]/y;
 }
 Float_t v=tmp+log(c[0]*ser/x);
 return v;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliMath::GamSer(Float_t a,Float_t x)
{
// Computation of the incomplete gamma function P(a,x)
// via its series representation.
//
// The algorithm is based on the formulas and code as denoted in
// Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
//
//--- Nve 14-nov-1998 UU-SAP Utrecht
 
 Int_t itmax=100;   // Maximum number of iterations
 Float_t eps=3.e-7; // Relative accuracy
 
 if (a<=0.)
 {
  cout << "*GamSer(a,x)* Invalid argument a = " << a << endl;
  return 0;
 }
 
 if (x<=0.)
 {
  if (x<0) cout << "*GamSer(a,x)* Invalid argument x = " << x << endl;
  return 0;
 }
 
 Float_t gln=LnGamma(a);
 Float_t ap=a;
 Float_t sum=1./a;
 Float_t del=sum;
 for (Int_t n=1; n<=itmax; n++)
 {
  ap+=1.;
  del=del*x/ap;
  sum+=del;
  if (fabs(del)<fabs(sum*eps)) break;
  if (n==itmax) cout << "*GamSer(a,x)* a too large or itmax too small" << endl;
 }
 Float_t v=sum*exp(-x+a*log(x)-gln);
 return v;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliMath::GamCf(Float_t a,Float_t x)
{
// Computation of the incomplete gamma function P(a,x)
// via its continued fraction representation.
//
// The algorithm is based on the formulas and code as denoted in
// Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
//
//--- Nve 14-nov-1998 UU-SAP Utrecht
 
 Int_t itmax=100;      // Maximum number of iterations
 Float_t eps=3.e-7;    // Relative accuracy
 Float_t fpmin=1.e-30; // Smallest Float_t value allowed here
 
 if (a<=0.)
 {
  cout << "*GamCf(a,x)* Invalid argument a = " << a << endl;
  return 0;
 }
 
 if (x<=0.)
 {
  if (x<0) cout << "*GamCf(a,x)* Invalid argument x = " << x << endl;
  return 0;
 }
 
 Float_t gln=LnGamma(a);
 Float_t b=x+1.-a;
 Float_t c=1./fpmin;
 Float_t d=1./b;
 Float_t h=d;
 Float_t an,del;
 for (Int_t i=1; i<=itmax; i++)
 {
  an=float(-i)*(float(i)-a);
  b+=2.;
  d=an*d+b;
  if (fabs(d)<fpmin) d=fpmin;
  c=b+an/c;
  if (fabs(c)<fpmin) c=fpmin;
  d=1./d;
  del=d*c;
  h=h*del;
  if (fabs(del-1.)<eps) break;
  if (i==itmax) cout << "*GamCf(a,x)* a too large or itmax too small" << endl;
 }
 Float_t v=exp(-x+a*log(x)-gln)*h;
 return (1.-v);
}
///////////////////////////////////////////////////////////////////////////
Float_t AliMath::Erf(Float_t x)
{
// Computation of the error function erf(x).
//
//--- NvE 14-nov-1998 UU-SAP Utrecht
 
 return (1.-Erfc(x));
}
///////////////////////////////////////////////////////////////////////////
Float_t AliMath::Erfc(Float_t x)
{
// Computation of the complementary error function erfc(x).
//
// The algorithm is based on a Chebyshev fit as denoted in
// Numerical Recipes 2nd ed. on p. 214 (W.H.Press et al.).
//
// The fractional error is always less than 1.2e-7.
//
//--- Nve 14-nov-1998 UU-SAP Utrecht
 
 // The parameters of the Chebyshev fit
 const Float_t a1=-1.26551223,  a2=1.00002368,
               a3= 0.37409196,  a4=0.09678418,
               a5=-0.18628806,  a6=0.27886807,
               a7=-1.13520398,  a8=1.48851587,
               a9=-0.82215223, a10=0.17087277;
 
 Float_t v=1.; // The return value
 
 Float_t z=fabs(x);
 
 if (z <= 0.) return v; // erfc(0)=1
 
 Float_t t=1./(1.+0.5*z);
 
 v=t*exp((-z*z)
         +a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10)))))))));
 
 if (x < 0.) v=2.-v; // erfc(-x)=2-erfc(x)
 
 return v;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliMath::Prob(Float_t chi2,Int_t ndf)
{
// Computation of the probability for a certain Chi-squared (chi2)
// and number of degrees of freedom (ndf).
//
// Calculations are based on the incomplete gamma function P(a,x),
// where a=ndf/2 and x=chi2/2.
//
// P(a,x) represents the probability that the observed Chi-squared
// for a correct model should be less than the value chi2.
//
// The returned probability corresponds to 1-P(a,x),
// which denotes the probability that an observed Chi-squared exceeds
// the value chi2 by chance, even for a correct model.
//
//--- NvE 14-nov-1998 UU-SAP Utrecht
 
 if (ndf <= 0) return 0; // Set CL to zero in case ndf<=0
 
 if (chi2 <= 0.)
 {
  if (chi2 < 0.)
  {
    return 0;
  }
  else
  {
   return 1;
  }
 }
 
// Alternative which is exact
// This code may be activated in case the gamma function gives problems
// if (ndf==1)
// {
//  Float_t v=1.-Erf(sqrt(chi2)/sqrt(2.));
//  return v;
// }
 
// Gaussian approximation for large ndf
// This code may be activated in case the gamma function shows a problem
// Float_t q=sqrt(2.*chi2)-sqrt(float(2*ndf-1));
// if (n>30 && q>0.)
// {
//  Float_t v=0.5*(1.-Erf(q/sqrt(2.)));
//  return v;
// }
 
 // Evaluate the incomplete gamma function
 Float_t a=float(ndf)/2.;
 Float_t x=chi2/2.;
 return (1.-Gamma(a,x));
}
///////////////////////////////////////////////////////////////////////////
