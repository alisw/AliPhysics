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
Revision 1.3  1999/11/03 14:23:17  fca
New version of RALICE introduced

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
//- Modified: NvE 12-mar-2000 UU-SAP Utrecht
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
Double_t AliMath::Gamma(Double_t z)
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
 
 Double_t v=LnGamma(z);
 return exp(v);
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::Gamma(Double_t a,Double_t x)
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
Double_t AliMath::LnGamma(Double_t z)
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
 Double_t v=tmp+log(c[0]*ser/x);
 return v;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::GamSer(Double_t a,Double_t x)
{
// Computation of the incomplete gamma function P(a,x)
// via its series representation.
//
// The algorithm is based on the formulas and code as denoted in
// Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
//
//--- Nve 14-nov-1998 UU-SAP Utrecht
 
 Int_t itmax=100;   // Maximum number of iterations
 Double_t eps=3.e-7; // Relative accuracy
 
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
 
 Double_t gln=LnGamma(a);
 Double_t ap=a;
 Double_t sum=1./a;
 Double_t del=sum;
 for (Int_t n=1; n<=itmax; n++)
 {
  ap+=1.;
  del=del*x/ap;
  sum+=del;
  if (fabs(del)<fabs(sum*eps)) break;
  if (n==itmax) cout << "*GamSer(a,x)* a too large or itmax too small" << endl;
 }
 Double_t v=sum*exp(-x+a*log(x)-gln);
 return v;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::GamCf(Double_t a,Double_t x)
{
// Computation of the incomplete gamma function P(a,x)
// via its continued fraction representation.
//
// The algorithm is based on the formulas and code as denoted in
// Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
//
//--- Nve 14-nov-1998 UU-SAP Utrecht
 
 Int_t itmax=100;      // Maximum number of iterations
 Double_t eps=3.e-7;    // Relative accuracy
 Double_t fpmin=1.e-30; // Smallest Double_t value allowed here
 
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
 
 Double_t gln=LnGamma(a);
 Double_t b=x+1.-a;
 Double_t c=1./fpmin;
 Double_t d=1./b;
 Double_t h=d;
 Double_t an,del;
 for (Int_t i=1; i<=itmax; i++)
 {
  an=double(-i)*(double(i)-a);
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
 Double_t v=exp(-x+a*log(x)-gln)*h;
 return (1.-v);
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::Erf(Double_t x)
{
// Computation of the error function erf(x).
//
//--- NvE 14-nov-1998 UU-SAP Utrecht
 
 return (1.-Erfc(x));
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::Erfc(Double_t x)
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
 const Double_t a1=-1.26551223,  a2=1.00002368,
                a3= 0.37409196,  a4=0.09678418,
                a5=-0.18628806,  a6=0.27886807,
                a7=-1.13520398,  a8=1.48851587,
                a9=-0.82215223, a10=0.17087277;
 
 Double_t v=1.; // The return value
 
 Double_t z=fabs(x);
 
 if (z <= 0.) return v; // erfc(0)=1
 
 Double_t t=1./(1.+0.5*z);
 
 v=t*exp((-z*z)
         +a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10)))))))));
 
 if (x < 0.) v=2.-v; // erfc(-x)=2-erfc(x)
 
 return v;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::Prob(Double_t chi2,Int_t ndf)
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
//  Double_t v=1.-Erf(sqrt(chi2)/sqrt(2.));
//  return v;
// }
 
// Gaussian approximation for large ndf
// This code may be activated in case the gamma function shows a problem
// Double_t q=sqrt(2.*chi2)-sqrt(double(2*ndf-1));
// if (n>30 && q>0.)
// {
//  Double_t v=0.5*(1.-Erf(q/sqrt(2.)));
//  return v;
// }
 
 // Evaluate the incomplete gamma function
 Double_t a=double(ndf)/2.;
 Double_t x=chi2/2.;
 return (1.-Gamma(a,x));
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselI0(Double_t x)
{
// Computation of the modified Bessel function I_0(x) for any real x.
//
// The algorithm is based on the article by Abramowitz and Stegun [1]
// as denoted in Numerical Recipes 2nd ed. on p. 230 (W.H.Press et al.).
//
// [1] M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
//     Applied Mathematics Series vol. 55 (1964), Washington.  
//
//--- NvE 12-mar-2000 UU-SAP Utrecht

 // Parameters of the polynomial approximation  
 const Double_t p1=1.0,          p2=3.5156229,    p3=3.0899424,
                p4=1.2067492,    p5=0.2659732,    p6=3.60768e-2,  p7=4.5813e-3;

 const Double_t q1= 0.39894228,  q2= 1.328592e-2, q3= 2.25319e-3,
                q4=-1.57565e-3,  q5= 9.16281e-3,  q6=-2.057706e-2,
                q7= 2.635537e-2, q8=-1.647633e-2, q9= 3.92377e-3; 

 Double_t ax=fabs(x);

 Double_t y=0,result=0;

 if (ax < 3.75)
 {
  y=pow(x/3.75,2);
  result=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
 }
 else
 {
  y=3.75/ax;
  result=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
 }

 return result;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselK0(Double_t x)
{
// Computation of the modified Bessel function K_0(x) for positive real x.
//
// The algorithm is based on the article by Abramowitz and Stegun [1]
// as denoted in Numerical Recipes 2nd ed. on p. 230 (W.H.Press et al.).
//
// [1] M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
//     Applied Mathematics Series vol. 55 (1964), Washington.  
//
//--- NvE 12-mar-2000 UU-SAP Utrecht

 // Parameters of the polynomial approximation  
 const Double_t p1=-0.57721566,  p2=0.42278420,   p3=0.23069756,
                p4= 3.488590e-2, p5=2.62698e-3,   p6=1.0750e-4,    p7=7.4e-5;

 const Double_t q1= 1.25331414,  q2=-7.832358e-2, q3= 2.189568e-2,
                q4=-1.062446e-2, q5= 5.87872e-3,  q6=-2.51540e-3,  q7=5.3208e-4;

 if (x <= 0)
 {
  cout << " *BesselK0* Invalid argument x = " << x << endl;
  return 0;
 }

 Double_t y=0,result=0;

 if (x <= 2)
 {
  y=x*x/4.;
  result=(-log(x/2.)*BesselI0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
 }
 else
 {
  y=2./x;
  result=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
 }

 return result;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselI1(Double_t x)
{
// Computation of the modified Bessel function I_1(x) for any real x.
//
// The algorithm is based on the article by Abramowitz and Stegun [1]
// as denoted in Numerical Recipes 2nd ed. on p. 230 (W.H.Press et al.).
//
// [1] M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
//     Applied Mathematics Series vol. 55 (1964), Washington.  
//
//--- NvE 12-mar-2000 UU-SAP Utrecht

 // Parameters of the polynomial approximation  
 const Double_t p1=0.5,          p2=0.87890594,   p3=0.51498869,
                p4=0.15084934,   p5=2.658733e-2,  p6=3.01532e-3,  p7=3.2411e-4;

 const Double_t q1= 0.39894228,  q2=-3.988024e-2, q3=-3.62018e-3,
                q4= 1.63801e-3,  q5=-1.031555e-2, q6= 2.282967e-2,
                q7=-2.895312e-2, q8= 1.787654e-2, q9=-4.20059e-3; 

 Double_t ax=fabs(x);

 Double_t y=0,result=0;

 if (ax < 3.75)
 {
  y=pow(x/3.75,2);
  result=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
 }
 else
 {
  y=3.75/ax;
  result=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
  if (x < 0) result=-result;
 }

 return result;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselK1(Double_t x)
{
// Computation of the modified Bessel function K_1(x) for positive real x.
//
// The algorithm is based on the article by Abramowitz and Stegun [1]
// as denoted in Numerical Recipes 2nd ed. on p. 230 (W.H.Press et al.).
//
// [1] M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
//     Applied Mathematics Series vol. 55 (1964), Washington.  
//
//--- NvE 12-mar-2000 UU-SAP Utrecht

 // Parameters of the polynomial approximation  
 const Double_t p1= 1.,          p2= 0.15443144,  p3=-0.67278579,
                p4=-0.18156897,  p5=-1.919402e-2, p6=-1.10404e-3,  p7=-4.686e-5;

 const Double_t q1= 1.25331414,  q2= 0.23498619,  q3=-3.655620e-2,
                q4= 1.504268e-2, q5=-7.80353e-3,  q6= 3.25614e-3,  q7=-6.8245e-4;

 if (x <= 0)
 {
  cout << " *BesselK1* Invalid argument x = " << x << endl;
  return 0;
 }

 Double_t y=0,result=0;

 if (x <= 2)
 {
  y=x*x/4.;
  result=(log(x/2.)*BesselI1(x))+(1./x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
 }
 else
 {
  y=2./x;
  result=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
 }

 return result;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselK(Int_t n,Double_t x)
{
// Computation of the Integer Order Modified Bessel function K_n(x)
// for n=0,1,2,... and positive real x.
//
// The algorithm uses the recurrence relation
//
//               K_n+1(x) = (2n/x)*K_n(x) + K_n-1(x) 
//
// as denoted in Numerical Recipes 2nd ed. on p. 232 (W.H.Press et al.).
//
//--- NvE 12-mar-2000 UU-SAP Utrecht

 if (x <= 0 || n < 0)
 {
  cout << " *BesselK* Invalid argument(s) (n,x) = (" << n << " , " << x << ")" << endl;
  return 0;
 }

 if (n==0) return BesselK0(x);

 if (n==1) return BesselK1(x);

 // Perform upward recurrence for all x
 Double_t tox=2./x;
 Double_t bkm=BesselK0(x);
 Double_t bk=BesselK1(x);
 Double_t bkp=0;
 for (Int_t j=1; j<n; j++)
 {
  bkp=bkm+double(j)*tox*bk;
  bkm=bk;
  bk=bkp;
 }

 return bk;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselI(Int_t n,Double_t x)
{
// Computation of the Integer Order Modified Bessel function I_n(x)
// for n=0,1,2,... and any real x.
//
// The algorithm uses the recurrence relation
//
//               I_n+1(x) = (-2n/x)*I_n(x) + I_n-1(x) 
//
// as denoted in Numerical Recipes 2nd ed. on p. 232 (W.H.Press et al.).
//
//--- NvE 12-mar-2000 UU-SAP Utrecht

 Int_t iacc=40; // Increase to enhance accuracy
 Double_t bigno=1.e10, bigni=1.e-10;

 if (n < 0)
 {
  cout << " *BesselI* Invalid argument (n,x) = (" << n << " , " << x << ")" << endl;
  return 0;
 }

 if (n==0) return BesselI0(x);

 if (n==1) return BesselI1(x);

 if (fabs(x) < 1.e-10) return 0;

 Double_t tox=2./fabs(x);
 Double_t bip=0,bim=0;
 Double_t bi=1;
 Double_t result=0;
 Int_t m=2*((n+int(sqrt(float(iacc*n))))); // Downward recurrence from even m
 for (Int_t j=m; j<=1; j--)
 {
  bim=bip+double(j)*tox*bi;
  bip=bi;
  bi=bim;
  if (fabs(bi) > bigno) // Renormalise to prevent overflows
  {
   result*=bigni;
   bi*=bigni;
   bip*=bigni;
  }
  if (j==n) result=bip;
 }

 result*=BesselI0(x)/bi; // Normalise with I0(x)
 if ((x < 0) && (n%2 == 1)) result=-result;

 return result;
}
///////////////////////////////////////////////////////////////////////////
