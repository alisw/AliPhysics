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

// $Id$

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
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliMath.h"
#include "Riostream.h"
 
ClassImp(AliMath) // Class implementation to enable ROOT I/O
 
AliMath::AliMath() : TObject()
{
// Default constructor
}
///////////////////////////////////////////////////////////////////////////
AliMath::~AliMath()
{
// Destructor
}
///////////////////////////////////////////////////////////////////////////
AliMath::AliMath(const AliMath& m) : TObject(m)
{
// Copy constructor
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::Gamma(Double_t z) const
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
Double_t AliMath::Gamma(Double_t a,Double_t x) const
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
Double_t AliMath::LnGamma(Double_t z) const
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
Double_t AliMath::GamSer(Double_t a,Double_t x) const
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
Double_t AliMath::GamCf(Double_t a,Double_t x) const
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
Double_t AliMath::Erf(Double_t x) const
{
// Computation of the error function erf(x).
//
//--- NvE 14-nov-1998 UU-SAP Utrecht
 
 return (1.-Erfc(x));
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::Erfc(Double_t x) const
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
 const Double_t ka1=-1.26551223,  ka2=1.00002368,
                ka3= 0.37409196,  ka4=0.09678418,
                ka5=-0.18628806,  ka6=0.27886807,
                ka7=-1.13520398,  ka8=1.48851587,
                ka9=-0.82215223, ka10=0.17087277;
 
 Double_t v=1.; // The return value
 
 Double_t z=fabs(x);
 
 if (z <= 0.) return v; // erfc(0)=1
 
 Double_t t=1./(1.+0.5*z);
 
 v=t*exp((-z*z)
   +ka1+t*(ka2+t*(ka3+t*(ka4+t*(ka5+t*(ka6+t*(ka7+t*(ka8+t*(ka9+t*ka10)))))))));
 
 if (x < 0.) v=2.-v; // erfc(-x)=2-erfc(x)
 
 return v;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::Prob(Double_t chi2,Int_t ndf,Int_t mode) const
{
// Computation of the probability for a certain Chi-squared (chi2)
// and number of degrees of freedom (ndf).
//
// A more clear and flexible facility is offered by Chi2Pvalue. 
//
// According to the value of the parameter "mode" various algorithms
// can be selected.
//
// mode = 0 : Calculations are based on the incomplete gamma function P(a,x),
//            where a=ndf/2 and x=chi2/2.
//
//        1 : Same as for mode=0. However, in case ndf=1 an exact expression
//            based on the error function Erf() is used.
//
//        2 : Same as for mode=0. However, in case ndf>30 a Gaussian approximation
//            is used instead of the gamma function.
//
// When invoked as Prob(chi2,ndf) the default mode=1 is used.
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

 Double_t v=-1.;

 switch (mode)
 {
  case 1: // Exact expression for ndf=1 as alternative for the gamma function
   if (ndf==1) v=1.-Erf(sqrt(chi2)/sqrt(2.));
   break;
 
  case 2: // Gaussian approximation for large ndf (i.e. ndf>30) as alternative for the gamma function
   if (ndf>30)
   {
    Double_t q=sqrt(2.*chi2)-sqrt(double(2*ndf-1));
    if (q>0.) v=0.5*(1.-Erf(q/sqrt(2.)));
   }
   break;
 }
 
 if (v<0.)
 {
  // Evaluate the incomplete gamma function
  Double_t a=double(ndf)/2.;
  Double_t x=chi2/2.;
  v=1.-Gamma(a,x);
 }

 return v;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselI0(Double_t x) const
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
 const Double_t kp1=1.0,          kp2=3.5156229,    kp3=3.0899424,
                kp4=1.2067492,    kp5=0.2659732,    kp6=3.60768e-2,  kp7=4.5813e-3;

 const Double_t kq1= 0.39894228,  kq2= 1.328592e-2, kq3= 2.25319e-3,
                kq4=-1.57565e-3,  kq5= 9.16281e-3,  kq6=-2.057706e-2,
                kq7= 2.635537e-2, kq8=-1.647633e-2, kq9= 3.92377e-3; 

 Double_t ax=fabs(x);

 Double_t y=0,result=0;

 if (ax < 3.75)
 {
  y=pow(x/3.75,2);
  result=kp1+y*(kp2+y*(kp3+y*(kp4+y*(kp5+y*(kp6+y*kp7)))));
 }
 else
 {
  y=3.75/ax;
  result=(exp(ax)/sqrt(ax))
         *(kq1+y*(kq2+y*(kq3+y*(kq4+y*(kq5+y*(kq6+y*(kq7+y*(kq8+y*kq9))))))));
 }

 return result;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselK0(Double_t x) const
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
 const Double_t kp1=-0.57721566,  kp2=0.42278420,   kp3=0.23069756,
                kp4= 3.488590e-2, kp5=2.62698e-3,   kp6=1.0750e-4,    kp7=7.4e-6;

 const Double_t kq1= 1.25331414,  kq2=-7.832358e-2, kq3= 2.189568e-2,
                kq4=-1.062446e-2, kq5= 5.87872e-3,  kq6=-2.51540e-3,  kq7=5.3208e-4;

 if (x <= 0)
 {
  cout << " *BesselK0* Invalid argument x = " << x << endl;
  return 0;
 }

 Double_t y=0,result=0;

 if (x <= 2)
 {
  y=x*x/4.;
  result=(-log(x/2.)*BesselI0(x))
         +(kp1+y*(kp2+y*(kp3+y*(kp4+y*(kp5+y*(kp6+y*kp7))))));
 }
 else
 {
  y=2./x;
  result=(exp(-x)/sqrt(x))
         *(kq1+y*(kq2+y*(kq3+y*(kq4+y*(kq5+y*(kq6+y*kq7))))));
 }

 return result;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselI1(Double_t x) const
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
 const Double_t kp1=0.5,          kp2=0.87890594,   kp3=0.51498869,
                kp4=0.15084934,   kp5=2.658733e-2,  kp6=3.01532e-3,  kp7=3.2411e-4;

 const Double_t kq1= 0.39894228,  kq2=-3.988024e-2, kq3=-3.62018e-3,
                kq4= 1.63801e-3,  kq5=-1.031555e-2, kq6= 2.282967e-2,
                kq7=-2.895312e-2, kq8= 1.787654e-2, kq9=-4.20059e-3; 

 Double_t ax=fabs(x);

 Double_t y=0,result=0;

 if (ax < 3.75)
 {
  y=pow(x/3.75,2);
  result=x*(kp1+y*(kp2+y*(kp3+y*(kp4+y*(kp5+y*(kp6+y*kp7))))));
 }
 else
 {
  y=3.75/ax;
  result=(exp(ax)/sqrt(ax))
         *(kq1+y*(kq2+y*(kq3+y*(kq4+y*(kq5+y*(kq6+y*(kq7+y*(kq8+y*kq9))))))));
  if (x < 0) result=-result;
 }

 return result;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselK1(Double_t x) const
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
 const Double_t kp1= 1.,          kp2= 0.15443144,  kp3=-0.67278579,
                kp4=-0.18156897,  kp5=-1.919402e-2, kp6=-1.10404e-3,  kp7=-4.686e-5;

 const Double_t kq1= 1.25331414,  kq2= 0.23498619,  kq3=-3.655620e-2,
                kq4= 1.504268e-2, kq5=-7.80353e-3,  kq6= 3.25614e-3,  kq7=-6.8245e-4;

 if (x <= 0)
 {
  cout << " *BesselK1* Invalid argument x = " << x << endl;
  return 0;
 }

 Double_t y=0,result=0;

 if (x <= 2)
 {
  y=x*x/4.;
  result=(log(x/2.)*BesselI1(x))
         +(1./x)*(kp1+y*(kp2+y*(kp3+y*(kp4+y*(kp5+y*(kp6+y*kp7))))));
 }
 else
 {
  y=2./x;
  result=(exp(-x)/sqrt(x))
         *(kq1+y*(kq2+y*(kq3+y*(kq4+y*(kq5+y*(kq6+y*kq7))))));
 }

 return result;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BesselK(Int_t n,Double_t x) const
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
Double_t AliMath::BesselI(Int_t n,Double_t x) const
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
TF1 AliMath::Chi2Dist(Int_t ndf) const
{
// Provide the Chi-squared distribution function corresponding to the
// specified ndf degrees of freedom.
//
// Details can be found in the excellent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <chi2>=ndf  Var(chi2)=2*ndf
 
 TF1 chi2dist("chi2dist","1./(TMath::Gamma([0]/2.)*pow(2,[0]/2.))*pow(x,[0]/2.-1.)*exp(-x/2.)");
 TString title="#chi^{2} distribution for ndf = ";
 title+=ndf; 
 chi2dist.SetTitle(title.Data());
 chi2dist.SetParName(0,"ndf");
 chi2dist.SetParameter(0,ndf);

 return chi2dist;
}
///////////////////////////////////////////////////////////////////////////
TF1 AliMath::StudentDist(Double_t ndf) const
{
// Provide the Student's T distribution function corresponding to the
// specified ndf degrees of freedom.
//
// In a frequentist approach, the Student's T distribution is particularly
// useful in making inferences about the mean of an underlying population
// based on the data from a random sample.
//
// In a Bayesian context it is used to characterise the posterior PDF
// for a particular state of information. 
//
// Note : ndf is not restricted to integer values
//
// Details can be found in the excellent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <T>=0  Var(T)=ndf/(ndf-2)
 
 TF1 tdist("tdist","(TMath::Gamma(([0]+1.)/2.)/(sqrt(pi*[0])*TMath::Gamma([0]/2.)))*pow(1.+x*x/[0],-([0]+1.)/2.)");
 TString title="Student's t distribution for ndf = ";
 title+=ndf; 
 tdist.SetTitle(title.Data());
 tdist.SetParName(0,"ndf");
 tdist.SetParameter(0,ndf);

 return tdist;
}
///////////////////////////////////////////////////////////////////////////
TF1 AliMath::FratioDist(Int_t ndf1,Int_t ndf2) const
{
// Provide the F (ratio) distribution function corresponding to the
// specified ndf1 and ndf2 degrees of freedom of the two samples.
//
// In a frequentist approach, the F (ratio) distribution is particularly useful
// in making inferences about the ratio of the variances of two underlying
// populations based on a measurement of the variance of a random sample taken
// from each one of the two populations.
// So the F test provides a means to investigate the degree of equality of
// two populations.
//
// Details can be found in the excellent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <F>=ndf2/(ndf2-2)  Var(F)=2*ndf2*ndf2*(ndf2+ndf1-2)/(ndf1*(ndf2-1)*(ndf2-1)*(ndf2-4))
 
 TF1 fdist("fdist",
 "(TMath::Gamma(([0]+[1])/2.)/(TMath::Gamma([0]/2.)*TMath::Gamma([1]/2.)))*pow([0]/[1],[0]/2.)*pow(x,([0]-2.)/2.)/pow(1.+x*[0]/[1],([0]+[1])/2.)");
 TString title="F (ratio) distribution for ndf1 = ";
 title+=ndf1;
 title+=" ndf2 = "; 
 title+=ndf2;
 fdist.SetTitle(title.Data());
 fdist.SetParName(0,"ndf1");
 fdist.SetParameter(0,ndf1);
 fdist.SetParName(1,"ndf2");
 fdist.SetParameter(1,ndf2);

 return fdist;
}
///////////////////////////////////////////////////////////////////////////
TF1 AliMath::BinomialDist(Int_t n,Double_t p) const
{
// Provide the Binomial distribution function corresponding to the
// specified number of trials n and probability p of success.
//
// p(k|n,p) = probability to obtain exactly k successes in n trials. 
//
// Details can be found in the excelent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <k>=n*p  Var(k)=n*p*(1-p)
 
 TF1 bindist("bindist","TMath::Binomial(int([0]),int(x))*pow([1],int(x))*pow(1.-[1],int([0])-int(x))");
 TString title="Binomial distribution for n = ";
 title+=n;
 title+=" p = "; 
 title+=p;
 bindist.SetTitle(title.Data());
 bindist.SetParName(0,"n");
 bindist.SetParameter(0,n);
 bindist.SetParName(1,"p");
 bindist.SetParameter(1,p);

 return bindist;
}
///////////////////////////////////////////////////////////////////////////
TF1 AliMath::NegBinomialDist(Int_t k,Double_t p) const
{
// Provide the Negative Binomial distribution function corresponding to the
// specified number of successes k and probability p of success.
//
// p(n|k,p) = probability to have reached k successes after n trials.
//
// Details can be found in the excelent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <n>=k*(1-p)/p  Var(n)=k*(1-p)/(p*p) 
 
 TF1 negbindist("negbindist","TMath::Binomial(int(x)-1,int([0])-1)*pow([1],int([0]))*pow(1.-[1],int(x)-int([0]))");
 TString title="Negative Binomial distribution for k = ";
 title+=k;
 title+=" p = "; 
 title+=p;
 negbindist.SetTitle(title.Data());
 negbindist.SetParName(0,"k");
 negbindist.SetParameter(0,k);
 negbindist.SetParName(1,"p");
 negbindist.SetParameter(1,p);

 return negbindist;
}
///////////////////////////////////////////////////////////////////////////
TF1 AliMath::PoissonDist(Double_t mu) const
{
// Provide the Poisson distribution function corresponding to the
// specified average rate (in time or space) mu of occurrences.
//
// p(n|mu) = probability for n occurrences in a certain time or space interval.
//
// Details can be found in the excelent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <n>=mu  Var(n)=mu
 
 TF1 poissdist("poissdist","exp(-[0])*pow([0],int(x))/TMath::Factorial(int(x))");
 TString title="Poisson distribution for mu = ";
 title+=mu;
 poissdist.SetTitle(title.Data());
 poissdist.SetParName(0,"mu");
 poissdist.SetParameter(0,mu);

 return poissdist;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::Chi2Pvalue(Double_t chi2,Int_t ndf,Int_t sides,Int_t sigma,Int_t mode) const
{
// Computation of the P-value for a certain specified Chi-squared (chi2) value 
// for a Chi-squared distribution with ndf degrees of freedom.
//
// The P-value for a certain Chi-squared value chi2 corresponds to the fraction of
// repeatedly drawn equivalent samples from a certain population, which is expected
// to yield a Chi-squared value exceeding (less than) the value chi2 for an
// upper (lower) tail test in case a certain hypothesis is true.
//
// Further details can be found in the excellent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <Chi2>=ndf  Var(Chi2)=2*ndf
// 
// With the "sides" parameter a one-sided or two-sided test can be selected
// using either the upper or lower tail contents.
// In case of automatic upper/lower selection the decision is made on basis
// of the location of the input chi2 value w.r.t. <Chi2> of the distribution. 
//
// sides =  1 : One-sided test using the upper tail contents
//          2 : Two-sided test using the upper tail contents
//         -1 : One-sided test using the lower tail contents
//         -2 : Two-sided test using the lower tail contents
//          0 : One-sided test using the auto-selected upper or lower tail contents
//          3 : Two-sided test using the auto-selected upper or lower tail contents
//
// The argument "sigma" allows for the following return values :
//
// sigma = 0 : P-value is returned as the above specified fraction
//         1 : The difference chi2-<Chi2> expressed in units of sigma
//             Note : This difference may be negative.
//  
// According to the value of the parameter "mode" various algorithms
// can be selected.
//
// mode = 0 : Calculations are based on the incomplete gamma function.
//
//        1 : Same as for mode=0. However, in case ndf=1 an exact expression
//            based on the error function Erf() is used.
//
//        2 : Same as for mode=0. However, in case ndf>30 a Gaussian approximation
//            is used instead of the gamma function.
//
// The default values are sides=0, sigma=0 and mode=1.
//
//--- NvE 21-may-2005 Utrecht University

 if (ndf<=0) return 0;

 Double_t mean=ndf;

 if (!sides) // Automatic one-sided test
 {
  sides=1;
  if (chi2<mean) sides=-1;
 }

 if (sides==3) // Automatic two-sided test
 {
  sides=2;
  if (chi2<mean) sides=-2;
 }

 Double_t val=0;
 if (sigma) // P-value in units of sigma
 {
  Double_t s=sqrt(double(2*ndf));
  val=(chi2-mean)/s;
 }
 else // P-value from tail contents
 {
  if (sides>0) // Upper tail
  { 
   val=Prob(chi2,ndf,mode);
  }
  else // Lower tail
  {
   val=1.-Prob(chi2,ndf,mode);
  }
 }

 if (abs(sides)==2) val=val*2.;

 return val;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::StudentPvalue(Double_t t,Double_t ndf,Int_t sides,Int_t sigma) const
{
// Computation of the P-value for a certain specified Student's t value 
// for a Student's T distribution with ndf degrees of freedom.
//
// In a frequentist approach, the Student's T distribution is particularly
// useful in making inferences about the mean of an underlying population
// based on the data from a random sample.
//
// The P-value for a certain t value corresponds to the fraction of
// repeatedly drawn equivalent samples from a certain population, which is expected
// to yield a T value exceeding (less than) the value t for an upper (lower)
// tail test in case a certain hypothesis is true.
//
// Further details can be found in the excellent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <T>=0  Var(T)=ndf/(ndf-2)
// 
// With the "sides" parameter a one-sided or two-sided test can be selected
// using either the upper or lower tail contents.
// In case of automatic upper/lower selection the decision is made on basis
// of the location of the input t value w.r.t. <T> of the distribution. 
//
// sides =  1 : One-sided test using the upper tail contents
//          2 : Two-sided test using the upper tail contents
//         -1 : One-sided test using the lower tail contents
//         -2 : Two-sided test using the lower tail contents
//          0 : One-sided test using the auto-selected upper or lower tail contents
//          3 : Two-sided test using the auto-selected upper or lower tail contents
//
// The argument "sigma" allows for the following return values :
//
// sigma = 0 : P-value is returned as the above specified fraction
//         1 : The difference t-<T> expressed in units of sigma
//             Note : This difference may be negative and sigma
//                    is only defined for ndf>2.
//  
// The default values are sides=0 and sigma=0.
//
//--- NvE 21-may-2005 Utrecht University

 if (ndf<=0) return 0;

 Double_t mean=0;

 if (!sides) // Automatic one-sided test
 {
  sides=1;
  if (t<mean) sides=-1;
 }

 if (sides==3) // Automatic two-sided test
 {
  sides=2;
  if (t<mean) sides=-2;
 }

 Double_t val=0;
 if (sigma) // P-value in units of sigma
 { 
  if (ndf>2) // Sigma is only defined for ndf>2
  {
   Double_t s=sqrt(ndf/(ndf-2.));
   val=t/s;
  }
 }
 else // P-value from tail contents
 {
  if (sides>0) // Upper tail
  {
   val=1.-TMath::StudentI(t,ndf);
  }
  else // Lower tail
  {
   val=TMath::StudentI(t,ndf);
  }
 }

 if (abs(sides)==2) val=val*2.;

 return val;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::FratioPvalue(Double_t f,Int_t ndf1,Int_t ndf2,Int_t sides,Int_t sigma) const
{
// Computation of the P-value for a certain specified F ratio f value 
// for an F (ratio) distribution with ndf1 and ndf2 degrees of freedom
// for the two samples X,Y respectively to be compared in the ratio X/Y.
//
// In a frequentist approach, the F (ratio) distribution is particularly useful
// in making inferences about the ratio of the variances of two underlying
// populations based on a measurement of the variance of a random sample taken
// from each one of the two populations.
// So the F test provides a means to investigate the degree of equality of
// two populations.
//
// The P-value for a certain f value corresponds to the fraction of
// repeatedly drawn equivalent samples from a certain population, which is expected
// to yield an F value exceeding (less than) the value f for an upper (lower)
// tail test in case a certain hypothesis is true.
//
// Further details can be found in the excellent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <F>=ndf2/(ndf2-2)  Var(F)=2*ndf2*ndf2*(ndf2+ndf1-2)/(ndf1*(ndf2-1)*(ndf2-1)*(ndf2-4))
// 
// With the "sides" parameter a one-sided or two-sided test can be selected
// using either the upper or lower tail contents.
// In case of automatic upper/lower selection the decision is made on basis
// of the location of the input f value w.r.t. <F> of the distribution. 
//
// sides =  1 : One-sided test using the upper tail contents
//          2 : Two-sided test using the upper tail contents
//         -1 : One-sided test using the lower tail contents
//         -2 : Two-sided test using the lower tail contents
//          0 : One-sided test using the auto-selected upper or lower tail contents
//          3 : Two-sided test using the auto-selected upper or lower tail contents
//
// The argument "sigma" allows for the following return values :
//
// sigma = 0 : P-value is returned as the above specified fraction
//         1 : The difference f-<F> expressed in units of sigma
//             Note : This difference may be negative and sigma
//                    is only defined for ndf2>4.
//  
// The default values are sides=0 and sigma=0.
//
//--- NvE 21-may-2005 Utrecht University

 if (ndf1<=0 || ndf2<=0 || f<=0) return 0;

 Double_t mean=double(ndf2/(ndf2-2));

 if (!sides) // Automatic one-sided test
 {
  sides=1;
  if (f<mean) sides=-1;
 }

 if (sides==3) // Automatic two-sided test
 {
  sides=2;
  if (f<mean) sides=-2;
 }

 Double_t val=0;
 if (sigma) // P-value in units of sigma
 { 
  // Sigma is only defined for ndf2>4
  if (ndf2>4)
  {
   Double_t s=sqrt(double(ndf2*ndf2*(2*ndf2+2*ndf1-4))/double(ndf1*pow(double(ndf2-1),2)*(ndf2-4)));
   val=(f-mean)/s;
  }
 }
 else // P-value from tail contents 
 {
  if (sides>0) // Upper tail
  {
   val=1.-TMath::FDistI(f,ndf1,ndf2);
  }
  else // Lower tail
  {
   val=TMath::FDistI(f,ndf1,ndf2);
  }
 }

 if (abs(sides)==2) val=val*2.;

 return val;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::BinomialPvalue(Int_t k,Int_t n,Double_t p,Int_t sides,Int_t sigma,Int_t mode) const
{
// Computation of the P-value for a certain specified number of successes k
// for a Binomial distribution with n trials and success probability p.
//
// The P-value for a certain number of successes k corresponds to the fraction of
// repeatedly drawn equivalent samples from a certain population, which is expected
// to yield a number of successes exceeding (less than) the value k for an
// upper (lower) tail test in case a certain hypothesis is true.
//
// Further details can be found in the excellent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <K>=n*p  Var(K)=n*p*(1-p)
// 
// With the "sides" parameter a one-sided or two-sided test can be selected
// using either the upper or lower tail contents.
// In case of automatic upper/lower selection the decision is made on basis
// of the location of the input k value w.r.t. <K> of the distribution. 
//
// sides =  1 : One-sided test using the upper tail contents
//          2 : Two-sided test using the upper tail contents
//         -1 : One-sided test using the lower tail contents
//         -2 : Two-sided test using the lower tail contents
//          0 : One-sided test using the auto-selected upper or lower tail contents
//          3 : Two-sided test using the auto-selected upper or lower tail contents
//
// The argument "sigma" allows for the following return values :
//
// sigma = 0 : P-value is returned as the above specified fraction
//         1 : The difference k-<K> expressed in units of sigma
//             Note : This difference may be negative.
//
// mode = 0 : Incomplete Beta function will be used to calculate the tail content.
//        1 : Straightforward summation of the Binomial terms is used.
//
// The Incomplete Beta function in general provides the most accurate values.
//
// The default values are sides=0, sigma=0 and mode=0.
//
//--- NvE 24-may-2005 Utrecht University

 Double_t mean=double(n)*p;

 if (!sides) // Automatic one-sided test
 {
  sides=1;
  if (k<mean) sides=-1;
 }

 if (sides==3) // Automatic two-sided test
 {
  sides=2;
  if (k<mean) sides=-2;
 }

 Double_t val=0;

 if (sigma) // P-value in units of sigma
 {
  Double_t s=sqrt(double(n)*p*(1.-p));
  val=(double(k)-mean)/s;
 }
 else // P-value from tail contents
 {
  if (sides>0)
  {
   if (!mode) // Use Incomplete Beta function
   {
    val=TMath::BetaIncomplete(p,k+1,n-k);
   }
   else // Use straightforward summation
   {
    for (Int_t i=k+1; i<=n; i++)
    {
     val+=TMath::Binomial(n,i)*pow(p,i)*pow(1.-p,n-i);
    }
   }
  }
  else
  {
   if (!mode) // Use Incomplete Beta function
   {
    val=1.-TMath::BetaIncomplete(p,k+1,n-k);
   }
   else
   {
    for (Int_t j=0; j<=k; j++)
    {
     val+=TMath::Binomial(n,j)*pow(p,j)*pow(1.-p,n-j);
    }
   }
  }
 }

 if (abs(sides)==2) val=val*2.;

 return val;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::PoissonPvalue(Int_t k,Double_t mu,Int_t sides,Int_t sigma) const
{
// Computation of the P-value for a certain specified number of occurrences k
// for a Poisson distribution with a specified average rate (in time or space)
// mu of occurrences.
//
// The P-value for a certain number of occurrences k corresponds to the fraction of
// repeatedly drawn equivalent samples from a certain population, which is expected
// to yield a number of occurrences exceeding (less than) the value k for an
// upper (lower) tail test in case a certain hypothesis is true.
//
// Further details can be found in the excellent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <K>=mu  Var(K)=mu
// 
// With the "sides" parameter a one-sided or two-sided test can be selected
// using either the upper or lower tail contents.
// In case of automatic upper/lower selection the decision is made on basis
// of the location of the input k value w.r.t. <K> of the distribution. 
//
// sides =  1 : One-sided test using the upper tail contents
//          2 : Two-sided test using the upper tail contents
//         -1 : One-sided test using the lower tail contents
//         -2 : Two-sided test using the lower tail contents
//          0 : One-sided test using the auto-selected upper or lower tail contents
//          3 : Two-sided test using the auto-selected upper or lower tail contents
//
// The argument "sigma" allows for the following return values :
//
// sigma = 0 : P-value is returned as the above specified fraction
//         1 : The difference k-<K> expressed in units of sigma
//             Note : This difference may be negative.
//
// The default values are sides=0 and sigma=0.
//
// Note : The tail contents are given by the incomplete Gamma function P(a,x).
//        Lower tail contents = 1-P(k,mu)
//        Upper tail contents = P(k,mu)
//
//--- NvE 24-may-2005 Utrecht University

 Double_t mean=mu;

 if (!sides) // Automatic one-sided test
 {
  sides=1;
  if (k<mean) sides=-1;
 }

 if (sides==3) // Automatic two-sided test
 {
  sides=2;
  if (k<mean) sides=-2;
 }

 Double_t val=0;

 if (sigma) // P-value in units of sigma
 {
  Double_t s=sqrt(mu);
  val=(double(k)-mean)/s;
 }
 else // P-value from tail contents
 {
  if (sides>0) // Upper tail
  {
   val=Gamma(k,mu);
  }
  else // Lower tail
  {
   val=1.-Gamma(k,mu);
  }
 }

 if (abs(sides)==2) val=val*2.;

 return val;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::NegBinomialPvalue(Int_t n,Int_t k,Double_t p,Int_t sides,Int_t sigma,Int_t mode) const
{
// Computation of the P-value for a certain specified number of trials n
// for a Negative Binomial distribution where exactly k successes are to
// be reached which have each a probability p.
//
// The P-value for a certain number of trials n corresponds to the fraction of
// repeatedly drawn equivalent samples from a certain population, which is expected
// to yield a number of trials exceeding (less than) the value n for an
// upper (lower) tail test in case a certain hypothesis is true.
//
// Further details can be found in the excellent textbook of Phil Gregory
// "Bayesian Logical Data Analysis for the Physical Sciences".
//
// Note : <N>=k*(1-p)/p  Var(N)=k*(1-p)/(p*p) 
// 
// With the "sides" parameter a one-sided or two-sided test can be selected
// using either the upper or lower tail contents.
// In case of automatic upper/lower selection the decision is made on basis
// of the location of the input n value w.r.t. <N> of the distribution. 
//
// sides =  1 : One-sided test using the upper tail contents
//          2 : Two-sided test using the upper tail contents
//         -1 : One-sided test using the lower tail contents
//         -2 : Two-sided test using the lower tail contents
//          0 : One-sided test using the auto-selected upper or lower tail contents
//          3 : Two-sided test using the auto-selected upper or lower tail contents
//
// The argument "sigma" allows for the following return values :
//
// sigma = 0 : P-value is returned as the above specified fraction
//         1 : The difference n-<N> expressed in units of sigma
//             Note : This difference may be negative.
//
// mode = 0 : Incomplete Beta function will be used to calculate the tail content.
//        1 : Straightforward summation of the Negative Binomial terms is used.
//
// The Incomplete Beta function in general provides the most accurate values.
//
// The default values are sides=0, sigma=0 and mode=0.
//
//--- NvE 24-may-2005 Utrecht University

 Double_t mean=double(k)*(1.-p)/p;

 if (!sides) // Automatic one-sided test
 {
  sides=1;
  if (n<mean) sides=-1;
 }

 if (sides==3) // Automatic two-sided test
 {
  sides=2;
  if (n<mean) sides=-2;
 }

 Double_t val=0;

 if (sigma) // P-value in units of sigma
 {
  Double_t s=sqrt(double(k)*(1.-p)/(p*p));
  val=(double(n)-mean)/s;
 }
 else // P-value from tail contents
 {
  if (sides>0) // Upper tail
  {
   if (!mode) // Use Incomplete Beta function
   {
    val=1.-TMath::BetaIncomplete(p,k,n-k);
   }
   else // Use straightforward summation
   {
    for (Int_t i=1; i<n; i++)
    {
     val+=TMath::Binomial(i-1,k-1)*pow(p,k)*pow(1.-p,i-k);
    }
    val=1.-val;
   }
  }
  else // Lower tail
  {
   if (!mode) // Use Incomplete Beta function
   {
    val=TMath::BetaIncomplete(p,k,n-k);
   }
   else
   {
    for (Int_t j=1; j<n; j++)
    {
     val+=TMath::Binomial(j-1,k-1)*pow(p,k)*pow(1.-p,j-k);
    }
   }
  }
 }

 if (abs(sides)==2) val=val*2.;

 return val;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::Nfac(Int_t n,Int_t mode) const
{
// Compute n!.
// The algorithm can be selected by the "mode" input argument.
//
// mode : 0 ==> Calculation by means of straightforward multiplication
//      : 1 ==> Calculation by means of Stirling's approximation
//      : 2 ==> Calculation by means of n!=Gamma(n+1)
//
// For large n the calculation modes 1 and 2 will in general be faster.
// By default mode=0 is used.
// For n<0 the value 0 will be returned.
//
// Note : Because of Double_t value overflow the maximum value is n=170.
//
//--- NvE 20-jan-2007 Utrecht University

 if (n<0) return 0;
 if (n==0) return 1;

 Double_t twopi=2.*acos(-1.);
 Double_t z=0;
 Double_t nfac=1;
 Int_t i=n;
 
 switch (mode)
 {
  case 0: // Straightforward multiplication
   while (i>1)
   {
    nfac*=Double_t(i);
    i--;
   }
   break;

  case 1: // Stirling's approximation 
   z=n;
   nfac=sqrt(twopi)*pow(z,z+0.5)*exp(-z)*(1.+1./(12.*z));
   break;

  case 2: // Use of Gamma(n+1)
   z=n+1;
   nfac=Gamma(z);
   break;

  default:
   nfac=0;
   break;
 }

 return nfac;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::LnNfac(Int_t n,Int_t mode) const
{
// Compute ln(n!).
// The algorithm can be selected by the "mode" input argument.
//
// mode : 0 ==>  Calculation via evaluation of n! followed by taking ln(n!)
//      : 1 ==>  Calculation via Stirling's approximation ln(n!)=0.5*ln(2*pi)+(n+0.5)*ln(n)-n+1/(12*n)
//      : 2 ==>  Calculation by means of ln(n!)=LnGamma(n+1)
//
// Note : Because of Double_t value overflow the maximum value is n=170 for mode=0.
//
// For mode=2 rather accurate results are obtained for both small and large n.
// By default mode=2 is used.
// For n<1 the value 0 will be returned.
//
//--- NvE 20-jan-2007 Utrecht University

 if (n<=1) return 0;

 Double_t twopi=2.*acos(-1.);
 Double_t z=0;
 Double_t lognfac=0;
 
 switch (mode)
 {
  case 0: // Straightforward ln(n!)
   z=Nfac(n);
   lognfac=log(z);
   break;

  case 1: // Stirling's approximation 
   z=n;
   lognfac=0.5*log(twopi)+(z+0.5)*log(z)-z+1./(12.*z);
   break;

  case 2: // Use of LnGamma(n+1)
   z=n+1;
   lognfac=LnGamma(z);
   break;

  default:
   lognfac=0;
   break;
 }

 return lognfac;
}
///////////////////////////////////////////////////////////////////////////
Double_t AliMath::LogNfac(Int_t n,Int_t mode) const
{
// Compute log_10(n!).
// First ln(n!) is evaluated via invokation of LnNfac(n,mode).
// Then the algorithm log_10(z)=ln(z)*log_10(e) is used.
//
//--- NvE 20-jan-2007 Utrecht University

 Double_t e=exp(1.);

 Double_t val=LnNfac(n,mode);
 val*=log10(e);

 return val;
}
///////////////////////////////////////////////////////////////////////////
