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
// Class AliRandom
// Generate universal random numbers on all common machines.
// Available distributions : Uniform, Gaussian, Poisson and
//                           User defined function
//
// Features :
// ----------
// 1) Period = 2**144
// 2) Same sequence of 24-bit real numbers on all common machines
//
// Reference :
// -----------
// G.Marsaglia and A.Zaman, FSU-SCRI-87-50, Florida State University, 1987.
//
// Coding example :
// ----------------
//
// Float_t rndm;          // Variable to hold a single random number
// const Int_t n=1000;
// Float_t rvec[n];       // Vector to hold n random numbers
//
// AliRandom r;           // Create a Random object with default sequence
//
// rndm=r.Uniform();      // Provide a uniform random number in <0,1>
// Float_t a=3.;
// Float_t b=5.;
// rndm=r.Uniform(a,b);   // Provide a uniform random number in <a,b>
// r.Uniform(rvec,n);     // Provide n uniform randoms in <0,1> in rvec
// r.Uniform(rvec,n,a,b); // Provide n uniform randoms in <a,b> in rvec
//
// rndm=r.Gauss();             // Provide a Gaussian random number with
//                             // mean=0 and sigma=1
// Float_t mean=25.;
// Float_t sigma=5.;
// rndm=r.Gauss(mean,sigma);   // Provide a Gaussian random number
//                             // with specified mean and sigma
// r.Gauss(rvec,n);            // n Gaussian randoms mean=0 sigma=1
// r.Gauss(rvec,n,mean,sigma); // n Gaussian randoms with specified
//                             //  mean and sigma
//
// rndm=r.Poisson(mean);  // Provide a Poisson random number with
//                        // specified mean
// r.Poisson(rvec,nmean); // n Poisson randoms with specified mean
//
// Int_t seed=1837724
// AliRandom p(seed);        // Create a Random object with specified seed.
//                           // The sequence is started from scratch.
// Int_t cnt1=25;
// Int_t cnt2=8;
// AliRandom q(seed,cnt1,cnt2); // Create a Random object with specified seed
//                              // The sequence is started at the location
//                              // denoted by the counters cnt1 and cnt2.
//
// q.Info();     // Print the current seed, cnt1 and cnt2 values.
// q.GetSeed();  // Provide the current seed value.
// q.GetCnt1();  // Provide the current cnt1 value.
// q.GetCnt2();  // Provide the current cnt2 value.
//
// Float_t udist(Float_t x) // A user defined distribution
// {
//  return x*x-4.*x;
// }
//
// Int_t nbins=100;
// q.SetUser(a,b,nbins,udist); // Initialise generator for udist distribution
// q.User(); // Provide a random number according to the udist distribution
// q.User(rvec,n); // Provide n randoms according to the udist distribution
//
// Float_t* x=new Float_t[nbins];
// Float_t* y=new Float_t[nbins];
//
// ... code to fill x[] and y[] ..
//
// AliRandom s;
// s.SetUser(x,y,nbins); // Initialise generator for (x[i],y[i]) distribution
// s.User(); // Provide a random number according to the user distribution
// s.User(rvec,n); // Provide n randoms according to the user distribution
//
// Notes :
// -------
// 1) Allowed seed values : 0 <= seed <= 921350143
//    Default seed = 53310452
// 2) To ensure a unique sequence for each run, one can automatically
//    construct a seed value by e.g. using the date and time.
// 3) Using the rvec facility saves a lot of CPU time for large n values.
//
//--- Author: Nick van Eijndhoven 11-oct-1997 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliRandom.h"
 
ClassImp(AliRandom) // Class implementation to enable ROOT I/O
 
AliRandom::AliRandom()
{
// Creation of an AliRandom object and default initialisation.
//
// A seed is used to create the initial u[97] table.
// This seed is converted into four startup parameters i j k and l
// (see member function "unpack").
//
// Suggested test values : i=12 j=34 k=56 l=78 (see article)
// which corresponds to  : seed = 53310452
 
 Int_t seed=53310452; // Default seed
 Start(seed,0,0);     // Start the sequence for this seed from scratch
}
///////////////////////////////////////////////////////////////////////////
AliRandom::AliRandom(Int_t seed)
{
// Creation of an AliRandom object and user defined initialisation
 
 Start(seed,0,0); // Start the sequence for this seed from scratch
}
///////////////////////////////////////////////////////////////////////////
AliRandom::AliRandom(Int_t seed,Int_t cnt1,Int_t cnt2)
{
// Creation of an AliRandom object and user defined initialisation
//
// seed is the seed to create the initial u[97] table.
// The range of the seed is : 0 <= seed <= 921350143
//
// cnt1 and cnt2 are the parameters for the counting system
// to enable a start of the sequence at a certain point.
// The current values of seed, cnt1 and cnt2 can be obtained
// via the member functions "GetSeed", "GetCnt1" and "GetCnt2" resp.
// To start from scratch one should select : cnt1=0 and cnt2=0
 
 Start(seed,cnt1,cnt2); // Start the sequence from a user defined point
}
///////////////////////////////////////////////////////////////////////////
AliRandom::~AliRandom()
{
// Destructor to delete memory allocated for the area function arrays
 if (fXa) delete [] fXa;
 fXa=0;
 if (fYa) delete [] fYa;
 fYa=0;
 if (fIbins) delete [] fIbins;
 fIbins=0;
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::Start(Int_t seed,Int_t cnt1,Int_t cnt2)
{
// Start a certain sequence from scratch or from a user defined point
//
// The algorithm to start from scratch is based on the routine RSTART
// as described in the report by G.Marsaglia and A.Zaman
// (FSU-SCRI-87-50 Florida State University 1987).
//
// seed is the seed to create the initial u[97] table.
// This seed is converted into four startup parameters i j k and l
// (see the member function "unpack").
//
// The range of the seed is : 0 <= seed <= 921350143
//
// Suggested test values : i=12 j=34 k=56 l=78 (see article)
// which corresponds to  : seed = 53310452
//
// cnt1 and cnt2 are the parameters for the counting system
// to enable a start of the sequence at a certain point.
// The current values of seed, cnt1 and cnt2 can be obtained
// via the member functions "GetSeed", "GetCnt1" and "GetCnt2" resp.
// To start from scratch one should select : cnt1=0 and cnt2=0
 
// Reset the area function
 fNa=0;
 fXa=0;
 fYa=0;
 fIbins=0;
 
// Clipping parameter to prevent overflow of the counting system
 fClip=1000000;
 
// Set the lags for the Fibonacci sequence of the first part
// The sequence is set to F(97,33,*) (see article)
 fI=97;
 fJ=33;
 
// Unpack the seed value and determine i, j, k and l
 fSeed=seed;
 Int_t i,j,k,l;
 Unpack(seed,i,j,k,l);
 
// Reset counters
 fCnt1=0;
 fCnt2=0;
 
// Fill the starting table for the first part of the combination
 Float_t s,t;
 Int_t m;
 for (Int_t ii=0; ii<97; ii++)
 {
  s=0.;
  t=0.5;
 
  for (Int_t jj=1; jj<=24; jj++)
  {
   m=(((i*j)%179)*k)%179;
   i=j;
   j=k;
   k=m;
   l=((53*l)+1)%169;
   if ((l*m)%64 >= 32) s+=t;
   t=0.5*t;
  }
  fU[ii]=s;
 }
 
// Initialise the second part of the combination
 fC=362436./16777216.;
 fCd=7654321./16777216.;
 fCm=16777213./16777216.;
 
// Generate random numbers upto the user selected starting point
// on basis of the counting system
 if (cnt1 > 0) Uniform(cnt1);
 if (cnt2 > 0)
 {
  for (Int_t n=1; n<=cnt2; n++)
  {
   Uniform(fClip);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::Unpack(Int_t seed,Int_t& i,Int_t& j,Int_t& k,Int_t& l)
{
// Unpack the seed into the four startup parameters i,j,k and l
//
// The range of the seed is : 0 <= seed <= 921350143
//
// From the article :
// The i,j and k values may be chosen in the interval : 1 <= n <= 178
// The l value may be in the interval : 0 <= l <= 168
//
// Taking into account the period of the 3-lagged Fibonacci sequence
// The following "bad" combinations have to be ruled out :
//
//   i    j    k    l   period
//   1    1    1    X      1
// 178    1    1    X      4
//   1  178    1    X      2
//   1    1  178    X      4
// 178  178    1    X      4
// 178    1  178    X      2
//   1  178  178    X      4
// 178  178  178    X      1
//
// To rule these "bad" combinations out all together, we choose
// the following allowed ranges :
// The i,j and k values may be chosen in the interval : 2 <= n <= 177
// The l value may be in the interval : 0 <= l <= 168
//
// and use the formula :
// seed = (i-2)*176*176*169 + (j-2)*176*169 + (k-2)*169 + l
 
 if ((seed >= 0) && (seed <= 921350143)) // Check seed value
 {
  Int_t idum=seed;
  Int_t imin2=idum/(176*176*169);
  idum=idum%(176*176*169);
  Int_t jmin2=idum/(176*169);
  idum=idum%(176*169);
  Int_t kmin2=idum/169;
 
  i=imin2+2;
  j=jmin2+2;
  k=kmin2+2;
  l=seed%169;
 }
 else
 {
  cout << " *AliRandom::unpack()* Unallowed seed value encountered."
       << " seed = " << seed << endl;
  cout << " Seed will be set to default value." << endl;
 
  seed=53310452;   // Default seed
  Start(seed,0,0); // Start the sequence for this seed from scratch
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliRandom::GetSeed()
{
// Provide the current seed value
 return fSeed;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliRandom::GetCnt1()
{
// Provide the current value of the counter cnt1
 return fCnt1;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliRandom::GetCnt2()
{
// Provide the current value of the counter cnt2
 return fCnt2;
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::Info()
{
// Print the current seed, cnt1 and cnt2 values
 cout << " *Random* seed = " << fSeed
      << " cnt1 = " << fCnt1 << " cnt2 = " << fCnt2 << endl;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliRandom::Uniform()
{
// Generate uniform random numbers in the interval <0,1>
//
// The algorithm is based on lagged Fibonacci sequences (first part)
// combined with a congruential method (second part)
// as described in the report by G.Marsaglia and A.Zaman
// (FSU-SCRI-87-50 Florida State University 1987).
//
// Features :
// 1) Period = 2**144
// 2) Same sequence of 24-bit real numbers on all common machines
 
// First part of the combination : F(97,33,*) (see article for explanation)
 Float_t unirnu=fU[fI-1]-fU[fJ-1];
 if (unirnu < 0) unirnu+=1.;
 fU[fI-1]=unirnu;
 fI-=1;
 if (fI == 0) fI=97;
 fJ-=1;
 if (fJ == 0) fJ=97;
 
// Second part of the combination (see article for explanation)
 fC-=fCd;
 if (fC < 0.) fC+=fCm;
 unirnu-=fC;
 if (unirnu < 0.) unirnu+=1.;
 
// Update the counting system to enable sequence continuation
// at an arbitrary starting position.
// Two counters have been introduced to avoid overflow
// fCnt1 : Counter which goes up to fClip
//         and is reset when fClip is reached
// fCnt2 : Counts the number of times fClip has been reached
 fCnt1+=1;
 if (fCnt1 >= fClip)
 {
  fCnt1=0;
  fCnt2+=1;
 }
 
 if (unirnu <= 0.) unirnu=Uniform(); // Exclude 0. from the range
 
 return unirnu;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliRandom::Uniform(Float_t a,Float_t b)
{
// Generate uniform random numbers in the interval <a,b>
 Float_t rmin=a;
 if (a > b) rmin=b;
 
 Float_t rndm=Uniform();
 rndm=rmin+fabs(rndm*(a-b));
 
 return rndm;
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::Uniform(Float_t* vec,Int_t n,Float_t a,Float_t b)
{
// Generate a vector of uniform random numbers in the interval <a,b>
// This saves lots of (member)function calls in case many random
// numbers are needed at once.
//
// n = The number of random numbers to be generated
//
// The algorithm is based on lagged Fibonacci sequences (first part)
// combined with a congruential method (second part)
// as described in the report by G.Marsaglia and A.Zaman
// (FSU-SCRI-87-50 Florida State University 1987).
//
// Features :
// 1) Period = 2**144
// 2) Same sequence of 24-bit real numbers on all common machines
 
// Determine the minimum of a and b
 Float_t rmin=a;
 if (a > b) rmin=b;
 
// First generate random numbers within <0,1>
 if (n > 0) // Check n value
 {
  for (Int_t jvec=0; jvec<n; jvec++)
  {
   // First part of the combination : F(97,33,*)
   Float_t unirnu=fU[fI-1]-fU[fJ-1];
   if (unirnu < 0) unirnu+=1.;
   fU[fI-1]=unirnu;
   fI-=1;
   if (fI == 0) fI=97;
   fJ-=1;
   if (fJ == 0) fJ=97;
 
   // Second part of the combination
   fC-=fCd;
   if (fC < 0.) fC+=fCm;
   unirnu-=fC;
   if (unirnu < 0.) unirnu+=1.;
 
   // Update the counting system to enable sequence continuation
   // at an arbitrary starting position.
   // Two counters have been introduced to avoid overflow
   // fCnt1 : Counter which goes up to fClip
   //         and is reset when fClip is reached
   // fCnt2 : Counts the number of times fClip has been reached
   fCnt1+=1;
   if (fCnt1 >= fClip)
   {
    fCnt1=0;
    fCnt2+=1;
   }
 
   if (unirnu <= 0.) unirnu=Uniform(); // Exclude 0. from the range
 
   // Fill the vector within the selected range
   vec[jvec]=rmin+fabs(unirnu*(a-b));
  }
 }
 else
 {
  cout << " *AliRandom::Uniform* Invalid value n = " << n << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::Uniform(Float_t* vec,Int_t n)
{
// Generate a vector of uniform random numbers in the interval <0,1>
// This saves lots of (member)function calls in case many random
// numbers are needed at once.
//
// n = The number of random numbers to be generated
 
 Uniform(vec,n,0.,1.);
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::Uniform(Int_t n)
{
// Generate n uniform random numbers in in one go.
// This saves lots of (member)function calls in case one needs to skip
// to a certain point in a sequence.
//
// n = The number of random numbers to be generated
//
// Note : No check is made here to exclude 0 from the range.
//        It's only the number of generated randoms that counts
//
// The algorithm is based on lagged Fibonacci sequences (first part)
// combined with a congruential method (second part)
// as described in the report by G.Marsaglia and A.Zaman
// (FSU-SCRI-87-50 Florida State University 1987).
//
// Features :
// 1) Period = 2**144
// 2) Same sequence of 24-bit real numbers on all common machines
 
 if (n > 0) // Check n value
 {
  for (Int_t jvec=0; jvec<n; jvec++)
  {
   // First part of the combination : F(97,33,*)
   Float_t unirnu=fU[fI-1]-fU[fJ-1];
   if (unirnu < 0) unirnu+=1.;
   fU[fI-1]=unirnu;
   fI-=1;
   if (fI == 0) fI=97;
   fJ-=1;
   if (fJ == 0) fJ=97;
 
   // Second part of the combination
   fC-=fCd;
   if (fC < 0.) fC+=fCm;
   unirnu-=fC;
   if (unirnu < 0.) unirnu+=1.;
 
   // Update the counting system to enable sequence continuation
   // at an arbitrary starting position.
   // Two counters have been introduced to avoid overflow
   // fCnt1 : Counter which goes up to fClip
   //         and is reset when fClip is reached
   // fCnt2 : Counts the number of times fClip has been reached
   fCnt1+=1;
   if (fCnt1 >= fClip)
   {
    fCnt1=0;
    fCnt2+=1;
   }
  }
 }
 else
 {
  cout << " *AliRandom::Uniform* Invalid value n = " << n << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliRandom::Gauss(Float_t mean,Float_t sigma)
{
// Generate gaussian distributed random numbers with certain mean and sigma
//
// Method :
// P(x) = The gaussian distribution function
// --> ln(P) provides an expression for (x-xmean)**2 from which
//     the following expression for x can be obtained
//
//           x = xmean +/- sigma * sqrt(-2*ln(q))
//
// in which q is an expression in terms of pi, sigma and p and lies within
// the interval <0,1>.
//
// The gaussian distributed x values are obtained as follows :
//
// 1) Two uniform random numbers q1 and q2 in <0,1> are generated.
// 2) q1 is in fact a uniform generated value for P which is substituted
//    directly in the formula above.
// 3) The value of q2 determines whether we use the + or - sign.
 
// Generate the two uniform random numbers q1 and q2 in <0,1>
 Float_t q1,q2;
 q1=Uniform();
 q2=Uniform();
 
// Use q1 and q2 to get the gaussian distributed random number
 Float_t pi=acos(-1.);
 Float_t unirng=mean+cos(2.*pi*q2)*sigma*sqrt(-2.*log(q1));
 
 return unirng;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliRandom::Gauss()
{
// Generate gaussian distributed random numbers with mean=0 and sigma=1
 
 return Gauss(0.,1.);
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::Gauss(Float_t* vec,Int_t n,Float_t mean,Float_t sigma)
{
// Generate a vector of gaussian random numbers with certain mean and sigma
// This saves lots of (member)function calls in case many random
// numbers are needed at once.
//
// n = The number of random numbers to be generated
 
 if (n > 0) // Check n value
 {
  // The vector to hold the q1 and q2 values.
  // Two subsequent q[] values are used for q1 and q2
  // in order to obtain identical random numbers in the vector
  // as when generating n single ones.
  Int_t m=2*n;
  Float_t* q=new Float_t[m];
  Uniform(q,m);
 
  // Fill the vector with gaussian random numbers
  Float_t pi=acos(-1.);
  Float_t q1,q2;
  for (Int_t jvec=0; jvec<n; jvec++)
  {
   q1=q[jvec*2];     // use two subsequent q[] values
   q2=q[(jvec*2)+1];
   vec[jvec]=mean+cos(2.*pi*q2)*sigma*sqrt(-2.*log(q1));
  }
  delete [] q;
 }
 else
 {
  cout << " *AliRandom::Gauss* Invalid value n = " << n << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::Gauss(Float_t* vec,Int_t n)
{
// Generate a vector of gaussian random numbers with mean=0 and sigma=1
// This saves lots of (member)function calls in case many random
// numbers are needed at once.
//
// n = The number of random numbers to be generated
 
 Gauss(vec,n,0.,1.);
}
///////////////////////////////////////////////////////////////////////////
Float_t AliRandom::Poisson(Float_t mean)
{
// Generate Poisson distributed random numbers with certain mean
//
// Method :
//
// P(n) = exp(-mean)*mean**n/n! is the Poisson distribution function
//
// with : n  = 0,1,2,3,...   and  mean > 0
//
// To generate the distribution, the "sum trick" is used as mentioned
// in "Formulae and Methods in Experimental data Evaluation Vol. 1"
 
 Float_t unirnp=0.; // Initialise the random number value
 
 if (mean <= 0.) // Check the mean value
 {
  cout << " *AliRandom::Poisson* Invalid value mean = " << mean
       << " Should be positive." << endl;
 }
 else
 {
  if (mean > 80.) // Use gaussian dist. for high mean values to save time
  {
   Float_t grndm=Gauss();
   Float_t rpois=mean+grndm*sqrt(mean);
   Int_t npois=int(rpois);
   if ((rpois-float(npois)) > 0.5) npois++;
   unirnp=float(npois);
  }
  else // Construct a Poisson random number from a uniform one
  {
   Int_t npois=-1;
   Float_t expxm=exp(-mean);
   Float_t poitst=1.;
   while (poitst > expxm)
   {
    Float_t rndm=Uniform();
    npois++;
    poitst=poitst*rndm;
   }
   unirnp=float(npois);
  } // end of check for using Gauss method
 }  // end of mean validity checkn
 return unirnp;
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::Poisson(Float_t* vec,Int_t n,Float_t mean)
{
// Generate a vector of Poisson dist. random numbers with certain mean
// This saves lots of (member)function calls in case many random
// numbers are needed at once.
//
// n = The number of random numbers to be generated
//
// Method :
//
// P(n) = exp(-mean)*mean**n/n! is the Poisson distribution function
//
// with : n  = 0,1,2,3,...   and  mean > 0
//
// To generate the distribution, the "sum trick" is used as mentioned
// in "Formulae and Methods in Experimental data Evaluation Vol. 1"
 
 if (n <= 0) // Check n value
 {
  cout << " *AliRandom::Poisson* Invalid value n = " << n << endl;
 }
 else
 {
  if (mean <= 0.) // Check the mean value
  {
   cout << " *AliRandom::Poisson* Invalid value mean = " << mean
        << " Should be positive." << endl;
  }
  else
  {
   if (mean > 80.) // Use gaussian dist. for high mean values to save time
   {
    Float_t* grndm=new Float_t[n];
    Gauss(grndm,n);
    Int_t npois;
    Float_t rpois;
    for (Int_t jvec=0; jvec<n; jvec++)
    {
     rpois=mean+grndm[jvec]*sqrt(mean);
     npois=int(rpois);
     if ((rpois-float(npois)) > 0.5) npois++;
     vec[jvec]=float(npois);
    }
    delete [] grndm;
   }
   else // Construct Poisson random numbers from a uniform ones
   {
    Float_t expxm=exp(-mean);
    Int_t npois;
    Float_t poitst;
    for (Int_t jvec=0; jvec<n; jvec++)
    {
     npois=-1;
     poitst=1.;
     while (poitst > expxm)
     {
      Float_t rndm=Uniform();
      npois++;
      poitst=poitst*rndm;
     }
     vec[jvec]=float(npois);
    }
   } // end of check for using Gauss method
  }  // end of mean validity check
 }   // end of n validity check
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::SetUser(Float_t a,Float_t b,Int_t n,Float_t (*f)(Float_t))
{
// Determine the area under f(x) as a function of x
// This is called the "area function" and serves as a basis to
// provide random numbers in [a,b] according to the user defined
// distribution f(x).
// The area function is normalised such that the most extreme
// value is 1 or -1.
 
 fNa=n+1;             // The number of bins for the area function
 fXa=new Float_t[fNa];  // The binned x values of the area function
 fYa=new Float_t[fNa];  // The corresponding summed f(x) values
 fIbins=new Int_t[fNa]; // The bin numbers of the random x candidates
 
 Float_t xmin=a;
 if (a > b) xmin=b;
 Float_t step=fabs(a-b)/float(n);
 
 Float_t x;
 Float_t extreme=0;
 for (Int_t i=0; i<fNa; i++) // Fill bins of the area function
 {
  x=xmin+float(i)*step;
  fXa[i]=x;
  fYa[i]=f(x);
  if (i > 0) fYa[i]+=fYa[i-1];
  if (fabs(fYa[i]) > extreme) extreme=fabs(fYa[i]);
 }
 fYamin=fYa[0]/extreme;
 fYamax=fYa[0]/extreme;
 for (Int_t j=0; j<fNa; j++) // Normalise the area function
 {
  fYa[j]=fYa[j]/extreme;
  if (fYa[j] < fYamin) fYamin=fYa[j];
  if (fYa[j] > fYamax) fYamax=fYa[j];
 }
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::SetUser(Float_t* x,Float_t* y,Int_t n)
{
// Determine the area under y[i] as a function of x[i]
// This is called the "area function" and serves as a basis to
// provide random numbers in x according to the user provided
// distribution (x[i],y[i]).
// The area function is normalised such that the most extreme
// value is 1 or -1.
 
 fNa=n;               // The number of bins for the area function
 fXa=new Float_t[fNa];  // The binned x values of the area function
 fYa=new Float_t[fNa];  // The corresponding summed y values
 fIbins=new Int_t[fNa]; // The bin numbers of the random x candidates
 
// Order input data with increasing x
 fXa[0]=x[0];
 fYa[0]=y[0];
 for (Int_t i=1; i<fNa; i++) // Loop over x[]
 {
  for (Int_t j=0; j<i; j++) // Loop over xa[]
  {
   if (x[i] < fXa[j])
   {
    for (Int_t k=i; k>=j; k--) // Create empty position
    {
     fXa[k+1]=fXa[k];
     fYa[k+1]=fYa[k];
    }
    fXa[j]=x[i]; // Put x[] value in empty position
    fYa[j]=y[i]; // Put y[] value in empty position
    break; // Go for next x[] value
   }
   if (j == i-1) // This x[] value is the largest so far
   {
    fXa[i]=x[i]; // Put x[] value at the end of x[]
    fYa[i]=y[i]; // Put y[] value at the end of y[]
   }
  } // End loop over fXa[]
 }  // End loop over x[]
 
 Float_t extreme=0;
 for (Int_t l=0; l<fNa; l++) // Fill area function
 {
  if (l > 0) fYa[l]+=fYa[l-1];
  if (fabs(fYa[l]) > extreme) extreme=fabs(fYa[l]);
 }
 fYamin=fYa[0]/extreme;
 fYamax=fYa[0]/extreme;
 for (Int_t m=0; m<fNa; m++) // Normalise the area function
 {
  fYa[m]=fYa[m]/extreme;
  if (fYa[m] < fYamin) fYamin=fYa[m];
  if (fYa[m] > fYamax) fYamax=fYa[m];
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliRandom::User()
{
// Provide a random number according to the user defined distribution
//
// Method :
// --------
// Select by a uniform random number a certain area fraction (from fYa[])
// of the area function.
// The required random number is given by the corresponding x value (fXa[])
// of the area function.
// In case of more than one x value candidate, select randomly one of them.
 
 Float_t unirnf=0;
 
 Float_t ra=Uniform(fYamin,fYamax);     // Random value of the area function
 Float_t dmin=100.*fabs(fYamax-fYamin); // Init. the min. distance encountered
 Int_t ncand=0;
 Float_t dist;
 for (Int_t i=0; i<fNa; i++) // Search for fYa[] value(s) closest to ra
 {
  dist=fabs(ra-fYa[i]);
  if (dist <= dmin) // fYa[i] within smallest distance --> extra candidate
  {
   ncand++;
   if (dist < dmin) ncand=1; // Smallest distance so far --> THE candidate
   dmin=dist;
   fIbins[ncand-1]=i+1;
  }
 }
 
 Int_t jbin=0; // The bin number to hold the required x value
 if (ncand == 1) jbin=fIbins[0];
 
 if (ncand > 1) // Multiple x value candidates --> pick one randomly
 {
  Float_t cand=Uniform(1.,float(ncand));
  Int_t jcand=int(cand);
  if ((cand-float(jcand)) > 0.5) jcand++;
  jbin=fIbins[jcand-1];
 }
 
 if (jbin > 0) // Pick randomly a value in this x-bin
 {
  Float_t xlow=fXa[jbin-1];
  if (jbin > 1) xlow=fXa[jbin-2];
  Float_t xup=fXa[jbin-1];
  if (jbin < fNa-1) xup=fXa[jbin];
  unirnf=Uniform(xlow,xup);
 }
 
 if (ncand == 0) cout << " *AliRandom::User* No candidate found." << endl;
 
 return unirnf;
}
///////////////////////////////////////////////////////////////////////////
void AliRandom::User(Float_t* vec,Int_t n)
{
// Generate a vector of random numbers according to a user defined dist.
// This saves lots of (member)function calls in case many random
// numbers are needed at once.
//
// n = The number of random numbers to be generated
//
// Method :
// --------
// Select by a uniform random number a certain area fraction (from fYa[])
// of the area function.
// The required random number is given by the corresponding x value (fXa[])
// of the area function.
// In case of more than one x value candidate, select randomly one of them.
 
 Float_t unirnf,ra,dmin,dist;
 Int_t ncand,jbin;
 for (Int_t jvec=0; jvec<n; jvec++)
 {
  unirnf=0;
  ra=Uniform(fYamin,fYamax);     // Random value of the area function
  dmin=100.*fabs(fYamax-fYamin); // Init. the min. distance encountered
  ncand=0;
  for (Int_t i=0; i<fNa; i++) // Search for fYa[] value(s) closest to ra
  {
   dist=fabs(ra-fYa[i]);
   if (dist <= dmin) // fYa[i] within smallest distance --> extra candidate
   {
    ncand++;
    if (dist < dmin) ncand=1; // Smallest distance so far --> THE candidate
    dmin=dist;
    fIbins[ncand-1]=i+1;
   }
  }
 
  jbin=0; // The bin number to hold the required x value
  if (ncand == 1) jbin=fIbins[0];
 
  if (ncand > 1) // Multiple x value candidates --> pick one randomly
  {
   Float_t cand=Uniform(1.,float(ncand));
   Int_t jcand=int(cand);
   if ((cand-float(jcand)) > 0.5) jcand++;
   jbin=fIbins[jcand-1];
  }
 
  if (jbin > 0) // Pick randomly a value in this x-bin
  {
   Float_t xlow=fXa[jbin-1];
   if (jbin > 1) xlow=fXa[jbin-2];
   Float_t xup=fXa[jbin-1];
   if (jbin < fNa-1) xup=fXa[jbin];
   unirnf=Uniform(xlow,xup);
  }
 
  if (ncand == 0) cout << " *AliRandom::User* No candidate found." << endl;
 
  vec[jvec]=unirnf;
 
 }
}
///////////////////////////////////////////////////////////////////////////
