#ifndef ALIRANDOM_H
#define ALIRANDOM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

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
//--- NvE 11-oct-1997 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////
 
#include <iostream.h>
#include <math.h>
 
#include "TObject.h"
 
class AliRandom : public TObject
{
 public:
  AliRandom();                                 // Constructor with default sequence
  AliRandom(Int_t seed);                       // Constructor with user defined seed
  AliRandom(Int_t seed,Int_t cnt1,Int_t cnt2); // User defined starting point
  ~AliRandom();                                // Destructor
  Int_t GetSeed();                             // Provide current seed value
  Int_t GetCnt1();                             // Provide current counter value cnt1
  Int_t GetCnt2();                             // Provide current counter value cnt2
  void Info();                                 // Print current seed, cnt1 and cnt2
  Float_t Uniform();                           // Uniform dist. within <0,1>
  Float_t Uniform(Float_t a,Float_t b);        // Uniform dist. within <a,b>
  void Uniform(Float_t* vec,Int_t n);          // n uniform randoms in <0,1>
  void Uniform(Float_t* vec,Int_t n,Float_t a,Float_t b); // see above
  Float_t Gauss();                             // Gaussian dist. with mean=0 sigma=1
  Float_t Gauss(Float_t mean,Float_t sigma);   // Gaussian dist. with mean and sigma
  void Gauss(Float_t* vec,Int_t n);            // n Gaussian randoms mean=0 sigma=1
  void Gauss(Float_t* vec,Int_t n,Float_t mean,Float_t sigma); // see above
  Float_t Poisson(Float_t mean);               // Poisson dist. with certain mean
  void Poisson(Float_t* vec,Int_t n,Float_t mean); // n Poisson randoms with mean
  void SetUser(Float_t a,Float_t b,Int_t n,Float_t (*f)(Float_t)); // User dist. f(x)
  void SetUser(Float_t* x,Float_t* y,Int_t n); // User dist. arrays
  Float_t User(); // Provide random in [a,b] according to user distribution
  void User(Float_t* vec,Int_t n); // n randoms in [a,b] from user dist.
 
 private:
  Int_t   fI,fJ,fSeed,fCnt1,fCnt2,fClip;                       // Indices, seed and counters
  Float_t fU[97],fC,fCd,fCm;                                   // The Fibonacci parameters
  void Start(Int_t seed,Int_t cnt1,Int_t cnt2);                // Start at certain point
  void Unpack(Int_t seed,Int_t& i,Int_t& j,Int_t& k,Int_t& l); // Unpack the seed
  void Uniform(Int_t n); // n uniform randoms for quick skipping
  Int_t fNa;             //! The number of bins of the area function
  Float_t* fXa;          //! The binned x values of the area function
  Float_t* fYa;          //! The corresponding y values of the area function
  Float_t fYamin,fYamax; //! The min. and max. y values of the area function
  Int_t* fIbins;         //! The bin numbers of the random x candidates
 
 ClassDef(AliRandom,1) // Class definition to enable ROOT I/O
};
#endif
