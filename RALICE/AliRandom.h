#ifndef ALIRANDOM_H
#define ALIRANDOM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

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
 
 ClassDef(AliRandom,1) // Generate universal random numbers on all common machines.
};
#endif
