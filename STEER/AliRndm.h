#ifndef ALIRNDM_H
#define ALIRNDM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Random Number Interface                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TRandom.h>

static TRandom *sRandom;

class AliRndm 
{
public:
  AliRndm() {SetRandom();}
  virtual ~AliRndm() {fRandom=sRandom=0;}
  
  // Random number generator bit
  virtual void SetRandom(TRandom *ran=0)
  {if(ran) fRandom=ran;
  else fRandom=sRandom=gRandom;}

  virtual TRandom* GetRandom() const {return fRandom;}
  virtual void Rndm(Float_t* array, const Int_t size) const; 
#ifdef CKNONE
  virtual Float_t Rndm() const {return fRandom->Rndm();}
#else
  virtual Float_t Rndm() const {
    Float_t r;
    do r=fRandom->Rndm(); while(0>=r || r>=1); return r;}
#endif
  virtual void WriteRandom(const char *filename) const;
  virtual void ReadRandom(const char *filename);

  protected:
  TRandom *fRandom;       // Pointer to the random number generator

  private:
  AliRndm(const AliRndm &) {}
  AliRndm & operator=(const AliRndm &) {return (*this);}

  ClassDef(AliRndm,1)  //Random Number generator wrapper
};

#endif 

