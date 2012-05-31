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

class AliRndm 
{
public:
  AliRndm();
  AliRndm(const AliRndm &rn);
  virtual ~AliRndm() {fRandom=0;}
  AliRndm & operator=(const AliRndm& rn) 
    {rn.Copy(*this); return (*this);}
  
  // Random number generator bit
  virtual void SetRandom(TRandom *ran=0)
  {if(ran) fRandom=ran;
  else fRandom=gRandom;}

  virtual TRandom* GetRandom() const {return fRandom;}
  virtual void Rndm(Float_t* array, Int_t size) const; 
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
  void Copy(AliRndm &rn) const;

  ClassDef(AliRndm,1)  //Random Number generator wrapper
};

#endif 

