// @(#) $Id$

#ifndef ALIL3TRANSBIT_H
#define ALIL3TRANSBIT_H

#include "AliL3RootTypes.h"

class AliL3TransBit {
 public:
  AliL3TransBit();
  virtual ~AliL3TransBit();
  inline Int_t Get0to1(Int_t val0);
  inline Int_t Get1to0(Int_t val1);
  Int_t GetBit0() {return fBit0;}
  Int_t GetBit1() {return fBit1;}
  Double_t GetX0() {return fX0;}
  void SetBits(Int_t bit0, Int_t bit1) {fBit0=bit0;fBit1=bit1;}
  void SetX0(Double_t x0) {fX0=x0;}
  virtual void Update()=0;
  virtual Double_t FindOptimumX0()=0;
 protected:
  Int_t  * fTable0; //! table
  Int_t  * fTable1; //! table
  Int_t fBit0;
  Int_t fBit1;
  Double_t fX0;

  ClassDef(AliL3TransBit,1)
};

class AliL3TransBit_v1 : public AliL3TransBit {
 public:
  virtual void Update();
  virtual Double_t FindOptimumX0();
 protected:
  
  ClassDef(AliL3TransBit_v1,1)
};

class AliL3TransBit_v2 : public AliL3TransBit {
 public:
  virtual void Update();
  virtual Double_t FindOptimumX0();
 protected:

  ClassDef(AliL3TransBit_v2,1)
};

Int_t AliL3TransBit::Get0to1(Int_t val0)
{
  //return compressed bit values
  return fTable0[val0];
}
 
Int_t AliL3TransBit::Get1to0(Int_t val1)
{
  //return uncompressed bit value
  return fTable1[val1];
}

#endif 

