// @(#) $Id$

#ifndef ALIL3TRANSBIT_H
#define ALIL3TRANSBIT_H

#include "AliHLTRootTypes.h"

class AliHLTTransBit {
 public:
  AliHLTTransBit();
  virtual ~AliHLTTransBit();
  inline Int_t Get0to1(Int_t val0) const;
  inline Int_t Get1to0(Int_t val1) const;
  Int_t GetBit0() const {return fBit0;}
  Int_t GetBit1() const {return fBit1;}
  Double_t GetX0() const {return fX0;}
  void SetBits(Int_t bit0, Int_t bit1) {fBit0=bit0;fBit1=bit1;}
  void SetX0(Double_t x0) {fX0=x0;}
  virtual void Update()=0;
  virtual Double_t FindOptimumX0()=0;
 protected:
  Int_t  * fTable0; //! table
  Int_t  * fTable1; //! table
  Int_t fBit0; // bit 0
  Int_t fBit1; // bit 1
  Double_t fX0; // optimal X value(?)

  ClassDef(AliHLTTransBit,1)
};

typedef AliHLTTransBit AliL3TransBit; // for backward compatibility

class AliHLTTransBitV1 : public AliHLTTransBit {
 public:
  virtual ~AliHLTTransBitV1(){}
  virtual void Update();
  virtual Double_t FindOptimumX0();
 protected:
  
  ClassDef(AliHLTTransBitV1,1)
};

typedef AliHLTTransBitV1 AliL3TransBitV1; // for backward compatibility

class AliHLTTransBitV2 : public AliHLTTransBit {
 public:
  virtual ~AliHLTTransBitV2(){}
  virtual void Update();
  virtual Double_t FindOptimumX0();
 protected:

  ClassDef(AliHLTTransBitV2,1)
};

typedef AliHLTTransBitV2 AliL3TransBitV2; // for backward compatibility

Int_t AliHLTTransBit::Get0to1(Int_t val0) const
{
  //return compressed bit values
  return fTable0[val0];
}
 
Int_t AliHLTTransBit::Get1to0(Int_t val1) const
{
  //return uncompressed bit value
  return fTable1[val1];
}

#endif 

