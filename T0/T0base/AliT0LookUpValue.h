#ifndef ALIT0LOOKUPVALUE_H
#define ALIT0LOOKUPVALUE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
#include "TObject.h"

class AliT0LookUpValue: public TObject
{
 public:
 
  AliT0LookUpValue(); 
  AliT0LookUpValue(Int_t trm, Int_t tdc, Int_t chain, Int_t channel );
 AliT0LookUpValue(const AliT0LookUpValue &o): TObject(),
    fTRM(-1),
    fTDC(-1),
    fChain(-1),
    fChannel(-1)
    { ((AliT0LookUpValue&) o).Copy(*this);} 
  AliT0LookUpValue& operator= (const AliT0LookUpValue & ) { return *this;};
  virtual ~AliT0LookUpValue() {};
 
  virtual Bool_t IsEqual(const TObject* obj) const ;
  virtual ULong_t Hash() const {return (ULong_t) 1000000*fTRM+100000*fTDC+1000*fChannel;}
  //   virtual ULong_t Hash(void) const {return TString::Hash(this, sizeof(*this));};
   Int_t GetTRM() const {return fTRM;};
   Int_t GetTDC() const {return fTDC;};
   Int_t GetChain() const {return fChain;};
   Int_t GetChannel() const {return fChannel;};
   void SetTRM(Int_t n) {fTRM=n;};
   void SetTDC(Int_t n) {fTDC=n;};
   void SetChain(Int_t n) {fChain=n;};
   void SetChannel(Int_t n) {fChannel=n;};
   virtual void Clear () {fTRM = -1; fTDC=-1; fChain=-1;  fChannel=-1;}
   void Clear(const Option_t*) {};
   void Print(Option_t* opt= "") const;

 protected:

   Int_t fTRM;       //#TRM
   Int_t fTDC;       //#TDC
   Int_t fChain;     //#chain 
   Int_t fChannel;   //#channel

   ClassDef(AliT0LookUpValue,1)  //Hits for detector T0
};

#endif
