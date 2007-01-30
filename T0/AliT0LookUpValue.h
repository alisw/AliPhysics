#ifndef ALIT0LOOKUPVALUE_H
#define ALIT0LOOKUPVALUE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
#include "TObject.h"
#include "TMap.h"

class AliT0LookUpValue: public TObject
{
 public:
 
  AliT0LookUpValue(); 
  AliT0LookUpValue(Int_t trm, Int_t tdc, Int_t chain, Int_t channel );
  virtual Bool_t IsEqual(const TObject* obj) const ;
  virtual ULong_t Hash(void) const;
   Int_t GetTRM(void) const {return fTRM;};
   Int_t GetTDC(void) const {return fTDC;};
   Int_t GetChain(void) const {return fChain;};
   Int_t GetChannel(void) const {return fChannel;};
   void SetTRM(Int_t n) {fTRM=n;};
   void SetTDC(Int_t n) {fTDC=n;};
   void SetChain(Int_t n) {fChain=n;};
   void SetChannel(Int_t n) {fChannel=n;};
   virtual void Clear (void) {fTRM = -1; fTDC=-1; fChain=-1;  fChannel=-1;}
//   void Clear();

 protected:

   Int_t fTRM;       //#TRM
   Int_t fTDC;       //#TDC
   Int_t fChain;     //#chain 
   Int_t fChannel;   //#channel

   ClassDef(AliT0LookUpValue,1)  //Hits for detector T0
};

class AliT0LookUpKey: public TObject
{
 public:
  AliT0LookUpKey();
  AliT0LookUpKey(Int_t key); 
  Int_t GetKey(void) const {return fKey;};
  void SetKey(Int_t n)  {fKey=n;};
 virtual Bool_t IsEqual(const TObject *obj) const;
  virtual ULong_t Hash(void) const {return fKey;};
//   virtual void Clear(void) {fKey=0;}
 protected:
  Int_t fKey;   //logical channel name
   ClassDef(AliT0LookUpKey,1)  //Hits for detector T0
};

#endif
