#ifndef ALIT0LOOKUPKEY_H
#define ALIT0LOOKUPKEY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"
#include "TString.h"

class AliT0LookUpKey: public TObject
{
 public:
  AliT0LookUpKey();
  AliT0LookUpKey(Int_t key); 
  AliT0LookUpKey(TString name); 
  AliT0LookUpKey& operator= (const AliT0LookUpKey &) { return *this;};
  AliT0LookUpKey(const AliT0LookUpKey &o);
  virtual ~AliT0LookUpKey() {};
  Int_t GetKey() const {return fKey;};
  void SetKey(Int_t n)  {fKey=n;};
  TString GetChannelName() {return fName;};
  void SetChannelName(TString name) {fName = name;};
  virtual Bool_t IsEqual(const TObject *obj) const;
  void Print(Option_t* opt= "") const;
  virtual ULong_t Hash() const {return 10000*fKey;}
  //    virtual ULong_t Hash(void) const {return TString::Hash(this, sizeof(*this));};
  //   virtual void Clear(void) {fKey=0;}
 protected:
  Int_t fKey;   //logical channel number
  TString fName; //logical channel name
  
   ClassDef(AliT0LookUpKey,1)  //Hits for detector T0
};


#endif
