#ifndef ALIALTRODATA_H
#define ALIALTRODATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class AliAltroBunch;

class AliAltroData: public TObject {
public:

  AliAltroData();
  ~ AliAltroData();
  Bool_t NextBunch(AliAltroBunch *altrobunch);

  Int_t  GetChannel() const;
  Int_t  GetChip() const;
  Int_t  GetCard() const;
  Int_t  GetBranch() const;
  void   Reset();

  Bool_t IsComplete()      const { return fIsComplete; }
  void   SetIsComplete(Bool_t iscomplete) { fIsComplete = iscomplete; }
  Int_t  GetHadd()         const { return fHadd; }
  Int_t  GetPrevHadd()     const { return fPrevHadd; }
  Bool_t IsNewHadd()       const { return (fHadd != fPrevHadd); }
  void   SetHadd(Int_t add)      { fPrevHadd = fHadd; fHadd = add; }
  const  UInt_t* GetData() const { return fData; }
  void   SetData(UInt_t *data)   { fData = data; }
  Int_t  GetDataSize()     const { return fDataSize; }
  void   SetDataSize(Int_t size) { fDataSize = size; }

private:

  AliAltroData& operator = (const AliAltroData& altrodata);
  AliAltroData(const AliAltroData& altrodata);

  UInt_t *fData;
  UInt_t *fBunchData;
  Int_t   fDataSize;
  Int_t   fWc;
  Int_t   fHadd;
  Int_t   fPrevHadd;
  Int_t   fBunchCounter;
  Bool_t  fIsComplete;

  ClassDef(AliAltroData, 0)  // container class for Altro payload

};

#endif

