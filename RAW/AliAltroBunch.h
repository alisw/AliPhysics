#ifndef ALIALTROBUNCH_H
#define ALIALTROBUNCH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class AliAltroBunch: public TObject {
public:

  AliAltroBunch();
  ~ AliAltroBunch();

  const  UInt_t* GetData() const { return fData; }
  void   SetData(UInt_t *data) { fData = data; }
  Int_t  GetBunchSize()    const { return fBunchSize; }
  void   SetBunchSize(Int_t size) { fBunchSize = size; }
  UInt_t GetEndTimeBin()   const { return fEndTimeBin; }
  void   SetEndTimeBin(UInt_t bin) { fEndTimeBin = bin; }
  UInt_t GetStartTimeBin() const { return fStartTimeBin; }
  void   SetStartTimeBin(UInt_t bin) { fStartTimeBin = bin; }

private:

  AliAltroBunch& operator = (const AliAltroBunch& bunch);
  AliAltroBunch(const AliAltroBunch& bunch);

  UInt_t *fData;
  Int_t   fBunchSize;
  UInt_t  fEndTimeBin;
  UInt_t  fStartTimeBin;

  ClassDef(AliAltroBunch,0) // container class for Altro bunches
};

#endif

