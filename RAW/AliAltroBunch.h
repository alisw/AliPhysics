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
  UInt_t GetStartTimeBin() const 
    { 
      return (fEndTimeBin - (fBunchSize -1)); 
    }

  void   SetStartTimeBin(UInt_t bin) { fStartTimeBin = bin; }

private:

  AliAltroBunch& operator = (const AliAltroBunch& bunch);
  AliAltroBunch(const AliAltroBunch& bunch);

  UInt_t *fData;          // pointer to data of current bunch
  Int_t   fBunchSize;     // total size of current bunch including timestamp and the size indicator (i.e a bunch with just one sample will have size 3)
  UInt_t  fEndTimeBin;    // Time stamp of the last sample in the bunch in entities of sample indexes
  UInt_t  fStartTimeBin;  // Time index of the first bin in the bunch 

  ClassDef(AliAltroBunch,0) // container class for Altro bunches
};

#endif

