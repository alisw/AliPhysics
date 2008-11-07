#ifndef ALIALTROBUNCH_H
#define ALIALTROBUNCH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

#define DECODERERROR -3

class AliAltroBunch: public TObject {
public:

  AliAltroBunch();
  ~ AliAltroBunch();
  inline const  UInt_t* GetData() const { return fData; }

  inline void   SetData(UInt_t *data) 
  {
    fData = data; 
  }
  
  inline Int_t  GetBunchSize()    const { return fBunchSize; }
  inline void SetBunchSize(Int_t size) 
  {
    fPrevBunchSize =  fBunchSize;  fBunchSize = size; 
  }
 
  inline int CheckConsistency()    
  {
    if( fPrevEndTimeBin <= fEndTimeBin + fPrevBunchSize )
      {
	//	printf("%s:%d, ERROR conistency check failed\n", __FILE__ , __LINE__ );
	return  DECODERERROR;
     }
    else
      {
	return kTRUE;
      }
  }

  UInt_t GetEndTimeBin()   const { return fEndTimeBin; }
  inline void   SetEndTimeBin(UInt_t bin) {fPrevEndTimeBin =  fEndTimeBin;  fEndTimeBin = bin; }
  UInt_t GetStartTimeBin() const { return (fEndTimeBin - (fBunchSize -1)); }
  void   SetStartTimeBin(UInt_t bin) { fStartTimeBin = bin; }

private:
  AliAltroBunch& operator = (const AliAltroBunch& bunch);
  AliAltroBunch(const AliAltroBunch& bunch);
  UInt_t *fData;          // pointer to data of current bunch
  Int_t   fBunchSize;     // total size of current bunch including timestamp and the size indicator (i.e a bunch with just one sample will have size 3)
  UInt_t  fEndTimeBin;    // Time stamp of the last sample in the bunch in entities of sample indexes
  UInt_t  fStartTimeBin;  // Time index of the first bin in the bunch 
  UInt_t  fPrevBunchSize; // Previous bunch size  
  UInt_t  fPrevEndTimeBin;// Previous end time bin
  // bool fIsFirstBunch;

  ClassDef(AliAltroBunch,0) // container class for Altro bunches
};

#endif

