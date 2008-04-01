#ifndef ALITRDDATAARRAYDIGITS_H
#define ALITRDDATAARRAYDIGITS_H

#include "AliTRDdataArrayS.h"

//#define DEBUG

class AliTRDdataArrayDigits;

class AliTRDdataArrayDigits : public AliTRDdataArrayS
{
 public:
  AliTRDdataArrayDigits(){};
  AliTRDdataArrayDigits(Int_t nrow, Int_t ncol, Int_t ntime);
  virtual ~AliTRDdataArrayDigits(){};

  Short_t GetDataUnchecked(Int_t row, Int_t col, Int_t time) const;
  Short_t GetData(Int_t row, Int_t col, Int_t time) const;
  Int_t   GetOverThreshold(Short_t threshold);  
  UChar_t GetPadStatus(Int_t row, Int_t col, Int_t time) const;
  void SetPadStatus(Int_t col, Int_t row, Int_t time, UChar_t status);
  Bool_t IsPadCorrupted(Int_t row, Int_t col, Int_t time);

  ClassDef(AliTRDdataArrayDigits, 1) // Container for TRD signals type of short taking pad masking into account
};
#endif
