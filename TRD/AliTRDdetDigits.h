#ifndef TRDdetDigits_H
#define TRDdetDigits_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////
//  Digits container for one TRD detector  //
/////////////////////////////////////////////

#include "AliDigits.h"

//_____________________________________________________________________________
class AliTRDdetDigits : public AliDigits {

 public:

  AliTRDdetDigits();
  ~AliTRDdetDigits() { };

  virtual void Allocate(Int_t nrow, Int_t ncol, Int_t ntime);
  virtual void SetDigit(Int_t row, Int_t col, Int_t time, Short_t value);

 protected:

  Int_t   fNrowTRD;                  // Number of rows in the TRD
  Int_t   fNcolTRD;                  // Number of colums in the TRD

  ClassDef(AliTRDdetDigits,1)        // Digits container for one TRD detector

};

#endif
