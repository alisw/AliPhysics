#ifndef ALITRDDATAARRAYDIGITS_H
#define ALITRDDATAARRAYDIGITS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDdataArrayDigits.h,v Exp $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Container for TRD signals type of short taking pad masking into account   //
//                                                                           //
// Author:                                                                   // 
//   Markus Fasel (markus.fasel@web.de)                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDdataArrayS.h"

class AliTRDdataArrayDigits : public AliTRDdataArrayS
{

 public:

  AliTRDdataArrayDigits(){};
  AliTRDdataArrayDigits(Int_t nrow, Int_t ncol, Int_t ntime);
  virtual ~AliTRDdataArrayDigits(){};

  virtual Short_t GetDataUnchecked(Int_t row, Int_t col, Int_t time) const;
  virtual Short_t GetData(Int_t row, Int_t col, Int_t time) const;
  virtual Int_t   GetOverThreshold(Short_t threshold);  
          UChar_t GetPadStatus(Int_t row, Int_t col, Int_t time) const;
          void    SetPadStatus(Int_t col, Int_t row, Int_t time, UChar_t status);
          Bool_t  IsPadCorrupted(Int_t row, Int_t col, Int_t time);

  ClassDef(AliTRDdataArrayDigits, 1) // Container for TRD signals type of short taking pad masking into account
};
#endif
