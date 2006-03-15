#ifndef ALITRDCALPADSTATUS_H
#define ALITRDCALPADSTATUS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for MCM status                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "AliTRDgeometry.h"
#include "AliTRDCalSingleChamberStatus.h"

class AliTRDCalPadStatus : public TNamed {

 public:
 
  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };
  enum { kMasked = 2, kPadBridgedLeft = 4,     kPadBridgedRight = 8    };

  AliTRDCalPadStatus();
  AliTRDCalPadStatus(const Text_t* name, const Text_t* title);
  AliTRDCalPadStatus(const AliTRDCalPadStatus &c);   
  virtual ~AliTRDCalPadStatus();
  AliTRDCalPadStatus &operator=(const AliTRDCalPadStatus &c);

  virtual void     Copy(TObject &c) const;

  AliTRDCalSingleChamberStatus *GetCalROC(Int_t d) const { return fROC[d]; };
  AliTRDCalSingleChamberStatus *GetCalROC(Int_t p, Int_t c, Int_t s) const
                                               { return fROC[AliTRDgeometry::GetDetector(p,c,s)]; };

  Bool_t IsMasked(Int_t d, Int_t col, Int_t row) const { return CheckStatus(d, col, row, kMasked); };
  Bool_t IsBridgedLeft(Int_t d, Int_t col, Int_t row) const { return CheckStatus(d, col, row, kPadBridgedLeft); };
  Bool_t IsBridgedRight(Int_t d, Int_t col, Int_t row) const { return CheckStatus(d, col, row, kPadBridgedRight); };
  inline Bool_t CheckStatus(Int_t d, Int_t col, Int_t row, Int_t bitMask) const;

 protected:

  AliTRDCalSingleChamberStatus *fROC[kNdet];          //  Array of ROC objects which contain the values per pad

  ClassDef(AliTRDCalPadStatus,1)                      //  TRD calibration class for MCM status

};

Bool_t AliTRDCalPadStatus::CheckStatus(Int_t d, Int_t col, Int_t row, Int_t bitMask) const
{
  AliTRDCalSingleChamberStatus* roc = GetCalROC(d);
  if (!roc)
    return kFALSE;

  return (roc->GetStatus(col, row) & bitMask) ? kTRUE : kFALSE;
}
          
#endif
