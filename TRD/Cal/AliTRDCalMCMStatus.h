#ifndef ALITRDCALMCMSTATUS_H
#define ALITRDCALMCMSTATUS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for MCM status                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTRDCalSingleChamberStatus;

class AliTRDCalMCMStatus : public TNamed {

 public:
 
  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };
  enum { kMasked = 2, kMCMTemperatureBit1 = 4, kMBMTemperatureBit2 = 8 };

  AliTRDCalMCMStatus();
  AliTRDCalMCMStatus(const Text_t* name, const Text_t* title);
  AliTRDCalMCMStatus(const AliTRDCalMCMStatus &c);   
  virtual ~AliTRDCalMCMStatus();
  AliTRDCalMCMStatus &operator=(const AliTRDCalMCMStatus &c);

  virtual void  Copy(TObject &c) const;


         Bool_t IsMasked(Int_t d, Int_t col, Int_t row) const 
                                                         { return CheckStatus(d, col, row, kMasked); };
         Bool_t CheckStatus(Int_t d, Int_t col, Int_t row, Int_t bitMask) const;
         Char_t GetStatus(Int_t d, Int_t col, Int_t row) const;

         AliTRDCalSingleChamberStatus *GetCalROC(Int_t p, Int_t c, Int_t s) const;
         AliTRDCalSingleChamberStatus *GetCalROC(Int_t d) const 
                                                         { return fROC[d]; };

 protected:

  AliTRDCalSingleChamberStatus *fROC[kNdet];          //  Array of ROC objects which contain the values per pad

  ClassDef(AliTRDCalMCMStatus,1)                      //  TRD calibration class for MCM status

};

#endif
