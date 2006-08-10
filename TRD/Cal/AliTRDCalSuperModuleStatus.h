#ifndef ALITRDCALSUPERMODULESTATUS_H
#define ALITRDCALSUPERMODULESTATUS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for status of supermodules                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTRDCalSuperModuleStatus : public TNamed {

 public:

  enum { kNdet = 540, kNstacks = 90, kNcham = 5, kNsect = 18 };
  enum { kInstalled = 1, kMasked = 2 };
  
  AliTRDCalSuperModuleStatus();
  AliTRDCalSuperModuleStatus(const Text_t* name, const Text_t* title);

  Char_t GetStatus(Int_t sm) const          { return fStatus[sm];   };
  void   SetStatus(Int_t sm, Char_t status) { fStatus[sm] = status; };

  Bool_t IsInstalled(Int_t sm) const        { return (GetStatus(sm) & kInstalled) ? kTRUE : kFALSE; }
  Bool_t IsMasked(Int_t sm) const           { return (GetStatus(sm) & kMasked)    ? kTRUE : kFALSE; }

 protected:

  Char_t fStatus[kNsect];                    //  Status byte

  ClassDef(AliTRDCalSuperModuleStatus,1)     //  Defines the status of a single super module

};

#endif
