#ifndef ALITRDPIDINFO_H
#define ALITRDPIDINFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliTRDpidInfo : public TObject
{
public:
  struct AliTRDpidData {
    AliTRDpidData();
    virtual ~AliTRDpidData(){}
    UChar_t fPLbin;   // momentum / layer bin
    Float_t fdEdx[8]; // dEdx array
    ClassDef(AliTRDpidData, 1)  // PID layer representation
  };

  AliTRDpidInfo();
  virtual ~AliTRDpidInfo();
  AliTRDpidData const* GetData() const { return fData;}
  Int_t   GetNtracklets() const        { return fNtracklets;}
  void    PushBack(Int_t ly, Int_t p, Float_t *dedx);
  void    Reset();

private:
  Int_t         fNtracklets;  // number of tracklets
  AliTRDpidData *fData;       //[fNtracklets] PID data array

  AliTRDpidInfo(const AliTRDpidInfo::AliTRDpidInfo& ref);
  AliTRDpidInfo& operator=(const AliTRDpidInfo::AliTRDpidInfo& ref);

  ClassDef(AliTRDpidInfo, 1)  // track PID data representation
};

#endif

