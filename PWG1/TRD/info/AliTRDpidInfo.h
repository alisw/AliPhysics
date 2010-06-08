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
  class AliTRDpidData {
  public:
    AliTRDpidData();
    virtual ~AliTRDpidData(){}
    Int_t   Layer() const    { return (fPLbin&0xf0)>>4;}
    Int_t   Momentum() const { return fPLbin&0xf;}

    UChar_t fPLbin;   // momentum / layer bin
    Float_t fdEdx[8]; // dEdx array
    ClassDef(AliTRDpidData, 1)  // PID layer representation
  };

  AliTRDpidInfo();
  AliTRDpidInfo(Int_t idx);
  virtual ~AliTRDpidInfo();
  inline AliTRDpidData const* GetData(Int_t itrklt) const;
  AliTRDpidData const* GetDataInLayer(Int_t ily) const;
  Int_t   GetNtracklets() const        { return fNtracklets;}
  Char_t  GetPID() const               { return fPID;}
  void    PushBack(Int_t ly, Int_t p, const Float_t *dedx);
  void    Reset();
  void    SetPID(Int_t idx)             { fPID = idx;}

private:
  Char_t        fPID;         // reference PID
  Int_t         fNtracklets;  // number of tracklets
  AliTRDpidData *fData;       //[fNtracklets] PID data array

  AliTRDpidInfo(const AliTRDpidInfo& ref);
  AliTRDpidInfo& operator=(const AliTRDpidInfo& ref);

  ClassDef(AliTRDpidInfo, 1)  // track PID data representation
};


//________________________________________________________________________
AliTRDpidInfo::AliTRDpidData const* AliTRDpidInfo::GetData(Int_t itrklt) const
{
// Retrive itrklt-th tracklet independent of layer
// For the layer specific getter see GetDataInLayer().

  if(!fData) return NULL;
  if(itrklt<0 || itrklt>=fNtracklets) return NULL;
  return &fData[itrklt];
}

#endif

