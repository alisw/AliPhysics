#ifndef ALITRDCALDET_H
#define ALITRDCALDET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for parameters which are saved per detector        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "AliTRDgeometry.h"

class AliTRDCalDet : public TNamed {

 public:
 
  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };

  AliTRDCalDet();
  AliTRDCalDet(const Text_t* name, const Text_t* title);
  AliTRDCalDet(const AliTRDCalDet &c);   
  virtual ~AliTRDCalDet();
  AliTRDCalDet &operator=(const AliTRDCalDet &c);

  virtual void     Copy(TObject &c) const;

  Float_t GetValue(Int_t d) { return fData[d]; };
  Float_t GetValue(Int_t p, Int_t c, Int_t s) { return fData[AliTRDgeometry::GetDetector(p,c,s)]; };

  void SetValue(Int_t d, Float_t value) { fData[d] = value; };
  void SetValue(Int_t p, Int_t c, Int_t s, Float_t value) { fData[AliTRDgeometry::GetDetector(p,c,s)] = value; };
  
  protected:

  Float_t  fData[kNdet];                          //[kNdet] Data

  ClassDef(AliTRDCalDet,1)                      //  TRD calibration class for parameters which are saved per detector

};

#endif
