#ifndef ALITPCCALDET_H
#define ALITPCCALDET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TPC calibration class for parameters which are saved per detector        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTPCCalDet : public TNamed {

 public:
 
  enum { kNdet = 72 };
  AliTPCCalDet();
  AliTPCCalDet(const Text_t* name, const Text_t* title);
  AliTPCCalDet(const AliTPCCalDet &c);   
  virtual ~AliTPCCalDet();
  AliTPCCalDet &operator=(const AliTPCCalDet &c);

  virtual void     Copy(TObject &c) const;
  Float_t GetValue(Int_t d) { return fData[d]; };
  void SetValue(Int_t d, Float_t value) { fData[d] = value; };
  
  protected:

  Float_t  fData[kNdet];                          //[kNdet] Data

  ClassDef(AliTPCCalDet,1)                      //  TPC calibration class for parameters which are saved per detector

};

#endif
