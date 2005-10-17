#ifndef ALITRDCALVDRIFT_H
#define ALITRDCALVDRIFT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for Vdrift                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTRDCalROCVdrift;

class AliTRDCalVdrift : public TNamed {

 public:
 
  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };

  AliTRDCalVdrift();
  AliTRDCalVdrift(const Text_t* name, const Text_t* title);
  AliTRDCalVdrift(const AliTRDCalVdrift &c);   
  virtual ~AliTRDCalVdrift();
  AliTRDCalVdrift &operator=(const AliTRDCalVdrift &c); 

  virtual void     Copy(TObject &c) const;

  Int_t               GetDet(Int_t p, Int_t c, Int_t s) { return p+c*kNplan+s*kNplan*kNcham; };

  AliTRDCalROCVdrift *GetCalROCVdrift(Int_t d) { return fROCVdrift[d]; };
  AliTRDCalROCVdrift *GetCalROCVdrift(Int_t p, Int_t c, Int_t s) 
                                               { return fROCVdrift[GetDet(p,c,s)]; };

 protected:

  AliTRDCalROCVdrift *fROCVdrift[kNdet];           //  Array of ROC Vdrift objects

  ClassDef(AliTRDCalVdrift,1)                      //  TRD calibration class for Vdrift

};

#endif
