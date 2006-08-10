#ifndef ALITRDCALPAD_H
#define ALITRDCALPAD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for parameters which are saved per pad                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTRDCalROC;
class AliTRDCalDet;

class AliTRDCalPad : public TNamed {

 public:
 
  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };

  AliTRDCalPad();
  AliTRDCalPad(const Text_t* name, const Text_t* title);
  AliTRDCalPad(const AliTRDCalPad &c);   
  virtual            ~AliTRDCalPad();
  AliTRDCalPad        &operator=(const AliTRDCalPad &c);

  virtual void        Copy(TObject &c) const;

  static inline Int_t GetDet(Int_t p, Int_t c, Int_t s) { return p+c*kNplan+s*kNplan*kNcham; };

  AliTRDCalROC       *GetCalROC(Int_t d) const          { return fROC[d]; };
  AliTRDCalROC       *GetCalROC(Int_t p, Int_t c, Int_t s) const
                                                        { return fROC[GetDet(p,c,s)]; };
  
  void                ScaleROCs(AliTRDCalDet* values);

 protected:

  AliTRDCalROC *fROC[kNdet];  //  Array of ROC objects which contain the values per pad

  ClassDef(AliTRDCalPad,1)    //  TRD calibration class for parameters which are saved per pad

};

#endif
