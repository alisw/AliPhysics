#ifndef ALITRDCALROC_H
#define ALITRDCALROC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalROC.h,v */

//////////////////////////////////////////////////
//                                              //
//  TRD calibration base class for one ROC      //
//                                              //
//////////////////////////////////////////////////

#include <TObject.h>

//_____________________________________________________________________________
class AliTRDCalROC : public TObject {

 public:

  AliTRDCalROC();
  AliTRDCalROC(Int_t p, Int_t c);
  AliTRDCalROC(const AliTRDCalROC &c);
  virtual           ~AliTRDCalROC();
  AliTRDCalROC      &operator=(const AliTRDCalROC &c);
  virtual void       Copy(TObject &c) const;

  Int_t    GetNrows() const                  { return fNrows; };
  Int_t    GetNcols() const                  { return fNcols; };

 protected:

  Int_t     fPla;             //  Plane number
  Int_t     fCha;             //  Chamber number

  Int_t     fNrows;           //  Number of rows
  Int_t     fNcols;           //  Number of columns

  ClassDef(AliTRDCalROC,1)    //  TRD ROC calibration class

};

#endif
