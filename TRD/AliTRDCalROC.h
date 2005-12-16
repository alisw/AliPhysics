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

  Int_t        GetChannel(Int_t c, Int_t r)     { return r+c*fNrows; };
  Int_t        GetNchannels()       const       { return fNchannels;   };
  Float_t      GetValue(Int_t ich) const       { return fData[ich]; };
  Float_t      GetValue(Int_t col, Int_t row)  { return fData[GetChannel(col,row)]; };

  void         SetValue(Int_t ich, Float_t vd) { fData[ich] = vd;   };
  void         SetValue(Int_t col, Int_t row, Float_t vd)
                                                { fData[GetChannel(col,row)] = vd; };

  void Scale(Float_t value);

 protected:

  Int_t     fPla;             //  Plane number
  Int_t     fCha;             //  Chamber number

  Int_t     fNrows;           //  Number of rows
  Int_t     fNcols;           //  Number of columns

  Int_t     fNchannels;             //  Number of channels
  Float_t  *fData;                //[fNchannels] Data

  ClassDef(AliTRDCalROC,1)    //  TRD ROC calibration class

};

#endif
