#ifndef ALITRDPRFINFO_H
#define ALITRDPRFINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDPrfInfo.h 27946 2008-08-13 15:26:24Z cblume $ */

//////////////////////////////////////////////////
//                                              //
//  TRD calibration base class for one ROC      //
//                                              //
//////////////////////////////////////////////////

#include <TObject.h>


//_____________________________________________________________________________
class AliTRDPrfInfo : public TObject 
{

 public:

  AliTRDPrfInfo();
  AliTRDPrfInfo(Int_t n);
  AliTRDPrfInfo(const AliTRDPrfInfo &c);
  virtual      ~AliTRDPrfInfo();
  AliTRDPrfInfo &operator=(const AliTRDPrfInfo &c);
  virtual void  Copy(TObject &c) const;

  Int_t         GetSize()  const                { return fSize;  };
  Float_t       At(Int_t bin) const             { return (Float_t) (fData[bin]/255.0); };
  
  void          AddAt(Float_t value, Int_t bin)   { fData[bin] = (UChar_t) (value*255.0); };
  void          SetSize(Int_t size);
  
 protected:
  
  Int_t     fSize;              //  Size
  UChar_t  *fData;              //[fSize] Data

  ClassDef(AliTRDPrfInfo, 2)   

};

#endif
