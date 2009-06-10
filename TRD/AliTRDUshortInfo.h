#ifndef ALITRDUSHORTINFO_H
#define ALITRDUSHORTINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDUshortInfo.h 27946 2008-08-13 15:26:24Z cblume $ */

//////////////////////////////////////////////////
//                                              //
//  TRD calibration base class for one ROC      //
//                                              //
//////////////////////////////////////////////////

#include <TObject.h>

//_____________________________________________________________________________
class AliTRDUshortInfo : public TObject 
{

 public:

  AliTRDUshortInfo();
  AliTRDUshortInfo(Int_t n);
  AliTRDUshortInfo(const AliTRDUshortInfo &c);
  virtual      ~AliTRDUshortInfo();
  AliTRDUshortInfo &operator=(const AliTRDUshortInfo &c);
  virtual void  Copy(TObject &c) const;

  Int_t             GetSize() const                 { return fSize; };
 
  void              SetSize(Int_t n);

 protected:

  Int_t     fSize;              // size
  UShort_t *fData;             //[fSize] Data

  ClassDef(AliTRDUshortInfo, 2)    

};

#endif
