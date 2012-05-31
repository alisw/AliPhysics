#ifndef ALITRDPHINFO_H
#define ALITRDPHINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDPhInfo.h 27946 2008-08-13 15:26:24Z cblume $ */

//////////////////////////////////////////////////
//                                              //
//  TRD calibration base class for one ROC      //
//                                              //
//////////////////////////////////////////////////

#include <AliTRDUshortInfo.h>
#include <TMath.h>

//_____________________________________________________________________________
class AliTRDPhInfo : public AliTRDUshortInfo
{

 public:

  AliTRDPhInfo();
  AliTRDPhInfo(Int_t n);
  AliTRDPhInfo(const AliTRDPhInfo &c);
  virtual      ~AliTRDPhInfo();
  AliTRDPhInfo &operator=(const AliTRDPhInfo &c);

  Float_t At(Int_t bin) const                  { return (Float_t) (fData[bin]*3000.0/65535.0);  };
  Float_t AtS(Int_t bin) const                 { return (Float_t) (fData[bin]*3000.0/65535.0*fData[bin]*3000.0/65535.0);  };

  void   AddAt(Float_t value, Int_t bin)       { fData[bin] = (UShort_t) (value*65535.0/3000.0);  };
  void   AddAtS(Float_t value, Int_t bin)      { fData[bin] = (UShort_t) (TMath::Sqrt(TMath::Abs(value))*65535.0/3000.0); };
  
  
 protected:

  ClassDef(AliTRDPhInfo, 2)    

};

#endif
