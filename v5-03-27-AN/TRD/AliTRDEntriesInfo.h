#ifndef ALITRDENTRIESINFO_H
#define ALITRDENTRIESINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDEntriesInfo.h 27946 2008-08-13 15:26:24Z cblume $ */

//////////////////////////////////////////////////
//                                              //
//  TRD calibration base class for one ROC      //
//                                              //
//////////////////////////////////////////////////

#include <AliTRDUshortInfo.h>

//_____________________________________________________________________________
class AliTRDEntriesInfo : public AliTRDUshortInfo
{

 public:

  AliTRDEntriesInfo();
  AliTRDEntriesInfo(Int_t n);
  AliTRDEntriesInfo(const AliTRDEntriesInfo &c);
  virtual      ~AliTRDEntriesInfo();
  AliTRDEntriesInfo &operator=(const AliTRDEntriesInfo &c);
  

  Int_t  At(Int_t bin) const                { return (Int_t) fData[bin]; };

  void   AddAt(Int_t value, Int_t bin)      { fData[bin] = (UShort_t) value;  };
  
  //
  // statistic
  //
  Int_t GetSum() const; 
  
  // algebra
  Bool_t TestAdd(const AliTRDEntriesInfo * info);
  void   Add(const AliTRDEntriesInfo *info);
  void   AddIf(const AliTRDEntriesInfo *info);
  
 protected:

  ClassDef(AliTRDEntriesInfo, 2)    

};

#endif
