#ifndef ALITRDTIMEBIN_H
#define ALITRDTIMEBIN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtimeBin.h,v */

//////////////////////////////////////////////////////////////////////
//                                                                  //
//  Hit compression class                                           //
//  Adapted from AliTPCTimeBin by Marian                            //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliTRDcluster;

//----------------------------------------------------------------- 
class AliTRDtimeBin : public TObject {

// Provides tools to address clusters which lay within one time bin
 
public: 

  AliTRDtimeBin();
  virtual ~AliTRDtimeBin() { };
  void InsertCluster(AliTRDcluster *c, UInt_t index);
 
  operator Int_t() const {return fN;}
  AliTRDcluster* operator[](Int_t i);
  UInt_t GetIndex(Int_t i) const {return fIndex[i];} 

  Int_t Find(Double_t y) const; 

protected:

  enum { kMaxClusterPerTimeBin=3500 };
 
   UInt_t         fN;                                 // ????
   AliTRDcluster *fClusters[kMaxClusterPerTimeBin];   // ????
   UInt_t         fIndex[kMaxClusterPerTimeBin];      // ????

  ClassDef(AliTRDtimeBin,1) // Provides tools to address clusters which lay within one time bin

}; 

#endif 

