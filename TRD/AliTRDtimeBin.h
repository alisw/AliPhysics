#ifndef ALITRDTIMEBIN_H
#define ALITRDTIMEBIN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtimeBin.h,v */

#include <TObject.h>

class AliTRDcluster;

const unsigned kMAX_CLUSTER_PER_TIME_BIN=3500; 

//----------------------------------------------------------------- 
class AliTRDtimeBin : public TObject {

// Provides tools to address clusters which lay within one time bin
 
public: 

  AliTRDtimeBin() {fN=0;}
  void InsertCluster(AliTRDcluster*,UInt_t);
 
  operator Int_t() const {return fN;}
  AliTRDcluster* operator[](Int_t i);
  UInt_t GetIndex(Int_t i) const {return fIndex[i];} 

  Int_t Find(Double_t y) const; 

protected:
 
   unsigned fN;
   AliTRDcluster *fClusters[kMAX_CLUSTER_PER_TIME_BIN];
   UInt_t fIndex[kMAX_CLUSTER_PER_TIME_BIN]; 

  ClassDef(AliTRDtimeBin,1) // Provides tools to address clusters which lay within one time bin

}; 

#endif 

