#ifndef ALITRDCLUSTERMI_H
#define ALITRDCLUSTERMI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliTRDcluster.h"  
#include "TMath.h"  

class AliTRDrecPoint;

class AliTRDclusterMI : public AliTRDcluster {

 public:
  AliTRDclusterMI();
  AliTRDclusterMI(AliTRDcluster&cl);
  AliTRDclusterMI(const AliTRDrecPoint &p);
  void SetRelPos(Float_t pos){fRelPos = TMath::Nint(pos*128.);}
  Float_t GetRelPos(){return float(fRelPos)/128.;}
  void SetNPads(Int_t npads){fNPads = npads;}
  Char_t GetNPads(){return fNPads;}
  Float_t fRmsY;
 protected:
  Char_t fNPads;
  Char_t fRelPos;		       	 
  ClassDef(AliTRDclusterMI,1) // ClusterMI for the TRD
 
};

#endif
