// @(#) $Id$

#ifndef ALIHLTTPCINTERMERGER_H
#define ALIHLTTPCINTERMERGER_H

#ifndef __CINT__ 
#include "AliHLTTPCMerger.h"
#endif

#include "AliHLTTPCRootTypes.h"

class AliHLTTPCInterMerger : public AliHLTTPCMerger {

 private:
  Int_t fPatch;
  Int_t fRowMin;
  Int_t fRowMax;
 public:
  AliHLTTPCInterMerger();
  virtual ~AliHLTTPCInterMerger();

  void Init(Int_t *row,Int_t p){fRowMin=row[0];fRowMax=row[1];fPatch=p;}
  void SlowMerge();
  Int_t Merge();
  void MMerge();  //Loop over tracks from different subsectors
  
  ClassDef(AliHLTTPCInterMerger,1) //Intermerging class
};

#endif
