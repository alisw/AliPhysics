// @(#) $Id$

#ifndef ALIL3INTERMERGER_H
#define ALIL3INTERMERGER_H

#ifndef __CINT__ 
#include "AliHLTMerger.h"
#endif

#include "AliHLTRootTypes.h"

class AliHLTInterMerger : public AliHLTMerger {

 private:
  Int_t fPatch;
  Int_t fRowMin;
  Int_t fRowMax;
 public:
  AliHLTInterMerger();
  virtual ~AliHLTInterMerger();

  void Init(Int_t *row,Int_t p){fRowMin=row[0];fRowMax=row[1];fPatch=p;}
  void SlowMerge();
  Int_t Merge();
  void MMerge();  //Loop over tracks from different subsectors
  
  ClassDef(AliHLTInterMerger,1) //Intermerging class
};

typedef AliHLTInterMerger AliL3InterMerger; // for backward compatibility

#endif
