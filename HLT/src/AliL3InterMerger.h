#ifndef ALIL3INTERMERGER_H
#define ALIL3INTERMERGER_H

#ifndef __CINT__ 
#include "AliL3Merger.h"
#endif

class AliL3InterMerger : public AliL3Merger {

 private:
  Int_t fPatch;
  Int_t fRowMin;
  Int_t fRowMax;
 public:
  AliL3InterMerger();
  virtual ~AliL3InterMerger();

  void Init(Int_t *row,Int_t p){fRowMin=row[0];fRowMax=row[1];fPatch=p;}
  void SlowMerge();
  Int_t Merge();
  void MMerge();  //Loop over tracks from different subsectors
  
  ClassDef(AliL3InterMerger,1) 
};

#endif
