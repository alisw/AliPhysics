// @(#) $Id$

#ifndef ALIL3_GLOBALMERGER_H
#define ALIL3_GLOBALMERGER_H

#ifndef  __CINT__
#include "AliL3Merger.h"
#endif

#include "AliL3RootTypes.h"

class AliL3GlobalMerger : public AliL3Merger{

 private:
  Int_t fNSlices;
  Int_t fFirst;
  Int_t fLast;

  Double_t CheckTracks(AliL3Track *innertrack,AliL3Track *outertrack,Int_t slice);
  
 public:
  AliL3GlobalMerger();
  virtual ~AliL3GlobalMerger();
  
  void Setup(Int_t first,Int_t last);
  void InitSlice(Int_t slice);
  void SlowMerge(Char_t *path="./");
  void Merge();  //Loop over tracks from different sectors

  ClassDef(AliL3GlobalMerger,1) //Slice merger
};

#endif
