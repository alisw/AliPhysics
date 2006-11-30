// @(#) $Id$

#ifndef ALIL3_GLOBALMERGER_H
#define ALIL3_GLOBALMERGER_H

#ifndef  __CINT__
#include "AliHLTMerger.h"
#endif

#include "AliHLTRootTypes.h"

class AliHLTGlobalMerger : public AliHLTMerger{

 private:
  Int_t fNSlices;
  Int_t fFirst;
  Int_t fLast;

  Double_t CheckTracks(AliHLTTrack *innertrack,AliHLTTrack *outertrack,Int_t slice);
  
 public:
  AliHLTGlobalMerger();
  virtual ~AliHLTGlobalMerger();
  
  void Setup(Int_t first,Int_t last);
  void InitSlice(Int_t slice);
  void SlowMerge(Char_t *path="./");
  void Merge();  //Loop over tracks from different sectors

  ClassDef(AliHLTGlobalMerger,1) //Slice merger
};

typedef AliHLTGlobalMerger AliL3GlobalMerger; // for backward compatibility

#endif
