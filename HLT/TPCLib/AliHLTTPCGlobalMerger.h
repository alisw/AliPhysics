// @(#) $Id$
// Original: AliHLTGlobalMerger.h,v 1.6 2004/02/02 15:00:34 loizides 

#ifndef ALIHLTTPC_GLOBALMERGER_H
#define ALIHLTTPC_GLOBALMERGER_H

#ifndef  __CINT__
#include "AliHLTTPCMerger.h"
#endif

#include "AliHLTTPCRootTypes.h"

class AliHLTTPCGlobalMerger : public AliHLTTPCMerger{

 private:
  Int_t fNSlices;
  Int_t fFirst;
  Int_t fLast;

  Double_t CheckTracks(AliHLTTPCTrack *innertrack,AliHLTTPCTrack *outertrack,Int_t slice);
  
 public:
  AliHLTTPCGlobalMerger();
  virtual ~AliHLTTPCGlobalMerger();
  
  void Setup(Int_t first,Int_t last);
  void InitSlice(Int_t slice);
  void SlowMerge(Char_t *path="./");
  void Merge();  //Loop over tracks from different sectors

  ClassDef(AliHLTTPCGlobalMerger,1) //Slice merger
};

#endif
