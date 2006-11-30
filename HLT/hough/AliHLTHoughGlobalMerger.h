// @(#) $Id$

#ifndef ALIL3_HOUGHGLOBALMERGER_H
#define ALIL3_HOUGHGLOBALMERGER_H

#include "AliHLTMerger.h"

class AliHLTTrackArray;
class AliHLTTrack;

class AliHLTHoughGlobalMerger {

 private:
  AliHLTTrackArray **fTracks; //!
  Int_t fNSlices;

 public:
  AliHLTHoughGlobalMerger();
  AliHLTHoughGlobalMerger(Int_t first,Int_t last);
  virtual ~AliHLTHoughGlobalMerger();
  
  void FillTracks(AliHLTTrackArray *tracks,Int_t i);
  void Merge();

  ClassDef(AliHLTHoughGlobalMerger,1) 
};

typedef AliHLTHoughGlobalMerger AliL3HoughGlobalMerger; // for backward compatibility

#endif
