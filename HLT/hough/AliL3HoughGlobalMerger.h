// @(#) $Id$

#ifndef ALIL3_HOUGHGLOBALMERGER_H
#define ALIL3_HOUGHGLOBALMERGER_H

#include "AliL3Merger.h"

class AliL3TrackArray;
class AliL3Track;

class AliL3HoughGlobalMerger {

 private:
  AliL3TrackArray **fTracks; //!
  Int_t fNSlices;

 public:
  AliL3HoughGlobalMerger();
  AliL3HoughGlobalMerger(Int_t first,Int_t last);
  virtual ~AliL3HoughGlobalMerger();
  
  void FillTracks(AliL3TrackArray *tracks,Int_t i);
  void Merge();

  ClassDef(AliL3HoughGlobalMerger,1) 
};

#endif
