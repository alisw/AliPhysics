#ifndef ALIL3_HOUGHGLOBALMERGER_H
#define ALIL3_HOUGHGLOBALMERGER_H

#include "AliL3Merger.h"

class AliL3TrackArray;
class AliL3Track;

class AliL3HoughGlobalMerger : public AliL3Merger {

 private:
  Int_t fNSlices;
  Int_t fFirst;
  Int_t fLast;
  Double_t fMaxY;
  Double_t fMaxZ;
  Double_t fMaxKappa;
  Double_t fMaxPsi;
  Double_t fMaxTgl;
  Double_t fMaxPhi0;

 public:
  AliL3HoughGlobalMerger();
  AliL3HoughGlobalMerger(Int_t first,Int_t last);
  virtual ~AliL3HoughGlobalMerger();

  Bool_t IsTrack(AliL3Track *innertrack,AliL3Track *outertrack);
  AliL3Track *MultiMerge(AliL3TrackArray *mergedtrack,AliL3Track **tracks, Int_t ntrack);
  void FillTracks(AliL3TrackArray *tracks,Int_t slice);
  void SlowMerge();
  void Merge();  //Loop over tracks from different sectors

  ClassDef(AliL3HoughGlobalMerger,1) 
};

#endif
