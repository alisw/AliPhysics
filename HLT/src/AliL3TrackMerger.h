#ifndef ALIL3TRACKMERGER_H
#define ALIL3TRACKMERGER_H


#ifndef __CINT__ 
#include "AliL3Merger.h"
#endif

class AliL3Merger;

class AliL3TrackMerger : public AliL3Merger {

 private:

  Int_t fSubSector;
  Int_t fNSubSector;
  Int_t fRowMin[5];
  Int_t fRowMax[5];
  Bool_t fSlow;
  void SlowMerge(AliL3TrackArray *mergedtrack,AliL3TrackArray *tracksin,AliL3TrackArray *tracksout,Double_t xval);
  Int_t Merge(AliL3TrackArray *mergedtrack,AliL3TrackArray *tracksin,AliL3TrackArray *tracksout);
 public:
  AliL3TrackMerger();
  AliL3TrackMerger(Int_t nsubsectors);
  virtual ~AliL3TrackMerger();

  void SetRows(Int_t *row);
  void InitSector(Int_t sector,Int_t subsector);
  void SlowMerge();
  void Merge();  //Loop over tracks from different subsectors
  void InterMerge();
  
  ClassDef(AliL3TrackMerger,1) 
};

#endif
