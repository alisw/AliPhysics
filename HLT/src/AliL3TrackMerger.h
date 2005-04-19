// @(#) $Id$

#ifndef ALIL3TRACKMERGER_H
#define ALIL3TRACKMERGER_H

//-------------------------------------------------------------------------
//                Class AliL3TrackMerger
//   This class is responsible for the merging of the HLT tracks
//   between TPC sectors and readout patches
//-------------------------------------------------------------------------

#ifndef __CINT__ 
#include "AliL3Merger.h"
#endif

#include "AliL3RootTypes.h"

class AliL3Merger;

class AliL3TrackMerger : public AliL3Merger {

 public:
  AliL3TrackMerger();
  AliL3TrackMerger(Int_t nsubsectors);
  virtual ~AliL3TrackMerger();

  void SetRows(Int_t *row);
  void InitSector(Int_t sector,Int_t subsector);
  void SlowMerge();
  void Merge();  //Loop over tracks from different subsectors
  void InterMerge();

 private:
  Int_t fSubSector;//Index of the readout patch inside given TPC sector
  Int_t fNSubSector;//Number of readout patches inside given TPC sector
  Int_t *fRowMin;//!
  Int_t *fRowMax;//!
  Bool_t fSlow;//Slow or fast merging
  void SlowMerge(AliL3TrackArray *mergedtrack,AliL3TrackArray *tracksin,AliL3TrackArray *tracksout,Double_t xval);
  Int_t Merge(AliL3TrackArray *mergedtrack,AliL3TrackArray *tracksin,AliL3TrackArray *tracksout);
  
  ClassDef(AliL3TrackMerger,1) //Track merging class 
};

#endif
