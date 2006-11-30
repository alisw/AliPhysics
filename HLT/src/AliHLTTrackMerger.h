// @(#) $Id$

#ifndef ALIL3TRACKMERGER_H
#define ALIL3TRACKMERGER_H

//-------------------------------------------------------------------------
//                Class AliHLTTrackMerger
//   This class is responsible for the merging of the HLT tracks
//   between TPC sectors and readout patches
//-------------------------------------------------------------------------

#ifndef __CINT__ 
#include "AliHLTMerger.h"
#endif

#include "AliHLTRootTypes.h"

class AliHLTMerger;

class AliHLTTrackMerger : public AliHLTMerger {

 public:
  AliHLTTrackMerger();
  AliHLTTrackMerger(Int_t nsubsectors);
  virtual ~AliHLTTrackMerger();

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
  void SlowMerge(AliHLTTrackArray *mergedtrack,AliHLTTrackArray *tracksin,AliHLTTrackArray *tracksout,Double_t xval);
  Int_t Merge(AliHLTTrackArray *mergedtrack,AliHLTTrackArray *tracksin,AliHLTTrackArray *tracksout);
  
  ClassDef(AliHLTTrackMerger,1) //Track merging class 
};

typedef AliHLTTrackMerger AliL3TrackMerger; // for backward compatibility

#endif
