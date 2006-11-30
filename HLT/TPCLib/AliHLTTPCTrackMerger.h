// @(#) $Id$
// Original: AliHLTTrackMerger.h,v 1.6 2005/04/19 04:29:01 cvetan 
#ifndef ALIHLTTPCTRACKMERGER_H
#define ALIHLTTPCTRACKMERGER_H

//-------------------------------------------------------------------------
//                Class AliHLTTPCTrackMerger
//   This class is responsible for the merging of the HLT tracks
//   between TPC sectors and readout patches
//-------------------------------------------------------------------------

#ifndef __CINT__ 
#include "AliHLTTPCMerger.h"
#endif

#include "AliHLTTPCRootTypes.h"

class AliHLTTPCMerger;

class AliHLTTPCTrackMerger : public AliHLTTPCMerger {

 public:
  AliHLTTPCTrackMerger();
  AliHLTTPCTrackMerger(Int_t nsubsectors);
  virtual ~AliHLTTPCTrackMerger();

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
  void SlowMerge(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrackArray *tracksin,AliHLTTPCTrackArray *tracksout,Double_t xval);
  Int_t Merge(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrackArray *tracksin,AliHLTTPCTrackArray *tracksout);
  
  ClassDef(AliHLTTPCTrackMerger,1) //Track merging class 
};

#endif
