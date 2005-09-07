// @(#) $Id$

#ifndef ALIHLTTPCTRACKMERGER_H
#define ALIHLTTPCTRACKMERGER_H

#ifndef __CINT__ 
#include "AliHLTTPCMerger.h"
#endif

class AliHLTTPCMerger;

class AliHLTTPCTrackMerger : public AliHLTTPCMerger {

 private:

  Int_t fSubSector;
  Int_t fNSubSector;
  Int_t *fRowMin;//!
  Int_t *fRowMax;//!
  Bool_t fSlow;
  void SlowMerge(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrackArray *tracksin,AliHLTTPCTrackArray *tracksout,Double_t xval);
  Int_t Merge(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrackArray *tracksin,AliHLTTPCTrackArray *tracksout);
 public:
  AliHLTTPCTrackMerger();
  AliHLTTPCTrackMerger(Int_t nsubsectors);
  virtual ~AliHLTTPCTrackMerger();

  void SetRows(Int_t *row);
  void InitSector(Int_t sector,Int_t subsector);
  void SlowMerge();
  void Merge();  //Loop over tracks from different subsectors
  void InterMerge();
  
  ClassDef(AliHLTTPCTrackMerger,1) //Track merging class 
};

#endif
