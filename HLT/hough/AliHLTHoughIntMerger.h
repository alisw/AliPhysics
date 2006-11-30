// @(#) $Id$

#ifndef ALIL3HOUGHINTMERGER_H
#define ALIL3HOUGHINTMERGER_H

#include "AliHLTMerger.h"

class AliHLTHoughTrack;
class AliHLTTrack;
class AliHLTTrackArray;

class AliHLTHoughIntMerger : public AliHLTMerger {
 
 public:
  AliHLTHoughIntMerger();
  virtual ~AliHLTHoughIntMerger();

  
  AliHLTTrack *MultiMerge(AliHLTTrackArray *mergedtrack,AliHLTTrack **tracks, Int_t ntrack);
  Bool_t IsTrack(AliHLTTrack *innertrack,AliHLTTrack *outertrack);
  void FillTracks(AliHLTTrackArray *tracks);
  void Init(Int_t *row,Int_t p){fRowMin=row[0];fRowMax=row[1];fPatch=p;}
  Int_t Merge();
  void MMerge();  //Loop over tracks from different subsectors
  void SetParameters(Double_t maxkappa=0.001, Double_t maxpsi=0.05, Double_t maxtgl=0.1);
  void SortTracks(AliHLTTrack **tracks, Int_t ntrack) const;
  void Print(AliHLTTrack **tracks);

 private:
  Int_t fPatch;//Index of the current patch
  Int_t fRowMin;//First padrow inside the patch
  Int_t fRowMax;//Last padrow inside the patch
  Double_t fMaxKappa;//Maximum track curvature
  Double_t fMaxPhi0;//Maximum phi0??
  Double_t fMaxTgl;//??

  ClassDef(AliHLTHoughIntMerger,1) 
};

typedef AliHLTHoughIntMerger AliL3HoughIntMerger; // for backward comaptibility

#endif
