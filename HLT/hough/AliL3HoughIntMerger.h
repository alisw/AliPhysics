// @(#) $Id$

#ifndef ALIL3HOUGHINTMERGER_H
#define ALIL3HOUGHINTMERGER_H

#include "AliL3Merger.h"

class AliL3HoughTrack;
class AliL3Track;
class AliL3TrackArray;

class AliL3HoughIntMerger : public AliL3Merger {
 
 public:
  AliL3HoughIntMerger();
  virtual ~AliL3HoughIntMerger();

  
  AliL3Track *MultiMerge(AliL3TrackArray *mergedtrack,AliL3Track **tracks, Int_t ntrack);
  Bool_t IsTrack(AliL3Track *innertrack,AliL3Track *outertrack);
  void FillTracks(AliL3TrackArray *tracks);
  void Init(Int_t *row,Int_t p){fRowMin=row[0];fRowMax=row[1];fPatch=p;}
  Int_t Merge();
  void MMerge();  //Loop over tracks from different subsectors
  void SetParameters(Double_t maxkappa=0.001, Double_t maxpsi=0.05, Double_t maxtgl=0.1);
  void SortTracks(AliL3Track **tracks, Int_t ntrack) const;
  void Print(AliL3Track **tracks);

 private:
  Int_t fPatch;//Index of the current patch
  Int_t fRowMin;//First padrow inside the patch
  Int_t fRowMax;//Last padrow inside the patch
  Double_t fMaxKappa;//Maximum track curvature
  Double_t fMaxPhi0;//Maximum phi0??
  Double_t fMaxTgl;//??

  ClassDef(AliL3HoughIntMerger,1) 
};

#endif
