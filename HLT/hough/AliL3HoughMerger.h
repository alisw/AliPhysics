// @(#) $Id$

#ifndef ALIL3_HOUGH_Merger
#define ALIL3_HOUGH_Merger

#include "AliL3RootTypes.h"
#include "AliL3Merger.h"

class AliL3TrackArray;
class AliL3Track;

class AliL3HoughMerger : public AliL3Merger {
  
 private:
  Double_t fMaxY;
  Double_t fMaxZ;
  Double_t fMaxKappa;
  Double_t fMaxPsi;
  Double_t fMaxTgl;
  Double_t fMaxPhi0;
  Bool_t fSlow;

  void Merge();
  Int_t Merge(AliL3TrackArray* mergedtrack,AliL3TrackArray *tracksin,AliL3TrackArray *tracksout);
  void SlowMerge(AliL3TrackArray *mergedtrack,AliL3TrackArray *tracksin,AliL3TrackArray *tracksout,Double_t xval);

 public:
  AliL3HoughMerger(); 
  AliL3HoughMerger(Int_t nsubsectors);
  virtual ~AliL3HoughMerger();
  
  virtual Bool_t IsTrack(AliL3Track *innertrack,AliL3Track *outertrack);
  virtual AliL3Track *MultiMerge(AliL3TrackArray *mergedtrack,AliL3Track **tracks, Int_t ntrack);
  virtual void AddTrack(AliL3TrackArray *mergedtrack,AliL3Track *track);
  void FillTracks(AliL3TrackArray *tracks,Int_t patch);
  
  void MergePatches(Bool_t slow=kTRUE);
  void MergeEtaSlices(Int_t /*patch*/) {};
  void MergeTracks(AliL3TrackArray */*intracks*/,Int_t /*i*/,Int_t /*j*/) {};
  void Print(AliL3Track **tracks);
  void SetParameters(Double_t maxkappa=0.001, Double_t maxpsi=0.05,Double_t maxphi0=0.1);
  
  ClassDef(AliL3HoughMerger,1) //Patch merger for Hough tracklets

};

#endif
