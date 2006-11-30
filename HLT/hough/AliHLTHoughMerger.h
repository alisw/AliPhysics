// @(#) $Id$

#ifndef ALIL3HOUGHMERGER_H
#define ALIL3HOUGHMERGER_H

#include "AliHLTRootTypes.h"
#include "AliHLTMerger.h"

class AliHLTTrackArray;
class AliHLTTrack;

class AliHLTHoughMerger : public AliHLTMerger {

 public:
  AliHLTHoughMerger(); 
  AliHLTHoughMerger(Int_t nsubsectors);
  virtual ~AliHLTHoughMerger();
  
  virtual Bool_t IsTrack(AliHLTTrack *innertrack,AliHLTTrack *outertrack);
  virtual AliHLTTrack *MultiMerge(AliHLTTrackArray *mergedtrack,AliHLTTrack **tracks, Int_t ntrack);
  virtual void AddTrack(AliHLTTrackArray *mergedtrack,AliHLTTrack *track);
  void FillTracks(AliHLTTrackArray *tracks,Int_t patch);
  
  void MergePatches(Bool_t slow=kTRUE);
  void MergeEtaSlices(Int_t /*patch*/) {};
  void MergeTracks(AliHLTTrackArray */*intracks*/,Int_t /*i*/,Int_t /*j*/) {};
  void Print(AliHLTTrack **tracks);
  void SetParameters(Double_t maxkappa=0.001, Double_t maxpsi=0.05,Double_t maxphi0=0.1);
  
 private:
  Double_t fMaxY;//maximum Y
  Double_t fMaxZ;//maximum Z
  Double_t fMaxKappa;//maximum track curvature
  Double_t fMaxPsi;//maximum track emission angle
  Double_t fMaxTgl;//??
  Double_t fMaxPhi0;//??
  Bool_t fSlow;//??

  void Merge();
  Int_t Merge(AliHLTTrackArray* mergedtrack,AliHLTTrackArray *tracksin,AliHLTTrackArray *tracksout);
  void SlowMerge(AliHLTTrackArray *mergedtrack,AliHLTTrackArray *tracksin,AliHLTTrackArray *tracksout,Double_t xval);
  
  ClassDef(AliHLTHoughMerger,1) //Patch merger for Hough tracklets

};

typedef AliHLTHoughMerger AliL3HoughMerger; // for backward compatibility

#endif
