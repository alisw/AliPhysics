#ifndef ALIL3_HOUGH_Merger
#define ALIL3_HOUGH_Merger

#include "AliL3RootTypes.h"

class AliL3TrackArray;

class AliL3HoughMerger {
  
 private:

  Int_t fNIn;
  AliL3TrackArray **fInTracks; //!
  AliL3TrackArray *fOutTracks; //!

  void SetArray();
  void DeleteArray();
  
 public:
  AliL3HoughMerger(); 
  AliL3HoughMerger(Int_t nsubsectors);
  virtual ~AliL3HoughMerger();
  
  void FillTracks(AliL3TrackArray *tracks,Int_t patch);
  void MergePatches();
  void MergeEtaSlices(Int_t patch);
  void MergeTracks(AliL3TrackArray *intracks,Int_t i,Int_t j);

  
  //Getters
  AliL3TrackArray *GetOutTracks() {return fOutTracks;}
  AliL3TrackArray *GetInTracks(Int_t i) {if(!fInTracks) return 0; if(!fInTracks[i]) return 0; return fInTracks[i];}

  ClassDef(AliL3HoughMerger,1)

};

#endif
