#ifndef ALIL3_HOUGH_Merge
#define ALIL3_HOUGH_Merge

#include "AliL3RootTypes.h"

class TClonesArray;
class TObjArray;
class TH2F;
class AliL3TrackArray;

class AliL3HoughMerge : public TObject {
  
 private:

  AliL3TrackArray **fInTracks; //!
  AliL3TrackArray *fOutTracks; //!
  
  Int_t fNPatches;
 public:
  AliL3HoughMerge(); 
  AliL3HoughMerge(Int_t slice,Int_t row_patches=5);
  virtual ~AliL3HoughMerge();

  void FillTracks(AliL3TrackArray *tracks,Int_t patch);
  void FillHisto(TH2F *merge_hist);
  
  ClassDef(AliL3HoughMerge,1)

};

#endif
