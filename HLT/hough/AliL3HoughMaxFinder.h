#ifndef ALIL3_HOUGH_MaxFinder
#define ALIL3_HOUGH_MaxFinder

#include "AliL3RootTypes.h"

class AliL3TrackArray;

class AliL3HoughMaxFinder : public TObject {
  
 private:

  Int_t fThreshold;
  //AliL3TrackArray *fTracks; //!

  Char_t fHistoType;

 public:
  AliL3HoughMaxFinder(); 
  AliL3HoughMaxFinder(Char_t *histotype);
  virtual ~AliL3HoughMaxFinder();

  AliL3TrackArray *FindMaxima(TH2F *hist,Int_t *rowrange=0,Int_t ref_row=0);
  AliL3TrackArray *FindPeak(TH2F *hist,Int_t t1,Double_t t2,Int_t t3);
  void SetThreshold(Int_t f) {fThreshold = f;}
  

  ClassDef(AliL3HoughMaxFinder,1)

};

#endif
