#ifndef ALIL3_HOUGH_MaxFinder
#define ALIL3_HOUGH_MaxFinder

#include "AliL3RootTypes.h"

class AliL3Histogram;
class AliL3TrackArray;
class AliL3HoughTrack;


class AliL3HoughMaxFinder : public TObject {
  
 private:

  Int_t fThreshold;
  AliL3Histogram *fCurrentHisto;

  Char_t fHistoType;

 public:
  AliL3HoughMaxFinder(); 
  AliL3HoughMaxFinder(Char_t *histotype,AliL3Histogram *hist=0);
  virtual ~AliL3HoughMaxFinder();

  Int_t *FindAbsMaxima();
  AliL3TrackArray *FindBigMaxima(AliL3Histogram *hist);
  AliL3TrackArray *FindMaxima(AliL3Histogram *hist,Int_t *rowrange=0,Int_t ref_row=0);
  AliL3TrackArray *LookForPeaks(AliL3Histogram *hist,Int_t nbins);
  AliL3TrackArray *LookInWindows(AliL3Histogram *hist,Int_t nbins,Int_t t1,Double_t t2,Int_t t3);
  Bool_t LocatePeak(AliL3Histogram *hist,AliL3HoughTrack *track,Int_t *xrange,Int_t *yrange,Int_t t1,Double_t t2,Int_t t3);
  AliL3HoughTrack *FindPeak(Int_t t1,Double_t t2,Int_t t3);
  AliL3HoughTrack *CalculatePeakInWindow(Int_t *maxbin,Int_t t0,Int_t t1,Double_t t2,Int_t t3);
  
  void SetThreshold(Int_t f) {fThreshold = f;}
  
  void SetHistogram(AliL3Histogram *hist) {fCurrentHisto = hist;}
  
  ClassDef(AliL3HoughMaxFinder,1)

};

#endif
