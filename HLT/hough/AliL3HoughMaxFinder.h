#ifndef ALIL3_HOUGH_MaxFinder
#define ALIL3_HOUGH_MaxFinder

#include "AliL3RootTypes.h"

class AliL3Histogram;
class AliL3TrackArray;
class AliL3HoughTrack;

struct AxisWindow
{
  Int_t ymin;
  Int_t ymax;
  Int_t xbin;
  Int_t weight;
};

class AliL3HoughMaxFinder {
  
 private:

  Int_t fThreshold;
  AliL3Histogram *fCurrentHisto;  //!

  Char_t fHistoType;

 public:
  AliL3HoughMaxFinder(); 
  AliL3HoughMaxFinder(Char_t *histotype,AliL3Histogram *hist=0);
  virtual ~AliL3HoughMaxFinder();

  void FindAbsMaxima(Int_t &max_xbin,Int_t &max_ybin);
  AliL3TrackArray *FindBigMaxima(AliL3Histogram *hist);
  void FindMaxima(Float_t *xpeaks,Float_t *ypeaks,Int_t *weight,Int_t &entries);
  AliL3TrackArray *LookForPeaks(AliL3Histogram *hist,Int_t nbins);
  
  void FindPeak(Int_t t1,Double_t t2,Int_t t3,Float_t &kappa,Float_t &phi0);
  AliL3HoughTrack *FindPeakLine(Double_t rho,Double_t theta);
  AliL3HoughTrack *CalculatePeakInWindow(Int_t *maxbin,Int_t t0,Int_t t1,Double_t t2,Int_t t3);
  void FindPeak1(Float_t *xpeaks,Float_t *ypeaks,Int_t *weight,Int_t &n,Int_t y_window=2,Int_t x_bin_sides=1);
  void SortPeaks(struct AxisWindow **a,Int_t first,Int_t last);
  Int_t PeakCompare(struct AxisWindow *a,struct AxisWindow *b);
  
  void SetThreshold(Int_t f) {fThreshold = f;}
  
  void SetHistogram(AliL3Histogram *hist) {fCurrentHisto = hist;}
  
  ClassDef(AliL3HoughMaxFinder,1) //Maximum finder class

};

#endif
