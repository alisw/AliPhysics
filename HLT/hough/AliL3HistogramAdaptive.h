// @(#) $Id$

#ifndef ALIL3_HISTOGRAMADAPTIVE
#define ALIL3_HISTOGRAMADAPTIVE

#include "AliL3RootTypes.h"
#include "AliL3Histogram.h"

class AliL3HistogramAdaptive : public AliL3Histogram {
  
 private:
  Double_t fPtres;
  Double_t fMinPt;
  Double_t fMaxPt;
  Double_t *fKappaBins; //!
  
  Int_t InitKappaBins();
  
 public:
  AliL3HistogramAdaptive();
  AliL3HistogramAdaptive(Char_t *name,Double_t minpt,Double_t maxpt,Double_t ptres,
			 Int_t nybins,Double_t ymin,Double_t ymax);
  ~AliL3HistogramAdaptive();

  void Fill(Double_t x,Double_t y,Int_t weight=1);
  Int_t FindBin(Double_t x,Double_t y) const;
  Int_t FindXbin(Double_t x) const;
  Int_t FindYbin(Double_t x) const;
  void Draw(Char_t *option = "hist");
  void Print();

  Double_t GetBinCenterX(Int_t xbin) const;
  Double_t GetBinCenterY(Int_t ybin) const;

  ClassDef(AliL3HistogramAdaptive,1) //2D histogram class
    
};

#endif
