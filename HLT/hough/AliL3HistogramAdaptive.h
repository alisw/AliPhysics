// @(#) $Id$

#ifndef ALIL3HISTOGRAMADAPTIVE_H
#define ALIL3HISTOGRAMADAPTIVE_H

#include "AliL3RootTypes.h"
#include "AliL3Histogram.h"

class AliL3HistogramAdaptive : public AliL3Histogram {
  
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
  void Print() const;

  Double_t GetBinCenterX(Int_t xbin) const;
  Double_t GetBinCenterY(Int_t ybin) const;

 private:
  Double_t fPtres;//The desired Pt resolution
  Double_t fMinPt;//Minimum Pt
  Double_t fMaxPt;//Maximum Pt
  Double_t *fKappaBins; //!
  
  Int_t InitKappaBins();
  
  ClassDef(AliL3HistogramAdaptive,1) //2D histogram class
    
};

#endif
