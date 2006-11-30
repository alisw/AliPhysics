// @(#) $Id$

#ifndef ALIL3HISTOGRAMADAPTIVE_H
#define ALIL3HISTOGRAMADAPTIVE_H

#include "AliHLTRootTypes.h"
#include "AliHLTHistogram.h"

class AliHLTHistogramAdaptive : public AliHLTHistogram {
  
 public:
  AliHLTHistogramAdaptive();
  AliHLTHistogramAdaptive(Char_t *name,Double_t minpt,Double_t maxpt,Double_t ptres,
			 Int_t nybins,Double_t ymin,Double_t ymax);
  ~AliHLTHistogramAdaptive();

  void Fill(Double_t x,Double_t y,Int_t weight=1);
  void Fill(Double_t x,Int_t ybin,Int_t weight=1) {
    AliHLTHistogram::Fill(x,ybin,weight);
  }
  void Fill(Int_t xbin,Double_t y,Int_t weight=1) {
    AliHLTHistogram::Fill(xbin,y,weight);
  }
  void Fill(Int_t xbin,Int_t ybin,Int_t weight=1) {
    AliHLTHistogram::Fill(xbin,ybin,weight);
  }
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
  
  ClassDef(AliHLTHistogramAdaptive,1) //2D histogram class
    
};

typedef AliHLTHistogramAdaptive AliL3HistogramAdaptive; // for backward comaatibility 

#endif
