#ifndef ALIL3_HISTOGRAMADAPTIVE
#define ALIL3_HISTOGRAMADAPTIVE

#include <stream.h>
#include "AliL3RootTypes.h"
#include "AliL3Histogram.h"

class AliL3HistogramAdaptive : public AliL3Histogram {
  
 private:
  Double_t *fKvalue; //!
  Double_t fMinPt;
  Double_t fMaxPt;
  Float_t fPtstep;
  
  Int_t InitPtBins(Int_t firstpatch,Int_t lastpatch,Double_t xyresolution);
  
 public:
  AliL3HistogramAdaptive();
  AliL3HistogramAdaptive(Char_t *name,Int_t firstpatch,Int_t lastpatch,Double_t minpt,Double_t maxpt,
			 Int_t nybins,Double_t ymin,Double_t ymax);
  ~AliL3HistogramAdaptive();

  void Fill(Double_t x,Double_t y,Int_t weight=1);
  Int_t FindBin(Double_t x,Double_t y);
  Int_t FindXbin(Double_t x);
  Int_t FindYbin(Double_t y);
  void Draw(Char_t *option = "hist");
  void Print();

  Double_t GetBinCenterX(Int_t xbin);
  Double_t GetBinCenterY(Int_t ybin);
  
  ClassDef(AliL3HistogramAdaptive,1) //2D histogram class
    
};

#endif
