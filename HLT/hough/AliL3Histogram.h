#ifndef ALIL3_HISTOGRAM
#define ALIL3_HISTOGRAM

#include "AliL3RootTypes.h"
#include <TH2.h>


class AliL3Histogram {
  
 private:
  
  Double_t *fContent; //!
  Char_t fName[100];
  Int_t fNxbins;
  Int_t fNybins;
  Int_t fNcells;
  Int_t fEntries;
  Int_t fFirstXbin;
  Int_t fFirstYbin;
  Int_t fLastXbin;
  Int_t fLastYbin;
  Int_t fThreshold;

  Double_t fXmin;
  Double_t fYmin;
  Double_t fXmax;
  Double_t fYmax;
  
  TH2F *fRootHisto;
  
 public:
  AliL3Histogram();
  AliL3Histogram(Char_t *name,Char_t *id,Int_t nxbin,Double_t xmin,Double_t xmax,Int_t nybin,Double_t ymin,Double_t ymax);
  virtual ~AliL3Histogram();
  
  void Reset();
  void Fill(Double_t x,Double_t y,Int_t weight=1);
  Int_t FindBin(Double_t x,Double_t y);
  Int_t FindXbin(Double_t x);
  Int_t FindYbin(Double_t y);
  Int_t GetBin(Int_t xbin,Int_t ybin);
  Double_t GetBinContent(Int_t bin);
  void SetBinContent(Int_t xbin,Int_t ybin,Int_t value);
  void SetBinContent(Int_t bin,Int_t value);
  void AddBinContent(Int_t xbin,Int_t ybin,Int_t weight);
  void AddBinContent(Int_t bin,Int_t weight);
  void Add(AliL3Histogram *h1,Double_t weight=1);
  void Draw(Char_t *option="hist");
  void SetThreshold(Int_t i) {fThreshold = i;}

  TH2F *GetRootHisto() {return fRootHisto;}
  Double_t GetXmin() {return fXmin;}
  Double_t GetXmax() {return fXmax;}
  Double_t GetYmin() {return fYmin;}
  Double_t GetYmax() {return fYmax;}
  Double_t GetBinCenterX(Int_t xbin);
  Double_t GetBinCenterY(Int_t ybin);
  Int_t GetFirstXbin() {return fFirstXbin;}
  Int_t GetLastXbin() {return fLastXbin;}
  Int_t GetFirstYbin() {return fFirstYbin;}
  Int_t GetLastYbin() {return fLastYbin;}
  Int_t GetNbinsX() {return fNxbins;}
  Int_t GetNbinsY() {return fNybins;}
  Int_t GetNEntries() {return fEntries;}
  
  
  ClassDef(AliL3Histogram,1)
    
};

#endif
