// @(#) $Id$

#ifndef ALIL3_HISTOGRAM1D
#define ALIL3_HISTOGRAM1D

#include "AliL3RootTypes.h"

#ifdef use_root
#include <TH1.h>
#endif

class AliL3Histogram1D {
  
 private:
  
  Double_t *fContent; //!
  Char_t fName[100];
  Int_t fNbins;
  Int_t fNcells;
  Int_t fEntries;

  Int_t fThreshold;
  Double_t fXmin;
  Double_t fXmax;

  
#ifdef use_root
  TH1F *fRootHisto;
#endif  

 public:
  AliL3Histogram1D();
  AliL3Histogram1D(Char_t *name,Char_t *id,Int_t nxbin,Double_t xmin,Double_t xmax);
  virtual ~AliL3Histogram1D();
  
  void Reset();
  void Fill(Double_t x,Int_t weight=1);
  void AddBinContent(Int_t bin,Int_t weight);
  Int_t GetMaximumBin();
  Int_t FindBin(Double_t x);
  Double_t GetBinContent(Int_t bin);
  Double_t GetBinCenter(Int_t bin);
  Int_t GetNEntries() {return fEntries;}
  
  void SetBinContent(Int_t bin,Int_t value);
  void SetThreshold(Int_t i) {fThreshold = i;}
  

#ifdef use_root
  void Draw(Char_t *option="hist");
  TH1F *GetRootHisto() {return fRootHisto;}
#endif
  
  ClassDef(AliL3Histogram1D,1) //1D histogram class
    
};

#endif
