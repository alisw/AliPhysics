// @(#) $Id$

#ifndef ALIL3HISTOGRAM1D_H
#define ALIL3HISTOGRAM1D_H

#include "AliHLTRootTypes.h"

#ifdef use_root
class TH1F;
#endif

class AliHLTHistogram1D {
  
 public:
  AliHLTHistogram1D();
  AliHLTHistogram1D(Char_t *name,Char_t *id,Int_t nxbin,Double_t xmin,Double_t xmax);
  virtual ~AliHLTHistogram1D();
  
  void Reset();
  void Fill(Double_t x,Int_t weight=1);
  void AddBinContent(Int_t bin,Int_t weight);
  Int_t GetMaximumBin() const;
  Int_t FindBin(Double_t x) const;
  Double_t GetBinContent(Int_t bin) const;
  Double_t GetBinCenter(Int_t bin) const;
  Int_t GetNEntries() const {return fEntries;}
  
  void SetBinContent(Int_t bin,Int_t value);
  void SetThreshold(Int_t i) {fThreshold = i;}
  

#ifdef use_root
  void Draw(Char_t *option="hist");
  TH1F *GetRootHisto() {return fRootHisto;}
#endif
  
 private:
  
  Double_t *fContent; //!
  Char_t fName[100];//Histogram title
  Int_t fNbins;//Number of bins
  Int_t fNcells;//Number of cells
  Int_t fEntries;//Number of entries

  Int_t fThreshold;//Bin content threshold
  Double_t fXmin;//Lower limit in X
  Double_t fXmax;//Upper limit in X

  
#ifdef use_root
  TH1F *fRootHisto;//The corresponding ROOT histogram
#endif  

  ClassDef(AliHLTHistogram1D,1) //1D histogram class
    
};

typedef AliHLTHistogram1D AliL3Histogram1D; // for backward compatibility

#endif
