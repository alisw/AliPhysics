#ifndef ALIL3_HISTOGRAM
#define ALIL3_HISTOGRAM

#include "AliL3RootTypes.h"
#include <TH2.h>

class AliL3Histogram : public TH2F {
  
 private:
  
  Double_t *fContent; //!
  Char_t fName[100];
  Int_t fNxbins;
  Int_t fNybins;
  Int_t fNcells;
  Int_t fEntries;
  
  Double_t fXmin;
  Double_t fYmin;
  Double_t fXmax;
  Double_t fYmax;
  
  
 public:
  AliL3Histogram();
  AliL3Histogram(Char_t *name,Char_t *id,Int_t nxbin,Double_t xmin,Double_t xmax,Int_t nybin,Double_t ymin,Double_t ymax);
  virtual ~AliL3Histogram();
  
  void Reset();
  void Fill(Double_t x,Double_t y,Int_t weight);
  Int_t FindBin(Double_t x,Double_t y);
  void AddBinContent(Int_t xbin,Int_t ybin,Int_t weight);
  void AddBinContent(Int_t bin,Int_t weight);
  void Draw();

  Double_t GetXmin() {return fXmin;}
  Double_t GetXmax() {return fXmax;}
  Double_t GetYmin() {return fYmin;}
  Double_t GetYmax() {return fXmax;}
  Double_t GetXBinCenter(Int_t xbin);
  Double_t GetYBinCenter(Int_t ybin);
  Int_t GetFirstXbin() {return 1 + (Int_t)(fNxbins*(fXmin-fXmin)/(fXmax-fXmin));}
  Int_t GetLastXbin() {return 1 + (Int_t)(fNxbins*(fXmax-fXmin)/(fXmax-fXmin));}
  Int_t GetFirstYbin() {return 1 + (Int_t)(fNxbins*(fXmin-fXmin)/(fXmax-fXmin));}
  Int_t GetLastYbin() {return 1 + (Int_t)(fNybins*(fYmax-fYmin)/(fYmax-fYmin));}
  Int_t GetNEntries() {return fEntries;}

  ClassDef(AliL3Histogram,1)
    
};

#endif
