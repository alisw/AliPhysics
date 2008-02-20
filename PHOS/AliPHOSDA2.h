#ifndef AliPHOSDA2_H
#define AliPHOSDA2_H

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

class AliPHOSDA2 : public TNamed {
  
 public:
  
  AliPHOSDA2(Int_t module);
  AliPHOSDA2(const AliPHOSDA2& );
  AliPHOSDA2& operator= (const AliPHOSDA2& );
  ~AliPHOSDA2();
  
  void  FillQualityHistograms(Float_t quality[64][56][2]);
  Int_t GetModule() { return fMod; }
  void  UpdateHistoFile();
  
 private:

  TFile* fHistoFile;            // root file to store histograms in
  TH1F* fHQuality[64][56][2];   // "quality" for high and low gains
  Int_t fMod;                   // PHOS module number (0..4)
  TH2F* fMaps[2];               // 2D quality map for low and high gains.
  
  ClassDef(AliPHOSDA2,1)

};

#endif
