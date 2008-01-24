#ifndef AliPHOSDA1_H
#define AliPHOSDA1_H

#include "TNamed.h"
#include "TH1.h"
#include "TH2F.h"
#include "TFile.h"

class AliPHOSDA1 : public TNamed {
  
 public:
  
  AliPHOSDA1(Int_t module);
  AliPHOSDA1(const AliPHOSDA1& );
  AliPHOSDA1& operator= (const AliPHOSDA1& );
  ~AliPHOSDA1();
  
  void  FillHistograms(Float_t e[64][56][2], Float_t t[64][56][2]);
  Int_t GetModule() { return fMod; }
  void  UpdateHistoFile();
  
 private:

  TFile* fHistoFile;            // root file to store histograms in
  TH1F* fHgLgRatio[64][56];     // high gain to low gain ratio  
  TH2F* fTimeEnergy[64][56][2]; // time and energy
  Int_t fMod;                   // PHOS module number (0..4)
  
  ClassDef(AliPHOSDA1,1)

};

#endif
