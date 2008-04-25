#ifndef AliPHOSDA1_H
#define AliPHOSDA1_H

#include "TNamed.h"
#include "TH1.h"
#include "TH2F.h"
#include "TFile.h"
#include "TObjArray.h"

class AliPHOSDA1 : public TNamed {
  
 public:
  
  AliPHOSDA1(Int_t module);
  AliPHOSDA1(Int_t module, TH2F oldTimeEnergy[64][56][2]);
  AliPHOSDA1(const AliPHOSDA1& );
  AliPHOSDA1& operator= (const AliPHOSDA1& );
  ~AliPHOSDA1();
  
  void  FillHistograms(Float_t e[64][56][2], Float_t t[64][56][2]);
  Int_t GetModule() { return fMod; }
  void  UpdateHistoFile();
  void  SetWriteToFile(Bool_t write);

  const TH2F* GetTimeEnergyHistogram(Int_t X, Int_t Z, Int_t gain) const 
  { return fTimeEnergy[X][Z][gain]; }
  const TH1F* GetHgLgRatioHistogram(Int_t X, Int_t Z) const
  { return fHgLgRatio[X][Z]; }

  const TObjArray* GetHistoContainer() const { return &fHistoArray; }
   
 private:

  TFile* fHistoFile;            // root file to store histograms in
  TH1F* fHgLgRatio[64][56];     // high gain to low gain ratio  
  TH2F* fTimeEnergy[64][56][2]; // time and energy
  Int_t fMod;                   // PHOS module number (0..4)
  Bool_t fWriteToFile;          // kTRUE to save histograms to ROOT file (default) 
  TObjArray fHistoArray;        // container for histograms
  
  ClassDef(AliPHOSDA1,1)

};

#endif
