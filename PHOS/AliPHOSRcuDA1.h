#ifndef AliPHOSRCUDA1_H
#define AliPHOSRCUDA1_H

#include "TNamed.h"
#include "TH1.h"
#include "TH2F.h"
#include "TFile.h"
#include "TObjArray.h"

class AliPHOSRcuDA1 : public TNamed {
  
 public:
  
  AliPHOSRcuDA1(Int_t module, Int_t rcu);
  AliPHOSRcuDA1(Int_t module, Int_t rcu, TObjArray* oldTimeEnergy);
  ~AliPHOSRcuDA1();
  
  void  FillHistograms(Float_t e[64][56][2], Float_t t[64][56][2]);
  Int_t GetModule() { return fMod; }
  Int_t GetRCU() { return fRCU; }
  void  UpdateHistoFile();
  void  SetWriteToFile(Bool_t write);

  const TH2F* GetTimeEnergyHistogram(Int_t X, Int_t Z, Int_t gain) const 
  { return fTimeEnergy[X][Z][gain]; }
  const TH1F* GetHgLgRatioHistogram(Int_t X, Int_t Z) const
  { return fHgLgRatio[X][Z]; }

  const TObjArray* GetHistoContainer() const { return &fHistoArray; }

 private:

  AliPHOSRcuDA1(const AliPHOSRcuDA1& );
  AliPHOSRcuDA1& operator= (const AliPHOSRcuDA1& );
  
 private:

  TFile* fHistoFile;            // root file to store histograms in
  TH1F* fHgLgRatio[64][56];     // high gain to low gain ratio  
  TH2F* fTimeEnergy[64][56][2]; // time and energy
  Int_t fMod;                   // PHOS module number (0..4)
  Int_t fRCU;                   // RCU number (0..3)
  Bool_t fWriteToFile;          // kTRUE to save histograms to ROOT file (default) 
  TObjArray fHistoArray;        // container for histograms
  
  ClassDef(AliPHOSRcuDA1,1)

};

#endif
