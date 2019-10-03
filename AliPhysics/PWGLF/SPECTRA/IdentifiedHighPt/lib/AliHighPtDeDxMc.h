#ifndef ALIHIGHPTDEDXMC_H
#define ALIHIGHPTDEDXMC_H

#include "AliHighPtDeDxBase.h"

class AliHighPtDeDxMc : public AliHighPtDeDxBase {
 public:
  AliHighPtDeDxMc(); // default constructor  
  AliHighPtDeDxMc(const char* name, const char* title); // named constructor  
  virtual ~AliHighPtDeDxMc(); // default destructor
  
  TH1D* GetPtSpectrum(); 

  virtual void SetTrackChargeMc(Int_t value)  { fTrackChargeMc = value; }
  virtual void SetTrackEtaMc(Double_t value)  { fTrackEtaMc = value; }
  virtual void SetTrackPtMc(Double_t value)   { fTrackPtMc = value; } 
  
  virtual void Init(Int_t nPtBins, Double_t* ptBins);
  virtual Bool_t TrackAcceptedMc();
  virtual void FillEventInfo();
  virtual void FillTrackInfoMc(Float_t weight=1);
  
  TH1D* GetHistNeventsMc()     { return hNeventsMc; };
  TH1D* GetHistNeventsMcTrig() { return hNeventsMcTrig; };
  TH1D* GetHistPtMc(Int_t pid, Int_t charge);
  TProfile* GetHistMeanPtMc()  { return hMeanPtMc; };

 protected:
  // Actual values for the event and track
  Int_t    fTrackChargeMc;       //! charge (+1 or -1)
  Double_t fTrackEtaMc;          //! eta
  Double_t fTrackPtMc;           //! pt

 private:
  void NormalizePt(TH1* hist);

  // histograms
  TH1D*     hVtxStatusMc;      // vtx status (all triggers) 
  TH1D*     hNeventsMc;        // Nevents (all triggers) - based on MC vtx 
  TH1D*     hNeventsMcTrig;    // Nevents (all triggers) - based on MC vtx + trig
  TH1D*     hPtMc;             // pt input distribution
  TH1D*     hPtMcNeg;          // pt input distribution (q<0)
  TH1D*     hPtMcPos;          // pt input distribution (q>0)
  TH1D*     hPtPiMc;           // pt input distribution
  TH1D*     hPtPiMcNeg;        // pt input distribution (q<0)
  TH1D*     hPtPiMcPos;        // pt input distribution (q>0)
  TH1D*     hPtKMc;            // pt input distribution
  TH1D*     hPtKMcNeg;         // pt input distribution (q<0)
  TH1D*     hPtKMcPos;         // pt input distribution (q>0)
  TH1D*     hPtPMc;            // pt input distribution
  TH1D*     hPtPMcNeg;         // pt input distribution (q<0)
  TH1D*     hPtPMcPos;         // pt input distribution (q>0)
  TProfile* hMeanPtMc;         // mean pt input distribution
    
  ClassDef(AliHighPtDeDxMc, 4)  // AliHighPtDeDxMc information
    };

#endif
	
