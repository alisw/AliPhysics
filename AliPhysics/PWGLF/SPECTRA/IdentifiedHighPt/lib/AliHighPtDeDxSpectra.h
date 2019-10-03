#ifndef ALIHIGHPTDEDXSPECTRA_H
#define ALIHIGHPTDEDXSPECTRA_H

#include <TNamed.h>
#include <TH1.h>

#include "AliHighPtDeDxData.h"
#include "AliHighPtDeDxMc.h"

class AliHighPtDeDxSpectra : public TNamed {
 public:
  AliHighPtDeDxSpectra(); // default constructor  
  AliHighPtDeDxSpectra(const char* name, const char* title); // named constructor  
  virtual ~AliHighPtDeDxSpectra(); // default destructor

  void GetCorrectedSpectra(AliHighPtDeDxData* data,
			   AliHighPtDeDxMc* mc);

  void SetUseMcNoVtxCorrection(Bool_t value) { fUseMcNoVtxCorrection = value; };
  void SetUseFittedEfficiency(Bool_t value) { fUseFittedEfficiency = value; };
  void SetUseBinCorrection(Bool_t value) { fUseBinCorrection = value; };
  void SetDebugLevel(Int_t value)  { fDebugLevel = value; };

  Double_t GetNevents()        { return fNevents; };
  TH1D* GetHistPt()            { return hPt; };
  TH1D* GetHistNevents()       { return hNevents; };
  TH1D* GetHistTriggerEfficiency() { return hTriggerEfficiency; };
  TH1D* GetHistEfficiency()    { return hEfficiency; };
  TH1D* GetHistBinCorrection() { return hBinCorrection; };
  TH1D* GetHistMcNoVtxCorrection() { return hMcNoVtxCorrection; };
  TH1D* GetHistMeanPt()        { return hMeanPt; };
  TF1*  GetFuncEfficiency()    { return fEfficiency; };
  TF1*  GetFuncBinFit()        { return fBinFit; };


 private:

  TH1D* ConstructBinCorrection(TH1D* histPt, TProfile* histMeanPt);
  TH1D* ConstructTriggerEfficiency(AliHighPtDeDxMc* mc);
  TH1D* GetEventCount(AliHighPtDeDxData* data, AliHighPtDeDxMc* mc);
  TH1D* ConstructTrackCorrection(AliHighPtDeDxMc* mc);
  TF1*  FitEfficiency(TH1D* histEff);
  void NormalizePt(TH1* hist);

  // members
  Int_t    fDebugLevel;     // debug level
  Bool_t   fUseMcNoVtxCorrection; 
  Bool_t   fUseFittedEfficiency; 
  Bool_t   fUseBinCorrection; 
  // use MC to correct down the novtx fraction outside the vtx cut
  // see AliHighPtDeDxMc::FillEventInfo() for details
  Double_t fNevents;        // corrected number of events

  // histograms
  TH1D* hPt;                // pt spectrum for 
  TH1D* hNevents;           // events
  TH1D* hTriggerEfficiency; // trigger efficency
  TH1D* hEfficiency;        // track correction
  TH1D* hMcNoVtxCorrection; // mc no vtx correction
  TH1D* hBinCorrection;     // bin correction
  TH1D* hMeanPt;            // mean pt of data
  TF1*  fEfficiency;        // fitted efficiency
  TF1*  fBinFit;            // fit used for bin correction
  
  ClassDef(AliHighPtDeDxSpectra, 1)  // AliHighPtDeDxSpectra information
    };

#endif
	
