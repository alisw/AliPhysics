#ifndef ALIHIGHPTDEDXCALIB_H
#define ALIHIGHPTDEDXCALIB_H

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include "AliHighPtDeDxBase.h"

class AliHighPtDeDxCalib : public AliHighPtDeDxBase {
 public:
  AliHighPtDeDxCalib(); // default constructor  
  AliHighPtDeDxCalib(const char* name, const char* title); // named constructor  
  virtual ~AliHighPtDeDxCalib(); // default destructor

  virtual void Init(Int_t nPtBins, Double_t* ptBins);
  virtual void Init(Int_t step, Int_t nPtBins, Double_t* ptBins);
  virtual void FillTrackInfo(Float_t weight);
  virtual void PerformEtaCal();
  virtual void PerformNclCal();

  TCanvas* DrawNclCal();
  TCanvas* DrawEta(Bool_t forMIP);
  TCanvas* DrawEtaCalibrated(Bool_t forMIP);
  TCanvas* DrawSelectionHistograms(Int_t step=3);
  
  virtual void SetStep(Int_t value) { fStep = value; }

  virtual void SetPMIPMin(Double_t value)    { fPMIPMin    = value; }
  virtual void SetPMIPMax(Double_t value)    { fPMIPMax    = value; }
  virtual void SetDeDxMIPMin(Double_t value) { fDeDxMIPMin = value; }
  virtual void SetDeDxMIPMax(Double_t value) { fDeDxMIPMax = value; }
  virtual void SetDeltaBeta(Double_t value)  { fDeltaBeta  = value; }

  TF1*  GetDeDxVsEtaNeg()  {return fDeDxVsEtaNeg;}
  TF1*  GetDeDxVsEtaPos()  {return fDeDxVsEtaPos;} 
  TF1*  GetDeDxVsNcl()     {return fDeDxVsNcl;}
  TProfile* GetHistMeanP() {return hMeanP;}
  TH1D* GetHistDeDx(Bool_t forMIP, Int_t etaBin);
  TH2D* GetHistDeDxVsNcl(Bool_t forMIP, Int_t etaBin);
  TH2D* GetHistDeDxVsP(Int_t pid);

  //  virtual void Init(Int_t nPtBins, Double_t* ptBins);
  //  virtual void FillTrackInfo(Double_t weight=1);
  
 private:

  Bool_t IsMIP();
  Bool_t IsElectron();

  Int_t fStep;             // Step 1 = Eta calibration, step 2 = dE/dx calibration
  Int_t fInit;             // Step 1 = Eta calibration, step 2 = dE/dx calibration
  
  Double_t fPMIPMin;       // Min P for MIP pion
  Double_t fPMIPMax;       // Max P for MIP pion
  Double_t fDeDxMIPMin;    // Min dE/dx for MIP pion
  Double_t fDeDxMIPMax;    // Max dE/dx for MIP pion
  Double_t fDeltaBeta;     // delta beta cut for electrons 

  TF1*   fDeDxPi;          // dE/dx vs p for pions
  TF1*   fSigmaDeDx;       // sigma dE/dx vs ncl

  // functions
  TF1* fDeDxVsEtaNeg;      // eta < 0 dE/dx calib
  TF1* fDeDxVsEtaPos;      // eta > 0 dE/dx calib
  TF1* fDeDxVsNcl;         // ncl dE/dx calib
  
  // histograms - step 0
  TH2D* hSelection1;            // selected region in p and dE/dx for pion MIPs
  TH2D* hSelection2;            // selected region in p and dE/dx for pion MIPs
  TH2D* hSelection3;            // selected region in p and dE/dx for pion MIPs
  TH2D* hSelectionElectrons2;   // selected region in p and dE/dx for electrons
  TH2D* hSelectionElectrons3;   // selected region in p and dE/dx for electrons
  TH2D* hDeDxVsEta;             // dE/dx vs eta uncalibrated (pion MIPs)
  TH2D* hDeDxVsEtaElectrons;    // dE/dx vs eta uncalibrated (electrons)
  TH2D* hNclVsEta;              // Ncl vs eta (pion MIPs)
  TH2D* hNclVsEtaElectrons;     // Ncl vs eta (electrons)

  // histograms - step 1
  TH2D* hDeDxVsEtaCal;          // dE/dx vs eta calibrated (pion MIPs)
  TH2D* hDeDxVsEtaCalElectrons; // dE/dx vs eta calibrated (electrons)
  TProfile* hMeanEta;           // <eta> in the 4 eta interval (pion MIPs)
  TProfile* hMeanEtaElectrons;  // <eta> in the 4 eta interval (electrons)

  TH1D* hDeDx;                  // dE/dx no eta cut    (pion MIPs)
  TH1D* hDeDx1;                 // dE/dx 0.0<|eta|<0.2 (pion MIPs)
  TH1D* hDeDx2;                 // dE/dx 0.2<|eta|<0.4 (pion MIPs)
  TH1D* hDeDx3;                 // dE/dx 0.4<|eta|<0.6 (pion MIPs)
  TH1D* hDeDx4;                 // dE/dx 0.6<|eta|<0.8 (pion MIPs)
  TH1D* hDeDxElectrons;         // dE/dx no eta cut    (electrons)
  TH1D* hDeDxElectrons1;        // dE/dx 0.0<|eta|<0.2 (electrons)
  TH1D* hDeDxElectrons2;        // dE/dx 0.2<|eta|<0.4 (electrons)
  TH1D* hDeDxElectrons3;        // dE/dx 0.4<|eta|<0.6 (electrons)
  TH1D* hDeDxElectrons4;        // dE/dx 0.6<|eta|<0.8 (electrons)
  
  TH2D* hDeDxVsNclBefore;       // dE/dx vs ncl for calib (step 2)

  TH2D* hDeDxVsNcl;             // dE/dx vs ncl no eta cut    (pion MIPs)
  TH2D* hDeDxVsNcl1;            // dE/dx vs ncl 0.0<|eta|<0.2 (pion MIPs)
  TH2D* hDeDxVsNcl2;            // dE/dx vs ncl 0.2<|eta|<0.4 (pion MIPs)
  TH2D* hDeDxVsNcl3;            // dE/dx vs ncl 0.4<|eta|<0.6 (pion MIPs)
  TH2D* hDeDxVsNcl4;            // dE/dx vs ncl 0.6<|eta|<0.8 (pion MIPs)
  TH2D* hDeDxVsNclElectrons;    // dE/dx vs ncl no eta cut    (electrons)
  TH2D* hDeDxVsNclElectrons1;   // dE/dx vs ncl 0.0<|eta|<0.2 (electrons)
  TH2D* hDeDxVsNclElectrons2;   // dE/dx vs ncl 0.2<|eta|<0.4 (electrons)
  TH2D* hDeDxVsNclElectrons3;   // dE/dx vs ncl 0.4<|eta|<0.6 (electrons)
  TH2D* hDeDxVsNclElectrons4;   // dE/dx vs ncl 0.6<|eta|<0.8 (electrons)
  
  /* TH1D* GetHistDeDx(Int_t bin = 0) { */
    

  /* } */

  // histograms - step 1 dE/dx

  TProfile* hMeanP;    // <p> vs p

  TH2D* hDeDxVsP;      // dE/dx vs p 
  TH2D* hDeDxVsPPiMc;  // dE/dx vs p for MC pions
  TH2D* hDeDxVsPKMc;   // dE/dx vs p for MC Kaons
  TH2D* hDeDxVsPPMc;   // dE/dx vs p for MC protons
  TH2D* hDeDxVsPEMc;   // dE/dx vs p for MC electrons
  TH2D* hDeDxVsPMuMc;  // dE/dx vs p for MC muons
  
  //  void Print(Option_t* option) const;
  
  ClassDef(AliHighPtDeDxCalib, 1)  // AliHighPtDeDxCalib information
    };

#endif
	
