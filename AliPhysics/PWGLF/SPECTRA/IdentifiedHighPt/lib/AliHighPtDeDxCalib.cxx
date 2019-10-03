#include "AliHighPtDeDxCalib.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif

using namespace std;

ClassImp(AliHighPtDeDxCalib);

//
// AliHighPtDeDxCalib class
//
// This class contains the AliHighPtDeDxCalib information 
//

//_________________________________________________________
AliHighPtDeDxCalib::AliHighPtDeDxCalib():
  AliHighPtDeDxBase(),
  fStep(1),             // Step 1 = Eta calibration, step 2 = dE/dx calibration, step = 3 ncl calib
  fInit(0),             // Step 1 = Eta calibration, step 2 = dE/dx calibration, step = 3 ncl calib
  fPMIPMin(0.4),        // Min P for MIP pion
  fPMIPMax(0.6),        // Max P for MIP pion
  fDeDxMIPMin(40),      // Min dE/dx for MIP pion
  fDeDxMIPMax(60),      // Max dE/dx for MIP pion
  fDeltaBeta(0.1),      // delta beta cut for electrons 
  fDeDxPi(0x0),         // dE/dx vs p for pions
  fSigmaDeDx(0x0),      // sigma dE/dx vs ncl
  fDeDxVsEtaNeg(0x0),   // eta < 0 dE/dx calib
  fDeDxVsEtaPos(0x0),   // eta > 0 dE/dx calib
  fDeDxVsNcl(0x0),      // ncl dE/dx calib
  hSelection1(0x0),             // selected region in p and dE/dx for pion MIPs
  hSelection2(0x0),             // selected region in p and dE/dx for pion MIPs
  hSelection3(0x0),             // selected region in p and dE/dx for pion MIPs
  hSelectionElectrons2(0x0),    // selected region in p and dE/dx for electrons
  hSelectionElectrons3(0x0),    // selected region in p and dE/dx for electrons
  hDeDxVsEta(0x0),             // dE/dx vs eta uncalibrated (pion MIPs)
  hDeDxVsEtaElectrons(0x0),    // dE/dx vs eta uncalibrated (electrons)
  hNclVsEta(0x0),              // Ncl vs eta (pion MIPs)
  hNclVsEtaElectrons(0x0),     // Ncl vs eta (electrons)
  hDeDxVsEtaCal(0x0),          // dE/dx vs eta calibrated (pion MIPs)
  hDeDxVsEtaCalElectrons(0x0), // dE/dx vs eta calibrated (electrons)
  hMeanEta(0x0),               // <eta> in the 4 eta interval (pion MIPs)
  hMeanEtaElectrons(0x0),      // <eta> in the 4 eta interval (electrons)
  hDeDx(0x0),                  // dE/dx no eta cut    (pion MIPs)
  hDeDx1(0x0),                 // dE/dx 0.0<|eta|<0.2 (pion MIPs)
  hDeDx2(0x0),                 // dE/dx 0.2<|eta|<0.4 (pion MIPs)
  hDeDx3(0x0),                 // dE/dx 0.4<|eta|<0.6 (pion MIPs)
  hDeDx4(0x0),                 // dE/dx 0.6<|eta|<0.8 (pion MIPs)
  hDeDxElectrons(0x0),         // dE/dx no eta cut    (electrons)
  hDeDxElectrons1(0x0),        // dE/dx 0.0<|eta|<0.2 (electrons)
  hDeDxElectrons2(0x0),        // dE/dx 0.2<|eta|<0.4 (electrons)
  hDeDxElectrons3(0x0),        // dE/dx 0.4<|eta|<0.6 (electrons)
  hDeDxElectrons4(0x0),        // dE/dx 0.6<|eta|<0.8 (electrons)
  hDeDxVsNclBefore(0x0),       // dE/dx vs ncl for step 2 calib (pion MIPs)
  hDeDxVsNcl(0x0),             // dE/dx vs ncl no eta cut    (pion MIPs)
  hDeDxVsNcl1(0x0),            // dE/dx vs ncl 0.0<|eta|<0.2 (pion MIPs)
  hDeDxVsNcl2(0x0),            // dE/dx vs ncl 0.2<|eta|<0.4 (pion MIPs)
  hDeDxVsNcl3(0x0),            // dE/dx vs ncl 0.4<|eta|<0.6 (pion MIPs)
  hDeDxVsNcl4(0x0),            // dE/dx vs ncl 0.6<|eta|<0.8 (pion MIPs)
  hDeDxVsNclElectrons(0x0),    // dE/dx vs ncl no eta cut    (electrons)
  hDeDxVsNclElectrons1(0x0),   // dE/dx vs ncl 0.0<|eta|<0.2 (electrons)
  hDeDxVsNclElectrons2(0x0),   // dE/dx vs ncl 0.2<|eta|<0.4 (electrons)
  hDeDxVsNclElectrons3(0x0),   // dE/dx vs ncl 0.4<|eta|<0.6 (electrons)
  hDeDxVsNclElectrons4(0x0),   // dE/dx vs ncl 0.6<|eta|<0.8 (electrons)
  hMeanP(0x0),        // <p> vs p
  hDeDxVsP(0x0),      // dE/dx vs p 
  hDeDxVsPPiMc(0x0),  // dE/dx vs p for MC pions
  hDeDxVsPKMc(0x0),   // dE/dx vs p for MC Kaons
  hDeDxVsPPMc(0x0),   // dE/dx vs p for MC protons
  hDeDxVsPEMc(0x0),   // dE/dx vs p for MC electrons
  hDeDxVsPMuMc(0x0)   // dE/dx vs p for MC muons
{
  // default constructor - do not use

}

//_________________________________________________________
AliHighPtDeDxCalib::AliHighPtDeDxCalib(const char* name, const char* title):
  AliHighPtDeDxBase(name, title),
  fStep(1),             // Step 1 = Eta calibration, step 2 = dE/dx calibration
  fInit(0),             // Step 1 = Eta calibration, step 2 = dE/dx calibration, step = 3 ncl calib
  fPMIPMin(0.4),        // Min P for MIP pion
  fPMIPMax(0.6),        // Max P for MIP pion
  fDeDxMIPMin(40),      // Min dE/dx for MIP pion
  fDeDxMIPMax(60),      // Max dE/dx for MIP pion
  fDeltaBeta(0.1),      // delta beta cut for electrons 
  fDeDxPi(0x0),         // dE/dx vs p for pions
  fSigmaDeDx(0x0),      // sigma dE/dx vs ncl
  fDeDxVsEtaNeg(0x0),   // eta < 0 dE/dx calib
  fDeDxVsEtaPos(0x0),   // eta > 0 dE/dx calib
  fDeDxVsNcl(0x0),      // ncl dE/dx calib
  hSelection1(0x0),             // selected region in p and dE/dx for pion MIPs
  hSelection2(0x0),             // selected region in p and dE/dx for pion MIPs
  hSelection3(0x0),             // selected region in p and dE/dx for pion MIPs
  hSelectionElectrons2(0x0),    // selected region in p and dE/dx for electrons
  hSelectionElectrons3(0x0),    // selected region in p and dE/dx for electrons
  hDeDxVsEta(0x0),             // dE/dx vs eta uncalibrated (pion MIPs)
  hDeDxVsEtaElectrons(0x0),    // dE/dx vs eta uncalibrated (electrons)
  hNclVsEta(0x0),              // Ncl vs eta (pion MIPs)
  hNclVsEtaElectrons(0x0),     // Ncl vs eta (electrons)
  hDeDxVsEtaCal(0x0),          // dE/dx vs eta calibrated (pion MIPs)
  hDeDxVsEtaCalElectrons(0x0), // dE/dx vs eta calibrated (electrons)
  hMeanEta(0x0),               // <eta> in the 4 eta interval (pion MIPs)
  hMeanEtaElectrons(0x0),      // <eta> in the 4 eta interval (electrons)
  hDeDx(0x0),                  // dE/dx no eta cut    (pion MIPs)
  hDeDx1(0x0),                 // dE/dx 0.0<|eta|<0.2 (pion MIPs)
  hDeDx2(0x0),                 // dE/dx 0.2<|eta|<0.4 (pion MIPs)
  hDeDx3(0x0),                 // dE/dx 0.4<|eta|<0.6 (pion MIPs)
  hDeDx4(0x0),                 // dE/dx 0.6<|eta|<0.8 (pion MIPs)
  hDeDxElectrons(0x0),         // dE/dx no eta cut    (electrons)
  hDeDxElectrons1(0x0),        // dE/dx 0.0<|eta|<0.2 (electrons)
  hDeDxElectrons2(0x0),        // dE/dx 0.2<|eta|<0.4 (electrons)
  hDeDxElectrons3(0x0),        // dE/dx 0.4<|eta|<0.6 (electrons)
  hDeDxElectrons4(0x0),        // dE/dx 0.6<|eta|<0.8 (electrons)
  hDeDxVsNclBefore(0x0),       // dE/dx vs ncl no eta cut    (pion MIPs)
  hDeDxVsNcl(0x0),             // dE/dx vs ncl no eta cut    (pion MIPs)
  hDeDxVsNcl1(0x0),            // dE/dx vs ncl 0.0<|eta|<0.2 (pion MIPs)
  hDeDxVsNcl2(0x0),            // dE/dx vs ncl 0.2<|eta|<0.4 (pion MIPs)
  hDeDxVsNcl3(0x0),            // dE/dx vs ncl 0.4<|eta|<0.6 (pion MIPs)
  hDeDxVsNcl4(0x0),            // dE/dx vs ncl 0.6<|eta|<0.8 (pion MIPs)
  hDeDxVsNclElectrons(0x0),    // dE/dx vs ncl no eta cut    (electrons)
  hDeDxVsNclElectrons1(0x0),   // dE/dx vs ncl 0.0<|eta|<0.2 (electrons)
  hDeDxVsNclElectrons2(0x0),   // dE/dx vs ncl 0.2<|eta|<0.4 (electrons)
  hDeDxVsNclElectrons3(0x0),   // dE/dx vs ncl 0.4<|eta|<0.6 (electrons)
  hDeDxVsNclElectrons4(0x0),   // dE/dx vs ncl 0.6<|eta|<0.8 (electrons)
  hMeanP(0x0),        // <p> vs p
  hDeDxVsP(0x0),      // dE/dx vs p 
  hDeDxVsPPiMc(0x0),  // dE/dx vs p for MC pions
  hDeDxVsPKMc(0x0),   // dE/dx vs p for MC Kaons
  hDeDxVsPPMc(0x0),   // dE/dx vs p for MC protons
  hDeDxVsPEMc(0x0),   // dE/dx vs p for MC electrons
  hDeDxVsPMuMc(0x0)   // dE/dx vs p for MC muons
{
  // named constructor
}

//_________________________________________________________
AliHighPtDeDxCalib::~AliHighPtDeDxCalib()
{
  delete fDeDxPi;
  delete fSigmaDeDx;
  delete fDeDxVsEtaNeg;
  delete fDeDxVsEtaPos;
  delete fDeDxVsNcl;
  delete hSelection1;
  delete hSelection2;
  delete hSelection3;
  delete hSelectionElectrons2;
  delete hSelectionElectrons3;
  delete hDeDxVsEta;
  delete hDeDxVsEtaElectrons;
  delete hNclVsEta;
  delete hNclVsEtaElectrons;
  delete hDeDxVsEtaCal;
  delete hDeDxVsEtaCalElectrons;
  delete hMeanEta;
  delete hMeanEtaElectrons;
  delete hDeDx;
  delete hDeDx1;
  delete hDeDx2;
  delete hDeDx3;
  delete hDeDx4;
  delete hDeDxElectrons;
  delete hDeDxElectrons1;
  delete hDeDxElectrons2;
  delete hDeDxElectrons3;
  delete hDeDxElectrons4;
  delete hDeDxVsNclBefore;
  delete hDeDxVsNcl;
  delete hDeDxVsNcl1;
  delete hDeDxVsNcl2;
  delete hDeDxVsNcl3;
  delete hDeDxVsNcl4;
  delete hDeDxVsNclElectrons;
  delete hDeDxVsNclElectrons1;
  delete hDeDxVsNclElectrons2;
  delete hDeDxVsNclElectrons3;
  delete hDeDxVsNclElectrons4;
  delete hMeanP;
  delete hDeDxVsP;
  delete hDeDxVsPPiMc;
  delete hDeDxVsPKMc;
  delete hDeDxVsPPMc;
  delete hDeDxVsPEMc;
  delete hDeDxVsPMuMc;
}

//_________________________________________________________
void AliHighPtDeDxCalib::Init(Int_t nPtBins, Double_t* ptBins)
{
  //
  // Create histograms and functions
  //

  //
  // init base class
  //
  AliHighPtDeDxBase::Init(nPtBins, ptBins);

  //
  // functions
  //
  fDeDxVsEtaNeg = new TF1("fDeDxVsEtaNeg", "pol3", -1,  0);
  fDeDxVsEtaNeg->SetParameters(50, 0, 0, 0);
  fDeDxVsEtaPos = new TF1("fDeDxVsEtaPos", "pol3",  0, +1);
  fDeDxVsEtaPos->SetParameters(50, 0, 0, 0);

  fDeDxVsNcl = new TF1("fDeDxVsNcl", "pol1",  0, 160);
  fDeDxVsNcl->SetParameters(50, 0);
}

void AliHighPtDeDxCalib::Init(Int_t step, Int_t nPtBins, Double_t* ptBins)
{
  if(step != fStep)
    cout << "Warning! Innit called for step " << step << " but next step is " << fStep << endl;

  if(fInit>=step) {
    cout << "Step " << step << " has been initialized!" << endl;
    return;
  }

  switch (step) {
  case 1:
    //
    // step 1 histograms - nCl calibration
    //
    hSelection1 = new TH2D("hSelection1", "dE/dx vs P (pion MIPs) step 1; P [GeV/c]; dE/dx",
			   100, 0.2, 0.8, 
			   100, 0, 100);
    hSelection1->SetDirectory(0);

    hDeDxVsNclBefore = new TH2D("hDeDxVsNclBefore", "dE/dx vs Ncl (pion MIPs); Ncl; dE/dx",
				18, 69.5, 159.5, 
				fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDxVsNclBefore->Sumw2();
    hDeDxVsNclBefore->SetDirectory(0);

    fInit = 1;
    
    break;

  case 2:
    //
    // step 2 histograms - eta calibration
    //
    
    hSelection2 = new TH2D("hSelection1", "dE/dx vs P (pion MIPs) step 2; P [GeV/c]; dE/dx",
			   100, 0.2, 0.8, 
			   100, 0, 100);
    hSelection2->SetDirectory(0);

    hSelectionElectrons2 = new TH2D("hSelectionElectrons3", "dE/dx vs P (electrons) step 2; P [GeV/c]; dE/dx",
				   100, 0.2, 0.8, 
				   100, 0, 100);
   hSelectionElectrons2->SetDirectory(0);

    hDeDxVsEta = new TH2D("hDeDxVsEta", "dE/dx vs #eta (pion MIPs); #eta; dE/dx",
			  50, -1.0, 1.0, 
			  fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDxVsEta->Sumw2();
    hDeDxVsEta->SetDirectory(0);
    
    hDeDxVsEtaElectrons = new TH2D("hDeDxVsEtaElectrons", "dE/dx vs #eta (electrons); #eta; dE/dx",
				 50, -1.0, 1.0, 
				   95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxVsEtaElectrons->Sumw2();
    hDeDxVsEtaElectrons->SetDirectory(0);

    hNclVsEta = new TH2D("hNclVsEta", "Ncl vs #eta (pion MIPs); #eta; Ncl",
			 50, -1.0, 1.0, 
			 18, 69.5, 159.5); 
    hNclVsEta->Sumw2();
    hNclVsEta->SetDirectory(0);
    hNclVsEtaElectrons = new TH2D("hNclVsEtaElectrons", "Ncl vs #eta (electrons); #eta; Ncl",
				  50, -1.0, 1.0, 
				  18, 69.5, 159.5); 
    hNclVsEtaElectrons->Sumw2();
    hNclVsEtaElectrons->SetDirectory(0);
    
    fInit = 2;

    break;

  case 3:

  //
  // step 3 histograms 
  //
    hSelection3 = new TH2D("hSelection3", "dE/dx vs P (pion MIPs) step 3; P [GeV/c]; dE/dx",
			   100, 0.2, 0.8, 
			   100, 0, 100);
    hSelection3->SetDirectory(0);

    hSelectionElectrons3 = new TH2D("hSelectionElectrons3", "dE/dx vs P (electrons) step 3; P [GeV/c]; dE/dx",
				   100, 0.2, 0.8, 
				   100, 0, 100);
    hSelectionElectrons3->SetDirectory(0);
    

    hMeanP = new TProfile("hMeanP", "mean p; p [GeV/c]; mean p",
			  nPtBins, ptBins);
    hMeanP->SetDirectory(0);
    
    hDeDxVsP = new TH2D("hDeDxVsP", "dE/dx vs P; p [GeV/c]; dE/dx",
			nPtBins, ptBins, 55, 40, 95);
    hDeDxVsP->Sumw2();
    hDeDxVsP->SetDirectory(0);
    // dE/dx vs p 
    if(fIsMc) {
      
      hDeDxVsPPiMc = new TH2D("hDeDxVsPPiMc", "dE/dx vs P; p [GeV/c]; dE/dx",
			      nPtBins, ptBins, 55, 40, 95);
      hDeDxVsPPiMc->Sumw2();
      hDeDxVsPPiMc->SetDirectory(0);
      
      hDeDxVsPKMc = new TH2D("hDeDxVsPKMc", "dE/dx vs P; p [GeV/c]; dE/dx",
			     nPtBins, ptBins, 55, 40, 95);
      hDeDxVsPKMc->Sumw2();
      hDeDxVsPKMc->SetDirectory(0);
      
      hDeDxVsPPMc = new TH2D("hDeDxVsPPMc", "dE/dx vs P; p [GeV/c]; dE/dx",
			     nPtBins, ptBins, 55, 40, 95);
      hDeDxVsPPMc->Sumw2();
      hDeDxVsPPMc->SetDirectory(0);
      
      hDeDxVsPEMc = new TH2D("hDeDxVsPEMc", "dE/dx vs P; p [GeV/c]; dE/dx",
			     nPtBins, ptBins, 55, 40, 95);
      hDeDxVsPEMc->Sumw2();
      hDeDxVsPEMc->SetDirectory(0);
      
      hDeDxVsPMuMc = new TH2D("hDeDxVsPMuMc", "dE/dx vs P; p [GeV/c]; dE/dx",
			      nPtBins, ptBins, 55, 40, 95);
      hDeDxVsPMuMc->Sumw2();
      hDeDxVsPMuMc->SetDirectory(0);
    }
    
    //
    // step 3 histograms reated to the eta calibration
    //
    
    hDeDxVsEtaCal = new TH2D("hDeDxVsEtaCal", "dE/dx vs #eta calibrated (pion MIPs); #eta; dE/dx",
			     50, -1.0, 1.0, 
			     fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDxVsEtaCal->Sumw2();
    hDeDxVsEtaCal->SetDirectory(0);
    
    hDeDxVsEtaCalElectrons = new TH2D("hDeDxVsEtaCalElectrons", "dE/dx vs #eta calibrated (electrons); #eta; dE/dx",
				      50, -1.0, 1.0, 
				      95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxVsEtaCalElectrons->Sumw2();
    hDeDxVsEtaCalElectrons->SetDirectory(0);
  
    //
    // step 3 histograms related to resolution
    //
    hMeanEta = new TProfile("hMeanEta", "<|#eta|> in |#eta| intervals (pion MIPs); |#eta|; <|#eta|>",
			    4, 0, 0.8);
    hMeanEta->SetDirectory(0);

    hMeanEtaElectrons = new TProfile("hMeanEtaElectrons", "<|#eta|> in |#eta| intervals (electrons); |#eta|; <|#eta|>",
				     4, 0, 0.8);
    hMeanEtaElectrons->SetDirectory(0);
    
    hDeDx = new TH1D("hDeDx", "dE/dx (pion MIPs); dE/dx; Counts",
		     fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDx->Sumw2();
    hDeDx->SetDirectory(0);
    hDeDx1 = new TH1D("hDeDx1", "dE/dx |#eta|<0.2 (pion MIPs); dE/dx; Counts",
		      fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDx1->Sumw2();
    hDeDx1->SetDirectory(0);
    hDeDx2 = new TH1D("hDeDx2", "dE/dx 0.2<|#eta|<0.4 (pion MIPs); dE/dx; Counts",
		      fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDx2->Sumw2();
    hDeDx2->SetDirectory(0);
    hDeDx3 = new TH1D("hDeDx3", "dE/dx 0.4<|#eta|<0.6 (pion MIPs); dE/dx; Counts",
		      fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDx3->Sumw2();
    hDeDx3->SetDirectory(0);
    hDeDx4 = new TH1D("hDeDx4", "dE/dx 0.6<|#eta|<0.8 (pion MIPs); dE/dx; Counts",
		      fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDx4->Sumw2();
    hDeDx4->SetDirectory(0);
    
    hDeDxElectrons = new TH1D("hDeDxElectrons", "dE/dx (electrons); dE/dx; Counts",
			      95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxElectrons->Sumw2();
    hDeDxElectrons->SetDirectory(0);
    hDeDxElectrons1 = new TH1D("hDeDxElectrons1", "dE/dx |#eta|<0.2 (electrons); dE/dx; Counts",
			       95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxElectrons1->Sumw2();
    hDeDxElectrons1->SetDirectory(0);
    hDeDxElectrons2 = new TH1D("hDeDxElectrons2", "dE/dx 0.2<|#eta|<0.4 (electrons); dE/dx; Counts",
			       95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxElectrons2->Sumw2();
    hDeDxElectrons2->SetDirectory(0);
    hDeDxElectrons3 = new TH1D("hDeDxElectrons3", "dE/dx 0.4<|#eta|<0.6 (electrons); dE/dx; Counts",
			       95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxElectrons3->Sumw2();
    hDeDxElectrons3->SetDirectory(0);
    hDeDxElectrons4 = new TH1D("hDeDxElectrons4", "dE/dx 0.6<|#eta|<0.8 (electrons); dE/dx; Counts",
			       95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxElectrons4->Sumw2();
    hDeDxElectrons4->SetDirectory(0);
    
    //
    // step 3 histograms for resolution
    //
    
    
    
    hDeDxVsNcl = new TH2D("hDeDxVsNcl", "dE/dx vs Ncl (pion MIPs); Ncl; dE/dx",
			  18, 69.5, 159.5, 
			  fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDxVsNcl->Sumw2();
    hDeDxVsNcl->SetDirectory(0);
    hDeDxVsNcl1 = new TH2D("hDeDxVsNcl1", "dE/dx vs Ncl |#eta|<0.2 (pion MIPs); Ncl; dE/dx",
			   18, 69.5, 159.5, 
			   fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDxVsNcl1->Sumw2();
    hDeDxVsNcl1->SetDirectory(0);
    hDeDxVsNcl2 = new TH2D("hDeDxVsNcl2", "dE/dx vs Ncl 0.2<|#eta|<0.4 (pion MIPs); Ncl; dE/dx",
			   18, 69.5, 159.5, 
			   fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDxVsNcl2->Sumw2();
    hDeDxVsNcl2->SetDirectory(0);
    hDeDxVsNcl3 = new TH2D("hDeDxVsNcl3", "dE/dx vs Ncl 0.4<|#eta|<0.6 (pion MIPs); Ncl; dE/dx",
			   18, 69.5, 159.5, 
			   fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDxVsNcl3->Sumw2();
    hDeDxVsNcl3->SetDirectory(0);
    hDeDxVsNcl4 = new TH2D("hDeDxVsNcl4", "dE/dx vs Ncl 0.6<|#eta|<0.8 (pion MIPs); Ncl; dE/dx",
			   18, 69.5, 159.5, 
			   fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
    hDeDxVsNcl4->Sumw2();
    hDeDxVsNcl4->SetDirectory(0);
    
    hDeDxVsNclElectrons = new TH2D("hDeDxVsNclElectrons", "dE/dx vs Ncl (electrons); Ncl; dE/dx",
				   18, 69.5, 159.5, 
				   95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxVsNclElectrons->Sumw2();
    hDeDxVsNclElectrons->SetDirectory(0);
    hDeDxVsNclElectrons1 = new TH2D("hDeDxVsNclElectrons1", "dE/dx vs Ncl |#eta|<0.2 (electrons); Ncl; dE/dx",
				    18, 69.5, 159.5, 
				    95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxVsNclElectrons1->Sumw2();
    hDeDxVsNclElectrons1->SetDirectory(0);
    hDeDxVsNclElectrons2 = new TH2D("hDeDxVsNclElectrons2", "dE/dx vs Ncl 0.2<|#eta|<0.4 (electrons); Ncl; dE/dx",
			 18, 69.5, 159.5, 
				    95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxVsNclElectrons2->Sumw2();
    hDeDxVsNclElectrons2->SetDirectory(0);
    hDeDxVsNclElectrons3 = new TH2D("hDeDxVsNclElectrons3", "dE/dx vs Ncl 0.4<|#eta|<0.6 (electrons); Ncl; dE/dx",
				    18, 69.5, 159.5, 
			 95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxVsNclElectrons3->Sumw2();
    hDeDxVsNclElectrons3->SetDirectory(0);
    hDeDxVsNclElectrons4 = new TH2D("hDeDxVsNclElectrons4", "dE/dx vs Ncl 0.6<|#eta|<0.8 (electrons); Ncl; dE/dx",
				    18, 69.5, 159.5, 
				    95-fDeDxMIPMax, fDeDxMIPMax, 95);
    hDeDxVsNclElectrons4->Sumw2();
    hDeDxVsNclElectrons4->SetDirectory(0);

    fInit = 3;

    break;
  default:
    cout << "No init implemented for step: " << step << endl;
  }
}

//_________________________________________________________
void AliHighPtDeDxCalib::FillTrackInfo(Float_t weight)
{
  if(fStep > 2) { // calibate eta dependence
    
    if(fTrackEta < 0) 
      fTrackDeDx *= 50.0 / fDeDxVsEtaNeg->Eval(fTrackEta);
    else
      fTrackDeDx *= 50.0 / fDeDxVsEtaPos->Eval(fTrackEta);
  }

  if(fStep > 1) { // calibate ncl dependence
    
    fTrackDeDx *= 50.0 / fDeDxVsNcl->Eval(fTrackNcl);
  }
  
  //
  // Fill information for MIP pions
  //
  if(IsMIP()) {
    
    if(fStep == 2) { // calibate eta
      
      hSelection2->Fill(fTrackP, fTrackDeDx);
      hDeDxVsEta->Fill(fTrackEta, fTrackDeDx);
      hNclVsEta->Fill(fTrackEta, fTrackNcl);
    } else if (fStep==1) { // calibrate nCl dependence
      
      hSelection1->Fill(fTrackP, fTrackDeDx);
      hDeDxVsNclBefore->Fill(fTrackNcl, fTrackDeDx);

    } else if (fStep==3) { // calibrate <dE/dx>pi and check eta calibration
      
      AliHighPtDeDxBase::FillTrackInfo(weight);

      hSelection3->Fill(fTrackP, fTrackDeDx);

      hDeDxVsEtaCal->Fill(fTrackEta, fTrackDeDx);
      
      hMeanEta->Fill(TMath::Abs(fTrackEta), TMath::Abs(fTrackEta));
      hDeDx->Fill(fTrackDeDx);
      hDeDxVsNcl->Fill(fTrackNcl, fTrackDeDx);
      if(TMath::Abs(fTrackEta)<0.2) {
	hDeDx1->Fill(fTrackDeDx);
	hDeDxVsNcl1->Fill(fTrackNcl, fTrackDeDx);
      } else if(TMath::Abs(fTrackEta)>=0.2&&TMath::Abs(fTrackEta)<0.4) {
	hDeDx2->Fill(fTrackDeDx);
	hDeDxVsNcl2->Fill(fTrackNcl, fTrackDeDx);
      } else if(TMath::Abs(fTrackEta)>=0.4&&TMath::Abs(fTrackEta)<0.6) {
	hDeDx3->Fill(fTrackDeDx);
	hDeDxVsNcl3->Fill(fTrackNcl, fTrackDeDx);
      }else if(TMath::Abs(fTrackEta)>=0.6&&TMath::Abs(fTrackEta)<0.8) { 
	hDeDx4->Fill(fTrackDeDx);
	hDeDxVsNcl4->Fill(fTrackNcl, fTrackDeDx);
      }
    }
  }

  //
  // Fill information for electrons with same p as MIP pions
  //
  // This is done to validate the calibration/assumptions about sigma for a
  // group of tracks with same topology but much higher dE/dx (plateau)
  //
  if(fStep > 1 && IsElectron()) {

    if(fStep == 2) { // calibate eta
      
      hSelectionElectrons2->Fill(fTrackP, fTrackDeDx);
      hDeDxVsEtaElectrons->Fill(fTrackEta, fTrackDeDx);
      hNclVsEtaElectrons->Fill(fTrackEta, fTrackNcl);
    } else if (fStep==3) { // calibrate <dE/dx>pi and check eta calibration
      
      hSelectionElectrons3->Fill(fTrackP, fTrackDeDx);

      hDeDxVsEtaCalElectrons->Fill(fTrackEta, fTrackDeDx);
      
      hMeanEtaElectrons->Fill(TMath::Abs(fTrackEta), TMath::Abs(fTrackEta));
      hDeDxElectrons->Fill(fTrackDeDx);
      hDeDxVsNclElectrons->Fill(fTrackNcl, fTrackDeDx);
      if(TMath::Abs(fTrackEta)<0.2) {
	hDeDxElectrons1->Fill(fTrackDeDx);
	hDeDxVsNclElectrons1->Fill(fTrackNcl, fTrackDeDx);
      } else if(TMath::Abs(fTrackEta)>=0.2&&TMath::Abs(fTrackEta)<0.4) {
	hDeDxElectrons2->Fill(fTrackDeDx);
	hDeDxVsNclElectrons2->Fill(fTrackNcl, fTrackDeDx);
      } else if(TMath::Abs(fTrackEta)>=0.4&&TMath::Abs(fTrackEta)<0.6) {
	hDeDxElectrons3->Fill(fTrackDeDx);
	hDeDxVsNclElectrons3->Fill(fTrackNcl, fTrackDeDx);
      }else if(TMath::Abs(fTrackEta)>=0.6&&TMath::Abs(fTrackEta)<0.8) { 
	hDeDxElectrons4->Fill(fTrackDeDx);
	hDeDxVsNclElectrons4->Fill(fTrackNcl, fTrackDeDx);
      }
    }
  }

  //
  // Fill information for high pT tracks (dE/dx vs p)
  //
  if(fStep==3) { // Fill dE/dx vs 
    
    hDeDxVsP->Fill(fTrackP, fTrackDeDx, weight);
    hMeanP->Fill(fTrackP, fTrackP);    

    if(fIsMc) {
      
      switch (fTrackPidMc) {
	
      case 1: // pion
	  hDeDxVsPPiMc->Fill(fTrackP, fTrackDeDx, weight);
	break;
      case 2: // kaon
	  hDeDxVsPKMc ->Fill(fTrackP, fTrackDeDx, weight);
	break;
      case 3: // proton
	  hDeDxVsPPMc ->Fill(fTrackP, fTrackDeDx, weight);
	break;
      case 4: // electron
	  hDeDxVsPEMc ->Fill(fTrackP, fTrackDeDx, weight);
	break;
      case 5: // muon
	  hDeDxVsPMuMc ->Fill(fTrackP, fTrackDeDx, weight);
	break;
      default:
	break;
      }
    }
  }
}
  

//___________________________________________________________________________
Bool_t AliHighPtDeDxCalib::IsMIP()
{
  if(fTrackP > fPMIPMin && fTrackP<fPMIPMax &&
     fTrackDeDx > fDeDxMIPMin && fTrackDeDx < fDeDxMIPMax)
    return kTRUE;
  
  return kFALSE;
}

//___________________________________________________________________________
Bool_t AliHighPtDeDxCalib::IsElectron()
{
  if(fTrackP > fPMIPMin && fTrackP<fPMIPMax &&
     fTrackDeDx > (fDeDxMIPMax+5) && TMath::Abs(fTrackBeta-1)<fDeltaBeta)
    return kTRUE;
  
  return kFALSE;
}


//___________________________________________________________________________
void AliHighPtDeDxCalib::PerformEtaCal()
{
  TProfile* hDeDxVsEtaProf = hDeDxVsEta->ProfileX();
  hDeDxVsEtaProf->SetMarkerStyle(29);
  hDeDxVsEtaProf->Fit(fDeDxVsEtaNeg, "0", "", -1, 0);
  hDeDxVsEtaProf->Fit(fDeDxVsEtaPos, "0", "",  0, 1);
  delete hDeDxVsEtaProf;
}

//___________________________________________________________________________
void AliHighPtDeDxCalib::PerformNclCal()
{
  TProfile* hDeDxVsNclProf = hDeDxVsNclBefore->ProfileX();
  hDeDxVsNclProf->SetMarkerStyle(29);
  hDeDxVsNclProf->Fit(fDeDxVsNcl, "0", "", 69.5, 159.5);
  delete hDeDxVsNclProf;
}

//___________________________________________________________________________
TCanvas* AliHighPtDeDxCalib::DrawNclCal()
{
  TCanvas* cNcl = new TCanvas("cNcl", "dE/dx vs Ncl for MIP pions", 600, 400);
  cNcl->Clear();
  cNcl->cd();
  hDeDxVsNclBefore->DrawCopy("COL");
  TProfile* hDeDxVsNclProf = hDeDxVsNclBefore->ProfileX();
  hDeDxVsNclProf->SetMarkerStyle(29);
  hDeDxVsNclProf->DrawCopy("SAME");
  fDeDxVsNcl->Draw("SAME");
  delete hDeDxVsNclProf;
  return cNcl;
}

//___________________________________________________________________________
TCanvas* AliHighPtDeDxCalib::DrawEta(Bool_t forMIP)
{
  // Draw dE/dx vs eta for step 1 (calibration)
  if(forMIP) {

    TCanvas* cEta = new TCanvas("cEta", "dE/dx vs Eta for MIP pions", 600, 400);
    cEta->Clear();
    cEta->cd();
    hDeDxVsEta->DrawCopy("COL");

    TProfile* hDeDxVsEtaProf = hDeDxVsEta->ProfileX();
    hDeDxVsEtaProf->SetMarkerStyle(29);
    hDeDxVsEtaProf->DrawCopy("SAME");
    fDeDxVsEtaNeg->DrawCopy("SAME");
    fDeDxVsEtaPos->DrawCopy("SAME");

    return cEta;
  }

  TCanvas* cEta = new TCanvas("cEtaElectrons", "dE/dx vs Eta for electrons", 600, 400);
  cEta->Clear();
  cEta->cd();
  hDeDxVsEtaElectrons->DrawCopy("COL");

  TProfile* hDeDxVsEtaProfElectrons = hDeDxVsEtaElectrons->ProfileX();
  hDeDxVsEtaProfElectrons->SetMarkerStyle(29);
  hDeDxVsEtaProfElectrons->DrawCopy("SAME");

  return cEta;
}     

//___________________________________________________________________________
TCanvas* AliHighPtDeDxCalib::DrawEtaCalibrated(Bool_t forMIP)
{
  // Draw dE/dx vs eta for step 2 (after calibration)
  if(forMIP) {

    TCanvas* cEtaCal = new TCanvas("cEtaCal", "dE/dx vs Eta for MIP pions (calibrated)", 600, 400);
    cEtaCal->Clear();
    cEtaCal->cd();
    hDeDxVsEtaCal->DrawCopy("COL");

    TProfile* hDeDxVsEtaCalProf = hDeDxVsEtaCal->ProfileX();
    hDeDxVsEtaCalProf->SetMarkerStyle(29);
    hDeDxVsEtaCalProf->DrawCopy("SAME");
    
    return cEtaCal;
  }

  TCanvas* cEtaCalElectrons = new TCanvas("cEtaCalElectrons", "dE/dx vs Eta (calibrated electrons)", 600, 400);
  cEtaCalElectrons->Clear();
  cEtaCalElectrons->cd();
  hDeDxVsEtaCalElectrons->DrawCopy("COL");
  TProfile* hDeDxVsEtaCalElectronsProf = 
    hDeDxVsEtaCalElectrons->ProfileX();
  hDeDxVsEtaCalElectronsProf->SetMarkerStyle(29);
  hDeDxVsEtaCalElectronsProf->DrawCopy("SAME");
  
  return cEtaCalElectrons;
}     

//___________________________________________________________________________
TCanvas* AliHighPtDeDxCalib::DrawSelectionHistograms(Int_t step)
{
  TCanvas* cSelection = new TCanvas("cSelection", "dE/dx vs p selection", 800, 800);
  cSelection->Clear();
  cSelection->Divide(2,2);

  cSelection->cd(1);
  if(step==1)
    hSelection1->DrawCopy("COL");
  else if(step==2)
    hSelection2->DrawCopy("COL");
  else if(step==3)
    hSelection3->DrawCopy("COL");

  cSelection->cd(2);
  if(step==2)
    hSelectionElectrons2->DrawCopy("COL");
  else if(step==3)
    hSelectionElectrons3->DrawCopy("COL");
  
  cSelection->cd(3);
  hNclVsEta->DrawCopy("COL");

  cSelection->cd(4);
  hNclVsEtaElectrons->DrawCopy("COL");
  
  return cSelection;
}

//___________________________________________________________________________
TH1D* AliHighPtDeDxCalib::GetHistDeDx(Bool_t forMIP, Int_t etaBin)
{
  switch (etaBin) {
	
  case 0:
    if(forMIP)
      return hDeDx;
    else
      return hDeDxElectrons;
    break;
  case 1:
    if(forMIP)
      return hDeDx1;
    else
      return hDeDxElectrons1;
    break;
  case 2:
    if(forMIP)
      return hDeDx2;
    else
      return hDeDxElectrons2;
    break;
  case 3:
    if(forMIP)
      return hDeDx3;
    else
      return hDeDxElectrons3;
    break;
  case 4:
    if(forMIP)
      return hDeDx4;
    else
      return hDeDxElectrons4;
    break;
  default:
    cout << "Eta bin: " << etaBin << " not found" << endl;
    break;
  }
  return 0;
}

//___________________________________________________________________________
TH2D* AliHighPtDeDxCalib::GetHistDeDxVsNcl(Bool_t forMIP, Int_t etaBin)
{
  switch (etaBin) {
	
  case 0:
    if(forMIP)
      return hDeDxVsNcl;
    else
      return hDeDxVsNclElectrons;
    break;
  case 1:
    if(forMIP)
      return hDeDxVsNcl1;
    else
      return hDeDxVsNclElectrons1;
    break;
  case 2:
    if(forMIP)
      return hDeDxVsNcl2;
    else
      return hDeDxVsNclElectrons2;
    break;
  case 3:
    if(forMIP)
      return hDeDxVsNcl3;
    else
      return hDeDxVsNclElectrons3;
    break;
  case 4:
    if(forMIP)
      return hDeDxVsNcl4;
    else
      return hDeDxVsNclElectrons4;
    break;
  default:
    cout << "Eta bin: " << etaBin << " not found" << endl;
    break;
  }
  return 0;
}

//___________________________________________________________________________
TH2D* AliHighPtDeDxCalib::GetHistDeDxVsP(Int_t pid)
{
  switch (pid) {
	
  case 0:
    return hDeDxVsP;
    break;
  case 1:
    return hDeDxVsPPiMc;
    break;
  case 2:
    return hDeDxVsPKMc;
    break;
  case 3:
    return hDeDxVsPPMc;
    break;
  case 4:
    return hDeDxVsPEMc;
    break;
  case 5:
    return hDeDxVsPMuMc;
    break;
  default:
    cout << "PID: " << pid << " not found" << endl;
    break;
  }
  return 0;
}

/*

      }

  }

*/
