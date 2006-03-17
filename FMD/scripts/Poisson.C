//____________________________________________________________________
//
// $Id$
//
// Script that contains a class to draw hits, using the
// AliFMDInputHits class in the util library. 
//
// It draws the energy loss versus the p/(mq^2).  It can be overlayed
// with the Bethe-Bloc curve to show how the simulation behaves
// relative to the expected. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <AliESDFMD.h>
#include <AliFMDInput.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TArrayF.h>
#include <iostream>

class Poisson : public AliFMDInput
{
private:
  TH2D*  fEmpty; // Histogram 
  TH2D*  fTotal; // Histogram 
  TH2D*  fMult;  // Histogram 
  TFile* fFile;  // File 
  Int_t  fEv;    // Event number
public:
  Poisson(Double_t threshold=.3,
	  Int_t nEta=120, Float_t minEta=-6, Float_t maxEta=6, 
	  Int_t nPhi=4,   Float_t minPhi=0,  Float_t maxPhi=2*TMath::Pi())
    : fFile(0), fEv(0)
  { 
    AddLoad(kESD);

    fEmpty = new TH2D("empty", "# of empty strips", nEta, minEta, maxEta, 
		      nPhi, minPhi, maxPhi);
    fTotal = new TH2D("total", "Total # of strips", nEta, minEta, maxEta, 
		      nPhi, minPhi, maxPhi);
    fMult = new TH2D("mult", "Multiplicity", nEta, minEta, maxEta, 
		      nPhi, minPhi, maxPhi);
    fEmpty->SetXTitle("#eta");  
    fEmpty->SetYTitle("#phi"); 
    fEmpty->SetZTitle("N");
    fTotal->SetXTitle("#eta");  
    fTotal->SetYTitle("#phi"); 
    fTotal->SetZTitle("N");
    fMult->SetXTitle("#eta");  
    fMult->SetYTitle("#phi"); 
    fMult->SetZTitle("<M_{ch}>");

  }
  Bool_t Init() 
  {
    if (!AliFMDInput::Begin(event)) return kFALSE;
    fFile = TFile::Open("poisson.root");
    if (!fFile) return kFALSE;
    return kTRUE;
  }
  Bool_t Begin(Int_t event) 
  {
    if (!AliFMDInput::Begin(event)) return kFALSE;
    fEv = event;
    fEmpty->Clear();
    fTotal->Clear();
    fMult->Clear();
    return kTRUE;
  }
  Bool_t ProcessESD(AliESDFMDHit* esd)
  {
    for (UShort_t det = 1; det <= esd->MaxDetector(); det++) {
      for (UShort_t rng = 0; rng < esd->MaxRing(); rng++) {
	Char_t ring = (rng == 0 ? 'I' : 'O');
	// Not covered channels 
	for (UShort_t sec = 0; sec < esd->MaxSector(); sec++) {
	  for (UShort_t str = 0; str < esd->MaxStrip(); str++) {
	    Float_t mult = esd->Multiplicity(det, ring, sec, str);
	    Float_t eta  = esd->Eta(det, ring, sec, str);
	    // Dead channels, or not covered. 
	    if (mult >= AliESDFMD::kInvalidMult) continue;
	    if (esd  >= AliESDFMD::kInvalidEta) continue;
	    Float_t phi;
	    switch (ring) { 
	    case 'I':  phi = (sec + .5) * 2 * TMath::Pi() / 20; break;
	    case 'O':  phi = (sec + .5) * 2 * TMath::Pi() / 40; break;
	    }
	    fTotal->Fill(eta, phi);
	    if (mult < threshold) fEmpty->Fill(eta, phi);
	  } // Loop over strips
	} // Loop over sectors
      } // Loop over rings
    } // Loop over detectors
  }
  Bool_t End() 
  {
    for (Int_t etaBin = 1; etaBin <= fEmpty->GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= fEmpty->GetNbinsY(); phiBin++) {
	Double_t empty  = fEmpty->GetBinContent(etaBin, phiBin);
	Double_t total  = fTotal->GetBinContent(etaBin, phiBin);
	Double_t lambda = (empty > 0 ? - TMath::Log(empty / nTotal) : 1);
	Double_t mult   = lambda * nTotal;
	fMult->SetBinContent(etaBin, phiBin, mult);
      }
    }
    fFile->cd();
    fMult->Write(Form("mult%03d", fEv));
    return AliFMDInput::End();
  }
  
  Bool_t Finish()
  {
    fFile->Write();
    fFile->Close();
    fFile = 0;
    return AliFMDInput::Finish();
  }
  ClassDef(Poisson,0);
};

//____________________________________________________________________
//
// EOF
//
