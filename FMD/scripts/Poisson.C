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
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
#include <iostream>

/** @class Poisson
    @brief Make a poisson reconstruction
    @code 
    Root> .L Compile.C
    Root> Compile("Poisson.C")
    Root> Poisson c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class Poisson : public AliFMDInput
{
protected:
  TH2D*  fEmpty; // Histogram 
  TH2D*  fTotal; // Histogram 
  TH2D*  fMult;  // Histogram 
  TFile* fFile;  // File 
  Int_t  fEv;    // Event number
  Double_t fThreshold;
public:
  /** Constructor 
      @param threshold Threshold
      @param nEta      # of @f$ \eta@f$ bins
      @param minEta    minimum @f$ \eta@f$
      @param maxEta    maximum @f$ \eta@f$
      @param nPhi      # of @f$ \eta@f$ bins
      @param minPhi    minimum @f$ \varphi@f$  
      @param maxPhi    maximum @f$ \varphi@f$ */
  Poisson(Double_t threshold=.3,
	  Int_t nEta=120, Float_t minEta=-6, Float_t maxEta=6, 
	  Int_t nPhi=4,   Float_t minPhi=0,  Float_t maxPhi=2*TMath::Pi())
    : fFile(0), fEv(0), fThreshold(threshold)
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
  /** Initialize the analyser. Opens the output file. 
      @return @c true on success. */
  virtual Bool_t Init() 
  {
    if (!AliFMDInput::Init()) return kFALSE;
    fFile = TFile::Open("poisson.root", "RECREATE");
    if (!fFile) return kFALSE;
    return kTRUE;
  }
  /** Begining of event
      @param event Event number
      @return @c false on error */
  virtual Bool_t Begin(Int_t event) 
  {
    if (!AliFMDInput::Begin(event)) return kFALSE;
    fEv = event;
    fEmpty->Clear();
    fTotal->Clear();
    fMult->Clear();
    return kTRUE;
  }
  /** Process ESD data.  For each strip, check if the
      psuedo-multiplicity is less than the threshold.  If it is, then
      count the strip as empty. 
      @param esd ESD data 
      @return @c true on success. */
  virtual Bool_t ProcessESD(AliESDFMD* esd)
  {
    for (UShort_t det = 1; det <= 3; det++) {
      for (UShort_t rng = 0; rng < 2; rng++) {
	Char_t ring = (rng == 0 ? 'I' : 'O');
	// Not covered channels 
	for (UShort_t sec = 0; sec < 40; sec++) {
	  for (UShort_t str = 0; str < 512; str++) {
	    Float_t mult = esd->Multiplicity(det, ring, sec, str);
	    Float_t eta  = esd->Eta(det, ring, sec, str);
	    // Dead channels, or not covered. 
	    if (mult >= AliESDFMD::kInvalidMult) continue;
	    if (eta  >= AliESDFMD::kInvalidEta) continue;
	    Float_t phi;
	    switch (ring) { 
	    case 'I':  phi = (sec + .5) * 2 * TMath::Pi() / 20; break;
	    case 'O':  phi = (sec + .5) * 2 * TMath::Pi() / 40; break;
	    }
	    fTotal->Fill(eta, phi);
	    if (mult < fThreshold) fEmpty->Fill(eta, phi);
	  } // Loop over strips
	} // Loop over sectors
      } // Loop over rings
    } // Loop over detectors
    return kTRUE;
  }
  /** For each bin, reconstruct the charge particle multiplicity as 
      @f[
      m = - N_{total} \log\left(\frac{N_{empty}}{N_{total}}\right)
      @f]
      where @f$ N_{total}@f$ is the total number of strips in the bin,
      and @f$ N_{empty}@f$ is the number of strips in the bin that did
      not fire. 
      @return @c true  */
  virtual Bool_t End() 
  {
    for (Int_t etaBin = 1; etaBin <= fEmpty->GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= fEmpty->GetNbinsY(); phiBin++) {
	Double_t empty  = fEmpty->GetBinContent(etaBin, phiBin);
	Double_t total  = fTotal->GetBinContent(etaBin, phiBin);
	Double_t lambda = (empty > 0 ? - TMath::Log(empty / total) : 1);
	Double_t mult   = lambda * total;
	fMult->SetBinContent(etaBin, phiBin, mult);
      }
    }
    fFile->cd();
    fMult->Write(Form("mult%03d", fEv));
    if (!gROOT->IsBatch()) { 
      gStyle->SetPalette(1);
      TCanvas* c = new TCanvas("poisson", "Poisson multiplicity");
      c->SetFillColor(0);
      fMult->Draw("colz");
    }
    return AliFMDInput::End();
  }
  
  /** At end of run.  Write and close output file @c poisson.root 
      @return @c true on success */
  virtual Bool_t Finish()
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
