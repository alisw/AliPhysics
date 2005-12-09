//
// $Id$
//
// Script that contains a class to draw eloss from hits, versus ADC
// counts from digits, using the AliFMDInputHits class in the util library. 
//
// It draws the energy loss versus the p/(mq^2).  It can be overlayed
// with the Bethe-Bloc curve to show how the simulation behaves
// relative to the expected. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <TH2D.h>
#include <AliFMDHit.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>
#include <AliFMDEdepMap.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>

class DrawHitsDigits : public AliFMDInputHits
{
private:
  TH2D* fElossVsAdc; // Histogram 
  AliFMDEdepMap fMap;
public:
  TArrayF MakeLogScale(Int_t n, Double_t min, Double_t max) 
  {
    TArrayF bins(n+1);
    Float_t dp   = n / TMath::Log10(max / min);
    Float_t pmin = TMath::Log10(min);
    bins[0]      = min;
    for (Int_t i = 1; i < n+1; i++) {
      Float_t p = pmin + i / dp;
      bins[i]   = TMath::Power(10, p);
    }
    return bins;
  }
  DrawHitsDigits(Int_t n=900, Double_t emin=1e-3, Double_t emax=10, 
		 Int_t m=1100, Double_t amin=-0.5, Double_t amax=1099.5) 
  { 
    AddLoad(kDigits);
    TArrayF eloss(MakeLogScale(n, emin, emax));
    TArrayF adcs(m+1);
    adcs[0] = amin;
    for (Int_t i = 1; i < m+1; i++) adcs[i] = adcs[i-1] + (amax-amin)/m;
    fElossVsAdc = new TH2D("bad", "#Delta E vs. ADC", 
			   eloss.fN-1, eloss.fArray, adcs.fN-1, adcs.fArray);
    fElossVsAdc->SetXTitle("#Delta E/#Delta x [MeV/cm]");
    fElossVsAdc->SetYTitle("ADC value");
  }
  Bool_t Begin(Int_t ev) 
  {
    fMap.Reset();
    return AliFMDInputHits::Begin(ev);
  }
  Bool_t ProcessHit(AliFMDHit* hit, TParticle*) 
  {
    // Cache the energy loss 
    if (!hit) return kFALSE;
    UShort_t det = hit->Detector();
    Char_t   rng = hit->Ring();
    UShort_t sec = hit->Sector();
    UShort_t str = hit->Strip();
    fMap(det, rng, sec, str).fEdep += hit->Edep();
    fMap(det, rng, sec, str).fN++;
    return kTRUE;
  }
  Bool_t Event() 
  {
    if (!AliFMDInputHits::Event()) return kFALSE;
    Int_t nEv = fTreeD->GetEntries();
    for (Int_t i = 0; i < nEv; i++) {
      Int_t digitRead  = fTreeD->GetEntry(i);
      if (digitRead <= 0) continue;
      Int_t nDigit = fArrayD->GetEntries();
      if (nDigit <= 0) continue;
      for (Int_t j = 0; j < nDigit; j++) {
	AliFMDDigit* digit = static_cast<AliFMDDigit*>(fArrayD->At(j));
	if (!digit) continue;
	UShort_t det = digit->Detector();
	Char_t   rng = digit->Ring();
	UShort_t sec = digit->Sector();
	UShort_t str = digit->Strip();
	fElossVsAdc->Fill(fMap(det, rng, sec, str).fEdep, digit->Counts());
      }    
    }
    return kTRUE;
  }
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    fElossVsAdc->SetStats(kFALSE);
    fElossVsAdc->Draw("COLZ");
    return kTRUE;
  }
};

//____________________________________________________________________
//
// EOF
//
