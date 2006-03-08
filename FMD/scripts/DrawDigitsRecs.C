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
#include <AliFMDUShortMap.h>
#include <AliFMDFloatMap.h>
#include <AliFMDMultStrip.h>
#include <AliFMDMultRegion.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <TCanvas.h>

class DrawDigitsRecs : public AliFMDInputDigits
{
private:
  TH2D* fAdcVsSingle; // Histogram 
  TH2D* fAdcVsRegion; // Histogram 
  TH2D* fSingleVsRegion; // Histogram 
  AliFMDUShortMap fMap;
  AliFMDFloatMap fEta;
  AliFMDFloatMap fPhi;
  AliFMDFloatMap fMult;
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
  DrawDigitsRecs(Int_t m=1100, Double_t amin=-0.5, Double_t amax=1099.5,
		 Int_t n=105, Double_t mmin=-0.5, Double_t mmax=20.5) 
  { 
    AddLoad(kRecPoints);
    fAdcVsSingle = new TH2D("adcVsSingle", "ADC vs. Multiplicity (strip)", 
			    m, amin, amax, n, mmin, mmax);
    fAdcVsSingle->SetXTitle("ADC value");
    fAdcVsSingle->SetYTitle("Strip Multiplicity");

    mmin *= 20;
    mmax *= 20;
    amin *= 50;
    amax *= 50;
    fAdcVsRegion = new TH2D("adcVsRegion", "ADC vs Muliplcity (region)", 
			    m, amin, amax, n, mmin, mmax);
    fAdcVsRegion->SetXTitle("ADC value");
    fAdcVsRegion->SetYTitle("Region Multiplicity");

    fSingleVsRegion = new TH2D("singleVsRegion", "Single vs Region", 
			       n, mmin, mmax, n, mmin, mmax);
    fSingleVsRegion->SetXTitle("Strip Multiplicity");
    fSingleVsRegion->SetYTitle("Region Multiplicity");
  }
  Bool_t Begin(Int_t ev) 
  {
    fMap.Reset();
    fEta.Reset();
    fPhi.Reset();
    return AliFMDInputDigits::Begin(ev);
  }
  Bool_t ProcessDigit(AliFMDDigit* digit) 
  {
    // Cache the energy loss 
    if (!digit) return kFALSE;
    UShort_t det = digit->Detector();
    Char_t   rng = digit->Ring();
    UShort_t sec = digit->Sector();
    UShort_t str = digit->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in digit", str));
      return kTRUE;
    }
    fMap(det, rng, sec, str) = digit->Counts();
    return kTRUE;
  }
  Bool_t Event() 
  {
    if (!AliFMDInputDigits::Event()) return kFALSE;
    Int_t nEv = fTreeR->GetEntries();
    for (Int_t i = 0; i < nEv; i++) {
      Int_t recoRead  = fTreeR->GetEntry(i);
      if (recoRead <= 0) continue;
      Int_t nSingle = fArrayN->GetEntries();
      if (nSingle <= 0) continue;
      for (Int_t j = 0; j < nSingle; j++) {
	AliFMDMultStrip* single = 
	  static_cast<AliFMDMultStrip*>(fArrayN->At(j));
	if (!single) continue;
	UShort_t det = single->Detector();
	Char_t   rng = single->Ring();
	UShort_t sec = single->Sector();
	UShort_t str = single->Strip();
	if (str > 511) {
	  AliWarning(Form("Bad strip number %d in single", str));
	  continue;
	}
	fAdcVsSingle->Fill(fMap(det, rng, sec, str), single->Particles());
	fEta(det, rng, sec, str)  = single->Eta();
	fPhi(det, rng, sec, str)  = single->Phi() * 180 / TMath::Pi();
	fMult(det, rng, sec, str) = single->Particles();
      }    

      Int_t nRegion = fArrayP->GetEntries();
      if (nRegion <= 0) continue;
      for (Int_t j = 0; j < nRegion; j++) {
	AliFMDMultRegion* region = 
	  static_cast<AliFMDMultRegion*>(fArrayP->At(j));
	if (!region) continue;
	UShort_t det  = region->Detector();
	Char_t   rng  = region->Ring();
	UInt_t   suma = 0;
	Float_t  summ = 0;
	
	for (size_t sec = 0; sec < AliFMDMap::kMaxSectors; sec++) {
	  Float_t phi = fPhi(det, rng, sec, 0);	    
	  Bool_t ok   = (phi <= region->MaxPhi() && phi >= region->MinPhi());
	  if (!ok) continue;
	  for (size_t str = 0; str < AliFMDMap::kMaxStrips; str++) {
	    Float_t eta  = fEta(det, rng, sec, str);
	    Float_t sign = eta < 0 ? -1. : 1.;
	    ok           = (sign * eta <= sign * region->MaxEta() 
			    && sign * eta >= sign * region->MinEta());
	    if (!ok) continue;
	    suma += fMap(det, rng, sec, str);
	    summ += fMult(det, rng, sec, str);
	  }
	}
	fAdcVsRegion->Fill(suma, region->Particles()); 
	fSingleVsRegion->Fill(summ, region->Particles());
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

    new TCanvas("c1", "ADC vs. Strip Multiplicity");
    fAdcVsSingle->SetStats(kFALSE);
    fAdcVsSingle->Draw("COLZ");

    new TCanvas("c2", "ADC vs. Region Multiplicity");
    fAdcVsRegion->SetStats(kFALSE);
    fAdcVsRegion->Draw("COLZ");

    new TCanvas("c3", "Single vs. Region Multiplicity");
    fSingleVsRegion->SetStats(kFALSE);
    fSingleVsRegion->Draw("COLZ");

    return kTRUE;
  }

  ClassDef(DrawDigitsRecs,0);
  
};

//____________________________________________________________________
//
// EOF
//
