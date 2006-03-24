//____________________________________________________________________
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
#include <AliFMDMultStrip.h>
#include <AliFMDMultRegion.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>
#include <AliFMDEdepMap.h>
#include <AliFMDFloatMap.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <AliStack.h>
#include <AliLog.h>

//____________________________________________________________________
/** @class DrawHitsRecs
    @brief Draw hit energy loss versus rec point mult
    @code 
    Root> .L Compile.C
    Root> Compile("DrawHitsRecs.C")
    Root> DrawHitsRecs c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawHitsRecs : public AliFMDInputHits
{
private:
  TH2D* fElossVsMult; // Histogram 
  TH2D* fHitsVsStrip;  // Histogram 
  TH2D* fHitsVsRegion;  // Histogram 
  AliFMDEdepMap fMap;
  AliFMDFloatMap fEta;
  AliFMDFloatMap fPhi;
  AliFMDFloatMap fMult;
  Bool_t fPrimary;
public:
  //__________________________________________________________________
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
  //__________________________________________________________________
  DrawHitsRecs(Bool_t primary=kFALSE,
	       Int_t n=900, Double_t emin=1e-3, Double_t emax=10, 
	       Int_t m=21, Double_t mmin=-0.5, Double_t mmax=20.5) 
  { 
    fPrimary = primary;
    AddLoad(kRecPoints);
    if (fPrimary) AddLoad(kKinematics);
    TArrayF eloss(MakeLogScale(n, emin, emax));
    TArrayF mults(m+1);
    mults[0] = mmin;
    for (Int_t i = 1; i < m+1; i++) mults[i] = mults[i-1] + (mmax-mmin)/m;
    fElossVsMult = new TH2D("elossVsMult", 
			    "#Delta E vs. Multiplicity (single)", 
			   eloss.fN-1, eloss.fArray, mults.fN-1, mults.fArray);
    fElossVsMult->SetXTitle("#Delta E/#Delta x [MeV/cm]");
    fElossVsMult->SetYTitle("Strip Multiplicity");

    Double_t omin = -.5;
    Double_t omax = 7.5;
    Int_t    o    = 8;
    fHitsVsStrip = new TH2D("hitsVsStrip", 
			    "# of Hits vs. Multiplicity (strip)",
			   o, omin, omax, m, mmin, mmax);
    fHitsVsStrip->SetXTitle("# of Hits");
    fHitsVsStrip->SetYTitle("Strip Multiplicity");
  }
  //__________________________________________________________________
  /** Begining of event
      @param ev Event number
      @return @c false on error */
  Bool_t Begin(Int_t ev) 
  {
    fMap.Reset();
    return AliFMDInputHits::Begin(ev);
  }
  //__________________________________________________________________
  Bool_t ProcessHit(AliFMDHit* hit, TParticle*) 
  {
    // Cache the energy loss 
    if (!hit) return kFALSE;
    UShort_t det = hit->Detector();
    Char_t   rng = hit->Ring();
    UShort_t sec = hit->Sector();
    UShort_t str = hit->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in hit", str));
      return kTRUE;
    }
    if (fPrimary) {
      TParticle* kine = fStack->Particle(hit->Track());
      if (!kine) return kTRUE;
      if (!kine->IsPrimary()) return kTRUE;
    }
    
    fMap(det, rng, sec, str).fEdep += hit->Edep();
    fMap(det, rng, sec, str).fN++;
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t ProcessRecPoint(AliFMDRecPoint* single) 
  {
    if (!single) continue;
    UShort_t det = single->Detector();
    Char_t   rng = single->Ring();
    UShort_t sec = single->Sector();
    UShort_t str = single->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in single", str));
      continue;
    }
    if (fMap(det, rng, sec, str).fEdep > 0) 
      fElossVsMult->Fill(fMap(det, rng, sec, str).fEdep, 
			 single->Particles());
    if (fMap(det, rng, sec, str).fN > 0) 
      fHitsVsStrip->Fill(fMap(det, rng, sec, str).fN, 
			 single->Particles());
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);

    new TCanvas("c1", "Energy loss vs. Strip Multiplicity");
    fElossVsMult->SetStats(kFALSE);
    fElossVsMult->Draw("COLZ");

    new TCanvas("c2", "# of Hits vs. Strip Multiplicity");
    fHitsVsStrip->SetStats(kFALSE);
    fHitsVsStrip->Draw("COLZ");

    return kTRUE;
  }

  ClassDef(DrawHitsRecs,0);
};

//____________________________________________________________________
//
// EOF
//
