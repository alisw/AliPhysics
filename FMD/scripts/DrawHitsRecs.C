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
#include <AliFMDDigit.h>
#include <AliFMDRecPoint.h>
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
#include <TF1.h>

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
class DrawHitsRecs : public AliFMDInput
{
private:
  TH2D* fHitEvsAdc;
  TH2D* fHitEvsRecM;  // Histogram 
  TH2D* fHitEvsRecE;  // Histogram 
  TH1D* fDiffE;       // Histogram 
  TH2D* fHitsVsRecM;  // Histogram 
  TH2D* fDiffM;       // Histogram 
  TH1*  fHitEloss;
  TH1*  fRecEloss;
  AliFMDEdepMap  fMap;
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
    AddLoad(kHits);
    AddLoad(kDigits);
    if (fPrimary) AddLoad(kKinematics);
    TArrayF eloss(MakeLogScale(n, emin, emax));
    TArrayF mults(m+1);
    mults[0] = mmin;
    for (Int_t i = 1; i < m+1; i++) mults[i] = mults[i-1] + (mmax-mmin)/m;

    fHitEvsAdc  = new TH2D("hitEvsAdc", "#Delta E_{sim} vs. ADC",
			   n, emin, emax, 1025, -.5, 1024.5);
    fHitEvsAdc->SetXTitle("#Delta E_{sim} [MeV]");
    fHitEvsAdc->SetYTitle("ADC");
    
    fHitEvsRecM = new TH2D("hitEvsRecM", "#Delta E_{sim} vs. M_{rec}", 
			   eloss.fN-1, eloss.fArray, mults.fN-1, mults.fArray);
    fHitEvsRecM->SetXTitle("#Delta E_{sim} [MeV]");
    fHitEvsRecM->SetYTitle("M_{rec}");

    fHitEvsRecE = new TH2D("hitEvsRecE", "#Delta E_{sim} vs. #Delta E_{rec}", 
			    n, emin, emax, n, emin, emax);
    fHitEvsRecE->SetXTitle("#Delta E_{sim} [MeV]");
    fHitEvsRecE->SetYTitle("#Delta E_{rec} [MeV]");


    fDiffE = new TH1D("diffE", 
		      "#frac{#Delta E_{sim}-#Delta E_{rec}}{#Delta E_{sim}}", 
		      1100, -1, 1.1);
    fDiffE->SetXTitle("#frac{#Delta E_{sim}-#Delta E_{rec}}{#Delta E_{sim}}");
    
    Double_t omin = -.5;
    Double_t omax = 7.5;
    Int_t    o    = 8;
    fHitsVsRecM = new TH2D("hitsVsStrip", "# of Hits vs. M_{rec}",
			   o, omin, omax, m, mmin, mmax);
    fHitsVsRecM->SetXTitle("# of Hits");
    fHitsVsRecM->SetYTitle("M_{rec}");

    fDiffM = new TH2D("diffM", "M_{sim} - M_{rec}", 
		      41, -20.5, 20.5, 70, 1.5, 5);
    // 36, -TMath::Pi(),TMath::Pi());
    fDiffM->SetXTitle("M_{sim} - M_{rec}");
    fDiffM->SetYTitle("|#eta|");
    // fDiffM->SetYTitle("Detector");

    fHitEloss = new TH1D("hitEloss", "#frac{#Delta E_{sim}}{#Delta x} (MeV/cm)", 
			 100, 0, 10);
    fHitEloss->SetFillColor(2);
    fHitEloss->SetFillStyle(3001);
    
    fRecEloss = new TH1D("recEloss", "#frac{#Delta E_{rec}}{#Delta x} (MeV/cm)", 
			 100, 0, 10);
    fRecEloss->SetFillColor(4);
    fRecEloss->SetFillStyle(3001);
  }
  //__________________________________________________________________
  /** Begining of event
      @param ev Event number
      @return @c false on error */
  Bool_t Begin(Int_t ev) 
  {
    fMap.Reset();
    return AliFMDInput::Begin(ev);
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
    
    if (hit->Edep()/hit->Length() > 0.1) fHitEloss->Fill(hit->Edep() / hit->Length());
    fMap(det, rng, sec, str).fEdep += hit->Edep();
    fMap(det, rng, sec, str).fN++;
    return kTRUE;
  }

  //__________________________________________________________________
  Bool_t ProcessDigit(AliFMDDigit* digit)
  {
    if (!digit) return kTRUE;

    UShort_t det = digit->Detector();
    Char_t   rng = digit->Ring();
    UShort_t sec = digit->Sector();
    UShort_t str = digit->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in digit", str));
      return kTRUE;
    }
    Double_t edep = fMap(det, rng, sec, str).fEdep;
    if (edep > 0) fHitEvsAdc->Fill(edep, digit->Counts());
      
    return kTRUE;
  }

  //__________________________________________________________________
  Bool_t ProcessRecPoint(AliFMDRecPoint* single) 
  {
    if (!single) return kTRUE;
    UShort_t det = single->Detector();
    Char_t   rng = single->Ring();
    UShort_t sec = single->Sector();
    UShort_t str = single->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in single", str));
      return kTRUE;
    }
    Double_t edep = fMap(det, rng, sec, str).fEdep;
    Int_t    nhit = fMap(det, rng, sec, str).fN;
    if (edep > 0) {
      fHitEvsRecM->Fill(edep, single->Particles());
      fHitEvsRecE->Fill(edep, single->Edep());
      fDiffE->Fill((single->Edep() - edep) / edep);
    }
    if (nhit > 0) fHitsVsRecM->Fill(nhit, single->Particles());
    fDiffM->Fill(nhit - single->Particles(), TMath::Abs(single->Eta()));
    if (single->Edep()/.03 > 0.1) fRecEloss->Fill(single->Edep() / 0.0300);
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
    TCanvas* c = 0;

    c = new TCanvas("c0", fHitEvsAdc->GetTitle());
    fHitEvsAdc->SetStats(kFALSE);
    fHitEvsAdc->Draw("COLZ");

    c = new TCanvas("c1", fHitEvsRecM->GetTitle());
    fHitEvsRecM->SetStats(kFALSE);
    fHitEvsRecM->Draw("COLZ");

    c = new TCanvas("c2", fHitEvsRecE->GetTitle());
    fHitEvsRecE->SetStats(kFALSE);
    fHitEvsRecE->Draw("COLZ");

    c = new TCanvas("c3", fDiffE->GetTitle());
    c->SetLogz();
    fDiffE->Draw();

    c = new TCanvas("c4", fHitsVsRecM->GetTitle());
    c->SetLogz();
    fHitsVsRecM->SetStats(kFALSE);
    fHitsVsRecM->Draw("COLZ");

    c = new TCanvas("c5", fDiffM->GetTitle());
    fDiffM->SetFillColor(2);
    fDiffM->SetFillStyle(3001);
    c->SetLogz();
    fDiffM->Draw("colz");

    c = new TCanvas("c6", "Hit Eloss, Reco Eloss");
    fRecEloss->Scale(1./fRecEloss->GetMaximum());
    fRecEloss->Draw();
    fRecEloss->Fit("landau", "", "SAME", 2, 4);
    TF1* recResp = new TF1(*fRecEloss->GetFunction("landau"));
    fHitEloss->Scale(1./fHitEloss->GetMaximum());
    fHitEloss->Draw("same");
    fHitEloss->Fit("landau", "", "SAME", 2, 10);
    TF1* hitResp = new TF1(*fHitEloss->GetFunction("landau"));
    std::cout << "Hit MPV,width: " << hitResp->GetParameter(1) << ","
	      << hitResp->GetParameter(2) << "\n"
	      << "Rec MPV,width: " << recResp->GetParameter(1) << ","
	      << recResp->GetParameter(2) << std::endl;
    c->SetLogy();

    return kTRUE;
  }

  ClassDef(DrawHitsRecs,0);
};

//____________________________________________________________________
//
// EOF
//
