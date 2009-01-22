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
#include <TH1D.h>
#include <AliFMDHit.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>
#include <AliFMDUShortMap.h>
#include <AliFMDFloatMap.h>
#include <AliFMDRecPoint.h>
#include <AliESDFMD.h>
#include <AliLog.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TLegend.h>
#include <TLine.h>

/** @class DrawESD
    @brief Draw digit ADC versus Rec point mult
    @code 
    Root> .L Compile.C
    Root> Compile("DrawESD.C")
    Root> DrawESD c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawESD : public AliFMDInput
{
private:
  TH1D* fMult; // Histogram 
  const Double_t fCorr;
public:
  //__________________________________________________________________
  DrawESD(Int_t n=1000, Double_t mmin=-0.5, Double_t mmax=20.5) 
    : fCorr(1) // 0.68377 / 1.1)
  { 
    AddLoad(kESD);
    fMult = new TH1D("mult", " Multiplicity (strip)", n, mmin, mmax);
    fMult->Sumw2();
    fMult->SetXTitle("Strip Multiplicity");
  }
  //__________________________________________________________________
  /** Begining of event
      @param ev Event number
      @return @c false on error */
  Bool_t Begin(Int_t ev) 
  {
    return AliFMDInput::Begin(ev);
  }
  //__________________________________________________________________
  Bool_t ProcessESD(UShort_t det, 
		    Char_t   rng, 
		    UShort_t sec, 
		    UShort_t str, 
		    Float_t  /* eta */, 
		    Float_t  mult)
  {
    // Cache the energy loss 
    if (mult > 0) 
      Info("ProcessESD", "FMD%d%c[%2d,%3d]=%f", det, rng, sec, str, mult);
    if (mult/fCorr > 0.001) fMult->Fill(mult/fCorr);
    return kTRUE;
  }
  //__________________________________________________________________
  TF1* FitPeak(Int_t n, TH1D* hist, Double_t min, Double_t& max)
  {
    if (TMath::Abs(max-min) < .25) return 0;
    std::cout << "Fit peack in range " << min << " to " << max << std::endl;
    TF1* l = new TF1(Form("l%02d", n), "landau", min, max);
    hist->GetXaxis()->SetRangeUser(0, 4);
    hist->Fit(l, "0Q", "", min, max);
    Double_t mpv   = l->GetParameter(1);
    Double_t empv  = l->GetParError(1);
    Double_t sigma = l->GetParameter(2);
    l->SetRange(mpv-empv, mpv+3*sigma);
    hist->Fit(l, "EMQ0", "", mpv-3*empv, mpv+3*sigma);
    std::cout << "Peak # " << n << " [" << min << "," << max << "]\n"
	      << " MPV: " << l->GetParameter(1) 
	      << " +/- "  << l->GetParError(1) 
	      << " Var: " << l->GetParameter(2) 
	      << " +/- "  << l->GetParError(2)
	      << " Chi^2/NDF: " << l->GetChisquare() / l->GetNDF() 
	      << std::endl;
    mpv   = l->GetParameter(1);
    sigma = l->GetParameter(2);
    min   = mpv - sigma * 2; // (n==1 ? 3 : 2);
    max   = mpv + sigma * 3;
    // l->SetRange(min, max);
    l->Draw("same");
    return  l;
  }
  //__________________________________________________________________
  void MaxInRange(TH1D* hist, Double_t min, Double_t& mean, Double_t& var)
  {
    hist->GetXaxis()->SetRangeUser(min, 4);
    mean = hist->GetMean();
    var  = hist->GetRMS();
  }

  //__________________________________________________________________
  Bool_t Finish()
  {
    Info("Finish", "Will draw results");
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1111111);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    
    if (fMult->GetEntries() <= 0) return kFALSE;
    
    TCanvas* c = new TCanvas("c", "C");
    c->cd();
    // c->SetLogy();
    fMult->GetXaxis()->SetRangeUser(0,4);
    fMult->Scale(1. / fMult->GetEntries());
    fMult->SetStats(kFALSE);
    fMult->SetFillColor(2);
    fMult->SetFillStyle(3001);
    fMult->Draw("hist e");
    
    return kTRUE;
    
    Double_t mean, rms;
    MaxInRange(fMult, 0.2, mean, rms);
    Double_t x1   = mean-rms/2; // .75; // .8;  // .65 / fCorr;
    Double_t x2   = mean+rms/2; // 1.3; // 1.7; // fCorr;
    TF1*     l1   = FitPeak(1, fMult, x1, x2);
    x2            = TMath::Max(mean+rms/2, x2);
    MaxInRange(fMult, x2, mean, rms);
    Double_t x3   = mean + rms;
    TF1*     l2   = FitPeak(2, fMult, x2, x3);
    TF1*     f    = 0;
    Double_t diff = 0;
    if (l2) {
      diff          = l2->GetParameter(1)-l1->GetParameter(1);
      f             = new TF1("user", "landau(0)+landau(3)", x1, x3);
      f->SetParNames("A_{1}", "Mpv_{1}", "#sigma_{1}",
		     "A_{2}", "Mpv_{2}", "#sigma_{2}",
		     "A_{3}", "Mpv_{3}", "#sigma_{3}");
      f->SetParameters(l1->GetParameter(0), 
		       l1->GetParameter(1), 
		       l1->GetParameter(2), 
		       l2 ? l2->GetParameter(0) : 0, 
		       l2 ? l2->GetParameter(1) : 0, 
		       l2 ? l2->GetParameter(2) : 0,
		       l2->GetParameter(0)/10, 
		       l2->GetParameter(1) + diff, 
		       l2->GetParameter(2));
    }
    else { 
      x3 = x2;
      f  = new TF1("usr", "landau", x1, x3);
    }
    
    std::cout << "Range: " << x1 << "-" << x3 << std::endl;
    
    fMult->GetXaxis()->SetRangeUser(0, 4);
    fMult->Fit(f, "0QR", "", x1, x3);
    fMult->Fit(f, "ME0R", "E1", x1, x3);
    fMult->DrawClone("same hist");

    l1->SetLineColor(3);
    l1->SetLineWidth(2);
    l1->SetRange(0, 4);
    l1->Draw("same");
    if (l2) {
      l2->SetLineColor(4);
      l2->SetLineWidth(2);
      l2->SetRange(0, 4);
      l2->Draw("same");
    }
    f->SetLineWidth(2);
    f->SetRange(0, 4);
    f->Draw("same");

    TLegend* l = new TLegend(0.6, 0.6, .89, .89);
    l->AddEntry(l1, "1 particle Landau", "l");
    if (l2) l->AddEntry(l2, "2 particle Landau", "l");
    l->AddEntry(f,  "1+2 particle Landau", "l");
    l->SetFillColor(0);
    l->Draw("same");


#if 0
    c = new TCanvas("c2", "Landaus");
    // c->SetLogy();
    fMult->DrawClone("axis");
    f->Draw("same");
    Double_t* p1 = f->GetParameters();
    Double_t* p2 = &(p1[3]);
    TF1* ll1 = new TF1("ll1", "landau", 0, 4);
    ll1->SetParameters(p1);
    ll1->SetLineColor(3);
    ll1->Draw("same");
    TF1* ll2 = new TF1("ll2", "landau", 0, 4);
    ll2->SetParameters(p2);
    ll2->SetLineColor(4);
    ll2->Draw("same");

    Double_t y1  = fMult->GetMinimum() * 1.1;
    Double_t y2  = fMult->GetMaximum() * .9;
    Double_t xc1 = p1[1]-3*p1[2];
    Double_t xc2 = p2[1]-2*p2[2];
    Double_t xc3 = p2[1]-2*p2[2]+diff;
    TLine* c1 = new TLine(xc1, y1, xc1, y2);
    c1->Draw("same");
    TLine* c2 = new TLine(xc2, y1, xc2, y2);
    c2->Draw("same");
    TLine* c3 = new TLine(xc3, y1, xc3, y2);
    c3->Draw("same");

    l = new TLegend(0.6, 0.6, .89, .89);
    l->AddEntry(ll1, "1 particle Landau", "l");
    l->AddEntry(ll2, "2 particle Landau", "l");
    l->AddEntry(f,  "1+2 particle Landau", "l");
    l->SetFillColor(0);
    l->Draw("same");
#endif


    return kTRUE;
  }

  ClassDef(DrawESD,0);
  
};

//____________________________________________________________________
//
// EOF
//
