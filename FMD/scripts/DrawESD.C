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
  DrawESD(Int_t n=420, Double_t mmin=-0.5, Double_t mmax=20.5) 
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
  Bool_t ProcessESD(UShort_t /* det */, 
		    Char_t   /* ring */, 
		    UShort_t /* sec */, 
		    UShort_t /* strip */, 
		    Float_t  /* eta */, 
		    Float_t  mult)
  {
    // Cache the energy loss 
    fMult->Fill(mult/fCorr);
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1111111);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    
    TCanvas* c = new TCanvas("c", "C");
    c->cd();
    c->SetLogy();
    // fMult->GetXaxis()->SetRangeUser(0.2,4);
    // fMult->Sumw2();
    // fMult->Scale(1. / fMult->GetMaximum());
    fMult->Scale(1. / fMult->GetEntries());
    fMult->SetStats(kFALSE);
    fMult->SetFillColor(2);
    fMult->SetFillStyle(3001);
    fMult->Draw();
    Double_t x1 = .75; // .8;  // .65 / fCorr;
    Double_t x2 = 1.3; // 1.7; // fCorr;
    Double_t x3 = 2.5; // 2.7; // fCorr;
    Double_t x4 = 3.7; // fCorr;
    TF1* l1 = new TF1("landau1", "landau", x1, x2);
    TF1* l2 = new TF1("landau2", "landau", x2, x3);
    // TF1* l3 = new TF1("landau3", "landau", x3, x4);
    TF1* f  = new TF1("user", "landau(0)+landau(3)", x1, x3);

    fMult->SetStats(kTRUE);
    fMult->GetXaxis()->SetRangeUser(0, 4);
    fMult->Fit(l1, "E0L", "", x1, x2);
    fMult->Fit(l2, "E0L", "", x2, x3);
    // fMult->Fit(l3, "E0L", "", x3, x4);
    f->SetParNames("A_{1}", "Mpv_{1}", "#sigma_{1}",
		   "A_{2}", "Mpv_{2}", "#sigma_{2}",
		   "A_{3}", "Mpv_{3}", "#sigma_{3}");
    f->SetParameters(l1->GetParameter(0), 
		     l1->GetParameter(1), 
		     l1->GetParameter(2), 
		     l2->GetParameter(0), 
		     l2->GetParameter(1), 
		     l2->GetParameter(2),
		     l1->GetParameter(0), 
		     l1->GetParameter(1) * 3, 
		     l1->GetParameter(2) * 3);
    fMult->Fit(f, "E0L", "", x1, x3);
    fMult->Fit(f, "MEL", "E1", x1, x3);

    TH1* h = static_cast<TH1*>(fMult->DrawClone("same hist"));

    // l1->SetLineStyle(2);
    l1->SetLineColor(3);
    l1->SetLineWidth(2);
    l1->SetRange(0, 4);
    l1->Draw("same");
    // l2->SetLineStyle(3);
    l2->SetLineColor(4);
    l2->SetLineWidth(2);
    l2->SetRange(0, 4);
    l2->Draw("same");
    // l3->SetLineStyle(3);
    // l3->SetLineWidth(2);
    // l3->SetRange(0, 5);
    // l3->Draw("same");
    f->SetLineWidth(2);
    f->SetRange(0, 4);
    f->Draw("same");

    TLegend* l = new TLegend(0.2, 0.6, .6, .89);
    l->AddEntry(l1, "1 particle Landau", "l");
    l->AddEntry(l2, "2 particle Landau", "l");
    l->AddEntry(f,  "1+2 particle Landau", "l");
    l->SetFillColor(0);
    l->Draw("same");


    c = new TCanvas("c2", "Landaus");
    c->SetLogy();
    fMult->DrawClone("axis");
    f->Draw("same");
    TF1* ll1 = new TF1("ll1", "landau", 0, 4);
    ll1->SetParameters(f->GetParameter(0), 
		       f->GetParameter(1), 
		       f->GetParameter(2));
    ll1->SetLineColor(3);
    ll1->Draw("same");
    TF1* ll2 = new TF1("ll2", "landau", 0, 4);
    ll2->SetParameters(f->GetParameter(3), 
		       f->GetParameter(4), 
		       f->GetParameter(5));
    ll2->SetLineColor(4);
    ll2->Draw("same");

    Double_t y1  = fMult->GetMinimum() * 1.1;
    Double_t y2  = fMult->GetMaximum() * .9;
    Double_t xc1 = f->GetParameter(1)-3*f->GetParameter(2);
    Double_t xc2 = f->GetParameter(4)-2*f->GetParameter(5);
    TLine* c1 = new TLine(xc1, y1, xc1, y2);
    c1->Draw("same");
    TLine* c2 = new TLine(xc2, y1, xc2, y2);
    c2->Draw("same");

    l = new TLegend(0.2, 0.6, .6, .89);
    l->AddEntry(ll1, "1 particle Landau", "l");
    l->AddEntry(ll2, "2 particle Landau", "l");
    l->AddEntry(f,  "1+2 particle Landau", "l");
    l->SetFillColor(0);
    l->Draw("same");

#if 0
    c = new TCanvas("s", "Spectrum");
    TSpectrum* s = new TSpectrum(16);
    h->GetXaxis()->SetRangeUser(0.3, 20);
    s->Search(h, .15, "", 0.01);
    c->Update();
    TH1* b = s->Background(h);
    b->SetFillColor(4);
    b->SetFillStyle(3001);
    b->Draw("same");
#endif

    return kTRUE;
  }

  ClassDef(DrawESD,0);
  
};

//____________________________________________________________________
//
// EOF
//
