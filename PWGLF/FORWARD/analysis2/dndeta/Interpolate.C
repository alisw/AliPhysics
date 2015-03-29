#include <THStack.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <TH1.h>
#include <TError.h>
#include <TMultiGraph.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPrincipal.h>

TH1* GetOne(UShort_t sNN, const TString& trigger)
{
  Long_t   p = gROOT->ProcessLine(Form("Drawer::GetStack(0, \"pp\", %d, "
				       "\"%s\", false, true)", 
				       sNN, trigger.Data()));
  THStack* s = (THStack*)p;
  TList*   l = s->GetHists();
  TH1*     h = 0;
  TIter    n(l);
  l->ls();
  while ((h = static_cast<TH1*>(n()))) {
    TString m(h->GetName());
    if (m.EqualTo("dndetaForward_all")) break;
  }

  if (h) {
    switch (sNN) { 
    case  900: h->SetTitle("900GeV");  h->SetMarkerColor(kRed+2);   break;
    case 2760: h->SetTitle("2.76TeV"); h->SetMarkerColor(kGreen+2); break;
    case 7000: h->SetTitle("7TeV");    h->SetMarkerColor(kBlue+2);  break;
    case 8000: h->SetTitle("8TeV");    h->SetMarkerColor(kBlack);   break;
    }
  }

  return h;
}
void OneBin(Int_t bin, 
	    TH1* h0900,
	    TH1* h2760, 
	    TH1* h7000, 
	    TH1* h8000, 
	    TMultiGraph* mg,
	    TNtuple* tuple,
	    Double_t sysErr=0.076)
{
  Info("OneBin", "Getting one bin %d,%p,%p,%p,%p,%p",
       bin, h0900,h2760,h7000,h8000,tuple);
  
  Double_t eta = h0900->GetXaxis()->GetBinCenter(bin);
  Double_t w   = h0900->GetXaxis()->GetBinWidth(bin);
  Info("", "Eta=%f +/- %f", eta, w);
  Double_t e[] = { 900.,  2760., 7000., 8000., 0 };
  TH1*     h[] = { h0900, h2760, h7000, h8000, 0 };
  Float_t  x[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  x[0] = eta;
  x[1] = w;

  TGraphErrors* g = new TGraphErrors(0);
  g->SetName(Form("eta%03d", bin));
  g->SetTitle(Form("%f", eta));
  g->SetMarkerStyle(bin % 10 + 20);
  g->SetMarkerColor(bin % 6 + 2);
  
  Double_t* pe = e;
  TH1**     ph = h;
  Int_t     i  = 0;
  Int_t     j  = 1;
  while (*pe && *ph) { 
    Double_t c = (*ph)->GetBinContent(bin);
    Double_t v = sysErr*c;
    if (c > 1e-6){
      g->SetPoint(i, *pe, c);
      g->SetPointError(i, w, v);
      x[Int_t(2*j+0)] = c;
      x[Int_t(2*j+1)] = v;
      i++;
    }
    j++;
    pe++;
    ph++;
  }
  if (tuple) tuple->Fill(x);
  if (i > 0) mg->Add(g);
  else delete g;
}


void Interpolate(const TString& trigger="INEL")
{
  if (gSystem->Getenv("FWD"))
    fwd = gSystem->Getenv("FWD");
  else 
    fwd = gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
  gROOT->SetMacroPath(Form("%s/dndeta:%s", gROOT->GetMacroPath(),fwd));

  if (!gROOT->GetClass("Drawer"))  gROOT->LoadMacro("Drawer.C+");

  TH1* h0900 = GetOne( 900, trigger);
  TH1* h2760 = GetOne(2760, trigger);
  TH1* h7000 = GetOne(7000, trigger);
  TH1* h8000 = GetOne(8000, trigger);
  Info("","900: %p 2760: %p 7000: %p 8000: %p", 
       h0900, h2760, h7000, h8000);
  Double_t e8000 = (trigger.EqualTo("INEL") ? 0.852 : 0.93);
  h8000->Scale(e8000);

  TFile* out = TFile::Open("trends.root", "RECREATE");
  THStack* sOrig = new THStack("orig", Form("pp - %s", trigger.Data()));
  sOrig->Add(h8000);
  sOrig->Add(h7000);
  sOrig->Add(h2760);
  sOrig->Add(h0900);

  TCanvas* cOrig = new TCanvas("cOrig", "Original", 1200, 1200);
  cOrig->SetTopMargin(0.01);
  cOrig->SetRightMargin(0.01);
  sOrig->Draw("nostack");
  sOrig->GetHistogram()->SetYTitle("1/#it{N} d#it{N}_{ch}/d#it{#eta}");
  sOrig->GetHistogram()->SetXTitle("#it{#eta}");
  sOrig->DrawClone("nostack");
  sOrig->Write();
  
  TLegend* l = cOrig->BuildLegend(.35, .2, .55, .6, "#sqrt{s}");
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);

  cOrig->Modified();
  cOrig->Update();
  cOrig->cd();
  cOrig->Write();
  Info("", "Wrote original"); 

  TCanvas* cG = new TCanvas("cG", "one", 1200, 1200);
  cG->SetTopMargin(0.01);
  cG->SetRightMargin(0.01);
  
  Info("","Creating tuple");
  TNtuple* tuple = new TNtuple("tuple", "Tuple", 
			       "eta:deta:"
			       "v0900:e0900:v2760:e2760:"
			       "v7000:e7000:v8000:e8000");
  TMultiGraph* mg = new TMultiGraph;
  Int_t n = h0900->GetNbinsX();
  Info("","Loop over bins %d", n);

  for (Int_t i = 1; i <= n; i++) {
    Info("", "Getting one bin %d,%p,%p,%p,%p,%p,%p",
	 i, h0900,h2760,h7000,h8000,mg,tuple);
    OneBin(i, h0900, h2760, h7000, h8000, mg, tuple);
  }
  mg->Draw("alp");

  cG->Modified();
  cG->Update();
  cG->cd();

  TPrincipal* p =tuple->Principal("v0900:v2760:v7000:v8000","eta<0", "npdhc");
}
