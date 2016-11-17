#ifndef __CINT__
#include "AliTrackletAODUtils.C"
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#else
class TCanvas;
class THStack;
class TVirtualPad;
class TH1;
#endif

struct Computer
{
#ifndef __CINT__
  typedef AliTrackletAODUtils U;
#endif
  
  Computer() {}

  void Run(const char* filename, Double_t c1, Double_t c2, Bool_t mid=true)
  {
    Printf("Opening file %s", filename);
    TString       cntN = U::CentName(c1,c2);
    TFile*        file = U::OpenFile(filename);
    U::Container* top  = U::GetC(file, "MidRapidityMCResults");
    U::Container* cent = U::GetC(top,  cntN);

    THStack*      prim = DoOne(cent, mid, "primaries");
    THStack*      seco = DoOne(cent, mid, "secondaries");
    THStack*      fake = DoOne(cent, mid, "combinatorics");
    THStack*      meas = MakeMeas(prim, seco, fake);
    THStack*      norm = Normalize(meas);
    THStack*      rati = ToOther(meas);
    
    prim->SetTitle("P: Primaries");
    seco->SetTitle("S: Secondaries");
    fake->SetTitle("C: Combinatorics");
    meas->SetTitle("M=P+S+C: Measured");
    norm->SetTitle("Normalized M");
    rati->SetTitle("Ratio to total");
    file->Close();
    
    Double_t yr = .93;
    TCanvas* canvas = new TCanvas(filename,filename, 1200, 1200);

    TLegend* l = new TLegend(0.01,yr+.01,.99,.99);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetNColumns(7);
    TIter next(norm->GetHists());
    TH1*  hist = 0;
    while ((hist = static_cast<TH1*>(next()))) {
      TLegendEntry* e = l->AddEntry("dummy",hist->GetTitle(), "p");
      e->SetMarkerStyle(hist->GetMarkerStyle());
      e->SetMarkerColor(hist->GetMarkerColor());
      e->SetMarkerSize (2*hist->GetMarkerSize());
    }
    canvas->cd();
    l->Draw();
    

    TPad* body = new TPad("body","body",0,0,1,yr);
    canvas->cd();
    body->Draw();
    body->SetTopMargin(0.01);
    body->SetRightMargin(0.01);
    body->Divide(2,1);
    
    TVirtualPad* q = body->cd(1);
    q->SetTopMargin(0.01);
    q->SetRightMargin(0.01);
    q->Divide(1,3,0,0);
    
    TVirtualPad* p = 0;
    DrawStack(q->cd(1),prim);
    DrawStack(q->cd(2),seco);
    DrawStack(q->cd(3),fake);

    q = body->cd(2);
    q->SetTopMargin(0.01);
    q->SetRightMargin(0.01);
    q->Divide(1,3,0,0);

    DrawStack(q->cd(1), meas);
    DrawStack(q->cd(2), norm);
    DrawStack(q->cd(3), rati, false, "Ratio");

    canvas->Modified();
    canvas->Update();
    canvas->cd();
    TString outName(filename); outName.ReplaceAll(".root","");
    outName.Prepend("deltaContrib");
    // outName.Append(".png");
    TFile* outRoot = TFile::Open(Form("%s.root", outName.Data()),"RECREATE");
    prim->Write();
    seco->Write();
    fake->Write();
    meas->Write();
    norm->Write();
    rati->Write();
    outRoot->Write();
    canvas->SaveAs(Form("%s.png",outName.Data()));
  }
  void DrawStack(TVirtualPad* p, THStack* hs, Bool_t logy=true,
		 const char* title="d#it{N}_{X}/d#it{#Delta}")
  {
    p->cd();
    p->SetLogy(logy);
    p->SetRightMargin(0.01);
    p->SetTicks();
    p->SetGridx();
    p->SetGridy();
    hs->Draw("nostack");
    hs->GetHistogram()->SetXTitle("#it{#Delta}");
    hs->GetHistogram()->SetYTitle(title);
    hs->GetHistogram()->GetYaxis()->SetTitleSize(0.025/p->GetHNDC());
    hs->GetHistogram()->GetXaxis()->SetTitleSize(0.025/p->GetHNDC());
    hs->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);
    hs->GetHistogram()->GetXaxis()->SetTitleOffset(0.6);
    hs->GetHistogram()->GetYaxis()->SetLabelSize(0.02/p->GetHNDC());
    hs->GetHistogram()->GetXaxis()->SetLabelSize(0.02/p->GetHNDC());
    hs->GetHistogram()->GetXaxis()->SetNdivisions(210);
    hs->GetHistogram()->GetXaxis()->SetLimits(-1,26);
    p->Modified();
  }
  THStack* DoOne(TList* c, Bool_t mid, const char* subName)
  {
    U::Container* sub  = U::GetC(c, subName);
    U::Container* spec = U::GetC(sub, "specie");
    U::Container* regi = U::GetC(spec, mid ? "mid" : "fwd");
    THStack*      all  = U::GetHS(regi, "all");
    if (!all) return 0;
    
    TIter         next(all->GetHists());
    TH1*          hist = 0;
    TH1*          k0s  = 0;
    TH1*          kpm  = 0;
    TH1*          lam  = 0;
    TH1*          xi   = 0;
    TH1*          sim  = 0;
    TH1*          sip  = 0;
    TH1*          tot  = 0;
    TH1*          oth  = 0;
    while ((hist = static_cast<TH1*>(next()))) {
      TString name(hist->GetName());
      if      (name.EqualTo("total")) tot = static_cast<TH1*>(hist->Clone());
      else if (name.EqualTo("321"))   kpm = static_cast<TH1*>(hist->Clone());
      else if (name.EqualTo("310"))   k0s = static_cast<TH1*>(hist->Clone());
      else if (name.EqualTo("3122"))  lam = static_cast<TH1*>(hist->Clone());
      else if (name.EqualTo("3112"))  sim = static_cast<TH1*>(hist->Clone());
      else if (name.EqualTo("3222"))  sip = static_cast<TH1*>(hist->Clone());
      else if (name.EqualTo("3312"))  xi  = static_cast<TH1*>(hist->Clone());
      else {
	if (!oth) oth = static_cast<TH1*>(hist->Clone("0"));
	else      oth->Add(hist);
      }
    }
    sim->Add(sip);
    tot->SetMarkerStyle(25);
    tot->SetMarkerSize(1.5);
    oth->SetMarkerColor(kRed+1);
    oth->SetLineColor(kRed+1);
    oth->SetFillColor(kRed+1);
    oth->SetMarkerStyle(20);
    oth->SetTitle("Other");
    oth->SetName("0");
#if 0
    TString       nam;
    Color_t       col;
    Style_t       sty;
    U::PdgAttr(321,  nam, col, sty); kpm->SetTitle(nam);
    U::PdgAttr(310,  nam, col, sty); k0s->SetTitle(nam);
    U::PdgAttr(3122, nam, col, sty); lam->SetTitle(nam);
    U::PdgAttr(3112, nam, col, sty); sim->SetTitle(nam);
    U::PdgAttr(3312, nam, col, sty); xi ->SetTitle(nam);
#endif 
    THStack*      ret  = new THStack(subName,subName);
    TH1*          a[]  = { k0s, kpm, lam, sim, xi, oth, tot };
    for (int i = 0; i < 7; i++) {
      if (!a[i]) continue;
      a[i]->SetDirectory(0);
      ret->Add(a[i]);
    }
    return ret;
  }
  TH1* FindHist(THStack* s, const char* name)
  {
    return U::GetH1(s->GetHists(), name, false);
  }
  TH1* MakeMeas(const char* name, THStack* prim, THStack* seco, THStack* fake)
  {
    TH1* r = 0;
    TH1* a[] = { FindHist(prim, name),
		 FindHist(seco, name),
		 FindHist(fake, name) };
    for (int i = 0; i < 3; i++) {
      if (!a[i]) continue;
      if (!r) {
	r = static_cast<TH1*>(a[i]->Clone());
	r->SetDirectory(0);
	// r->Sumw2();
	r->Reset();
      }
      r->Add(a[i]);
    }
    return r;
  }
  THStack* MakeMeas(THStack* prim, THStack* seco, THStack* fake)
  {
    THStack* ret = new THStack("measured", "measured");
    ret->Add(MakeMeas("310",   prim, seco, fake));
    ret->Add(MakeMeas("321",   prim, seco, fake));
    ret->Add(MakeMeas("3122",  prim, seco, fake));
    ret->Add(MakeMeas("3112",  prim, seco, fake));
    ret->Add(MakeMeas("3312",  prim, seco, fake));
    ret->Add(MakeMeas("0",     prim, seco, fake));
    ret->Add(MakeMeas("total", prim, seco, fake));
    return ret;
  }
  THStack* Normalize(THStack* inp)
  {
    THStack* ret  = new THStack(Form("norm%s",inp->GetName()),
				Form("%s normalised", inp->GetTitle()));
    TIter    next(inp->GetHists());
    TH1*     hist = 0;
    while ((hist = static_cast<TH1*>(next()))) {
      Double_t err = 0;
      Double_t ing = hist->IntegralAndError(1,hist->GetNbinsX(), err);
      TH1*     hir = static_cast<TH1*>(hist->Clone());
      hir->SetDirectory(0);
      // hir->Sumw2();
      Printf("%10s integral %f +/- %f", hist->GetTitle(), ing, err);
      // U::Scale(hir, 1/ing, err);
      hir->Scale(1/ing);
      ret->Add(hir);
    }
    return ret;
  }
  THStack* ToOther(THStack* inp)
  {
    TH1*     o    = FindHist(inp, "0");
    TH1*     t    = FindHist(inp, "total");
    THStack* ret  = new THStack(Form("rat%s",inp->GetName()),
				Form("%s ratios", inp->GetTitle()));
    TIter    next(inp->GetHists());
    TH1*     hist = 0;
    while ((hist = static_cast<TH1*>(next()))) {
      // if (hist == o) continue;
      if (hist == t) {
#if 0
	Printf("%10s: integral %f \t-> %f",
	       hist->GetTitle(),
	       hist->Integral(1,hist->GetNbinsX()),
	       t->Integral(1,t->GetNbinsX()));
#endif 
	continue;
      }
      TH1*     hir = static_cast<TH1*>(hist->Clone());
      hir->Reset();
      hir->Divide(hist,t);
      hir->SetDirectory(0);
      hir->Scale(25./hir->Integral(1,hir->GetNbinsX()), "width");
      ret->Add(hir);
#if 0
      Printf("%10s: integral %f \t-> %f",
	     hist->GetTitle(),
	     hist->Integral(1,hist->GetNbinsX()),
	     hir->Integral(1,hir->GetNbinsX(),"width"));
#endif 
    }
    return ret;

  }
  
};


void DeltaShapes(const char* what="hijing.root", Double_t c1=48, Double_t c2=52)
{
  if (!gROOT->GetClass("AliTrackletAODUtils")) {
    Printf("Loading utilities");    
    gROOT->LoadMacro("$ANA_SRC/dndeta/tracklets3/AliTrackletAODUtils.C+g");
  }

  Computer c;
  c.Run(what, c1, c2);
}
