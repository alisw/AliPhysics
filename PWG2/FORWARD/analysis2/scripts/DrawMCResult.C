THStack*
GetStack(const TList& forward,  const char* sub, const char* name)
{
  TList* lsub = static_cast<TList*>(forward.FindObject(sub));
  if (!lsub) { 
    Warning("GetStack", "Sub list %s not found in %s", sub, forward.GetName());
    return 0;
  }
  THStack* ret = static_cast<THStack*>(lsub->FindObject(name));
  if (!ret) 
    Warning("GetStack" "Stack %s not found in %s", name, sub);
  return ret;
}

TH1* 
Rebin(TH1* h, Int_t rebin)
{
  if (rebin <= 1) return h;
  h->Rebin(rebin);
  h->Scale(1. / rebin);
  return h;
}

void
DrawMCResult(const char* filename="forward.root", Int_t rebin=1)
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawMCResult", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get("Forward"));
  if (!forward) { 
    Error("DrawMCResult", "List Forward not found in %s", filename);
    return;
  }
  THStack* res    = GetStack(*forward, "ringResults", "all");
  THStack* mcRes  = GetStack(*forward, "mcRingResults", "all");
  THStack* deltas = GetStack(*forward, "fmdSharingFilter", "sums");
  THStack* nchs   = GetStack(*forward, "fmdDensityCalculator", "sums");
  THStack* prims  = GetStack(*forward, "fmdCorrector", "sums");
  
  TH1* sumEta = static_cast<TH1*>(forward->FindObject("mcSumEta"));
  if (!sumEta) { 
    Warning("DrawMCResults", "mcSumEta not found in Forward");
  }



  TCanvas* c = new TCanvas("c", "C", 900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleStyle(0);
  // gStyle->SetTitleColor(kBlack);

  THStack* all = new THStack("all", "Analysis steps");
  if (res) {
    TIter next(res->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) all->Add(Rebin(h,rebin));
  }
  if (mcRes) {
    TIter next(mcRes->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) all->Add(Rebin(h,rebin));
  }
  if (deltas) {
    TIter next(deltas->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) { 
      h->SetMarkerStyle(25);
      all->Add(Rebin(h,rebin));
    }
  }
  if (nchs) {
    TIter next(nchs->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) { 
      h->SetMarkerStyle(21);
      all->Add(Rebin(h,rebin));
    }
  }
  if (prims) {
    TIter next(prims->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) { 
      h->SetMarkerStyle(22);
      all->Add(Rebin(h,rebin));
    }
  }
  if (sumEta) all->Add(sumEta);
  all->Draw("nostack");
  all->GetHistogram()->SetXTitle("#eta");
  all->GetHistogram()->SetYTitle("signal");
  all->GetHistogram()->GetXaxis()->SetLabelFont(132);
  all->GetHistogram()->GetXaxis()->SetTitleFont(132);
  all->GetHistogram()->GetYaxis()->SetLabelFont(132);
  all->GetHistogram()->GetYaxis()->SetTitleFont(132);
  c->SetGridx();

  TLegend* l = new TLegend(.33, .2, .53, .9);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetNColumns(1);
  l->SetTextFont(132);
  Int_t i = 0;
  if (res) { 
    TIter next(res->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      TLegendEntry* e = l->AddEntry(Form("dummy%02d", i++),h->GetTitle(),"pl");
      e->SetMarkerStyle(20);
      e->SetMarkerColor(h->GetMarkerColor());
    }
    TLegendEntry* e1 = l->AddEntry(Form("dummy%02d", i++),
				   "dN_{ch}/d#eta","pl");
    e1->SetMarkerStyle(20);
    e1->SetMarkerColor(kBlack);
  }
  if (deltas) { 
    TLegendEntry* e1 = l->AddEntry(Form("dummy%02d", i++),
				   "#sum_{} #Delta/#Delta_{mip}","pl");
    TH1* h = static_cast<TH1*>(deltas->GetHists()->At(0));
    e1->SetMarkerStyle(h->GetMarkerStyle());
    e1->SetMarkerColor(kBlack);
  }
  if (nchs) { 
    TLegendEntry* e1 = l->AddEntry(Form("dummy%02d", i++),
				   "#sum_{} N_{ch,incl}","pl");
    TH1* h = static_cast<TH1*>(nchs->GetHists()->At(0));
    e1->SetMarkerStyle(h->GetMarkerStyle());
    e1->SetMarkerColor(kBlack);
  }
  if (prims) { 
    TLegendEntry* e1 = l->AddEntry(Form("dummy%02d", i++),
				   "#sum_{} N_{ch,primary}","pl");
    TH1* h = static_cast<TH1*>(prims->GetHists()->At(0));
    e1->SetMarkerStyle(h->GetMarkerStyle());
    e1->SetMarkerColor(kBlack);
  }

  if (mcRes) { 
    TLegendEntry* e1 = l->AddEntry(Form("dummy%02d", i++),"Track-Refs", "pl");
    TH1* h = static_cast<TH1*>(mcRes->GetHists()->At(0));
    e1->SetMarkerStyle(h->GetMarkerStyle());
    e1->SetMarkerColor(kBlack);
  }

  if (sumEta) l->AddEntry(sumEta);
  l->Draw();

}

  
  
 
