void
DrawFits(const char* fname="AnalysisResults.root")
{
  TFile* file = TFile::Open(fname, "READ");
  if (!file) {
    Error("DrawFits", "Couldn't open %s", fname);
    return;
  }
    
  TList* forward = 
    static_cast<TList*>(file->Get("PWG2forwardDnDeta/Forward"));
  if (!forward) { 
    Error("DrawFits", "Couldn't get forward list from %s", fname);
    return;
  }

  TList* fitter = 
    static_cast<TList*>(forward->FindObject("fmdEnergyFitter"));
  if (!fitter) { 
    Error("DrawFits", "Couldn't get fitter folder");
    return;
  }
  
  TList stacks;
  stacks.Add(fitter->FindObject("chi2"));
  stacks.Add(fitter->FindObject("c"));
  stacks.Add(fitter->FindObject("mpv"));
  stacks.Add(fitter->FindObject("w"));
  stacks.Add(fitter->FindObject("n"));
  Int_t i=2;
  while (true) { 
    TObject* o = static_cast<THStack*>(fitter->FindObject(Form("a%d",i++)));
    if (!o) break;
    Info("DrawFits", "Adding %s", o->GetName());
    stacks.Add(o);
  }
  // stacks.ls();
  Int_t nMax = stacks.GetEntries();
  for (Int_t i = nMax-1; i > 4; i--) { 
    THStack* stack   = static_cast<THStack*>(stacks.At(i));
    TIter    nextH(stack->GetHists());
    TH1*     hist    = 0;
    Bool_t   hasData = kFALSE;
    while ((hist = static_cast<TH1*>(nextH()))) 
      if (hist->Integral() > 0) hasData = kTRUE;
    if (!hasData) nMax--;
  }

  gStyle->SetOptTitle(0);
#if 1
  TCanvas* c = new TCanvas("c", "C",800,800);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(0);
  c->Divide(1, nMax,0,0);

  TIter next(&stacks);
  THStack* stack = 0;
  i = 1;
  while ((stack = static_cast<THStack*>(next()))) {
    if (i > nMax) break;
    TVirtualPad* p = c->cd(i);
    p->SetFillColor(0);
    p->SetFillStyle(0);
    p->SetGridx();
    stack->Draw("nostack");
    stack->GetHistogram()->SetYTitle(stack->GetTitle());
    stack->GetHistogram()->SetXTitle("#eta");
    TAxis* yaxis = stack->GetHistogram()->GetYaxis();
    yaxis->SetTitleSize(0.15);
    yaxis->SetLabelSize(0.08);
    yaxis->SetTitleOffset(0.35);
    yaxis->SetNdivisions(10);
    TAxis* xaxis = stack->GetHistogram()->GetXaxis();
    xaxis->SetTitleSize(0.15);
    xaxis->SetLabelSize(0.08);
    xaxis->SetTitleOffset(0.35);
    xaxis->SetNdivisions(320);
    i++;
    p->cd();
  }
#endif
    
  gStyle->SetOptFit(111111);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  TCanvas* c1 = new TCanvas("c1", "c1", 800,800);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->Divide(1, 5,0,0);

  const char* dets[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3I", "FMD3O", 0 };
  for (Int_t i = 0; i < 5; i++) { 
    TVirtualPad* p = c1->cd(i+1);
    p->SetGridx();
    p->SetFillColor(0);
    p->SetFillStyle(0);
    TList* d = static_cast<TList*>(fitter->FindObject(dets[i]));
    if (!d) { 
      Warning("DrawFits", "List %s not found", dets[i]);
      continue;
    }
    TH1* edist = static_cast<TH1*>(d->FindObject(Form("%s_edist", dets[i])));
    if (!edist) {
      Warning("DrawFits", "Histogram %s_edist not found", dets[i]);
      continue;
    }
    edist->Draw();
    TF1*   f = 0;
    TIter  nextF(edist->GetListOfFunctions());
    while ((f = static_cast<TF1*>(nextF()))) {
      Double_t chi2 = f->GetChisquare();
      Int_t    ndf  = f->GetNDF();
      Printf("%s %s:\n  Range: %f-%f\n" 
             "chi^2/ndf= %f / %d = %f", 
	     edist->GetName(), f->GetName(), 
             f->GetXmin(), f->GetXmax(), chi2, ndf, 
	     (ndf > 0) ? chi2/ndf : 0);
      for (Int_t j = 0; j < f->GetNpar(); j++) { 
	Printf("  %-20s : %9.4f +/- %9.4f", 
	       f->GetParName(j), f->GetParameter(j), f->GetParError(j));
      }
    }
    p->cd();
  }
  c1->cd();
}
