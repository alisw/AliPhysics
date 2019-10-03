/** 
 * Read in a single file and display Delta distributions for selected
 * centrality class.
 * 
 * @param filename File to read 
 * @param c1       Centrality class 
 * @param c2       Centrality class 
 */
TCanvas* DrawDeltas2(const char* filename, Double_t c1=0, Double_t c2=5)
{
  TFile* file = TFile::Open(filename, "READ");
  if (!file) {
    Warning("DrawDeltas2", "File %s couldn't be opened", filename);
    return 0;
  }

  TString dname; dname.Form("cent%06.2f_%06.2f", c1, c2);
  dname.ReplaceAll(".", "d");
  TDirectory* d = file->GetDirectory(dname);
  if (!d) {
    Warning("DrawDeltas2", "Directory %s not found in %s",
	    dname.Data(), file->GetName());
    return 0;
  }

  TDirectory* det = d->GetDirectory("details");
  if (!det) {
    Warning("DrawDeltas2", "Directory details not found in %s",
	    d->GetName());
    d->ls();
    return 0;
  }

  TObject* o = det->Get("deltas");
  if (!o) {
    Warning("DrawDeltas2", "Object deltas not found in %s",
	    det->GetName());
    return 0;
  }

  if (!o->IsA()->InheritsFrom(THStack::Class())) {
    Warning("DrawDeltas2", "Object %s is not a THStack, but a %s",
	    o->GetName(), o->ClassName());
    return 0;
  }
  THStack* s = static_cast<THStack*>(o);

  TString t(filename);
  t.ReplaceAll("results/", "");
  t.ReplaceAll("combine_","");
  t.ReplaceAll("_0x3.root", "");
  TString nm(filename);
  nm.ReplaceAll(".root", "_");
  nm.ReplaceAll("results/", "plots/");
  nm.ReplaceAll("combine", "deltas");
  if (t.Contains("none")) 
    t.ReplaceAll("none", "No weights,");
  else
    t.Append(" weights,");
  if (det->Get("scalar")) {
    t.Append(" #it{k}(#eta),"); nm.Append("keta_");
  }
  else {
    t.Append(" #it{k}#equiv1,"); nm.Append("kunit_");
  }
  t.Append(Form(" %5.2f - %5.2f%%", c1, c2));
  nm.Append(dname);
  
  Int_t cW = 1200;
  Int_t cH =  800;
  TCanvas* c = new TCanvas(nm,t,cW, cH);
  c->SetTopMargin(0.10);
  c->SetRightMargin(0.01);
  c->Divide(1,2,0,0);

  TLatex* tit = new TLatex(0.55, 0.99, t);
  tit->SetTextFont(42);
  tit->SetTextAlign(23);
  tit->SetTextSize(0.03);
  tit->Draw();
  
  TVirtualPad* p = c->cd(1);
  p->SetRightMargin(0.01);
  p->SetLogx();
  p->SetLogy();
  p->SetGridx();
  p->SetGridy();
  p->SetTicks();
  s->Draw("nostack");

  TH1*     h   = s->GetHistogram();
  Double_t min = h->GetMinimum();
  Double_t max = h->GetMaximum();
  h->SetXTitle("#Delta");
  h->SetYTitle("Tracklets/event");

  TLine* ll = new TLine(1.5, min, 1.5, max);
  ll->SetLineColor(kGreen+2);
  ll->SetLineWidth(3);
  ll->Draw();

  TLine* lh = new TLine(5, min, 5, max);
  lh->SetLineColor(kRed+2);
  lh->SetLineWidth(3);
  lh->Draw();

  // TLegend* l = p->BuildLegend(.6,.75,.99,.99);
  TLegend* l = p->BuildLegend(.11,.0,.7,.35);
  l->SetNColumns(2);
  Int_t    n = l->GetListOfPrimitives()->GetEntries();
  ((TNamed*)(l->GetListOfPrimitives()->At(n-2)))->SetTitle("Signal cut");
  ((TNamed*)(l->GetListOfPrimitives()->At(n-1)))->SetTitle("Background cut");
  // l->SetBorderSize(0);
  // l->SetFillStyle(0);
  p->Modified();
  // s->GetHists()->Print();

  TH1* realMeas = static_cast<TH1*>(s->GetHists()->At(0));
  TH1* simMeas  = static_cast<TH1*>(s->GetHists()->At(1));
  TH1* realInj  = static_cast<TH1*>(s->GetHists()->At(2));
  TH1* simInj   = static_cast<TH1*>(s->GetHists()->At(3));
  TH1* simComb  = static_cast<TH1*>(s->GetHists()->At(4));

  TH1* measRat = static_cast<TH1*>(simMeas->Clone("ratioMeas"));
  measRat->Divide(realMeas);
  measRat->SetTitle("Sim./Real Measured");

  
  TH1* simInjRat = static_cast<TH1*>(simInj->Clone("ratioSimInj"));
  simInjRat->Divide(simComb);
  simInjRat->SetDirectory(0); 
  simInjRat->SetTitle("Sim. Inj./Combinatorics");

  TH1* realInjRat = static_cast<TH1*>(realInj->Clone("ratioRealInj"));
  realInjRat->Divide(simComb);
  realInjRat->SetDirectory(0); 
  realInjRat->SetTitle("Real Inj./Combinatorics");
  
  THStack* ratios = new THStack("ratios", "");
  ratios->Add(measRat);
  ratios->Add(realInjRat);
  ratios->Add(simInjRat);

  p = c->cd(2);
  p->SetRightMargin(0.01);
  p->SetLogx();
  // p->SetLogy();
  p->SetGridx();
  p->SetGridy();
  p->SetTicks();
  ratios->SetMinimum(0.5);
  ratios->SetMaximum(1.5);
  ratios->Draw("nostack");
  h   = ratios->GetHistogram();
  h->SetXTitle(realMeas->GetXaxis()->GetTitle());
  h->SetYTitle("Ratio");

  l = p->BuildLegend(.11,.11,.7,.35);
  // l->SetNColumns(2);

  p->Modified();
  c->Modified();
  c->Update();
  c->cd();

  return c;
}

    
void DrawDeltas2(UShort_t flags, const char* which)
{
  TCanvas* c =
    DrawDeltas2(Form("results/combine_%s_0x%x.root", which, flags), 0, 5);
  c->SaveAs(Form("plots/deltas_%s_0x%x.png", which, flags));
}

