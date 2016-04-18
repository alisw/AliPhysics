

void DrawDeltas2(const char* filename, Double_t c1, Double_t c2)
{
  TFile* file = TFile::Open(filename, "READ");
  if (!file) {
    Warning("DrawDeltas2", "File %s couldn't be opened", filename);
    return;
  }

  TString dname; dname.Form("cent%06.2f_%06.2f", c1, c2);
  dname.ReplaceAll(".", "d");
  TDirectory* d = file->GetDirectory(dname);
  if (!d) {
    Warning("DrawDeltas2", "Directory %s not found in %s",
	    dname.Data(), file->GetName());
    return;
  }

  TDirectory* det = d->GetDirectory("details");
  if (!det) {
    Warning("DrawDeltas2", "Directory details not found in %s",
	    d->GetName());
    d->ls();
    return;
  }

  TObject* o = det->Get("deltas");
  if (!o) {
    Warning("DrawDeltas2", "Object deltas not found in %s",
	    det->GetName());
    return;
  }

  if (!o->IsA()->InheritsFrom(THStack::Class())) {
    Warning("DrawDeltas2", "Object %s is not a THStack, but a %s",
	    o->GetName(), o->ClassName());
    return;
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
  if (det->Get("deltaInt")) {
    t.Append(" #it{k}(#eta),"); nm.Append("keta_");
  }
  else {
    t.Append(" #it{k}#equiv1,"); nm.Append("kunit_");
  }
  t.Append(Form(" %5.2f - %5.2f%%", c1, c2));
  nm.Append(dname);
  
  TCanvas* c = new TCanvas(nm,t,1000, 800);
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

  TH1* realMeas = static_cast<TH1*>(s->GetHists()->At(0));
  TH1* simMeas  = static_cast<TH1*>(s->GetHists()->At(1));
  TH1* simSum   = static_cast<TH1*>(s->GetHists()->At(n-3));

  TH1* measRat = static_cast<TH1*>(simMeas->Clone("simMeas"));
  measRat->Divide(realMeas);
  measRat->SetTitle(Form("%s/%s", simMeas->GetTitle(), realMeas->GetTitle()));

  TH1* sumRat = static_cast<TH1*>(simSum->Clone("simSum"));
  sumRat->Divide(realMeas);
  sumRat->SetTitle(Form("%s/%s", simSum->GetTitle(), realMeas->GetTitle()));

  THStack* ratios = new THStack("ratios", "");
  ratios->Add(measRat);
  ratios->Add(sumRat);

  p = c->cd(2);
  p->SetRightMargin(0.01);
  p->SetLogx();
  // p->SetLogy();
  p->SetGridx();
  p->SetGridy();
  p->SetTicks();
  ratios->SetMinimum(0.8);
  ratios->SetMaximum(1.2);
  ratios->Draw("nostack");
  h   = ratios->GetHistogram();
  h->SetXTitle("#Delta");
  h->SetYTitle("Ratio to measured");

  l = p->BuildLegend(.11,.11,.7,.35);
  // l->SetNColumns(2);

  p->Modified();
  c->Modified();
  c->Update();
  c->cd();
  c->SaveAs(Form("%s.png", nm.Data()));
}

    
