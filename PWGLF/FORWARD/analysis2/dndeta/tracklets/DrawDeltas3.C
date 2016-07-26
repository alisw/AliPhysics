

THStack* GetDeltas(const char* filename, Double_t c1, Double_t c2,
		   Bool_t& k)
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
  k = (det->Get("scalar") != 0);
  
  TObject* o = det->Get("deltas");
  /// file->Close();
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

  return s;
}

void DrawStack(TVirtualPad* p, THStack* s, Double_t min, Double_t max)
{
  p->SetLogx();
  p->SetLogy();
  p->SetGridx();
  p->SetGridy();
  p->SetTicks();

  s->SetMaximum(max);
  s->SetMinimum(min);
  s->Draw("nostack");

  TH1*     h   = s->GetHistogram();
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
  TLegend* l = p->BuildLegend(p->GetLeftMargin(),
			      p->GetBottomMargin(),
			      .7,
			       p->GetBottomMargin()+.35);
  l->SetNColumns(2);
  Int_t    n = l->GetListOfPrimitives()->GetEntries();
  ((TNamed*)(l->GetListOfPrimitives()->At(n-2)))->SetTitle("Signal cut");
  ((TNamed*)(l->GetListOfPrimitives()->At(n-1)))->SetTitle("Background cut");
  // l->SetBorderSize(0);
  // l->SetFillStyle(0);
  p->Modified();

  // s->GetHists()->Print();
  
}
void DrawRatios(TVirtualPad* p, THStack* s, Double_t min=0.8, Double_t max=1.22)
{
  p->SetLogx();
  // p->SetLogy();
  p->SetGridx();
  p->SetGridy();
  p->SetTicks();

  TH1* realMeas = static_cast<TH1*>(s->GetHists()->At(0));
  TH1* simMeas  = static_cast<TH1*>(s->GetHists()->At(1));
  TH1* realInj  = static_cast<TH1*>(s->GetHists()->At(2));
  TH1* simInj   = static_cast<TH1*>(s->GetHists()->At(3));

  TH1* ratioMeas = static_cast<TH1*>(simMeas->Clone("ratioMeas"));
  TH1* ratioInj  = static_cast<TH1*>(simInj ->Clone("ratioInj"));
  ratioMeas->SetDirectory(0);
  ratioInj ->SetDirectory(0);
  ratioMeas->SetTitle(Form("%s/%s",simMeas->GetTitle(),realMeas->GetTitle()));
  ratioInj ->SetTitle(Form("%s/%s",simInj ->GetTitle(),realInj ->GetTitle()));
  ratioMeas->SetTitle("Measured");
  ratioInj ->SetTitle("Injection");
  ratioInj ->SetMarkerColor(realMeas->GetMarkerColor());
  ratioMeas->Divide(realMeas); 
  ratioInj ->Divide(realInj );
  ratioMeas->SetYTitle("Sim./Real");
  ratioInj ->SetYTitle("Sim./Real");

  THStack* r = new THStack("ratios", "");
  r->Add(ratioMeas);
  r->Add(ratioInj);
  r->SetMinimum(min);
  r->SetMaximum(max);

  r->Draw("nostack");
  TH1* h = r->GetHistogram();
  h->SetXTitle(ratioMeas->GetXaxis()->GetTitle());
  h->SetYTitle(ratioMeas->GetYaxis()->GetTitle());

  TLegend* l = p->BuildLegend(p->GetLeftMargin(),
			      p->GetBottomMargin(),
			      1-p->GetRightMargin(),
			      p->GetBottomMargin()+.15);
  l->SetNColumns(2);
  
  p->Modified();
}
  
void DrawTitle(TVirtualPad* c,
	       TVirtualPad* p,
	       const char*  file,
	       Double_t     c1,
	       Double_t     c2,
	       Bool_t       k)
{
  TString t(file);
  t.ReplaceAll("results/", "");
  t.ReplaceAll("combine_","");
  t.ReplaceAll("_0x3.root", "");
  t.Prepend("Weights: ");
  if (k) t.Append(", #it{k}(#eta),");
  else   t.Append(", #it{k}#equiv1,");
  t.Append(Form(" %4.1f - %4.1f%%", c1, c2));

  c->cd();
  Double_t xm = p->GetXlowNDC() + p->GetWNDC()/2;
  TLatex*  tit = new TLatex(xm, 0.99, t);
  tit->SetTextFont(42);
  tit->SetTextAlign(23);
  tit->SetTextSize(0.03);
  tit->Draw();

  c->Modified();
}
    
  

TCanvas* DrawDeltas3(const char* file1, const char* file2,
		     Double_t    c1=0,  Double_t    c2=5)
{
  Bool_t   k1  = false;
  Bool_t   k2  = false;
  THStack* s1  = GetDeltas(file1, c1, c2, k1);
  THStack* s2  = GetDeltas(file2, c1, c2, k2);
  Double_t min = TMath::Min(s1->GetMinimum("nostack"),
			    s2->GetMinimum("nostack"));
  Double_t max = TMath::Max(s1->GetMaximum("nostack"),
			    s2->GetMaximum("nostack"));

  Int_t cW = 1200;
  Int_t cH =  900;
  TCanvas* c = new TCanvas("c", "C", cW, cH);
  c->SetTopMargin(0.1);
  c->SetRightMargin(0.01);
  c->Divide(2,2,0,0);

  TVirtualPad* p = c->cd(1);
  DrawStack(p, s1, min, max);

  p = c->cd(2);
  p->SetRightMargin(0.01);
  DrawStack(p, s2, min, max);

  p = c->cd(3);
  DrawRatios(p, s1);

  p = c->cd(4);
  p->SetRightMargin(0.01);
  DrawRatios(p, s2);

  DrawTitle(c, c->GetPad(1), file1, c1, c2, k1);
  DrawTitle(c, c->GetPad(2), file2, c1, c2, k2);

  return c;
}

void DrawDeltas3(UShort_t    flags,
		 const char* var, 
		 Double_t    c1=0,
		 Double_t    c2=5)
{
  TCanvas* c = DrawDeltas3(Form("results/combine_none_0x%x.root",flags),
			   Form("results/combine_%s_0x%x.root",var,flags),
			   c1, c2);
  c->Modified();
  c->Update();
  c->cd();
  c->SaveAs(Form("plots/deltas2_%s_0x%x.png", var, flags));
}
