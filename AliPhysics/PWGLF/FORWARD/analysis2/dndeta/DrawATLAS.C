//====================================================================
TCanvas* MakeCanvas(const char* name, const char* title,
		    const char* img=0, Double_t r=0)
{
  TImage* dimg = 0;
  if (img) {
    const char* full = Form("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta/%s",
			    img);
    TImage* dimg = TImage::Open(gSystem->ExpandPathName(full));
    if (!dimg) {
      Error("Failed to open image file: %s", img);
      return 0;
    }
  }
  UInt_t w = (dimg ? dimg->GetWidth() : 1000);
  UInt_t h = (dimg ? dimg->GetHeight() : 900);
  Printf("Canvas size: (%dx%d)", w, h);
  TCanvas* c = new TCanvas(name, title, w, h);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);

  if (dimg) {
    c->SetTopMargin(0);
    c->SetRightMargin(0);
    c->SetLeftMargin(0);
    c->SetRightMargin(0);
    c->cd();
    // dimg->Scale(c->GetWw(), c->GetWh());
    dimg->Draw("x");
    // dimg->SetDrawOption("X");

    TPad* p = new TPad("overlay", "Overlay", 0, 0, 1, 1, 0, 0, 0);
    p->SetFillStyle(0);
    p->SetFillColor(0);
    p->SetLineColor(0);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetNumber(10);
    c->cd();
    p->Draw();
    
    c->Modified();
    c->Update();
    c->cd();
    return c;
  }

  if (r < 1e-3) return c;

  c->Divide(1,2,0,0);
  TVirtualPad* p = c->cd(1);
  p->SetPad(0, r , 1, 1);
  p->SetBottomMargin(0.001);
  p->SetTopMargin(0.01);
  p->SetRightMargin(0.01);
  p->Modified();
  p->Update();
  p->cd();

  p = c->cd(2);
  p->SetPad(0, 0, 1, r);
  Double_t yp = r / (1-r);
  p->SetRightMargin(0.01);
  p->SetBottomMargin(p->GetBottomMargin()/yp);
  p->Modified();
  p->Update();
  p->cd();

  c->Modified();
  c->Update();
  c->cd();

  return c;
}

//====================================================================
TPair* GetUs(TLegend*& dl, Bool_t sumLast=false)
{
  if (!dl){
    dl = new TLegend(.12, .65, .35, .93);
    dl->SetTextAlign(32);
    dl->SetNColumns(1);
    dl->SetFillColor(0);
    dl->SetFillStyle(0);
    dl->SetBorderSize(0);
  }
  e = dl->AddEntry("dummy", "PWGLF/GEO", "p");
  e->SetMarkerStyle(20);
  e->SetMarkerSize(1.5);
  e = dl->AddEntry("dummy", "PWGPP/GLOB", "p");
  e->SetMarkerStyle(21);
  e->SetMarkerSize(1.5);
  
  TObjArray unique;
  const char* trigs[] = { "CENTV0A", 0 };
  const char* exps[]  = { "ALICE", "WIP", 0 };
  TPair* pair = Drawer::GetDataOther(dl, unique, "pPb", 5023, trigs,
				     exps, "e2", true);
  if (!pair || !pair->Key()) {
    Error("", "Didn't get the results");
    return;
  }
  THStack*     data  = static_cast<THStack*>(pair->Key());
  TMultiGraph* other = static_cast<TMultiGraph*>(pair->Value());
  TString      tit = Form("pPb %s %s", data->GetTitle(), trigs[0]);
  other->GetListOfGraphs()->Sort();
  Drawer::AddSystematics(data, 0.076);
  data->SetTitle(tit);

  if (!sumLast) return pair;

  TList*        hl    = data->GetHists();
  TList*        ol    = other->GetListOfGraphs();
  TH1*          h6080 = static_cast<TH1*>(hl->At(hl->GetEntries()-1));
  TH1*          h8099 = static_cast<TH1*>(hl->At(hl->GetEntries()-2));
  TGraphErrors* o6080 = static_cast<TGraphErrors*>(ol->At(ol->GetEntries()-1));
  TGraphErrors* o8099 = static_cast<TGraphErrors*>(ol->At(ol->GetEntries()-2));

  hl->Remove(h6080);
  hl->Remove(h8099);
  ol->Remove(o6080);
  ol->Remove(o8099);
  TH1*          h6099 = Average("h6099", "(*) 60%-90%", h6080, .25, h8099, .75);
  TGraphErrors* o6099 = Average("o6099", "60%-90%", o6080, .25, o8099, .75);
  hl->Add(h6099);
  ol->Add(o6099);

  TList* el = dl->GetListOfPrimitives();
  el->Remove(el->Last());
  el->Remove(el->Last());
  dl->AddEntry(h6099, h6099->GetTitle(), "pl");
  
  return pair;
}

void DrawUs(TPair* pair, TLegend* dl, Double_t max=105)
{
  THStack*     data  = static_cast<THStack*>(pair->Key());
  TMultiGraph* other = static_cast<TMultiGraph*>(pair->Value());
  
  data->Draw("nostack");
  data->GetHistogram()->SetXTitle("#it{#eta}");
  data->GetHistogram()->GetYaxis()->SetTitleOffset(1.1);
  data->GetHistogram()->SetYTitle("1/#it{N} d#it{N}_{ch}/d#it{#eta}");
  if (other) other->Draw("p");
  dl->Draw();

  data->SetMaximum(max);

  return pair;
}


//====================================================================  
void
DrawOrigATLAS(Bool_t overlay=false)
{
  TCanvas* c  = MakeCanvas("orig", "ATLAS & ALICE",
			   (overlay ? "ATLAS_pPb_dndeta.png" : 0));

  TVirtualPad* p = (overlay ? c->GetPad(10) : c);
  if (!overlay)
    p->SetLeftMargin(0.15);
  else {
    c->SetTopMargin(0.2);
    p->SetLeftMargin(0.117);
    p->SetRightMargin(0.016);
    p->SetTopMargin(0.01);
    p->SetBottomMargin(0.115);
    
  }
  p->cd();

  TLegend* dl = 0;
  TPair* pair = GetUs(dl);
  DrawUs(pair, dl, overlay ? 79.5 : 105);
  
  // al->AddEntry("dummy","(read-off arXiv/1403.5738v1)", "");
  TLegend* al = new TLegend(.79, .65, .97, .93);
  al->SetTextAlign(32);
  al->SetNColumns(1);
  al->SetFillColor(0);
  al->SetFillStyle(0);
  al->SetBorderSize(0);
  e = al->AddEntry("dummy", "ATLAS", "p");
  e->SetMarkerStyle(22);
  e->SetMarkerSize(1.5);
  
  TLatex* ltx = new TLatex(.97, .64, "(read off arXiv/1403.5738v1)");
  ltx->SetTextAlign(33);
  ltx->SetTextColor(kGray+2);
  ltx->SetNDC();
  ltx->SetTextSize(0.015);
  ltx->Draw();

  TMultiGraph* mg = ATLASpPb(false, al);
  mg->Draw("p");
  al->Draw();

  if (overlay) {
    THStack* s = static_cast<THStack*>(pair->Key());
    // mg->GetHistogram()->SetMaximum(79);
    s->GetHistogram()->SetMinimum(0);
    s->GetHistogram()->GetXaxis()->SetRangeUser(-3,3);
    TFrame* f =
      static_cast<TFrame*>(p->GetListOfPrimitives()->FindObject("TFrame"));
    if (!f) {
      Warning("", "No frame");
    }
    else {
      f->SetFillStyle(0);
    }
    p->Modified();
    p->Update();
    p->cd();
  }
  
  if (!overlay) c->Print("plots/alice_atlas_v0a.pdf");
  c->Print("plots/alice_atlas_v0a.png");
  c->Print("plots/alice_atlas_v0a.C");
}

//====================================================================  
TMultiGraph* Ratio2Peripheral(TMultiGraph* mg, const char* func=0)
{
  TGraph*      perip = static_cast<TGraph*>(mg->GetListOfGraphs()->Last());
  if (!perip) {
    Error("", "No periphiral bin found!");
    return;
  }
  TMultiGraph* rmg = new TMultiGraph;
  rmg->SetName(mg->GetName());
  rmg->SetTitle(mg->GetTitle());

  TGraph* g = 0;
  TIter next(mg->GetListOfGraphs());
  // TList funcs;
  while ((g = static_cast<TGraph*>(next()))) {
    if (g == perip) continue;
    TGraph* r = Drawer::GOverG(g, perip, 0, 0, 0);
    r->SetFillStyle(3002);
    r->SetLineColor(g->GetMarkerColor());
    rmg->Add(r, "p 1 z");
    if (!func) continue;
    TFitResultPtr fit = r->Fit(func, "Q0S+");
    TF1* f = r->GetListOfFunctions()->FindObject(func)->Clone();
    Printf("%s | %f | %f +/- %f | %f +/- %f | %f +/- %f | %p",
	   r->GetName(),
	   fit->Chi2()/fit->Ndf(),
	   fit->Parameter(0), fit->ParError(0),
	   fit->Parameter(1), fit->ParError(1),
	   fit->Parameter(2), fit->ParError(2), f);
    rmg->GetListOfFunctions()->Add(f);
  }
  return rmg;
}
//____________________________________________________________________  
THStack* Ratio2Peripheral(THStack* mg, const char* func=0, TList* l=0)
{
  TH1*  perip = static_cast<TH1*>(mg->GetHists()->Last());
  if (!perip) {
    Error("", "No periphiral bin found!");
    return;
  }
  THStack* rmg = new THStack(mg->GetName(), mg->GetTitle());

  TH1* g = 0;
  TIter next(mg->GetHists());
  // TList funcs;
  while ((g = static_cast<TH1*>(next()))) {
    if (g == perip) continue;
    TH1* r = static_cast<TH1*>(g->Clone());
    r->SetDirectory(0);
    r->Divide(perip);
    r->SetFillStyle(3002);
    r->SetLineColor(g->GetMarkerColor());
    rmg->Add(r);
    if (!func) continue;

    TFitResultPtr fit = r->Fit(func, "Q0S+");
    TF1* f = r->GetListOfFunctions()->FindObject(func); // ->Clone();
    Printf("%s | %f | %f +/- %f | %f +/- %f | %f +/- %f | %p",
	   r->GetName(),
	   fit->Chi2()/fit->Ndf(),
	   fit->Parameter(0), fit->ParError(0),
	   fit->Parameter(1), fit->ParError(1),
	   fit->Parameter(2), fit->ParError(2), f);
    if (l) l->Add(f->Clone());
    // rmg->GetListOfFunctions()->Add(f);
  }
  return rmg;
}

//====================================================================  
void DrawOrigATLASRatio(Bool_t overlay=false)
{
  TCanvas* c = MakeCanvas("ratio", "ATLAS dN/deta|cp",
			  (overlay ? "ATLAS_pPb_dndeta_cp.png" : 0));
  TVirtualPad* p = (overlay ? c->GetPad(10) : c);
  if (!overlay)
    p->SetLeftMargin(0.15);
  else {
    p->SetLeftMargin(0.123);
    p->SetRightMargin(0.0175);
    p->SetTopMargin(0.005);
    p->SetBottomMargin(0.11);
    
  }
  p->cd();
  TLegend* al = new TLegend(.18, .65, .4, .93);
  al->SetTextAlign(32);
  al->SetNColumns(1);
  al->SetFillColor(0);
  al->SetFillStyle(0);
  al->SetBorderSize(0);
  e = al->AddEntry("dummy", "ATLAS", "p");
  e->SetMarkerStyle(22);
  e->SetMarkerSize(1.5);
  al->GetListOfPrimitives()->Remove(al->GetListOfPrimitives()->Last());
  
  TMultiGraph* mg  = ATLASpPb(false, al);
  TMultiGraph* rmg = Ratio2Peripheral(mg, "pol2");
  rmg->SetTitle(Form("%s pPb ratio to peripheral", mg->GetTitle()));
  rmg->Draw("ap");
  rmg->GetListOfFunctions()->Draw("same");

  const char* dndeta = "1/#it{N} d#it{N}_{ch}/d#it{#eta}";
  rmg->GetHistogram()->GetXaxis()->SetTitle("#it{#eta}");
  rmg->GetHistogram()->SetYTitle(Form("#frac{%s|_{#it{c}}}{%s|_{60%%-90%%}}",
				      dndeta, dndeta));
  rmg->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);

  al->GetListOfPrimitives()->Remove(al->GetListOfPrimitives()->Last());
  al->Draw();
  
  TLatex* ltx = new TLatex(.97, .14, "(read off arXiv/1403.5738v1)");
  ltx->SetTextAlign(33);
  ltx->SetTextColor(kGray+2);
  ltx->SetNDC();
  ltx->SetTextSize(0.015);
  ltx->Draw();

  if (overlay) {
    rmg->GetHistogram()->SetMaximum(9.9);
    rmg->GetHistogram()->SetMinimum(0);
    TFrame* f =
      static_cast<TFrame*>(p->GetListOfPrimitives()->FindObject("TFrame"));
    if (!f) {
      Warning("", "No frame");
    }
    else {
      f->SetFillStyle(0);
    }
    p->Modified();
    p->Update();
    p->cd();
  }
  c->Modified();
  c->Update();
  c->cd();
  
  
  if (!overlay) c->Print("plots/alice_atlas_cp_v0a.pdf");
  c->Print("plots/alice_atlas_cp_v0a.png");
  c->Print("plots/alice_atlas_cp_v0a.C");
}

//====================================================================  
void
Add2Ratio(TMultiGraph* mg, TH1* num, TGraphErrors* den)
{
  TGraph* low = 0;
  TGraph* upp = 0;
  Drawer::ErrorGraphs(den, low, upp);

  TGraph* gnu = Drawer::H2G(num, 0, 1);

  TGraph* rat = Drawer::GOverG(gnu, den, low, upp);
  rat->SetMarkerStyle(num->GetMarkerStyle());
  rat->SetMarkerSize(num->GetMarkerSize());
  mg->Add(rat);
}
//____________________________________________________________________  
void
Add2Ratio(TMultiGraph* mg, TGraphErrors* num, TGraphErrors* den)
{
  TGraph* low = 0;
  TGraph* upp = 0;
  Drawer::ErrorGraphs(den, low, upp);

  TGraph* rat = Drawer::GOverG(num, den, low, upp);
  rat->SetMarkerStyle(num->GetMarkerStyle());
  rat->SetMarkerSize(num->GetMarkerSize());
  mg->Add(rat);
}
//____________________________________________________________________  
void 
Add2Ratio(TMultiGraph* mg, TGraphErrors* g, TList* hl, Int_t idx)
{
  Add2Ratio(mg, static_cast<TH1*>(hl->At(idx)), g);

}
//____________________________________________________________________  
void 
Add2Ratio(TMultiGraph* mg, TGraphErrors* g, TList* gl, const char* name)
{
  Add2Ratio(mg, static_cast<TGraphErrors*>(gl->FindObject(name)),g);
}
  
//____________________________________________________________________  
void DrawAvgATLAS()
{
  Double_t yp = .3;
  TCanvas* c  = MakeCanvas("averaged", "ATLAS & ALICE (averaged)", 0, yp);
  c->cd(1);
  
  TLegend*     dl    = new TLegend(.12, .6, .35, .93);
  dl->SetTextAlign(32);
  dl->SetNColumns(1);
  dl->SetFillColor(0);
  dl->SetFillStyle(0);
  dl->SetBorderSize(0);
  TLegendEntry* e = dl->AddEntry("dummy", "ATLAS (resummed)", "p");
  e->SetMarkerStyle(22);
  e->SetMarkerSize(1.7);
  
  TPair*       pair  = GetUs(dl, true);
  THStack*     data  = static_cast<THStack*>(pair->Key());
  TMultiGraph* other = static_cast<TMultiGraph*>(pair->Value());
  DrawUs(pair, dl, 81);
  
  TMultiGraph* mg = ATLASpPbResum(false);
  mg->Draw("p");
  TIter         nextG(mg->GetListOfGraphs());
  TIter         nextO(other->GetListOfGraphs());
  TIter         nextH(data->GetHists());
  TGraphErrors* g      = 0;
  TGraphErrors* o      = 0;
  TH1*          h      = 0;  
  TMultiGraph* ratios = new TMultiGraph();

  while ((g = static_cast<TGraphErrors*>(nextG())) &&
	 (o = static_cast<TGraphErrors*>(nextO())) &&
	 (h = static_cast<TH1*>(nextH()))) {
    g->SetMarkerColor(h->GetMarkerColor());
    g->SetLineColor(h->GetLineColor());
    Add2Ratio(ratios, h, g);
    Add2Ratio(ratios, o, g);
    Printf("%s %s %s", g->GetName(), o->GetName(), h->GetName());
  }

  TLatex* ltx = new TLatex(.15, .6, 
			   "#splitline{(ATLAS read off arXiv/1403.5738v1}"
			   "{Averaged w/centrality bin width)}");
  ltx->SetTextAlign(13);
  ltx->SetTextColor(kGray+2);
  ltx->SetNDC();
  ltx->SetTextSize(0.015);
  ltx->Draw();

  ltx->DrawLatex(.15, .56,
		 "#splitline{(*) For ALICE: Weighted sum}"
		 "{of 60-80% and 80-100%}");
  // l->Draw();
  c->cd(2);
  yp /= (1-yp);
  
  TGraphErrors* band = new TGraphErrors(2);
  band->SetPoint(0,-3.95,1);
  band->SetPoint(1,5.95,1);
  band->SetPointError(0,0,.05);
  band->SetPointError(1,0,.05);
  band->SetFillColor(kYellow-10);
  band->SetFillStyle(1001);
  band->SetLineStyle(2);
  ratios->GetListOfGraphs()->AddFirst(band, "3l");

  ratios->Draw("ap");
  TAxis* yaxis = ratios->GetHistogram()->GetYaxis();
  TAxis* xaxis = ratios->GetHistogram()->GetXaxis();
  yaxis->SetTitle("ALICE/ATLAS");
  yaxis->SetLabelSize(yaxis->GetLabelSize()/yp);
  yaxis->SetTitleSize(yaxis->GetTitleSize()/yp);
  yaxis->SetTitleOffset(yaxis->GetTitleOffset()*yp);
  xaxis->SetTitle("#eta");
  xaxis->Set(200, -4, 6);
  xaxis->SetLabelSize(xaxis->GetLabelSize()/yp);
  xaxis->SetTitleSize(xaxis->GetTitleSize()/yp);
  band->Draw("same l");

  c->Modified();
  c->Update();
  c->cd();

  c->Print("plots/alice_watlas_v0a.pdf");
  c->Print("plots/alice_watlas_v0a.C");
  c->Print("plots/alice_watlas_v0a.png");
}

//====================================================================  
void DrawAvgATLASRatio()
{
  Double_t yp = .3;
  TCanvas* c = MakeCanvas("ratioboth", "ATLAS & ALICE dN/deta|cp", 0, yp);
  c->SetLeftMargin(0.15);
  c->cd(1);

  TLegend*     dl    = new TLegend(.12, .6, .35, .93);
  dl->SetTextAlign(32);
  dl->SetNColumns(1);
  dl->SetFillColor(0);
  dl->SetFillStyle(0);
  dl->SetBorderSize(0);
  TLegendEntry* e = dl->AddEntry("dummy", "ATLAS (resummed)", "p");
  e->SetMarkerStyle(22);
  e->SetMarkerSize(1.7);
  TPair*       pair   = GetUs(dl, true);
  THStack*     data   = static_cast<THStack*>(pair->Key());
  TMultiGraph* other  = static_cast<TMultiGraph*>(pair->Value());
  TList        fdata;
  THStack*     rdata  = Ratio2Peripheral(data,  0); // "pol1", &fdata);
  TMultiGraph* rother = Ratio2Peripheral(other, 0); // "pol1");
  
  TMultiGraph* mg    = ATLASpPbResum(false);
  TMultiGraph* rmg   = Ratio2Peripheral(mg, 0); // "pol2");
  rmg->SetTitle(Form("%s pPb ratio to peripheral", mg->GetTitle()));

  rdata->Draw("nostack");
  fdata.Draw("same");
  rother->Draw("p");
  rother->GetListOfFunctions()->Draw("same");
  rmg->Draw("p");
  dl->GetListOfPrimitives()->Remove(dl->GetListOfPrimitives()->Last());
  dl->Draw();
  rdata->SetMaximum(9.1);
  
  const char* dndeta = "1/#it{N} d#it{N}_{ch}/d#it{#eta}";
  rdata->GetHistogram()->GetXaxis()->SetTitle("#it{#eta}");
  rdata->GetHistogram()->SetYTitle(Form("#frac{%s|_{#it{c}}}"
					"{%s|_{60%%-90%% (*)}}",
					dndeta, dndeta));
  rdata->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
  
  TLatex* ltx = new TLatex(.15, .6, 
			   "#splitline{(ATLAS read off arXiv/1403.5738v1}"
			   "{Averaged w/centrality bin width)}");
  ltx->SetTextAlign(13);
  ltx->SetTextColor(kGray+2);
  ltx->SetNDC();
  ltx->SetTextSize(0.015);
  ltx->Draw();

  ltx->DrawLatex(.15, .56,
		 "#splitline{(*) For ALICE: Weighted sum}"
		 "{of 60-80% and 80-100%}");

  TIter         nextG(rmg->GetListOfGraphs());
  TIter         nextO(rother->GetListOfGraphs());
  TIter         nextH(rdata->GetHists());
  TGraphErrors* g      = 0;
  TGraphErrors* o      = 0;
  TH1*          h      = 0;  
  TMultiGraph* ratios = new TMultiGraph();

  while ((g = static_cast<TGraphErrors*>(nextG())) &&
	 (o = static_cast<TGraphErrors*>(nextO())) &&
	 (h = static_cast<TH1*>(nextH()))) {
    g->SetMarkerColor(h->GetMarkerColor());
    g->SetLineColor(h->GetLineColor());
    Add2Ratio(ratios, h, g);
    Add2Ratio(ratios, o, g);
    Printf("%s %s %s", g->GetName(), o->GetName(), h->GetName());
  }

  c->cd(2);
  yp /= (1-yp);
  
  TGraphErrors* band = new TGraphErrors(2);
  band->SetPoint(0,-3.95,1);
  band->SetPoint(1,5.95,1);
  band->SetPointError(0,0,.05);
  band->SetPointError(1,0,.05);
  band->SetFillColor(kYellow-10);
  band->SetFillStyle(1001);
  band->SetLineStyle(2);
  ratios->GetListOfGraphs()->AddFirst(band, "3l");

  ratios->Draw("ap");
  TAxis* yaxis = ratios->GetHistogram()->GetYaxis();
  TAxis* xaxis = ratios->GetHistogram()->GetXaxis();
  yaxis->SetTitle("ALICE/ATLAS");
  yaxis->SetLabelSize(yaxis->GetLabelSize()/yp);
  yaxis->SetTitleSize(yaxis->GetTitleSize()/yp);
  yaxis->SetTitleOffset(yaxis->GetTitleOffset()*yp);
  xaxis->SetTitle("#eta");
  xaxis->Set(200, -4, 6);
  xaxis->SetLabelSize(xaxis->GetLabelSize()/yp);
  xaxis->SetTitleSize(xaxis->GetTitleSize()/yp);
  band->Draw("same l");
  
  c->Modified();
  c->Update();
  c->cd();  

  c->Print("plots/alice_watlas_cp_v0a.pdf");
  c->Print("plots/alice_watlas_cp_v0a.C");
  c->Print("plots/alice_watlas_cp_v0a.png");
}

//====================================================================
void
DrawATLAS(Bool_t overlay=false)
{
  gROOT->SetMacroPath(Form("%s:%s", gROOT->GetMacroPath(),
			   "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/dndeta"));
  gROOT->LoadMacro("Drawer.C+g");
  gROOT->LoadMacro("ATLASpPb.C");

  DrawOrigATLAS(overlay);
  DrawAvgATLAS();
  
  DrawOrigATLASRatio(overlay);
  DrawAvgATLASRatio();
}
//====================================================================
//
// EOF
// 
