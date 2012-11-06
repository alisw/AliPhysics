void 
DrawMG(TObject* o, TLegend* l, TLegend* cent=0)
{
  TMultiGraph* mg = static_cast<TMultiGraph*>(o);
  TIter next(mg->GetListOfGraphs());
  TObject* g = 0;
  Int_t i = 0;
  while ((g = next())) {
    g->Draw("p same");
    TGraph* gg = static_cast<TGraph*>(g);
    TLegendEntry* e = 0;
    if (i == 0) {
      e = l->AddEntry("dummy", mg->GetTitle(), "pl");
      e->SetMarkerColor(kBlack);
      e->SetMarkerStyle(gg->GetMarkerStyle());
      e->SetFillStyle(0);
    }
    i++;
    if (!cent) continue;
    e = cent->AddEntry("dummy", g->GetTitle(), "f");
    e->SetMarkerColor(kBlack);
    e->SetFillStyle(1001);
    e->SetFillColor(gg->GetMarkerColor());
  }

}

void 
DrawEmpirical(const char* corrs="EmpiricalCorrection.root")
{
  TFile* fcorr = TFile::Open(corrs, "READ");
  if (!fcorr) return;

  TObject* odfmdnom  = fcorr->Get("nominal");
  TObject* odfmdsat  = fcorr->Get("fmdfmd/satellite");
  TObject* odfullsat = fcorr->Get("fmdfull/satellite");
  TObject* ocfmdfmd  = fcorr->Get("fmdfmd/average");
  TObject* ocfmdfull = fcorr->Get("fmdfull/average");
  TObject* oraverage = fcorr->Get("ratio");

  Int_t size = 800;
  TCanvas* c = new TCanvas("c", "c", size / TMath::Sqrt(2), size);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);
  c->Divide(1, 3, 0, 0);

  TVirtualPad* p =  c->cd(1);
  p->SetRightMargin(0.01);
  p->SetGridx();
  TH1* h = new TH1D("h1", "h", 200, -4, 6);
  h->SetStats(0);
  h->SetXTitle("#eta");
  h->SetYTitle("1/N_{ev} dN_{ch}/d#eta");
  h->SetTitle("");
  h->SetMinimum(1);
  h->SetMaximum(2100);
  h->DrawClone();
  TLegend* l  = new TLegend(.2,  .05, .45, .25, "dN/d#eta");
  TLegend* lc = new TLegend(.55, .05, .8,  .25, "Centrality");
  DrawMG(odfmdnom, l);
  DrawMG(odfmdsat, l);
  DrawMG(odfullsat, l, lc);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->Draw();
  lc->SetFillStyle(0);
  lc->SetBorderSize(0);
  lc->Draw();

  p = c->cd(2);
  p->SetGridx();
  p->SetRightMargin(0.01);
  h->SetMinimum(.91);
  h->SetMaximum(1.21);
  h->SetYTitle("#LTdN_{ch}/d#eta|_{nominal}/dN_{ch}/d#eta|_{satelitte}#GT");
  h->DrawClone();
  // DrawMG(ocfmdfmd);
  // DrawMG(ocfmdfull);
  ocfmdfmd->Draw("p same");
  ocfmdfull->Draw("p same");
  l = new TLegend(.2, .05, .8, .2);
  l->SetNColumns(2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(ocfmdfmd,  "FMD/FMD",  "pl");
  l->AddEntry(ocfmdfull, "Full/FMD", "pl");
  l->Draw();

  p = c->cd(3);
  p->SetGridx();
  p->SetRightMargin(0.01);
  h->SetMinimum(-.035);
  h->SetMaximum(+.065);
  h->SetYTitle("#frac{#LTFull/FMD#GT}{#LTFMD/FMD#GT}");
  h->DrawClone();
  oraverage->Draw("p same");
  TGraph* g = static_cast<TGraph*>(oraverage);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.1);
  // DrawMG(orfmdfmd);
  // DrawMG(orfmdfull);

  c->SaveAs("EmpericalCorrection.png"); 

}

