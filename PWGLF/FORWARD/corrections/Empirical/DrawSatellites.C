
TGraphErrors*
GetGraph(TDirectory* dir, Int_t low, Int_t high, TGraphErrors* avg=0)
{
  TObject* o = dir->Get(Form("graphErrors_cent_%d_%d", low, high));
  if (!o) return 0;
  TGraphErrors* g = static_cast<TGraphErrors*>(o);
  g->SetName(Form("g%02d%02d", low, high));
  g->SetTitle(Form("%2d-%2d%%", low, high));
  g->SetMarkerSize(2.2);

  if (!avg) return g;
  for (Int_t i = 0; i < g->GetN(); i++) {
    Double_t x   = g->GetX()[i];
    Double_t ex  = g->GetEX()[i];
    Double_t y   = g->GetY()[i];
    Double_t ey  = g->GetEY()[i];
    Double_t ay  = avg->GetY()[i];
    Double_t eay = avg->GetEY()[i];
    Double_t w   = (high-low); // 1/ey/ey;
    Printf("%3d: x=%7f y=%7f +/- %8f -> %9f",
	   i, x, y, ey, y*w);
    avg->SetPoint(i, x, y*w + ay);
    avg->SetPointError(i, ex, eay + w*ey);
  }
  
  return g;
}
TF1*
GetFunc(TDirectory* dir, Int_t low, Int_t high)
{
  TObject* o = dir->Get(Form("fitfunc_cent_%d_%d", low, high));
  if (!o) return 0;
  TF1* g = static_cast<TF1*>(o);
  g->SetName(Form("f%02d%02d", low, high));
  g->SetTitle(Form("%2d-%2d%%", low, high));
  return g;
}

TGraph*
GetRel(TGraphErrors* g) 
{
  TGraph* ret = new TGraph(g->GetN());
  ret->SetTitle(g->GetTitle());
  ret->SetName(Form("e%s", g->GetName()));
  ret->SetMarkerColor(g->GetMarkerColor());
  ret->SetMarkerStyle(g->GetMarkerStyle());
  ret->SetMarkerSize(g->GetMarkerSize());
  ret->SetFillColor(g->GetFillColor());
  ret->SetFillStyle(g->GetFillStyle());
  ret->SetLineColor(g->GetLineColor());
  ret->SetLineStyle(g->GetLineStyle());
  ret->SetLineWidth(g->GetLineWidth());

  for (Int_t i = 0; i < g->GetN(); i++) {
    Double_t x  = g->GetX()[i];
    Double_t ex = g->GetEX()[i];
    Double_t y  = g->GetY()[i];
    Double_t ey = g->GetEY()[i];
    Double_t r  = 100*ey/y;
    ret->SetPoint(i, x, r);
  }
  return ret;
}
	       

void
DrawSatellites(Bool_t doAvg=false)
{
  TFile *_file0 = TFile::Open("dndetaPreliminaryQM12.root");

  TCanvas* c = new TCanvas("c","c",1400,900);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);
  c->SetLeftMargin(0.08);
  c->Divide(1,2,0,0);

  TVirtualPad* p = c->cd(1);
  p->SetRightMargin(0.01);
  TMultiGraph* mg = new TMultiGraph;

  TGraphErrors* avg = 0;
  if (doAvg) {
    avg = new TGraphErrors(42);
    avg->SetTitle("Mean");
    avg->SetLineColor(kGray+2);
    avg->SetLineWidth(3);
    avg->SetLineStyle(2);
    avg->SetFillColor(0);
    avg->SetFillStyle(0);
  }
  
  mg->Add(GetGraph(_file0,  0,  5, avg));
  mg->Add(GetGraph(_file0,  5, 10, avg));
  mg->Add(GetGraph(_file0, 10, 20, avg));
  mg->Add(GetGraph(_file0, 20, 30, avg));

  if (avg) {
    for (Int_t i = 0; i < avg->GetN(); i++) {
      Double_t y  = avg->GetY()[i];
      Double_t ey = avg->GetEY()[i];
      avg->SetPoint(i, avg->GetX()[i], y/30); // ey);
      avg->SetPointError(i, avg->GetEX()[i], ey/30); // 1/TMath::Sqrt(ey));
    }
    mg->Add(avg, "c");
  }
  
  mg->Draw("a p 3");
  
  TLegend* l = p->BuildLegend(0.85,0.5,0.97,0.97);
  l->SetBorderSize(0);
  l->SetFillColor(0);
  
  mg->Draw("p");
  mg->GetHistogram()->SetMaximum(2100);
  mg->GetHistogram()->SetXTitle("#it{#eta}");
  mg->GetHistogram()->SetYTitle("1/#it{N} d#it{N}_{ch}/d#it{#eta}");
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  mg->GetHistogram()->GetYaxis()->SetNdivisions(10);
  
  TList* funcs = new TList();
  funcs->Add(GetFunc(_file0, 0, 5));
  funcs->Add(GetFunc(_file0, 5, 10));
  funcs->Add(GetFunc(_file0, 10, 20));
  funcs->Add(GetFunc(_file0, 20, 30));
  funcs->Draw("same");
  
  TLatex* ltx = new TLatex(.15, .97,
			   "Pb-Pb @ #sqrt{#it{s}_{NN}} = "
			   "2.76 TeV (satellites)");
  ltx->SetTextAlign(13);
  ltx->SetTextFont(42);
  ltx->SetNDC();
  ltx->Draw();
  
  p = c->cd(2);
  p->SetRightMargin(0.01);
  
  mg = new TMultiGraph;
  mg->Add(GetRel(GetGraph(_file0, 0, 5)));
  mg->Add(GetRel(GetGraph(_file0, 5, 10)));
  mg->Add(GetRel(GetGraph(_file0, 10, 20)));
  mg->Add(GetRel(GetGraph(_file0, 20, 30)));
  TGraph* eavg = 0;
  if (avg) mg->Add(eavg = GetRel(avg), "C");
  mg->Draw("a p");

  mg->GetHistogram()->SetMinimum(0);
  mg->GetHistogram()->SetMaximum(8.5);
  mg->GetHistogram()->SetXTitle("#it{#eta}");
  mg->GetHistogram()->SetYTitle("#delta(1/#it{N} d#it{N}_{ch}/d#it{#eta}) [%]");
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
  mg->GetHistogram()->GetYaxis()->SetNdivisions(10);
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(0.6);
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  mg->GetHistogram()->GetXaxis()->SetLabelSize(0.05);

  if (eavg) {
    Double_t min = 1000;
    Double_t max = 0;
    Double_t sum = 0;
    Int_t    cnt = 0;
    for (Int_t i = 0; i < eavg->GetN(); i++) {
      Double_t x = eavg->GetX()[i];
      if (x < -3.5 || x > 5 || TMath::Abs(x) < 2) continue;
      Double_t y = eavg->GetY()[i];
      min        = TMath::Min(min,y);
      max        = TMath::Max(max,y);
      sum        += y;
      cnt++;
    }	
    Printf("%f - %f %f", min, max, sum/cnt);
  }
  c->Modified();
  c->Update();
  c->cd();
   
  
  c->Print("dndeta_satellites.pdf");
}

