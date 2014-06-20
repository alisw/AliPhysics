THStack* DrawOne(TVirtualPad* p, 
		 Double_t     yr,
		 Bool_t       top,
		 TDirectory*  dir, 
		 const char*  name)
{
  p->cd();
  p->SetFillColor(0);
  p->SetFillStyle(0);
  p->SetLineColor(0);
  p->SetRightMargin(0.01);
  p->SetLeftMargin(0.12);
  p->SetGridx();
  if (top) p->SetBottomMargin(0.001);
  else     p->SetBottomMargin(0.2);
  if (top) p->SetTopMargin(0.02);
  else     p->SetTopMargin(0.0001);
  
  
  THStack* s = static_cast<THStack*>(dir->Get(name));
  s->Draw("nostack");
  Double_t sc = (top ? 1-yr : yr);
  TAxis* ya = s->GetHistogram()->GetYaxis();
  ya->SetLabelSize(1/sc*ya->GetLabelSize());
  ya->SetTitleSize(1/sc*ya->GetTitleSize());
  ya->SetTitleOffset(sc*(ya->GetTitleOffset()+.5));
  ya->SetTitleFont(42);
  ya->SetLabelFont(42);
  TAxis* xa = s->GetHistogram()->GetXaxis();
  xa->SetLabelSize(!top ? 1/yr*xa->GetLabelSize() : 0);
  xa->SetTitleSize(!top ? 1/yr*xa->GetTitleSize() : 0);
  xa->SetTitleOffset(yr*(xa->GetTitleOffset()+2));
  xa->SetTitleFont(42);
  xa->SetLabelFont(42);

  p->Modified();
  p->Update();
  p->cd();

  return s;
}

void DrawEmpirical(const char* filename="Empirical.root", 
		   Bool_t fmd=true)
{
  gStyle->SetOptTitle(0);

  TFile* file = TFile::Open(filename, "READ");
  if (!file) return;

  Double_t yr = 0.3;
  TCanvas* c  = new TCanvas("c","c", 1000,1000);
  TPad*    p1 = new TPad("p1","p1",0,0,1,yr);
  TPad*    p2 = new TPad("p2","p2",0,yr,1,1);
  c->cd(); p1->Draw();
  c->cd(); p2->Draw();
  
  gDirectory->cd("Forward");
  THStack* r = DrawOne(p1, yr, false, gDirectory, "ratios");  
  THStack* e = DrawOne(p2, yr, true, gDirectory, "empirical");

  r->SetMinimum(0.945);
  r->SetMaximum(1.055);
  r->GetXaxis()->SetTitle("#it{#eta}");
  r->GetYaxis()->SetTitle("Ratio to mean");
  e->SetMinimum(0.005);
  e->GetYaxis()->SetTitle("#it{E_{c}}(#it{#eta})");
  TIter nextE(e->GetHists());
  TIter nextR(r->GetHists());
  TH1*  hist = 0;
  Color_t cols[]  = { kRed+2, kGreen+2, kBlue+2, kMagenta+2, 0 };
  Color_t *ptr    = cols;
  Style_t stys[]  = { 20, 21, 22, 23 };
  Style_t* sty    = stys;
  while (*ptr) { 
    hist = static_cast<TH1*>(nextE()); 
    hist->SetMarkerColor(*ptr);
    hist->SetMarkerSize(2);
    hist->SetMarkerStyle(*sty);
    hist = static_cast<TH1*>(nextR()); 
    hist->SetMarkerColor(*ptr);
    hist->SetMarkerSize(2);
    hist->SetMarkerStyle(*sty);
    ptr++;
    sty++;
  }


  TLegend* l = p2->BuildLegend(0.35, .2, .65, .8);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);

  c->Modified();
  c->Update();
  c->cd();
  c->Print("empirical.png");
  
}
  
