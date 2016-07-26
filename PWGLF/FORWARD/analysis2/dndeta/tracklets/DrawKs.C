const Color_t cc[] = { kMagenta+2,
		       kBlue+2,
		       kAzure-1, // 10,
		       kCyan+2,
		       kGreen+1,
		       kSpring+5,//+10,
		       kYellow+1,
		       kOrange+5,//+10,
		       kRed+1,
		       kPink+5,//+10,
		       kBlack };

TH1* GetCentK(TDirectory* top, Double_t c1, Double_t c2, Int_t s,
	      TLegend* l)
{
  TString dname; dname.Form("cent%06.2f_%06.2f", c1, c2);
  dname.ReplaceAll(".", "d");
  TDirectory* d = top->GetDirectory(dname);
  if (!d) {
    Warning("GetCetnK", "Directory %s not found in %s",
	    dname.Data(), top->GetName());
    return;
  }

  TDirectory* det = d->GetDirectory("details");
  if (!det) {
    Warning("GetCetnK", "Directory details not found in %s",
	    d->GetName());
    d->ls();
    return;
  }

  TObject* o = det->Get("scalar");
  if (!o) {
    Warning("GetCetnK", "Object scalar not found in %s",
	    det->GetName());
    return;
  }

  if (!o->IsA()->InheritsFrom(TH1::Class())) {
    Warning("GetCetnK", "Object %s is not a TH1, but a %s",
	    o->GetName(), o->ClassName());
    return;
  }
  TH1* h = static_cast<TH1*>(o->Clone());
  Color_t col = cc[(s-1)%10];
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->SetFillColor(col);
  h->SetFillStyle(1001);
  // h->SetTitle(Form("%5.2f-%5.2f%% #times %d", c1, c2, s));
  h->SetTitle(Form("%2.0f-%2.0f%% + %d", c1, c2, s-1));
  TF1* f = new TF1("", "[0]",-2.2,2.2);
  f->SetParameter(0,s-1);
  f->SetLineColor(col);
  f->SetLineStyle(7);
  f->SetLineWidth(1);
  // h->Scale(s);
  h->Add(f);
  h->GetListOfFunctions()->Add(f);
  f->SetParameter(0,s);
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    if (TMath::Abs(h->GetBinCenter(i)) > 2) {
      h->SetBinContent(i,0);
      h->SetBinError(i,0);
    }
  }
  
  TLegendEntry* e = l->AddEntry(h, h->GetTitle(), "f");
  e->SetFillColor(col);
  e->SetFillStyle(1001);
  e->SetLineColor(col);

  return h;
}

TCanvas* DrawKs(const char* filename)
{

  TFile* file = TFile::Open(filename, "READ");
  if (!file) {
    Warning("DrawKs", "File %s couldn't be opened", filename);
    return 0;
  }

  TH1* cent = static_cast<TH1*>(file->Get("cent"));
  if (!cent) {
    Warning("DrawKs", "Failed to find cent in %s", file->GetName());
    return 0;
  }

  TString t(filename);
  t.ReplaceAll("results/", "");
  t.ReplaceAll("combine_","");
  t.ReplaceAll("_0x3.root", "");
  TString nm(filename);
  nm.ReplaceAll(".root", "");
  nm.ReplaceAll("results/", "plots/");
  nm.ReplaceAll("combine", "ks");
  if (t.Contains("none")) 
    t.ReplaceAll("none", "No weights");
  else
    t.Append(" weights");
  
  Int_t cW = 1200;
  Int_t cH =  800;
  TCanvas* c = new TCanvas(nm,t,cW, cH);
  c->SetTopMargin(0.07);
  c->SetRightMargin(0.20);
  c->SetTicks();
  
  TLegend* l = new TLegend(1-c->GetRightMargin(),
			   c->GetBottomMargin(), 
			   .99,
			   1-c->GetTopMargin());
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  
  THStack* s = new THStack("ks", "");
  Int_t nCent = cent->GetXaxis()->GetNbins();
  for (Int_t i = 1; i <= nCent; i++) {
    Double_t c1 = cent->GetXaxis()->GetBinLowEdge(i);
    Double_t c2 = cent->GetXaxis()->GetBinUpEdge(i);

    TH1*     h  = GetCentK(file, c1, c2, nCent-i+1-2, l);
    if (!h) continue;

    s->Add(h);
  }
  s->Draw("nostack");
  TH1* f = s->GetHistogram();
  if (f) {
    f->SetXTitle("#eta");
    f->SetYTitle("#it{k}(#eta)");
  }
  
  TLatex* tit = new TLatex(0.55, 0.99, t);
  tit->SetTextFont(42);
  tit->SetTextAlign(23);
  tit->SetTextSize(0.03);
  tit->SetNDC();
  tit->Draw();
  l->SetBorderSize(0);
  l->Draw();

  c->Modified();
  c->Update();
  c->cd();

  return c;
}

void DrawKs(UShort_t flags, const char* var)
{
  TCanvas* c = DrawKs(Form("results/combine_%s_0x%x.root", var, flags));
  c->SaveAs(Form("plots/ks_%s_0x%x.png", var, flags));
}
