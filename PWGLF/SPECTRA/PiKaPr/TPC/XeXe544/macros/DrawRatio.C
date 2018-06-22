//
//calculates pos-particle to neg-particle ratio
//
TH1* Get(const Char_t* particleType, Int_t cent)
{
  TVirtualPad* currpad = gPad;
  TFile f(Form("FinalSpectra_forQM/canvPtDist%s.root", particleType), "READ");
  //Printf("%s",f.GetName());
  TCanvas* c = (TCanvas*)f.Get(Form("canvPtDist%s", particleType));
  //Printf("%s",c->GetName());
  TList* l = c->GetListOfPrimitives();
  l->ls();
  TH1* h = (TH1*)l->At(cent)->Clone(Form("%s%i", particleType, cent));
  h->SetDirectory(0);
  f.Close();
  currpad->cd();
  return h;
}

void SetColor(TH1* h, Int_t cent)
{
  const TString col[10] = { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a" }; //http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10
  h->SetLineColor(TColor::GetColor(col[cent]));
  h->SetMarkerColor(TColor::GetColor(col[cent]));
}

const Char_t* GetCent(Int_t cent)
{
  Float_t min = 0, width = 0;
  for (Int_t i = 0; i <= cent; i++) {
    if (i < 2)
      width = 5;
    if (i >= 2 && i <= 7)
      width = 10;
    else if (i > 7)
      width = 20;
    if (cent == i)
      break;
    min += width;
  }

  return Form("%.0f-%.0f", min, min + width);
}

TH1* DoRatio(const Char_t* particleType, Int_t cent)
{
  TH1* hnum = Get(particleType, cent);
  //Printf("%s", hnum->GetName());
  TH1* hden = Get(Form("Anti%s", particleType), cent);
  //Printf("%s", hden->GetName());
  hnum->Divide(hden);
  hnum->SetMaximum(1.3);
  hnum->SetMinimum(0.7);
  hnum->SetDirectory(0);
  hnum->SetName(Form("Ratio%s%i", particleType, cent));
  hnum->SetTitle(Form("%s%%;#it{p}_{T} (GeV/#it{c});%s/Anti%s", GetCent(cent), particleType, particleType));
  SetColor(hnum, cent);
  return hnum;
}

TCanvas* DrawRatioPart(const Char_t* particleType, Int_t min = 0, Int_t max = 9)
{
  TCanvas* canv = new TCanvas(Form("c%s", particleType), particleType);
  for (Int_t i = min; i < max; i++) {
    TH1* h= Get(particleType, i);
    h->SetLineColor(kRed);
    h->Draw(i == 0 ? "" : "same");
    h= Get(Form("Anti%s", particleType), i);
    h->SetLineColor(kBlue);
    h->Draw("EPsame");
  }
  
  canv->BuildLegend();
  TCanvas* canvRatio = new TCanvas(Form("cRatio%s", particleType), particleType);
  for (Int_t i = min; i < max; i++)
    DoRatio(particleType, i)->Draw("EPsame");
  canvRatio->BuildLegend(0.75,0.6,0.9,0.9);
  TLatex* label = new TLatex(.8, .91, particleType);
  label->SetNDC();
  label->Draw();
  canvRatio->SetGridy();
  canvRatio->SaveAs(Form("Part_antiPart_Ratios/part_antipart_ratio_%s.png", particleType));
  canvRatio->SaveAs(Form("Part_antiPart_Ratios/part_antipart_ratio_%s.root", particleType));
  return canvRatio;
}

void DrawRatio()
{
  gStyle->SetOptTitle(0);
  DrawRatioPart("Pion");
  DrawRatioPart("Kaon");
  DrawRatioPart("Proton");
}
