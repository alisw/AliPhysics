// -*- C++ -*-

TGraph* MakeGraph(TH1* h) {
    TGraphErrors *g = new TGraphErrors;
    g->SetMarkerStyle(kFullDotMedium);
    g->SetMarkerColor(h->GetLineColor());
    g->SetLineColor(h->GetLineColor());
    for (Int_t i=0,j=0,n=hrA->GetNbinsX();i<n; ++i) {
      j=(n+i-1)%n;
      g->SetPoint(i,
                  (3564 + atoi(h->GetXaxis()->GetBinLabel(1+i)) - atoi(h->GetXaxis()->GetBinLabel(1+j)))%3564,
                  h->GetBinContent(1+i));
      g->SetPointError(i,
                       0,
                       h->GetBinError(1+i));
    }
    return g;
}

void PlotSingle(TString rn, TString title) {
  TFile::Open(rn);
  TGraph *gA = MakeGraph(hrA);
  TGraph *gC = MakeGraph(hrC);
  TCanvas *c1 = new TCanvas;
  TH1 *hf = new TH1D("hf", title+" one-arm ratios;distance to last collision (BC)", 3564,0,3564);
  hf->SetStats(0);
  hf->SetMinimum(0);
  hf->SetMaximum(1.1*TMath::Max(TMath::MaxElement(gA->GetN(), gA->GetY()),
                                TMath::MaxElement(gC->GetN(), gC->GetY())));
  hf->Draw();
  gPad->SetLogx();
  gA->Draw("PE");
  gC->Draw("PE");
  TLegend *leg = new TLegend(0.75, title.Contains("V0") ? 0.65 : 0.25,
                             0.9,  title.Contains("V0") ? 0.80 : 0.4);
  leg->AddEntry(gA, "r_{A}", "PE");
  leg->AddEntry(gC, "r_{C}", "PE");
  leg->Draw();

  TString pn = rn;
  pn.ReplaceAll("root", "png");
  c1->SaveAs(pn);
}

void plotRatioVsDistance()
{
  PlotSingle("pileup_fill5533_c2VBAandVBC.root", "fill 5533 V0");
  PlotSingle("pileup_fill5533_c2UBAandUBC.root", "fill 5533 AD");
  PlotSingle("pileup_fill5533_c2TVX.root", "fill 5533 T0");

  PlotSingle("pileup_fill5568_c2VBAandVBC.root", "fill 5568 V0");
  PlotSingle("pileup_fill5568_c2UBAandUBC.root", "fill 5568 AD");
  PlotSingle("pileup_fill5568_c2TVX.root", "fill 5568 T0");
}
