void PresentationStyle() {
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(101);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetPaperSize(TStyle::kA4);
  gStyle->SetTitleOffset(1.2);

  Int_t error = 0;
  gROOT->Macro("TDR_style.C", &error);
  gStyle->SetOptDate(0);
  gStyle->SetTitleColor(1);
  gStyle->SetMarkerStyle(kFullCircle);
}
