{
  TCanvas *canvasT = new TCanvas("canvasT", "canvasT", 900, 500);
  canvasT->Divide(1,1);
  canvasT->cd(1);
  TH1F *his1 = new TH1F("his1.class(Error).style(line-color:#f30000;)", "his1", 100, -5, 5);
  his1->FillRandom("gaus", 15000);
  his1->SetStats(kFALSE);
  TH1D *his2 = new TH1D("his2.class(Error)", "his2", 100, -5, 5);
  his2->FillRandom("gaus", 10000);
  TH1I *his3 = new TH1I("his3.class(Error)", "his3", 100, -5, 5);
  his3->FillRandom("gaus", 5000);
  his1->Draw();
  his2->Draw("same");
  his3->Draw("same");
  gStyle->SetOptTitle(0);
  TPaveText *pt = new TPaveText(-0.438183, 694.575009, 0.438183, 740.053135);
  pt->AddText("Example of styling with using AliDrawStyle");
  pt->SetTextSize(0.04);
  pt->SetShadowColor(0);
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->Draw();
  gPad->BuildLegend();
}